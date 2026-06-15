// Standalone Python binding for slim decompression: vendored slim C++ + libc
// only (no libactpol/zzip/getdata). Decodes straight from an in-memory buffer
// via fmemopen, no temp files.

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <unistd.h>
#include <sys/mman.h>
#include <cstring>
#include <cerrno>
#include <cstdlib>
#include <vector>
#include <thread>
#include <atomic>

#include "_vendor_slim/slim.h"

// Decode a slim payload directly into a caller-provided buffer (no malloc, no
// extra copy). Python/GIL-free so it can run on a worker thread.
static bool decode_into(const char *data, Py_ssize_t len,
                        char *dst, size_t dst_cap,
                        size_t *outlen, const char **errmsg) {
    *outlen = 0;

    FILE *fp = fmemopen((void *)data, (size_t)len, "rb");
    if (fp == NULL) { *errmsg = "fmemopen failed"; return false; }

    size_t written = 0;
    bool ok = true;
    // slim_expander_t fcloses fp in its destructor; do not close it here.
    try {
        slim_expander_t exp(fp, (size_t)len);
        if (!exp.is_open()) {
            *errmsg = "could not open slim stream";
            ok = false;
        } else {
            size_t raw_size = exp.get_rawsize();
            if (raw_size != 0 && raw_size > dst_cap) {
                *errmsg = "destination buffer too small for raw size";
                ok = false;
            } else {
                size_t want = raw_size ? raw_size : dst_cap;
                while (written < want) {
                    size_t got = exp.read((unsigned char *)(dst + written),
                                          want - written);
                    if (got == 0) break;
                    written += got;
                }
            }
        }
    } catch (const char *s) {
        *errmsg = s; ok = false;
    } catch (...) {
        *errmsg = "unknown slim decode error"; ok = false;
    }

    if (!ok) return false;
    *outlen = written;
    return true;
}

// Like decode_into, but mallocs its own output (for decompress / many).
static bool decode_one(const char *data, Py_ssize_t len,
                       char **out, size_t *outlen, const char **errmsg) {
    *out = NULL;
    *outlen = 0;

    FILE *fp = fmemopen((void *)data, (size_t)len, "rb");
    if (fp == NULL) { *errmsg = "fmemopen failed"; return false; }
    char *outbuf = NULL;
    size_t raw_size = 0;
    bool ok = true;
    try {
        slim_expander_t exp(fp, (size_t)len);
        if (!exp.is_open()) { *errmsg = "could not open slim stream"; ok = false; }
        else {
            raw_size = exp.get_rawsize();
            if (raw_size == 0) {
                size_t cap = 1 << 20, used = 0;
                outbuf = (char *)malloc(cap);
                if (!outbuf) { *errmsg = "out of memory"; ok = false; }
                while (ok) {
                    if (used == cap) {
                        cap *= 2;
                        char *nb = (char *)realloc(outbuf, cap);
                        if (!nb) { *errmsg="out of memory"; ok=false; break; }
                        outbuf = nb;
                    }
                    size_t got = exp.read((unsigned char *)(outbuf + used), cap - used);
                    if (got == 0) break;
                    used += got;
                }
                raw_size = used;
            } else {
                outbuf = (char *)malloc(raw_size ? raw_size : 1);
                if (!outbuf) { *errmsg = "out of memory"; ok = false; }
                size_t total = 0;
                while (ok && total < raw_size) {
                    size_t got = exp.read((unsigned char *)(outbuf + total),
                                          raw_size - total);
                    if (got == 0) break;
                    total += got;
                }
                raw_size = total;
            }
        }
    } catch (const char *s) { *errmsg = s; ok = false; }
      catch (...) { *errmsg = "unknown slim decode error"; ok = false; }

    if (!ok) { free(outbuf); return false; }
    *out = outbuf;
    *outlen = raw_size;
    return true;
}

static PyObject *actslim_decompress(PyObject *self, PyObject *args) {
    Py_buffer in;
    if (!PyArg_ParseTuple(args, "y*", &in))
        return NULL;

    char *outbuf = NULL;
    size_t outlen = 0;
    const char *errmsg = NULL;
    bool ok;

    Py_BEGIN_ALLOW_THREADS
    ok = decode_one((const char *)in.buf, in.len, &outbuf, &outlen, &errmsg);
    Py_END_ALLOW_THREADS

    PyBuffer_Release(&in);
    if (!ok) {
        PyErr_SetString(PyExc_ValueError, errmsg ? errmsg : "slim decode failed");
        return NULL;
    }
    PyObject *result = PyBytes_FromStringAndSize(outbuf, (Py_ssize_t)outlen);
    free(outbuf);
    return result;
}

// decompress_many(payloads, nthreads=0) -> list[bytes]
static PyObject *actslim_decompress_many(PyObject *self, PyObject *args) {
    PyObject *seq;
    int nthreads = 0;
    if (!PyArg_ParseTuple(args, "O|i", &seq, &nthreads))
        return NULL;

    PyObject *fast = PySequence_Fast(seq, "payloads must be a sequence");
    if (!fast) return NULL;
    Py_ssize_t n = PySequence_Fast_GET_SIZE(fast);

    // Acquire buffer views for every input (under the GIL).
    std::vector<Py_buffer> bufs((size_t)n);
    Py_ssize_t acquired = 0;
    for (; acquired < n; acquired++) {
        PyObject *item = PySequence_Fast_GET_ITEM(fast, acquired);
        if (PyObject_GetBuffer(item, &bufs[(size_t)acquired], PyBUF_SIMPLE) != 0)
            break;
    }
    if (acquired != n) {
        for (Py_ssize_t i = 0; i < acquired; i++) PyBuffer_Release(&bufs[(size_t)i]);
        Py_DECREF(fast);
        return NULL;
    }

    std::vector<char *> outs((size_t)n, NULL);
    std::vector<size_t> outlens((size_t)n, 0);
    std::vector<const char *> errs((size_t)n, NULL);

    if (nthreads <= 0) {
        unsigned hc = std::thread::hardware_concurrency();
        nthreads = hc ? (int)hc : 1;
    }
    if (nthreads > n) nthreads = (int)(n > 0 ? n : 1);

    std::atomic<Py_ssize_t> next(0);

    Py_BEGIN_ALLOW_THREADS
    auto worker = [&]() {
        for (;;) {
            Py_ssize_t i = next.fetch_add(1);
            if (i >= n) break;
            decode_one((const char *)bufs[(size_t)i].buf, bufs[(size_t)i].len,
                       &outs[(size_t)i], &outlens[(size_t)i], &errs[(size_t)i]);
        }
    };
    if (nthreads <= 1) {
        worker();
    } else {
        std::vector<std::thread> pool;
        pool.reserve((size_t)nthreads);
        for (int t = 0; t < nthreads; t++) pool.emplace_back(worker);
        for (auto &th : pool) th.join();
    }
    Py_END_ALLOW_THREADS

    for (Py_ssize_t i = 0; i < n; i++) PyBuffer_Release(&bufs[(size_t)i]);
    Py_DECREF(fast);

    PyObject *result = PyList_New(n);
    if (!result) { for (Py_ssize_t i=0;i<n;i++) free(outs[(size_t)i]); return NULL; }
    for (Py_ssize_t i = 0; i < n; i++) {
        if (outs[(size_t)i] == NULL && errs[(size_t)i] != NULL) {
            Py_DECREF(result);
            PyErr_Format(PyExc_ValueError, "slim decode failed on item %zd: %s",
                         (Py_ssize_t)i, errs[(size_t)i]);
            for (Py_ssize_t j = 0; j < n; j++) free(outs[(size_t)j]);
            return NULL;
        }
        PyObject *b = PyBytes_FromStringAndSize(outs[(size_t)i],
                                                (Py_ssize_t)outlens[(size_t)i]);
        free(outs[(size_t)i]);
        outs[(size_t)i] = NULL;
        if (!b) { Py_DECREF(result); for (Py_ssize_t j=i+1;j<n;j++) free(outs[(size_t)j]); return NULL; }
        PyList_SET_ITEM(result, i, b);
    }
    return result;
}

// decompress_into(payloads, out, item_nbytes, nthreads=0) -> None
// Decode each payload into out[i*item_nbytes:], in parallel. The hot path: no
// per-channel allocation, scales like the original OpenMP decode.
static PyObject *actslim_decompress_into(PyObject *self, PyObject *args) {
    PyObject *seq;
    Py_buffer out;
    Py_ssize_t item_nbytes;
    int nthreads = 0;
    if (!PyArg_ParseTuple(args, "Ow*n|i", &seq, &out, &item_nbytes, &nthreads))
        return NULL;

    PyObject *fast = PySequence_Fast(seq, "payloads must be a sequence");
    if (!fast) { PyBuffer_Release(&out); return NULL; }
    Py_ssize_t n = PySequence_Fast_GET_SIZE(fast);

    if (out.len < n * item_nbytes) {
        PyErr_SetString(PyExc_ValueError, "output buffer too small");
        Py_DECREF(fast); PyBuffer_Release(&out); return NULL;
    }

    std::vector<Py_buffer> bufs((size_t)n);
    Py_ssize_t acquired = 0;
    for (; acquired < n; acquired++) {
        PyObject *item = PySequence_Fast_GET_ITEM(fast, acquired);
        if (PyObject_GetBuffer(item, &bufs[(size_t)acquired], PyBUF_SIMPLE) != 0)
            break;
    }
    if (acquired != n) {
        for (Py_ssize_t i = 0; i < acquired; i++) PyBuffer_Release(&bufs[(size_t)i]);
        Py_DECREF(fast); PyBuffer_Release(&out);
        return NULL;
    }

    std::vector<const char *> errs((size_t)n, NULL);
    std::vector<size_t> wrote((size_t)n, 0);
    char *base = (char *)out.buf;

    if (nthreads <= 0) {
        unsigned hc = std::thread::hardware_concurrency();
        nthreads = hc ? (int)hc : 1;
    }
    if (nthreads > n) nthreads = (int)(n > 0 ? n : 1);

    std::atomic<Py_ssize_t> next(0);
    Py_BEGIN_ALLOW_THREADS
    auto worker = [&]() {
        for (;;) {
            Py_ssize_t i = next.fetch_add(1);
            if (i >= n) break;
            decode_into((const char *)bufs[(size_t)i].buf, bufs[(size_t)i].len,
                        base + i * item_nbytes, (size_t)item_nbytes,
                        &wrote[(size_t)i], &errs[(size_t)i]);
        }
    };
    if (nthreads <= 1) worker();
    else {
        std::vector<std::thread> pool;
        pool.reserve((size_t)nthreads);
        for (int t = 0; t < nthreads; t++) pool.emplace_back(worker);
        for (auto &th : pool) th.join();
    }
    Py_END_ALLOW_THREADS

    for (Py_ssize_t i = 0; i < n; i++) PyBuffer_Release(&bufs[(size_t)i]);
    Py_DECREF(fast);
    PyBuffer_Release(&out);

    for (Py_ssize_t i = 0; i < n; i++) {
        if (errs[(size_t)i] != NULL) {
            PyErr_Format(PyExc_ValueError, "slim decode failed on item %zd: %s",
                         (Py_ssize_t)i, errs[(size_t)i]);
            return NULL;
        }
        if (wrote[(size_t)i] != (size_t)item_nbytes) {
            PyErr_Format(PyExc_ValueError,
                         "item %zd expanded to %zu bytes, expected %zd",
                         (Py_ssize_t)i, wrote[(size_t)i], (Py_ssize_t)item_nbytes);
            return NULL;
        }
    }
    Py_RETURN_NONE;
}

static PyMethodDef methods[] = {
    {"decompress", actslim_decompress, METH_VARARGS,
     "decompress(slm_bytes) -> bytes. Expand raw slim-compressed bytes."},
    {"decompress_many", actslim_decompress_many, METH_VARARGS,
     "decompress_many(payloads, nthreads=0) -> list[bytes]. Decode many slim\n"
     "payloads in parallel using an internal C thread pool (GIL released)."},
    {"decompress_into", actslim_decompress_into, METH_VARARGS,
     "decompress_into(payloads, out, item_nbytes, nthreads=0) -> None.\n"
     "Decode payloads directly into one preallocated writable buffer."},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_actslim",
    "Standalone slim decompressor (no libactpol/zzip).", -1, methods,
};

PyMODINIT_FUNC PyInit__actslim(void) {
    return PyModule_Create(&moduledef);
}
