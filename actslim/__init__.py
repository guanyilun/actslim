"""actslim -- standalone reader for ACT slim-compressed dirfile TODs.

No libactpol, no zzip, no getdata, no autotools.  Just the vendored slim
codec (compiled into ``actslim._actslim``) plus the Python standard library.

Typical use::

    import actslim
    tod = actslim.read_zip("/path/to/1572374891....ar6.zip",
                           channels=["tesdatar00c01", ...])  # or all tes by default
    tod["tesdatar00c01"]   # -> numpy int32 array

The ``.slm`` payload expands to the *original* raw channel bytes; the numpy
dtype for each channel comes from the dirfile ``format`` text file, exactly as
getdata's RAW entries describe it.
"""

import ctypes
import ctypes.util
import mmap
import struct
import threading
import zipfile

import numpy as np

from ._actslim import decompress, decompress_many, decompress_into

__all__ = ["decompress", "decompress_many", "decompress_into", "read_zip",
           "read_dirfile", "read_zip_array", "read_fields", "parse_format",
           "parse_lincom", "FORMAT_DTYPES", "TOD", "read_tod"]

# dirfile RAW type code -> numpy dtype.  Authoritative source: libactpol
# getdata.c (sizes at lines 337-347, signedness at 824-855):
#   c=1B  s=int16  u=uint16  S=i=int32(signed,4B)  U=uint32  f=float32  d=float64
FORMAT_DTYPES = {
    "c": np.uint8,
    "s": np.int16,    "u": np.uint16,
    "S": np.int32,    "U": np.uint32,    "i": np.int32,
    "f": np.float32,  "d": np.float64,
}


def parse_format(text):
    """Parse a dirfile ``format`` file.

    Returns ``{field_name: (spf, numpy_dtype)}`` for every RAW field.
    Non-RAW entries (bit, lincom, ...) are ignored -- they are derived fields,
    not stored channels.
    """
    fields = {}
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        # RAW form:  <name> RAW <typecode> <spf>
        if len(parts) >= 4 and parts[1].upper() == "RAW":
            name, _, tcode, spf = parts[0], parts[1], parts[2], parts[3]
            dt = FORMAT_DTYPES.get(tcode)
            if dt is None:
                continue
            try:
                spf = int(spf)
            except ValueError:
                continue
            fields[name] = (spf, dt)
    return fields


def parse_lincom(text):
    """Parse LINCOM (derived) fields from a dirfile ``format`` file.

    Returns ``{name: [(in_field, scale, offset), ...]}``.  A LINCOM field is
    ``sum_i (scale_i * in_field_i + offset_i)`` -- this is how the dirfile
    encodes e.g. encoder counts -> degrees and cpu_s/cpu_us -> ctime.
    """
    out = {}
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        parts = line.split()
        if len(parts) >= 3 and parts[1].upper() == "LINCOM":
            try:
                nterm = int(parts[2])
                terms = []
                for k in range(nterm):
                    f = parts[3 + 3 * k]
                    m = float(parts[4 + 3 * k])
                    b = float(parts[5 + 3 * k])
                    terms.append((f, m, b))
            except (ValueError, IndexError):
                continue
            out[parts[0]] = terms
    return out


_LFH = struct.Struct("<IHHHHHIIIHH")  # local file header (zip spec, 30 bytes)

# Sequential prefetch tuning.  The reader thread pulls the file into the page
# cache in pieces of this size.  A decode batch waits for at least _BATCH
# channels, then also takes every later channel that has already been read.
_PREFETCH_PIECE = 32 << 20
_BATCH = 64
# Prefetch only when the requested entries make up at least this fraction of
# the byte range that spans them; sparser reads are cheaper as random access.
_PREFETCH_MIN_FRACTION = 0.25
# Skip the prefetch when at least this fraction of the range is already in
# the page cache: there is nothing to gain and copying it costs ~10%.
_CACHED_FRACTION = 0.9


def _cached_fraction(fd, start, stop):
    """Fraction of ``[start, stop)`` of ``fd`` resident in the page cache.

    Uses mincore(2) on a throwaway mapping; returns None where that isn't
    available (non-Linux, odd filesystems), so callers just prefetch.
    """
    try:
        libc = ctypes.CDLL(ctypes.util.find_library("c"), use_errno=True)
        libc.mmap.restype = ctypes.c_void_p
        libc.mmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int,
                              ctypes.c_int, ctypes.c_int, ctypes.c_long]
        libc.mincore.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p]
        libc.munmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
        page = mmap.PAGESIZE
        off = start - start % page
        length = stop - off
        if length <= 0:
            return None
        addr = libc.mmap(None, length, mmap.PROT_READ, mmap.MAP_SHARED, fd, off)
        if addr in (None, ctypes.c_void_p(-1).value):
            return None
        try:
            vec = np.zeros((length + page - 1) // page, dtype=np.uint8)
            if libc.mincore(addr, length, vec.ctypes.data) != 0:
                return None
            return float(np.count_nonzero(vec & 1)) / len(vec)
        finally:
            libc.munmap(addr, length)
    except Exception:
        return None


class _Prefetcher:
    """Read ``[start, stop)`` of a file front to back in a background thread.

    The bytes are discarded: the point is to fill the page cache with one
    sequential pass, so decoding from the mmap never waits on a seek.  On a
    spinning disk this runs at ~2x the speed of letting 40 decode threads fault
    pages in concurrently (which turns into random I/O), and it lets decoding
    of early channels overlap with reading of later ones.

    ``wait(pos)`` blocks until everything before ``pos`` has been read.  If the
    reader stops early (error, short file) ``wait`` returns anyway -- the mmap
    still works, it just faults the remaining pages in itself.
    """

    def __init__(self, path, start, stop):
        self._pos = start
        self._done = False
        self._cancel = False
        self._cv = threading.Condition()
        self._thread = threading.Thread(target=self._run,
                                        args=(path, start, stop), daemon=True)
        self._thread.start()

    def _run(self, path, start, stop):
        try:
            buf = memoryview(bytearray(min(_PREFETCH_PIECE, max(stop - start, 1))))
            with open(path, "rb", buffering=0) as fh:
                fh.seek(start)
                pos = start
                while pos < stop and not self._cancel:
                    n = fh.readinto(buf[:min(len(buf), stop - pos)])
                    if not n:
                        break
                    pos += n
                    with self._cv:
                        self._pos = pos
                        self._cv.notify_all()
        except Exception:
            pass
        finally:
            with self._cv:
                self._done = True
                self._cv.notify_all()

    def wait(self, pos):
        with self._cv:
            self._cv.wait_for(lambda: self._done or self._pos >= pos)

    def ready(self, pos):
        """True if ``wait(pos)`` would return without blocking."""
        with self._cv:
            return self._done or self._pos >= pos

    def close(self):
        self._cancel = True
        self._thread.join()


class _Payloads:
    """The ``.slm`` payloads of the requested channels, in channel order.

    For *Stored* entries (the ACT case) payloads are zero-copy memoryview
    slices of an mmap of the zip; otherwise they are ``bytes`` from
    ``z.read``.  Use :meth:`batches` to decode while a prefetch is running,
    or :meth:`all` to get every payload at once.
    """

    def __init__(self, mm, infos, prefetch, payloads=None):
        self._mm, self._infos, self._prefetch = mm, infos, prefetch
        self._payloads = payloads
        self._views = []
        self._spans = {}

    def __len__(self):
        return len(self._infos)

    def _span(self, i):
        """``(start, stop)`` of payload ``i`` in the file, from its local header."""
        span = self._spans.get(i)
        if span is None:
            info = self._infos[i]
            off = info.header_offset
            if self._prefetch is not None:
                self._prefetch.wait(off + 30)
            sig, _, _, _, _, _, _, _, _, nlen, elen = _LFH.unpack(self._mm[off:off + 30])
            if sig != 0x04034b50:
                raise ValueError("bad local header for %s" % info.filename)
            start = off + 30 + nlen + elen
            span = self._spans[i] = (start, start + info.compress_size)
        return span

    def _view(self, i):
        """Payload ``i``, waiting for the prefetch to reach it first."""
        span = self._span(i)
        if self._prefetch is not None:
            self._prefetch.wait(span[1])
        v = memoryview(self._mm)[span[0]:span[1]]
        self._views.append(v)
        return v

    def all(self):
        if self._payloads is not None:
            return self._payloads
        return [self._view(i) for i in range(len(self._infos))]

    def _ready(self, i):
        """True if payload ``i`` can be sliced without waiting on the prefetch."""
        pf = self._prefetch
        if pf is None:
            return True
        if i not in self._spans and not pf.ready(self._infos[i].header_offset + 30):
            return False
        return pf.ready(self._span(i)[1])

    def batches(self):
        """Yield ``(first_index, payloads)`` for runs of consecutive channels.

        Runs follow file order, so each batch can be decoded as soon as the
        prefetch has read it.  A batch is at least ``_BATCH`` channels (or the
        rest of the run) and grows to include every following channel already
        read, so a cached file decodes in one batch.  If file order doesn't
        follow channel order (never seen in ACT data) everything is yielded as
        one batch once it has all been read.
        """
        n = len(self._infos)
        if self._payloads is not None or n == 0:
            yield 0, self.all()
            return
        order = sorted(range(n), key=lambda i: self._infos[i].header_offset)
        runs = []
        for i in order:
            if runs and i == runs[-1][-1] + 1:
                runs[-1].append(i)
            else:
                runs.append([i])
        if len(runs) > max(1, n // 8):
            yield 0, self.all()
            return
        for run in runs:
            k = 0
            while k < len(run):
                end = min(len(run), k + _BATCH)
                while end < len(run) and self._ready(run[end]):
                    end += 1
                yield run[k], [self._view(i) for i in run[k:end]]
                k = end

    def release(self):
        for v in self._views:
            try:
                v.release()
            except Exception:
                pass
        self._views = []


class _open_payloads:
    """Context manager yielding ``(channels, payloads, fields)``.

    ``payloads`` is a :class:`_Payloads`.  The zip is mmapped and, when the
    requested entries cover most of their byte range (a full-TOD read), a
    :class:`_Prefetcher` reads it sequentially in the background.  ``prefetch``
    forces that on (True) or off (False); ``None`` decides automatically.
    """

    def __init__(self, path, channels, fields, prefetch=None):
        self._path, self._channels, self._fields = path, channels, fields
        self._want_prefetch = prefetch
        self._z = self._fh = self._mm = self._pf = self._payloads = None

    def __enter__(self):
        self._z = zipfile.ZipFile(self._path)
        names = set(self._z.namelist())
        fields = self._fields
        if fields is None:
            fields = (parse_format(self._z.read("format").decode("latin-1"))
                      if "format" in names else {})
        channels = self._channels
        if channels is None:
            slm = {n[:-4] for n in names if n.endswith(".slm")}
            channels = sorted(c for c in slm if c.startswith("tesdata"))
        arcnames = [ch + ".slm" if ch + ".slm" in names else ch for ch in channels]
        infos = [self._z.getinfo(a) for a in arcnames]

        if not all(i.compress_type == zipfile.ZIP_STORED for i in infos):
            # Fallback: ordinary (CRC-checked) reads.
            self._payloads = _Payloads(None, infos, None,
                                       payloads=[self._z.read(a) for a in arcnames])
            return channels, self._payloads, fields

        self._fh = open(self._path, "rb")
        self._mm = mmap.mmap(self._fh.fileno(), 0, access=mmap.ACCESS_READ)
        if infos:
            lo = min(i.header_offset for i in infos)
            # The local header can carry an extra field the central directory
            # doesn't; 64 KiB covers the largest one the zip format allows.
            hi = min(len(self._mm),
                     max(i.header_offset + 30 + len(i.filename) + i.compress_size
                         for i in infos) + 65536)
            want = self._want_prefetch
            if want is None:
                dense = sum(i.compress_size for i in infos) >= _PREFETCH_MIN_FRACTION * (hi - lo)
                cached = _cached_fraction(self._fh.fileno(), lo, hi) if dense else None
                want = dense and (cached is None or cached < _CACHED_FRACTION)
            if want:
                self._pf = _Prefetcher(self._path, lo, hi)
        self._payloads = _Payloads(self._mm, infos, self._pf)
        return channels, self._payloads, fields

    def __exit__(self, *exc):
        # Release exported memoryviews before closing the mmap (decoding is done
        # by the time the with-block exits, so any data the caller keeps must be
        # its own copy -- read_zip_array/read_zip both return fresh arrays).
        if self._pf is not None:
            self._pf.close()
        if self._payloads is not None:
            self._payloads.release()
        if self._mm is not None:
            self._mm.close()
        if self._fh is not None:
            self._fh.close()
        if self._z is not None:
            self._z.close()
        return False


def read_zip(path, channels=None, fields=None, workers=None, prefetch=None):
    """Read channels from a zipped slim dirfile into numpy arrays.

    Parameters
    ----------
    path : str
        Path to the ``*.zip`` dirfile.
    channels : list[str] or None
        Channel base names (without the ``.slm`` suffix).  ``None`` reads all
        ``tesdata*`` detector channels.
    fields : dict or None
        Pre-parsed format mapping (see :func:`parse_format`).  Parsed from the
        archive if omitted.
    workers : int or None
        Number of decode threads.  ``decompress`` releases the GIL, so this
        scales across cores.  ``None`` (default) uses ``os.cpu_count()``; pass
        ``1`` to force serial.
    prefetch : bool or None
        Read the zip front to back in a background thread while decoding.
        ``None`` (default) does so when the channels make up most of the file
        (e.g. a full TOD) and it isn't already in the page cache; on a cold
        spinning disk that is ~2x faster than random access.  Small subsets
        are read by random access instead.

    Returns
    -------
    dict[str, numpy.ndarray]
    """
    if workers is None:
        workers = 0  # let the C layer pick hardware_concurrency()

    with _open_payloads(path, channels, fields, prefetch) as (channels, payloads, fields):
        dtypes = [fields.get(ch, (1, np.int32))[1] for ch in channels]

        # Uniform-dtype hot path: decode into one preallocated 2D array.
        if channels and len(set(dtypes)) == 1:
            arr = _decode_uniform(payloads, dtypes[0], int(workers))
            if arr is not None:
                return {ch: arr[i] for i, ch in enumerate(channels)}

        decoded = decompress_many(payloads.all(), int(workers))
        return {ch: np.frombuffer(raw, dtype=dt)
                for ch, raw, dt in zip(channels, decoded, dtypes)}


def _decode_uniform(payloads, dtype, workers):
    """Decode same-dtype payloads into one (n_chan, n_samp) array.

    Decodes batch by batch in file order, so it overlaps with a running
    prefetch.  Returns the 2D array, or None if the first channel's size isn't
    a whole number of samples (caller then uses the general path).
    """
    dtype = np.dtype(dtype)
    n = len(payloads)
    arr = flat = None
    item_nbytes = 0
    for first, batch in payloads.batches():
        if arr is None:
            # Decode one channel to learn the per-channel byte size, then
            # decode everything straight into one preallocated array.
            item_nbytes = len(decompress(batch[0]))
            if item_nbytes % dtype.itemsize:
                return None
            arr = np.empty((n, item_nbytes // dtype.itemsize), dtype=dtype)
            flat = arr.reshape(-1).view(np.uint8)
        decompress_into(batch, flat[first * item_nbytes:(first + len(batch)) * item_nbytes],
                        item_nbytes, workers)
    return arr


def read_zip_array(path, channels=None, fields=None, workers=None, prefetch=None):
    """Like :func:`read_zip` but return ``(channels, 2d_array)``.

    Only valid when all selected channels share a dtype and sample count
    (the normal TOD detector case).  This is the most efficient form: one
    contiguous array, decoded in parallel with no per-channel allocation.
    """
    if workers is None:
        workers = 0
    with _open_payloads(path, channels, fields, prefetch) as (channels, payloads, fields):
        dtypes = {fields.get(ch, (1, np.int32))[1] for ch in channels}
        if len(dtypes) != 1:
            raise ValueError(
                "read_zip_array requires a single dtype across channels")
        arr = _decode_uniform(payloads, dtypes.pop(), int(workers))
    return channels, arr


def read_fields(path, names, workers=None):
    """Read fields by name, resolving RAW channels and LINCOM derived fields.

    ``names`` may include derived fields such as ``"Enc_Az_Deg"`` or
    ``"C_Time"``; their LINCOM definitions in the dirfile ``format`` are applied
    automatically (returned as float64).  RAW fields are returned at native
    dtype.  Missing fields are omitted from the result.

    Returns ``dict[str, numpy.ndarray]``.
    """
    if workers is None:
        workers = 0
    with zipfile.ZipFile(path) as z:
        names_in = set(z.namelist())
        text = z.read("format").decode("latin-1") if "format" in names_in else ""
    raw_fmt = parse_format(text)
    lincom = parse_lincom(text)

    # Figure out which RAW channels we actually need to decode.
    need = set()
    plan = {}  # requested name -> ("raw", ch) or ("lincom", terms)
    for nm in names:
        if nm in raw_fmt:
            plan[nm] = ("raw", nm)
            need.add(nm)
        elif nm in lincom:
            plan[nm] = ("lincom", lincom[nm])
            for f, _, _ in lincom[nm]:
                if f in raw_fmt:
                    need.add(f)
    need = [c for c in need]
    if not need:
        return {}

    raw = read_zip(path, channels=need, fields=raw_fmt, workers=workers)

    out = {}
    for nm, (kind, spec) in plan.items():
        if kind == "raw":
            out[nm] = raw[spec]
        else:
            acc = None
            ok = True
            for f, m, b in spec:
                if f not in raw:
                    ok = False
                    break
                term = raw[f].astype(np.float64) * m + b
                acc = term if acc is None else acc + term
            if ok and acc is not None:
                out[nm] = acc
    return out


def read_dirfile(path, channels=None):
    """Read channels from an *unzipped* slim dirfile directory."""
    import os
    fmt_path = os.path.join(path, "format")
    fields = {}
    if os.path.exists(fmt_path):
        with open(fmt_path, "rb") as f:
            fields = parse_format(f.read().decode("latin-1"))

    if channels is None:
        channels = sorted(
            n[:-4] for n in os.listdir(path)
            if n.endswith(".slm") and n.startswith("tesdata")
        )

    out = {}
    for ch in channels:
        p = os.path.join(path, ch + ".slm")
        if not os.path.exists(p):
            p = os.path.join(path, ch)
        with open(p, "rb") as f:
            raw = decompress(f.read())
        dt = fields.get(ch, (1, np.int32))[1]
        out[ch] = np.frombuffer(raw, dtype=dt)
    return out


# High-level data object (imported last: tod.py depends on the names above).
from .tod import TOD, read_tod  # noqa: E402
