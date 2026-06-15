# actslim — standalone reader for ACT slim-compressed dirfile TODs

Reads ACT (ACTPol/`merlin`) time-ordered data straight from the legacy
**slim**-compressed, zipped *dirfile* archives — **without libactpol, zzip,
getdata, or the autotools build chain**. The only things it needs are a C++
compiler at install time and `numpy` at run time, so it works in any modern
Python (tested on 3.8 and 3.12 / numpy 2.x).

## Why

The legacy path is `libactpol → getdata → zzip → libslim`, built with autotools
and only importable from a pinned old Python with a hand-set `LD_LIBRARY_PATH`.
But for *reading*, the only piece that does real work is the slim codec — the
zip entries are *Stored* (uncompressed), so the zip is just a table of offsets.
`actslim` vendors the proven slim **decode** source and pairs it with the Python
standard library (`zipfile`, `mmap`):

```
zip (Stored)  --mmap zero-copy-->  .slm bytes  --vendored slim decode-->  raw bytes  --np.frombuffer-->  ndarray
```

## Install

```bash
pip install actslim       # prebuilt wheels (Linux x86_64, macOS x86_64/arm64)
```

From source (needs a C++ compiler):

```bash
pip install -e .
```

Wheels for CPython 3.8–3.13 are built with `cibuildwheel` in GitHub Actions and
published to PyPI on a `v*` tag via Trusted Publishing (see
`.github/workflows/wheels.yml`).

## Use

### Detector timeseries as an `(ndet, nsample)` array

```python
import actslim

dets, data = actslim.read_zip_array("1572374891.1572382965.ar6.zip")
# dets: list of 1760 channel names; data: (1760, 259864) int32
```

### A full TOD object (detectors + pointing)

```python
tod = actslim.read_tod("1572374891.1572382965.ar6.zip")

tod.data        # (ndet, nsample) int32 detector array  (alias: tod.signal)
tod.det_uid     # integer detector indices
tod.az          # boresight azimuth,  radians  (n_sample)
tod.alt         # boresight altitude, radians  (n_sample)
tod.ctime       # unix timestamp,     seconds  (n_sample)
tod.enc_flags   # encoder validity flags
tod["tesdatar00c01"]          # one detector row by name
tod = actslim.read_tod("....zip", aux=["enc_status", "data_rate"])  # extra channels
```

`TOD` is a plain dataclass of numpy arrays + metadata (modelled on moby2's TOD),
so it is trivial to convert/wrap for other frameworks. `az`/`alt`/`ctime` are
derived from the dirfile's LINCOM fields (`Enc_Az_Deg`, `Enc_El_Deg`, `C_Time`)
automatically.

### Lower level

```python
# dict {channel: ndarray}, optional subset:
d = actslim.read_zip("....zip", channels=["tesdatar00c01", "az", "el"])

# resolve RAW or derived (LINCOM) fields by name -> {name: ndarray}:
p = actslim.read_fields("....zip", ["Enc_Az_Deg", "Enc_El_Deg", "C_Time"])

# unzipped dirfile directory, or raw bytes:
d = actslim.read_dirfile("/path/to/dirfile")
raw = actslim.decompress(open("tesdatar00c01.slm", "rb").read())
```

`workers=N` controls decode threads (default: all cores; `decompress` releases
the GIL and decoding runs in an internal C thread pool).

## Correctness

Validated **byte-for-byte against libactpol** across all 1760 channels of a real
season-7 TOD (`tests/test_against_libactpol.py`). Channel dtypes come from the
dirfile `format` file using libactpol's getdata type codes
(`S`=int32, `U`=uint32, `s`=int16, `u`=uint16, `f`=float32, `d`=float64, `c`=1 byte).

## Performance

Full TOD (1760 channels × 259864 samples), warm cache, 40-core Xeon:

| reader | end-to-end |
|---|---|
| legacy `DirfileManager.load_channels` (libactpol, 40 OMP threads) | ~0.84 s |
| `actslim.read_zip_array` (40 threads) | **~0.54–0.64 s** |

Comparable-to-better than the original, with none of the legacy dependencies.
The single biggest real-world lever for the legacy reader — making sure
`OMP_NUM_THREADS` is high — applies here too via `workers`.

### Why it's faster than well-optimized legacy code

We did **not** make slim decoding faster — both readers run the same codec at
the same speed (~0.78 s for 1760 channels at 40 threads, measured head to head).
The wins are entirely *around* the decode, by removing work the layered
`libactpol → getdata → zzip` stack can't avoid without a rewrite:

1. **No data copy in.** The legacy path reads each `.slm` through `zzip`, which
   copies bytes out of the zip before slim sees them. We `mmap` the zip and let
   slim decode *straight from the mapped file* (`fmemopen` over a memoryview
   slice) — ~600 MB that used to be copied is never copied. This is only
   possible because the zip is *Stored* (no deflate), so each `.slm` is a
   contiguous slice of the file.
2. **No CRC pass.** Python's `zipfile.read()` computes a CRC-32 over every entry
   (~0.9 s serial on 600 MB). The mmap path skips it. (An earlier version that
   used `zipfile.read()` was slower for exactly this reason.)
3. **No per-channel allocation.** `decompress_into` decodes every channel into
   one preallocated 2-D array. An earlier version that `malloc`'d per channel
   plateaued at ~10 threads from glibc malloc-arena contention; decoding into
   one array scales to all cores.

Takeaway: "optimized" usually means the inner loop is fast. A layered design can
still be slower than the sum of its tuned parts, because each layer boundary
forces a copy.

**Caveats.** These numbers are **warm cache** (data already in the page cache);
cold from disk both readers are I/O-bound and roughly equal. mmap zero-copy
shines on **local** filesystems — on a network FS with poor random access,
staging the zip to `/dev/shm` first (the legacy `MOBY2_TOD_STAGING_PATH` trick)
still helps and applies to both.

## Layout

- `actslim/_actslim.cpp` — CPython extension: `decompress`, `decompress_many`,
  `decompress_into` (decode straight into a preallocated array).
- `actslim/_vendor_slim/` — vendored slim **decode** source (GPLv3, J. Fowler);
  the only local change is an added `slim_expander_t(FILE*)` constructor so we
  can decode from an in-memory buffer via `fmemopen`.
- `actslim/__init__.py` — zip/mmap container handling and `format` parsing.
