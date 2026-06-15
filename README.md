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

```python
import actslim

# All detector channels as one (n_chan, n_samp) array (fastest):
channels, data = actslim.read_zip_array("1572374891.1572382965.ar6.zip")

# Or a dict {channel: ndarray}, optionally a subset:
tod = actslim.read_zip("....zip", channels=["tesdatar00c01", "tesdatar00c02"])

# Unzipped dirfile directory:
tod = actslim.read_dirfile("/path/to/dirfile")

# Low level:
raw = actslim.decompress(open("tesdatar00c01.slm","rb").read())   # bytes
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

## Layout

- `actslim/_actslim.cpp` — CPython extension: `decompress`, `decompress_many`,
  `decompress_into` (decode straight into a preallocated array).
- `actslim/_vendor_slim/` — vendored slim **decode** source (GPLv3, J. Fowler);
  the only local change is an added `slim_expander_t(FILE*)` constructor so we
  can decode from an in-memory buffer via `fmemopen`.
- `actslim/__init__.py` — zip/mmap container handling and `format` parsing.
