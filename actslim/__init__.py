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

import mmap
import struct
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


def _mmap_stored_payloads(z, mm, arcnames):
    """Return zero-copy memoryview slices for *Stored* zip entries.

    Reads each entry's local file header to find where its data begins, then
    slices the mmap directly -- no decompression, no CRC, no copy.  Falls back
    to ``None`` for any entry that is not Stored (caller uses ``z.read``).
    """
    views = []
    for name in arcnames:
        info = z.getinfo(name)
        if info.compress_type != zipfile.ZIP_STORED:
            return None
        off = info.header_offset
        sig, ver, flag, method, t, d, crc, csz, usz, nlen, elen = \
            _LFH.unpack(mm[off:off + 30])
        if sig != 0x04034b50:
            return None
        start = off + 30 + nlen + elen
        views.append(memoryview(mm)[start:start + info.compress_size])
    return views


def read_zip(path, channels=None, fields=None, workers=None):
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

    Returns
    -------
    dict[str, numpy.ndarray]
    """
    if workers is None:
        workers = 0  # let the C layer pick hardware_concurrency()

    with _open_payloads(path, channels, fields) as (channels, payloads, fields):
        dtypes = [fields.get(ch, (1, np.int32))[1] for ch in channels]

        # Uniform-dtype hot path: decode into one preallocated 2D array.
        if channels and len(set(dtypes)) == 1:
            arr = _decode_uniform(payloads, dtypes[0], int(workers))
            if arr is not None:
                return {ch: arr[i] for i, ch in enumerate(channels)}

        decoded = decompress_many(payloads, int(workers))
        return {ch: np.frombuffer(raw, dtype=dt)
                for ch, raw, dt in zip(channels, decoded, dtypes)}


class _open_payloads:
    """Context manager yielding ``(channels, payloads, fields)``.

    ``payloads`` are zero-copy mmap memoryviews for Stored entries (the ACT
    case), or plain ``bytes`` from ``z.read`` otherwise.  The mmap is held open
    for the duration of the ``with`` block.
    """

    def __init__(self, path, channels, fields):
        self._path, self._channels, self._fields = path, channels, fields
        self._z = self._fh = self._mm = None
        self._views = None

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

        payloads = None
        try:
            self._fh = open(self._path, "rb")
            self._mm = mmap.mmap(self._fh.fileno(), 0, access=mmap.ACCESS_READ)
            payloads = _mmap_stored_payloads(self._z, self._mm, arcnames)
        except Exception:
            payloads = None
        if payloads is None:
            # Fallback: ordinary (CRC-checked) reads.
            payloads = [self._z.read(a) for a in arcnames]
        else:
            self._views = payloads  # keep refs so we can release before close
        return channels, payloads, fields

    def __exit__(self, *exc):
        # Release exported memoryviews before closing the mmap (decoding is done
        # by the time the with-block exits, so any data the caller keeps must be
        # its own copy -- read_zip_array/read_zip both return fresh arrays).
        if self._views is not None:
            for v in self._views:
                try:
                    v.release()
                except Exception:
                    pass
            self._views = None
        if self._mm is not None:
            self._mm.close()
        if self._fh is not None:
            self._fh.close()
        if self._z is not None:
            self._z.close()
        return False


def _decode_uniform(payloads, dtype, workers):
    """Decode same-dtype payloads into one (n_chan, n_samp) array.

    Returns the 2D array, or None if the channels don't share a sample count
    (caller then uses the general path).
    """
    dtype = np.dtype(dtype)
    # Decode the first channel to learn the per-channel byte size, then decode
    # the rest straight into one preallocated array.
    item_nbytes = len(decompress(payloads[0]))
    if item_nbytes % dtype.itemsize:
        return None
    n_samp = item_nbytes // dtype.itemsize
    n = len(payloads)
    arr = np.empty((n, n_samp), dtype=dtype)
    decompress_into(payloads, arr.reshape(-1).view(np.uint8), item_nbytes, workers)
    return arr


def read_zip_array(path, channels=None, fields=None, workers=None):
    """Like :func:`read_zip` but return ``(channels, 2d_array)``.

    Only valid when all selected channels share a dtype and sample count
    (the normal TOD detector case).  This is the most efficient form: one
    contiguous array, decoded in parallel with no per-channel allocation.
    """
    if workers is None:
        workers = 0
    with _open_payloads(path, channels, fields) as (channels, payloads, fields):
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
