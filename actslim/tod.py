"""A small, dependency-light data object for ACT slim TODs.

A generic container -- numpy arrays plus metadata -- modelled on moby2's TOD
(``data``/``det_uid``/``az``/``alt``/``ctime``) but without dragging in moby2
or any legacy types, so it is easy to hand to modern code.
"""

import os
import zipfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from . import _decode_uniform, parse_format, read_fields, _open_payloads

__all__ = ["TOD", "read_tod"]

# Derived pointing fields, mapped to TOD attributes.  Names follow the standard
# ACT dirfile LINCOM definitions; angles are converted from degrees to radians
# to match moby2 (boresight az/alt in radians, ctime in unix seconds).
_POINTING = [
    ("az", "Enc_Az_Deg", np.deg2rad),
    ("alt", "Enc_El_Deg", np.deg2rad),
    ("ctime", "C_Time", None),
    ("enc_flags", "enc_flags", None),
]


@dataclass
class TOD:
    """Time-ordered data from one ACT slim dirfile.

    Attributes
    ----------
    name : str
        Source basename (e.g. ``"1572374891.1572349526.ar5"``).
    dets : list[str]
        Detector channel names, in row order of ``data``.
    data : numpy.ndarray
        ``(n_dets, n_samp)`` detector array (int32 raw MCE counts).
    det_uid : numpy.ndarray
        Integer detector indices ``arange(n_dets)`` (row index into ``data``).
    nsamps : int
        Number of samples.
    az, alt : numpy.ndarray or None
        Boresight azimuth / altitude in **radians** (n_samp), if available.
    ctime : numpy.ndarray or None
        Unix timestamp per sample (n_samp), if available.
    enc_flags : numpy.ndarray or None
        Encoder/pointing validity flags (n_samp), if available.
    aux : dict[str, numpy.ndarray]
        Any other channels that were requested, at their native sample rate.
    spf : dict[str, int]
        Samples-per-frame per field (key ``"dets"`` for the detector block).
    fmt : dict[str, tuple[int, numpy.dtype]]
        The full parsed dirfile ``format`` (every RAW field), for reference.
    """

    name: str
    dets: List[str]
    data: np.ndarray
    det_uid: np.ndarray = None
    nsamps: int = 0
    az: Optional[np.ndarray] = None
    alt: Optional[np.ndarray] = None
    ctime: Optional[np.ndarray] = None
    enc_flags: Optional[np.ndarray] = None
    aux: Dict[str, np.ndarray] = field(default_factory=dict)
    spf: Dict[str, int] = field(default_factory=dict)
    fmt: Dict[str, tuple] = field(default_factory=dict)

    def __post_init__(self):
        if self.det_uid is None:
            self.det_uid = np.arange(len(self.dets))
        if not self.nsamps:
            self.nsamps = 0 if self.data is None else self.data.shape[1]

    # -- moby2-compatible alias ------------------------------------------
    @property
    def signal(self) -> np.ndarray:
        return self.data

    # -- sizes -----------------------------------------------------------
    @property
    def ndets(self) -> int:
        return len(self.dets)

    @property
    def nframes(self) -> int:
        s = self.spf.get("dets", 1)
        return self.nsamps // s if s else 0

    def __len__(self) -> int:
        return self.ndets

    # -- access ----------------------------------------------------------
    def det_index(self, name: str) -> int:
        return self.dets.index(name)

    def __getitem__(self, key):
        """``tod[i]`` / ``tod["tesdatar00c01"]`` -> a detector row;
        ``tod["enc_status"]`` -> an auxiliary channel."""
        if isinstance(key, str):
            if key in self.aux:
                return self.aux[key]
            return self.data[self.det_index(key)]
        return self.data[key]

    def get(self, name: str, default=None):
        if name in self.aux:
            return self.aux[name]
        if name in self.dets:
            return self.data[self.det_index(name)]
        return default

    def __contains__(self, name: str) -> bool:
        return name in self.aux or name in self.dets

    def __repr__(self) -> str:
        has = ",".join(k for k in ("az", "alt", "ctime")
                       if getattr(self, k) is not None)
        return ("TOD(name=%r, ndets=%d, nsamps=%d, pointing=[%s], aux=%d)"
                % (self.name, self.ndets, self.nsamps, has, len(self.aux)))

    # -- export ----------------------------------------------------------
    def to_hdf5(self, path, compression="gzip", compression_opts=4):
        """Write this TOD to an HDF5 file (requires ``h5py``)."""
        import h5py

        with h5py.File(path, "w") as h:
            h.attrs["name"] = self.name
            h.attrs["nsamps"] = self.nsamps
            h.attrs["spf_dets"] = self.spf.get("dets", 0)
            h.create_dataset("dets", data=np.array(self.dets, dtype="S"))
            h.create_dataset("det_uid", data=self.det_uid)
            h.create_dataset("data", data=self.data, compression=compression,
                             compression_opts=compression_opts)
            for nm in ("az", "alt", "ctime", "enc_flags"):
                v = getattr(self, nm)
                if v is not None:
                    h.create_dataset(nm, data=v, compression=compression,
                                     compression_opts=compression_opts)
            g = h.create_group("aux")
            for k, v in self.aux.items():
                d = g.create_dataset(k, data=v, compression=compression,
                                     compression_opts=compression_opts)
                d.attrs["spf"] = self.spf.get(k, 0)
        return path


def read_tod(path, dets=None, aux=(), pointing=True, workers=None):
    """Read an ACT slim zipped dirfile into a :class:`TOD`.

    Parameters
    ----------
    path : str
        Path to the ``*.zip`` dirfile.
    dets : list[str] or None
        Detector channels for the 2D ``data`` block.  ``None`` loads all
        ``tesdata*`` channels (they must share a dtype and sample count).
    aux : iterable[str]
        Extra channel names to load into ``TOD.aux`` at native rate.
    pointing : bool
        If True (default), populate ``az``/``alt`` (radians), ``ctime`` (unix s)
        and ``enc_flags`` from the dirfile's derived (LINCOM) fields when present.
    workers : int or None
        Decode threads (default: all cores).

    Returns
    -------
    TOD
    """
    aux = list(aux)
    name = os.path.basename(path)
    if name.endswith(".zip"):
        name = name[:-4]
    if workers is None:
        workers = 0

    with zipfile.ZipFile(path) as z:
        names = set(z.namelist())
        fmt = (parse_format(z.read("format").decode("latin-1"))
               if "format" in names else {})
        if dets is None:
            slm = {n[:-4] for n in names if n.endswith(".slm")}
            dets = sorted(c for c in slm if c.startswith("tesdata"))

    # Detector block: one preallocated (n_det, n_samp) array.
    with _open_payloads(path, list(dets), fmt) as (dets, payloads, fmt):
        det_dtypes = {fmt.get(c, (1, np.int32))[1] for c in dets}
        if len(det_dtypes) != 1:
            raise ValueError("detector channels must share a single dtype")
        data = _decode_uniform(payloads, det_dtypes.pop(), int(workers))
        spf = {"dets": fmt.get(dets[0], (1, None))[0] if dets else 1}

    # Auxiliary channels (mixed dtype / rate).
    aux_data = {}
    if aux:
        got = read_fields(path, aux, workers=workers)
        for cname in aux:
            if cname in got:
                aux_data[cname] = got[cname]
                spf[cname] = fmt.get(cname, (1, None))[0]

    tod = TOD(name=name, dets=list(dets), data=data,
              nsamps=data.shape[1] if data is not None else 0,
              aux=aux_data, spf=spf, fmt=fmt)

    if pointing:
        want = [src for _, src, _ in _POINTING]
        got = read_fields(path, want, workers=workers)
        for attr, src, conv in _POINTING:
            if src in got:
                v = got[src]
                setattr(tod, attr, conv(v) if conv is not None else v)

    return tod
