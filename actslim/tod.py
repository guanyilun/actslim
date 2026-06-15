"""A small, dependency-light data object for ACT slim TODs.

The goal is a *generic* container -- just numpy arrays plus metadata -- that
does not drag in moby2 or any legacy types, so it is easy to hand to modern
code (numpy/scipy/xarray/HDF5/...).
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from . import _decode_uniform, decompress_many, parse_format, _open_payloads

__all__ = ["TOD", "read_tod"]


@dataclass
class TOD:
    """Time-ordered data from one ACT slim dirfile.

    Attributes
    ----------
    name : str
        Source basename (e.g. ``"1572374891.1572382965.ar6"``).
    dets : list[str]
        Detector channel names, in row order of ``signal``.
    signal : numpy.ndarray
        ``(n_dets, n_samp)`` detector array (int32 for raw MCE counts).
    aux : dict[str, numpy.ndarray]
        Any other channels that were read (housekeeping, encoders, flags...),
        each at its own native sample rate -- *not* resampled.
    spf : dict[str, int]
        Samples-per-frame for every field, including the key ``"dets"``.
        Channels with different ``spf`` are sampled at different rates.
    fmt : dict[str, tuple[int, numpy.dtype]]
        The full parsed dirfile ``format`` (every RAW field), for reference.
    """

    name: str
    dets: List[str]
    signal: np.ndarray
    aux: Dict[str, np.ndarray] = field(default_factory=dict)
    spf: Dict[str, int] = field(default_factory=dict)
    fmt: Dict[str, tuple] = field(default_factory=dict)

    # -- sizes -----------------------------------------------------------
    @property
    def ndets(self) -> int:
        return len(self.dets)

    @property
    def nsamp(self) -> int:
        return 0 if self.signal is None else self.signal.shape[1]

    @property
    def nframes(self) -> int:
        s = self.spf.get("dets", 1)
        return self.nsamp // s if s else 0

    def __len__(self) -> int:
        return self.ndets

    # -- access ----------------------------------------------------------
    def det_index(self, name: str) -> int:
        return self.dets.index(name)

    def __getitem__(self, key):
        """``tod[i]`` or ``tod["tesdatar00c01"]`` -> a detector row;
        ``tod["enc_flags"]`` -> an auxiliary channel."""
        if isinstance(key, str):
            if key in self.aux:
                return self.aux[key]
            return self.signal[self.det_index(key)]
        return self.signal[key]

    def get(self, name: str, default=None):
        if name in self.aux:
            return self.aux[name]
        if name in self.dets:
            return self.signal[self.det_index(name)]
        return default

    def __contains__(self, name: str) -> bool:
        return name in self.aux or name in self.dets

    def __repr__(self) -> str:
        return ("TOD(name=%r, ndets=%d, nsamp=%d, nframes=%d, aux=%d fields)"
                % (self.name, self.ndets, self.nsamp, self.nframes, len(self.aux)))

    # -- export (Option B: convert to a modern format) -------------------
    def to_hdf5(self, path, compression="gzip", compression_opts=4):
        """Write this TOD to an HDF5 file (requires ``h5py``).

        Detectors go to ``/signal`` (with ``/dets`` names); each aux field to
        ``/aux/<name>``; metadata to attributes.  This makes a clean one-time
        conversion away from the legacy slim format.
        """
        import h5py

        with h5py.File(path, "w") as h:
            h.attrs["name"] = self.name
            h.attrs["ndets"] = self.ndets
            h.attrs["nsamp"] = self.nsamp
            h.attrs["spf_dets"] = self.spf.get("dets", 0)
            h.create_dataset("dets",
                             data=np.array(self.dets, dtype="S"))
            h.create_dataset("signal", data=self.signal,
                             compression=compression,
                             compression_opts=compression_opts)
            g = h.create_group("aux")
            for k, v in self.aux.items():
                d = g.create_dataset(k, data=v, compression=compression,
                                     compression_opts=compression_opts)
                d.attrs["spf"] = self.spf.get(k, 0)
        return path


def read_tod(path, dets=None, aux=(), workers=None):
    """Read an ACT slim zipped dirfile into a :class:`TOD`.

    Parameters
    ----------
    path : str
        Path to the ``*.zip`` dirfile.
    dets : list[str] or None
        Detector channels to load as the 2D ``signal`` block.  ``None`` loads
        all ``tesdata*`` channels.  All selected detectors must share a dtype
        and sample count (the normal case).
    aux : iterable[str]
        Additional channel names to load into ``TOD.aux`` at native rate
        (e.g. ``["enc_flags", "data_rate"]``).
    workers : int or None
        Decode threads (default: all cores).

    Returns
    -------
    TOD
    """
    import os
    import zipfile

    aux = list(aux)
    name = os.path.basename(path)
    if name.endswith(".zip"):
        name = name[:-4]

    if workers is None:
        workers = 0

    # Resolve the detector list and format up front (cheap: central directory).
    with zipfile.ZipFile(path) as z:
        names = set(z.namelist())
        fmt = (parse_format(z.read("format").decode("latin-1"))
               if "format" in names else {})
        if dets is None:
            slm = {n[:-4] for n in names if n.endswith(".slm")}
            dets = sorted(c for c in slm if c.startswith("tesdata"))

    allch = list(dets) + aux
    n_det = len(dets)

    with _open_payloads(path, allch, fmt) as (allch, payloads, fmt):
        det_payloads = payloads[:n_det]
        aux_payloads = payloads[n_det:]

        det_dtypes = {fmt.get(c, (1, np.int32))[1] for c in dets}
        if len(det_dtypes) != 1:
            raise ValueError("detector channels must share a single dtype")
        signal = _decode_uniform(det_payloads, det_dtypes.pop(), int(workers))

        aux_data = {}
        spf = {"dets": fmt.get(dets[0], (1, None))[0] if dets else 1}
        if aux_payloads:
            decoded = decompress_many(aux_payloads, int(workers))
            for cname, raw in zip(aux, decoded):
                spf_c, dt = fmt.get(cname, (1, np.int32))
                aux_data[cname] = np.frombuffer(raw, dtype=dt)
                spf[cname] = spf_c

    return TOD(name=name, dets=list(dets), signal=signal,
               aux=aux_data, spf=spf, fmt=fmt)
