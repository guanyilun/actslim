"""Byte-for-byte validation of actslim against the libactpol oracle.

This needs the legacy moby2/libactpol stack to produce reference data, so run
it under that environment, e.g.:

    OMP_NUM_THREADS=40 \
    PYTHONPATH=/home/yguan/shared/software/moby2/lib/python3.8/site-packages \
    LD_LIBRARY_PATH=/home/yguan/shared/software/moby2/lib \
    python tests/test_against_libactpol.py [TOD.zip]

actslim itself needs none of that; the oracle does.
"""
import sys
import zipfile

import numpy as np

sys.path.insert(0, "/home/yguan/software/actslim")
import actslim

DEFAULT_TOD = ("/home/data/act/tod/actpol/tod/season7/merlin/15723/"
               "1572374891.1572382965.ar6.zip")


def main(path=DEFAULT_TOD):
    tes = sorted(n[:-4] for n in zipfile.ZipFile(path).namelist()
                 if n.startswith("tesdatar") and n.endswith(".slm"))

    from moby2.util.dirfile import DirfileManager
    oracle = np.asarray(DirfileManager(path).load_channels(tes, dtype="int32"))

    _, got = actslim.read_zip_array(path, channels=tes)

    assert got.shape == oracle.shape, (got.shape, oracle.shape)
    bad = int((got != oracle).any(axis=1).sum())
    print("channels: %d  mismatches: %d" % (len(tes), bad))
    assert bad == 0, "actslim output differs from libactpol"
    print("OK: byte-for-byte identical to libactpol")


if __name__ == "__main__":
    main(*sys.argv[1:])
