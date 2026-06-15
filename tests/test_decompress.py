"""Self-contained decode tests (no ACT data, no libactpol).

Validates the vendored slim decoder against a committed fixture: a real .slm
payload and its expected raw bytes.
"""
import json
import os

import numpy as np

import actslim

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")


def _load_fixture():
    meta = json.load(open(os.path.join(FIX, "sample.json")))
    slm = open(os.path.join(FIX, "sample.slm"), "rb").read()
    raw = open(os.path.join(FIX, "sample.raw"), "rb").read()
    return meta, slm, raw


def test_decompress_matches_fixture():
    _, slm, raw = _load_fixture()
    assert actslim.decompress(slm) == raw


def test_decompress_many_matches_single():
    _, slm, raw = _load_fixture()
    out = actslim.decompress_many([slm] * 5, 4)
    assert out == [raw] * 5


def test_decompress_into_preallocated():
    _, slm, raw = _load_fixture()
    n, item = 5, len(raw)
    buf = np.empty(n * item, dtype=np.uint8)
    actslim.decompress_into([slm] * n, buf, item, 4)
    for i in range(n):
        assert buf[i * item:(i + 1) * item].tobytes() == raw


def test_dtype_roundtrip():
    meta, slm, raw = _load_fixture()
    arr = np.frombuffer(actslim.decompress(slm), dtype=np.dtype(meta["dtype"]))
    assert arr.nbytes == len(raw)
