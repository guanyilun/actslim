"""Tests for the sequential-prefetch read path (no ACT data needed).

Builds small zipped dirfiles from committed real ``.slm`` payloads and checks
that every read path -- prefetch on, off, automatic; file order matching
channel order or not; Stored or deflated -- puts each channel in the right row.
"""
import glob
import os
import random
import zipfile

import numpy as np
import pytest

import actslim

HERE = os.path.dirname(os.path.abspath(__file__))
ROWS = sorted(glob.glob(os.path.join(HERE, "fixtures", "rows", "*.slm")))


def _payloads():
    return [open(p, "rb").read() for p in ROWS]


def _make_zip(path, n_chan, order="sorted", compression=zipfile.ZIP_STORED):
    """Write ``n_chan`` tesdata channels cycling through the fixture payloads.

    Returns ``(channels, expected)`` with channels sorted by name and expected
    the ``(n_chan, n_samp)`` uint32 array they should decode to.
    """
    slms = _payloads()
    raws = [np.frombuffer(actslim.decompress(s), dtype=np.uint32) for s in slms]
    chans = ["tesdatar%02dc%02d" % divmod(i, 64) for i in range(n_chan)]
    idx = list(range(n_chan))
    if order == "reversed":
        idx.reverse()
    elif order == "shuffled":
        random.Random(0).shuffle(idx)
    fmt = "".join("%s RAW U 1\n" % c for c in chans)
    with zipfile.ZipFile(path, "w", compression) as z:
        z.writestr("format", fmt)
        for i in idx:
            z.writestr(chans[i] + ".slm", slms[i % len(slms)])
    expected = np.stack([raws[i % len(raws)] for i in range(n_chan)])
    return chans, expected


def test_fixture_payloads_are_distinct():
    raws = {actslim.decompress(s) for s in _payloads()}
    assert len(ROWS) >= 4 and len(raws) == len(ROWS)


@pytest.mark.parametrize("order", ["sorted", "reversed", "shuffled"])
@pytest.mark.parametrize("prefetch", [None, True, False])
def test_read_zip_array_rows(tmp_path, order, prefetch):
    # 300 channels > _BATCH, so the sorted case decodes in several batches.
    path = str(tmp_path / "tod.zip")
    chans, expected = _make_zip(path, 300, order)
    got_chans, arr = actslim.read_zip_array(path, prefetch=prefetch, workers=4)
    assert got_chans == chans
    assert arr.dtype == np.uint32
    np.testing.assert_array_equal(arr, expected)


@pytest.mark.parametrize("prefetch", [None, True, False])
def test_read_zip_dict(tmp_path, prefetch):
    path = str(tmp_path / "tod.zip")
    chans, expected = _make_zip(path, 200)
    got = actslim.read_zip(path, prefetch=prefetch, workers=4)
    assert sorted(got) == chans
    for i, c in enumerate(chans):
        np.testing.assert_array_equal(got[c], expected[i])


def _count_prefetches(monkeypatch):
    started = []
    real = actslim._Prefetcher
    monkeypatch.setattr(actslim, "_Prefetcher",
                        lambda *a: started.append(a) or real(*a))
    return started


def test_subset_uses_random_access(tmp_path, monkeypatch):
    path = str(tmp_path / "tod.zip")
    chans, expected = _make_zip(path, 200)
    monkeypatch.setattr(actslim, "_cached_fraction", lambda *a: 0.0)  # "cold"
    started = _count_prefetches(monkeypatch)
    pick = [chans[3], chans[150]]
    got = actslim.read_zip(path, channels=pick)
    assert not started  # 2 of 200 channels: too sparse to prefetch
    np.testing.assert_array_equal(got[pick[0]], expected[3])
    np.testing.assert_array_equal(got[pick[1]], expected[150])
    actslim.read_zip_array(path)
    assert len(started) == 1  # full read: prefetched


def test_cached_file_skips_prefetch(tmp_path, monkeypatch):
    path = str(tmp_path / "tod.zip")
    chans, expected = _make_zip(path, 200)
    started = _count_prefetches(monkeypatch)
    actslim.read_zip_array(path)  # just written, so it is in the page cache
    fd = os.open(path, os.O_RDONLY)
    try:
        mincore_works = actslim._cached_fraction(fd, 0, 1) is not None
    finally:
        os.close(fd)
    if mincore_works:
        assert not started
    monkeypatch.setattr(actslim, "_cached_fraction", lambda *a: None)  # unknown
    _, arr = actslim.read_zip_array(path)
    assert len(started) == 1
    np.testing.assert_array_equal(arr, expected)


def test_deflated_zip_fallback(tmp_path):
    path = str(tmp_path / "tod.zip")
    chans, expected = _make_zip(path, 50, compression=zipfile.ZIP_DEFLATED)
    got_chans, arr = actslim.read_zip_array(path)
    assert got_chans == chans
    np.testing.assert_array_equal(arr, expected)


def test_read_tod_prefetch(tmp_path):
    path = str(tmp_path / "tod.zip")
    chans, expected = _make_zip(path, 130)
    for prefetch in (None, True, False):
        tod = actslim.read_tod(path, prefetch=prefetch, pointing=False)
        assert tod.dets == chans
        np.testing.assert_array_equal(tod.data, expected)


def test_prefetcher_stops_at_eof(tmp_path):
    # A stop past the end of the file must not leave wait() hanging.
    path = tmp_path / "f.bin"
    path.write_bytes(b"x" * 1000)
    pf = actslim._Prefetcher(str(path), 0, 10 ** 9)
    pf.wait(10 ** 9)
    pf.close()
