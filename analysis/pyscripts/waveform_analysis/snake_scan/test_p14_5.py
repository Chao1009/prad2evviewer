#!/usr/bin/env python3
"""P14.5 slice-5 tests (pytest or `python test_p14_5.py`).  tracks_rate and the peak-sample half of mip_mpv decode the
first 30 k events of 23744 file 00000 (about 20 s); the integral MPV, file consistency and group counts read
out/mip_23744.tsv (made by `python cosmic_mip.py <3 files> --out out/mip_23744.tsv`)."""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cosmic_mip as CM  # noqa: E402

EVIO = Path("/Users/L00338853/work/PRad/data/evio/prad_023744.evio.00000")
TSV = HERE / "out" / "mip_23744.tsv"
_R = {}


def _scan30k():
    if "r" not in _R:
        _R["r"] = CM.scan_file(EVIO, CM.SSN.DERIVED_MAP, max_events=30000)
    return _R["r"]


def _table():
    lines = [ln for ln in TSV.read_text().splitlines() if not ln.startswith("#")]
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"))) for ln in lines[1:]]


def test_tracks_rate():
    n = _scan30k()["counts"]["n_tracks"]
    assert 150 <= n <= 300, n


def test_mip_mpv():
    h = _scan30k()["hits"]
    pk = float(np.median(h["peak_adc"]))
    assert 25 <= pk <= 40, pk
    rows = {r["group"]: r for r in _table()}
    a = rows["all"]
    mpv = float(a["mpv"])
    assert 150 <= mpv <= 260, mpv
    files = [k[4:] for k in a if k.startswith("mpv_")]
    assert len(files) >= 2, files
    for i, f in enumerate(files):
        for g in files[i + 1:]:
            d = abs(float(a["mpv_" + f]) - float(a["mpv_" + g]))
            se = np.hypot(float(a["err_" + f]), float(a["err_" + g]))
            assert d <= 2 * se, (f, g, d, se)


def test_mip_groups():
    rows = {r["group"]: r for r in _table()}
    for g in CM.GROUPS:
        assert int(rows[g]["n_modules"]) >= 20 and int(rows[g]["n_hits"]) >= 300, (g, rows[g]["n_modules"], rows[g]["n_hits"])


if __name__ == "__main__":
    ok = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print("PASS", name)
            except AssertionError as e:
                ok = 1
                print("FAIL", name, e)
    sys.exit(ok)
