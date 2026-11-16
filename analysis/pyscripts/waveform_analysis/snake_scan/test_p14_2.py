#!/usr/bin/env python3
"""P14.2 slice-2 tests (pytest or `python test_p14_2.py`).  Read out/edge_hist/edge_<module>.npz and out/edge_23701.tsv
(made by `snake_scan.py edge --dwell out/dwell_23701.tsv <9 files>`).  Criteria are the plan's, unmodified."""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import edge_pass as EP  # noqa: E402

TSV = HERE / "out" / "edge_23701.tsv"
HIST = HERE / "out" / "edge_hist"
W = EP.load_calib(EP.resolve_calib(23701))


def _hist(m):
    d = np.load(HIST / f"edge_{m}.npz")
    return {k: d[k] for k in d.files}


def _table():
    rows = {}
    for ln in TSV.read_text().splitlines():
        if ln.startswith("#") or ln.startswith("module\t"):
            continue
        v = ln.split("\t")
        rows[v[0]] = dict(zip(EP.COLS, [v[0]] + [float(x) for x in v[1:]]))
    return rows


def test_edge_w547():
    E, sig, *_ = EP.fit_edge(_hist("W547"))
    exp = 18500.0 * W["W547"]
    print(f"E_edge={E:.1f} MeV ({E / W['W547']:.0f} ADC) expected {exp:.1f}; sigma={sig:.1f} ({100 * sig / E:.2f} %)")
    assert abs(E / exp - 1) < 0.03, (E, exp)
    assert sig < 0.05 * E


def test_edge_pileup():
    h = _hist("W547")
    f = float(h["n_pileup"]) / (float(h["n_single"]) + float(h["n_pileup"]))
    print(f"n_pileup fraction {100 * f:.1f} %")
    assert 0.02 <= f <= 0.20, f


def test_edge_windows():
    sp = EP.fit_edge(_hist("W547"))[5]
    print(f"window spread {100 * sp:.3f} %")
    assert sp < 0.005, sp


def test_edge_consistency():
    t = _table()
    s = np.array([t[m]["s"] for m in ("W547", "W552", "W580", "W599")])
    print("s:", dict(zip(("W547", "W552", "W580", "W599"), np.round(s, 4))))
    assert (s.max() - s.min()) / s.mean() < 0.01, s


if __name__ == "__main__":
    bad = 0
    for k, f in sorted(globals().items()):
        if k.startswith("test_"):
            try:
                f()
                print("PASS", k)
            except AssertionError as e:
                bad += 1
                print("FAIL", k, e)
    sys.exit(bad)
