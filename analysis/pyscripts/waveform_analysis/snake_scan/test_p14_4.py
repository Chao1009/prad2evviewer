#!/usr/bin/env python3
"""P14.4 slice-4 tests (pytest or `python test_p14_4.py`).  transfer_units is synthetic; the other two read the tables
in out/ and the V3 ring/gains.tsv (control_phase2_v3_20261114/v3_defaults); nothing is written to out/.
transfer_spread compares the verdict with ER/improvement/P14_snake_scan.tsv (primary L_scan = 23701 alone, pool group)."""
import math
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import transfer as T  # noqa: E402

RTS = ("2p24GeV_H2", "3p49GeV_H2")
STUDY = HERE.parent / "energy_resolution" / "improvement" / "P14_snake_scan.tsv"


def _synthetic(s=1.0, Ls=5000.0, Lp=5000.0, g=0.97, n=4):
    mods = [f"W{i}" for i in range(n)]
    edge = {m: (s, 0.005) for m in mods}
    ls = {m: (Ls, 10.0, "23701") for m in mods}
    lp = {m: (Lp, 5.0, "x") for m in mods}
    gp = {m: (g, 0.002, 1) for m in mods}
    return T.transfer(edge, ls, lp, gp, {m: "pool" for m in mods})


def test_transfer_units():
    """s = 1, L_scan = L_prod: g_T = N (and N = g_pool)."""
    res = _synthetic()
    assert abs(res["N"] - 0.97) < 1e-12
    for r in res["rows"]:
        assert abs(r["g_T"] - res["N"]) < 1e-12 and abs(r["r"] - 1.0) < 1e-12
    # a module with s = 1.05 and the same L: g_T = N / 1.05; an LMS ratio of 1.1 multiplies g_T by 1.1
    edge = {"A": (1.0, 0.0), "B": (1.05, 0.0)}
    ls = {"A": (1000.0, 1.0, ""), "B": (1100.0, 1.0, "")}
    lp = {"A": (1000.0, 1.0, ""), "B": (1000.0, 1.0, "")}
    out = T.transfer(edge, ls, lp, {"A": (0.97, 0.001, 1)}, {"A": "pool", "B": "O"})
    b = [r for r in out["rows"] if r["module"] == "B"][0]
    assert abs(b["g_T"] - 0.97 * 1.1 / 1.05) < 1e-12 and math.isnan(b["r"])


def test_transfer_pool_median():
    """median r over the pool modules = 1 by construction (both run types, both L_scan variants)."""
    reg = T.region_of()
    edge = T.read_edge(T.OUT / "edge_23701.tsv")
    for rt in RTS:
        gp = T.read_pool_gains(T.V3 / rt / "ring" / "gains.tsv")
        lp = T.read_L(T.OUT / f"lms_prod_{rt}.tsv")
        for ls in (T.read_L(T.OUT / "lms_scan.tsv"), T.read_L(T.OUT / "lms_scan_merged.tsv")):
            res = T.transfer(edge, ls, lp, gp, reg)
            r = [x["r"] for x in res["rows"] if x["region"] == "pool" and np.isfinite(x["r"])]
            assert len(r) == res["n_pool"] and len(r) >= 3, (rt, len(r))
            assert abs(np.median(r) - 1.0) < 1e-12, (rt, np.median(r))


def test_transfer_spread():
    """spread reported with its expected value; gate decided (pass/fail with >= 3 pool modules, else cannot be assessed)."""
    study = {}
    for ln in STUDY.read_text().splitlines():
        if ln and not ln.startswith("#") and not ln.startswith("rt\t"):
            v = ln.split("\t")
            study[(v[0], v[1], v[2])] = (float(v[5]), v[8])      # spread, verdict
    for rt in RTS:
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            for f in ("edge_23701.tsv", "lms_scan.tsv", f"lms_prod_{rt}.tsv"):
                shutil.copy(T.OUT / f, td / f)
            res = T.build(rt, out_dir=td, prior=False)
        sp, verdict = study[(rt, "primary", "pool")]
        assert res["verdict"] == verdict and abs(res["spread"] - sp) < 1e-5, (rt, res["verdict"], verdict, res["spread"], sp)
        assert np.isfinite(res["spread"]) and np.isfinite(res["expected"]) and res["expected"] > 0
        assert res["verdict"] == ("pass" if res["spread"] <= T.GATE else "fail"), res["verdict"]
        print(f"  {rt}: n={res['n_pool']} spread={res['spread']:.4f} expected={res['expected']:.4f} gate {res['verdict']}")
    assert _synthetic(n=2)["verdict"] == "cannot be assessed"


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
