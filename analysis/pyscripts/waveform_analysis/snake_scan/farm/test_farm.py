#!/usr/bin/env python3
"""
test_farm.py PRODUCTS_DIR MERGED_DIR   -- compare farm products / farm_merge output of local files with the wave-1 results in SS/out
(the local 23701 splits 00000-00080 step 10 and the 23744 splits 00000/00030/00060).  Prints one line per check:
    PASS|FAIL name: numbers (tolerance)
MERGED_DIR must come from `farm_merge.py --in PRODUCTS_DIR --out MERGED_DIR --lms-runs 23701` (primary L_scan) with products made
by `farm_job.sh ... --keep-every 1` (the exact configuration; the default keep-every 2 changes the histograms by sampling).
"""
import sys
from pathlib import Path

import numpy as np

FARM = Path(__file__).resolve().parent
SS = FARM.parent
for p in (str(FARM), str(SS), str(SS.parents[1])):
    if p not in sys.path:
        sys.path.insert(0, p)
import farm_merge as FM        # noqa: E402
import snake_scan as SSN       # noqa: E402
import transfer as TR          # noqa: E402

PROD, MER = Path(sys.argv[1]), Path(sys.argv[2])
OUT = SS / "out"
LMS_NPZ = SS.parents[1] / "output" / "snake_scan" / "lms_npz"
res = []


def rep(name, ok, msg):
    res.append(ok)
    print(("PASS " if ok else "FAIL ") + f"{name}: {msg}")


def table(path):
    rows, _ = TR.read_tsv(path)
    return rows


prods = [p for p in FM.index(PROD) if p.mode == "scan"]
# 1 dwell rows, exact
w1 = SSN.read_dwell_tsv(OUT / "dwell_23701.tsv")
bad, n = 0, 0
for p in prods:
    mine = SSN.blocks_to_rows(p.name, p.blocks(), p.epics())
    ref = [r for r in w1 if r["file"] == p.name]
    n += len(ref)
    ok = len(mine) == len(ref) and all(
        all(a[k] == b[k] for k in ("ev_first", "ev_last", "module", "n_events", "kind", "ts_first", "ts_last")) and
        abs(a["share"] - b["share"]) < 5e-5 and abs(a["ptrans_x"] - b["ptrans_x"]) < 5e-5 and abs(a["ptrans_y"] - b["ptrans_y"]) < 5e-5
        for a, b in zip(mine, ref))
    bad += not ok
rep("dwell rows", bad == 0, f"{len(prods)} files, {n} wave-1 rows, mismatching files {bad} (exact event ranges/modules/kinds, share/ptrans 5e-5)")
# 2 LMS/alpha series, exact against the wave-1 per-split cache
bad, tot = 0, 0
for p in prods:
    split = p.name.split(".")[-1]
    ref = np.load(LMS_NPZ / f"lms_{p.run}_{split}.npz", allow_pickle=False)
    for k in ("lms_I", "lms_ref", "alpha_ref", "lms_ts", "alpha_ts", "lms_nch", "lms_pk_med", "ped"):
        tot += 1
        bad += not np.array_equal(p["L_" + k], ref[k], equal_nan=True)
rep("LMS/alpha/ped series", bad == 0, f"{tot} arrays compared with output/snake_scan/lms_npz, differing: {bad} (exact)")
# 3 edge table
a, b = table(MER / "edge_scan.tsv"), table(OUT / "edge_23701.tsv")
A = {r["module"]: r for r in a}
worst = dict(s=0, s_err=0, n=0, pile=0, sc=0)
for r in b:
    m = A.get(r["module"])
    if m is None:
        worst = None
        break
    worst["s"] = max(worst["s"], abs(float(m["s"]) - float(r["s"])))
    worst["s_err"] = max(worst["s_err"], abs(float(m["s_err"]) - float(r["s_err"])))
    worst["sc"] = max(worst["sc"], abs(float(m["s_central"]) - float(r["s_central"])))
    worst["n"] = max(worst["n"], abs(int(m["n_single"]) - int(r["n_single"])) / int(r["n_single"]))
    worst["pile"] = max(worst["pile"], abs(int(m["n_pileup"]) - int(r["n_pileup"])) / int(r["n_pileup"]))
rep("edge table", worst is not None and worst["s"] < 5e-5 and worst["s_err"] < 2e-5 and worst["sc"] < 1e-4 and worst["n"] < 1e-3 and worst["pile"] < 1e-3,
    "modules missing" if worst is None else f"{len(b)} modules: max|d s| {worst['s']:.2e}, max|d s_err| {worst['s_err']:.2e}, "
    f"max|d s_central| {worst['sc']:.2e}, max rel d n_single {worst['n']:.1e}, n_pileup {worst['pile']:.1e} "
    "(tolerances s 5e-5, s_err 2e-5, s_central 1e-4, counts 1e-3: the first ~100 large-peak events of a file/dwell are not used by the running-module rule)")
# 4 L_scan (23701 alone) and transfer
for nm in ("lms_scan.tsv",):
    a, b = table(MER / nm), table(OUT / nm)
    A = {r["module"]: r for r in a}
    d = max(abs(float(A[r["module"]]["L"]) - float(r["L"])) / float(r["L"]) for r in b)
    e = max(abs(float(A[r["module"]]["err"]) - float(r["err"])) / float(r["err"]) for r in b)
    rep(nm, len(a) == len(b) and d < 1e-6 and e < 1e-6, f"{len(b)} modules, max rel d L {d:.1e}, d err {e:.1e} (1e-6)")
for rt in ("2p24GeV_H2", "3p49GeV_H2"):
    f = MER / f"transfer_{rt}.tsv"
    if not f.exists():
        rep(f"transfer_{rt}", False, "not produced")
        continue
    a, b = table(f), table(OUT / f"transfer_{rt}.tsv")
    cols = ("s", "s_err", "L_scan", "L_prod", "g_T", "sigma_T", "g_pool", "r")
    worst = {c: 0.0 for c in cols}
    for ra, rb in zip(a, b):
        assert ra["module"] == rb["module"]
        for c in cols:
            x, y = float(ra[c]), float(rb[c])
            if np.isfinite(y):
                worst[c] = max(worst[c], abs(x - y) if c in ("s", "s_err") else abs(x - y) / max(abs(y), 1e-12))
    ok = len(a) == len(b) and worst["s"] < 5e-5 and worst["s_err"] < 2e-5 and all(worst[c] < (2e-3 if c == "sigma_T" else 1e-4) for c in cols[2:])
    rep(f"transfer_{rt}", ok, f"{len(b)} rows, max |d| of s, s_err: {worst['s']:.1e}, {worst['s_err']:.1e}; max rel d of L_scan, L_prod, g_T, "
        f"sigma_T, g_pool, r: " + ", ".join(f"{worst[c]:.1e}" for c in cols[2:]) + " (abs 5e-5/2e-5, rel 1e-4; sigma_T 2e-3 because s_err is printed with 5 decimals)")
    pa, pb = table(MER / f"prior_{rt}.tsv"), table(OUT / f"prior_{rt}.tsv")
    wp = max((abs(float(x["g"]) - float(y["g"])) for x, y in zip(pa, pb)), default=0)
    rep(f"prior_{rt}", len(pa) == len(pb) and wp < 1e-4, f"{len(pb)} rows, max |d g| {wp:.1e}")
# 5 cosmic
cos = [p for p in FM.index(PROD) if p.mode == "cosmic"]
if cos and (MER / "mip_23744.tsv").exists():
    a, b = table(MER / "mip_23744.tsv"), table(OUT / "mip_23744.tsv")
    worst = 0.0
    for ra, rb in zip(a, b):
        for c in rb:
            if c == "group":
                continue
            x, y = float(ra[c]), float(rb[c])
            worst = max(worst, abs(x - y) / max(abs(y), 1e-12))
    rep("mip_23744.tsv", len(a) == len(b) and worst < 1e-3, f"{len(b)} groups, max rel diff of all numeric columns {worst:.1e} (1e-3: the 4-digit printout)")
    za, zb = np.load(MER / "lms_23744.npz"), np.load(OUT / "lms_23744.npz")
    bad = [k for k in zb.files if k != "files" and not (np.array_equal(za[k], zb[k], equal_nan=True) if za[k].dtype.kind in "fc"
                                                          else np.array_equal(za[k], zb[k]))]
    rep("lms_23744.npz", not bad, f"{len(zb.files)} arrays compared, differing: {bad or 'none'} (exact)")
print(f"{sum(res)}/{len(res)} checks passed")
sys.exit(0 if all(res) else 1)
