#!/usr/bin/env python3
"""p1b_arms.py <out tsv> <arm dir root> <control arm> <arm>...: AMP P1b arm table on one common set of elastic clusters.

Per arm (report dir <root>/<arm>/2p24GeV_H2: summary.tsv, thr_clusters.npz of the control-gate pool clusters, default arm in sample): Crystal Ball
sigma/E of the report (g16, in sample; one split, no standard error), and on the clusters present in every arm (same split, event and seed) with
|e / E' - median| < 7 % in the control: robust sigma (1.4826 MAD) of e / E', its difference to the control, the difference of the variances, and the
standard errors of both by resampling 40 blocks of consecutive events (paired: the same blocks for every arm)."""
import csv
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
import energy_resolution_benchmark as bm        # noqa: E402


def rsig(v):
    return 1.4826 * np.median(np.abs(v - np.median(v)))


def load(p):
    z = np.load(p)
    ep = bm.e_elastic(np.degrees(np.arctan(np.hypot(z["x"], z["y"]) / float(z["hz"]))), float(z["ebeam"]), bm.TARGET_MASS["H2"])
    return {"key": (z["split"].astype(np.int64) << 48) + (z["ev"].astype(np.int64) << 16) + z["seed"].astype(np.int64), "ev": z["ev"], "r": z["e"] / ep, "used": z["used"]}


def cb(root, arm, path="g16"):
    r = {x["path"]: x for x in csv.DictReader(open(Path(root) / arm / "2p24GeV_H2" / "summary.tsv"), delimiter="\t")}
    return float(r[path]["res_cb"]), float(r[path]["yield_pct"])


def main(out, root, ctl, *arms):
    D = {a: load(Path(root) / a / "2p24GeV_H2" / "thr_clusters.npz") for a in (ctl,) + arms}
    keys = D[ctl]["key"]
    for a in arms:
        keys = np.intersect1d(keys, D[a]["key"])
    sel = {}
    for a, d in D.items():
        _, i, _ = np.intersect1d(d["key"], keys, return_indices=True)
        sel[a] = d["r"][i]
        if a == ctl:
            ev = d["ev"][i]
            used = d["used"][i]
    ok = used & (np.abs(sel[ctl] - np.median(sel[ctl][used])) < 0.07)
    ev = ev[ok]
    blk = np.minimum(((ev - ev.min()) * 40 // (ev.max() - ev.min() + 1)).astype(int), 39)
    order = [np.flatnonzero(blk == b) for b in range(40)]
    rng = np.random.default_rng(1)
    boots = [np.concatenate([order[b] for b in rng.integers(0, 40, 40)]) for _ in range(200)]
    s0 = {a: rsig(sel[a][ok]) for a in D}
    rows = []
    for a in (ctl,) + arms:
        sd, sv = [], []
        for ix in boots:
            sc, sa = rsig(sel[ctl][ok][ix]), rsig(sel[a][ok][ix])
            sd.append(sa - sc)
            sv.append(sa ** 2 - sc ** 2)
        res, yl = cb(root, a)
        rows.append([a, res, yl, int(ok.sum()), 100 * s0[a], 100 * (s0[a] - s0[ctl]), 100 * np.std(sd), 1e4 * (s0[a] ** 2 - s0[ctl] ** 2), 1e4 * np.std(sv),
                     100 * np.sqrt(max(s0[a] ** 2 - s0[ctl] ** 2, 0.0))])
    CON = [("phase table (vs integral)", "rawpk", "rawpk@corrhephase"), ("amplitude table after phase table", "rawpk@corrhephase", "rawpk@corrhephase+amp"),
           ("phase table (vs h_a)", "rawpk", "rawpk@corrphase"), ("raw -> 5-sample sum, no tables", "rawpk", "sum5"),
           ("raw -> 5-sample sum, both with phase + amplitude tables", "rawpk@corrhephase+amp", "sum5@corrhephase+amp"),
           ("raw -> 7-sample sum, both with tables", "rawpk@corrhephase+amp", "sum7@corrhephase+amp"),
           ("raw -> 4-sample sum, both with tables", "rawpk@corrhephase+amp", "sum4@corrhephase+amp"),
           ("raw -> CFD-anchored 7-sample sum, both with tables", "rawpk@corrhephase+amp", "cfd7@corrhephase+amp"),
           ("5-sample peak-anchored -> CFD-anchored 7-sample sum, both with tables", "sum5@corrhephase+amp", "cfd7@corrhephase+amp"),
           ("pile-up hits only (raw in kn)", "integral", "rawpk@kn"), ("single-pulse hits only (raw in k1)", "integral", "rawpk@k1"),
           ("hits >= 400 MeV only", "integral", "rawpk@above400"), ("hits < 400 MeV only", "integral", "rawpk@below400"),
           ("hits >= 1000 MeV only", "integral", "rawpk@above1000"), ("hits < 1000 MeV only", "integral", "rawpk@below1000"),
           ("integral -> raw (total)", "integral", "rawpk")]
    crow = []
    for lab, a, b in CON:
        if a not in sel or b not in sel:
            continue
        v = [rsig(sel[b][ok][ix]) ** 2 - rsig(sel[a][ok][ix]) ** 2 for ix in boots]
        crow.append([lab, a, b, 1e4 * (s0[b] ** 2 - s0[a] ** 2), 1e4 * np.std(v)])
    bm.write_tsv(str(out).replace(".tsv", "_contrasts.tsv"), ("contrast", "from", "to", "d_var_pct2", "d_var_se_pct2"), crow)
    print(open(str(out).replace(".tsv", "_contrasts.tsv")).read())
    bm.write_tsv(out, ("arm", "res_cb_g16_pct", "yield_pct", "n_clusters", "sig_robust_pct", "d_sig_pct", "d_sig_se_pct", "d_var_pct2", "d_var_se_pct2", "excess_quadrature_pct"), rows)
    print(open(out).read())


if __name__ == "__main__":
    main(*sys.argv[1:])
