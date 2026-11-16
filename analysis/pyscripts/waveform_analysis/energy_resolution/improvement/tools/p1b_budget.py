#!/usr/bin/env python3
"""p1b_budget.py <out dir> <scan npz> <fit_coef.npy> <control thr_clusters.npz> [label=thr_clusters.npz ...]: AMP P1b cluster-level budget.

Control = the integral arm (thr_clusters.npz of a report with THR_DUMP, the control-gate pool clusters, default arm in sample).  Elastic clusters:
used by the control's gain solve and 0.93 < e / E' < 1.07 in both arms (E' from the cluster position).  Tables in <out dir>:
  pairs.tsv      per arm: matched clusters (split, ev, seed), robust sigma (1.4826 MAD) of e / E' for the control and the arm, of the paired
                 e_arm / e_control - 1, and the quadrature excess sqrt(s_arm^2 - s_control^2)
  classes.tsv    sigma of the paired difference by cluster class (pile-up flag, members, seed fraction, summed energy of the three largest members)
  predict.tsv    the hit-level variance v(A) = (n / A)^2 + a / A + c^2 of the raw peak (fit_coef.npy row 0) propagated member by member to the
                 elastic clusters of the control: sqrt(mean_cluster sum_m [(c_m n)^2 + e_m c_m a + e_m^2 c^2]) / E', by component
Member amplitude A_m = e_m / s_m with s_m the module median of h_e / h_a (single-pulse hits, h_a > 200).
"""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
import energy_resolution_benchmark as bm        # noqa: E402

WIN = 0.07


def rsig(v):
    return 1.4826 * np.median(np.abs(v - np.median(v)))


def load(p):
    z = np.load(p)
    n = z["e"].size
    off = np.r_[0, np.cumsum(z["nmem"])]
    ep = bm.e_elastic(np.degrees(np.arctan(np.hypot(z["x"], z["y"]) / float(z["hz"]))), float(z["ebeam"]), bm.TARGET_MASS["H2"])
    return {"key": (z["split"].astype(np.int64) << 48) + (z["ev"].astype(np.int64) << 16) + z["seed"].astype(np.int64), "e": z["e"], "r": z["e"] / ep, "pu": z["pu"],
            "used": z["used"], "mi": z["mi"], "me": z["me"], "off": off, "ep": ep, "n": n}


def main(out, scan, coef, ctl, *arms):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    C = load(ctl)
    inwin = C["used"] & (np.abs(C["r"] - np.median(C["r"][C["used"]])) < WIN)
    rows, cls = [], []
    for a in arms:
        lab, p = a.split("=", 1)
        B = load(p)
        _, ia, ib = np.intersect1d(C["key"], B["key"], return_indices=True)
        ok = inwin[ia] & (np.abs(B["r"][ib] - np.median(B["r"][B["used"]])) < WIN)
        ia, ib = ia[ok], ib[ok]
        d = B["e"][ib] / C["e"][ia] - 1
        sc, sb = rsig(C["r"][ia]), rsig(B["r"][ib])
        rows.append([lab, ia.size, 100 * sc, 100 * sb, 100 * np.sqrt(max(sb ** 2 - sc ** 2, 0)), 100 * rsig(d), 100 * np.std(d[np.abs(d) < 0.2]), 100 * np.median(d),
                     np.corrcoef(C["r"][ia], d)[0, 1] if ia.size > 10 else np.nan])
        nm = np.diff(C["off"])[ia]
        top3 = np.array([np.sort(C["me"][C["off"][i]:C["off"][i + 1]])[-1] / C["e"][i] for i in ia])
        for cname, sel in (("pu flag off", ~C["pu"][ia]), ("pu flag on", C["pu"][ia]), ("members <= 9", nm <= 9), ("members 10-14", (nm > 9) & (nm <= 14)),
                           ("members >= 15", nm >= 15), ("seed fraction < 0.6", top3 < 0.6), ("seed fraction 0.6-0.8", (top3 >= 0.6) & (top3 < 0.8)),
                           ("seed fraction >= 0.8", top3 >= 0.8)):
            if sel.sum() > 100:
                cls.append([lab, cname, int(sel.sum()), 100 * rsig(d[sel]), 100 * rsig(B["r"][ib][sel]), 100 * rsig(C["r"][ia][sel])])
    bm.write_tsv(out / "pairs.tsv", ("arm", "n", "sig_control_pct", "sig_arm_pct", "excess_quadrature_pct", "sig_pair_diff_pct", "rms_pair_diff_pct", "median_diff_pct",
                                    "corr_control_vs_diff"), rows)
    bm.write_tsv(out / "classes.tsv", ("arm", "class", "n", "sig_pair_diff_pct", "sig_arm_pct", "sig_control_pct"), cls)
    # prediction
    z = np.load(scan)
    m = z["h_mod"].astype(int)
    ha, he = z["h_a"].astype(np.float64), z["h_e"].astype(np.float64)
    one = (z["h_k"] == 1) & ~z["h_extra"] & (ha > 200)
    s_m = np.full(z["mod_name"].size, np.nan)
    o = np.argsort(m[one], kind="stable")
    mm, rr = m[one][o], (he[one] / ha[one])[o]
    b = np.r_[0, np.cumsum(np.bincount(mm, minlength=s_m.size))]
    for k in range(s_m.size):
        if b[k + 1] - b[k] > 30:
            s_m[k] = np.median(rr[b[k]:b[k + 1]])
    s_m = np.where(np.isfinite(s_m), s_m, np.nanmedian(s_m))
    n, a, c = (float(v) for v in (np.sqrt(np.load(coef)[0][0]), np.load(coef)[0][1], np.sqrt(np.load(coef)[0][2])))
    idx = np.flatnonzero(inwin)
    comp = np.zeros((idx.size, 3))
    for q, i in enumerate(idx):
        lo, hi = C["off"][i], C["off"][i + 1]
        me, sm = C["me"][lo:hi], s_m[C["mi"][lo:hi]]
        comp[q] = [np.sum((sm * n) ** 2), np.sum(me * sm * a), np.sum(me ** 2) * c ** 2]
    E2 = (C["ep"][idx]) ** 2
    pred = [["white noise n/A", n, 100 * np.sqrt(np.mean(comp[:, 0] / E2))], ["stochastic a/A", a, 100 * np.sqrt(np.mean(comp[:, 1] / E2))],
            ["constant c", c, 100 * np.sqrt(np.mean(comp[:, 2] / E2))], ["total", np.nan, 100 * np.sqrt(np.mean(comp.sum(1) / E2))]]
    # analytic single-sample noise from the measured pedestal statistics (noise_summary.tsv next to the fit coefficients)
    ns = {r["quantity"]: float(r["value"]) for r in bm.tsv_rows(Path(coef).parent / "noise_summary.tsv")}
    psd, exc = ns["psd_adc_median"], ns["sigma_excess_white_adc (psd sqrt(1-1/sq))"]
    seed_i = np.array([np.argmax(C["me"][C["off"][i]:C["off"][i + 1]]) + C["off"][i] for i in idx])
    A_seed = C["me"][seed_i] / s_m[C["mi"][seed_i]]
    for lab, sg in (("raw sample noise (pedestal std)", psd), ("noise in excess of the matched filter (psd sqrt(1 - 1/sumsq))", exc)):
        allm = np.array([np.sum((s_m[C["mi"][C["off"][i]:C["off"][i + 1]]] * sg) ** 2) for i in idx])
        seedm = (s_m[C["mi"][seed_i]] * sg) ** 2
        pred += [[lab + ", all members", sg, 100 * np.sqrt(np.mean(allm / E2))], [lab + ", seed only", sg, 100 * np.sqrt(np.mean(seedm / E2))]]
    pred += [["median seed amplitude (ADC)", np.median(A_seed), np.nan], ["median seed sigma_noise / A (percent)", psd, 100 * np.median(psd / A_seed)]]
    bm.write_tsv(out / "predict.tsv", ("component", "hit_parameter", "cluster_sigma_pct"), pred, comment="# raw-peak excess over the matched filter propagated to the elastic clusters")
    for f in ("pairs.tsv", "classes.tsv", "predict.tsv"):
        print(open(out / f).read())


if __name__ == "__main__":
    main(*sys.argv[1:])
