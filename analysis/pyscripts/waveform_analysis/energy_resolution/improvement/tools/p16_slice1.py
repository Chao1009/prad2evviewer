#!/usr/bin/env python3
"""P16 slice 1: does the readout threshold of the modules around a cluster change its energy response and width?

python p16_slice1.py <report root> [--rts 0p73GeV_H2 ...] [--tsv out.tsv]

<report root>/<run type>/thr_clusters.npz from a report with --set THR_DUMP=True (control-gate path, default arm).
Per cluster T_w = sum_m w(d_m) T_m / sum_m w(d_m) over the usable PbWO4 modules m of the 5 x 5 window (seed excluded), T_m = TET x k_m (MeV),
w(d) = P(1 - P) with P the share of window modules at distance d (module sizes, bins of 0.125) that have a stored hit (peak clusters of
the run type): the weight of the modules that sit near threshold.  Classes: lo = T_w below its 25 % quantile, hi = above its 75 %, quantiles within 25 mm rings of the seed radius.
(a) per split and class, fit_path in sample (e), out of sample (e_oos), before gains (er: mu_none) and with gains only (eg = sum g me, no position
correction), pool / edge / centre; T_m in calibrated MeV (x the run's gain of the module); hi - lo and
class - all paired per split, mean +- SE (D3).  (b) within-seed slope: per split, r = E / E'(seed) - 1 (%) of in-peak clusters regressed on
T_w after subtracting each seed's mean (r, T_w and the nuisance u, v, u^2, v^2, uv), mean +- SE over splits; none = er, gain = eg, corr = e; *_thr2 / *_near2: joint fit with T_near, the same mean over the
always-read modules (P >= 0.95) weighted by their mean energy fraction (calibration-scale confound).  Mechanism rows: within-seed
slopes of the recorded energy of the modules at d >= 1.25 (ring, MeV) and of the number of absent modules with w_thr > 0.05 on T_w.
"""
import argparse, math, sys
import numpy as np
import p16_common as pc
bm = pc.bm

Q = ("res_cb", "mu_ep", "tail_lo_f", "tail_hi_f", "tail_lo3", "tail_hi3", "eff_cut", "leak_cut")
DB = np.arange(0.0, 3.6, 0.125)


def weights(z, d, T, ep):
    """Per d bin: w_thr = P (1 - P) (P = share of window modules with a stored hit) and w_near = mean energy fraction of the module where
    P >= 0.95 (0 elsewhere), from the peak clusters; also P and the robust width s of E / E'."""
    r = z["e"] / ep - 1
    s = 1.4826 * np.median(np.abs(r[z["used"]] - np.median(r[z["used"]])))
    pk = z["used"] & (np.abs(r - np.median(r[z["used"]])) < 2 * s)
    ok = (z["win_m"] >= 0) & np.isfinite(T) & pk[:, None]
    ib = np.digitize(d, DB) - 1
    P, F = np.full(DB.size - 1, np.nan), np.zeros(DB.size - 1)
    me = member_energy(z)
    for b in range(DB.size - 1):
        k = ok & (ib == b)
        if k.sum() > 100:
            P[b] = z["win_p"][k].mean()
            F[b] = np.mean(me[k] / z["er"][np.nonzero(k)[0]])
    return np.nan_to_num(P * (1 - P)), np.where(np.nan_to_num(P) >= 0.95, F, 0.0), P, s


def member_energy(z):
    """[n, 24] raw energy of each window module in its cluster (0 when not a member)."""
    out = np.zeros(z["win_m"].shape)
    for c, (mi, me) in enumerate(zip(z["mi"], z["me"])):
        hit = np.nonzero(z["win_m"][c][:, None] == mi[None, :])
        out[c, hit[0]] = me[hit[1]]
    return out


def wmean(W, Tm):
    return np.sum(W * np.nan_to_num(Tm), 1) / np.maximum(W.sum(1), 1e-12)


def demean(seed, X):
    """X [n, k] minus the mean of its seed."""
    u, inv = np.unique(seed, return_inverse=True)
    cnt = np.bincount(inv)
    return X - (np.stack([np.bincount(inv, weights=X[:, j]) for j in range(X.shape[1])], 1) / cnt[:, None])[inv]


def slope(z, sel, regs, ep, key, u, v):
    """Within-seed OLS coefficients (%/MeV) of r = 100 (E / E' - 1) on the regressors regs (list of [n] arrays), nuisance u, v, u^2, v^2, uv."""
    r = 100 * (z[key][sel] / ep[sel] - 1)
    X = np.column_stack([x[sel] for x in regs] + [u[sel], v[sel], u[sel] ** 2, v[sel] ** 2, u[sel] * v[sel]])
    D = demean(z["seed"][sel], np.column_stack([r, X]))
    return np.linalg.lstsq(D[:, 1:], D[:, 0], rcond=None)[0][:len(regs)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--rts", nargs="+", default=["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "3p49GeV_H2_old"])
    ap.add_argument("--tsv")
    a = ap.parse_args()
    out = []
    for rt in a.rts:
        z, geo = pc.load(f"{a.root}/{rt}/thr_clusters.npz"), pc.geometry(rt)
        T = pc.thresholds(rt)
        G = pc.gains(a.root, rt, geo["x"].size, list(geo["name"]))
        z["eg"] = pc.egain(z, G)
        m = z["win_m"].astype(np.int64)
        mm = np.where(m >= 0, m, 0)
        gm = np.empty(m.shape)
        for run, g in G.items():
            s_ = z["run"] == run
            gm[s_] = g[mm[s_]]
        Tm = np.where(m >= 0, T[mm] * gm, np.nan)      # threshold in calibrated MeV: TET x k x gain
        d, _ = pc.window_d(z, geo)
        ep_mod = bm.e_elastic(bm.theta_deg(geo["x"], geo["y"], geo["hz"]), geo["ebeam"], pc.consts(geo)[2])
        ep = ep_mod[z["seed"].astype(np.int64)]
        w, wn, P, s = weights(z, d, Tm, ep)
        ib = np.clip(np.digitize(d, DB) - 1, 0, DB.size - 2)
        Tw = wmean(np.where(np.isfinite(Tm), w[ib], 0.0), Tm)            # near-threshold modules
        Tn = wmean(np.where(np.isfinite(Tm), wn[ib], 0.0), Tm)           # always-read neighbours, energy weighted (scale confound)
        ring = np.digitize(np.hypot(geo["x"], geo["y"])[z["seed"].astype(np.int64)], np.arange(25.0, 200.0, 25.0))
        lo, hi = np.zeros(Tw.size, bool), np.zeros(Tw.size, bool)        # quartiles within 25 mm rings of the seed radius (matched in radius)
        for rg_ in np.unique(ring):
            k = ring == rg_
            q25, q75 = np.percentile(Tw[k & z["used"]], (25, 75)) if (k & z["used"]).sum() > 100 else (np.nan, np.nan)
            lo |= k & (Tw < q25)
            hi |= k & (Tw > q75)
        cls = dict(all=np.ones(Tw.size, bool), lo=lo, hi=hi)
        q25, q75 = np.percentile(Tw[z["used"]], (25, 75))
        far = (d >= 1.25) & np.isfinite(Tm)                              # the ring where modules sit near threshold
        z["ring"] = np.where(far, member_energy(z) * gm, 0.0).sum(1)
        z["nabs"] = (far & ~z["win_p"] & (w[ib] > 0.05)).sum(1).astype(float)
        u, v = bm.uv(z, geo)
        edge = (np.abs(u) > bm.EDGE_UV) | (np.abs(v) > bm.EDGE_UV)
        reg = dict(pool=np.ones(Tw.size, bool), edge=edge, centre=~edge)
        print(f"== {rt}: {z['e'].size} clusters, T_w quartiles {q25:.2f} / {q75:.2f} MeV, mean lo {Tw[cls['lo'] & z['used']].mean():.2f} "
              f"hi {Tw[cls['hi'] & z['used']].mean():.2f}; P(d) {np.round(P[4:24:2], 2)}", flush=True)
        out.append(dict(rt=rt, region="-", cls="-", q="Tw_lo_hi", mean=float(Tw[cls['lo'] & z['used']].mean()),
                        se=float(Tw[cls['hi'] & z['used']].mean()), n=0))
        val = {}                                                   # (region, class, quantity) -> {split key: value}
        sl = {}
        keys = sorted(set(zip(z["run"].tolist(), z["split"].tolist())))
        for k in keys:
            sp = (z["run"] == k[0]) & (z["split"] == k[1])
            for rg, rm in reg.items():
                for c, cm in cls.items():
                    sel = sp & rm & cm
                    for tag, ek, oo in (("", "e", False), ("_oos", "e_oos", True), ("_none", "er", False), ("_gain", "eg", False)):
                        if oo and not np.isfinite(z["e_oos"][sel]).any():
                            continue
                        f = pc.fit(z, sel, geo, rt, ek, oo)
                        if f is None or f["n_range"] < bm.REGION_MIN_FIT:
                            continue
                        for q in (Q if tag == "" else ("res_cb", "mu_ep")):
                            val.setdefault((rg, c, q + tag), {})[k] = f[q]
            r = z["e"] / ep - 1
            pk = sp & z["used"] & (np.abs(r - np.median(r[z["used"]])) < 3 * s)
            for key in ("ring", "nabs"):                                  # mechanism: recorded ring energy (MeV) and absent near-threshold modules per MeV of T_w
                X = np.column_stack([Tw[pk], u[pk], v[pk], u[pk] ** 2, v[pk] ** 2, u[pk] * v[pk]])
                D = demean(z["seed"][pk], np.column_stack([z[key][pk], X]))
                sl.setdefault(key, []).append(np.linalg.lstsq(D[:, 1:], D[:, 0], rcond=None)[0][0])
            for key in ("er", "eg", "e"):
                sl.setdefault(key, []).append(slope(z, pk, [Tw], ep, key, u, v)[0])
                b2 = slope(z, pk, [Tw, Tn], ep, key, u, v)
                sl.setdefault(key + "_thr2", []).append(b2[0])
                sl.setdefault(key + "_near2", []).append(b2[1])
        for (rg, c, q), d_ in sorted(val.items()):
            mn, se, n = pc.mse(list(d_.values()))
            out.append(dict(rt=rt, region=rg, cls=c, q=q, mean=mn, se=se, n=n))
            if c == "hi" and (rg, "lo", q) in val:
                lo = val[rg, "lo", q]
                mn, se, n = pc.mse([d_[k] - lo[k] for k in d_ if k in lo])
                out.append(dict(rt=rt, region=rg, cls="hi-lo", q=q, mean=mn, se=se, n=n))
            if c in ("lo", "hi") and (rg, "all", q) in val:
                al = val[rg, "all", q]
                mn, se, n = pc.mse([d_[k] - al[k] for k in d_ if k in al])
                out.append(dict(rt=rt, region=rg, cls=f"{c}-all", q=q, mean=mn, se=se, n=n))
        for key, v_ in sl.items():
            mn, se, n = pc.mse(v_)
            if key in ("ring", "nabs"):
                out.append(dict(rt=rt, region="pool", cls="within_seed", q=dict(ring="ring_MeV_per_MeV", nabs="n_absent_per_MeV")[key], mean=mn, se=se, n=n))
                continue
            out.append(dict(rt=rt, region="pool", cls="within_seed", q=f"slope_{dict(er='none', eg='gain', e='corr')[key.split('_')[0]]}{('_' + key.split('_')[1]) if '_' in key else ''}_pct_per_MeV", mean=mn, se=se, n=n))
        for r_ in out:
            if r_["rt"] == rt and (r_["cls"] in ("hi-lo", "within_seed")) and r_["region"] in ("pool", "edge", "centre"):
                print(f"  {r_['region']:6s} {r_['cls']:11s} {r_['q']:24s} {r_['mean']:+.4f} +- {r_['se']:.4f} (n {r_['n']})", flush=True)
    if a.tsv:
        with open(a.tsv, "w") as f:
            f.write("# P16 slice 1 (p16_slice1.py); Tw_lo_hi row: mean = <T_w> lo class, se column = <T_w> hi class (MeV)\n")
            f.write("run_type\tregion\tclass\tquantity\tmean\tse\tn\n")
            for r_ in out:
                f.write(f"{r_['rt']}\t{r_['region']}\t{r_['cls']}\t{r_['q']}\t{r_['mean']:.6g}\t{r_['se']:.6g}\t{r_['n']}\n")


if __name__ == "__main__":
    main()
