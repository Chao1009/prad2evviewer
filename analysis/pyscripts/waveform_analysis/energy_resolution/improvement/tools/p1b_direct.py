#!/usr/bin/env python3
"""p1b_direct.py <out dir> <scan npz> <control thr_clusters.npz> <rawpk corr tables npz (column@h_e)> [col=tablenpz ...]: AMP P1b direct cluster-level budget of a peak column.

Gain-solve free.  For every member of an elastic control cluster (integral arm, used by the gain solve, |e / E' - median| < 7 %) the hit (event, module,
energy) is looked up in the scan; rho = column / h_e (ADC per integral MeV), u = rho / (module median of rho over single-pulse hits with h_a > 200) - 1,
and D = sum_m e_m u_m / sum_m e_m is the relative energy shift the column would give the cluster against the integral (the integral's own error is
common to both).  var(D) (mean removed) is the excess variance of the column over the integral at cluster level, as the report arms measure it with the gain solve in
between.  Variants of u: phase table (u' = rho / (module median x curve[module, phase])), phase and amplitude tables, hits with a pile-up pulse
(h_k >= 2 or extra) set to u = 0, members below / above E MeV set to 0, and the independent-member prediction sum_m e_m^2 sigma^2(A_m) / E^2 from the binned
hit scatter (sigma(A) robust, single-pulse hits, phase curve removed or not).  Standard errors: 40 event blocks, 200 resamples.
Output: direct.tsv (variance in percent^2, mean over clusters of D^2 with |D| < 10 % kept; robust variants in the md), clusters.npz (D per variant)."""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
import energy_resolution_benchmark as bm        # noqa: E402

COLS = {"rawpk": "h_apk_raw", "sum3": "h_s3", "sum4": "h_s4", "sum5": "h_s5", "sum7": "h_s7"}


def rsig(v):
    return 1.4826 * np.median(np.abs(v - np.median(v)))


def main(out, scan, ctl, raw_tab, *extra):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    tabs = {"rawpk": raw_tab} | dict(e.split("=", 1) for e in extra)
    z, c = np.load(scan), np.load(ctl)
    names = z["mod_name"]
    t50 = bm.template_columns(bm.PYDIR / str(z["bank"]), names)[0]
    hev, hm = z["h_ev"].astype(np.int64), z["h_mod"].astype(np.int64)
    he, ha, hk, hx = z["h_e"].astype(np.float64), z["h_a"].astype(np.float64), z["h_k"], z["h_extra"]
    x = ((z["h_t"].astype(np.float64) - t50[hm]) % bm.DT) / bm.DT
    key = hev * 4096 + hm
    order = np.argsort(key, kind="stable")
    ks = key[order]
    # elastic control clusters and their members
    ep = bm.e_elastic(np.degrees(np.arctan(np.hypot(c["x"], c["y"]) / float(c["hz"]))), float(c["ebeam"]), bm.TARGET_MASS["H2"])
    r = c["e"] / ep
    used = c["used"]
    ok = used & (np.abs(r - np.median(r[used])) < 0.07)
    off = np.r_[0, np.cumsum(c["nmem"])]
    cid = np.repeat(np.arange(c["e"].size), c["nmem"])
    mkey = c["ev"].astype(np.int64)[cid] * 4096 + c["mi"].astype(np.int64)
    sel = ok[cid]
    cid, mkey, me = cid[sel], mkey[sel], c["me"][sel]
    lo, hi = np.searchsorted(ks, mkey, "left"), np.searchsorted(ks, mkey, "right")
    hit = np.full(me.size, -1)
    best = np.full(me.size, np.inf)
    for d in range(int((hi - lo).max())):
        i = np.flatnonzero(lo + d < hi)
        cand = order[lo[i] + d]
        dd = np.abs(he[cand] - me[i])
        b = dd < best[i]
        hit[i[b]], best[i[b]] = cand[b], dd[b]
    good = (hit >= 0) & (best < 1e-3 * np.maximum(me, 1.0))
    print(f"members {me.size}, matched {good.mean():.4f} (clusters {ok.sum()})")
    h = np.where(good, hit, 0)
    mod = hm[h]
    pile = (hk[h] >= 2) | hx[h]
    one = (hk == 1) & ~hx & (ha > 200)
    nb = int(ok.sum())
    ids = np.cumsum(ok) - 1
    cl = ids[cid]                                           # cluster row 0..nb-1 of every member
    E = np.bincount(cl, weights=me, minlength=nb)
    rows, dump = [], {}
    ev = c["ev"][ok]
    blk = np.minimum(((ev - ev.min()) * 40 // (ev.max() - ev.min() + 1)).astype(int), 39)
    order_b = [np.flatnonzero(blk == b) for b in range(40)]
    rng = np.random.default_rng(1)
    boots = [np.concatenate([order_b[b] for b in rng.integers(0, 40, 40)]) for _ in range(200)]

    def stat(D, label, comp=""):
        D = np.where(np.abs(D) < 0.1, D, 0.0)
        v = 1e4 * np.var(D)                                 # variance: a global offset is absorbed by the gain solve and is not a resolution term
        se = np.std([1e4 * np.var(D[ix]) for ix in boots])
        rows.append([label, v, se, 100 * rsig(D), f"mean {100 * D.mean():.3f} %"])
        dump[label] = D
        return v

    for name, col in COLS.items():
        if name not in tabs:
            continue
        T = np.load(tabs[name])
        cv = z[col].astype(np.float64)
        rho = np.where(np.isfinite(cv) & (he > 0), cv / he, np.nan)
        med = np.full(names.size, np.nan)
        o = np.argsort(hm[one & np.isfinite(rho)], kind="stable")
        mm, rr = hm[one & np.isfinite(rho)][o], rho[one & np.isfinite(rho)][o]
        bb = np.r_[0, np.cumsum(np.bincount(mm, minlength=names.size))]
        for k in range(names.size):
            if bb[k + 1] - bb[k] > 30:
                med[k] = np.median(rr[bb[k]:bb[k + 1]])
        med = np.where(np.isfinite(med), med, np.nanmedian(med))
        u = rho[h] / med[mod] - 1
        fin = good & np.isfinite(u)
        ph = T["curve"][mod, np.minimum((x[h] * 16).astype(int), 15)]
        ga = np.log(T["g_a"])
        gA = np.interp(np.log(np.maximum(cv[h], 1e-3)), ga, T["g"])

        def D_of(uu, mask=None):
            w = np.where(fin & (True if mask is None else mask), uu, 0.0) * me
            return np.bincount(cl, weights=w, minlength=nb) / E
        def parts(uu, label):
            w = np.where(fin, uu, 0.0) * me
            D = np.bincount(cl, weights=w, minlength=nb) / E
            diag = np.bincount(cl, weights=w ** 2, minlength=nb) / E ** 2
            ok2 = np.abs(D) < 0.1
            Dm = np.where(ok2, D - D[ok2].mean(), 0.0)
            rows.append([f"{label}: independent part (mean of the sum of squared member terms; holds member-level mean biases)", 1e4 * np.mean(np.where(ok2, diag, 0.0)),
                         np.std([1e4 * np.mean(np.where(ok2, diag, 0.0)[ix]) for ix in boots]), np.nan, ""])
            cr = np.where(ok2, Dm ** 2 - diag, 0.0)
            rows.append([f"{label}: cross terms = variance of D minus the independent part", 1e4 * np.mean(cr), np.std([1e4 * np.mean(cr[ix]) for ix in boots]), np.nan, ""])
            # residual first and second sampling-phase harmonics of the cluster shift, member weighted (as the P10 fit)
            th = 2 * np.pi * x[h]
            F = np.stack([np.bincount(cl, weights=np.where(fin, me, 0.0) * g(k * th), minlength=nb) / E for k in (1, 2) for g in (np.cos, np.sin)], axis=1)
            A = F[ok2]
            cf, *_ = np.linalg.lstsq(A, D[ok2], rcond=None)
            rows.append([f"{label}: removable by a member-weighted k = 1, 2 phase harmonic fit (in sample)", 1e4 * (np.var(D[ok2]) - np.var(D[ok2] - A @ cf)), np.nan, np.nan, ""])
        stat(D_of(u), f"{name}: all members")
        parts(u, f"{name}: all members")
        stat(D_of((u + 1) / ph - 1), f"{name}: phase table", "")
        stat(D_of((u + 1) / (ph * gA) - 1), f"{name}: phase + amplitude tables")
        parts((u + 1) / (ph * gA) - 1, f"{name}: phase + amplitude tables")
        if name == "rawpk":
            stat(D_of(u, ~pile), f"{name}: single-pulse members only (pile-up members u=0)")
            stat(D_of((u + 1) / (ph * gA) - 1, ~pile), f"{name}: phase + amplitude tables, single-pulse members only")
            for e0 in (100.0, 400.0, 1000.0):
                stat(D_of(u, me >= e0), f"{name}: members >= {e0:g} MeV only")
                stat(D_of(u, me < e0), f"{name}: members < {e0:g} MeV only")
            # independent-member prediction from the binned hit scatter of the single-pulse hits
            edges = np.array([30, 40, 50, 65, 80, 100, 130, 170, 220, 300, 400, 550, 750, 1000, 1400, 2000, 4000.0])
            s2 = []
            for a_lo, a_hi in zip(edges[:-1], edges[1:]):
                s = one & np.isfinite(rho) & (cv / np.nanmedian(med[hm]) >= a_lo) & (cv / np.nanmedian(med[hm]) < a_hi)
                s2.append(rsig(rho[s] / med[hm[s]] - 1) ** 2 if s.sum() > 500 else np.nan)
            cen = np.sqrt(edges[:-1] * edges[1:])
            s2 = np.array(s2)
            Am = cv[h] / np.nanmedian(rho)
            sig2 = np.interp(np.log(np.maximum(Am, 5.0)), np.log(cen)[np.isfinite(s2)], s2[np.isfinite(s2)])
            ind = np.bincount(cl, weights=np.where(fin & ~pile, (me * 1.0) ** 2 * sig2, 0.0), minlength=nb) / E ** 2
            rows.append(["rawpk: independent-member prediction (single-pulse members, hit scatter with phase)", 1e4 * ind.mean(), 1e4 * ind.std() / np.sqrt(nb), np.nan, ""])
    bm.write_tsv(out / "direct.tsv", ("variant", "var_pct2", "se_pct2", "robust_sigma_pct", "note"), rows)
    np.savez_compressed(out / "clusters.npz", ev=ev, **{k.replace(" ", "_"): v for k, v in dump.items()})
    print(open(out / "direct.tsv").read())


if __name__ == "__main__":
    main(*sys.argv[1:])
