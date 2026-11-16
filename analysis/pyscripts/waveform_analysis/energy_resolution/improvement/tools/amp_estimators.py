#!/usr/bin/env python3
"""amp_estimators.py <scan dir> <out dir> [label=phase_clusters.npz ...]: P1 estimator study on scans written with PEAK_COLUMNS (amplitude track).

Tables (tsv) in <out dir>:
  estimator_table.tsv  per column (integral = h_e, rawpk, parpk, parpk0 = parabola without the U2 subtraction): global scale s = median(h_e / column),
                       the column / h_e ratio (module medians: median and IQR over modules), the first harmonic of the sampling-phase curve
                       (column / h_a, each module divided by its median; pooled curve and the median over modules with >= PH_MIN hits)
  phase_curves.tsv     the pooled 16-bin curves behind it
  energy_bins.tsv      column / h_e in 8 energy bins (h_e / kE)
  overlap.tsv          h_k >= 2 primary hits against the signed time to the module's largest extra pulse: size of the U2 subtraction and the
                       ratio column / h_a (single-pulse median of the module = 1) with and without it
  cfd_vs_mf.tsv        h_tcfd - h_t per amplitude class after removing each module's median
  phase_true.tsv      (with label=dump pairs, PHASE_STUDY dumps of the reports) the phase dependence against the true energy as in P10: clusters with
                       0.95 < E/E' < 1.05, E/E' / median - 1 = sum over the three largest members f_j [c cos(k th_j) + s sin(k th_j)], th = 2 pi x / 256,
                       k = 1 (per unit f, permille), amplitude and standard error; the same with k = 1, 2 (rms of the fitted curve)
Sampling phase: x = ((h_t - t50) mod 4 ns) / 4 ns (the onset phase of the scan, mirrored).  Reference of the curves: h_a, the template-fit amplitude."""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
import energy_resolution_benchmark as bm        # noqa: E402

A_MIN, PH_MIN, NB = 200.0, 2000, 16             # ADC amplitude floor of the curves, hits per module for a per-module harmonic, phase bins
COLS = {"integral": "h_e", "rawpk": "h_apk_raw", "parpk": "h_apk_par", "parpk0": "h_apk_par0"}


def tsv(path, head, rows):
    bm.write_tsv(path, head, rows)


def harm(r, x):
    """First-harmonic amplitude (fraction), peak-to-peak of the binned curve, curve (bin medians) of the values r at phase x."""
    b = np.minimum((x * NB).astype(int), NB - 1)
    cur = np.array([np.median(r[b == k]) if (b == k).any() else np.nan for k in range(NB)])
    th = 2 * np.pi * (np.arange(NB) + 0.5) / NB
    ok = np.isfinite(cur)
    c = cur[ok] - cur[ok].mean()
    return 2 * abs(np.sum(c * np.exp(1j * th[ok]))) / ok.sum(), np.nanmax(cur) - np.nanmin(cur), cur


def by_module(m, n_mod):
    """Index arrays per module (sorted order of m)."""
    o = np.argsort(m, kind="stable")
    b = np.r_[0, np.cumsum(np.bincount(m, minlength=n_mod))]
    return [o[b[i]:b[i + 1]] for i in range(n_mod)]


def mod_median(v, groups, nmin=30):
    """Median of the finite values of v per module group; NaN below nmin values."""
    out = np.full(len(groups), np.nan)
    for i, g in enumerate(groups):
        w = v[g]
        w = w[np.isfinite(w)]
        if w.size > nmin:
            out[i] = np.median(w)
    return out


def phase_true(label, path):
    """Row of phase_true.tsv for a PHASE_STUDY dump."""
    z = np.load(path)
    ok = (z["r"] > 0.95) & (z["r"] < 1.05) & np.isfinite(z["x0"])
    y = z["r"][ok] / np.median(z["r"][ok]) - 1.0
    th = [2 * np.pi * np.nan_to_num(z[f"x{j}"][ok]) / 256 for j in range(3)]
    f = [np.where(np.isfinite(z[f"x{j}"][ok]), z[f"f{j}"][ok], 0.0) for j in range(3)]
    res = [label, int(ok.sum()), float(np.std(y)) * 1e3]
    for K in (1, 2):
        A = np.stack([sum(f[j] * g(k * th[j]) for j in range(3)) for k in range(1, K + 1) for g in (np.cos, np.sin)], axis=1)
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        se = np.sqrt(np.diag(np.linalg.inv(A.T @ A)) * np.var(y - A @ c))
        res += [1e3 * c[0], 1e3 * c[1], 1e3 * np.hypot(c[0], c[1]), 1e3 * se[0]] if K == 1 else [1e3 * np.std(A @ c)]
    return res


def main(src, out, *dumps):
    zs = [np.load(f) for f in sorted(Path(src).glob("*.npz"))]
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    d = {q: np.concatenate([z[q] for z in zs]) for q in ("h_mod", "h_t", "h_a", "h_k", "h_extra", "h_tcfd", *COLS.values())}
    names = zs[0]["mod_name"]
    t50 = bm.template_columns(bm.PYDIR / str(zs[0]["bank"]), names)[0]
    kE = float(zs[0]["beam_mev"]) / bm.E_REF
    m = d["h_mod"].astype(int)
    grp = by_module(m, len(names))
    x = ((d["h_t"].astype(np.float64) - t50[m]) % bm.DT) / bm.DT
    ha, he = d["h_a"].astype(np.float64), d["h_e"].astype(np.float64)
    one = (d["h_k"] == 1) & ~d["h_extra"] & (ha > A_MIN)
    pool = ~d["h_extra"] & (d["h_k"] > 0) & (he >= bm.SEED_MIN * kE)
    rows, curves, ebins = [], [], []
    edges = np.array([30, 100, 200, 400, 700, 1000, 1500, 2000, 1e9]) * kE
    for name, q in COLS.items():
        c = d[q].astype(np.float64)
        ok = pool & np.isfinite(c) & (c > 0)
        s = float(np.median(he[ok] / c[ok]))
        with np.errstate(invalid="ignore", divide="ignore"):
            rm = mod_median(np.where(ok, c / he, np.nan), grp)
            r = np.where(one, c / ha, np.nan)
        rn = r / mod_median(r, grp)[m]
        good = np.isfinite(rn)
        a1, pp, cur = harm(rn[good], x[good])
        per = np.array([harm(rn[g][np.isfinite(rn[g])], x[g][np.isfinite(rn[g])])[:2] for g in grp if np.isfinite(rn[g]).sum() >= PH_MIN])
        q1, q3 = np.nanpercentile(rm, [25, 75])
        rows.append([name, q, s, np.nanmedian(rm), q1, q3, a1, pp, np.median(per[:, 0]), np.median(per[:, 1]), len(per), int(good.sum())])
        curves += [[name, k, (k + 0.5) / NB, cur[k]] for k in range(NB)]
        for lo, hi in zip(edges[:-1], edges[1:]):
            sel = ~d["h_extra"] & (d["h_k"] > 0) & np.isfinite(c) & (c > 0) & (he >= lo) & (he < hi)
            ebins.append([name, lo / kE, hi / kE, int(sel.sum()), np.median(c[sel] / he[sel]) if sel.any() else np.nan])
    tsv(out / "estimator_table.tsv", ("estimator", "column", "scale_MeV_per_ADC", "ratio_to_h_e_median", "ratio_q1", "ratio_q3", "harm1_pooled",
                                     "ptp_pooled", "harm1_module_median", "ptp_module_median", "n_modules", "n_hits"), rows)
    tsv(out / "phase_curves.tsv", ("estimator", "bin", "phase", "median_ratio"), curves)
    tsv(out / "energy_bins.tsv", ("estimator", "e_lo_kE", "e_hi_kE", "n", "median_col_over_h_e"), ebins)
    # U2: the primary hit of a module and event with h_k >= 2 against its largest extra pulse
    ev = np.concatenate([z["h_ev"].astype(np.int64) + (i << 40) for i, z in enumerate(zs)])
    key = ev * 4096 + m
    first = np.flatnonzero(~d["h_extra"])
    kk = key[first]
    ex = np.flatnonzero(d["h_extra"])
    ex = ex[np.argsort(-he[ex], kind="stable")]
    _, u = np.unique(key[ex], return_index=True)
    ex = ex[u]                                                         # largest extra pulse per (event, module)
    pi = np.searchsorted(kk, key[ex])
    pr = first[pi]
    dt = d["h_t"][ex].astype(np.float64) - d["h_t"][pr]
    ok = (d["h_k"][pr] >= 2) & (ha[pr] > A_MIN)
    pr, dt = pr[ok], dt[ok]
    ref = {}
    for name in ("parpk", "parpk0", "rawpk"):
        c = d[COLS[name]].astype(np.float64)
        ref[name] = mod_median(np.where(one, c / ha, np.nan), grp)
    bins = [-60, -40, -30, -20, -14, -8, 0, 8, 14, 20, 30, 40, 60]
    rows = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        sel = (dt >= lo) & (dt < hi)
        i = pr[sel]
        row = [lo, hi, int(sel.sum())]
        sub = (d["h_apk_par"][i] / d["h_apk_par0"][i]).astype(np.float64)
        row += [np.nanmedian(sub) - 1 if sel.any() else np.nan, np.nanmean(np.abs(sub - 1) > 0.01) if sel.any() else np.nan]
        for name in ("rawpk", "parpk", "parpk0"):
            rr = d[COLS[name]][i] / ha[i] / ref[name][m[i]]
            row.append(np.nanmedian(rr) if sel.any() else np.nan)
        rows.append(row)
    tsv(out / "overlap.tsv", ("dt_lo_ns", "dt_hi_ns", "n", "par_over_par0_minus1_median", "frac_sub_gt_1pct", "rawpk_over_h_a_rel", "parpk_over_h_a_rel",
                              "parpk0_over_h_a_rel"), rows)
    # CFD against MF time
    dtc = d["h_tcfd"].astype(np.float64) - d["h_t"].astype(np.float64)
    okc = np.isfinite(dtc) & (d["h_k"] > 0)
    off = mod_median(np.where(okc, dtc, np.nan), grp)
    res = dtc - off[m]
    okc &= np.isfinite(res)
    rows = []
    for lab, sel in (("k1 a 50-200", (d["h_k"] == 1) & (ha > 50) & (ha <= 200)), ("k1 a 200-1000", (d["h_k"] == 1) & (ha > 200) & (ha <= 1000)),
                     ("k1 a >1000", (d["h_k"] == 1) & (ha > 1000)), ("k>=2 primary", (d["h_k"] >= 2) & ~d["h_extra"] & (ha > 200)),
                     ("extra pulses a>200", d["h_extra"] & (ha > 200))):
        v = res[sel & okc]
        qs = np.percentile(v, [1, 5, 50, 95, 99])
        rows.append([lab, int((sel & okc).sum()), 1.4826 * np.median(np.abs(v - np.median(v))), *qs, np.mean(np.abs(v) > 1.0), np.mean(np.abs(v) > 4.0),
                     1 - (sel & okc).sum() / max(sel.sum(), 1)])
    tsv(out / "cfd_vs_mf.tsv", ("class", "n", "sigma_mad_ns", "p01", "p05", "p50", "p95", "p99", "frac_gt_1ns", "frac_gt_4ns", "frac_no_cfd"), rows)
    if dumps:
        tsv(out / "phase_true.tsv", ("arm", "n_clusters", "std_E_over_Ep_permille", "cos1", "sin1", "amp1", "se1", "rms_fit_k12_permille"),
            [phase_true(*d_.split("=", 1)) for d_ in dumps])
        print(open(out / "phase_true.tsv").read())
    print(f"module median offset tcfd - t: median {np.nanmedian(off):.3f} ns, rms across modules {np.nanstd(off):.3f} ns, range "
          f"{np.nanmin(off):.2f} .. {np.nanmax(off):.2f}")
    print(open(out / "estimator_table.tsv").read())


if __name__ == "__main__":
    main(*sys.argv[1:])
