#!/usr/bin/env python3
"""p1b_hits.py <p1b scan npz> <out dir> [thr_clusters.npz]: AMP P1b hit-level diagnostics of the raw peak sample against the template amplitude.

Reference: h_a, the matched-filter amplitude.  For an efficient (matched) estimator the error of any other unbiased estimator splits into the
efficient one's error plus an orthogonal part, so var(col / h_a - 1) is the variance the column carries in excess of the matched filter (exact for
white noise and a unit-peak template, where cov(raw, h_a) = var(h_a)).  Hits: single pulse (h_k == 1, not extra), finite column, h_a >= A_LO.
Tables (tsv) in <out dir>:
  noise_modules.tsv   per module: pedestal sigma (median of the per-waveform clipped std, ADC), lag-1 correlation, benchmark sigma reference, template sum of squares
  noise_summary.tsv   medians over modules, the white-noise sigma_a of the matched filter (sigma / sqrt(h_sq)) and of the raw sample, per hit class
  scatter_<col>.tsv   robust sigma of col / h_a / module phase curve - 1 against h_a (bins), with and without the per-module phase curve removed
  fit.tsv             v(A) = (n / A)^2 + a / A + c^2 fitted to the binned variances: n (ADC), a (ADC), c, with the white-noise prediction n_pred
  classes.tsv         scatter at fixed amplitude for h_k == 1 against h_k >= 2 primaries and a pile-up neighbour split
  cluster_budget.tsv  (with thr_clusters.npz of the integral arm) propagation of v(A) to the elastic clusters, member by member, by component
"""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
import energy_resolution_benchmark as bm        # noqa: E402

A_LO, NB = 30.0, 16
EDGES = np.array([30, 40, 50, 65, 80, 100, 130, 170, 220, 300, 400, 550, 750, 1000, 1400, 2000, 4000.0])
COLS = {"rawpk": "h_apk_raw", "sum3": "h_s3", "sum4": "h_s4", "sum5": "h_s5", "sum7": "h_s7", "cfd5": "h_c5", "cfd7": "h_c7"}


def rsig(v):
    return 1.4826 * np.median(np.abs(v - np.median(v)))


def by_module(m, n_mod):
    o = np.argsort(m, kind="stable")
    b = np.r_[0, np.cumsum(np.bincount(m, minlength=n_mod))]
    return [o[b[i]:b[i + 1]] for i in range(n_mod)]


def residuals(r, mm, x, A, n_mod, hi=200.0):
    """r = col / h_a.  Returns r / module median (module = hits with A > hi), the same divided by the module's 16-bin phase curve (pooled curve
    for modules with < 2000 hits)."""
    b = np.minimum((x * NB).astype(int), NB - 1)
    grp = by_module(mm, n_mod)
    med = np.full(n_mod, np.nan)
    for k, g in enumerate(grp):
        s = g[A[g] > hi]
        if s.size > 200:
            med[k] = np.median(r[s])
    rn = r / med[mm]
    pooled = np.array([np.nanmedian(rn[(b == q) & (A > hi)]) for q in range(NB)])
    cur = np.tile(pooled, (n_mod, 1))
    for k, g in enumerate(grp):
        if g.size >= 2000:
            for q in range(NB):
                t = g[(b[g] == q) & (A[g] > hi)]
                if t.size > 20:
                    cur[k, q] = np.median(rn[t])
    return rn, rn / cur[mm, b]


def fit_v(Ac, v, w):
    """v = n^2 / A^2 + a / A + c^2 by weighted linear least squares on (1/A^2, 1/A, 1); non-negative terms by active-set."""
    X = np.stack([1 / Ac ** 2, 1 / Ac, np.ones_like(Ac)], axis=1)
    best = None
    for mask in ([1, 1, 1], [0, 1, 1], [1, 0, 1], [1, 1, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]):
        k = np.flatnonzero(mask)
        c, *_ = np.linalg.lstsq(X[:, k] * w[:, None], v * w, rcond=None)
        if (c >= 0).all():
            full = np.zeros(3)
            full[k] = c
            chi = np.sum(((X @ full - v) * w) ** 2)
            if best is None or chi < best[1]:
                best = (full, chi)
    return best[0]


def main(src, out, thr=None):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    z = np.load(src)
    names, n_mod = z["mod_name"], z["mod_name"].size
    t50 = bm.template_columns(bm.PYDIR / str(z["bank"]), names)[0]
    kE = float(z["beam_mev"]) / bm.E_REF
    have = "h_psd" in z.files                      # a plain PEAK_COLUMNS scan has only the raw column and no noise statistics
    cols = {k: v for k, v in COLS.items() if v in z.files}
    d = {q: z[q] for q in ("h_mod", "h_t", "h_a", "h_e", "h_k", "h_extra", *cols.values())}
    for q in ("h_psd", "h_c1", "h_sref", "h_sq", "h_st", "h_dp"):
        d[q] = z[q] if have else np.full(z["h_mod"].size, np.nan, np.float32)
    m = d["h_mod"].astype(int)
    ha, he = d["h_a"].astype(np.float64), d["h_e"].astype(np.float64)
    x = ((d["h_t"].astype(np.float64) - t50[m]) % bm.DT) / bm.DT
    fin = np.isfinite(ha) & (d["h_k"] > 0)
    one = fin & (d["h_k"] == 1) & ~d["h_extra"] & (ha >= A_LO)
    grp = by_module(m[fin], n_mod)
    idx = np.flatnonzero(fin)
    rows = []
    for k, g in enumerate(grp if have else []):
        s = idx[g]
        if s.size >= 200:
            rows.append([names[k], s.size, np.median(d["h_psd"][s]), np.median(d["h_c1"][s]), np.median(d["h_sref"][s]), np.nanmedian(d["h_sq"][s]) if np.isfinite(d["h_sq"][s]).any() else np.nan, np.nanmedian(d["h_st"][s]) if np.isfinite(d["h_st"][s]).any() else np.nan])
    if have:
        bm.write_tsv(out / "noise_modules.tsv", ("module", "n", "psd_adc", "lag1", "sig_ref_adc", "sq", "st"), rows)
    R = np.array([r[2:] for r in rows], dtype=float) if have else np.full((1, 5), np.nan)
    psd_med = np.median(R[:, 0])
    sq_med = np.nanmedian(R[:, 3])
    summ = [["modules", len(rows)], ["psd_adc_median", psd_med], ["psd_adc_q10", np.percentile(R[:, 0], 10)], ["psd_adc_q90", np.percentile(R[:, 0], 90)],
            ["lag1_median", np.median(R[:, 1])], ["sig_ref_median", np.median(R[:, 2])], ["template_sumsq_median", sq_med], ["template_sum_median", np.nanmedian(R[:, 4])],
            ["sigma_a_MF_adc (psd/sqrt(sq))", psd_med / np.sqrt(sq_med)], ["sigma_raw_adc (psd)", psd_med],
            ["sigma_excess_white_adc (psd sqrt(1-1/sq))", psd_med * np.sqrt(1 - 1 / sq_med)]]
    if have:
        bm.write_tsv(out / "noise_summary.tsv", ("quantity", "value"), summ)
    n_pred_mod = np.array([r[2] for r in rows]) * np.sqrt(1 - 1 / np.array([r[5] for r in rows]))
    n_pred_mod = dict(zip([r[0] for r in rows], n_pred_mod))
    npred_hit = np.array([n_pred_mod.get(str(names[k]), np.nan) for k in range(n_mod)])[m]
    fits, cres = [], {}
    A = ha[one]
    for name, q in cols.items():
        c = d[q].astype(np.float64)
        sel = one & np.isfinite(c)
        r = c[sel] / ha[sel]
        rn, res = residuals(r, m[sel], x[sel], ha[sel], n_mod)
        Asel = ha[sel]
        good = np.isfinite(res)
        rows = []
        for lo, hi in zip(EDGES[:-1], EDGES[1:]):
            s = good & (Asel >= lo) & (Asel < hi)
            if s.sum() < 500:
                continue
            v1, v0 = res[s] - 1, rn[s] / np.median(rn[s]) - 1
            v1, v0 = v1[np.abs(v1) < 0.5], v0[np.abs(v0) < 0.5]
            rows.append([lo, hi, int(s.sum()), np.median(Asel[s]), rsig(v1), rsig(v0), np.median(rn[s]) - 1 if True else 0,
                         np.median(npred_hit[sel][s] / Asel[s]) if name == "rawpk" else np.nan])
        bm.write_tsv(out / f"scatter_{name}.tsv", ("a_lo", "a_hi", "n", "a_median", "sigma_resid", "sigma_with_phase", "median_ratio_minus_1", "white_pred"), rows)
        T = np.array(rows, dtype=float)
        Ac, v = T[:, 3], T[:, 4] ** 2
        w = 1.0 / np.maximum(v, 1e-12)
        coef = fit_v(Ac, v, w * 0 + w)       # relative weights (equal relative accuracy per bin)
        cres[name] = coef
        npred = float(np.nanmedian(npred_hit[sel]))
        fits.append([name, int(good.sum()), np.sqrt(coef[0]), coef[1], np.sqrt(coef[2]), npred if name == "rawpk" else np.nan,
                     np.sqrt(np.sum(((coef[0] / Ac ** 2 + coef[1] / Ac + coef[2]) / v - 1) ** 2) / len(v))])
    bm.write_tsv(out / "fit.tsv", ("column", "n", "n_adc", "a_adc", "c", "n_pred_white_adc", "rel_fit_rms"), fits)
    # classes: fixed amplitude, raw only
    c = d["h_apk_raw"].astype(np.float64)
    allp = fin & ~d["h_extra"] & np.isfinite(c) & (ha >= A_LO)
    rn, res = residuals(np.where(allp, c / np.where(ha > 0, ha, 1), np.nan)[allp], m[allp], x[allp], ha[allp], n_mod)
    Aa, kk = ha[allp], d["h_k"][allp]
    rows = []
    for lo, hi in ((50, 100), (100, 200), (200, 500), (500, 1000), (1000, 2500)):
        for lab, s0 in (("h_k=1", kk == 1), ("h_k>=2 primary", kk >= 2)):
            s = s0 & (Aa >= lo) & (Aa < hi) & np.isfinite(res)
            if s.sum() > 200:
                v = res[s] - 1
                rows.append([lo, hi, lab, int(s.sum()), rsig(v), np.median(v), np.mean(np.abs(v) > 0.1)])
    bm.write_tsv(out / "classes.tsv", ("a_lo", "a_hi", "class", "n", "sigma_resid", "median", "frac_gt_10pct"), rows)
    if have:   # peak-sample position against the template peak (samples): is the scatter carried by hits whose largest sample is not the expected one?
        dp = d["h_dp"][allp].astype(np.float64)
        rows = []
        for lo, hi in ((50, 200), (200, 1000), (1000, 4000)):
            for k_lab, sel_dp in (("dp <= -1", dp <= -1), ("dp = -0.5 .. 0.5 (|dp| < 0.5)", np.abs(dp) < 0.5), ("0.5 <= dp < 1.5", (dp >= 0.5) & (dp < 1.5)), ("dp >= 1.5", dp >= 1.5)):
                s1 = sel_dp & (Aa >= lo) & (Aa < hi) & np.isfinite(res)
                if s1.sum() > 200:
                    v = res[s1] - 1
                    rows.append([lo, hi, k_lab, int(s1.sum()), float(s1.sum() / ((Aa >= lo) & (Aa < hi) & np.isfinite(res)).sum()), rsig(v), np.median(v)])
        bm.write_tsv(out / "peakpos.tsv", ("a_lo", "a_hi", "peak_sample_offset_from_template_peak", "n", "fraction", "sigma_resid", "median"), rows)
    np.save(out / "fit_coef.npy", np.array([cres[k] for k in cols]))
    print(open(out / "fit.tsv").read())
    if have:
        print(open(out / "noise_summary.tsv").read())


if __name__ == "__main__":
    main(*sys.argv[1:])
