#!/usr/bin/env python3
"""integral_vs_rawsum.py -- quick check: for single-pulse (k == 1) waveforms with a monotonic main peak, is the fitted integral A * I_m consistent with the raw sum?

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/integral_vs_rawsum.py [--cache <resid cache npz>] [--out <dir>] [--module-class PbWO4|LG|all]

Data: the S1 cache of resid_common.py (output/resid_study/resid_cache_N60000_P1500.npz; default fit, MONO_EDGE/NOISE_FLAG off; k == 1, no overflow sample,
fitted A 50-3000 ADC; raw samples Y, fitted pedestal, residual stored; first 1500 events (PbWO4) / 60000 (LG) of each EVIO split of the 13-run monoedge sample).
The scans store no raw sum, so the cache is the only source of raw samples: no new pass.  Overflow (nclip > 0) k == 1 channels are counted from the sample scans.
Monotonic main peak (offline, same definition as MF.protected / PILEUP_MONO_EDGE, no overflow): P = argmax of the raw samples, d_j = y[j+1] - y[j],
rise-monotonic: all d_j, j = P-3..P-1, > -2 sqrt2 sigma_ref;  fall-monotonic: all d_j, j = P..P+9 (cut at the waveform end), < +2 sqrt2 sigma_ref;  monotonic = both.
Fitted integral = A * I_m, I_m = "iint_wa" of the module's bank row (for lead glass and modules without a valid row the average template's integral, which the fit uses as their template): the phase-averaged
WaveAnalyzer walk-out sum (ADC x samples) of the noiseless unit template, triangle-smoothed (order 2), walk outward from the smoothed peak, a sample is "below" at < 0.1 x peak,
2 consecutive below samples end the walk.
Raw sum, def 1 (template window) = sum(Y - ped) over the samples that walk keeps on the fitted pulse (A x shape at the fitted onset, sampled, same smoothing and walk), contiguous around its peak.
Raw sum, def 2 (fixed window) = sum(Y - ped) over samples with t_i in [t0 - 20, t0 + 200] ns (t_i = 4 i, clipped to the waveform).
Ratio = A I_m / raw sum.  chi2/dof is the fit's.  Also kept: model window sum (the fitted pulse summed over the def-1 window) / (A I_m), the phase/window self-consistency of I_m.
-> integral_vs_rawsum_20261117.pdf / .tsv (table) / _waveforms.tsv (per waveform, capped 200k rows)
--module-class: PbWO4 (default, modules W*), LG (lead glass, modules G*) or all; every page, table and row is restricted to it, and the summary page states it.
"""
from __future__ import annotations

import sys
import json
import argparse
import textwrap
from pathlib import Path

HERE = Path(__file__).resolve().parent
PYDIR = HERE.parent.parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(PYDIR))
import numpy as np                                      # noqa: E402
import matplotlib                                       # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                         # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages    # noqa: E402
from mf_pulse import smooth_tri                         # noqa: E402

STAMP = "20261117"
CACHE = PYDIR / "output" / "resid_study" / "resid_cache_N60000_P1500.npz"
SAMPLE = PYDIR / "output" / "pileup_rates" / "monoedge_sample" / "old"
BANK_ROOT = PYDIR / "waveform_analysis" / "pulse_templates"
RTS = ["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17"]
DT, NS = 4.0, 100
MONO_NSIG, RISE_N, FALL_N = 2.0 * np.sqrt(2.0), 3, 10
CHI_EDGES = np.array([0, 1, 2, 5, 20, 100, np.inf])
CHI_LAB = ["<1", "1-2", "2-5", "5-20", "20-100", ">100"]
A_EDGES = np.array([50, 100, 200, 400, 800, 1600, 3000.0])
FIX_LO, FIX_HI = -20.0, 200.0
MIN_N = 20
CAP_ROWS = 200000
N_EX = 20
SEED = 20261117
COL = dict(zip(RTS, ["tab:blue", "tab:orange", "tab:green", "tab:red"]))
CLS_LAB = {"PbWO4": "PbWO4", "LG": "lead glass (LG)", "all": "PbWO4 + LG"}


def mono_test(Y, sig):
    """Monotonic main peak (no overflow rows): returns rise, fall [n] bool and the peak index."""
    n, ns = Y.shape
    y = Y.astype(np.float64)
    p = np.argmax(y, axis=1)
    d = y[:, 1:] - y[:, :-1]
    j = np.arange(ns - 1)[None, :]
    thr = (MONO_NSIG * sig.astype(np.float64))[:, None]
    rise = ~(((j >= (p - RISE_N)[:, None]) & (j <= (p - 1)[:, None])) & (d <= -thr)).any(axis=1)
    fall = ~(((j >= p[:, None]) & (j <= (p + FALL_N - 1)[:, None])) & (d >= thr)).any(axis=1)
    return rise, fall, p


def walk_window(v, ratio=0.1, nbreak=2):
    """Samples kept by the WaveAnalyzer outward walk on smoothed rows v [n, N] (mf_pulse.walk_integral): (first, last) sample index [n]."""
    n, N = v.shape
    ar = np.arange(n)
    pk = np.argmax(v, axis=1)
    cut = ratio * v[ar, pk]
    ext = []
    for side in (-1, 1):
        run = np.zeros(n, np.int64)
        act = np.ones(n, bool)
        e = np.zeros(n, np.int64)
        for dd in range(1, N):
            pos = pk + side * dd
            inr = act & (pos >= 0) & (pos < N)
            if not inr.any():
                break
            vv = v[ar, np.clip(pos, 0, N - 1)]
            below = (vv < cut) | (vv < 0)
            mb, ma = inr & below, inr & ~below
            run = np.where(mb, run + 1, np.where(ma, 0, run))
            act &= ~(mb & (run >= nbreak))
            e = np.where(ma, dd, e)
        ext.append(e)
    return pk - ext[0], pk + ext[1]


def load_iint(names, pwo4, rt):
    """I_m per module name for a run type: iint_wa of valid bank rows, median for the others (as the fit's fallback)."""
    z = np.load(BANK_ROOT / rt / f"templates_{rt}.npz", allow_pickle=False)
    d = {str(n): float(v) for n, v, ok in zip(z["names"], z["iint_wa"], z["valid"]) if ok > 0}
    inv = z["iint_wa"][z["valid"] == 0]                  # invalid rows carry the average template's integral (spline_template_extract: own.get(midx, avg_i)), the template the fit uses for them
    med = float(inv[0]) if inv.size else float(np.median(list(d.values())))
    return np.array([d.get(str(n), med) if p else med for n, p in zip(names, pwo4)]), len(d)


def derive(c):
    n = c["mi"].size
    Y = c["Y"].astype(np.float64)
    ped = c["ped"].astype(np.float64)
    y = Y - ped[:, None]
    model = y - c["resid"].astype(np.float64)
    rise, fall, p = mono_test(Y, c["sig"])
    A = c["a"].astype(np.float64)
    im = np.zeros(n)
    for rt in RTS:
        s = c["rt"] == rt
        tab, _ = load_iint(c["names"], c["pwo4"], rt)
        im[s] = tab[c["mi"][s]]
    fit = A * im
    f, l = walk_window(smooth_tri(model, 2))
    idx = np.arange(NS)[None, :]
    wm = (idx >= f[:, None]) & (idx <= l[:, None])
    t = DT * idx
    t0 = c["t0"].astype(np.float64)[:, None]
    fm = (t >= t0 + FIX_LO) & (t <= t0 + FIX_HI)
    r = dict(A=A, im=im, fit=fit, rise=rise, fall=fall, mono=rise & fall, f=f, l=l,
             raw1=(y * wm).sum(axis=1), raw2=(y * fm).sum(axis=1), mod1=(model * wm).sum(axis=1), chi2=c["chi2"].astype(np.float64))
    with np.errstate(divide="ignore", invalid="ignore"):
        r["rat1"], r["rat2"], r["modrat"] = fit / r["raw1"], fit / r["raw2"], r["mod1"] / fit
    r["cb"] = np.digitize(r["chi2"], CHI_EDGES[1:-1])
    r["cb"][~np.isfinite(r["chi2"])] = -1
    r["ab"] = np.clip(np.digitize(A, A_EDGES[1:-1]), 0, len(A_EDGES) - 2)
    return r


def overflow_counts(cls):
    out = {}
    for fn in sorted(SAMPLE.glob("run_*.npz")):
        z = np.load(fn)
        rt = str(z["run_type"])
        k1 = (z["ch_k"] == 1) & class_mask(z["mod_pwo4"][z["ch_midx"]], cls)
        o = out.setdefault(rt, [0, 0])
        o[0] += int((k1 & (z["ch_nclip"] > 0)).sum())
        o[1] += int(k1.sum())
    return out


def class_mask(pw, cls):
    """Rows of module class cls: PbWO4 = W* modules, LG = G* lead glass, all = every row (pw = PbWO4 flag per row)."""
    return np.ones(pw.shape, bool) if cls == "all" else (pw if cls == "PbWO4" else ~pw)


def q(x, p):
    return np.percentile(x, p) if x.size else np.nan


def stats(x):
    if x.size == 0:
        return dict(n=0, med=np.nan, rms=np.nan, lo=np.nan, hi=np.nan, f2=np.nan)
    return dict(n=int(x.size), med=float(np.median(x)), rms=float(x.std()), lo=q(x, 16), hi=q(x, 84), f2=float((np.abs(x - 1) <= 0.02).mean()))


def make_rows(d, r, ok, okall):
    rows = []
    for cls, mm in (("mono", r["mono"]), ("nonmono", ~r["mono"]), ("trunc_mono", r["mono"])):
        okc = okall & ~r["inside"] if cls == "trunc_mono" else ok
        for rt in RTS + ["all"]:
            srt = (d["rt"] == rt) if rt != "all" else np.ones(d["rt"].size, bool)
            for cb in range(len(CHI_LAB)):
                s = okc & mm & srt & (r["cb"] == cb)
                row = dict(sel=cls, run_type=rt, chi2_bin=CHI_LAB[cb], n=int(s.sum()))
                for k, tag in (("rat1", "tmpl"), ("rat2", "fixed")):
                    st = stats(r[k][s])
                    row.update({f"med_{tag}": st["med"], f"rms_{tag}": st["rms"], f"p16_{tag}": st["lo"], f"p84_{tag}": st["hi"], f"frac2pct_{tag}": st["f2"]})
                row["med_modrat"] = float(np.median(r["modrat"][s])) if s.any() else np.nan
                rows.append(row)
    return rows


def write_tsv(path, rows):
    keys = list(rows[0])
    with open(path, "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[k]:.5g}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")


def text_page(pdf, title, lines, size=9):
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.04, 0.96, title, fontsize=12, weight="bold", va="top")
    fig.text(0.04, 0.92, "\n".join(textwrap.fill(l, 150) for l in lines), fontsize=size, va="top", family="monospace" if size < 9 else None)
    pdf.savefig(fig)
    plt.close(fig)


def page_scatter(pdf, d, r, ok, key, lab, rng, mcls):
    fig, axs = plt.subplots(2, 2, figsize=(11, 8.5), constrained_layout=True)
    for ax, rt in zip(axs.flat, RTS):
        s = ok & r["mono"] & (d["rt"] == rt)
        x, y = r[key][s], r["fit"][s]
        g = (x > 0) & (y > 0)
        lo, hi = 200.0, 4e4
        _, _, _, im = ax.hist2d(x[g], y[g], bins=[np.geomspace(lo, hi, 151)] * 2, cmin=1, norm=matplotlib.colors.LogNorm(vmin=1), cmap="viridis")   # 150 x 150 bins, counts >= 1
        fig.colorbar(im, ax=ax, label="waveforms per bin")
        hi20 = g & (r["chi2"][s] > 20)
        sel = np.nonzero(hi20)[0]
        if sel.size > 1500:
            sel = rng.choice(sel, 1500, replace=False)
        ax.plot(x[sel], y[sel], "o", mfc="tab:red", mec="white", mew=0.4, ms=4, label=f"chi2/dof > 20 (n = {int(hi20.sum())}, {sel.size} shown)")
        ax.plot([lo, hi], [lo, hi], "b-", lw=0.8, label="y = x")
        ax.set_xscale("log"), ax.set_yscale("log")
        ax.set_xlabel(f"raw sum, {lab} (ADC x samples)"), ax.set_ylabel("fitted integral A I_m")
        ax.set_title(f"{rt}: monotonic k = 1, n = {int(s.sum())}", fontsize=10)
        ax.legend(fontsize=7, loc="upper left")
    fig.suptitle(f"(1) Fitted integral vs raw sum ({lab}), monotonic single peaks, module class {CLS_LAB[mcls]}", fontsize=11)
    pdf.savefig(fig)
    plt.close(fig)


def profile(ax, xb, x, y, col, ls, label, nmin=MIN_N):
    xs, m, lo, hi = [], [], [], []
    for b in range(len(xb)):
        s = x == b
        if s.sum() >= nmin:
            xs.append(b), m.append(np.median(y[s])), lo.append(q(y[s], 16)), hi.append(q(y[s], 84))
    xs = np.array(xs)
    if xs.size:
        ax.plot(xs, m, ls, color=col, lw=1.3, marker="o", ms=3, label=label)
        ax.fill_between(xs, lo, hi, color=col, alpha=0.12 if ls == "-" else 0.06)


def page_profile(pdf, d, r, ok, xkey, xlab, labs, n_b, tag, mcls):
    fig, axs = plt.subplots(2, 2, figsize=(11, 8.5), constrained_layout=True, sharey=True)
    pw = d["pwo4"][d["mi"]]
    for row, (rk, rl) in enumerate((("rat1", "def 1 template window"), ("rat2", "def 2 fixed window"))):
        for colx, (cls, mm) in enumerate((("monotonic k = 1", r["mono"]), ("non-monotonic k = 1 (comparison)", ~r["mono"]))):
            ax = axs[row, colx]
            for rt in RTS:
                for cl, ls in ((True, "-"), (False, "--")):
                    s = ok & mm & (d["rt"] == rt) & (pw == cl) & np.isfinite(r[rk])
                    profile(ax, range(n_b), r[xkey][s], r[rk][s], COL[rt], ls, f"{rt} {'PbWO4' if cl else 'LG'}")
            ax.axhline(1, color="k", lw=0.6)
            ax.set_xticks(range(n_b)), ax.set_xticklabels(labs, fontsize=8)
            ax.set_xlabel(xlab), ax.set_ylabel("A I_m / raw sum")
            ax.set_title(f"{cls}, {rl}", fontsize=9)
            ax.set_ylim(0.8, 1.3)
            if row == 0 and colx == 0:
                ax.legend(fontsize=6, ncol=2)
    fig.suptitle(f"({tag}) Ratio fitted / raw, module class {CLS_LAB[mcls]}: median and 16-84 % band (bins with >= {MIN_N} waveforms{'; dashed = lead glass' if mcls == 'all' else ''})", fontsize=11)
    pdf.savefig(fig)
    plt.close(fig)


def page_hist(pdf, d, r, ok):
    fig, axs = plt.subplots(2, 4, figsize=(14, 7), constrained_layout=True)
    bins = np.linspace(0.8, 1.3, 101)
    for row, (rk, rl) in enumerate((("rat1", "def 1"), ("rat2", "def 2"))):
        for ax, rt in zip(axs[row], RTS):
            for (lo, hi), col, lab in (((0, 5), "tab:blue", "chi2/dof < 5"), ((20, np.inf), "tab:red", "chi2/dof > 20")):
                s = ok & r["mono"] & (d["rt"] == rt) & (r["chi2"] >= lo) & (r["chi2"] < hi) & np.isfinite(r[rk])
                if s.sum():
                    ax.hist(np.clip(r[rk][s], 0.8, 1.3), bins=bins, density=True, histtype="step", color=col, lw=1.3, label=f"{lab} (n {int(s.sum())}, med {np.median(r[rk][s]):.3f})")
            ax.axvline(1, color="k", lw=0.6)
            ax.set_title(f"{rt} {rl}", fontsize=9), ax.set_xlabel("A I_m / raw sum (clipped to 0.8-1.3)"), ax.legend(fontsize=6)
    fig.suptitle("(4) Ratio fitted / raw, monotonic k = 1: chi2/dof < 5 vs > 20 (area normalised)", fontsize=11)
    pdf.savefig(fig)
    plt.close(fig)


def page_table(pdf, rows, rt_list):
    for sel in ("mono", "nonmono", "trunc_mono"):
        fig = plt.figure(figsize=(11, 8.5))
        fig.text(0.03, 0.97, f"(5) Table, {dict(mono='monotonic', nonmono='NON-monotonic (comparison)', trunc_mono='monotonic but TRUNCATED fit window (late pulse, A I_m extrapolates past the waveform end)')[sel]} k = 1: n, median, rms of the ratio A I_m / raw sum, fraction within +-2 %  (def 1 template window | def 2 fixed window)", fontsize=9.5, va="top")
        L = [f"{'run type':16s}{'chi2/dof':>8s}{'n':>8s} | {'med1':>7s}{'rms1':>7s}{'f2%1':>6s} | {'med2':>7s}{'rms2':>7s}{'f2%2':>6s}"]
        for rt in rt_list:
            for r in rows:
                if r["sel"] == sel and r["run_type"] == rt:
                    f = lambda v, w, p=3: f"{v:{w}.{p}f}" if np.isfinite(v) else " " * (w - 2) + "NA"   # noqa: E731
                    L.append(f"{rt:16s}{r['chi2_bin']:>8s}{r['n']:8d} | {f(r['med_tmpl'], 7)}{f(r['rms_tmpl'], 7)}{f(100 * r['frac2pct_tmpl'], 6, 1)} | {f(r['med_fixed'], 7)}{f(r['rms_fixed'], 7)}{f(100 * r['frac2pct_fixed'], 6, 1)}")
            L.append("")
        fig.text(0.03, 0.93, "\n".join(L), fontsize=7, family="monospace", va="top")
        pdf.savefig(fig)
        plt.close(fig)


def page_examples(pdf, c, d, r, ok, rng, trunc=False):
    cand = np.nonzero(ok & r["mono"] & (r["chi2"] > 20))[0]
    pick = []
    per = N_EX // len(RTS)
    for rt in RTS:
        cr = cand[d["rt"][cand] == rt]
        pick += list(rng.choice(cr, min(per, cr.size), replace=False)) if cr.size else []
    rest = np.setdiff1d(cand, pick)
    if len(pick) < N_EX and rest.size:                    # run types with fewer than 5: fill with random others
        pick += list(rng.choice(rest, min(N_EX - len(pick), rest.size), replace=False))
    fig, axs = plt.subplots(5, 4, figsize=(15, 13), constrained_layout=True)
    tt = DT * np.arange(NS)
    for ax, i in zip(axs.flat, pick):
        ped = float(c["ped"][i])
        y = c["Y"][i]
        ax.plot(tt, y, "k.", ms=3.5)
        ax.axhline(ped, color="0.6", lw=0.7)
        ax.plot(tt, y - c["resid"][i], "r-", lw=1.0)
        ax.axvspan(DT * r["f"][i], DT * r["l"][i], color="tab:blue", alpha=0.12)
        ax.axvline(c["t0"][i] + FIX_LO, color="tab:green", ls="--", lw=0.8)
        ax.axvline(min(c["t0"][i] + FIX_HI, tt[-1]), color="tab:green", ls="--", lw=0.8)
        ax.set_title(f"{d['rt'][i]} {c['names'][c['mi'][i]]} r{int(c['run'][i])} e{int(c['ev'][i])}  A {r['A'][i]:.0f}  chi2 {r['chi2'][i]:.0f}\n"
                     f"fit int {r['fit'][i]:.0f}, raw1 {r['raw1'][i]:.0f} ({r['rat1'][i]:.3f}), raw2 {r['raw2'][i]:.0f} ({r['rat2'][i]:.3f})", fontsize=7, loc="left", pad=2)
        ax.tick_params(labelsize=6)
    fig.suptitle(f"(6{'b, TRUNCATED fit window' if trunc else ''}) {len(pick)} monotonic k = 1 waveforms with chi2/dof > 20 ({'fit window beyond the waveform end; ' if trunc else ''}up to {per} per run type, rest random, seed {SEED}).  Black data, red fit (Y - resid), grey pedestal; "
                 "blue band = def 1 template window, green dashed = def 2 fixed window [t0 - 20, t0 + 200] ns", fontsize=9)
    pdf.savefig(fig)
    plt.close(fig)
    return pick


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cache", default=str(CACHE))
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--module-class", choices=["PbWO4", "LG", "all"], default="PbWO4", help="restrict every plot, table and example to this module class (PbWO4 = W* modules, LG = G* lead glass, all = both)")
    a = ap.parse_args()
    c = dict(np.load(a.cache))
    print(f"[cache] {a.cache}: {c['mi'].size} S1 waveforms", flush=True)
    cm = class_mask(c["pwo4"][c["mi"]], a.module_class)                 # per waveform: module class (pwo4 flag of its module)
    # S1 definition of resid_common.load_s1: the 32-sample fit window (from ceil(t0/4) - 3) must lie inside the waveform.  The truncated ones (late pulses whose A I_m extrapolates beyond the waveform end) stay out of pages 1-4 and 6 and form their own class in the table (trunc_mono).
    w0 = np.ceil(c["t0"] / DT).astype(np.int64) - 3
    inside = (w0 >= 0) & (w0 + 32 <= NS)
    trunc = {}
    for rt in RTS:
        s = (c["rt"] == rt) & cm
        trunc[rt] = (int((s & ~inside).sum()), int((s & ~inside & (c["chi2"] > 20)).sum()), int((s & inside & (c["chi2"] > 20)).sum()))
    r = derive(c)
    d = dict(rt=c["rt"], mi=c["mi"], pwo4=c["pwo4"])
    okall = np.isfinite(r["chi2"]) & np.isfinite(r["rat1"]) & np.isfinite(r["rat2"]) & (r["raw1"] > 0) & (r["raw2"] > 0)
    n_bad = int((~okall & cm).sum())
    okall &= cm
    r["inside"] = inside
    ok = okall & inside
    ovf = overflow_counts(a.module_class)
    rng = np.random.default_rng(SEED)
    rows = make_rows(d, r, ok, okall)
    out = Path(a.out)
    write_tsv(out / f"integral_vs_rawsum_{STAMP}.tsv", rows)
    # per-waveform TSV (capped)
    idx = np.nonzero(cm)[0]
    if idx.size > CAP_ROWS:
        idx = np.sort(rng.choice(idx, CAP_ROWS, replace=False))
    keys = ["run_type", "run", "event", "module", "class", "A", "t0", "chi2dof", "mono", "inside", "rise_mono", "fall_mono", "I_m", "fit_int", "win1_lo_ns", "win1_hi_ns", "raw_tmpl", "raw_fixed", "ratio_tmpl", "ratio_fixed", "model_over_fit_tmpl"]
    with open(out / f"integral_vs_rawsum_{STAMP}_waveforms.tsv", "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for i in idx:
            fh.write("\t".join(str(v) for v in (c["rt"][i], c["run"][i], c["ev"][i], c["names"][c["mi"][i]], "PbWO4" if c["pwo4"][c["mi"][i]] else "LG", f"{r['A'][i]:.2f}", f"{c['t0'][i]:.2f}",
                                                f"{r['chi2'][i]:.4g}", int(r["mono"][i]), int(r["inside"][i]), int(r["rise"][i]), int(r["fall"][i]), f"{r['im'][i]:.4f}", f"{r['fit'][i]:.2f}", DT * r["f"][i], DT * r["l"][i],
                                                f"{r['raw1'][i]:.2f}", f"{r['raw2'][i]:.2f}", f"{r['rat1'][i]:.5f}", f"{r['rat2'][i]:.5f}", f"{r['modrat'][i]:.5f}")) + "\n")
    # summary lines
    L = [f"Module class: {CLS_LAB[a.module_class]} (--module-class {a.module_class}; PbWO4 = modules W*, lead glass = modules G*). Every plot, table, example and TSV row is restricted to this class.",
         f"Source: S1 cache {Path(a.cache).name}: default fit (all switches off), k == 1, no overflow sample, fitted A 50-3000 ADC, fit window inside the waveform (as load_s1), {int(cm.sum())} waveforms of this class ({n_bad} dropped: nan chi2 or raw sum <= 0).",
         "Removed as truncated (fit window beyond the waveform end; n, of which chi2/dof > 20, vs chi2/dof > 20 kept): " + "; ".join(f"{rt} {v[0]}, {v[1]} vs {v[2]}" for rt, v in trunc.items()),
         "Waveforms per run type (monotonic / non-monotonic, class split PbWO4 / LG):"]
    for rt in RTS:
        s = ok & (d["rt"] == rt)
        pw = d["pwo4"][d["mi"]]
        L.append(f"  {rt:16s} total {int(s.sum()):6d}  mono {int((s & r['mono']).sum()):6d} (PbWO4 {int((s & r['mono'] & pw).sum())}, LG {int((s & r['mono'] & ~pw).sum())})  non-mono {int((s & ~r['mono']).sum()):6d}"
                 f"   |  overflow k == 1 channels in the 2-file scans (excluded from the cache): {ovf.get(rt, [0, 0])[0]} of {ovf.get(rt, [0, 0])[1]} k == 1 channels")
    L += ["", "Definitions (see script docstring): fitted integral = A x I_m (iint_wa, phase-averaged WaveAnalyzer walk-out sum of the smoothed unit template, ADC x samples; 0.1 x peak, 2 below samples end the walk).",
          "raw sum def 1 = sum(Y - ped) over the samples the same walk keeps on the fitted pulse (A x shape at the fitted onset, sampled); def 2 = sum(Y - ped) for t in [t0 - 20, t0 + 200] ns.",
          "Monotonic = rise (3 differences up to the peak > -2 sqrt2 sigma_ref) and fall (10 differences after the peak < +2 sqrt2 sigma_ref); PILEUP_MONO_EDGE's test, applied offline to the raw samples.",
          "Limits: A in 50-3000 ADC; PbWO4 only from the first 1500 events of each split (LG: 60000), so PbWO4 is a sample of the sample; cache run set = 13 runs x 2 splits.", "",
          "Headline (median ratio, rms; template window | fixed window), monotonic:"]
    for rt in RTS:
        for lab, lo, hi in (("chi2 < 5", 0, 5), ("chi2 > 20", 20, np.inf)):
            s = ok & r["mono"] & (d["rt"] == rt) & (r["chi2"] >= lo) & (r["chi2"] < hi)
            a1, a2 = stats(r["rat1"][s]), stats(r["rat2"][s])
            L.append(f"  {rt:16s} {lab:10s} n {a1['n']:6d}  def 1 med {a1['med']:.4f} rms {a1['rms']:.4f} [16-84: {a1['lo']:.3f}-{a1['hi']:.3f}] in +-2 % {100 * a1['f2']:.1f} %"
                     f"  |  def 2 med {a2['med']:.4f} rms {a2['rms']:.4f} in +-2 % {100 * a2['f2']:.1f} %")
    L.append("Monotonic, TRUNCATED fit window (not in the pages above):")
    for rt in RTS:
        for lab, lo, hi in (("chi2 < 5", 0, 5), ("chi2 > 20", 20, np.inf)):
            s = okall & ~inside & r["mono"] & (d["rt"] == rt) & (r["chi2"] >= lo) & (r["chi2"] < hi)
            a1 = stats(r["rat1"][s])
            L.append(f"  {rt:16s} {lab:10s} n {a1['n']:6d}  def 1 med {a1['med']:.4f} rms {a1['rms']:.4f} [16-84: {a1['lo']:.3f}-{a1['hi']:.3f}]")
    with PdfPages(out / f"integral_vs_rawsum_{STAMP}.pdf") as pdf:
        text_page(pdf, "Fitted integral A I_m vs raw sum for monotonic single peaks (quick check, 2026-11-17)", L, 8)
        page_scatter(pdf, d, r, ok, "raw1", "def 1 template window", rng, a.module_class)
        page_scatter(pdf, d, r, ok, "raw2", "def 2 fixed window", rng, a.module_class)
        page_profile(pdf, d, r, ok, "cb", "chi2/dof bin", CHI_LAB, len(CHI_LAB), "2", a.module_class)
        page_profile(pdf, d, r, ok, "ab", "fitted amplitude A (ADC)", [f"{int(lo)}-{int(hi)}" for lo, hi in zip(A_EDGES[:-1], A_EDGES[1:])], len(A_EDGES) - 1, "3", a.module_class)
        page_hist(pdf, d, r, ok)
        page_table(pdf, rows, RTS + ["all"])
        pick = page_examples(pdf, c, d, r, ok, rng)
        page_examples(pdf, c, d, r, okall & ~r["inside"], rng, trunc=True)
    # amplitude drift summary + example summary, printed
    print("\n".join(L))
    print(f"\nAmplitude drift (monotonic, module class {CLS_LAB[a.module_class]}, def 1 median ratio per A bin, chi2<5 | chi2>20):")
    for rt in RTS:
        parts = []
        for b in range(len(A_EDGES) - 1):
            v = []
            for lo, hi in ((0, 5), (20, np.inf)):
                s = ok & r["mono"] & (d["rt"] == rt) & (r["ab"] == b) & (r["chi2"] >= lo) & (r["chi2"] < hi)
                v.append(f"{np.median(r['rat1'][s]):.3f}(n{int(s.sum())})" if s.sum() >= MIN_N else "NA")
            parts.append(f"A{int(A_EDGES[b])}-{int(A_EDGES[b + 1])}: {v[0]} | {v[1]}")
        print(f"  {rt}: " + "; ".join(parts))
    print("\nExamples:")
    for i in pick:
        print(f"  {c['rt'][i]} {c['names'][c['mi'][i]]} run {int(c['run'][i])} ev {int(c['ev'][i])} A {r['A'][i]:.0f} chi2 {r['chi2'][i]:.0f} ratio1 {r['rat1'][i]:.3f} ratio2 {r['rat2'][i]:.3f} modrat {r['modrat'][i]:.3f}")
    print(f"[model window sum / fit integral] median over ok: {np.median(r['modrat'][ok]):.4f}")


if __name__ == "__main__":
    main()
