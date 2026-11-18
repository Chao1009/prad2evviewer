#!/usr/bin/env python3
"""integral_vs_rawsum_full.py -- fitted integral A I_m vs raw sums for monotonic single peaks, from a PILEUP_RAWSUM=1 scan directory (PbWO4 modules W* only).

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/integral_vs_rawsum_full.py <scan dir> [--stamp <tag>] [--out <dir>] [--nproc N]
    e.g. output/pileup_rates/rawsum_sample  ->  integral_vs_rawsum_full_20261117.pdf / .tsv / _counts.tsv
         output/pileup_rates/rawsum_full    --stamp all_<date>  ->  integral_vs_rawsum_full_all_<date>.pdf

Data: run_*.npz of pileup_rates.py scan with PILEUP_RAWSUM=1 (ch_rs_raw1, ch_rs_raw2, ch_rs_rise, ch_rs_fall, ch_rs_phase per kept channel; main-pulse quantities).
Class: PbWO4 channels (mod_pwo4).  Selection ("used"): k == 1, no overflow sample (nclip == 0), monotonic (rise and fall flags = MF.protected / PILEUP_MONO_EDGE test),
main-pulse fit window inside the waveform (not trunc_main), fitted A >= A_MIN ADC, finite chi2/dof, both raw sums > 0.
Fitted integral = A * I_m, I_m = iint_wa of the module's bank row (integral_vs_rawsum.load_iint).  Raw sum def 1 = sum(Y - ped) over the template walk-out samples placed at the fitted
onset, def 2 = sum(Y - ped) over [t0 - 20, t0 + 200] ns.  Ratio = A I_m / raw sum.
All-cases pages (all_page, after fit_page): any k, monotonic or not, overflow and truncated kept, A_main >= A_MIN with A_main the largest pulse (as rawsum_vars);
x1: A_main I_m vs raw sum def 1; x2: sum over the fitted pulses of A_i I_m vs raw sum def 2 (per-pulse amplitudes ch_a, columns < k).
Every file is reduced to histograms (ratio bins 0.002 on 0.5-1.5, clipped into the end bins), so the full scan needs no per-waveform memory; medians and 16-84 % are
read off the histograms (resolution 0.002), the rms and the +-2 % fraction are exact.  Examples: the scan's chi2 > 10 refit reservoirs (qx_chi2), k == 1, chi2/dof > 20, monotonic, refitted shape rebuilt from the bank.
"""
from __future__ import annotations

import os
import sys
import argparse
import textwrap
from pathlib import Path
from multiprocessing import Pool

HERE = Path(__file__).resolve().parent
PYDIR = HERE.parent.parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(PYDIR))
sys.path.insert(0, str(HERE))
import numpy as np                                      # noqa: E402
import matplotlib                                       # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                         # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages    # noqa: E402
from matplotlib.colors import LogNorm                   # noqa: E402
import pileup_rates as pr                               # noqa: E402  (its import sets a plot style; reset below)
from integral_vs_rawsum import load_iint                # noqa: E402

plt.rcdefaults()
RTS = pr.RUN_TYPES
A_MIN = 50.0
PNG_DIR = os.environ.get("RAWSUM_PNG_DIR", "")
SEED = 20261117
N_EX = 20
RB = np.linspace(0.5, 1.5, 501)                         # ratio bins (0.002)
NR = RB.size - 1
RAW2D = np.geomspace(100.0, 1e5, 161)                   # (a): raw sum and fitted integral
RAWB = np.geomspace(100.0, 1e5, 41)                     # (b), (f): raw sum bins
PHB = np.linspace(0.0, 4.0, 41)                         # (c): onset phase (ns)
CHB = np.geomspace(0.1, 1000.0, 61)                     # (d): chi2/dof
CHI_EDGES = np.array([0, 1, 2, 5, 20, 100, np.inf])     # table bins
CHI_LAB = ["<1", "1-2", "2-5", "5-20", "20-100", ">100"]
MIN_MOD = 50                                            # (e): modules with >= this many used waveforms
MIN_BIN = 50                                            # median / band drawn for bins with >= this many waveforms
COL = dict(zip(RTS, ["tab:blue", "tab:orange", "tab:green", "tab:red"]))


def idx(x, e):
    return np.clip(np.searchsorted(e, x, side="right") - 1, 0, e.size - 2)


def h2(x, y, xe, ye):
    nx, ny = xe.size - 1, ye.size - 1
    return np.bincount(idx(x, xe) * ny + idx(y, ye), minlength=nx * ny).reshape(nx, ny).astype(np.float64)


def quant(H, p):
    """Quantile p of each row of the ratio histogram H [nx, NR] (linear within the bin); nan where the row is empty."""
    c = np.cumsum(H, axis=1)
    n = c[:, -1]
    tgt = p * n
    i = np.minimum((c < tgt[:, None]).sum(axis=1), NR - 1)
    ar = np.arange(H.shape[0])
    prev = np.where(i > 0, c[ar, np.maximum(i - 1, 0)], 0.0)
    f = (tgt - prev) / np.maximum(H[ar, i], 1e-12)
    out = RB[i] + np.clip(f, 0, 1) * (RB[i + 1] - RB[i])
    return np.where(n > 0, out, np.nan)


def new_acc(nmod):
    z = lambda *s: np.zeros(s)                          # noqa: E731
    d = dict(n=np.zeros(8), xn=np.zeros(2))                             # chans, k1, k1 no overflow, mono, trunc, A < A_MIN, bad value, used
    for t in ("1", "2"):
        d.update({f"a{t}": z(RAW2D.size - 1, RAW2D.size - 1), f"b{t}": z(RAWB.size - 1, NR), f"c{t}": z(PHB.size - 1, NR), f"d{t}": z(CHB.size - 1, NR), f"m{t}": z(nmod, NR),
                  f"tab{t}": z(len(CHI_LAB), NR), f"s{t}": z(len(CHI_LAB), 3), f"s{t}_f2": z(len(CHI_LAB)), f"mod{t}": z(nmod, 3), f"x{t}": z(RAW2D.size - 1, RAW2D.size - 1), f"xr{t}": z(NR)})   # s: n, sum, sum2 ; per module: n, sum, sum2
    return d


def add_acc(a, b):
    for k in a:
        a[k] += b[k]


def reduce_file(path):
    z = np.load(path)
    rt = str(z["run_type"])
    names, pw = z["mod_name"], z["mod_pwo4"]
    im = load_iint(names, pw, rt)[0]
    mi = z["ch_midx"].astype(np.int64)
    sel = pw[mi]
    acc = new_acc(names.size)
    acc["n"][0] = sel.sum()
    k1 = sel & (z["ch_k"] == 1)
    k1o = k1 & (z["ch_nclip"] == 0)
    mono = k1o & z["ch_rs_rise"] & z["ch_rs_fall"]
    tr = mono & z["ch_trunc_main"]
    a = z["ch_a"][:, 0].astype(np.float64)
    low = mono & ~tr & (a < A_MIN)
    raw1, raw2 = z["ch_rs_raw1"].astype(np.float64), z["ch_rs_raw2"].astype(np.float64)
    chi = z["ch_chi2"].astype(np.float64)
    ok = np.isfinite(chi) & (raw1 > 0) & (raw2 > 0)
    used = mono & ~tr & (a >= A_MIN) & ok
    acc["n"][1:7] = [k1.sum(), k1o.sum(), mono.sum(), tr.sum(), low.sum(), (mono & ~tr & (a >= A_MIN) & ~ok).sum()]
    acc["n"][7] = used.sum()
    u = np.nonzero(used)[0]
    fit = a[u] * im[mi[u]]
    ph, ch, m = z["ch_rs_phase"][u].astype(np.float64), chi[u], mi[u]
    cb = np.clip(np.digitize(ch, CHI_EDGES[1:-1]), 0, len(CHI_LAB) - 1)
    for t, raw in (("1", raw1[u]), ("2", raw2[u])):
        r = fit / raw
        rbin = idx(r, RB)
        acc[f"a{t}"] += h2(raw, fit, RAW2D, RAW2D)
        acc[f"b{t}"] += h2(raw, r, RAWB, RB)
        acc[f"c{t}"] += h2(ph, r, PHB, RB)
        acc[f"d{t}"] += h2(ch, r, CHB, RB)
        acc[f"m{t}"] += np.bincount(m * NR + rbin, minlength=names.size * NR).reshape(names.size, NR)
        acc[f"tab{t}"] += np.bincount(cb * NR + rbin, minlength=len(CHI_LAB) * NR).reshape(len(CHI_LAB), NR)
        w2 = (np.abs(r - 1) <= 0.02).astype(np.float64)
        for j, v in enumerate((np.ones_like(r), r, r * r)):
            acc[f"s{t}"][:, j] += np.bincount(cb, weights=v, minlength=len(CHI_LAB))
            acc[f"mod{t}"][:, j] += np.bincount(m, weights=v, minlength=names.size)
        acc[f"s{t}_f2"] += np.bincount(cb, weights=w2, minlength=len(CHI_LAB))
    # all cases (pages x1, x2): any k, monotonic or not, overflow and truncated kept; A_main = largest pulse (as rawsum_vars); x2 total = sum over the fitted pulses
    va = np.where(np.arange(z["ch_a"].shape[1]) < z["ch_k"][:, None], z["ch_a"].astype(np.float64), 0.0)
    ua = np.nonzero(sel & (va.max(axis=1) >= A_MIN) & ok)[0]
    acc["xn"] += [ua.size, (z["ch_k"][ua] >= 2).sum()]
    for t, raw, fit in (("1", raw1[ua], va[ua].max(axis=1) * im[mi[ua]]), ("2", raw2[ua], va[ua].sum(axis=1) * im[mi[ua]])):
        acc[f"x{t}"] += h2(raw, fit, RAW2D, RAW2D)
        acc[f"xr{t}"] += np.bincount(idx(fit / raw, RB), minlength=NR)
    return rt, acc, str(z["run"]), names, pw, z["n_seen"], z["n_pass"]


def load_examples(scan_dir, rng):
    """chi2 > 20 monotonic k == 1 PbWO4 examples from the qx_chi2 reservoirs: dicts with samples, ped, model, windows."""
    pool = {rt: [] for rt in RTS}
    for fn in sorted(scan_dir.glob("run_*.npz")):
        z = np.load(fn)
        rt = str(z["run_type"])
        if "qx_chi2_samples" not in z.files:
            continue
        names, pwo = z["mod_name"], z["mod_pwo4"]
        Y, mi = z["qx_chi2_samples"], z["qx_chi2_midx"].astype(np.int64)
        ped, sig, t0, a, chi, k = (z[f"qx_chi2_{q}"].astype(np.float64) for q in ("ped", "sig", "t", "a", "chi2", "k"))
        t0, a = t0[:, 0], a[:, 0]                        # k == 1 rows: the pulse is column 0
        w = (Y < pr.OVERFLOW).astype(pr.F32)
        rise, fall, _ = pr.MF.mono_flags((Y.astype(pr.F32) - ped.astype(pr.F32)[:, None]) * w, w, sig.astype(pr.F32))
        sel = pwo[mi] & (k == 1) & (chi > 20) & rise & fall & (w > 0).all(axis=1) & (a >= A_MIN)
        for i in np.nonzero(sel)[0]:
            pool[rt].append((fn.stem, names[mi[i]], int(mi[i]), Y[i], ped[i], t0[i], a[i], chi[i], z["qx_chi2_evnum"][i] if "qx_chi2_evnum" in z.files else -1))
    pick = []
    per = N_EX // len(RTS)
    rest = []
    for rt in RTS:
        p = pool[rt]
        o = rng.permutation(len(p))
        pick += [(rt, p[j]) for j in o[:per]]
        rest += [(rt, p[j]) for j in o[per:]]
    if len(pick) < N_EX and rest:
        o = rng.permutation(len(rest))
        pick += [rest[j] for j in o[:N_EX - len(pick)]]
    return pick


def example_model(rt, names, pwo, scan_dir, ent, banks):
    if rt not in banks:
        banks[rt] = (pr.Bank(pr.BANK_ROOT / rt / f"templates_{rt}.npz", names, pwo), load_iint(names, pwo, rt)[0])
    b, im = banks[rt]
    run, nm, mi, Y, ped, t0, a, chi, ev = ent
    mf = pr.MF(b)
    m = np.array([b.row[mi]])
    model = a * mf.basis(m, np.array([t0], pr.F32), Y.size)[0].astype(np.float64)
    y = Y.astype(np.float64) - ped
    r1, r2, f, l = pr.raw_sums(y[None], model[None], np.array([t0]))
    fit = a * im[mi]
    return model + ped, f[0], l[0], fit, r1[0], r2[0]


def stats_row(s, f2, tab):
    n = s[:, 0]
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = s[:, 1] / n
        rms = np.sqrt(np.maximum(s[:, 2] / n - mean ** 2, 0))
        return n, quant(tab, 0.5), rms, f2 / n, quant(tab, 0.16), quant(tab, 0.84)


def sawtooth(H):
    med = quant(H, 0.5)
    ok = np.isfinite(med) & (H.sum(axis=1) >= MIN_BIN)
    return (np.nanmax(med[ok]) - np.nanmin(med[ok])) if ok.sum() > 1 else np.nan, med


def density(ax, H, xe, ylo, yhi, label=None):
    j0, j1 = np.searchsorted(RB, ylo), np.searchsorted(RB, yhi)
    Z = H[:, j0:j1].T
    q = ax.pcolormesh(xe, RB[j0:j1 + 1], np.where(Z > 0, Z, np.nan), norm=LogNorm(vmin=1), cmap="viridis", rasterized=True)
    plt.colorbar(q, ax=ax, label="waveforms per bin")


def median_band(ax, H, xe, log=True, col="r"):
    x = np.sqrt(xe[1:] * xe[:-1]) if log else 0.5 * (xe[1:] + xe[:-1])
    ok = H.sum(axis=1) >= MIN_BIN
    for p, ls in ((0.16, ":"), (0.5, "-"), (0.84, ":")):
        ax.plot(x[ok], quant(H, p)[ok], ls, color=col, lw=1.6 if ls == "-" else 1.0)
    ax.axhline(1, color="w", lw=0.6)


def rt_page(pdf, rt, A, n, names, ttl):
    fig, axs = plt.subplots(2, 3, figsize=(16, 9.5), constrained_layout=True)
    n1, med1, rms1, f21, p161, p841 = stats_row(A["s1"].sum(axis=0, keepdims=True), np.array([A["s1_f2"].sum()]), A["tab1"].sum(axis=0, keepdims=True))
    n2, med2, rms2, f22, p162, p842 = stats_row(A["s2"].sum(axis=0, keepdims=True), np.array([A["s2_f2"].sum()]), A["tab2"].sum(axis=0, keepdims=True))
    # (a)
    ax = axs[0, 0]
    q = ax.pcolormesh(RAW2D, RAW2D, np.where(A["a1"].T > 0, A["a1"].T, np.nan), norm=LogNorm(vmin=1), cmap="viridis", rasterized=True)
    plt.colorbar(q, ax=ax, label="waveforms per bin")
    ax.plot([100, 1e5], [100, 1e5], "w-", lw=0.8)
    ax.set_xscale("log"), ax.set_yscale("log")
    ax.set_xlabel("raw sum, def 1 template window (ADC x sample)"), ax.set_ylabel("fitted integral A I_m (ADC x sample)")
    ax.set_title(f"(a) fitted integral vs raw sum, n = {int(n1[0])}", fontsize=10)
    # (b), (f)
    for ax, t, tag, med, p16, p84 in ((axs[0, 1], "1", "(b) def 1", med1, p161, p841), (axs[1, 2], "2", "(f) def 2", med2, p162, p842)):
        density(ax, A[f"b{t}"], RAWB, *((0.8, 1.2) if t == "1" else (0.7, 1.1)))
        median_band(ax, A[f"b{t}"], RAWB)
        ax.set_xscale("log")
        ax.set_ylim(0.8, 1.2) if t == "1" else ax.set_ylim(0.7, 1.1)
        ax.set_xlabel(f"raw sum, def {t} (ADC x sample)"), ax.set_ylabel("A I_m / raw sum")
        ax.set_title(f"{tag}: median {med[0]:.4f}, 16-84 % {p16[0]:.3f}-{p84[0]:.3f}", fontsize=10)
    # (c)
    ax = axs[0, 2]
    pc = 0.5 * (PHB[1:] + PHB[:-1])
    txt = []
    for t, col in (("1", "tab:blue"), ("2", "tab:orange")):
        sw, mc = sawtooth(A[f"c{t}"])
        ok = A[f"c{t}"].sum(axis=1) >= MIN_BIN
        ax.plot(pc[ok], mc[ok], "o-", color=col, ms=3, label=f"def {t}: peak-to-peak {sw:.4f}")
        txt.append(sw)
    ax.set_xlabel("onset phase t0 mod 4 ns"), ax.set_ylabel("median A I_m / raw sum")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_title("(c) ratio vs onset phase (median per 0.1 ns bin)", fontsize=10)
    # (d)
    ax = axs[1, 0]
    density(ax, A["d1"], CHB, 0.8, 1.2)
    median_band(ax, A["d1"], CHB)
    ax.set_xscale("log"), ax.set_ylim(0.8, 1.2)
    ax.set_xlabel("chi2/dof"), ax.set_ylabel("A I_m / raw sum (def 1)")
    ax.set_title("(d) ratio vs chi2/dof, median and 16-84 %", fontsize=10)
    # (e)
    ax = axs[1, 1]
    nm = A["mod1"][:, 0]
    okm = nm >= MIN_MOD
    mm = quant(A["m1"], 0.5)
    c0 = float(np.nanmedian(mm[okm])) if okm.any() else 1.0
    ax.hist(np.clip(mm[okm], c0 - 0.0999, c0 + 0.0999), bins=np.linspace(c0 - 0.1, c0 + 0.1, 81), color="tab:orange")
    sd = float(np.std(mm[okm])) if okm.any() else np.nan
    ax.set_xlabel(f"module median A I_m / raw sum (def 1), modules with >= {MIN_MOD} waveforms"), ax.set_ylabel("modules")
    ax.set_title(f"(e) {int(okm.sum())} modules: median {np.nanmedian(mm[okm]):.4f}, sd {sd:.4f}", fontsize=10)
    ins = ax.inset_axes([0.58, 0.55, 0.4, 0.4])
    wi = np.nonzero(okm)[0]
    ins.plot(wi, mm[wi], ".", ms=2)
    ins.set_xlabel("module index", fontsize=7), ins.set_ylabel("median ratio", fontsize=7)
    ins.tick_params(labelsize=6)
    ins.set_ylim(c0 - 0.1, c0 + 0.1)
    fig.suptitle(ttl, fontsize=12)
    if PNG_DIR:                                          # quick look: also write the page as a PNG
        fig.savefig(Path(PNG_DIR) / f"page_{rt}.png", dpi=70)
    pdf.savefig(fig)
    plt.close(fig)
    return dict(med1=med1[0], p16_1=p161[0], p84_1=p841[0], med2=med2[0], p16_2=p162[0], p84_2=p842[0], saw1=txt[0], saw2=txt[1], nmod=int(okm.sum()), mod_med=float(np.nanmedian(mm[okm])), mod_sd=sd,
                rms1=rms1[0], rms2=rms2[0], f2_1=f21[0], f2_2=f22[0])


def fit_page(pdf, acc):
    """(a) for every run type on one 2 x 2 page: fitted integral A I_m vs raw sum def 1, same axes in every panel."""
    fig, axs = plt.subplots(2, 2, figsize=(12, 11), constrained_layout=True)
    for ax, rt in zip(axs.flat, RTS):
        A = acc[rt]
        n = int(A["n"][7])
        if n > 0:
            q = ax.pcolormesh(RAW2D, RAW2D, np.where(A["a1"].T > 0, A["a1"].T, np.nan), norm=LogNorm(vmin=1), cmap="viridis", rasterized=True)
            plt.colorbar(q, ax=ax, label="waveforms per bin")
        else:
            ax.text(0.5, 0.5, "no used waveforms", ha="center", va="center", transform=ax.transAxes)
        ax.plot([100, 1e5], [100, 1e5], "--", color="tab:red", lw=0.9)
        ax.set_xscale("log"), ax.set_yscale("log")
        ax.set_xlim(100, 1e5), ax.set_ylim(100, 1e5)
        ax.set_xlabel("raw sum, def 1 template window (ADC x sample)"), ax.set_ylabel("fitted integral A I_m (ADC x sample)")
        ax.set_title(f"(a) {rt}, n = {n}", fontsize=10)
    fig.suptitle(f"(a) fitted integral A I_m vs raw sum def 1, all run types: PbWO4 W*, monotonic k = 1, no overflow, A >= {A_MIN:g} ADC\n"
                 "same axes in every panel, log colour scale, dashed red line y = x", fontsize=11)
    if PNG_DIR:                                          # quick look: also write the page as a PNG
        fig.savefig(Path(PNG_DIR) / "page_fit_vs_raw.png", dpi=70)
    pdf.savefig(fig)
    plt.close(fig)


def all_page(pdf, acc, t):
    """All cases on one 2 x 2 page, after fit_page.  t = '1': fitted A_main I_m vs raw sum def 1; t = '2': total fitted integral sum_i A_i I_m vs raw sum def 2 (fixed window).
    Any k, monotonic or not, overflow and truncated kept; only A_main >= A_MIN.  Same axes as fit_page; each panel has the fraction with k >= 2 and the median ratio."""
    fig, axs = plt.subplots(2, 2, figsize=(12, 11), constrained_layout=True)
    what = "fitted integral A_main I_m vs raw sum def 1" if t == "1" else "total fitted integral sum_i A_i I_m vs raw sum def 2"
    for ax, rt in zip(axs.flat, RTS):
        A = acc[rt]
        n = int(A["xn"][0])
        if n > 0:
            q = ax.pcolormesh(RAW2D, RAW2D, np.where(A[f"x{t}"].T > 0, A[f"x{t}"].T, np.nan), norm=LogNorm(vmin=1), cmap="viridis", rasterized=True)
            plt.colorbar(q, ax=ax, label="waveforms per bin")
            med = quant(A[f"xr{t}"][None], 0.5)[0]
            ax.text(0.03, 0.97, f"k >= 2: {100 * A['xn'][1] / n:.1f} %\nmedian ratio: {med:.4f}", transform=ax.transAxes, va="top", fontsize=9,
                    bbox=dict(fc="white", ec="0.5", lw=0.5))
        else:
            ax.text(0.5, 0.5, "no waveforms", ha="center", va="center", transform=ax.transAxes)
        ax.plot([100, 1e5], [100, 1e5], "--", color="tab:red", lw=0.9)
        ax.set_xscale("log"), ax.set_yscale("log")
        ax.set_xlim(100, 1e5), ax.set_ylim(100, 1e5)
        ax.set_xlabel(f"raw sum, def {t} {'template window' if t == '1' else 'fixed window'} (ADC x sample)")
        ax.set_ylabel("fitted integral A_main I_m (ADC x sample)" if t == "1" else "fitted integral sum_i A_i I_m (ADC x sample)")
        ax.set_title(f"(a) all cases, {rt}, n = {n}", fontsize=10)
    fig.suptitle(f"(a) all cases: {what}; all run types: PbWO4 W*, A_main >= {A_MIN:g} ADC\n"
                 "any k, monotonic or not, overflow and truncated kept; A_main = largest pulse of the channel\n"
                 "same axes in every panel, log colour scale, dashed red line y = x", fontsize=11)
    if PNG_DIR:                                          # quick look: also write the page as a PNG
        fig.savefig(Path(PNG_DIR) / f"page_all_{t}.png", dpi=70)
    pdf.savefig(fig)
    plt.close(fig)


def text_page(pdf, title, lines, size=8):
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.03, 0.97, title, fontsize=12, weight="bold", va="top")
    fig.text(0.03, 0.93, "\n".join(textwrap.fill(l, 170, subsequent_indent="    ") if not l.startswith(" ") else l for l in lines), fontsize=size, va="top", family="monospace")
    pdf.savefig(fig)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("scan_dir")
    ap.add_argument("--stamp", default="20261117")
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--nproc", type=int, default=8)
    a = ap.parse_args()
    sd = Path(a.scan_dir)
    files = sorted(sd.glob("run_*.npz"))
    if not files:
        sys.exit(f"[ERROR] no run_*.npz in {sd}")
    with Pool(min(a.nproc, len(files))) as p:
        res = p.map(reduce_file, files, chunksize=1)
    nmod = res[0][3].size
    acc = {rt: new_acc(nmod) for rt in RTS}
    runs = {rt: [] for rt in RTS}
    evn = {rt: [0, 0] for rt in RTS}
    for rt, ac, run, names, pw, ns, npass in res:
        add_acc(acc[rt], ac)
        runs[rt].append(run)
        evn[rt][0] += int(ns)
        evn[rt][1] += int(npass)
    names = res[0][3]
    rng = np.random.default_rng(SEED)
    out = Path(a.out)
    stem = f"integral_vs_rawsum_full_{a.stamp}"
    # counts, table
    cl = ["n_channels", "k1", "k1_no_overflow", "monotonic", "trunc_main_dropped", "A_below_min_dropped", "bad_value_dropped", "used"]
    with open(out / f"{stem}_counts.tsv", "w") as fh:
        fh.write("run_type\truns\tevents_seen\tevents_after_trigger\t" + "\t".join(cl) + "\n")
        for rt in RTS:
            fh.write(f"{rt}\t{','.join(runs[rt])}\t{evn[rt][0]}\t{evn[rt][1]}\t" + "\t".join(str(int(v)) for v in acc[rt]["n"]) + "\n")
    rows = []
    for rt in RTS + ["all"]:
        A = {k: (sum(acc[r][k] for r in RTS) if rt == "all" else acc[rt][k]) for k in ("s1", "s2", "tab1", "tab2", "s1_f2", "s2_f2")}
        r1, r2 = stats_row(A["s1"], A["s1_f2"], A["tab1"]), stats_row(A["s2"], A["s2_f2"], A["tab2"])
        for i, lab in enumerate(CHI_LAB):
            rows.append(dict(run_type=rt, chi2_bin=lab, n=int(r1[0][i]), med_def1=r1[1][i], rms_def1=r1[2][i], within2pct_def1=r1[3][i], med_def2=r2[1][i], rms_def2=r2[2][i], within2pct_def2=r2[3][i]))
    with open(out / f"{stem}.tsv", "w") as fh:
        keys = list(rows[0])
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[k]:.5g}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")
    pick = load_examples(sd, rng)
    ex_names = {}
    for fn in files[:1]:
        z = np.load(fn)
        ex_names = (z["mod_name"], z["mod_pwo4"])
    with PdfPages(out / f"{stem}.pdf") as pdf:
        # title page
        tot = np.sum([acc[rt]["n"] for rt in RTS], axis=0)
        L = [f"Sample: {sd} ({len(files)} run files; runs " + "; ".join(f"{rt}: {','.join(runs[rt])}" for rt in RTS) + ").",
             "Scan: pileup_rates.py with PILEUP_RAWSUM=1 (default fit, monotonic-edge rule and noisy flag off).  Class: PbWO4 modules (W*) only.",
             f"Selection: k == 1, no overflow sample, monotonic main peak (rise and fall flags of MF.protected, raw samples), main-pulse fit window inside the waveform, fitted A >= {A_MIN:g} ADC, finite chi2/dof, raw sums > 0.",
             "Fitted integral = A I_m (I_m = iint_wa of the module's bank row).  Raw def 1 = sum(Y - ped) over the template walk-out samples at the fitted onset; def 2 = sum(Y - ped) over [t0 - 20, t0 + 200] ns.  Ratio = A I_m / raw sum.",
             "", f"{'run type':16s}{'events':>10s}{'passing':>10s}" + "".join(f"{c[:14]:>15s}" for c in cl)]
        for rt in RTS:
            L.append(f"{rt:16s}{evn[rt][0]:10d}{evn[rt][1]:10d}" + "".join(f"{int(v):15d}" for v in acc[rt]["n"]))
        L.append(f"{'all':16s}{sum(e[0] for e in evn.values()):10d}{sum(e[1] for e in evn.values()):10d}" + "".join(f"{int(v):15d}" for v in tot))
        L += ["", "Medians and 16-84 % come from 0.002-wide ratio histograms; bins with fewer than " + str(MIN_BIN) + " waveforms are not drawn; modules need >= " + str(MIN_MOD) + " waveforms in (e)."]
        summ = {}
        text_page(pdf, "Fitted integral A I_m vs raw sum: PbWO4 monotonic single peaks (" + a.stamp + ")", L)
        for rt in [r for r in RTS if acc[r]["n"][7] > 0]:
            ttl = f"{rt}: PbWO4 W*, monotonic k = 1, no overflow, A >= {A_MIN:g} ADC; n = {int(acc[rt]['n'][7])} waveforms (runs {','.join(runs[rt])}); ratio = A I_m / raw sum"
            summ[rt] = rt_page(pdf, rt, acc[rt], 0, names, ttl)
        fit_page(pdf, acc)
        all_page(pdf, acc, "1")
        all_page(pdf, acc, "2")
        # table pages
        lines = [f"{'run type':16s}{'chi2/dof':>9s}{'n':>11s} | {'med1':>7s}{'rms1':>7s}{'+-2%1':>7s} | {'med2':>7s}{'rms2':>7s}{'+-2%2':>7s}"]
        f = lambda v, w, p=3: f"{v:{w}.{p}f}" if np.isfinite(v) else " " * (w - 2) + "NA"   # noqa: E731
        for rt in RTS + ["all"]:
            for r in rows:
                if r["run_type"] == rt:
                    lines.append(f"{rt:16s}{r['chi2_bin']:>9s}{r['n']:11d} | {f(r['med_def1'], 7)}{f(r['rms_def1'], 7)}{f(100 * r['within2pct_def1'], 7, 1)} | {f(r['med_def2'], 7)}{f(r['rms_def2'], 7)}{f(100 * r['within2pct_def2'], 7, 1)}")
            lines.append("")
        text_page(pdf, "Table: ratio A I_m / raw sum per run type and chi2/dof bin (def 1 template window | def 2 fixed window; +-2 % in percent)", lines, 7)
        # examples
        banks = {}
        fig, axs = plt.subplots(5, 4, figsize=(15, 13), constrained_layout=True)
        tt = pr.DT * np.arange(pr.NS_EXPECT)
        shown = 0
        for ax, (rt, ent) in zip(axs.flat, pick):
            run, nm, mi, Y, ped, t0, amp, chi, ev = ent
            mod, fl, ll, fit, r1, r2 = example_model(rt, ex_names[0], ex_names[1], sd, ent, banks)
            ax.plot(tt, Y, "k.", ms=3.5)
            ax.axhline(ped, color="0.6", lw=0.7)
            ax.plot(tt, mod, "r-", lw=1.0)
            ax.axvspan(pr.DT * fl, pr.DT * ll, color="tab:blue", alpha=0.12)
            ax.axvline(t0 + pr.RS_FIX_LO, color="tab:green", ls="--", lw=0.8)
            ax.axvline(min(t0 + pr.RS_FIX_HI, tt[-1]), color="tab:green", ls="--", lw=0.8)
            ax.set_title(f"{rt} {nm} {run} A {amp:.0f} chi2 {chi:.0f}\nratio def1 {fit / r1:.3f}, def2 {fit / r2:.3f} (fit int {fit:.0f}, raw {r1:.0f} | {r2:.0f})", fontsize=7, loc="left", pad=2)
            ax.tick_params(labelsize=6)
            shown += 1
        fig.suptitle(f"{shown} monotonic k = 1 PbWO4 waveforms with chi2/dof > 20 (scan reservoir qx_chi2, up to {N_EX // len(RTS)} per run type, seed {SEED}).  Black data, red refitted pulse, grey pedestal; "
                     "blue = def 1 window, green dashed = def 2 window.  Reservoir rows carry no trunc_main flag, so a late pulse may appear.  Values outside the axis ranges are clipped into the end bins of the density panels.", fontsize=9)
        pdf.savefig(fig)
        plt.close(fig)
    # findings to stdout
    print("\n".join(L))
    for rt in summ:
        s = summ[rt]
        print(f"{rt}: def1 median {s['med1']:.4f} [16-84 {s['p16_1']:.3f}-{s['p84_1']:.3f}] rms {s['rms1']:.4f} +-2% {100 * s['f2_1']:.1f}% | def2 median {s['med2']:.4f} [{s['p16_2']:.3f}-{s['p84_2']:.3f}] rms {s['rms2']:.4f} +-2% {100 * s['f2_2']:.1f}% | "
              f"phase peak-to-peak def1 {s['saw1']:.4f} def2 {s['saw2']:.4f} | modules {s['nmod']} median {s['mod_med']:.4f} sd {s['mod_sd']:.4f}")
    for rt in RTS:
        A = acc[rt]
        print(f"all cases {rt}: n {int(A['xn'][0])}, k >= 2 {100 * A['xn'][1] / max(A['xn'][0], 1):.1f} %, median A_main I_m / raw def 1 {quant(A['xr1'][None], 0.5)[0]:.4f}, "
              f"median sum A_i I_m / raw def 2 {quant(A['xr2'][None], 0.5)[0]:.4f}")
    print(f"[written] {out / stem}.pdf, .tsv, _counts.tsv")


if __name__ == "__main__":
    main()
