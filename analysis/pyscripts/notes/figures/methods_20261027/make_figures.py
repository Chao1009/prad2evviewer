#!/usr/bin/env python3
"""make_figures.py -- figure set for the PbWO4 (HyCal, PRad-II) energy-resolution methods note, 2026-10-27.

Run from anywhere (reads existing outputs only, no replay, nothing under database/):

    MPLBACKEND=Agg python3 analysis/pyscripts/notes/figures/methods_20261027/make_figures.py [fig01 fig02 ...]

Writes PNGs next to this file.  New plots (this script):
  fig01_dataflow           matplotlib box diagram (same chain as dataflow.mmd; mmdc is not installed)
  fig02_waveform           one real waveform: pedestal, WaveAnalyzer walk-out range, linear CFD 0.5, MF fit
  fig03_template_banks     one module, production / best-50 / Wassim single-pulse bank
  fig04a_pileup_example    pseudo-data pile-up (dt_peak = 20 ns, ratio 1), MF decomposition
  fig04b_pileup_eff_cells  60-cell detection-efficiency heat maps, one panel per bank
  fig06_position_corr      position-correction map: bilinear run-type map vs hybrid per-module poly8
  fig07_resolution_vs_E    sigma_E/E vs E with the published 2-term curve
  fig08_module_sigma_map   composite of two existing plot_module_resolution.py PNGs (PIL paste)
Existing figure copied unchanged (shell, from analysis/pyscripts):
  fig05_elastic_spectrum.png = output/scan_20261025/pos_hybrid/scan_025320/split_00050/result_calibrated.png
      cp <that file> notes/figures/methods_20261027/fig05_elastic_spectrum.png
Reads python modules read-only: mf_pulse.py, benchmark_hycal_timing_test_20261021.py (poscorr_load/eval), position_corr.py.
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                    # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle           # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.normpath(os.path.join(HERE, "..", "..", ".."))        # analysis/pyscripts
sys.path.insert(0, P)
os.chdir(P)

BLUE, ORANGE, GREEN, GREY, RED = "#2a78d6", "#eb6834", "#1baf7a", "#6b6b6b", "#c8302b"
INK = "#222222"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#888888", "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": INK, "ytick.color": INK, "figure.dpi": 100})


def save(fig, name, dpi=130):
    path = os.path.join(HERE, name)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    kb = os.path.getsize(path) / 1024
    if kb > 280:                                   # keep under the 300 kB limit
        fig2 = None
        from PIL import Image
        im = Image.open(path)
        im = im.convert("P", palette=Image.ADAPTIVE, colors=128)
        im.save(path, optimize=True)
        kb = os.path.getsize(path) / 1024
    print(f"{name}: {kb:.0f} kB")


# ---------------------------------------------------------------- banks -----
BANKS = {
    "prod": ("output/templates_spline/2p24GeV_H2/templates_2p24GeV_H2.npz",
             "production (first-come 150 / class, R3)", BLUE),
    "b50": ("output/templates_spline_best50_20261026/2p24GeV_H2/templates_2p24GeV_H2.npz",
            "best-50 (50 highest-energy pulses)", ORANGE),
    "wass": ("output/templates_wassim_20261027/2p24GeV_H2/templates_2p24GeV_H2.npz",
             "Wassim single pulse", GREEN),
}


def load_bank(path):
    import mf_pulse as M
    z = np.load(path)
    names = [str(n) for n in z["names"]]
    bank = M.SampledTemplateBank.from_npz(path, names, [True] * len(names), smooth_order=2,
                                          int_tail_ratio=0.1, tail_break_n=2)
    return bank, names, z


# ---------------------------------------------------------------- fig01 -----
def fig01():
    fig, ax = plt.subplots(figsize=(11.5, 6.4))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 65)
    ax.axis("off")

    def box(x, y, w, h, title, sub="", col="#eaf1fb", edge=BLUE):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=0.8",
                                    fc=col, ec=edge, lw=1.2))
        ax.text(x + w / 2, y + h - 1.6, title, ha="center", va="top", fontsize=9, weight="bold")
        if sub:
            ax.text(x + w / 2, y + h - 4.4, sub, ha="center", va="top", fontsize=7.4, color="#333333",
                    linespacing=1.25)

    def arrow(x0, y0, x1, y1, txt=""):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1.2, shrinkA=0, shrinkB=0))
        if txt:
            ax.text((x0 + x1) / 2, (y0 + y1) / 2 + 0.9, txt, ha="center", fontsize=7, color=GREY)

    W, H = 20, 12.5
    y1, y2, y3 = 44, 24, 4
    # row 1
    box(1, y1, W, H, "Raw FADC waveform", "100 samples x 4 ns\n(250 MHz), one per\nchannel and event")
    box(26, y1, W, H, "Pedestal", "first 30 samples,\niterative clip\n-> mean, rms")
    box(51, y1 + 7, 24, 8.5, "WaveAnalyzer (C++)", "smooth; local max; walk-out\nintegral; linear CFD 0.5\n+ log-normal time", col="#fdf0e8", edge=ORANGE)
    box(51, y1 - 4, 24, 8.5, "Matched filter MF (Python)", "sampled template bank,\nup to 4 pulses: A, onset t0", col="#e6f6ef", edge=GREEN)
    arrow(21.4, y1 + 6, 25.6, y1 + 6)
    arrow(46.4, y1 + 8, 50.6, y1 + 11)
    arrow(46.4, y1 + 4, 50.6, y1 + 0)
    box(80, y1, 19, H, "Pulse list", "integral (rawsum) or\nA x I_int x f(A) x kappa\n(MF); time per pulse")
    arrow(75.4, y1 + 11, 79.6, y1 + 8)
    arrow(75.4, y1 + 0, 79.6, y1 + 4)
    # row 2
    box(1, y2, W, H, "Module energy", "calibration file -> MeV;\nLMS gain correction;\nper-crystal gains")
    box(26, y2, W, H, "Clustering", "legacy: largest pulse in\npre-window, no time cut\ngated: |t - t_seed| <= W")
    box(51, y2, W, H, "Cluster energy", "HyCalCluster island,\nseed + neighbours\n(single-cluster events)")
    box(76, y2, 23, H, "Position correction", "E -> E / R(u, v), u, v =\nimpact in seed crystal;\nbilinear or hybrid poly8")
    arrow(90, y1 - 0.4, 12, y2 + H + 0.4)
    arrow(21.4, y2 + 6, 25.6, y2 + 6)
    arrow(46.4, y2 + 6, 50.6, y2 + 6)
    arrow(71.4, y2 + 6, 75.6, y2 + 6)
    # row 3
    box(76, y3, 23, H, "Elastic selection", "seed r < 200 mm, seed E cut,\nrow/col exclusion, elastic\nwindow; seed-indexed gains")
    box(51, y3, W, H, "Resolution fit", "post-refinement peak:\nGaussian (iterative +-2 sigma)\nand Crystal Ball")
    box(26, y3, W, H, "Per-split numbers", "sigma_E/E per split\n(Gaussian, CB),\nyields")
    box(1, y3, W, H, "Per-module and headline", "module medians;\nrun-type mean over splits;\n2-term fit a, b", col="#f3f3f3", edge=GREY)
    arrow(87.5, y2 - 0.4, 87.5, y3 + H + 0.4)
    arrow(75.6, y3 + 6, 71.4, y3 + 6)
    arrow(50.6, y3 + 6, 46.4, y3 + 6)
    arrow(25.6, y3 + 6, 21.4, y3 + 6)
    ax.text(1, 64.5, "Analysis chain per EVIO split (benchmark_hycal_timing_test_20261021.py, driven by benchmark_hycal_timing_scan.py)",
            fontsize=9.5, weight="bold", va="top")
    save(fig, "fig01_dataflow.png")


# ---------------------------------------------------------------- fig02 -----
def fig02():
    import mf_pulse as M
    from prad2py import dec
    d = np.load("output/templates_wassim_20261027/2p24GeV_H2/pulses/025320_00050.npz", allow_pickle=True)
    ok = (d["wa_npeaks"] == 1) & (d["tier"] == 2) & (d["amp"] > 1000) & (d["amp"] < 1500)
    j = int(np.flatnonzero(ok)[0])
    raw = d["samples"][j].astype(float)
    name = str(d["name"][j])
    wa = dec.WaveAnalyzer(dec.WaveConfig())
    ped, rms, pks = wa.analyze(d["samples"][j])
    pk = list(pks)[0]
    nsamp = len(raw)
    t = 4.0 * np.arange(nsamp)
    # linear CFD at 0.5 on the smoothed pedestal-subtracted samples (WaveAnalyzer.cpp 526-545)
    buf = M.smooth_tri(raw[None, :], 2)[0] - ped
    height = raw[pk.pos] - ped
    thr = 0.5 * height
    tcfd = None
    for jj in range(pk.left - 1, pk.pos):
        v0, v1 = buf[jj], buf[jj + 1]
        if v0 < thr <= v1 and (v1 - v0) > 1.5 * rms:
            tcfd = 4.0 * (jj + (thr - v0) / (v1 - v0))
    bank, names, _ = load_bank(BANKS["prod"][0])
    rec = M.MFReconstructor(bank, nsigma=5.0, min_amp=5.0, kmax=4, r_model=0.3, ped_nsamples=30)
    m = names.index(name)
    res = rec.run(d["samples"][j:j + 1], np.array([m]))
    a, t0 = float(res.p_amp[0]), float(res.p_t0[0])
    tf = np.arange(0, 4.0 * nsamp, 0.25)
    model = res.ped[0] + a * bank.shape(m, tf - t0, dtype=np.float64)
    model_s = res.ped[0] + a * bank.shape(m, t - t0, dtype=np.float64)

    fig = plt.figure(figsize=(10.5, 4.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 1], hspace=0.08, wspace=0.22)
    for col, (lo, hi, ttl) in enumerate([(0, 400, "full readout window"), (120, 260, "pulse region")]):
        ax = fig.add_subplot(gs[0, col])
        axr = fig.add_subplot(gs[1, col], sharex=ax)
        ax.plot(t, raw, "o", ms=3.2, color=INK, label="raw samples (ADC)")
        ax.plot(tf, model, "-", color=GREEN, lw=1.6, label=f"MF fit: A = {a:.0f} ADC, onset t0 = {t0:.1f} ns")
        ax.axhline(ped, color=GREY, lw=1, ls=":")
        ax.axvspan(0, 4.0 * 30, color="#dddddd", alpha=0.45, lw=0)
        ax.axvspan(4.0 * pk.left, 4.0 * pk.right, color=ORANGE, alpha=0.18, lw=0)
        ax.plot([4.0 * pk.left, 4.0 * pk.right], [ped + 0.5 * height * 0 - 25] * 2, "|", color=ORANGE, ms=8)
        ax.axhline(ped + thr, color=BLUE, lw=0.8, ls="--")
        if tcfd is not None:
            ax.axvline(tcfd, color=BLUE, lw=1.3)
        ax.set_xlim(lo, hi)
        ax.set_ylim(ped - 0.06 * height, ped + 1.12 * height)
        ax.set_title(f"{name}, run 25320 split 50, {ttl}", fontsize=9)
        ax.tick_params(labelbottom=False)
        ax.set_ylabel("ADC")
        if col == 0:
            ax.text(60, ped + 0.25 * height, "pedestal\nsamples 0-29", ha="center", fontsize=7.5, color=GREY)
            ax.text(ped * 0 + 260, ped + 0.12 * height, f"pedestal {ped:.1f} ADC, rms {rms:.2f}", fontsize=7.5, color=GREY)
        else:
            ax.text(hi - 3, ped + thr + 0.02 * height, "0.5 x height", ha="right", fontsize=7.5, color=BLUE)
        axr.axhline(0, color=GREY, lw=0.8)
        axr.plot(t, raw - model_s, "o", ms=2.8, color=INK)
        axr.set_xlabel("time (ns)")
        axr.set_ylabel("data - fit")
        axr.axvspan(4.0 * pk.left, 4.0 * pk.right, color=ORANGE, alpha=0.18, lw=0)
    ax = fig.axes[0]
    h = [ax.lines[0], ax.lines[1],
         Rectangle((0, 0), 1, 1, fc="#dddddd", alpha=0.6),
         Rectangle((0, 0), 1, 1, fc=ORANGE, alpha=0.25),
         plt.Line2D([], [], color=BLUE, lw=1.3)]
    lab = ["raw samples", f"MF fit (A = {a:.0f} ADC, t0 = {t0:.1f} ns)", "pedestal window (first 30 samples)",
           f"WaveAnalyzer walk-out range, samples {pk.left}-{pk.right} (inclusive bounds; integral {pk.integral:.0f} ADC.sample)",
           f"linear CFD 0.5 crossing, t = {tcfd:.1f} ns (pk.time incl. log-normal refit: {pk.time:.1f} ns)"]
    fig.legend(h, lab, loc="lower center", ncol=1, fontsize=7.8, frameon=False, bbox_to_anchor=(0.5, -0.17))
    save(fig, "fig02_waveform.png")
    print("fig02:", name, "amp", d["amp"][j], "pk", pk.pos, pk.left, pk.right, pk.time, "cfd", tcfd, "MF", a, t0)


# ---------------------------------------------------------------- fig03 -----
def fig03(module=None):
    banks = {k: load_bank(v[0]) for k, v in BANKS.items()}
    if module is None:
        cand = []
        for nm in banks["prod"][1]:
            ok = True
            for k, (b, names, z) in banks.items():
                r = names.index(nm)
                ok &= (z["fallback"][r] == 0) and (z["valid"][r] > 0)
            if ok:
                cand.append(nm)
        module = "W458" if "W458" in cand else cand[len(cand) // 2]
    tt = np.arange(-8, 160, 0.25)
    fig, (ax, axd) = plt.subplots(2, 1, figsize=(7.2, 5.4), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.3], hspace=0.07))
    ref = None
    for k, (b, names, z) in banks.items():
        r = names.index(module)
        y = b.shape(r, tt, dtype=np.float64)
        n = int(z["n_pool"][r])
        lab = f"{BANKS[k][1]}, N = {n}"
        ax.plot(tt, y, color=BANKS[k][2], lw=1.5, label=lab)
        if ref is None:
            ref = y
        else:
            axd.plot(tt, 100 * (y - ref), color=BANKS[k][2], lw=1.2)
    ax.set_ylabel("unit-peak template")
    ax.set_title(f"{module}, 2.24 GeV H2 banks (t = 0 at the 0.5 % foot of the template)", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    axd.axhline(0, color=GREY, lw=0.8)
    axd.set_ylabel("bank - production\n(% of peak)")
    axd.set_xlabel("time after onset (ns)")
    save(fig, "fig03_template_banks.png")
    print("fig03 module", module)


# ---------------------------------------------------------------- fig04a ----
def fig04a():
    import mf_pulse as M
    bank, names, _ = load_bank(BANKS["prod"][0])
    rec = M.MFReconstructor(bank, nsigma=5.0, min_amp=5.0, kmax=4, r_model=0.3, ped_nsamples=30)
    d = np.load("output/pileup_eff_20261024_fix/2p24GeV_H2/synth/pileup_PbWO4_dt020_ratio1000.npz", allow_pickle=True)
    nm = [str(x) for x in d["base_channel_names"]]
    idx = np.array([names.index(n) for n in nm])
    W = d["waveforms"]
    res = rec.run(W, idx)
    pk_of = {}
    for r_, c in enumerate(res.p_ch):
        pk_of.setdefault(int(c), []).append(r_)
    resolved = next(i for i in range(len(nm)) if len(pk_of.get(i, [])) == 2)
    unres = next(i for i in range(len(nm)) if len(pk_of.get(i, [])) == 1)
    n2 = sum(1 for i in range(len(nm)) if len(pk_of.get(i, [])) >= 2)
    tf = np.arange(0, 400, 0.25)
    t = 4.0 * np.arange(100)
    fig, axs = plt.subplots(1, 2, figsize=(10.2, 3.9), sharey=False)
    for ax, i, ttl in ((axs[0], resolved, "two MF pulses found"), (axs[1], unres, "one MF pulse found")):
        m = idx[i]
        ped = float(res.ped[i])
        ax.plot(t, W[i], "o", ms=3.2, color=INK, label="pseudo-data samples")
        tot = np.full_like(tf, ped)
        for r_ in pk_of[i]:
            a, t0 = float(res.p_amp[r_]), float(res.p_t0[r_])
            comp = a * bank.shape(m, tf - t0, dtype=np.float64)
            tot += comp
            ax.plot(tf, ped + comp, "--", color=GREEN, lw=1.0)
        ax.plot(tf, tot, "-", color=GREEN, lw=1.6, label="MF total")
        for k in range(2):
            ax.axvline(float(d["truth_onset_sampled_ns"][i, k]), color=ORANGE, lw=1.0, ls=":",
                       label="truth onsets" if k == 0 else None)
        info = ", ".join(f"A = {float(res.p_amp[r_]):.0f}, t0 = {float(res.p_t0[r_]):.1f} ns" for r_ in pk_of[i])
        ax.set_title(f"{nm[i]}: {ttl}\ntruth heights {d['truth_heights_adc'][i, 0]:.0f} + {d['truth_heights_adc'][i, 1]:.0f} ADC; MF: {info}", fontsize=8)
        ax.set_xlim(110, 260)
        ax.set_xlabel("time (ns)")
        ax.set_ylabel("ADC")
    axs[0].legend(frameon=False, fontsize=7.5, loc="center left")
    save(fig, "fig04a_pileup_example.png")
    print(f"fig04a: {n2}/{len(nm)} events with >= 2 MF pulses in this file; resolved {nm[resolved]} ev {resolved}, unresolved {nm[unres]} ev {unres}")


# ---------------------------------------------------------------- fig04b ----
def fig04b():
    import csv
    rows = list(csv.DictReader(open("output/templates_wassim_20261027/compare/pileup_3way_cells.tsv"), delimiter="\t"))
    blk = "r_model 0.3 (headline)"
    dts = sorted({float(r["dt_ns"]) for r in rows})
    rts = sorted({float(r["ratio"]) for r in rows})
    grid = {k: np.zeros((len(rts), len(dts))) for k in ("R3", "B50", "W")}
    cnt = np.zeros((len(rts), len(dts)))
    for r in rows:
        if r["block"] != blk:
            continue
        i, j = rts.index(float(r["ratio"])), dts.index(float(r["dt_ns"]))
        if r["arm"] == "B":
            grid["R3"][i, j] += float(r["eff_ref"])
            grid["B50"][i, j] += float(r["eff_new"])
            cnt[i, j] += 1
        else:
            grid["W"][i, j] += float(r["eff_new"])
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.6), sharey=True)
    for ax, (k, ttl) in zip(axs, (("R3", "production (R3)"), ("B50", "best-50"), ("W", "Wassim single pulse"))):
        g = grid[k] / cnt
        im = ax.imshow(g, origin="lower", cmap="Blues", vmin=0, vmax=100, aspect="auto")
        for i in range(len(rts)):
            for j in range(len(dts)):
                ax.text(j, i, f"{g[i, j]:.0f}", ha="center", va="center", fontsize=6.5,
                        color="white" if g[i, j] > 60 else INK)
        ax.set_xticks(range(len(dts)))
        ax.set_xticklabels([f"{x:g}" for x in dts], fontsize=7)
        ax.set_yticks(range(len(rts)))
        ax.set_yticklabels([f"{x:g}" for x in rts])
        ax.set_xlabel("peak separation dt (ns)")
        ax.set_title(ttl, fontsize=9)
        ax.spines[:].set_visible(False)
    axs[0].set_ylabel("amplitude ratio")
    cb = fig.colorbar(im, ax=axs, fraction=0.025, pad=0.02)
    cb.set_label("detection efficiency (%)")
    save(fig, "fig04b_pileup_eff_cells.png")


# ---------------------------------------------------------------- fig06 -----
def fig06():
    import benchmark_hycal_timing_test_20261021 as B
    base = "output/poscorr_20261025"
    bil = B.poscorr_load(f"{base}/2p24GeV_H2/poscorr_2p24GeV_H2_all.json")
    hyb = B.poscorr_load(f"{base}/hybrid/2p24GeV_H2/poscorr_hybrid_2p24GeV_H2_all.json")
    hj = json.load(open(f"{base}/hybrid/2p24GeV_H2/poscorr_hybrid_2p24GeV_H2_all.json"))
    own = hj["modules"]
    best = max(own, key=lambda k: own[k]["n"])
    g = np.linspace(-0.6, 0.6, 121)
    U, V = np.meshgrid(g, g)                       # rows = v, cols = u
    Rb, _ = B.poscorr_eval(bil, U.ravel(), V.ravel())
    Rb = Rb.reshape(U.shape)
    acc = np.zeros(U.shape)
    for nm_ in own:
        Rm, _ = B.poscorr_eval(hyb, U.ravel(), V.ravel(), seed=np.array([nm_] * U.size))
        acc += Rm.reshape(U.shape)
    Rh_avg = acc / len(own)
    Rm, _ = B.poscorr_eval(hyb, U.ravel(), V.ravel(), seed=np.array([best] * U.size))
    Rh_one = Rm.reshape(U.shape)
    lim = 1.6
    fig, axs = plt.subplots(1, 3, figsize=(11.2, 3.7), sharey=True)
    for ax, R, ttl in ((axs[0], Rb, "bilinear run-type map\n(12 x 12 bin grid, all seeds)"),
                       (axs[1], Rh_avg, f"hybrid: mean of the {len(own)} own-module poly8 maps"),
                       (axs[2], Rh_one, f"hybrid: module {best} only (N = {own[best]['n']})")):
        im = ax.pcolormesh(g, g, 100 * (R - 1), cmap="RdBu_r", vmin=-lim, vmax=lim, shading="auto", rasterized=True)
        ax.add_patch(Rectangle((-0.5, -0.5), 1, 1, fill=False, ec="#444444", lw=0.8, ls="--"))
        ax.set_aspect("equal")
        ax.set_title(ttl, fontsize=8.5)
        ax.set_xlabel("u = (x_cl - x_seed) / size")
    axs[0].set_ylabel("v = (y_cl - y_seed) / size")
    cb = fig.colorbar(im, ax=axs, fraction=0.025, pad=0.02)
    cb.set_label("R(u, v) - 1  (%)")
    fig.suptitle("Position correction E -> E / R(u, v), 2.24 GeV H2 (dashed box: seed crystal)", fontsize=9, y=1.02)
    save(fig, "fig06_position_corr.png")
    print("fig06 best module", best, own[best]["n"], "range avg", Rh_avg.min(), Rh_avg.max())


# ---------------------------------------------------------------- fig07 -----
def fig07():
    # summary.md section 3.6, row "R3, legacy, CB for all": sigma/E 2.974 / 1.806 / 1.492 %, a = 2.49, b = 0.68, chi2 3.4
    # (same numbers: 20261027 section 2.3, row "R3, CB for all": a = 2.49 +- 0.02, b = 0.68 +- 0.03)
    E = np.array([0.727, 2.236, 3.484])
    s = np.array([2.974, 1.806, 1.492])
    a, b = 2.49, 0.68
    x = np.linspace(0.4, 4.0, 200)
    fig, ax = plt.subplots(figsize=(5.8, 4.1))
    ax.plot(x, np.hypot(a / np.sqrt(x), b), color=BLUE, lw=1.8, label="a / sqrt(E) (+) b")
    ax.plot(x, a / np.sqrt(x), color=BLUE, lw=0.9, ls="--")
    ax.axhline(b, color=BLUE, lw=0.9, ls=":")
    ax.plot(E, s, "o", ms=7, color=ORANGE, mec="white", mew=1.2, label="R3 legacy, CB for all (H2)", zorder=5)
    # Ta (X17) runs are never combined: three separate per-run R3 points (legacy, CB, mean +- SE over splits),
    # output/resolution_fit_20261027_best50/compare/r4_vs_r3_splits.tsv column R3_CB_post; not in the a, b fit
    for k_, (run_, v_, e_) in enumerate((("25800", 1.986, 0.028), ("26138", 1.950, 0.114), ("26230", 1.983, 0.066))):
        xe = 2.236 + (-0.28, 0.0, 0.28)[k_]
        ax.errorbar(xe, v_, yerr=e_, fmt="s", ms=5, mfc="white", color=GREY, capsize=2, lw=1,
                    label="Ta runs 25800 / 26138 / 26230, separate, all at 2.24 GeV (x offset for display; not fitted)" if k_ == 0 else None)
        ax.annotate(run_, (xe, v_), textcoords="offset points", xytext=(0, 12 if k_ != 1 else 22), ha="center", fontsize=6.5, color=GREY)
    for e, v in zip(E, s):
        dx, dy = ((6, 6), (-8, -18), (6, 8))[list(E).index(e)]
        ax.annotate(f"{v:.3f} %", (e, v), textcoords="offset points", xytext=(dx, dy), fontsize=8)
    ax.text(3.0, a / np.sqrt(3.0) - 0.2, "a / sqrt(E)", ha="center", fontsize=7.5, color=BLUE)
    ax.text(3.9, b + 0.07, "b", ha="right", fontsize=7.5, color=BLUE)
    ax.set_xlabel("beam energy E (GeV)")
    ax.set_ylabel("sigma_E / E  (%)")
    ax.set_xlim(0.3, 4.0)
    ax.set_ylim(0, 3.9)
    ax.set_title("a = 2.49 +- 0.02 % sqrt(GeV), b = 0.68 +- 0.03 %, chi2 = 3.4 (1 dof)", fontsize=8.5)
    ax.legend(frameon=False, fontsize=8)
    save(fig, "fig07_resolution_vs_E.png")


# ---------------------------------------------------------------- fig08 -----
def fig08():
    from PIL import Image
    base = "output/resolution_fit_20261027_best50/module_resolution/r3_ref"
    a = Image.open(f"{base}/2p24_H2/2p24_H2_sigma_post.png").convert("RGB")
    b = Image.open(f"{base}/3p49_H2/3p49_H2_sigma_post.png").convert("RGB")
    h = min(a.height, b.height)
    out = Image.new("RGB", (a.width + b.width, h), "white")
    out.paste(a.crop((0, 0, a.width, h)), (0, 0))
    out.paste(b.crop((0, 0, b.width, h)), (a.width, 0))
    path = os.path.join(HERE, "fig08_module_sigma_map.png")
    out.convert("P", palette=Image.ADAPTIVE, colors=200).save(path, optimize=True)
    print(f"fig08_module_sigma_map.png: {os.path.getsize(path) / 1024:.0f} kB")


ALL = dict(fig01=fig01, fig02=fig02, fig03=fig03, fig04a=fig04a, fig04b=fig04b, fig06=fig06, fig07=fig07, fig08=fig08)

if __name__ == "__main__":
    for k in (sys.argv[1:] or ALL):
        ALL[k]()
