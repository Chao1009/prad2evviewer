#!/usr/bin/env python3
"""calib_figs.py -- five figures for the calibration / energy-resolution section, read only from existing outputs.

Run from analysis/pyscripts:

    python3 waveform_analysis/writeup/make_figs/calib_figs.py

PRAD2_DATABASE_DIR defaults to ~/work/PRad/prad2evviewer/database (beam energies from runinfo/general.json via
energy_resolution_fit.beam_gev; time offsets from hycal_time_offsets/time_calib.json).

Inputs (paths relative to analysis/pyscripts; ER = waveform_analysis/energy_resolution, IMP = ER/improvement,
V3 = output/energy_resolution/control_phase2_v3_20261114):
  V3/v3/3p49GeV_H2/spectra.npz                 E/E' histograms of the v3 chain (runs 24698 + 24911): key 'legacy' = the
                                                ungated reference path of the v3 chain, 'g16' = the 16 ns control gate
  output/energy_resolution/baseline_20261112/pre/3p49GeV_H2/spectra.npz
                                                key 'legacy' = reference path of the legacy chain (runs 24698, 24911, 24937),
                                                drawn as a light grey context line
  IMP/control_phase2_v3_final.tsv              oos res_cb per run type and stage (legacy, v1, v2, v3)
  V3/resfit.out                                v3 resolution fits (cross-check of the refit done here)
  ER/resolution_fit.tsv                        previous three-term fit (older reports)
  ER/<run type>/summary.tsv                    raw-time gate scan (legacy chain), res_cb_oos(_se)
  IMP/P15_tgate_calibrated.tsv                 calibrated-time gate scan (legacy chain), res_cb_oos(_se)
  IMP/control_phase2_v3_gates.tsv              v3 chain near the control gate, res_cb_oos(_se)
  IMP/P5_latest.tsv                            legacy -> best calibrated gate, paired d_oos(_se); best_oos
  IMP/control_phase2.tsv                       v1 control, res_cb_oos
  IMP/control_phase2_v2_delta.tsv              v2 - v1 paired d_res_cb_oos(_se) (P17 on Ta, P9 on 2.24 / 3.49 H2)
  IMP/control_phase2_v3_attr.tsv               pool steps v2 -> +P16 -> +M1, paired d_res_cb_oos(_se)
  crosstalk_study/cache/<run type>/*_legacy.npz  seed and member times (first 3 splits per dir), inputs of
                                                IMP/tools/p15_dt_hist.py, whose figure is redrawn here in the house style

Outputs (vector PDF, waveform_analysis/writeup/figs/):
  calib_elastic_peak_3p49gev_h2.pdf   pooled E/E' at 3.49 GeV, v3 reference path vs v3 16 ns gate, single Crystal Ball fits
  calib_resolution_vs_energy.pdf      oos sigma/E vs beam energy, v3 chain, three-term and two-term fits
  calib_seed_member_dt.pdf            seed-minus-member time, raw vs offset-corrected (redraw of P15_seed_member_dt)
  calib_gate_scan.pdf                 oos sigma/E vs gate width, raw vs calibrated times, one panel per run type
  calib_chain_comparison.pdf          (a) oos sigma/E per chain stage; (b) per-lever paired deltas
"""
import csv
import glob
import json
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "figs"
WA = HERE.parent.parent                                  # waveform_analysis/
PY = WA.parent                                           # analysis/pyscripts/
ER = WA / "energy_resolution"
IMP = ER / "improvement"
V3 = PY / "output" / "energy_resolution" / "control_phase2_v3_20261114"
PRE = PY / "output" / "energy_resolution" / "baseline_20261112" / "pre"
CACHE = PY / "crosstalk_study" / "cache"
os.environ.setdefault("PRAD2_DATABASE_DIR", str(Path.home() / "work/PRad/prad2evviewer/database"))
DB = Path(os.environ["PRAD2_DATABASE_DIR"])
sys.path.insert(0, str(WA))
import energy_resolution_fit as F                        # noqa: E402  (beam_gev, model, model2)

RTS = ["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "3p49GeV_H2_old", "2p24GeV_Ta_X17"]
LABEL = {"0p73GeV_H2": "0.73 GeV H$_2$", "2p24GeV_H2": "2.24 GeV H$_2$", "3p49GeV_H2": "3.49 GeV H$_2$",
         "3p49GeV_H2_old": "3.49 GeV H$_2$ (old runs)", "2p24GeV_Ta_X17": "2.24 GeV Ta (X17)"}
RUN = {"0p73GeV_H2": 24590, "2p24GeV_H2": 25228, "3p49GeV_H2": 24698, "3p49GeV_H2_old": 24386, "2p24GeV_Ta_X17": 25800}
# Okabe-Ito hues (black replaced by a dark ochre), as make_deconv_performance.py
BLUE, ORANGE, GREEN, VERM, SKY, PINK, OCHRE = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7", "#8C6D1F"
RED = "#C0392B"
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5, "legend.fontsize": 7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.6,
                     "lines.linewidth": 1.0, "lines.markersize": 3.2, "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.labelweight": "normal"})   # energy_resolution_fit sets bold labels on import


def tsv(path):
    with open(path) as f:
        return list(csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"))


def grid(ax, axis="y"):
    ax.grid(axis=axis, color="0.88", lw=0.5)
    ax.set_axisbelow(True)


FINAL = {(r["rt"], r["stage"]): r for r in tsv(IMP / "control_phase2_v3_final.tsv")}


def final(rt, stage):
    r = FINAL[(rt, stage)]
    return float(r["res_cb_oos"]), float(r["res_cb_oos_se"]), r["gate"]


# ---------------------------------------------------------------- 1. elastic peak
def crystal_ball(x, A, mu, sigma, alpha, n):
    """Crystal Ball with a left power-law tail (same form as energy_resolution_benchmark.crystal_ball; that module needs
    prad2py and is not importable here)."""
    z = (np.asarray(x, dtype=float) - mu) / sigma
    a = abs(alpha)
    arg = n / a - a - z
    with np.errstate(over="ignore", invalid="ignore"):
        tail = (n / a) ** n * np.exp(-0.5 * a * a) * np.where(arg > 0, arg, 1e-30) ** (-n)
    return A * np.where(z > -alpha, np.exp(-0.5 * z * z), tail)


CB_WIN = 0.18      # fit range [mu0 - CB_WIN, mu0 + CB_WIN / 2] in E/E', as cb_single (400 MeV at 2.24 GeV) of the benchmark


def fit_cb(c, y):
    m0 = (c > 0.85) & (c < 1.12)
    mu0 = float(c[m0][np.argmax(y[m0])])
    m = (c > mu0 - CB_WIN) & (c < mu0 + 0.5 * CB_WIN)
    p, cov = curve_fit(crystal_ball, c[m], y[m].astype(float), p0=[float(y[m].max()), mu0, 0.02, 1.5, 2.0],
                       bounds=([0, mu0 - 0.1, 1e-3, 0.1, 1.01], [np.inf, mu0 + 0.1, CB_WIN, 20.0, 50.0]),
                       sigma=np.sqrt(np.maximum(y[m], 1.0)), maxfev=20000)
    return p, (c > p[1] - CB_WIN) & (c < p[1] + 0.5 * CB_WIN)


def elastic_peak():
    v3 = np.load(V3 / "v3" / "3p49GeV_H2" / "spectra.npz")
    pre = np.load(PRE / "3p49GeV_H2" / "spectra.npz")
    e = v3["edges"]
    assert np.allclose(e, pre["edges"])
    c = 0.5 * (e[:-1] + e[1:])
    curves = [(pre["legacy"], "0.65", "legacy chain, reference path", 0.8, True),
              (v3["legacy"], BLUE, "v3 chain, reference path", 1.2, False),
              (v3["g16"], ORANGE, "v3 chain, 16 ns gate", 1.2, False)]
    print("[1] calib_elastic_peak_3p49gev_h2: pooled single Crystal Ball fits (E/E', in sample)")
    print("    v3 spectra.npz keys:", list(v3.files), "; 'legacy' there = v3 chain without time gate (reference path)")
    print("    legacy-chain reference from baseline_20261112/pre/3p49GeV_H2/spectra.npz (runs 24698, 24911, 24937)")
    fits = []
    for y, col, lab, lw, ctx in curves:
        p, m = fit_cb(c, y)
        fits.append(p)
        print(f"    {lab:30s} N={int(y.sum()):8d}  mu={p[1]:.4f}  sigma/E={100 * p[2] / p[1]:.3f} %  alpha={p[3]:.2f} n={p[4]:.1f}")
    fig, ax = plt.subplots(figsize=(5.4, 3.6), layout="constrained")
    mu_r, sg_r = fits[1][1], fits[1][2]
    for a in (ax, ax.inset_axes([0.6, 0.42, 0.38, 0.52])):
        a.axvspan(mu_r - 2 * sg_r, mu_r + 2 * sg_r, color="0.9", lw=0, zorder=0)
        for (y, col, lab, lw, ctx), p in zip(curves, fits):
            a.stairs(y, e, color=col, lw=lw, zorder=1 if ctx else 2,
                     label=f"{lab}" if a is ax else None)
            xx = np.linspace(p[1] - CB_WIN, p[1] + 0.5 * CB_WIN, 400)
            if not ctx:
                a.plot(xx, crystal_ball(xx, *p), ls="--", color=col, lw=0.9, zorder=3,
                       label=f"CB fit, $\\sigma_E/E$ = {100 * p[2] / p[1]:.2f} %" if a is ax else None)
        if a is not ax:
            a.set_xlim(0.85, 1.12)
            a.set_ylim(0, 1.08 * max(y.max() for y, *_ in curves))
            a.tick_params(labelsize=6)
            a.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
            a.yaxis.get_offset_text().set_fontsize(6)
            a.set_title("linear, 0.85-1.12", fontsize=6.5, pad=2)
    ax.set_yscale("log")
    ax.set_xlim(e[0], e[-1])
    ax.set_ylim(50, 5 * max(y.max() for y, *_ in curves))
    ax.set_xlabel("$E/E'$")
    ax.set_ylabel("clusters / 0.004")
    grid(ax)
    h, l = ax.get_legend_handles_labels()
    order = [0, 1, 3, 2, 4]               # legacy line, then v3 ref + its fit, v3 gate + its fit
    ax.legend([h[i] for i in order], [l[i] for i in order], frameon=False, loc="upper left")
    fig.savefig(OUT / "calib_elastic_peak_3p49gev_h2.pdf")
    plt.close(fig)
    print(f"    grey band: v3 reference-path peak +- 2 sigma = [{mu_r - 2 * sg_r:.4f}, {mu_r + 2 * sg_r:.4f}]")


# ---------------------------------------------------------------- 2. resolution vs energy
def resolution_vs_energy():
    H2 = ["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2"]
    E = np.array([F.beam_gev(RUN[rt]) for rt in H2])
    R, S = np.array([final(rt, "v3")[:2] for rt in H2]).T
    q, cov = curve_fit(F.model, E, R, p0=[1, 1, 1], sigma=S, absolute_sigma=True, bounds=(0, np.inf), maxfev=20000)
    q2, cov2 = curve_fit(F.model2, E, R, p0=[1, 1], sigma=S, absolute_sigma=True, bounds=(0, np.inf), maxfev=20000)
    chi2 = float(np.sum(((R - F.model2(E, *q2)) / S) ** 2))
    # cross-check against resfit.out (first block = the three H2 points)
    txt = (V3 / "resfit.out").read_text().split("points")[1]
    a, b, c_ = (float(re.search(rf"three-term.*?{k} ([\d.]+)", txt).group(1)) for k in ("a", "b", "c"))
    a2, b2 = (float(re.search(rf"two-term.*?{k} ([\d.]+)", txt).group(1)) for k in ("a", "b"))
    chi_out = float(re.search(r"two-term.*?chi2 ([\d.]+)", txt).group(1))
    assert np.allclose([round(x, 3) for x in (*q, *q2)], [a, b, c_, a2, b2]) and round(chi2, 2) == chi_out, "refit != resfit.out"
    old = tsv(ER / "resolution_fit.tsv")[0]
    qo = [float(old[k]) for k in ("a", "b", "c")]
    Eo, (Ro, So) = F.beam_gev(RUN["3p49GeV_H2_old"]), final("3p49GeV_H2_old", "v3")[:2]
    Et, (Rt, St) = F.beam_gev(RUN["2p24GeV_Ta_X17"]), final("2p24GeV_Ta_X17", "v3")[:2]
    RL, SL = np.array([final(rt, "legacy")[:2] for rt in H2]).T
    print("[2] calib_resolution_vs_energy (oos res_cb %, v3 control gate)")
    for e_, r_, s_, rl, rt in zip(E, R, S, RL, H2):
        print(f"    {rt:15s} E={e_:.4f} GeV  v3 {r_:.3f} +- {s_:.3f}   legacy {rl:.3f}")
    print(f"    3p49GeV_H2_old  E={Eo:.4f} GeV  v3 {Ro:.3f} +- {So:.3f} (not fitted)")
    print(f"    2p24GeV_Ta_X17  E={Et:.4f} GeV  v3 {Rt:.3f} +- {St:.3f} (not fitted)")
    print(f"    three-term refit a={q[0]:.3f} b={q[1]:.3f} c={q[2]:.3f} (= resfit.out); two-term a={q2[0]:.3f} b={q2[1]:.3f} "
          f"chi2={chi2:.1f}/1; previous three-term (resolution_fit.tsv) {qo[0]:.2f} / {qo[1]:.2f} / {qo[2]:.2f}")
    fig, ax = plt.subplots(figsize=(4.6, 3.4), layout="constrained")
    x = np.linspace(0.5, 4.0, 300)
    ax.plot(x, F.model(x, *q), color=RED, lw=1.3,
            label=f"three-term fit: {q[0]:.2f} $\\oplus$ {q[1]:.2f}/$\\sqrt{{E}}$ $\\oplus$ {q[2]:.2f}/$E$")
    ax.plot(x, F.model(x, *qo), color="0.45", ls="--", lw=1.0,
            label=f"previous three-term: {qo[0]:.2f} $\\oplus$ {qo[1]:.2f}/$\\sqrt{{E}}$ $\\oplus$ {qo[2]:.2f}/$E$")
    ax.plot(x, F.model2(x, *q2), color=GREEN, ls=":", lw=1.4,
            label=f"two-term fit: {q2[0]:.2f} $\\oplus$ {q2[1]:.2f}/$\\sqrt{{E}}$ ($\\chi^2$ = {chi2:.1f}/1)")
    ax.errorbar(E, RL, yerr=SL, fmt="o", color="0.7", ms=4, capsize=1.5, label="H$_2$, legacy chain")
    ax.errorbar(E, R, yerr=S, fmt="o", color=BLUE, ms=4, capsize=1.5, label="H$_2$, v3 chain (fitted)")
    ax.errorbar([Eo], [Ro], yerr=[So], fmt="s", mfc="white", mec=BLUE, color=BLUE, ms=4, capsize=1.5,
                label="3.49 GeV old runs, v3 (not fitted)")
    ax.errorbar([Et], [Rt], yerr=[St], fmt="*", color=VERM, ms=8, capsize=1.5, label="2.24 GeV Ta X17, v3 (not fitted)")
    ax.set_xlim(0.5, 4.0)
    ax.set_ylim(1.2, 3.8)
    ax.set_xlabel("beam energy $E$ (GeV)")
    ax.set_ylabel("out-of-sample $\\sigma_E/E$ (%)")
    grid(ax, "both")
    ax.legend(frameon=False, fontsize=6.5, loc="upper right")
    fig.savefig(OUT / "calib_resolution_vs_energy.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- 3. seed-member dt (redraw of tools/p15_dt_hist.py)
DT_PANELS = [("0p73GeV_H2", ["0p73GeV_H2"], 12), ("2p24GeV_H2", ["2p24GeV_H2"], 16),
             ("2p24GeV_Ta_X17", ["2p24GeV_Ta_X17_25800", "2p24GeV_Ta_X17_26138", "2p24GeV_Ta_X17_26230"], 4),
             ("3p49GeV_H2", ["3p49GeV_H2"], 12)]     # P15 one-SE gates, as the producer
NSPLIT = 3
DT_BINS = np.arange(-20, 30.01, 0.5)


def seed_member_dt():
    off = {m["name"]: float(m["offset_ns"]) for m in json.load(open(DB / "hycal_time_offsets" / "time_calib.json"))["modules"]
           if m.get("resolved")}
    rsig = lambda v: 0.5 * np.subtract(*np.percentile(v, [84.135, 15.865]))   # noqa: E731
    print("[3] calib_seed_member_dt (redraw of improvement/figures/P15_seed_member_dt by tools/p15_dt_hist.py; same inputs)")
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 5.0), sharex=True, sharey=True, layout="constrained")
    for ax, (rt, dirs, gate) in zip(axs.flat, DT_PANELS):
        raw, cal = [], []
        for d in dirs:
            for f in sorted(glob.glob(str(CACHE / d / "*_legacy.npz")))[:NSPLIT]:
                z = np.load(f)
                o = np.array([off.get(str(n), 0.0) for n in z["module_names"]])
                mi, si = z["mod_i"].astype(np.int64), z["seed_mod_i"].astype(np.int64)
                ts, t = z["seed_t"].astype(float), z["t"].astype(float)
                raw.append(ts - t)
                cal.append((ts - o[si]) - (t - o[mi]))
        raw, cal = np.concatenate(raw), np.concatenate(cal)
        sr, sc = rsig(raw), rsig(cal)
        fr = lambda v, lo, hi: 100 * np.mean((v > lo) & (v <= hi))                # noqa: E731
        print(f"    {rt:15s} N={raw.size / 1e6:.2f} M  robust sigma raw {sr:.2f} cal {sc:.2f} ns  "
              f"dt<-8: raw {fr(raw, -1e9, -8):.2f} % cal {fr(cal, -1e9, -8):.2f} %  dt>8: raw {fr(raw, 8, 1e9):.2f} % cal {fr(cal, 8, 1e9):.2f} %  "
              f"7<dt<=10: raw {fr(raw, 7, 10):.2f} % cal {fr(cal, 7, 10):.2f} %  |dt|>{gate}: raw {100 * np.mean(abs(raw) > gate):.2f} % "
              f"cal {100 * np.mean(abs(cal) > gate):.2f} %")
        hr, hc = np.histogram(raw, DT_BINS)[0], np.histogram(cal, DT_BINS)[0]
        ax.stairs(hr, DT_BINS, color="0.5", lw=1.0, label=f"raw, robust $\\sigma$ = {sr:.2f} ns")
        ax.stairs(hc, DT_BINS, color=BLUE, lw=1.1, label=f"offset-corrected, robust $\\sigma$ = {sc:.2f} ns")
        ax.axvline(gate, color=VERM, ls="--", lw=0.9, label=f"P15 one-SE gate {gate} ns")
        ax.set_yscale("log")
        ax.set_xlim(-20, 30)
        ax.set_ylim(0.5, max(hr.max(), hc.max()) * 60)
        ax.set_title(f"{LABEL[rt]} ({raw.size / 1e6:.2f} M members)")
        grid(ax, "both")
        ax.legend(frameon=False, loc="upper left", fontsize=6.5)
    for ax in axs[1]:
        ax.set_xlabel("$\\Delta t = t_\\mathrm{seed} - t_\\mathrm{member}$ (ns)")
    for ax in axs[:, 0]:
        ax.set_ylabel("members / 0.5 ns")
    fig.savefig(OUT / "calib_seed_member_dt.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- 4. gate scan
def gscan(rows, key="path"):
    out = {}
    for r in rows:
        if r[key].startswith("g") and np.isfinite(float(r["res_cb_oos"])):
            out[float(r[key][1:])] = (float(r["res_cb_oos"]), float(r["res_cb_oos_se"]))
    return out


def gate_scan():
    p15 = tsv(IMP / "P15_tgate_calibrated.tsv")
    v3g = tsv(IMP / "control_phase2_v3_gates.tsv")
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.6), layout="constrained")
    print("[4] calib_gate_scan (oos res_cb %; raw = ER/<RT>/summary.tsv, cal = P15_tgate_calibrated.tsv, both legacy chain)")
    print(f"    {'run type':15s} {'gate':>4s}  {'raw':>12s}  {'calibrated':>12s}  {'v3 chain':>12s}")
    for ax, rt in zip(axs.flat, RTS):
        S = tsv(ER / rt / "summary.tsv")
        leg = next(r for r in S if r["path"] == "legacy")
        L, Ls = float(leg["res_cb_oos"]), float(leg["res_cb_oos_se"])
        raw = gscan(S)
        cal = gscan([r for r in p15 if r["run_type"] == rt])
        v3 = {float(r["path"][1:]): (float(r["res_cb_oos"]), float(r["res_cb_oos_se"])) for r in v3g if r["run_type"] == rt}
        ctrl = float(final(rt, "v3")[2][1:])
        for g in (4, 8, 12, 16):
            f = lambda d: f"{d[g][0]:.3f}+-{d[g][1]:.3f}" if g in d else "-"     # noqa: E731
            print(f"    {rt:15s} {g:4d}  {f(raw):>12s}  {f(cal):>12s}  {f(v3):>12s}")
        print(f"    {rt:15s} legacy (ungated) {L:.3f}+-{Ls:.3f}; control gate {ctrl:g} ns")
        vis = [v for d in (cal, v3) for g, (v, _) in d.items() if g >= 3] + [L]
        lo, hi = min(vis), max(vis)
        top = hi + 0.25 * (hi - lo)
        ax.axhspan(L - Ls, L + Ls, color="0.9", lw=0)
        ax.axhline(L, color="0.4", lw=0.8, ls="-", label="no time gate (legacy chain)")
        ax.axvline(ctrl, color=VERM, ls="--", lw=0.8, label="control gate")
        for d, col, mk, lab in ((raw, "0.45", "o", "raw times"), (cal, BLUE, "s", "calibrated times"),
                                (v3, GREEN, "D", "v3 chain, calibrated")):
            g = np.array(sorted(d))
            v = np.array([d[x][0] for x in g])
            s = np.array([d[x][1] for x in g])
            inn = v <= top
            ms = 2.2 if d is raw else 3.0
            ax.errorbar(g[inn], v[inn], yerr=s[inn], fmt=mk, ms=ms, color=col, lw=0.8, capsize=0, label=lab,
                        mfc="white" if d is v3 else col)
            if d is not v3:
                ax.plot(g, np.minimum(v, top), color=col, lw=0.7)
                if (~inn).any():
                    ax.plot(g[~inn], np.full((~inn).sum(), top), marker="^", ls="none", ms=3, color=col, clip_on=False)
        ax.set_ylim(lo - 0.1 * (hi - lo), top)
        ax.set_xlim(0, 25)
        ax.set_xticks([0, 4, 8, 12, 16, 20, 24])
        ax.set_title(LABEL[rt])
        grid(ax)
    for ax in axs[1, :2]:
        ax.set_xlabel("time gate (ns)")
    axs[0, 2].set_xlabel("time gate (ns)")
    for ax in axs[:, 0]:
        ax.set_ylabel("out-of-sample $\\sigma_E/E$ (%)")
    h, l = axs[0, 0].get_legend_handles_labels()
    axs[1, 2].axis("off")
    axs[1, 2].legend(h, l, frameon=False, loc="center", fontsize=7.5,
                     title="triangles at the top edge:\nvalues above the axis", title_fontsize=6.5)
    fig.savefig(OUT / "calib_gate_scan.pdf")
    plt.close(fig)


# ---------------------------------------------------------------- 5. chain comparison and lever attribution
def chain_comparison():
    stages = [("legacy", "0.6", "o"), ("v1", SKY, "s"), ("v2", ORANGE, "^"), ("v3", BLUE, "D")]
    P5 = {r["run_type"]: r for r in tsv(IMP / "P5_latest.tsv")}
    C1 = {r["run_type"]: r for r in tsv(IMP / "control_phase2.tsv")}
    D2 = {r["run_type"]: r for r in tsv(IMP / "control_phase2_v2_delta.tsv")}
    AT = {(r["rt"], r["step_to"].split(" ")[0]): r for r in tsv(IMP / "control_phase2_v3_attr.tsv")
          if r["region"] == "pool" and r["step_from"].split(" ")[0] in ("v2", "+P16")}
    lev = {}                                            # (lever, rt) -> (delta, se or nan)
    for rt in RTS:
        lev["calibrated gate (P15)", rt] = (float(P5[rt]["d_oos"]), float(P5[rt]["d_oos_se"]))
        if rt not in ("0p73GeV_H2", "2p24GeV_Ta_X17"):  # log weight applied at 2.24 / 3.49 only; unpaired difference
            lev["log weight (P4b/P5)", rt] = (float(C1[rt]["res_cb_oos"]) - float(P5[rt]["best_oos"]), np.nan)
        d = (float(D2[rt]["d_res_cb_oos"]), float(D2[rt]["d_res_cb_oos_se"]))
        if rt == "2p24GeV_Ta_X17":
            lev["Ta 4 ns gate (P17)", rt] = d
        elif rt != "0p73GeV_H2":
            lev["shape factor (P9)", rt] = d
        lev["below-threshold (P16)", rt] = (float(AT[rt, "+P16"]["d_res_cb_oos"]), float(AT[rt, "+P16"]["d_res_cb_oos_se"]))
        lev["phase residual (M1)", rt] = (float(AT[rt, "+M1"]["d_res_cb_oos"]), float(AT[rt, "+M1"]["d_res_cb_oos_se"]))
    levers = ["calibrated gate (P15)", "log weight (P4b/P5)", "Ta 4 ns gate (P17)", "shape factor (P9)",
              "below-threshold (P16)", "phase residual (M1)"]
    lcol = [BLUE, SKY, VERM, ORANGE, GREEN, PINK]
    print("[5] calib_chain_comparison: oos res_cb % per stage (control_phase2_v3_final.tsv)")
    for rt in RTS:
        print("    " + f"{rt:15s} " + "  ".join(f"{s} {final(rt, s)[0]:.3f}+-{final(rt, s)[1]:.3f} ({final(rt, s)[2]})" for s, *_ in stages))
    print("    lever deltas (pp, oos, negative = better); closure = sum of levers vs v3 - legacy")
    for rt in RTS:
        ds = {k: v for (k, r), v in lev.items() if r == rt}
        tot = final(rt, "v3")[0] - final(rt, "legacy")[0]
        print(f"    {rt:15s} " + "  ".join(f"{k.split(' (')[1][:-1]} {v[0]:+.4f}" + (f"+-{v[1]:.4f}" if np.isfinite(v[1]) else "(no SE)")
                                         for k, v in ds.items())
              + f"  | sum {sum(v[0] for v in ds.values()):+.4f} vs total {tot:+.4f}")
    fig = plt.figure(figsize=(7.2, 5.6), layout="constrained")
    sf = fig.subfigures(2, 1, height_ratios=[1, 1.15])
    top = sf[0].subplots(1, 5)
    for ax, rt in zip(top, RTS):
        for i, (s, col, mk) in enumerate(stages):
            v, e, g = final(rt, s)
            ax.errorbar([i], [v], yerr=[e], fmt=mk, color=col, ms=4, capsize=1.5, lw=0.9)
        ax.plot(range(4), [final(rt, s)[0] for s, *_ in stages], color="0.75", lw=0.6, zorder=0)
        ax.set_xticks(range(4), [s for s, *_ in stages], fontsize=6.5)
        ax.set_xlim(-0.5, 3.5)
        ax.set_title(LABEL[rt], fontsize=7.5)
        grid(ax)
    top[0].set_ylabel("out-of-sample $\\sigma_E/E$ (%)")
    sf[0].suptitle("(a) chain stages: legacy $\\to$ v1 (gate, log weight) $\\to$ v2 (P17, P9) $\\to$ v3 (P16, M1)", fontsize=8)
    ax = sf[1].subplots()
    x = np.arange(len(RTS))
    w = 0.13
    for j, (k, col) in enumerate(zip(levers, lcol)):
        xs, vs, es = [], [], []
        for i, rt in enumerate(RTS):
            if (k, rt) in lev:
                xs.append(i + (j - 2.5) * w)
                vs.append(lev[k, rt][0])
                es.append(lev[k, rt][1] if np.isfinite(lev[k, rt][1]) else 0.0)
        bars = ax.bar(xs, vs, w * 0.92, color=col, label=k, hatch="////" if "log weight" in k else None,
                      edgecolor="white" if "log weight" in k else col, lw=0)
        ax.errorbar(xs, vs, yerr=es, fmt="none", ecolor="0.2", lw=0.6, capsize=1)
    ax.axhline(0, color="0.3", lw=0.6)
    ax.set_ylim(-0.36, 0.02)
    ax.set_xticks(x, [LABEL[rt] for rt in RTS])
    ax.set_ylabel("$\\Delta$ out-of-sample $\\sigma_E/E$ (pp)")
    ax.set_title("(b) per-lever paired change (negative = better); hatched: unpaired difference, no SE", fontsize=8)
    grid(ax)
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 0.02), fontsize=6.5)
    fig.savefig(OUT / "calib_chain_comparison.pdf")
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    elastic_peak()
    resolution_vs_energy()
    seed_member_dt()
    gate_scan()
    chain_comparison()
