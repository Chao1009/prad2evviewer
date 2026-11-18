#!/usr/bin/env python3
"""t1x_highamp_residual.py -- T1 extension: coherent fit residual of single pulses up to the top of the ADC range (PbWO4), per run type.

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/t1x_highamp_residual.py [--hi-cap N] [--rebuild-hi] [--out <dir>]

Sample: PbWO4 S1 waveforms (default fit, k == 1) from the T1 cache (resid_common, A 50-2000 ADC used here) plus the high-amplitude pass of resid_common.build_hi
(A 2000-4095 ADC, drawn from the 2-file sample scans with at most --hi-cap rows per run type and amplitude bin, refit).  Waveforms with overflow samples (4095) are kept in
the high bins; the overflow samples are masked in every sum (the fit gives them weight 0).  x = (data - ped - fit) / A at tau = t_j - t0, 4 ns bins -20..200 ns, phase pooled.
Amplitude bins: 50-200, 200-800, 800-2000, 2000-3000, 3000-4095 ADC.  0p73 GeV has no pulses above ~800 ADC.
-> T1x_highamp_residual_20261116.pdf / .tsv (profiles) and _peakint.tsv (peak / integral per amplitude bin)
"""
from __future__ import annotations

import sys
import argparse
import textwrap
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
import numpy as np                                      # noqa: E402
import matplotlib                                       # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                         # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages    # noqa: E402
import resid_common as rc                               # noqa: E402

TAU = np.arange(-20.0, 200.1, 4.0)
NT = TAU.size
WIN = (-12.0, 116.0)
AE = np.array([50.0, 200.0, 800.0, 2000.0, 3000.0, 4096.0])
AL = ["50-200", "200-800", "800-2000", "2000-3000", "3000-4095"]
NA = len(AL)
RTS = ["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17"]
COL = dict(zip(RTS, ["tab:green", "tab:red", "tab:blue", "tab:orange"]))
ACOL = plt.cm.viridis(np.linspace(0.0, 0.9, NA))
FE = np.array([50, 100, 200, 400, 800, 1200, 1600, 2000, 2500, 3000, 3500, 4095.0001])      # peak / integral bins (as T1x_chi2_vs_amplitude)
FL = [f"{int(a)}-{int(b) if b < 4096 else 4095}" for a, b in zip(FE[:-1], FE[1:])]
OVF = 4095
MINN = 100          # samples per tau bin for a curve point
SEED = 20261116


def assemble(a):
    """PbWO4 waveforms: cache with A < 2000 + high pass.  Dict of arrays with 'ovf' (waveform has an overflow sample) and 'src'."""
    c = rc.load_s1(a)
    h = rc.load_hi(a)
    pc = (c["cls"] == 0) & (c["a"] < AE[3])
    ph = np.asarray(h["a"] >= AE[3])
    keys = ("run", "ev", "mi", "Y", "ped", "sig", "a", "t0", "chi2", "resid", "var", "rt")
    d = {k: np.concatenate([c[k][pc], h[k][ph]]) for k in keys}
    d["ovf"] = np.concatenate([np.zeros(int(pc.sum()), bool), h["ovf"][ph]])
    d["src"] = np.concatenate([np.zeros(int(pc.sum()), np.int8), np.ones(int(ph.sum()), np.int8)])
    d["names"] = c["names"]
    d["ab"] = (np.digitize(d["a"], AE[1:-1])).astype(np.int8)
    d["rti"] = np.array([RTS.index(r) for r in d["rt"]], np.int8)
    d["tr"] = (c, h)
    return d


def accumulate(d, use):
    """Sums [rt, ab, tau] over the samples of the waveforms in `use`, overflow samples masked: n, sx, sxx, sv (V / A^2), sm (model / A)."""
    S = {k: np.zeros(len(RTS) * NA * NT) for k in ("n", "sx", "sxx", "sv", "sm")}
    idx = np.nonzero(use)[0]
    for i0 in range(0, idx.size, 20000):
        s = idx[i0:i0 + 20000]
        A = d["a"][s].astype(np.float64)[:, None]
        x = d["resid"][s] / A
        v = d["var"][s] / A ** 2
        mod = (d["Y"][s].astype(np.float64) - d["ped"][s][:, None] - d["resid"][s]) / A
        tau = rc.TT[None, :] - d["t0"][s][:, None]
        ib = np.floor((tau + 22.0) / 4.0).astype(np.int64)
        ok = (ib >= 0) & (ib < NT) & (d["Y"][s] < OVF)
        g = (((d["rti"][s].astype(np.int64) * NA + d["ab"][s])[:, None]) * NT + ib)[ok]
        for key, w in (("n", None), ("sx", x[ok]), ("sxx", x[ok] ** 2), ("sv", v[ok]), ("sm", mod[ok])):
            S[key] += np.bincount(g, weights=w, minlength=S[key].size)
    return {k: v.reshape(len(RTS), NA, NT) for k, v in S.items()}


def stats(S):
    with np.errstate(invalid="ignore", divide="ignore"):
        n = S["n"]
        mean = S["sx"] / n
        std = np.sqrt(np.maximum(S["sxx"] / n - mean ** 2, 0.0))
        return mean, std / np.sqrt(n), std, np.sqrt(S["sv"] / n), S["sm"] / n


def coherent(S, r, b):
    """Numbers of one (run type, amp bin): tau of the model peak, mean r/A there (+- sem), the deepest mean after the peak (fall: tau_pk + 8 .. + 60 ns), rms / assumed."""
    mean, sem, rms, asg, mod = stats(S)
    n = S["n"][r, b]
    ok = n >= MINN
    if ok.sum() < 5:
        return None
    ip = int(np.nanargmax(np.where(ok, mod[r, b], -np.inf)))
    fall = np.nonzero(ok & (TAU >= TAU[ip] + 8) & (TAU <= TAU[ip] + 60))[0]
    jf = fall[int(np.nanargmin(mean[r, b][fall]))]
    rise = np.nonzero(ok & (TAU >= TAU[ip] - 28) & (TAU <= TAU[ip] - 4))[0]
    jr = rise[int(np.nanargmax(mean[r, b][rise]))] if rise.size else ip
    inwin = ok & (TAU >= WIN[0]) & (TAU <= WIN[1])
    return dict(tau_pk=TAU[ip], pk=100 * mean[r, b, ip], pk_sem=100 * sem[r, b, ip], mod_pk=mod[r, b, ip],
                tau_rise=TAU[jr], rise=100 * mean[r, b, jr], rise_sem=100 * sem[r, b, jr],
                tau_fall=TAU[jf], fall=100 * mean[r, b, jf], fall_sem=100 * sem[r, b, jf],
                ratio=float(np.nanmedian(rms[r, b][inwin] / asg[r, b][inwin])), rms_win=100 * float(np.nanmean(rms[r, b][inwin])), asg_win=100 * float(np.nanmean(asg[r, b][inwin])),
                n_mid=int(n[ip]))


def peak_integral(d):
    """Per (run type, finer A bin): mean peak / integral of data and of the fit over the fit window (32 samples), waveforms without overflow.  Returns dict arrays [rt, bin]."""
    use = np.nonzero(~d["ovf"])[0]
    out = {k: np.full((len(RTS), len(FL)), np.nan) for k in ("n", "pd", "pm", "sd", "sm", "rd", "rm", "rsem")}
    rat_d, rat_m, fb, ri = [], [], [], []
    for i0 in range(0, use.size, 20000):
        s = use[i0:i0 + 20000]
        wm = rc.window_mask(d["t0"][s])
        dat = d["Y"][s].astype(np.float64) - d["ped"][s][:, None]
        mod = dat - d["resid"][s]
        pk = lambda z: np.where(wm, z, -np.inf).max(axis=1)                         # noqa: E731
        it = lambda z: np.where(wm, z, 0.0).sum(axis=1)                             # noqa: E731
        rat_d.append(pk(dat) / it(dat))
        rat_m.append(pk(mod) / it(mod))
        fb.append(np.digitize(d["a"][s], FE) - 1)
        ri.append(d["rti"][s])
    rat_d, rat_m, fb, ri = map(np.concatenate, (rat_d, rat_m, fb, ri))
    for r in range(len(RTS)):
        for b in range(len(FL)):
            m = (ri == r) & (fb == b)
            if m.sum() < 30:
                continue
            out["n"][r, b] = m.sum()
            out["rd"][r, b], out["rm"][r, b] = rat_d[m].mean(), rat_m[m].mean()
            out["rsem"][r, b] = (rat_d[m] - rat_m[m]).std() / np.sqrt(m.sum())
    return out


def write_tsv(path, S):
    mean, sem, rms, asg, mod = stats(S)
    with open(path, "w") as f:
        f.write("run_type\tamp_bin\ttau\tmean\tsem\trms\tassumed_sigma\tmodel_mean\tn\n")
        for r, rt in enumerate(RTS):
            for b in range(NA):
                for i in range(NT):
                    n = S["n"][r, b, i]
                    if n > 0:
                        f.write(f"{rt}\t{AL[b]}\t{TAU[i]:.0f}\t{mean[r, b, i]:.6g}\t{sem[r, b, i]:.3g}\t{rms[r, b, i]:.6g}\t{asg[r, b, i]:.6g}\t{mod[r, b, i]:.6g}\t{int(n)}\n")


def write_pi_tsv(path, PI):
    with open(path, "w") as f:
        f.write("run_type\tA_bin\tn\tpeak_over_integral_data\tpeak_over_integral_fit\tdata_minus_fit\tsem_of_difference\n")
        for r, rt in enumerate(RTS):
            for b in range(len(FL)):
                if not np.isnan(PI["n"][r, b]):
                    f.write(f"{rt}\t{FL[b]}\t{int(PI['n'][r, b])}\t{PI['rd'][r, b]:.6g}\t{PI['rm'][r, b]:.6g}\t{PI['rd'][r, b] - PI['rm'][r, b]:.4g}\t{PI['rsem'][r, b]:.3g}\n")


def table_lines(d, S, S0):
    L = [f"{'run type':<15}{'A bin':<10}{'nwf':>7}{'nmod':>6}{'tau_pk':>7}{'mean@peak %A':>17}{'rise max %A':>16}{'fall min %A (tau)':>22}{'nclip0 only: peak / fall':>26}{'rms/asg':>8}{'med chi2':>9}"]
    for r, rt in enumerate(RTS):
        for b in range(NA):
            sel = (d["rti"] == r) & (d["ab"] == b)
            if sel.sum() == 0:
                continue
            c, c0 = coherent(S, r, b), coherent(S0, r, b)
            if c is None:
                L.append(f"{rt:<15}{AL[b]:<10}{int(sel.sum()):>7}  too few samples")
                continue
            s0 = f"{c0['pk']:+6.2f} / {c0['fall']:+6.2f}" if c0 else "n/a"
            L.append(f"{rt:<15}{AL[b]:<10}{int(sel.sum()):>7}{len(np.unique(d['mi'][sel])):>6}{c['tau_pk']:>7.0f}{c['pk']:>+9.2f} +-{c['pk_sem']:5.2f}{c['rise']:>+9.2f} +-{c['rise_sem']:5.2f}"
                     f"{c['fall']:>+9.2f} +-{c['fall_sem']:5.2f} ({c['tau_fall']:3.0f}){s0:>26}{c['ratio']:>8.2f}{np.median(d['chi2'][sel]):>9.2f}")
    return L


def decision_lines(S, d, PI):
    L = ["Decision input (PbWO4, coherent mean of r / A, overflow samples masked; peak = tau bin of the largest fitted pulse model, fall = deepest mean 8..60 ns after it):"]
    for r, rt in enumerate(RTS):
        cs = [(b, coherent(S, r, b)) for b in range(NA)]
        cs = [(b, c) for b, c in cs if c]
        if len(cs) < 2:
            continue
        L.append(f"  {rt}: " + "; ".join(f"A {AL[b]}: peak {c['pk']:+.2f} %, fall {c['fall']:+.2f} %" for b, c in cs))
    L.append("  Equal-A comparison (peak / fall, % of A):")
    for b in range(NA):
        parts = []
        for r, rt in enumerate(RTS):
            c = coherent(S, r, b)
            if c and S["n"][r, b].max() >= 2000:
                parts.append(f"{rt} {c['pk']:+.2f}+-{c['pk_sem']:.2f} / {c['fall']:+.2f}+-{c['fall_sem']:.2f}")
        L.append(f"   A {AL[b]}: " + "; ".join(parts))
    L.append("  Peak / integral over the fit window, data minus fit (waveforms without overflow), by run type:")
    for r, rt in enumerate(RTS):
        ok = ~np.isnan(PI["n"][r])
        if ok.sum() == 0:
            continue
        L.append(f"   {rt}: " + "; ".join(f"{FL[b]}: {PI['rd'][r, b] - PI['rm'][r, b]:+.4f}+-{PI['rsem'][r, b]:.4f}" for b in np.nonzero(ok)[0] if PI["n"][r, b] >= 30))
    return L


def page_text(pdf, d, S, S0, PI, hi, a):
    fig = plt.figure(figsize=(17, 10))
    fig.text(0.02, 0.97, "T1x: coherent residual of single pulses vs amplitude up to 4095 ADC (PbWO4)", fontsize=16, fontweight="bold")
    poolk = [f"{k.split(':')[0]} {('2000-3000', '3000-4095')[int(k.split(':')[1])]}: pool {int(p)}, drawn {int(q)}" for k, p, q in zip(hi["pool_keys"], hi["pool_n"], hi["pick_n"])]
    nr = {rt: int((d['rti'] == i).sum()) for i, rt in enumerate(RTS)}
    head = (f"Sample: T1 cache (A < 2000, first-N-events pass, default fit) + high-amplitude pass ({int(hi['hi_cap'])} rows max per run type and bin; pool / drawn rows: " + "; ".join(poolk) + f"; "
            f"{int(hi['wall'])} s wall; {hi['n_trunc']} waveforms dropped with the fit window outside the record).  Waveforms: " + ", ".join(f"{k} {v}" for k, v in nr.items()) +
            f";  {int(d['ovf'].sum())} with overflow samples (masked).\nCommand: " + " ".join(sys.argv))
    dec = "\n".join(textwrap.fill(ln, 215, subsequent_indent="      ") for ln in decision_lines(S, d, PI))
    txt = textwrap.fill(head, 215) + "\n\n" + "\n".join(table_lines(d, S, S0)) + "\n\n" + dec
    fig.text(0.02, 0.94, txt, fontsize=6.3, va="top", family="monospace", linespacing=1.25)
    pdf.savefig(fig)
    plt.close(fig)


def page_runtype(pdf, S, r):
    mean, sem, rms, asg, mod = stats(S)
    rt = RTS[r]
    fig, ax = plt.subplots(4, 1, figsize=(11, 13), sharex=True, gridspec_kw=dict(height_ratios=[3, 3, 1.6, 1.6]))
    fig.suptitle(f"PbWO4, {rt}: residual / A vs tau per amplitude bin (overflow samples masked)", fontsize=12, fontweight="bold")
    for b in range(NA):
        ok = S["n"][r, b] >= MINN
        if ok.sum() < 5:
            continue
        lab = f"A {AL[b]}  (n wf-samples at peak {int(S['n'][r, b].max())})"
        ax[0].plot(TAU[ok], 100 * mean[r, b][ok], ".-", color=ACOL[b], lw=1, label=lab)
        ax[0].fill_between(TAU[ok], 100 * (mean[r, b] - 2 * sem[r, b])[ok], 100 * (mean[r, b] + 2 * sem[r, b])[ok], color=ACOL[b], alpha=0.2, lw=0)
        ax[1].plot(TAU[ok], 100 * rms[r, b][ok], "-", color=ACOL[b], lw=1.2)
        ax[1].plot(TAU[ok], 100 * asg[r, b][ok], "--", color=ACOL[b], lw=1.0)
        ax[2].plot(TAU[ok], (rms / asg)[r, b][ok], ".-", color=ACOL[b], lw=1)
        ax[3].plot(TAU[ok], mod[r, b][ok], "-", color=ACOL[b], lw=1)
    for x in ax:
        x.axvspan(*WIN, color="0.94", lw=0, zorder=0)
        x.axvline(0, color="tab:blue", ls=":", lw=0.8)
        x.grid(alpha=0.25)
    ax[0].axhline(0, color="r", lw=0.7)
    ax[0].set_ylabel("mean r / A  (% of A; band +- 2 sem)")
    ax[0].legend(fontsize=8)
    ax[1].set_ylabel("sigma / A (%): solid rms, dashed assumed")
    ax[1].set_yscale("log")
    ax[2].axhline(1, color="r", lw=0.7)
    ax[2].set_ylabel("rms / assumed")
    ax[2].set_ylim(0.3, 10)
    ax[2].set_yscale("log")
    ax[3].set_ylabel("fitted pulse / A")
    ax[3].set_xlabel("tau = t_j - t0 (ns)   (grey: nominal fit window)")
    ax[3].set_xlim(-22, 202)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    pdf.savefig(fig)
    plt.close(fig)


def page_overlay(pdf, S, kind):
    mean, sem, rms, asg, mod = stats(S)
    fig, axs = plt.subplots(3, 2, figsize=(15, 11), sharex=True)
    fig.suptitle("Equal amplitude, run types overlaid (PbWO4): " + ("mean r / A (%), band +- 2 sem" if kind == "mean" else "rms / assumed sigma"), fontsize=12, fontweight="bold")
    for b in range(NA):
        ax = axs[b // 2, b % 2]
        for r, rt in enumerate(RTS):
            ok = S["n"][r, b] >= MINN
            if ok.sum() < 5:
                continue
            if kind == "mean":
                ax.plot(TAU[ok], 100 * mean[r, b][ok], ".-", color=COL[rt], lw=1, label=f"{rt} (n {int(S['n'][r, b].max())})")
                ax.fill_between(TAU[ok], 100 * (mean[r, b] - 2 * sem[r, b])[ok], 100 * (mean[r, b] + 2 * sem[r, b])[ok], color=COL[rt], alpha=0.2, lw=0)
            else:
                ax.plot(TAU[ok], (rms / asg)[r, b][ok], ".-", color=COL[rt], lw=1, label=rt)
        ax.axvspan(*WIN, color="0.94", lw=0, zorder=0)
        ax.axhline(0 if kind == "mean" else 1, color="r", lw=0.6)
        ax.set_title(f"A {AL[b]} ADC", fontsize=10, loc="left")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7)
        if kind != "mean":
            ax.set_ylim(0.3, 10)
            ax.set_yscale("log")
        if b >= 3:
            ax.set_xlabel("tau (ns)")
    axs[2, 1].axis("off")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)


def page_peakint(pdf, PI):
    fig, axs = plt.subplots(1, 3, figsize=(17, 6))
    x = np.arange(len(FL))
    fig.suptitle("PbWO4, waveforms without overflow: peak / integral over the 32-sample fit window, data vs fitted template", fontsize=12, fontweight="bold")
    for r, rt in enumerate(RTS):
        ok = ~np.isnan(PI["n"][r])
        if not ok.any():
            continue
        axs[0].plot(x[ok], PI["rd"][r][ok], "o-", color=COL[rt], label=f"{rt} data")
        axs[0].plot(x[ok], PI["rm"][r][ok], "s--", color=COL[rt], mfc="none", label=f"{rt} fit")
        axs[1].errorbar(x[ok], (PI["rd"][r] - PI["rm"][r])[ok], yerr=PI["rsem"][r][ok], color=COL[rt], marker="o", capsize=2, label=rt)
        axs[2].plot(x[ok], PI["rd"][r][ok] / PI["rm"][r][ok] - 1, "o-", color=COL[rt], label=rt)
    axs[1].axhline(0, color="k", lw=0.7)
    axs[2].axhline(0, color="k", lw=0.7)
    for ax, yl in zip(axs, ("peak / integral", "data - fit", "data / fit - 1")):
        ax.set_xticks(x)
        ax.set_xticklabels(FL, rotation=40, fontsize=8)
        ax.set_ylabel(yl)
        ax.set_xlabel("fitted A (ADC)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    pdf.savefig(fig)
    plt.close(fig)


def pages_examples(pdf, d, S):
    rng = np.random.default_rng(SEED)
    mean, sem, rms, asg, mod = stats(S)
    quota = {("2p24GeV_H2", 3): 4, ("2p24GeV_H2", 4): 4, ("2p24GeV_Ta_X17", 3): 3, ("3p49GeV_H2", 3): 3, ("3p49GeV_H2", 4): 3, ("2p24GeV_Ta_X17", 4): 0}
    picks = []
    for (rt, b), q in quota.items():
        pool = np.nonzero((d["rt"] == rt) & (d["ab"] == b))[0]
        picks += [int(i) for i in rng.choice(pool, size=min(q, pool.size), replace=False)]
    have = set(picks)
    if len(picks) < 20:                               # fill from the pooled high bins
        pool = np.array([i for i in np.nonzero(d["a"] >= 2000)[0] if int(i) not in have])
        picks += [int(i) for i in rng.choice(pool, size=min(20 - len(picks), pool.size), replace=False)]
    picks = picks[:20]
    PER = 6
    for p0 in range(0, len(picks), PER):
        fig = plt.figure(figsize=(16, 12.5))
        fig.text(0.03, 0.985, va="top", s=textwrap.fill(f"Example waveforms above 2000 ADC ({p0 + 1}-{min(p0 + PER, len(picks))} of {len(picks)}), random per run type and amplitude bin (seed {SEED}); lower strip = data - fit (open red circles: "
                              "overflow samples, weight 0); grey band = measured +- rms of the run type / amplitude bin (x A), red dashed = assumed +- sqrt(V)", 190), fontsize=9.5)
        outer = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.12, left=0.05, right=0.99, top=0.92, bottom=0.04)
        for q, i in enumerate(picks[p0:p0 + PER]):
            gs = outer[q // 2, q % 2].subgridspec(2, 1, height_ratios=[2.2, 1.4], hspace=0.08)
            ax, ar = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
            r, b = int(d["rti"][i]), int(d["ab"][i])
            A, t0, ped = float(d["a"][i]), float(d["t0"][i]), float(d["ped"][i])
            model = d["Y"][i] - ped - d["resid"][i]
            rc.panel(ax, d["Y"][i], ped, model, t0, f"{d['rt'][i]} {d['names'][d['mi'][i]]} run {d['run'][i]} ev {d['ev'][i]}  A {A:.0f}  t0 {t0:.1f}  chi2/dof {d['chi2'][i]:.2f}" + ("  [overflow]" if d["ovf"][i] else ""))
            ax.tick_params(labelbottom=False)
            band = np.interp(rc.TT - t0, TAU, rms[r, b], left=np.nan, right=np.nan) * A
            sg = np.sqrt(d["var"][i])
            ar.fill_between(rc.TT, -band, band, color="0.8", label="measured +- rms")
            ar.plot(rc.TT, sg, "r--", lw=1, label="assumed +- sqrt(V)")
            ar.plot(rc.TT, -sg, "r--", lw=1)
            ov = d["Y"][i] >= OVF
            ar.plot(rc.TT[~ov], d["resid"][i][~ov], "k.-", lw=0.6, ms=3)
            ar.plot(rc.TT[ov], d["resid"][i][ov], "o", mfc="none", mec="r", ms=5)
            ar.axhline(0, color="r", lw=0.5)
            ar.axvline(t0, color="tab:blue", ls=":", lw=0.8)
            ar.set_xlim(0, rc.TT[-1])
            ar.tick_params(labelsize=7)
            ar.set_xlabel("t (ns)", fontsize=7)
            if q == 0:
                ar.legend(fontsize=6, loc="upper right")
        pdf.savefig(fig)
        plt.close(fig)
    return len(picks)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    rc.add_common_args(ap)
    rc.add_hi_args(ap)
    a = ap.parse_args()
    d = assemble(a)
    c, hi = d.pop("tr")
    S = accumulate(d, np.ones(d["a"].size, bool))
    S0 = accumulate(d, ~d["ovf"])
    PI = peak_integral(d)
    lines = table_lines(d, S, S0) + [""] + decision_lines(S, d, PI)
    print("\n".join(lines))
    stem = Path(a.out) / f"T1x_highamp_residual_{rc.STAMP}"
    write_tsv(str(stem) + ".tsv", S)
    write_pi_tsv(str(stem) + "_peakint.tsv", PI)
    with PdfPages(str(stem) + ".pdf") as pdf:
        page_text(pdf, d, S, S0, PI, hi, a)
        for r in (1, 3, 2, 0):
            page_runtype(pdf, S, r)
        page_overlay(pdf, S, "mean")
        page_overlay(pdf, S, "ratio")
        page_peakint(pdf, PI)
        ne = pages_examples(pdf, d, S)
    print(f"[done] {stem}.pdf / .tsv / _peakint.tsv, {ne} examples")


if __name__ == "__main__":
    main()
