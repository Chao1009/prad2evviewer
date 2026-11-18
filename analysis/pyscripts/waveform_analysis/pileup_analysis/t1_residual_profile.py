#!/usr/bin/env python3
"""t1_residual_profile.py -- T1: fit residual r_j / A of single-pulse waveforms (S1) vs tau = t_j - t0, measured and assumed sigma(t).

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/t1_residual_profile.py [--n-events N] [--n-events-pwo N] [--cache-dir <dir>] [--out <dir>] [--rebuild]

S1, the shared pass and the cache: resid_common.py (default fit; k == 1, no overflow, A 50-3000 ADC; PbWO4 W* vs lead glass G*).  Classes x amplitude bins
(50-200-800-3000 ADC) x phase bins (t0 mod 4 ns, 1 ns each).  tau bins 4 ns, -20 .. +200 ns (bin centres), all 100 samples of the waveform, also outside the fit window.
Per bin: mean of r/A, sem = std / sqrt(n), rms = std of r/A about the mean, assumed sigma = sqrt(mean(V / A^2)) with V the MF.variance used in the chi2
(sig_ref^2 + 1/12 + template-shape term), n samples.  Per-module spread: std over modules (>= 200 waveforms in the amp bin) of the per-module mean profile.
-> T1_residual_profile_20261116.pdf / .tsv (pages: summary, profiles per class and amp bin, phase overlays, per-module spread, 20 example waveforms).
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
NPH = 4
MIN_MOD = 200
WIN = (-12.0, 116.0)       # nominal fit window in tau: 32 samples from 3 samples before the onset (varies with the phase by < 4 ns)
SEED = 20261116


def accumulate(d):
    """Sums per (class, amp bin, phase bin, tau bin) and per (module, amp bin, tau bin): n, sum x, sum x^2, sum V / A^2 (x = r / A)."""
    nm = len(d["names"])
    G = {k: np.zeros(2 * 3 * NPH * NT) for k in ("n", "sx", "sxx", "sv")}
    M = {k: np.zeros(nm * 3 * NT) for k in ("n", "sx", "sxx")}
    nwf = np.bincount(d["mi"].astype(np.int64) * 3 + d["ab"], minlength=nm * 3)
    ph = np.clip(np.floor(np.mod(d["t0"], 4.0)).astype(np.int64), 0, NPH - 1)
    N = d["a"].size
    for i0 in range(0, N, 40000):
        s = slice(i0, i0 + 40000)
        A = d["a"][s].astype(np.float64)[:, None]
        x = d["resid"][s] / A
        v = d["var"][s] / A ** 2
        tau = rc.TT[None, :] - d["t0"][s][:, None]
        ib = np.floor((tau + 22.0) / 4.0).astype(np.int64)
        ok = (ib >= 0) & (ib < NT)
        g = (((d["cls"][s].astype(np.int64) * 3 + d["ab"][s]) * NPH + ph[s])[:, None] * NT + ib)[ok]
        m = ((d["mi"][s].astype(np.int64) * 3 + d["ab"][s])[:, None] * NT + ib)[ok]
        xs, vs = x[ok], v[ok]
        for key, w in (("n", None), ("sx", xs), ("sxx", xs * xs), ("sv", vs)):
            G[key] += np.bincount(g, weights=w, minlength=G[key].size)
        for key, w in (("n", None), ("sx", xs), ("sxx", xs * xs)):
            M[key] += np.bincount(m, weights=w, minlength=M[key].size)
    G = {k: v.reshape(2, 3, NPH, NT) for k, v in G.items()}
    M = {k: v.reshape(nm, 3, NT) for k, v in M.items()}
    return G, M, nwf.reshape(nm, 3)


def stats(n, sx, sxx, sv):
    """mean, sem, rms (std about the mean), assumed sigma from the sums (arrays over the last axis = tau)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = sx / n
        std = np.sqrt(np.maximum(sxx / n - mean ** 2, 0.0))
        return mean, std / np.sqrt(n), std, np.sqrt(sv / n)


def runs_where(flag, minlen=3):
    """[(i0, i1)] of the contiguous True stretches of flag of at least minlen entries."""
    out, i = [], 0
    while i < len(flag):
        if flag[i]:
            j = i
            while j + 1 < len(flag) and flag[j + 1]:
                j += 1
            if j - i + 1 >= minlen:
                out.append((i, j))
            i = j + 1
        else:
            i += 1
    return out


def pooled(G):
    """Phase-pooled sums [2, 3, NT]."""
    return {k: v.sum(axis=2) for k, v in G.items()}


def decisions(P, mod_spread):
    """Text lines: coherent mean stretches and sigma ratio stretches per class / amp bin (phase pooled)."""
    mean, sem, rms, asg = stats(P["n"], P["sx"], P["sxx"], P["sv"])
    lines = []
    for c, cn in enumerate(rc.CLASSES):
        for a, al in enumerate(rc.AMP_LABELS):
            tag = f"{cn} A {al}"
            n = P["n"][c, a]
            ok = n > 100
            coh = runs_where(ok & (np.abs(mean[c, a]) > 2 * sem[c, a]))
            sig = "; ".join(f"tau {TAU[i]:.0f}..{TAU[j]:.0f} ns (peak {100 * mean[c, a, i:j + 1][np.argmax(np.abs(mean[c, a, i:j + 1]))]:+.2f} % A at {TAU[i + int(np.argmax(np.abs(mean[c, a, i:j + 1])))]:.0f} ns)" for i, j in coh)
            big = runs_where(ok & (np.abs(mean[c, a]) > 2 * sem[c, a]) & (np.abs(mean[c, a]) > 0.005))
            bigs = "; ".join(f"tau {TAU[i]:.0f}..{TAU[j]:.0f} ns (peak {100 * mean[c, a, i:j + 1][np.argmax(np.abs(mean[c, a, i:j + 1]))]:+.2f} % A)" for i, j in big)
            ratio = rms[c, a] / asg[c, a]
            wide = runs_where(ok & (ratio > 1.1))
            narrow = runs_where(ok & (ratio < 0.9))
            inwin = (TAU >= WIN[0]) & (TAU <= WIN[1])
            fm = lambda r: "; ".join(f"tau {TAU[i]:.0f}..{TAU[j]:.0f} ns (ratio {ratio[i:j + 1].min():.2f}-{ratio[i:j + 1].max():.2f}, max at {TAU[i + int(np.argmax(ratio[i:j + 1]))]:.0f})" for i, j in r) or "none"    # noqa: E731
            lines.append(f"{tag}: n(waveforms) ~ {int(n[np.argmax(n)] if n.max() > 0 else 0)}")
            lines.append(f"   coherent mean (|mean| > 2 sem, >= 3 bins): {sig or 'none'}")
            lines.append(f"   ... and |mean| > 0.5 % A: {bigs or 'none'}")
            lines.append(f"   rms / assumed sigma: in fit window median {np.nanmedian(ratio[inwin & ok]):.2f}, outside median {np.nanmedian(ratio[~inwin & ok]):.2f}; wider (> 1.1): {fm(wide)}; narrower (< 0.9): {fm(narrow)}")
    lines.append("")
    lines.append("Per-module spread (std over modules of the per-module mean profile, class mean subtracted; rms about zero also given) vs the pure-noise expectation sqrt(mean sem_m^2):")
    for (c, a), (sp, ex, nmod, r0) in mod_spread.items():
        inwin = (TAU >= WIN[0]) & (TAU <= WIN[1])
        lines.append(f"   {rc.CLASSES[c]} A {rc.AMP_LABELS[a]}: {nmod} modules; std over modules (class mean subtracted) / expected in window median {np.nanmedian((sp / ex)[inwin]):.2f} (max {np.nanmax((sp / ex)[inwin]):.2f}); "
                     f"std in window {100 * np.nanmean(sp[inwin]):.2f} % A, expected {100 * np.nanmean(ex[inwin]):.2f} % A;  rms about zero {100 * np.nanmean(r0[inwin]):.2f} % A (ratio to expected median {np.nanmedian((r0 / ex)[inwin]):.2f})")
    return lines


def module_spread(d, M, nwf):
    """{(class, amp bin): (spread [NT], expected noise-only [NT], n modules, rms about zero [NT])} for modules with >= MIN_MOD waveforms.
    spread = std over modules of the per-module mean (class mean subtracted); rms about zero = sqrt(mean over modules of mean^2)."""
    out, curves = {}, {}
    cls = np.where(d["pwo4"], 0, 1)
    for c in range(2):
        for a in range(3):
            sel = np.nonzero((cls == c) & (nwf[:, a] >= MIN_MOD))[0]
            if sel.size < 2:
                continue
            mean, sem, _, _ = stats(M["n"][sel, a], M["sx"][sel, a], M["sxx"][sel, a], M["sx"][sel, a] * 0)
            out[(c, a)] = (np.nanstd(mean, axis=0), np.sqrt(np.nanmean(sem ** 2, axis=0)), sel.size, np.sqrt(np.nanmean(mean ** 2, axis=0)))
            curves[(c, a)] = mean
    return out, curves


def write_tsv(path, G, P):
    with open(path, "w") as f:
        f.write("class\tamp_bin\tphase_bin\ttau\tmean\tsem\trms\tassumed_sigma\tn\n")
        for c, cn in enumerate(rc.CLASSES):
            for a, al in enumerate(rc.AMP_LABELS):
                for p in list(range(NPH)) + ["all"]:
                    S = P if p == "all" else G
                    sl = (c, a) if p == "all" else (c, a, p)
                    mean, sem, rms, asg = stats(S["n"][sl], S["sx"][sl], S["sxx"][sl], S["sv"][sl])
                    for i in range(NT):
                        n = S["n"][sl][i]
                        f.write(f"{cn}\t{al}\t{p}\t{TAU[i]:.0f}\t{mean[i]:.6g}\t{sem[i]:.3g}\t{rms[i]:.6g}\t{asg[i]:.6g}\t{int(n)}\n")


def page_summary(pdf, d, lines, a):
    fig = plt.figure(figsize=(16, 9))
    fig.text(0.03, 0.96, "T1: fit residual profile of single-pulse waveforms", fontsize=18, fontweight="bold")
    txt = (rc.run_runs_note(d) + "\n" + "\n".join(rc.counts_text(d)) + "\n"
           "Pass: default fit (MONO_EDGE, NOISE_FLAG off, SHAPE_SCALE 1).  S1 = fit k == 1, no sample at 4095, fitted A 50-3000 ADC.  tau = t_j - t0, 4 ns bins -20..200 ns; x = (data - ped - fit) / A.\n"
           "mean +- sem (sem = std / sqrt n); rms = std of x about its mean; assumed = sqrt(mean V / A^2), V = sigma_ref^2 + 1/12 + template-shape term (MF.variance).  Fit window ~ tau -12..116 ns.\n"
           "Command: " + " ".join(sys.argv) + "\n"
           "Coherent mean: >= 3 contiguous tau bins with |mean| > 2 sem (phase pooled).  Ratio = rms / assumed.\n\n" + "\n".join(lines))
    fig.text(0.02, 0.92, txt, fontsize=6.8, va="top", family="monospace", linespacing=1.25)
    pdf.savefig(fig)
    plt.close(fig)


def page_profile(pdf, c, a, P, d):
    mean, sem, rms, asg = stats(P["n"][c, a], P["sx"][c, a], P["sxx"][c, a], P["sv"][c, a])
    fig, ax = plt.subplots(3, 1, figsize=(11, 10), sharex=True, gridspec_kw=dict(height_ratios=[3, 3, 1.5]))
    fig.suptitle(f"{rc.CLASSES[c]}, fitted A {rc.AMP_LABELS[a]} ADC: S1 residual / A vs tau   (n waveforms {int(P['n'][c, a].max())})", fontsize=13, fontweight="bold")
    for x in ax:
        x.axvspan(*WIN, color="0.93", lw=0)
        x.axvline(0, color="tab:blue", ls=":", lw=0.8)
    ax[0].plot(TAU, 100 * mean, "k.-", lw=0.8, label="mean")
    ax[0].fill_between(TAU, 100 * (mean - 2 * sem), 100 * (mean + 2 * sem), color="tab:red", alpha=0.35, label="+- 2 sem")
    ax[0].axhline(0, color="r", lw=0.7)
    ax[0].set_ylabel("mean of r / A  (% of A)")
    ax[0].legend(fontsize=8)
    ax[1].plot(TAU, 100 * rms, "k.-", lw=0.8, label="measured rms (std about mean)")
    ax[1].plot(TAU, 100 * asg, "r.-", lw=0.8, label="assumed sqrt(V) / A")
    ax[1].set_ylabel("sigma / A  (% of A)")
    ax[1].legend(fontsize=8)
    ax[2].plot(TAU, rms / asg, "k.-", lw=0.8)
    ax[2].axhline(1, color="r", lw=0.7)
    ax[2].set_ylabel("rms / assumed")
    ax[2].set_xlabel("tau = t_j - t0 (ns)    (grey: nominal fit window)")
    ax[2].set_ylim(0.5, 2.0)
    for x in ax:
        x.set_xlim(-22, 202)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    pdf.savefig(fig)
    plt.close(fig)


def page_phase(pdf, G, kind):
    fig, axs = plt.subplots(3, 2, figsize=(15, 10), sharex=True)
    fig.suptitle("Phase overlay (t0 mod 4 ns bins 0-1, 1-2, 2-3, 3-4 ns): " + ("mean of r / A (%)" if kind == "mean" else "measured rms / assumed sigma"), fontsize=13, fontweight="bold")
    for c in range(2):
        for a in range(3):
            ax = axs[a, c]
            for p in range(NPH):
                mean, sem, rms, asg = stats(G["n"][c, a, p], G["sx"][c, a, p], G["sxx"][c, a, p], G["sv"][c, a, p])
                ax.plot(TAU, 100 * mean if kind == "mean" else rms / asg, ".-", lw=0.8, label=f"phase {p}", color=plt.cm.viridis(p / 3.2))
            ax.axvspan(*WIN, color="0.93", lw=0)
            ax.axhline(0 if kind == "mean" else 1, color="r", lw=0.6)
            ax.set_title(f"{rc.CLASSES[c]}  A {rc.AMP_LABELS[a]}", fontsize=10, loc="left")
            if kind != "mean":
                ax.set_ylim(0.5, 2.0)
            if a == 2:
                ax.set_xlabel("tau (ns)")
    axs[0, 0].legend(fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)


def page_spread(pdf, spread, curves):
    fig, axs = plt.subplots(3, 2, figsize=(15, 10), sharex=True)
    fig.suptitle(f"Per-module spread of the mean profile (modules with >= {MIN_MOD} waveforms in the amp bin): thin = per-module means (% A), thick black = std over modules, green = rms about zero, red = noise-only expectation",
                 fontsize=11, fontweight="bold")
    for c in range(2):
        for a in range(3):
            ax = axs[a, c]
            ax.set_title(f"{rc.CLASSES[c]}  A {rc.AMP_LABELS[a]}", fontsize=10, loc="left")
            if (c, a) not in spread:
                ax.text(0.5, 0.5, "fewer than 2 modules", transform=ax.transAxes, ha="center")
                continue
            sp, ex, nmod, r0 = spread[(c, a)]
            cur = curves[(c, a)]
            for k in range(min(cur.shape[0], 150)):
                ax.plot(TAU, 100 * cur[k], color="tab:blue", lw=0.3, alpha=0.35)
            ax.plot(TAU, 100 * r0, "g-", lw=1.2, label="rms about zero")
            ax.plot(TAU, -100 * r0, "g-", lw=1.2)
            ax.plot(TAU, 100 * sp, "k-", lw=1.8, label=f"std over {nmod} modules (class mean subtracted)")
            ax.plot(TAU, -100 * sp, "k-", lw=1.8)
            ax.plot(TAU, 100 * ex, "r--", lw=1.2, label="noise-only expectation")
            ax.plot(TAU, -100 * ex, "r--", lw=1.2)
            ax.axvspan(*WIN, color="0.93", lw=0, zorder=0)
            lim = 3 * np.nanmax(sp[(TAU > -10) & (TAU < 150)]) * 100
            ax.set_ylim(-max(lim, 0.1), max(lim, 0.1))
            if a == 2:
                ax.set_xlabel("tau (ns)")
    axs[0, 0].legend(fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    pdf.savefig(fig)
    plt.close(fig)


def pages_examples(pdf, d, P):
    rng = np.random.default_rng(SEED)
    picks = []
    for c, quota in ((0, (3, 3, 4)), (1, (3, 3, 4))):
        for a in range(3):
            pool = np.nonzero((d["cls"] == c) & (d["ab"] == a))[0]
            picks += [(c, a, int(i)) for i in rng.choice(pool, size=min(quota[a], pool.size), replace=False)]
        short = 10 - sum(1 for p in picks if p[0] == c)           # amp bins with too few waveforms (lead glass): fill from the lowest bin
        if short > 0:
            have = {p[2] for p in picks}
            pool = np.array([i for i in np.nonzero((d["cls"] == c) & (d["ab"] == 0))[0][:5000] if int(i) not in have])
            picks += [(c, 0, int(i)) for i in rng.choice(pool, size=min(short, pool.size), replace=False)]
    PER = 6
    for p0 in range(0, len(picks), PER):
        fig = plt.figure(figsize=(16, 12.5))
        fig.text(0.03, 0.985, va="top", s=textwrap.fill(f"Example S1 waveforms ({p0 + 1}-{min(p0 + PER, len(picks))} of {len(picks)}): random per class / amp bin (seed {SEED}); lower strip = data - fit; "
                              "grey band = measured +- rms of the class/amp bin (x A), red dashed = assumed +- sqrt(V)", 190), fontsize=9.5)
        outer = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.12, left=0.05, right=0.99, top=0.95, bottom=0.04)
        for q, (c, a, i) in enumerate(picks[p0:p0 + PER]):
            gs = outer[q // 2, q % 2].subgridspec(2, 1, height_ratios=[2.2, 1.4], hspace=0.08)
            ax, ar = fig.add_subplot(gs[0]), fig.add_subplot(gs[1], sharex=None)
            A, t0, ped = float(d["a"][i]), float(d["t0"][i]), float(d["ped"][i])
            model = d["Y"][i] - ped - d["resid"][i]
            rc.panel(ax, d["Y"][i], ped, model, t0, f"{rc.CLASSES[c]} {d['modname'][i]}  run {d['run'][i]} ev {d['ev'][i]}  A {A:.0f}  t0 {t0:.1f} ns  chi2/dof {d['chi2'][i]:.2f}")
            ax.tick_params(labelbottom=False)
            _, _, rms, _ = stats(P["n"][c, a], P["sx"][c, a], P["sxx"][c, a], P["sv"][c, a])
            band = np.interp(rc.TT - t0, TAU, rms, left=np.nan, right=np.nan) * A
            sg = np.sqrt(d["var"][i])
            ar.fill_between(rc.TT, -band, band, color="0.8", label="measured +- rms")
            ar.plot(rc.TT, sg, "r--", lw=1, label="assumed +- sqrt(V)")
            ar.plot(rc.TT, -sg, "r--", lw=1)
            ar.plot(rc.TT, d["resid"][i], "k.-", lw=0.6, ms=3)
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
    a = ap.parse_args()
    d = rc.load_s1(a)
    G, M, nwf = accumulate(d)
    P = pooled(G)
    spread, curves = module_spread(d, M, nwf)
    lines = decisions(P, spread)
    print(rc.run_runs_note(d))
    print("\n".join(rc.counts_text(d)))
    print("\n".join(lines))
    stem = Path(a.out) / f"T1_residual_profile_{rc.STAMP}"
    write_tsv(str(stem) + ".tsv", G, P)
    with PdfPages(str(stem) + ".pdf") as pdf:
        page_summary(pdf, d, lines, a)
        for c in range(2):
            for b in range(3):
                page_profile(pdf, c, b, P, d)
        page_phase(pdf, G, "mean")
        page_phase(pdf, G, "ratio")
        page_spread(pdf, spread, curves)
        ne = pages_examples(pdf, d, P)
    print(f"[done] {stem}.pdf / .tsv, {ne} examples")


if __name__ == "__main__":
    main()
