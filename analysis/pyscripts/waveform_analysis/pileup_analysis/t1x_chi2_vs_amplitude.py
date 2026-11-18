#!/usr/bin/env python3
"""t1x_chi2_vs_amplitude.py -- T1 extension: fraction of single-pulse fits with chi2/dof > 20 (and > 5) and the median chi2/dof vs fitted amplitude, per run type.

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/t1x_chi2_vs_amplitude.py [--scan-dir <dir>] [--out <dir>]

No refit: reads the stored scan rows of output/pileup_rates/scan/run_*.npz (full scans, all splits, default fit): ch_k == 1, ch_a[:, 0] = A, ch_chi2 = chi2/dof,
ch_nclip = number of overflow samples (4095).  Class from mod_pwo4 (PbWO4 / lead glass).  Rows with NaN chi2 are counted in the table but are not in any fraction.
Two sets per curve: all k == 1 rows, and rows without overflow samples (nclip == 0).
-> T1x_chi2_vs_amplitude_20261116.pdf (curves, ratio page, count tables) and .tsv
"""
from __future__ import annotations

import sys
import argparse
import textwrap
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
import numpy as np                                      # noqa: E402
import matplotlib                                       # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                         # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages    # noqa: E402

EDGES = np.array([50, 100, 200, 400, 800, 1200, 1600, 2000, 2500, 3000, 3500, 4095.0001])
LABELS = [f"{int(a)}-{int(b) if b < 4096 else 4095}" for a, b in zip(EDGES[:-1], EDGES[1:])]
RTS = ["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17"]
COL = dict(zip(RTS, ["tab:green", "tab:red", "tab:blue", "tab:orange"]))
CLASSES = ["PbWO4", "LG"]
NB = len(LABELS)
STAMP = "20261116"


def load(scan_dir):
    """cnt[rt][class][set][metric] arrays [NB] and median chi2 lists; sets 0 = all, 1 = nclip == 0.  Chi2 values kept per bin for the median."""
    res = {rt: {c: {s: dict(n=np.zeros(NB, np.int64), nnan=np.zeros(NB, np.int64), g20=np.zeros(NB, np.int64), g5=np.zeros(NB, np.int64), chi=[[] for _ in range(NB)])
                    for s in (0, 1)} for c in range(2)} for rt in RTS}
    runs = {rt: [] for rt in RTS}
    over = {rt: [0, 0] for rt in RTS}
    for fn in sorted(Path(scan_dir).glob("run_*.npz")):
        z = np.load(fn)
        rt = str(z["run_type"])
        runs[rt].append((int(z["run"]), len(z["splits"])))
        k = z["ch_k"]
        sel = np.nonzero(k == 1)[0]
        a = z["ch_a"][sel, 0].astype(np.float64)
        chi = z["ch_chi2"][sel].astype(np.float64)
        ncl = z["ch_nclip"][sel]
        cls = np.where(z["mod_pwo4"][z["ch_midx"][sel]], 0, 1)
        ib = np.digitize(a, EDGES) - 1
        okb = (ib >= 0) & (ib < NB)
        for c in range(2):
            over[rt][c] += int(((a >= EDGES[-1]) & (cls == c)).sum())
        for c in range(2):
            for s in (0, 1):
                m = okb & (cls == c) & ((ncl == 0) if s else True)
                R = res[rt][c][s]
                nan = np.isnan(chi)
                R["n"] += np.bincount(ib[m & ~nan], minlength=NB)
                R["nnan"] += np.bincount(ib[m & nan], minlength=NB)
                R["g20"] += np.bincount(ib[m & ~nan & (chi > 20)], minlength=NB)
                R["g5"] += np.bincount(ib[m & ~nan & (chi > 5)], minlength=NB)
                for b in range(NB):
                    mb = m & ~nan & (ib == b)
                    if mb.any():
                        R["chi"][b].append(chi[mb].astype(np.float32))
        print(f"[load] {fn.name} {rt}: {sel.size} k==1 rows", flush=True)
    for rt in RTS:
        for c in range(2):
            for s in (0, 1):
                R = res[rt][c][s]
                R["med"] = np.array([np.median(np.concatenate(v)) if v else np.nan for v in R["chi"]])
                R["q75"] = np.array([np.percentile(np.concatenate(v), 75) if v else np.nan for v in R["chi"]])
                del R["chi"]
    res["over"] = over
    return res, runs


NMIN = 100        # rows per bin needed to draw a point or enter a ratio


def frac(R, key, nmin=NMIN):
    n = R["n"].astype(float)
    f = np.where(n >= nmin, R[key] / np.maximum(n, 1), np.nan)
    e = np.sqrt(np.maximum(f * (1 - f), 0) / np.maximum(n, 1))
    return f, e


def verdict(res):
    """Text lines: ratio of the 2p24 H2 fraction to the pooled other three at equal A (PbWO4, nclip == 0), per bin."""
    L = []
    for key, nm in (("g20", "chi2/dof > 20"), ("g5", "chi2/dof > 5")):
        L.append(f"PbWO4, nclip == 0, fraction with {nm}: 2p24 H2 vs the other three (pooled) at equal A")
        for b in range(NB):
            o = [res[rt][0][1] for rt in RTS if rt != "2p24GeV_H2"]
            no, ko = sum(R["n"][b] for R in o), sum(R[key][b] for R in o)
            h = res["2p24GeV_H2"][0][1]
            if h["n"][b] < NMIN or no < NMIN:
                L.append(f"   A {LABELS[b]:>10}: n(2p24 H2) {int(h['n'][b])}, n(others) {int(no)}: too few")
                continue
            fh, fo = h[key][b] / h["n"][b], ko / no
            se = np.sqrt(fh * (1 - fh) / h["n"][b] + fo * (1 - fo) / no)
            spread = [R[key][b] / R["n"][b] for R in o if R["n"][b] >= NMIN]
            L.append(f"   A {LABELS[b]:>10}: 2p24 H2 {100 * fh:6.2f} %  others {100 * fo:6.2f} % (per type {', '.join(f'{100 * s:.2f}' for s in spread)})  ratio {fh / fo if fo > 0 else np.nan:5.2f}  diff {100 * (fh - fo):+.2f} +- {100 * se:.2f} %")
    L.append("")
    L.append("PbWO4, nclip == 0: 2p24 H2 / 3p49 H2 (same target, different beam energy), fraction > 20 (n 2p24 / n 3p49):")
    for b in range(NB):
        h, t = res["2p24GeV_H2"][0][1], res["3p49GeV_H2"][0][1]
        if h["n"][b] >= NMIN and t["n"][b] >= NMIN and t["g20"][b] > 0:
            L.append(f"   A {LABELS[b]:>10}: {100 * h['g20'][b] / h['n'][b]:6.3f} % / {100 * t['g20'][b] / t['n'][b]:6.3f} % = {(h['g20'][b] / h['n'][b]) / (t['g20'][b] / t['n'][b]):6.2f}  ({int(h['n'][b])} / {int(t['n'][b])})")
    return L


def page_curves(pdf, res, c, runs):
    fig, axs = plt.subplots(3, 1, figsize=(11, 12.5), sharex=True)
    x = np.arange(NB)
    fig.suptitle(f"{CLASSES[c]}: single-pulse (k == 1) fit quality vs fitted amplitude, default fit, full scans", fontsize=13, fontweight="bold")
    for ax, key, yl in ((axs[0], "g20", "fraction with chi2/dof > 20"), (axs[1], "g5", "fraction with chi2/dof > 5"), (axs[2], "med", "median chi2/dof")):
        for rt in RTS:
            for s, ls, lab in ((1, "-", "no overflow sample"), (0, "--", "all (incl. nclip > 0)")):
                R = res[rt][c][s]
                if key == "med":
                    y, e = np.where(R["n"] >= NMIN, R["med"], np.nan), None
                else:
                    y, e = frac(R, key)
                ax.errorbar(x + (0.06 if s == 0 else -0.06) * (RTS.index(rt) - 1.5), y, yerr=e, color=COL[rt], ls=ls, marker="o" if s else "s", ms=4 if s else 3, lw=1.4 if s else 0.9,
                            mfc=COL[rt] if s else "none", capsize=2, label=f"{rt}, {lab}")
        ax.set_ylabel(yl)
        ax.grid(alpha=0.3)
        ax.set_yscale("log")
        if key == "med":
            ax.axhline(1, color="0.5", lw=0.7)
    axs[0].legend(fontsize=7, ncol=2)
    axs[2].set_xticks(x)
    axs[2].set_xticklabels(LABELS, rotation=30, fontsize=8)
    axs[2].set_xlabel("fitted amplitude A (ADC); bins with < 100 rows are not drawn; error bars: binomial")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    pdf.savefig(fig)
    plt.close(fig)


def page_text(pdf, runs, res, L):
    fig = plt.figure(figsize=(16, 9.5))
    fig.text(0.03, 0.96, "T1x: chi2/dof of single-pulse fits vs amplitude per run type (no refit)", fontsize=16, fontweight="bold")
    rl = "\n".join(f"{rt}: runs " + ", ".join(f"{r} ({n} splits)" for r, n in sorted(runs[rt])) for rt in RTS)
    txt = ("Source: output/pileup_rates/scan/run_*.npz (full scans, all splits, default fit); k == 1 rows, A = ch_a[:, 0], chi2/dof = ch_chi2, overflow = ch_nclip > 0.\n" + rl +
           "\nk == 1 rows with A > 4095 (outside all bins, almost all with overflow samples; PbWO4 / LG): " + "; ".join(f"{rt} {res['over'][rt][0]} / {res['over'][rt][1]}" for rt in RTS) + "\nCommand: " + " ".join(sys.argv) + "\nPage 1-2: curves (solid = nclip == 0, dashed = all).  Page 3: 2p24 H2 minus others.  Then count tables.\n\n" + "\n".join(L))
    fig.text(0.02, 0.92, txt, fontsize=6.6, va="top", family="monospace", linespacing=1.2)
    pdf.savefig(fig)
    plt.close(fig)


def page_ratio(pdf, res):
    fig, axs = plt.subplots(1, 3, figsize=(16, 5.5))
    x = np.arange(NB)
    fig.suptitle("PbWO4, nclip == 0: each run type divided by the pooled other three at equal A (1 = same curve; 2p24 H2 above 1 = worse than the rest)", fontsize=11, fontweight="bold")
    for ax, key, yl in zip(axs, ("g20", "g5", "med"), ("frac chi2>20 ratio", "frac chi2>5 ratio", "median chi2/dof ratio")):
        for rt in RTS:
            y = []
            for b in range(NB):
                o = [res[r][0][1] for r in RTS if r != rt]
                R = res[rt][0][1]
                if key == "med":
                    ov = [Ro["med"][b] for Ro in o if Ro["n"][b] >= NMIN]
                    y.append(R["med"][b] / np.mean(ov) if R["n"][b] >= NMIN and ov else np.nan)
                else:
                    no, ko = sum(Ro["n"][b] for Ro in o), sum(Ro[key][b] for Ro in o)
                    y.append((R[key][b] / R["n"][b]) / (ko / no) if R["n"][b] >= NMIN and no >= NMIN and ko > 0 and R[key][b] > 0 else np.nan)
            ax.plot(x, y, "o-", color=COL[rt], label=rt)
        ax.axhline(1, color="k", lw=0.7)
        ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels(LABELS, rotation=40, fontsize=7)
        ax.set_ylabel(yl)
        ax.grid(alpha=0.3)
    axs[0].legend(fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    pdf.savefig(fig)
    plt.close(fig)


def page_table(pdf, res, c, s):
    fig = plt.figure(figsize=(16, 9.5))
    fig.text(0.03, 0.96, f"Counts per amplitude bin: {CLASSES[c]}, {'no overflow sample (nclip == 0)' if s else 'all k == 1 rows'}", fontsize=14, fontweight="bold")
    hdr = f"{'A bin':>10} | " + " | ".join(f"{rt:^34}" for rt in RTS)
    sub = f"{'':>10} | " + " | ".join(f"{'n':>9} {'>20 %':>7} {'>5 %':>7} {'med':>6}" for _ in RTS)
    lines = [hdr, sub, "-" * len(hdr)]
    for b in range(NB):
        cells = []
        for rt in RTS:
            R = res[rt][c][s]
            n = int(R["n"][b])
            cells.append(f"{n:>9d} {100 * R['g20'][b] / n if n else np.nan:7.2f} {100 * R['g5'][b] / n if n else np.nan:7.2f} {R['med'][b]:6.2f}")
        lines.append(f"{LABELS[b]:>10} | " + " | ".join(cells))
    lines.append("")
    lines.append("g20 / g5 counts (rows with chi2/dof > 20 / > 5) and NaN-chi2 rows are in the .tsv.")
    fig.text(0.02, 0.92, "\n".join(lines), fontsize=10, va="top", family="monospace", linespacing=1.5)
    pdf.savefig(fig)
    plt.close(fig)


def write_tsv(path, res):
    with open(path, "w") as f:
        f.write("run_type\tclass\toverflow_set\tA_lo\tA_hi\tn\tn_nan_chi2\tn_gt20\tn_gt5\tfrac_gt20\tfrac_gt5\tmedian_chi2\tq75_chi2\n")
        for rt in RTS:
            for c in range(2):
                for s, sn in ((0, "all"), (1, "nclip0")):
                    R = res[rt][c][s]
                    for b in range(NB):
                        n = int(R["n"][b])
                        f.write(f"{rt}\t{CLASSES[c]}\t{sn}\t{EDGES[b]:.0f}\t{min(EDGES[b + 1], 4095):.0f}\t{n}\t{int(R['nnan'][b])}\t{int(R['g20'][b])}\t{int(R['g5'][b])}\t"
                                f"{R['g20'][b] / n if n else np.nan:.6g}\t{R['g5'][b] / n if n else np.nan:.6g}\t{R['med'][b]:.4g}\t{R['q75'][b]:.4g}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scan-dir", default=str(HERE.parent.parent / "output" / "pileup_rates" / "scan"))
    ap.add_argument("--out", default=str(HERE))
    a = ap.parse_args()
    res, runs = load(a.scan_dir)
    L = verdict(res)
    print("\n".join(L))
    stem = Path(a.out) / f"T1x_chi2_vs_amplitude_{STAMP}"
    write_tsv(str(stem) + ".tsv", res)
    with PdfPages(str(stem) + ".pdf") as pdf:
        page_text(pdf, runs, res, L)
        for c in range(2):
            page_curves(pdf, res, c, runs)
        page_ratio(pdf, res)
        for c in range(2):
            for s in (1, 0):
                page_table(pdf, res, c, s)
    print(f"[done] {stem}.pdf / .tsv")


if __name__ == "__main__":
    main()
