#!/usr/bin/env python3
"""make_deconv_performance.py -- performance plots of the pile-up deconvolution study for the writeup.

Read-only: reads waveform_analysis/pileup_analysis/efficiency.tsv (written by pileup_deconvolution_study.py)
and writes two vector PDFs into waveform_analysis/writeup/figs/:

  deconv_efficiency_vs_dt.pdf   secondary-pulse efficiency eff_ta (timing + amplitude match) vs dt, one line per
                                r = A2/A1, one panel per A1; pooled over the 5 modules (found counts summed,
                                eff_ta x n, divided by the summed n); error bar = binomial sqrt(e (1 - e) / N).
  deconv_main_pull_vs_dt.pdf    main-pulse amplitude pull (A1_fit - A1)/A1 vs dt: top row the median (pull_med),
                                bottom row the width 1.4826 MAD (pull_mad), each the plain mean over the 5 modules
                                of the per-cell TSV values over the module cells with >= 100 found events (a point
                                is left out when no module cell qualifies); dashed grey line = the same for the
                                r = 0 control (its pull uses the one-pulse control events).

    python3 waveform_analysis/writeup/make_figs/make_deconv_performance.py
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WA = os.path.dirname(os.path.dirname(HERE))                       # waveform_analysis/
TSV = os.path.join(WA, "pileup_analysis", "efficiency.tsv")
OUT = os.path.join(os.path.dirname(HERE), "figs")

# Okabe-Ito hues (black replaced by a dark ochre); markers give a second, colour-free encoding
COLORS = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7", "#8C6D1F"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X"]
MIN_FOUND = 100       # pull points need at least this many found events in a (module, A1, r, dt) cell

plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5, "legend.fontsize": 7,
                     "xtick.labelsize": 6.5, "ytick.labelsize": 7, "axes.linewidth": 0.6,
                     "lines.linewidth": 1.0, "lines.markersize": 3.2, "pdf.fonttype": 42})


def load():
    with open(TSV) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    for r in rows:
        for k in r:
            if k != "module":
                r[k] = float(r[k])
    return rows


def cell(rows, a1, ratio, dt):
    return [r for r in rows if r["a1"] == a1 and r["ratio"] == ratio and r["dt_ns"] == dt]


def style(ax, dts):
    ax.set_xticks(range(len(dts)), [f"{d:g}" for d in dts])
    ax.grid(axis="y", color="0.88", lw=0.5)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def efficiency(rows, ratios, dts, a1s):
    fig, axs = plt.subplots(len(a1s), 1, figsize=(6.5, 7.2), sharex=True)
    for ax, a1 in zip(axs, a1s):
        for i, rt in enumerate(ratios):
            x, e, s = [], [], []
            for j, d in enumerate(dts):
                c = cell(rows, a1, rt, d)
                if not c:                                 # combination dropped by the generator (no clipping)
                    continue
                n = sum(r["n"] for r in c)
                k = sum(r["eff_ta"] * r["n"] for r in c)
                x.append(j + 0.06 * (i - 3))
                e.append(100 * k / n)
                s.append(100 * np.sqrt(k / n * (1 - k / n) / n))
            ax.errorbar(x, e, s, color=COLORS[i], marker=MARKERS[i], lw=1.0, capsize=0, label=f"r = {rt:g}")
        ax.set_ylim(-3, 103)
        ax.set_ylabel("efficiency (%)")
        ax.set_title(f"A1 = {a1:g} ADC (pooled over 5 modules, 50 000 events per point)", loc="left")
        style(ax, dts)
    axs[0].legend(ncol=7, loc="lower center", bbox_to_anchor=(0.5, 1.12), frameon=False, handletextpad=0.3,
                  columnspacing=0.9)
    axs[-1].set_xlabel("time gap dt = t2 - t1 (ns), grid points (not to scale)")
    fig.tight_layout(h_pad=0.6)
    fig.savefig(os.path.join(OUT, "deconv_efficiency_vs_dt.pdf"))
    plt.close(fig)


def pull(rows, ratios, dts, a1s):
    fig, axs = plt.subplots(2, len(a1s), figsize=(7.0, 4.6), sharex=True)
    for c, a1 in enumerate(a1s):
        ctl = cell(rows, a1, 0.0, 0.0)
        for row, key in enumerate(("pull_med", "pull_mad")):
            ax = axs[row, c]
            ax.axhline(100 * np.mean([r[key] for r in ctl]), color="0.45", ls="--", lw=0.8,
                       label="r = 0 control" if (row, c) == (0, 0) else None)
            for i, rt in enumerate(ratios):
                x, y = [], []
                for j, d in enumerate(dts):                # cells with < MIN_FOUND found events are left out
                    v = [r[key] for r in cell(rows, a1, rt, d) if np.isfinite(r[key]) and r["eff_ta"] * r["n"] >= MIN_FOUND]
                    x.append(j)
                    y.append(100 * np.mean(v) if v else np.nan)   # NaN breaks the line at a missing cell
                ax.plot(x, y, color=COLORS[i], marker=MARKERS[i], ms=2.4, lw=0.8,
                        label=f"r = {rt:g}" if (row, c) == (0, 0) else None)
            style(ax, dts)
            ax.tick_params(axis="x", labelrotation=90, labelsize=5.5)
            if row == 0:
                ax.set_title(f"A1 = {a1:g} ADC")
            else:
                ax.set_xlabel("dt (ns), grid points")
                ax.set_ylim(bottom=0)
        axs[0, 0].set_ylabel("median pull (%)")
        axs[1, 0].set_ylabel("pull width 1.4826 MAD (%)")
    fig.legend(ncol=8, loc="upper center", frameon=False, handletextpad=0.3, columnspacing=0.8)
    fig.tight_layout(rect=(0, 0, 1, 0.93), w_pad=0.8, h_pad=0.5)
    fig.savefig(os.path.join(OUT, "deconv_main_pull_vs_dt.pdf"))
    plt.close(fig)


def main():
    rows = load()
    ratios = sorted({r["ratio"] for r in rows if r["ratio"] > 0})
    dts = sorted({r["dt_ns"] for r in rows if r["ratio"] > 0})
    a1s = sorted({r["a1"] for r in rows})
    os.makedirs(OUT, exist_ok=True)
    efficiency(rows, ratios, dts, a1s)
    pull(rows, ratios, dts, a1s)
    print("wrote", os.path.join(OUT, "deconv_efficiency_vs_dt.pdf"), os.path.join(OUT, "deconv_main_pull_vs_dt.pdf"))


if __name__ == "__main__":
    main()
