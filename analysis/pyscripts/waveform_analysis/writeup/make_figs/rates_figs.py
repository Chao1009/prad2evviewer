#!/usr/bin/env python3
"""rates_figs.py -- two summary figures for the pile-up rate section, read only from existing outputs.

Reads waveform_analysis/pileup_analysis/pileup_rates_out/{pileup_modules_<RT>.tsv, pileup_rates.tsv, pileup_slices.tsv}
and $PRAD2_DATABASE_DIR/hycal_map.json (default ~/work/PRad/prad2evviewer/database); writes into ../figs/:
  rates_map_h2_vs_ta.pdf  per-module channel cut rate (rate_cut) for 2p24GeV_H2 and 2p24GeV_Ta_X17, one log colour scale
  rates_summary.pdf       (a) channel cut, cluster cut (all / in-window extras) and raw k >= 2 rates per run type;
                          (b) pooled cut rate per E_main / E_beam slice (pileup_slices.tsv)
"""
import csv
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "figs"
SRC = HERE.parent.parent / "pileup_analysis" / "pileup_rates_out"
DB = Path(os.environ.get("PRAD2_DATABASE_DIR", str(Path.home() / "work/PRad/prad2evviewer/database")))
RTS = ["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17"]
LABEL = {"0p73GeV_H2": "0.73 GeV H$_2$", "2p24GeV_H2": "2.24 GeV H$_2$", "3p49GeV_H2": "3.49 GeV H$_2$",
         "2p24GeV_Ta_X17": "2.24 GeV Ta (X17)"}
OKABE = {"0p73GeV_H2": "#E69F00", "2p24GeV_H2": "#56B4E9", "3p49GeV_H2": "#009E73", "2p24GeV_Ta_X17": "#D55E00"}
MARK = {"0p73GeV_H2": "o", "2p24GeV_H2": "s", "3p49GeV_H2": "^", "2p24GeV_Ta_X17": "D"}
plt.rcParams.update({"font.size": 10, "axes.labelsize": 10, "legend.fontsize": 8.5, "axes.linewidth": 0.8,
                     "axes.spines.top": False, "axes.spines.right": False})


def tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def rate_map():
    geo = {m["n"]: m["geo"] for m in json.load(open(DB / "hycal_map.json"))}
    vmin, vmax = 0.05, 50.0                                    # % ; zero and smaller rates drawn at vmin
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.7), layout="constrained")
    for ax, rt in zip(axs, ("2p24GeV_H2", "2p24GeV_Ta_X17")):
        rows = [r for r in tsv(SRC / f"pileup_modules_{rt}.tsv") if r["is_pbwo4"] == "1"]
        rect = lambda r: Rectangle((geo[r["module"]]["x"] - geo[r["module"]]["sx"] / 2, geo[r["module"]]["y"] - geo[r["module"]]["sy"] / 2),
                                   geo[r["module"]]["sx"], geo[r["module"]]["sy"])   # noqa: E731
        has = [r for r in rows if r["rate_cut"] != "" and r["dead"] == "0"]
        few = [r for r in rows if r["rate_cut"] == "" and r["dead"] == "0"]
        dead = [r for r in rows if r["dead"] == "1"]
        ax.add_collection(PatchCollection([rect(r) for r in few], facecolor="0.85", edgecolor="0.7", lw=0.15))
        ax.add_collection(PatchCollection([rect(r) for r in dead], facecolor="white", edgecolor="0.3", lw=0.3, hatch="////"))
        pc = PatchCollection([rect(r) for r in has], cmap="viridis", edgecolor="0.35", lw=0.1,
                             norm=LogNorm(vmin, vmax))
        pc.set_array(np.clip([100 * float(r["rate_cut"]) for r in has], vmin, vmax))
        ax.add_collection(pc)
        ax.set_xlim(-360, 360)
        ax.set_ylim(-360, 360)
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)")
        ax.set_title(f"{LABEL[rt]}: {len(has)} modules", fontsize=10)
    axs[0].set_ylabel("y (mm)")
    cb = fig.colorbar(pc, ax=axs, shrink=0.85, extend="both")
    cb.set_label("channel cut rate per module (%)")
    fig.savefig(OUT / "rates_map_h2_vs_ta.pdf")
    plt.close(fig)


def summary():
    s = {r["run_type"]: r for r in tsv(SRC / "pileup_rates.tsv")}
    sl = tsv(SRC / "pileup_slices.tsv")
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.5), layout="constrained")
    qs = [("ch_cut", "channel\ncut"), ("cl_cut", "cluster cut\n(all)"), ("cl_cutin", "cluster cut\n(in-window)"),
          ("ch", "channel raw\n$k \\geq 2$")]
    x = np.arange(len(qs))
    for j, rt in enumerate(RTS):
        v = np.array([100 * float(s[rt]["rate_" + q]) for q, _ in qs])
        e = np.array([100 * float(s[rt]["err_" + q]) for q, _ in qs])
        a.errorbar(x + (j - 1.5) * 0.15, v, yerr=e, fmt=MARK[rt], color=OKABE[rt], ms=6, lw=1, capsize=2, label=LABEL[rt])
    a.set_yscale("log")
    a.set_xticks(x, [lab for _, lab in qs], fontsize=8.5)
    a.set_ylabel("rate (%)")
    a.set_title("(a) rates per run type", fontsize=10)
    a.grid(axis="y", which="major", color="0.88", lw=0.6)
    for rt in RTS:
        rr = [r for r in sl if r["run_type"] == rt]
        c = np.array([(float(r["slice_lo"]) + float(r["slice_hi"])) / 2 for r in rr])
        n = np.array([int(r["n_cut"]) for r in rr])
        N = np.array([int(r["n"]) for r in rr])
        b.plot(c, 100 * n / N, MARK[rt] + "-", color=OKABE[rt], ms=6, lw=1.5)   # binomial errors: smaller than the markers except at n_cut <= 12
    b.set_yscale("log")
    b.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0, 1.2], ["0.1-0.3", "0.3-0.5", "0.5-0.7", "0.7-0.9", "0.9-1.1", "1.1-1.3"], fontsize=8.5)
    b.set_xlabel("$E_{main}/E_{beam}$ slice")
    b.set_ylabel("N(cut in slice) / N(kept) (%)")
    b.set_title("(b) pooled cut rate per energy slice", fontsize=10)
    b.grid(axis="y", which="major", color="0.88", lw=0.6)
    fig.legend(*a.get_legend_handles_labels(), loc="outside upper center", ncol=4, frameon=False)
    fig.savefig(OUT / "rates_summary.pdf")
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    rate_map()
    summary()
    print("wrote", OUT / "rates_map_h2_vs_ta.pdf", OUT / "rates_summary.pdf")
