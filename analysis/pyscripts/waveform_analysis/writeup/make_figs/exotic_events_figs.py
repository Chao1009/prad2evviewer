#!/usr/bin/env python3
"""Figures of the exotic-waveform chapter (Section "Exotic waveforms" of the write-up).

Reads  (read-only) ../../pileup_analysis/exotic_events/
         exotic_P1_quantities_20261117.pdf            page 1, the 12 user examples (crops)
         exotic_P4_module_rates_20261117.tsv          per-module exclusive rates (maps)
         exotic_P6_amp_vs_w50_20261118.pdf            page 1, fitted A_max vs w50 (crop)
         exotic_P7_saturated_seeds_20261118.tsv       saturated seeds in the elastic window (scatter)
         exotic_P8_energy_error_20261118.tsv          energy-scale shift per ring (table scale)
         exotic_P12_rates_vs_energy_20261118.tsv      rates per energy slice (tables chan, cluster)
         exotic_P13_drop_flagged_20261118.tsv         energy removed per sample under the drop policy (table p8_scale)
       $PRAD2_DATABASE_DIR/hycal_map.json (default ~/work/PRad/prad2evviewer/database)
Writes ../figs/exotic_events/
         exotic_types_examples.pdf    one example per type (1a, 2c, 3a, 4a of P1), 2 x 2
         exotic_module_maps.pdf       f1, f3 and saturated (f4 + f5) rate per module, 2.24 GeV H2 and Ta X17
         exotic_amp_vs_w50.pdf        fitted A_max vs w50 for the four run types, f2 rows marked
         exotic_saturated_seeds.pdf   seed energy fraction vs cluster size of the saturated elastic seeds
         exotic_energy_scale.pdf      energy-scale shift S per seed-radius ring, f1 and f3 corrections (P8)
         exotic_rates_vs_energy.pdf   per-channel category rate vs energy slice, cluster any-flag fraction dashed (P12)
         exotic_energy_removed.pdf    energy removed per calibration sample if every flagged channel is dropped (P13)
Needs pypdf, matplotlib, numpy.  Run from anywhere; paths are relative to this file.
"""
import csv
import io
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PatchCollection  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from pypdf import PageObject, PdfReader, PdfWriter, Transformation  # noqa: E402

HERE = Path(__file__).resolve().parent
SRC = HERE.parent.parent / "pileup_analysis" / "exotic_events"
OUT = HERE.parent / "figs" / "exotic_events"
DB = Path(os.environ.get("PRAD2_DATABASE_DIR", str(Path.home() / "work/PRad/prad2evviewer/database")))
LABEL = {"0p73GeV_H2": "0.73 GeV H$_2$", "2p24GeV_H2": "2.24 GeV H$_2$", "3p49GeV_H2": "3.49 GeV H$_2$",
         "2p24GeV_Ta_X17": "2.24 GeV Ta"}
plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "legend.fontsize": 8, "axes.linewidth": 0.8,
                     "axes.spines.top": False, "axes.spines.right": False})


def tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def crop(page, box):
    """Return a new page holding the (left, bottom, right, top) box of `page`, in points."""
    left, bottom, right, top = box
    new = PageObject.create_blank_page(width=right - left, height=top - bottom)
    new.merge_transformed_page(page, Transformation().translate(-left, -bottom), expand=False)
    return new


def white_mask(w, h, x0, y1):
    """A w x h pt page that is transparent except for a white box x >= x0, y <= y1."""
    fig = plt.figure(figsize=(w / 72, h / 72))
    fig.patch.set_alpha(0)
    fig.patches.append(Rectangle((x0 / w, 0), 1 - x0 / w, y1 / h, transform=fig.transFigure, facecolor="white", edgecolor="none"))
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", transparent=True)
    plt.close(fig)
    buf.seek(0)
    return PdfReader(buf).pages[0]


def types_examples():
    """2 x 2 crops of P1 page 1 (3168 x 1152 pt, gridspec 2 x 6): examples 1a, 2c, 3a, 4a."""
    page = PdfReader(SRC / "exotic_P1_quantities_20261117.pdf").pages[0]
    W, H = 3168.0, 1152.0
    pitch, width = 0.16376 * W, 0.14621 * W          # column pitch and axis width of the 2 x 6 grid
    cols = {j: (0.03 * W + j * pitch - 40, 0.03 * W + j * pitch + width + 50) for j in range(6)}   # +50: titles overflow the axis
    rows = {0: (0.595 * H, 0.905 * H), 1: (0.115 * H, 0.425 * H)}   # title + waveform + residual strip
    axtop = {0: 0.86 * H, 1: 0.38207 * H}                        # top of the waveform axis of each row
    cells = [(0, 0), (1, 0), (1, 1), (1, 4)]                      # 1a, 2c, 3a, 4a
    panels = []
    for i, j in cells:
        box = (cols[j][0], rows[i][0], cols[j][1], rows[i][1])
        panel = crop(page, box)
        # white mask over the tick labels of the neighbouring column (right strip, below the title)
        panel.merge_page(white_mask(box[2] - box[0], box[3] - box[1], x0=40 + width + 14, y1=axtop[i] - box[1] + 2))
        panels.append(panel)
    pw, ph = float(panels[0].mediabox.width), float(panels[0].mediabox.height)
    gap = 12.0
    sheet = PageObject.create_blank_page(width=2 * pw + gap, height=2 * ph + gap)
    for n, p in enumerate(panels):
        dx, dy = (n % 2) * (pw + gap), (1 - n // 2) * (ph + gap)
        sheet.merge_transformed_page(p, Transformation().translate(dx, dy), expand=False)
    w = PdfWriter()
    w.add_page(sheet)
    w.write(OUT / "exotic_types_examples.pdf")


def module_maps():
    geo = {m["n"]: m["geo"] for m in json.load(open(DB / "hycal_map.json"))}
    rows = [r for r in tsv(SRC / "exotic_P4_module_rates_20261117.tsv")
            if r["threshold_set"] == "uniform" and r["basis"] == "excl"]
    cats = [("f1", "f1 (tail)"), ("f3", "f3 (smooth-wide)"), ("sat", "f4 + f5 (saturated)")]
    rts = ["2p24GeV_H2", "2p24GeV_Ta_X17"]
    vmin, vmax = 1e-3, 3.0                                     # % ; log scale shared by all panels
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 5.0), layout="constrained")
    rect = lambda m: Rectangle((geo[m]["x"] - geo[m]["sx"] / 2, geo[m]["y"] - geo[m]["sy"] / 2),  # noqa: E731
                               geo[m]["sx"], geo[m]["sy"])
    for i, rt in enumerate(rts):
        for j, (flag, title) in enumerate(cats):
            ax = axs[i, j]
            sel = [r for r in rows if r["run_type"] == rt and r["flag"] == flag]
            few = [r for r in sel if int(r["n_kept"]) < 30]
            has = [r for r in sel if int(r["n_kept"]) >= 30 and int(r["n_flag"]) > 0]
            zero = [r for r in sel if int(r["n_kept"]) >= 30 and int(r["n_flag"]) == 0]
            ax.add_collection(PatchCollection([rect(r["module"]) for r in few], facecolor="0.85", edgecolor="0.7", lw=0.15))
            ax.add_collection(PatchCollection([rect(r["module"]) for r in zero], facecolor="white", edgecolor="0.7", lw=0.15))
            pc = PatchCollection([rect(r["module"]) for r in has], cmap="viridis", edgecolor="0.35", lw=0.1, norm=LogNorm(vmin, vmax))
            pc.set_array(np.clip([100 * float(r["rate"]) for r in has], vmin, vmax))
            ax.add_collection(pc)
            nf, nk = sum(int(r["n_flag"]) for r in sel), sum(int(r["n_kept"]) for r in sel)
            ax.set_xlim(-370, 370)
            ax.set_ylim(-370, 370)
            ax.set_aspect("equal")
            pw = (nf + 0.5) / (nk + 1)                         # Wilson 68 % (z = 1) interval of the total fraction
            hw = np.sqrt(nf / nk * (1 - nf / nk) / nk + 0.25 / nk ** 2) * nk / (nk + 1)
            ax.set_title(f"{LABEL[rt]}, {title}\n{nf} / {nk} = {100 * nf / nk:.3g} %\n[{100 * max(pw - hw, 0):.3g}, {100 * (pw + hw):.3g}] (68 % Wilson)", fontsize=7)
            if i == 1:
                ax.set_xlabel("x (mm)")
            if j == 0:
                ax.set_ylabel("y (mm)")
    cb = fig.colorbar(pc, ax=axs, shrink=0.8, extend="both")
    cb.set_label("rate per module (% of kept channels)")
    fig.savefig(OUT / "exotic_module_maps.pdf")
    plt.close(fig)


def amp_vs_w50():
    """P6 page 1 (1152 x 828 pt) without its page title line."""
    page = PdfReader(SRC / "exotic_P6_amp_vs_w50_20261118.pdf").pages[0]
    w = PdfWriter()
    w.add_page(crop(page, (0, 0, 1152, 790)))
    w.write(OUT / "exotic_amp_vs_w50.pdf")


def saturated_seeds():
    rows = tsv(SRC / "exotic_P7_saturated_seeds_20261118.tsv")
    fig, ax = plt.subplots(figsize=(5.2, 3.4), layout="constrained")
    mark = {"2p24GeV_H2": ("s", "2.24 GeV H$_2$"), "3p49GeV_H2": ("^", "3.49 GeV H$_2$")}
    rng = np.random.default_rng(1)
    for rt, (m, lab) in mark.items():
        sel = [r for r in rows if r["dataset"] == rt]
        x = np.array([int(r["n_members"]) for r in sel]) + rng.uniform(-0.18, 0.18, len(sel))
        y = [float(r["seed_frac"]) for r in sel]
        c = [float(r["chi2_dof"]) for r in sel]
        sc = ax.scatter(x, y, c=c, marker=m, s=34, cmap="plasma", norm=LogNorm(1, 1500), edgecolor="0.3", lw=0.4, label=lab)
    ax.axhline(0.95, color="0.4", ls="--", lw=0.8)
    ax.text(10.5, 0.955, "V3: seed fraction > 0.95", fontsize=8, ha="right", va="bottom", color="0.3")
    ax.annotate("W398 (4)", xy=(16, 0.79), xytext=(12.5, 0.86), fontsize=8, arrowprops=dict(arrowstyle="-", lw=0.6, color="0.3"))
    ax.annotate("W457 (17)", xy=(1.5, 0.997), xytext=(3.2, 0.905), fontsize=8, arrowprops=dict(arrowstyle="-", lw=0.6, color="0.3"))
    ax.set_xlabel("cluster members")
    ax.set_ylabel("seed energy fraction")
    ax.set_xlim(0.3, 20.5)
    ax.set_ylim(0.75, 1.01)
    cb = fig.colorbar(sc, ax=ax, shrink=0.9)
    cb.set_label(r"$\chi^2/\mathrm{dof}$ of the seed fit")
    ax.legend(loc="lower left", frameon=False)
    fig.savefig(OUT / "exotic_saturated_seeds.pdf")
    plt.close(fig)


def energy_scale():
    """S = sum dE / sum E (%) per seed-radius ring from the P8 scale table: Ta f1 (bs) and f3 (st), elastic and Moller; 2.24 and 3.49 GeV H2 f1 (bs), elastic (0.73 GeV has no f1 member)."""
    rows = [r for r in tsv(SRC / "exotic_P8_energy_error_20261118.tsv") if r["table"] == "scale"]
    rings = ["0-120", "120-170", "170-230", "230-320", "320-450", "450-700"]
    series = [("2p24GeV_Ta_X17", "elastic", "f1 tail-baseline", "bs", "Ta elastic, f1 sloped baseline", "#D55E00", "o", "full"),
              ("2p24GeV_Ta_X17", "elastic", "f3 smooth-wide", "st", "Ta elastic, f3 stretched single", "#0072B2", "o", "full"),
              ("2p24GeV_Ta_X17", "moller", "f1 tail-baseline", "bs", "Ta M\u00f8ller, f1", "#D55E00", "s", "none"),
              ("2p24GeV_Ta_X17", "moller", "f3 smooth-wide", "st", "Ta M\u00f8ller, f3", "#0072B2", "s", "none"),
              ("2p24GeV_H2", "elastic", "f1 tail-baseline", "bs", "2.24 GeV H$_2$ elastic, f1", "#56B4E9", "^", "full"),
              ("3p49GeV_H2", "elastic", "f1 tail-baseline", "bs", "3.49 GeV H$_2$ elastic, f1", "#009E73", "^", "full")]
    fig, ax = plt.subplots(figsize=(5.6, 3.3), layout="constrained")
    x0 = np.arange(len(rings))
    for n, (rt, pop, cat, var, lab, col, m, fill) in enumerate(series):
        sel = {r["ring"]: r for r in rows if r["dataset"] == rt and r["population"] == pop and r["category"] == cat and r["variant"] == var}
        x = [x0[i] + (n - 2.5) * 0.09 for i, rg in enumerate(rings) if rg in sel and sel[rg]["S"] not in ("", "nan") and int(sel[rg]["n_flag"]) > 0]
        y = [100 * float(sel[rg]["S"]) for rg in rings if rg in sel and sel[rg]["S"] not in ("", "nan") and int(sel[rg]["n_flag"]) > 0]
        e = [100 * float(sel[rg]["S_err"]) for rg in rings if rg in sel and sel[rg]["S"] not in ("", "nan") and int(sel[rg]["n_flag"]) > 0]
        ax.errorbar(x, y, yerr=e, fmt=m, ms=5, color=col, mfc=col if fill == "full" else "white", lw=0.8, capsize=2, label=lab)
    ax.axhline(0, color="0.5", lw=0.6)
    ax.set_xticks(x0)
    ax.set_xticklabels(rings)
    ax.set_xlabel("seed module radius (mm)")
    ax.set_ylabel("energy-scale shift $S$ (%)")
    ax.set_ylim(-0.022, 0.006)
    ax.legend(loc="lower right", ncol=2, frameon=False, fontsize=7)
    fig.savefig(OUT / "exotic_energy_scale.pdf")
    plt.close(fig)


def energy_removed():
    """-S1 = sum E_flagged / sum E_cluster (%) per run type and calibration sample from the P13 p8_scale table (ratio-estimator error)."""
    rows = [r for r in tsv(SRC / "exotic_P13_drop_flagged_20261118.tsv") if r["table"] == "p8_scale"]
    samples = [("all", "all clusters", "#0072B2", "o"), ("hiE", "$E > 0.5\\,E_b$", "#56B4E9", "^"),
               ("elastic", "elastic window", "#D55E00", "s"), ("moller", "M\u00f8ller pairs", "#009E73", "D")]
    fig, ax = plt.subplots(figsize=(5.0, 3.0), layout="constrained")
    x0 = np.arange(len(LABEL))
    for n, (smp, lab, col, m) in enumerate(samples):
        sel = {r["run_type"]: r for r in rows if r["sample"] == smp}
        x = [x0[i] + (n - 1.5) * 0.15 for i, rt in enumerate(LABEL) if rt in sel]
        y = [-100 * float(sel[rt]["S1"]) for rt in LABEL if rt in sel]
        e = [100 * float(sel[rt]["S1_err"]) for rt in LABEL if rt in sel]
        ax.errorbar(x, y, yerr=e, fmt=m, ms=5, color=col, lw=0.8, capsize=2, label=lab)
    ax.set_yscale("log")
    ax.set_xticks(x0)
    ax.set_xticklabels([LABEL[rt] for rt in LABEL])
    ax.set_ylabel("energy removed $-S$ (%)")
    ax.grid(axis="y", alpha=0.25, lw=0.4)
    ax.legend(loc="upper left", frameon=False, fontsize=7, ncol=2)
    fig.savefig(OUT / "exotic_energy_removed.pdf")
    plt.close(fig)


def rates_vs_energy():
    rows = tsv(SRC / "exotic_P12_rates_vs_energy_20261118.tsv")
    cats = [("f4 clipped-multi", "f4", "#D55E00", "o"), ("f5 saturated single", "f5", "#E69F00", "v"), ("f1 tail-baseline", "f1", "#0072B2", "s"),
            ("f2 wide", "f2", "#009E73", "^"), ("f3 smooth-wide", "f3", "#CC79A7", "D")]
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 4.6), sharey=True, layout="constrained")
    for ax, rt in zip(axs.ravel(), LABEL):
        sl = []
        for r in rows:
            if r["table"] == "chan" and r["run_type"] == rt and r["slice"] not in sl:
                sl.append(r["slice"])
        x = np.arange(len(sl))
        pick = lambda tab, cat, role=None: [next(r for r in rows if r["table"] == tab and r["run_type"] == rt and r["slice"] == s and r["category"] == cat and (role is None or r["role"] == role)) for s in sl]  # noqa: E731
        for cat, lab, col, mk in cats:
            rr = pick("chan", cat)
            ok = np.array([int(r["k"]) > 0 for r in rr])
            y = np.array([100 * float(r["rate"]) for r in rr])
            lo = np.array([100 * (float(r["rate"]) - float(r["rate_lo"])) for r in rr])
            hi = np.array([100 * (float(r["rate_hi"]) - float(r["rate"])) for r in rr])
            ax.errorbar(x[ok], y[ok], yerr=[lo[ok], hi[ok]], fmt="-" + mk, ms=3.5, lw=0.8, capsize=1.2, color=col, label=lab)
        cl = pick("cluster", "any", "any")
        okc = np.array([int(r["k"]) > 0 for r in cl])
        yc = np.array([100 * float(r["rate"]) for r in cl])
        lc = np.array([100 * (float(r["rate"]) - float(r["rate_lo"])) for r in cl])
        hc = np.array([100 * (float(r["rate_hi"]) - float(r["rate"])) for r in cl])
        ax.errorbar(x[okc], yc[okc], yerr=[lc[okc], hc[okc]], fmt="k--", lw=1.1, capsize=1.2, label="clusters, any flag (E_cluster)")
        ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels(sl, rotation=45, ha="right", fontsize=6.5)
        ax.set_title(LABEL[rt], fontsize=9)
        ax.grid(alpha=0.25, lw=0.4)
    for ax in axs[1]:
        ax.set_xlabel("energy slice (MeV): channel E (solid), E_cluster (dashed)", fontsize=8)
    for ax in axs[:, 0]:
        ax.set_ylabel("flagged fraction (%, Wilson 68 %)")
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=7, loc="outside upper center", ncol=3, columnspacing=1.0, handletextpad=0.4)
    fig.savefig(OUT / "exotic_rates_vs_energy.pdf")
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    types_examples()
    module_maps()
    amp_vs_w50()
    saturated_seeds()
    energy_scale()
    rates_vs_energy()
    energy_removed()
    print("wrote", *sorted(p.name for p in OUT.glob("exotic_*.pdf")))
