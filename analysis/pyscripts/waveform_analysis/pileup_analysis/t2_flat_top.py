#!/usr/bin/env python3
"""t2_flat_top.py -- T2: flat-top width of single-pulse waveforms (S1) vs amplitude, peak/integral ratio vs integral, saturation amplitude A_sat per module.

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/t2_flat_top.py [--n-events N] [--n-events-pwo N] [--cache-dir <dir>] [--out <dir>] [--rebuild]

S1, the shared pass and the cache: resid_common.py (default fit; k == 1, no overflow, A 50-3000 ADC; PbWO4 W* vs lead glass G*).
y = Y - fitted pedestal; template prediction = the fitted pulse (A x shape of the module's bank template at the fitted onset, i.e. same amplitude and phase), sampled like the data.
Width def 1 = number of consecutive samples (the run containing the maximum) with y >= 0.98 max;  def 2 = consecutive samples with y >= max - sigma_ref (module noise reference).
The same definitions on the template prediction give the expected width.  Ratio = max(y) / integral of y over the fit window (32 samples from 3 before the onset).
A_sat (per definition) = lower edge of the lowest amplitude bin from which the median of (data width - template width) is >= 1 sample in that and every higher bin
with >= MIN_BIN waveforms (amplitude bins SAT_EDGES); NA if none.
-> T2_flat_top_20261116.pdf / .tsv (pages: summary, width vs amplitude per class, ratio vs integral, A_sat table, 20 widest flat tops at A > 800 and 10 at A < 200).
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

SAT_EDGES = np.array([50, 100, 200, 400, 800, 1200, 1800, 3000.0])
MIN_BIN = 20
MIN_MOD = 200
N_WIDE_HI, N_WIDE_LO = 20, 10
FRAC = 0.02


def run_width(y, thr):
    """Samples of the run of y >= thr that contains the maximum, per row: (width, first, last) [n]."""
    ok = y >= thr[:, None]
    p = np.argmax(y, axis=1)
    ar = np.arange(y.shape[0])
    lab = np.cumsum(~ok, axis=1)
    same = (lab == lab[ar, p][:, None]) & ok
    idx = np.arange(y.shape[1])[None, :]
    first = np.where(same, idx, 10 ** 6).min(axis=1)
    last = np.where(same, idx, -1).max(axis=1)
    return same.sum(axis=1), first, last


def derive(d):
    y = d["Y"].astype(np.float64) - d["ped"][:, None]
    tm = y - d["resid"]                                   # template prediction (fitted pulse)
    mx, tmx = y.max(axis=1), tm.max(axis=1)
    sig = d["sig"].astype(np.float64)
    r = {}
    r["w1"], r["f1"], r["l1"] = run_width(y, (1 - FRAC) * mx)
    r["w2"], r["f2"], r["l2"] = run_width(y, mx - sig)
    r["t1"], _, _ = run_width(tm, (1 - FRAC) * tmx)
    r["t2"], _, _ = run_width(tm, tmx - sig)
    win = rc.window_mask(d["t0"])
    r["integ"] = (y * win).sum(axis=1)
    r["ratio"] = mx / r["integ"]
    r["tratio"] = tmx / (tm * win).sum(axis=1)
    r["mx"] = mx
    return r


def a_sat(a, diff):
    """Lower edge of the lowest SAT_EDGES bin from which median(diff) >= 1 in every populated (>= MIN_BIN) bin; None if none."""
    b = np.digitize(a, SAT_EDGES[1:-1])
    med, ok = np.full(len(SAT_EDGES) - 1, np.nan), np.zeros(len(SAT_EDGES) - 1, bool)
    for k in range(len(med)):
        s = b == k
        if s.sum() >= MIN_BIN:
            med[k], ok[k] = np.median(diff[s]), True
    good = ok & (med >= 1)
    for k in range(len(med)):
        if good[k] and np.all(good[k:][ok[k:]]):
            return float(SAT_EDGES[k])
    return None


def module_table(d, r):
    rows = []
    for m in np.unique(d["mi"]):
        s = d["mi"] == m
        a = d["a"][s]
        row = dict(module=str(d["names"][m]), cls=rc.CLASSES[0 if d["pwo4"][m] else 1], n=int(s.sum()))
        ab = d["ab"][s]
        for k in range(3):
            t = ab == k
            row[f"n{k}"] = int(t.sum())
            for key in ("w1", "w2", "t1", "t2", "ratio"):
                row[f"{key}_{k}"] = float(np.median(r[key][s][t])) if t.sum() else np.nan
        row["sat1"] = a_sat(a, (r["w1"] - r["t1"])[s]) if s.sum() >= MIN_BIN else None
        row["sat2"] = a_sat(a, (r["w2"] - r["t2"])[s]) if s.sum() >= MIN_BIN else None
        rows.append(row)
    return rows


def write_tsv(path, rows):
    f = lambda x: "NA" if (x is None or (isinstance(x, float) and np.isnan(x))) else (f"{x:.4g}" if isinstance(x, float) else str(x))    # noqa: E731
    hdr = ["module", "class", "n"] + [f"n_A{al}" for al in rc.AMP_LABELS]
    for nm in ("median_width_def1", "median_width_def2", "template_width_def1", "template_width_def2", "median_peak_over_integral"):
        hdr += [f"{nm}_A{al}" for al in rc.AMP_LABELS]
    hdr += ["A_sat_def1", "A_sat_def2"]
    with open(path, "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            vals = [r["module"], r["cls"], r["n"]] + [r[f"n{k}"] for k in range(3)]
            for key in ("w1", "w2", "t1", "t2", "ratio"):
                vals += [r[f"{key}_{k}"] for k in range(3)]
            vals += [r["sat1"], r["sat2"]]
            fh.write("\t".join(f(v) for v in vals) + "\n")


def prof(x, y, edges, q=(50, 16, 84)):
    """Per-bin percentiles of y in bins of x: [len(q), nbins] (nan where < 5)."""
    out = np.full((len(q), len(edges) - 1), np.nan)
    b = np.digitize(x, edges[1:-1])
    for k in range(len(edges) - 1):
        s = b == k
        if s.sum() >= 5:
            out[:, k] = np.percentile(y[s], q)
    return out


def page_summary(pdf, d, r, rows):
    fig = plt.figure(figsize=(16, 9))
    fig.text(0.03, 0.96, "T2: flat-top width of single-pulse waveforms", fontsize=18, fontweight="bold")
    L = [rc.run_runs_note(d)] + rc.counts_text(d) + ["Command: " + " ".join(sys.argv), 
        "Width def 1: consecutive samples (run containing the maximum) within 2 % of the maximum;  def 2: within sigma_ref (the module's ADC noise) of the maximum.  y = data - fitted pedestal.",
        "Template prediction: the fitted pulse (same amplitude and phase) with the same definitions.  Ratio = peak / integral of y over the fit window.",
        f"A_sat: lower edge of the lowest amplitude bin (edges {SAT_EDGES.astype(int).tolist()}) from which the median(data width - template width) >= 1 sample in that and all higher bins with >= {MIN_BIN} waveforms.", ""]
    for c, cn in enumerate(rc.CLASSES):
        for k, al in enumerate(rc.AMP_LABELS):
            s = (d["cls"] == c) & (d["ab"] == k)
            if s.sum() == 0:
                continue
            L.append(f"{cn} A {al} (n {int(s.sum())}): median width data/template def 1 {np.median(r['w1'][s]):.0f}/{np.median(r['t1'][s]):.0f}, def 2 {np.median(r['w2'][s]):.0f}/{np.median(r['t2'][s]):.0f};  "
                     f"mean width def 1 {r['w1'][s].mean():.2f}/{r['t1'][s].mean():.2f};  median peak/integral data {np.median(r['ratio'][s]):.4f}, template {np.median(r['tratio'][s]):.4f}")
    for c, cn in enumerate(rc.CLASSES):
        ms = [x for x in rows if x["cls"] == cn]
        L.append(f"{cn}: {len(ms)} modules, A_sat found (def 1) in {sum(x['sat1'] is not None for x in ms)}, (def 2) in {sum(x['sat2'] is not None for x in ms)}"
                 + (f"; def 1 values {sorted(x['sat1'] for x in ms if x['sat1'] is not None)}" if 0 < sum(x['sat1'] is not None for x in ms) <= 12 else ""))
    fig.text(0.02, 0.92, "\n".join(L), fontsize=8.5, va="top", family="monospace", linespacing=1.35)
    pdf.savefig(fig)
    plt.close(fig)


def page_width(pdf, d, r, c):
    rng = np.random.default_rng(1)
    s = np.nonzero(d["cls"] == c)[0]
    sub = rng.choice(s, size=min(30000, s.size), replace=False)
    edges = np.geomspace(rc.A_LO, rc.A_HI, 16)
    ctr = np.sqrt(edges[1:] * edges[:-1])
    fig, axs = plt.subplots(1, 2, figsize=(15, 6.5))
    fig.suptitle(f"{rc.CLASSES[c]}: flat-top width vs fitted amplitude ({s.size} S1 waveforms; points: random {sub.size}, jittered)", fontsize=13, fontweight="bold")
    for ax, (k, lab) in zip(axs, (("1", "def 1: within 2 % of max"), ("2", "def 2: within sigma_ref of max"))):
        ax.plot(d["a"][sub], r["w" + k][sub] + rng.uniform(-0.3, 0.3, sub.size), ".", ms=1.5, alpha=0.15, color="0.3")
        pd, pt = prof(d["a"][s], r["w" + k][s], edges), prof(d["a"][s], r["t" + k][s], edges)
        ax.fill_between(ctr, pd[1], pd[2], color="tab:blue", alpha=0.25, label="data 16-84 %")
        ax.plot(ctr, pd[0], "b.-", lw=1.5, label="data median")
        ax.plot(ctr, pt[0], "r.-", lw=1.5, label="template median (same A, phase)")
        ax.plot(ctr, pt[1], "r:", lw=0.8)
        ax.plot(ctr, pt[2], "r:", lw=0.8, label="template 16-84 %")
        for m in [m for m in np.unique(d["mi"][s]) if (d["mi"][s] == m).sum() >= MIN_MOD][:80]:
            t = s[d["mi"][s] == m]
            pm = prof(d["a"][t], r["w" + k][t], edges, q=(50,))[0]
            ax.plot(ctr, pm, color="tab:green", lw=0.3, alpha=0.4)
        ax.plot([], [], color="tab:green", lw=0.8, label="per-module medians (<= 80 modules)")
        ax.set_xscale("log")
        ax.set_ylim(0, min(25, np.nanmax(pd[2]) + 2))
        ax.set_xlabel("fitted amplitude A (ADC)")
        ax.set_ylabel("flat-top width (samples)")
        ax.set_title(lab, fontsize=10, loc="left")
        ax.legend(fontsize=7)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    pdf.savefig(fig)
    plt.close(fig)


def page_ratio(pdf, d, r):
    rng = np.random.default_rng(2)
    fig, axs = plt.subplots(1, 2, figsize=(15, 6.5))
    fig.suptitle("Peak / integral (fit window) vs integral: data (black/blue), template (red); thin green = per-module medians", fontsize=13, fontweight="bold")
    for c, ax in enumerate(axs):
        s = np.nonzero(d["cls"] == c)[0]
        sub = rng.choice(s, size=min(30000, s.size), replace=False)
        edges = np.geomspace(max(np.percentile(r["integ"][s], 0.5), 50), r["integ"][s].max(), 16)
        ctr = np.sqrt(edges[1:] * edges[:-1])
        ax.plot(r["integ"][sub], r["ratio"][sub], ".", ms=1.5, alpha=0.15, color="0.3")
        pd, pt = prof(r["integ"][s], r["ratio"][s], edges), prof(r["integ"][s], r["tratio"][s], edges)
        ax.fill_between(ctr, pd[1], pd[2], color="tab:blue", alpha=0.25)
        ax.plot(ctr, pd[0], "b.-", lw=1.5, label="data median, 16-84 %")
        ax.plot(ctr, pt[0], "r.-", lw=1.5, label="template median")
        for m in [m for m in np.unique(d["mi"][s]) if (d["mi"][s] == m).sum() >= MIN_MOD][:80]:
            t = s[d["mi"][s] == m]
            ax.plot(ctr, prof(r["integ"][t], r["ratio"][t], edges, q=(50,))[0], color="tab:green", lw=0.3, alpha=0.4)
        ax.set_xscale("log")
        lo, hi = np.nanmin(pd[1]), np.nanmax(pd[2])
        ax.set_ylim(lo - 0.3 * (hi - lo), hi + 0.3 * (hi - lo))
        ax.set_title(rc.CLASSES[c], fontsize=10, loc="left")
        ax.set_xlabel("integral of y over the fit window (ADC x samples)")
        ax.set_ylabel("max(y) / integral")
        ax.legend(fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    pdf.savefig(fig)
    plt.close(fig)


def page_table(pdf, rows):
    sat = [x for x in rows if x["sat1"] is not None or x["sat2"] is not None]
    per = 45
    pages = max(1, -(-len(sat) // per))
    for p in range(pages):
        fig = plt.figure(figsize=(11, 8.5))
        fig.text(0.04, 0.97, textwrap.fill(f"A_sat per module (page {p + 1}/{pages}): modules with a flat-top width above the template prediction by >= 1 sample (median) above A_sat", 100), fontsize=11, fontweight="bold", va="top")
        if not sat:
            fig.text(0.04, 0.88, textwrap.fill("Table empty: no module has A_sat (the median flat-top width never exceeds the template prediction by >= 1 sample above any amplitude, in either definition).", 110), fontsize=10, va="top")
        else:
            L = ["module  class   n      A_sat def 1   A_sat def 2   median width A800-3000 (data/template def 1)"]
            for x in sat[p * per:(p + 1) * per]:
                L.append(f"{x['module']:7s} {x['cls']:6s} {x['n']:6d}  {str(x['sat1']):>10s}   {str(x['sat2']):>10s}      {x['w1_2']:.1f} / {x['t1_2']:.1f}")
            fig.text(0.04, 0.93, "\n".join(L), fontsize=8, va="top", family="monospace")
        pdf.savefig(fig)
        plt.close(fig)


def pages_examples(pdf, d, r):
    hi = np.nonzero(d["a"] > 800)[0]
    lo = np.nonzero(d["a"] < 200)[0]
    key = lambda i: np.lexsort((-d["a"][i], -r["w2"][i], -r["w1"][i]))    # noqa: E731
    picks = [("A > 800", int(i)) for i in hi[key(hi)][:N_WIDE_HI]] + [("A < 200", int(i)) for i in lo[key(lo)][:N_WIDE_LO]]
    PER = 6
    for p0 in range(0, len(picks), PER):
        fig = plt.figure(figsize=(16, 12.5))
        fig.text(0.03, 0.985, va="top", s=textwrap.fill(f"Widest flat tops (def 1 width, then def 2): {N_WIDE_HI} at A > 800 and {N_WIDE_LO} at A < 200 ADC; examples {p0 + 1}-{min(p0 + PER, len(picks))} of {len(picks)}.  "
                              "Black dots data, red single-pulse fit; shaded blue = def 1 flat stretch (2 % of max), orange hatched = def 2 (within sigma_ref).  Strip: data - fit.", 190), fontsize=9.5)
        outer = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.12, left=0.05, right=0.99, top=0.95, bottom=0.04)
        for q, (grp, i) in enumerate(picks[p0:p0 + PER]):
            gs = outer[q // 2, q % 2].subgridspec(2, 1, height_ratios=[2.2, 1.0], hspace=0.08)
            ax, ar = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
            ped = float(d["ped"][i])
            model = d["Y"][i] - ped - d["resid"][i]
            rc.panel(ax, d["Y"][i], ped, model, float(d["t0"][i]),
                     f"{grp}: {rc.CLASSES[int(d['cls'][i])]} {d['modname'][i]} run {d['run'][i]} ev {d['ev'][i]}  A {d['a'][i]:.0f}  width def1 {r['w1'][i]} (tmpl {r['t1'][i]}), def2 {r['w2'][i]} (tmpl {r['t2'][i]})")
            for f_, l_, kw in ((r["f1"][i], r["l1"][i], dict(color="tab:blue", alpha=0.2)), (r["f2"][i], r["l2"][i], dict(facecolor="none", edgecolor="tab:orange", hatch="//", lw=1))):
                ax.axvspan(rc.TT[f_] - 2, rc.TT[l_] + 2, **kw)
            ax.tick_params(labelbottom=False)
            ar.plot(rc.TT, d["resid"][i], "k.-", lw=0.6, ms=3)
            ar.axhline(0, color="r", lw=0.5)
            ar.axvline(float(d["t0"][i]), color="tab:blue", ls=":", lw=0.8)
            ar.set_xlim(0, rc.TT[-1])
            ar.tick_params(labelsize=7)
            ar.set_xlabel("t (ns)", fontsize=7)
        pdf.savefig(fig)
        plt.close(fig)
    return len(picks)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    rc.add_common_args(ap)
    a = ap.parse_args()
    d = rc.load_s1(a)
    r = derive(d)
    rows = module_table(d, r)
    print(rc.run_runs_note(d))
    print("\n".join(rc.counts_text(d)))
    stem = Path(a.out) / f"T2_flat_top_{rc.STAMP}"
    write_tsv(str(stem) + ".tsv", rows)
    with PdfPages(str(stem) + ".pdf") as pdf:
        page_summary(pdf, d, r, rows)
        for c in range(2):
            page_width(pdf, d, r, c)
        page_ratio(pdf, d, r)
        page_table(pdf, rows)
        ne = pages_examples(pdf, d, r)
    print(f"[done] {stem}.pdf / .tsv, {len(rows)} modules, {ne} examples")


if __name__ == "__main__":
    main()
