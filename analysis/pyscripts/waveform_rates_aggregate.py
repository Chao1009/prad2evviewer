#!/usr/bin/env python3
"""
waveform_rates_aggregate.py - aggregate the per-file pickles of
`waveform_good_bad.py rates` (output/waveform_rates_20261029/pkl) into
tables/rates.md and per-module HyCal heat maps (figs/).

usage: python3 waveform_rates_aggregate.py [OUTDIR]
"""
from __future__ import annotations

import glob
import json
import pickle
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "output/waveform_rates_20261029"
CLS = ("good", "bad", "mid", "nofit")
CHI2_EDGES = [0, .5, 1, 1.5, 2, 2.5, 3, 4, 6, 10, "inf"]
RUN_TYPE = {24386: "3.49 GeV H2", 24503: "3.49 GeV H2", 24512: "3.49 GeV H2",
            24590: "0.73 GeV H2", 24603: "0.73 GeV H2",
            24578: "0.73 GeV H2 (24578, non-standard trigger)",
            25228: "2.24 GeV H2", 25257: "2.24 GeV H2", 25260: "2.24 GeV H2",
            25299: "2.24 GeV H2", 25320: "2.24 GeV H2", 25308: "2.24 GeV C",
            25800: "2.24 GeV Ta/X17", 26138: "2.24 GeV Ta/X17",
            26230: "2.24 GeV Ta/X17"}


def binom(k, n):
    p = k / n if n else float("nan")
    return p, (np.sqrt(p * (1 - p) / n) if n else float("nan"))


def fmt(k, n):
    p, e = binom(k, n)
    return f"{100*p:.2f} +- {100*e:.2f}"


def main():
    files = sorted(glob.glob(str(OUT / "pkl/*.pkl")))
    recs = []
    for f in files:
        with open(f, "rb") as fh:
            d = pickle.load(fh)
        t = d["tot"]
        ab = np.zeros((3, 4), np.int64)
        for e in d["mods"].values():
            ab += e["ab"]
        recs.append(dict(run=int(d["run"]), split=int(d["split"]), trig=d["n_trigger"],
                         cl=d["cnt"]["cluster"], c=np.array([t[k] for k in CLS]),
                         ab=ab, mods=d["mods"]))
    recs.sort(key=lambda r: (r["run"], r["split"]))
    for r in recs:
        r["n"] = int(r["c"].sum())
        r["type"] = RUN_TYPE.get(r["run"], "unknown")
    L = []
    w = L.append
    runs = sorted({r["run"] for r in recs})
    types = []
    for r in recs:
        if r["type"] not in types:
            types.append(r["type"])
    w("# Good/bad MF fits of sub-50-ADC pulses: rates\n")
    w("Population: RawSum trigger, 1 HyCal cluster, < 20 modules fired, one clean peak "
      "20 <= h < 50 ADC, h >= 10 x ped rms (waveform_good_bad.py). good chi2/dof < 2, "
      "bad > 3, between 2-3, nofit = MF found no pulse. Fractions are of qualifying "
      "pulses; uncertainty = binomial (stat only, pulses in one event/file are "
      "correlated, so see the file-to-file columns for the realistic spread).\n")
    w(f"Files: {len(recs)}; runs: {len(runs)}\n")

    def block(title, groups):
        w(f"\n## {title}\n")
        w("| group | files | triggered events | qualifying pulses | good % | bad % | between % | nofit % |")
        w("|---|---|---|---|---|---|---|---|")
        for name, rs in groups:
            c = sum(r["c"] for r in rs)
            n = int(c.sum())
            w(f"| {name} | {len(rs)} | {sum(r['trig'] for r in rs)} | {n} | "
              + " | ".join(fmt(int(c[i]), n) for i in range(4)) + " |")

    block("Overall and per run type",
          [("ALL (14 R5 runs)", [r for r in recs if r["run"] != 24578])]
          + [(t, [r for r in recs if r["type"] == t]) for t in types])
    block("Per run", [(f"{ru} ({RUN_TYPE.get(ru)})", [r for r in recs if r["run"] == ru])
                      for ru in runs])

    # detector
    w("\n## Detector (PbWO4 = W*, lead glass = G*), R5 runs\n")
    w("| type | detector | qualifying | good % | bad % | between % | nofit % |")
    w("|---|---|---|---|---|---|---|")
    for t in ["ALL"] + types:
        for det in ("W", "G"):
            c = np.zeros(4, np.int64)
            for r in recs:
                if (t == "ALL" and r["run"] == 24578) or (t not in ("ALL", r["type"])):
                    continue
                for nme, e in r["mods"].items():
                    if nme[0] == det:
                        c += [e[k] for k in CLS]
            n = int(c.sum())
            w(f"| {t} | {'PbWO4' if det == 'W' else 'PbGlass'} | {n} | "
              + " | ".join(fmt(int(c[i]), n) for i in range(4)) + " |")

    # amplitude
    w("\n## Amplitude dependence (fractions within each bin)\n")
    w("| group | amp bin (ADC) | qualifying | good % | bad % | between % | nofit % |")
    w("|---|---|---|---|---|---|---|")
    for t in ["ALL"] + types:
        rs = [r for r in recs if (t == "ALL" and r["run"] != 24578) or r["type"] == t]
        ab = sum(r["ab"] for r in rs)
        for j, lab in enumerate(("20-30", "30-40", "40-50")):
            n = int(ab[j].sum())
            w(f"| {t} | {lab} | {n} | "
              + " | ".join(fmt(int(ab[j][i]), n) for i in range(4)) + " |")

    # chi2 hist
    w("\n## chi2/dof histogram (fitted pulses, R5 runs)\n")
    h = np.zeros(10, np.int64)
    for r in recs:
        if r["run"] != 24578:
            for e in r["mods"].values():
                h += e["hist"]
    w("| bin | " + " | ".join(f"{CHI2_EDGES[i]}-{CHI2_EDGES[i+1]}" for i in range(10)) + " |")
    w("|---|" + "---|" * 10)
    w("| fraction % | " + " | ".join(f"{100*x/h.sum():.2f}" for x in h) + " |")

    # per-file stability
    w("\n## File-to-file stability (bad fraction)\n")
    w("RMS = std of per-file bad fraction in the run; trend = OLS slope of per-file "
      "bad fraction vs file number, in percentage points per 100 files "
      "(10 files = ~1/10 of a run), with |t| its significance.\n")
    w("| run | type | files | mean bad % | RMS % | RMS/mean | trend pp/100 files | t | "
      "good RMS % | qualifying/event |")
    w("|---|---|---|---|---|---|---|---|---|---|")
    tr_sig = 0
    for ru in runs:
        rs = [r for r in recs if r["run"] == ru]
        bad = np.array([r["c"][1] / max(r["n"], 1) for r in rs])
        good = np.array([r["c"][0] / max(r["n"], 1) for r in rs])
        x = np.array([r["split"] for r in rs], float)
        slope = tstat = float("nan")
        if len(rs) >= 4:
            A = np.vstack([x, np.ones_like(x)]).T
            coef, res, *_ = np.linalg.lstsq(A, bad, rcond=None)
            resid = bad - A @ coef
            s2 = (resid ** 2).sum() / (len(x) - 2)
            se = np.sqrt(s2 / ((x - x.mean()) ** 2).sum())
            slope, tstat = coef[0] * 100 * 100, coef[0] / se if se > 0 else float("nan")
            tr_sig += int(abs(tstat) > 3)
        c = sum(r["c"] for r in rs)
        w(f"| {ru} | {RUN_TYPE[ru]} | {len(rs)} | {100*bad.mean():.2f} | "
          f"{100*bad.std(ddof=1) if len(rs) > 1 else float('nan'):.2f} | "
          f"{bad.std(ddof=1)/bad.mean() if len(rs) > 1 else float('nan'):.3f} | "
          f"{slope:.2f} | {tstat:.1f} | "
          f"{100*good.std(ddof=1) if len(rs) > 1 else float('nan'):.2f} | "
          f"{c.sum()/sum(r['trig'] for r in rs):.3f} |")
    w(f"\nRuns with trend |t| > 3: {tr_sig} of {len(runs)}.\n")
    w("\nPer-type file-to-file RMS of bad fraction (R5 files) and standard error of the "
      "mean over files:\n")
    w("| type | files | mean bad % | RMS % | SEM % | binomial % |")
    w("|---|---|---|---|---|---|")
    for t in types:
        rs = [r for r in recs if r["type"] == t]
        bad = np.array([r["c"][1] / max(r["n"], 1) for r in rs])
        c = sum(r["c"] for r in rs)
        w(f"| {t} | {len(rs)} | {100*bad.mean():.2f} | {100*bad.std(ddof=1):.2f} | "
          f"{100*bad.std(ddof=1)/np.sqrt(len(rs)):.2f} | {100*binom(int(c[1]), int(c.sum()))[1]:.3f} |")

    # per event
    w("\n## Rates per triggered event\n")
    w("Triggered event = RawSum-bit physics event (before the < 20 fired and 1-cluster cuts).\n")
    w("| group | triggered events | 1-cluster events | qualifying/trig event | bad fits/trig event | "
      "qualifying/1-cluster event | bad/1-cluster event |")
    w("|---|---|---|---|---|---|---|")
    for name, rs in [("ALL R5", [r for r in recs if r["run"] != 24578])] + \
                    [(t, [r for r in recs if r["type"] == t]) for t in types]:
        tr = sum(r["trig"] for r in rs)
        cl = sum(r["cl"] for r in rs)
        c = sum(r["c"] for r in rs)
        nq, nb = int(c.sum()), int(c[1])
        w(f"| {name} | {tr} | {cl} | {nq/tr:.4f} +- {np.sqrt(nq)/tr:.4f} | "
          f"{nb/tr:.4f} +- {np.sqrt(nb)/tr:.4f} | {nq/max(cl,1):.3f} | {nb/max(cl,1):.3f} |")

    # modules
    mod = {}
    for r in recs:
        if r["run"] == 24578:
            continue
        for nme, e in r["mods"].items():
            m = mod.setdefault(nme, np.zeros(4, np.int64))
            m += [e[k] for k in CLS]
    rows = [(n_, c) for n_, c in mod.items() if c.sum() >= 500]
    rows.sort(key=lambda x: -x[1][1] / x[1].sum())
    w(f"\n## Notable modules (R5 runs pooled, >= 500 qualifying pulses; {len(rows)} of "
      f"{len(mod)} modules)\n")
    fr = np.array([c[1] / c.sum() for _, c in rows])
    w(f"Module bad fraction: median {100*np.median(fr):.1f} %, 5-95 % range "
      f"{100*np.percentile(fr,5):.1f}-{100*np.percentile(fr,95):.1f} %.\n")
    w("| module | qualifying | bad % | good % | nofit % |")
    w("|---|---|---|---|---|")
    for lab, sel in (("highest bad", rows[:15]), ("lowest bad", rows[-8:])):
        w(f"| *{lab}* | | | | |")
        for n_, c in sel:
            w(f"| {n_} | {c.sum()} | {fmt(int(c[1]), int(c.sum()))} | "
              f"{100*c[0]/c.sum():.1f} | {100*c[3]/c.sum():.2f} |")

    (OUT / "tables").mkdir(exist_ok=True)
    (OUT / "tables/rates.md").write_text("\n".join(L) + "\n")
    print("wrote", OUT / "tables/rates.md")

    # heat maps
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle
    geo = {m["n"]: m for m in json.load(open(HERE.parent.parent / "database/hycal_map.json"))
           if "geo" in m}
    (OUT / "figs").mkdir(exist_ok=True)
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(OUT / "figs/bad_fraction_maps.pdf") as pdf:
        for t in types:
            rs = [r for r in recs if r["type"] == t]
            agg = {}
            for r in rs:
                for nme, e in r["mods"].items():
                    a = agg.setdefault(nme, np.zeros(4, np.int64))
                    a += [e[k] for k in CLS]
            fig, axes = plt.subplots(1, 2, figsize=(15, 7.5))
            for ax, det, lab in zip(axes, ("W", "G"), ("PbWO4", "PbGlass")):
                patches, vals = [], []
                for nme, a in agg.items():
                    if nme[0] != det or a.sum() < 30 or nme not in geo:
                        continue
                    g = geo[nme]["geo"]
                    patches.append(Rectangle((g["x"] - g["sx"] / 2, g["y"] - g["sy"] / 2),
                                             g["sx"], g["sy"]))
                    vals.append(a[1] / a.sum())
                allp = [geo[n_]["geo"] for n_ in geo if n_[0] == det]
                pc = PatchCollection(patches, cmap="viridis", ec="0.5", lw=0.2)
                pc.set_array(np.array(vals))
                ax.add_collection(pc)
                ax.set_xlim(min(g["x"] - g["sx"] for g in allp), max(g["x"] + g["sx"] for g in allp))
                ax.set_ylim(min(g["y"] - g["sy"] for g in allp), max(g["y"] + g["sy"] for g in allp))
                ax.set_aspect("equal")
                ax.set_title(f"{lab}: bad fraction ({len(vals)} modules, >= 30 pulses)")
                ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
                fig.colorbar(pc, ax=ax, shrink=0.8, label="bad (chi2/dof > 3) fraction")
            fig.suptitle(f"{t}: {len(rs)} files, sub-50-ADC pulses")
            fig.tight_layout()
            tag = t.split(" (")[0].replace(" ", "_").replace("/", "_")
            fig.savefig(OUT / f"figs/bad_fraction_map_{tag}.png", dpi=130)
            pdf.savefig(fig)
            plt.close(fig)
    print("wrote figs")


if __name__ == "__main__":
    main()
