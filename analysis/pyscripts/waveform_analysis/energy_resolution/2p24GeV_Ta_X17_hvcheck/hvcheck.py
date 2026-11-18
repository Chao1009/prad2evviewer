#!/usr/bin/env python3
"""HV-window check of the X17 Ta runs (25800, 26138, 26230) from the stored v3 outputs.

    source analysis/pyscripts/setup_env.sh; python3 hvcheck.py        (run from anywhere; one process)

Inputs : output/energy_resolution/control_phase2_v3_20261114/v3_defaults/2p24GeV_Ta_X17/{splits,gains}.tsv,
         output/energy_resolution/scan/2p24GeV_Ta_X17/run_<run>.npz (hit times), database/calibration/x17_*.json,
         the first EVIO split of each run (example waveforms only).
Outputs: summary.tsv, report.pdf, waveform_examples.tsv in this directory.
"""
import csv
import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

HERE = Path(__file__).resolve().parent
PYDIR = HERE.parents[2]
sys.path.insert(0, str(HERE.parents[1]))
V3 = PYDIR / "output/energy_resolution/control_phase2_v3_20261114/v3_defaults/2p24GeV_Ta_X17"
SCAN = PYDIR / "output/energy_resolution/scan/2p24GeV_Ta_X17"
DB = Path(os.environ["PRAD2_DATABASE_DIR"])
EVIO = Path("/Users/L00338853/work/PRad/data/evio")
RUNS = (25800, 26138, 26230)
CALIB = {25800: "x17_25825", 26138: "x17_26094", 26230: "x17_26214"}
DATE = {25800: "~07-11 to 07-13 (no entry; between 25783 on 07-10 and 25833 on 07-14/15)", 26138: "08-06", 26230: "08-12"}
WINDOW = {25800: "before", 26138: "before", 26230: "inside (08-10 to 08-14)"}
PATH = "g4"                                    # the 1.987 % arm of the v3 defaults
T_WIN = (140.0, 156.0)
E_MIN = 200.0                                  # MeV, hits used for the time distribution
OUTER = 5                                      # reference: PbWO4 layers >= OUTER

geo = [e for e in json.load(open(DB / "hycal_map.json")) if e["t"] == "PbWO4"]
LAY = {e["n"]: max(round(abs(e["geo"]["x"]) / e["geo"]["sx"] + 0.5), round(abs(e["geo"]["y"]) / e["geo"]["sy"] + 0.5)) - 1 for e in geo}
XY = {e["n"]: (e["geo"]["x"], e["geo"]["y"]) for e in geo}

rows = list(csv.DictReader(open(V3 / "splits.tsv"), delimiter="\t"))
gains = list(csv.DictReader(open(V3 / "gains.tsv"), delimiter="\t"))
res = {}
for run in RUNS:
    for path in (PATH, "legacy"):
        s = [r for r in rows if int(r["run"]) == run and r["path"] == path and r["variant"] == "all"]
        v = np.array([float(r["res_cb_oos"]) for r in s])
        res[run, path] = dict(n_splits=len(s), res=v.mean(), se=v.std(ddof=1) / np.sqrt(len(v)),
                              n_ep=sum(int(r["n_used"]) for r in s), mu=np.mean([float(r["mu_ep_oos"]) for r in s]))

fac = {}
for run, c in CALIB.items():
    fac[run] = {e["name"]: e["factor"] for e in json.load(open(DB / "calibration" / f"{c}.json"))}


def ratio(d):
    inner = [v for m, v in d.items() if m in LAY and LAY[m] <= 2 and v > 0]
    outer = [v for m, v in d.items() if m in LAY and LAY[m] >= OUTER and v > 0]
    return np.median(inner) / np.median(outer), len(inner), len(outer)


out, tdist, gmap = [], {}, {}
for run in RUNS:
    z = np.load(SCAN / f"run_{run}.npz")
    names = z["mod_name"]
    lay = np.array([LAY.get(n, 0) for n in names])
    ok = (z["h_k"] >= 1) & ~z["h_extra"] & (z["h_e"] >= E_MIN) & (lay[z["h_mod"]] > 0) & (z["h_t"] > 100) & (z["h_t"] < 200)
    t, ly = z["h_t"][ok], lay[z["h_mod"][ok]]
    tdist[run] = (t, ly)
    g = {r["module"]: float(r["g"]) for r in gains if int(r["run"]) == run}
    gmap[run] = g
    gr, gi, go = ratio(g)
    fr, fi, fo = ratio({m: fac[run].get(m, 0.0) for m in LAY})
    eff = {m: g[m] * fac[run][m] for m in g if m in fac[run]}
    er, ei, eo = ratio(eff)
    a, b = res[run, PATH], res[run, "legacy"]
    out.append(dict(run=run, date=DATE[run], hv_window=WINDOW[run], calib=CALIB[run], n_seen=int(z["n_seen"]), n_physics=int(z["n_pass"]),
                    n_splits=a["n_splits"], n_ep_g4=a["n_ep"], oos_cb_g4=a["res"], oos_cb_g4_se=a["se"],
                    oos_cb_legacy=b["res"], oos_cb_legacy_se=b["se"], mu_ep_oos_g4=a["mu"],
                    gain_ratio_inner_outer=gr, n_inner_gain=gi, n_outer_gain=go,
                    calib_factor_ratio_inner_outer=fr, effective_ratio=er, n_inner_eff=ei, n_outer_eff=eo,
                    t_median=np.median(t), t_mean=t.mean(), t_std=t.std(), frac_in_window=np.mean((t > T_WIN[0]) & (t < T_WIN[1])),
                    t_median_inner=np.median(t[ly <= 2]), t_median_outer=np.median(t[ly >= OUTER]), n_t=len(t)))
    print(run, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out[-1].items()}, flush=True)
    z.close()

with open(HERE / "summary.tsv", "w") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(out[0].keys())
    for o in out:
        w.writerow([f"{v:.5g}" if isinstance(v, float) else v for v in o.values()])

# ---- example waveforms: innermost layer, in-time, first split of each run
import energy_resolution_benchmark as B
from prad2py import dec
EX = {}
for run in RUNS:
    z = np.load(SCAN / f"run_{run}.npz")
    names = z["mod_name"]
    lay = np.array([LAY.get(n, 0) for n in names])
    sp0 = int(z["ev_split"].min())
    evn = z["ev_number"]
    ok = (z["h_k"] >= 1) & ~z["h_extra"] & (lay[z["h_mod"]] == 1) & (z["h_t"] > T_WIN[0]) & (z["h_t"] < T_WIN[1]) \
        & (z["h_e"] > 150) & (z["h_e"] < 2000) & (z["ev_split"][z["h_ev"]] == sp0)
    idx = np.nonzero(ok)[0]
    idx = idx[np.argsort(z["h_e"][idx])]
    pick = idx[np.linspace(0, len(idx) - 1, 10).astype(int)]
    want = {(int(evn[z["h_ev"][i]]), int(z["h_mod"][i])): (float(z["h_t"][i]), float(z["h_e"][i])) for i in pick}
    cfg = B.run_config(run)
    core = B.build(run, str(B.DB_DIR / cfg["calib"]))
    hy = core.hycal
    crate_map = {int(r.tag): int(r.crate) for r in core.daq_cfg.roc_tags}
    ch = dec.EvChannel()
    ch.set_config(core.daq_cfg)
    ch.open_auto(str(EVIO / f"prad_0{run}.evio.{sp0:05d}"))
    got = {}
    want_ev, max_en = {k[0] for k in want}, max(k[0] for k in want)
    while ch.read() == dec.Status.success and len(got) < len(want) and "stop" not in got:
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            d = ch.decode_event(i, with_ssp=False)
            if not d["ok"]:
                continue
            fe = d["event"]
            en = int(fe.info.event_number)
            if en > max_en:
                got["stop"] = None
                break
            if en not in want_ev:
                continue
            for ri in range(fe.nrocs):
                roc = fe.roc(ri)
                if not roc.present or roc.tag not in crate_map:
                    continue
                for sl in roc.present_slots():
                    slot = roc.slot(sl)
                    for c in slot.present_channels():
                        m = hy.module_by_daq(crate_map[roc.tag], sl, c)
                        if m is not None and (en, m.index) in want and slot.channel(c).nsamples == B.NS_EXPECT:
                            got[en, m.index] = np.asarray(slot.channel(c).samples, dtype=float)
    EX[run] = [(names[mi], want[en, mi][0], want[en, mi][1], got[en, mi]) for (en, mi) in got if en != 'stop' and (en, mi) in want]
    print(run, "waveforms", len(EX[run]), flush=True)
    z.close()
with open(HERE / "waveform_examples.tsv", "w") as f:
    f.write("run\tmodule\thit_t_ns\thit_E_MeV\tsamples_4ns\n")
    for run, L in EX.items():
        for n, t, e, y in L:
            f.write(f"{run}\t{n}\t{t:.1f}\t{e:.0f}\t{' '.join(f'{v:.0f}' for v in y)}\n")

# ---- report
col = {25800: "#1b9e77", 26138: "#7570b3", 26230: "#d95f02"}
with PdfPages(HERE / "report.pdf") as pdf:
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
    a = ax[0, 0]
    for k, o in enumerate(out):
        a.errorbar(k, o["oos_cb_g4"], o["oos_cb_g4_se"], fmt="o", color=col[o["run"]], capsize=4)
        a.annotate(f'{o["oos_cb_g4"]:.2f}±{o["oos_cb_g4_se"]:.2f}', (k, o["oos_cb_g4"]), textcoords="offset points", xytext=(8, 4), fontsize=8)
    a.axvspan(1.5, 2.5, color="0.9"); a.set_xticks(range(3)); a.set_xticklabels([f'{o["run"]}\n{o["hv_window"].split(" ")[0]}' for o in out])
    a.set_ylabel("out-of-sample CB σ/E (%), gate 4 ns"); a.set_title("Resolution per run (grey: raised-HV window)")
    a = ax[0, 1]
    for k, o in enumerate(out):
        a.plot(k, o["gain_ratio_inner_outer"], "o", color=col[o["run"]], label=f'fitted g {o["run"]}')
        a.plot(k + 0.15, o["calib_factor_ratio_inner_outer"], "s", color=col[o["run"]], mfc="none", label=None)
        a.plot(k + 0.3, o["effective_ratio"], "^", color=col[o["run"]], label=None)
    a.plot([], [], "ko", label="fitted g (inner/outer)"); a.plot([], [], "ks", mfc="none", label="calibration file factor"); a.plot([], [], "k^", label="factor x g")
    a.set_yscale("log"); a.axvspan(1.5, 2.5, color="0.9"); a.set_xticks(range(3)); a.set_xticklabels([str(o["run"]) for o in out]); a.legend(fontsize=7)
    a.set_title(f"Median ratio, layers 1-2 over layers >= {OUTER}")
    a = ax[1, 0]
    bins = np.arange(100, 200.1, 1.0)
    for run in RUNS:
        t, ly = tdist[run]
        a.hist(t, bins=bins, histtype="step", density=True, color=col[run], label=f"{run} all (median {np.median(t):.1f} ns)")
        a.hist(t[ly <= 2], bins=bins, histtype="step", density=True, color=col[run], ls="--", label=f"{run} layers 1-2")
    a.axvspan(*T_WIN, color="y", alpha=0.25); a.set_xlabel("hit time t0 + t50 (ns), E >= 200 MeV"); a.legend(fontsize=6); a.set_title("Hit time, window 140-156 ns shaded")
    a = ax[1, 1]; a.axis("off")
    txt = "run   date        n_ep(g4)  σ/E   ratio g  factor  t_med  in win\n" + "\n".join(
        f'{o["run"]}  {o["date"][:12]:12s} {o["n_ep_g4"]:7d}  {o["oos_cb_g4"]:.3f}  {o["gain_ratio_inner_outer"]:.3f}  {o["calib_factor_ratio_inner_outer"]:.3f}  {o["t_median"]:.1f}  {o["frac_in_window"]:.2f}' for o in out)
    a.text(0, 0.9, txt, family="monospace", fontsize=7, va="top")
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
    fig, ax = plt.subplots(1, 3, figsize=(11, 4.2))
    for a, run in zip(ax, RUNS):
        g = gmap[run]
        x = [XY[m][0] for m in g if m in XY]; y = [XY[m][1] for m in g if m in XY]; c = [g[m] for m in g if m in XY]
        s = a.scatter(x, y, c=c, s=12, cmap="viridis", vmin=0.85, vmax=1.15); a.set_aspect("equal"); a.set_title(f"{run}: fitted gain g"); plt.colorbar(s, ax=a, shrink=0.7)
    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
    for run in RUNS:
        fig, ax = plt.subplots(2, 5, figsize=(11, 5), sharex=True)
        for a, (n, t, e, y) in zip(ax.ravel(), EX[run]):
            a.plot(4 * np.arange(len(y)), y, ".-", ms=3, lw=0.7, color=col[run]); a.axvspan(*T_WIN, color="y", alpha=0.25); a.axvline(t, color="k", lw=0.6)
            a.set_title(f"{n}  t={t:.0f} ns  E={e:.0f} MeV", fontsize=7); a.tick_params(labelsize=6)
        fig.suptitle(f"Run {run}: innermost layer, in-time hits (window shaded, line = hit time), raw ADC"); fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
print("done")
