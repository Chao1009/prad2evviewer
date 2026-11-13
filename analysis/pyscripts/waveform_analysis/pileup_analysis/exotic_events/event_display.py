#!/usr/bin/env python3
"""event_display.py -- one HyCal event: energy map, waveforms with the pile-up fit, same-board coherence check.

    source setup_env.sh      # PRAD2_DATABASE_DIR, prad2py
    python3 waveform_analysis/pileup_analysis/exotic_events/event_display.py \
        --run 24590 --event 4084847 2069998 --module W121 W493 --out <dir> [--persistence W121]

One PDF per event, <out>/event_display_run<run>_ev<event>_<module>.pdf, plus <out>/summary.txt (all events of the call; with
--keep-png also one PNG per page).  Event = the EVIO event number (fe.info.event_number), the number pileup_rates.py stores
as `evnum` and prints in its example titles.

Page 1  HyCal map, colour = module energy E_main (MeV, logarithmic) + linear 9 x 9 zoom with the energy in every cell, reconstructed clusters.
Page 2  highlighted module and its 8 neighbours: pedestal-subtracted waveform with the pile-up fit (pileup_rates.draw_panel's overlay).
Page 3  all channels of the highlighted module's fADC board (same y range), event-wide mean waveform of the channels without a fitted pulse
        (common-mode estimate), correlations.
Page 4  (only with --persistence <module>) the module's baseline noise over a stride sample of the event's EVIO split.

Reused by import, not copied: pileup_rates.py (Bank, MF, fit, run_config, module_gain, noise_reference, pedestal, bad_waveform,
neighbours, constants) -- the fit of pileup_rates.scan is rerun on the event with the same per-file noise reference (first N_REF_EVENTS
decoded events of the event's EVIO file, as pileup_rates.read_reference).  Replicated here because the originals are closures inside
pileup_rates.scan / energy_resolution_benchmark.py (being edited): the event locating, the noise-reference reading loop and the HyCalCluster set-up
(energy_resolution_benchmark.clusterers: builder configuration, MIN_CENTER_E 30 MeV, seed_time_window -1 = legacy, no time gate).
Energy: E_main = A_main x gain, A_main the largest fitted pulse amplitude of the channel (ADC, pedestal-subtracted), gain = calibration factor
x I_int of the run type's template bank (pileup_rates.module_gain).  Channels removed by the chain (bad waveform, A_main < KEEP_ADC, truncated main
pulse) have no energy.  Clustering: per module the largest pulse with onset in the pre-window (PRE_LO, PRE_HI) = (130, 200) ns, energy A x gain.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import glob
import pickle
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True                         # no __pycache__ residue next to pileup_rates.py
sys.path.insert(0, str(HERE.parent))
if not os.environ.get("PRAD2_DATABASE_DIR"):
    sys.exit("[ERROR] PRAD2_DATABASE_DIR is not set: run `source setup_env.sh` (analysis/pyscripts) first")
import pileup_rates as pr                              # noqa: E402  (module import only sets rcParams; main is guarded)
from prad2py import dec, det                           # noqa: E402

matplotlib.rcdefaults()
plt.rcParams.update({"font.size": 11, "axes.titlesize": 11, "axes.labelsize": 12, "axes.titleweight": "bold",
                     "axes.spines.top": False, "axes.spines.right": False, "figure.titlesize": 15,
                     "figure.titleweight": "bold", "axes.edgecolor": "0.25", "text.color": "0.12",
                     "axes.labelcolor": "0.12", "xtick.color": "0.25", "ytick.color": "0.25"})

CMAP = plt.get_cmap("viridis")          # perceptually uniform, sequential, one scale for the magnitude
GREY, DEAD_GREY = "#dcdcdc", "#f2f2f2"
PULSE_COLORS = ("tab:blue", "tab:orange", "tab:green", "tab:purple", "tab:brown", "tab:pink")   # as pileup_rates.draw_panel
NS, DT = pr.NS_EXPECT, pr.DT
TT = DT * np.arange(NS)
MIN_CENTER_E = 30.0                     # MeV, seed threshold (energy_resolution_benchmark.MIN_CENTER_E)
CLUSTER_FIELDS = ("min_module_energy", "min_center_energy", "min_cluster_energy", "min_cluster_size",
                  "corner_conn", "split_iter", "least_split", "log_weight_thres")
N_REF = pr.N_REF_EVENTS


# ---------------------------------------------------------------------------------------------------------------- setup
class Ctx:
    """Run context: pipeline core, module tables, bank, gains, config."""

    def __init__(self, run):
        self.run = run
        self.cfg = pr.run_config(run)
        self.rt = pr.RUN_TYPE_OF.get(run)
        if self.rt is None:
            sys.exit(f"[ERROR] run {run} is not in pileup_rates.RUN_TYPE_OF")
        b = det.PipelineBuilder()
        b.set_hycal_calib(self.cfg["calib"])
        b.set_run_number(run)
        self.core = b.build()
        hy = self.hycal = self.core.hycal
        self.mods = [hy.module(i) for i in range(hy.module_count())]
        self.names = np.array([m.name for m in self.mods])
        self.pwo4 = np.array([bool(m.is_pwo4()) for m in self.mods])
        self.x = np.array([m.x for m in self.mods])
        self.y = np.array([m.y for m in self.mods])
        self.sx = np.array([m.size_x for m in self.mods])
        self.sy = np.array([m.size_y for m in self.mods])
        self.index = {n: i for i, n in enumerate(self.names)}
        self.dead = np.array([n in set(self.cfg["dead"]) for n in self.names])
        self.crate_of_tag = {int(r.tag): int(r.crate) for r in self.core.daq_cfg.roc_tags}
        self.daq = {i: (int(m.daq.crate), int(m.daq.slot), int(m.daq.channel)) for i, m in enumerate(self.mods) if m.is_hycal()}
        self.by_daq = {v: k for k, v in self.daq.items()}
        self.bank = pr.Bank(pr.BANK_ROOT / self.rt / f"templates_{self.rt}.npz", self.names, self.pwo4)
        self.mf = pr.MF(self.bank)
        self.gain = pr.module_gain(self.cfg, self.rt, self.names)
        self.beam = self.cfg["beam"]
        self.cache = {}

    def module_of(self, tag, crate, sl, c):
        """HyCal module index of a DAQ channel, -1 if none or dead (as pileup_rates.scan.module_of)."""
        k = (tag, sl, c)
        mi = self.cache.get(k)
        if mi is None:
            m = self.hycal.module_by_daq(crate, sl, c)
            mi = self.cache[k] = m.index if (m is not None and m.is_hycal() and m.name not in set(self.cfg["dead"])) else -1
        return mi

    def channels(self, fe):
        """[(module index, samples)] of the HyCal channels with NS samples of a decoded event (dead modules left out as in the scan)."""
        out = []
        for ri in range(fe.nrocs):
            roc = fe.roc(ri)
            crate = self.crate_of_tag.get(roc.tag) if roc.present else None
            if crate is None:
                continue
            for sl in roc.present_slots():
                slot = roc.slot(sl)
                for c in slot.present_channels():
                    mi = self.module_of(roc.tag, crate, sl, c)
                    cd = slot.channel(c)
                    if mi >= 0 and cd.nsamples == NS:
                        out.append((mi, np.asarray(cd.samples, dtype=np.uint16)))
        return out

    def open(self, f):
        ch = dec.EvChannel()
        ch.set_config(self.core.daq_cfg)
        if ch.open_auto(str(f)) != dec.Status.success:
            sys.exit(f"[ERROR] cannot open {f}")
        return ch


def physics_blocks(ch):
    """Yield the number of sub-events of each physics block read from an open channel."""
    while ch.read() == dec.Status.success:
        if ch.scan() and ch.get_event_type() == dec.EventType.Physics:
            yield ch.get_n_events()


# ------------------------------------------------------------------------------------------------------ event locating
def first_event_number(ctx, f):
    ch = ctx.open(f)
    for n in physics_blocks(ch):
        ch.select_event(0)
        en = int(ch.info().event_number)
        ch.close()
        return en
    ch.close()
    return None


def locate(ctx, evio_dir, evnum):
    """EVIO split holding event number evnum: the last split of the run whose first physics event number is <= evnum.
    Returns (path, split ordinal in the sorted file list, [paths])."""
    files = sorted(glob.glob(str(Path(evio_dir) / f"prad_{ctx.run:06d}.evio.*")))
    if not files:
        sys.exit(f"[ERROR] no prad_{ctx.run:06d}.evio.* in {evio_dir}")
    firsts = [first_event_number(ctx, f) for f in files]
    ok = [i for i, e in enumerate(firsts) if e is not None and e <= evnum]
    if not ok:
        sys.exit(f"[ERROR] event {evnum} is before the first event of run {ctx.run}")
    return files[ok[-1]], ok[-1], files


def find_event(ctx, f, evnum):
    """Decode event number evnum of EVIO file f: skip block by block, decode only the matching event.  Returns fe or None."""
    ch = ctx.open(f)
    for n in physics_blocks(ch):
        for i in range(n):
            ch.select_event(i)
            if int(ch.info().event_number) == evnum:
                d = ch.decode_event(i, with_ssp=False)
                ch.close()
                return d["event"] if d["ok"] else None
    ch.close()
    return None


def noise_reference(ctx, f):
    """Per-module (n_ref, b_ref, sigma_ref) of EVIO file f, as pileup_rates.scan.read_reference (first N_REF decoded events, all modules)."""
    ch = ctx.open(f)
    X, M, rows, mis = [], [], [], []

    def take():
        if rows:
            Y, mi = np.stack(rows), np.array(mis)
            rows.clear(), mis.clear()
            good = ~pr.bad_waveform(Y)
            X.append(Y[good, :pr.PED_NSAMPLES])
            M.append(mi[good])

    n = 0
    while n < N_REF and ch.read() == dec.Status.success:
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            d = ch.decode_event(i, with_ssp=False)
            if not d["ok"]:
                continue
            n += 1
            for mi, s in ctx.channels(d["event"]):
                rows.append(s)
                mis.append(mi)
            if n % 200 == 0:
                take()
            if n >= N_REF:
                break
    take()
    ch.close()
    return pr.noise_reference(np.concatenate(X), np.concatenate(M), len(ctx.mods))


# ------------------------------------------------------------------------------------------------------------ event fit
def fit_event(ctx, fe, b_ref, sig_ref):
    """Fit every HyCal channel of the event as pileup_rates.scan.flush does.  Returns dict of per-channel arrays (row order = ch['mi'])."""
    chs = ctx.channels(fe)
    mi = np.array([c[0] for c in chs])
    Y = np.stack([c[1] for c in chs])
    n = len(mi)
    bad = pr.bad_waveform(Y)
    ped0 = pr.pedestal(Y)
    x = np.median(Y[:, :pr.PED_NSAMPLES], axis=1)
    sel = np.nonzero((Y.max(axis=1) - x >= pr.PRECUT_ADC) & ~bad)[0]
    K = pr.KMAX
    r = dict(mi=mi, Y=Y, bad=bad, ped30=ped0.astype(np.float64), ped=ped0.astype(np.float64), k=np.zeros(n, int), t=np.full((n, K), np.nan), a=np.zeros((n, K)),
             comp=np.zeros((n, K, NS)), chi2=np.full(n, np.nan), trunc=np.zeros(n, bool), nclip=pr.clipped(Y), fitted=np.zeros(n, bool))
    if sel.size:
        f = pr.fit(ctx.mf, Y[sel], ctx.bank.row[mi[sel]], sig_ref[mi[sel]], b_ref[mi[sel]])
        for key, v in (("ped", f["ped"]), ("k", f["k"]), ("t", f["t"]), ("a", f["a"]), ("comp", f["comp"]), ("chi2", f["chi2"]), ("trunc", f["trunc_main"])):
            r[key][sel] = v
        r["fitted"][sel] = True
    r["row"] = {int(m): i for i, m in enumerate(mi)}
    amain = r["a"].max(axis=1)
    r["amain"] = amain
    keep = (amain >= pr.KEEP_ADC) & ~bad
    r["keep"] = keep
    r["emain"] = np.where(keep & ~r["trunc"], amain * ctx.gain[mi], np.nan)       # chain: truncated main pulses have no E_main
    return r


def pedestal_stats(y):
    """rms (std of the residuals) and slope (ADC per 100 ns, least squares) of the first PED_NSAMPLES samples."""
    s = y[:pr.PED_NSAMPLES].astype(float)
    t = DT * np.arange(s.size)
    slope = np.polyfit(t, s, 1)[0] * 100.0
    ped = pr.pedestal(y[None, :])[0]
    ok = np.abs(s - np.median(s)) <= max(2.0 * 1.4826 * np.median(np.abs(s - np.median(s))), 1.0)
    return (float(np.sqrt(np.mean((s[ok] - ped) ** 2))) if ok.sum() > 2 else float(np.std(s))), float(slope), float(np.std(s))


# ----------------------------------------------------------------------------------------------------------- clustering
def make_clusterer(ctx):
    """HyCalCluster as energy_resolution_benchmark.clusterers builds the legacy path."""
    base = ctx.core.hycal_cluster_cfg
    cfg = det.HyCalClusterConfig()
    for f in CLUSTER_FIELDS:
        setattr(cfg, f, getattr(base, f))
    cfg.min_center_energy = MIN_CENTER_E
    cfg.seed_time_window = -1.0
    cl = det.HyCalCluster(ctx.hycal)
    cl.set_config(base)
    cl.set_config(cfg)
    return cl


def cluster_event(ctx, r):
    """Clusters of the event (list of dicts) from the largest in-window pulse of each kept module."""
    cl = make_clusterer(ctx)
    cl.clear()
    hits = {}
    for i in np.nonzero(r["keep"] & ~r["trunc"])[0]:
        k = int(r["k"][i])
        t, a = r["t"][i][:k], r["a"][i][:k]
        win = (t > pr.PRE_LO) & (t < pr.PRE_HI)
        if not win.any() or not np.isfinite(ctx.gain[r["mi"][i]]):
            continue
        j = np.where(win, a, -1).argmax()
        e = a[j] * ctx.gain[r["mi"][i]]
        hits[int(r["mi"][i])] = (e, t[j])
        cl.add_hit(int(r["mi"][i]), float(e), float(t[j]))
    cl.form_clusters()
    out = []
    for c in cl.reconstruct_matched():
        mi = [h.index for h in c.cluster.hits]
        me = [h.energy for h in c.cluster.hits]
        out.append(dict(seed=int(c.cluster.center.index), seed_e=float(c.cluster.center.energy), e=float(c.hit.energy), esum=float(np.sum(me)),
                        x=float(c.hit.x), y=float(c.hit.y), mi=mi, me=me, n=len(mi)))
    out.sort(key=lambda c: -c["e"])
    return out, hits


def chain_cluster(ctx, r, nbr):
    """pileup_rates' own leading cluster: seed = kept channel with the largest in-window pulse, members within CL_RADIUS x mean size with
    an in-window pulse.  Returns (seed index, [member indices]) or None."""
    cand = {}
    for i in np.nonzero(r["keep"])[0]:                  # as pileup_rates.clusters: every kept channel, truncated ones included
        k = int(r["k"][i])
        t, a = r["t"][i][:k], r["a"][i][:k]
        win = (t > pr.PRE_LO) & (t < pr.PRE_HI)
        if win.any():
            cand[int(r["mi"][i])] = float(a[win].max())
    if not cand:
        return None
    seed = max(cand, key=cand.get)
    return seed, [m for m in cand if nbr[seed, m]]


# ---------------------------------------------------------------------------------------------------------- geometry
def neighbours3x3(ctx, m):
    """3 x 3 grid [row (top = larger y)][col (left = smaller x)] of module indices around module m (None where no module)."""
    grid = []
    for dj in (1, 0, -1):
        row = []
        for di in (-1, 0, 1):
            px, py = ctx.x[m] + di * ctx.sx[m], ctx.y[m] + dj * ctx.sy[m]
            hit = np.nonzero((np.abs(ctx.x - px) <= ctx.sx / 2 + 0.01) & (np.abs(ctx.y - py) <= ctx.sy / 2 + 0.01))[0]
            row.append(int(hit[0]) if hit.size else None)
        grid.append(row)
    return grid


def rect(ctx, i, shrink=0.0):
    return Rectangle((ctx.x[i] - ctx.sx[i] / 2 + shrink, ctx.y[i] - ctx.sy[i] / 2 + shrink), ctx.sx[i] - 2 * shrink, ctx.sy[i] - 2 * shrink)


def text_color(rgba):
    return "white" if 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2] < 0.5 else "0.1"


def draw_modules(ax, ctx, idx, energy, norm, lw=0.4):
    """Module rectangles with outlines, colour = energy (nan: grey), dead modules hatched."""
    for i in idx:
        e = energy.get(int(i), np.nan)
        fc = CMAP(norm(e)) if np.isfinite(e) and e > 0 else (DEAD_GREY if ctx.dead[i] else GREY)
        p = rect(ctx, i)
        p.set(facecolor=fc, edgecolor="0.35", linewidth=lw, hatch="////" if ctx.dead[i] else None)
        ax.add_patch(p)


def draw_clusters(ax, ctx, clusters, labels=True, fs=9):
    for n, c in clusters:
        s = c["seed"]
        ax.add_patch(Rectangle((ctx.x[s] - ctx.sx[s] / 2, ctx.y[s] - ctx.sy[s] / 2), ctx.sx[s], ctx.sy[s], fill=False,
                               edgecolor="white", linewidth=2.6, zorder=5))
        ax.add_patch(Rectangle((ctx.x[c["seed"]] - ctx.sx[c["seed"]] / 2, ctx.y[c["seed"]] - ctx.sy[c["seed"]] / 2),
                               ctx.sx[c["seed"]], ctx.sy[c["seed"]], fill=False, edgecolor="black", linewidth=1.0, linestyle="--", zorder=6))
        ax.plot(c["x"], c["y"], marker="*", ms=14 if fs >= 9 else 10, mfc="#ffd400", mec="black", mew=1.0, ls="none", zorder=7)
        if labels:
            ax.annotate(f"C{n + 1}: {c['e']:.0f} MeV", (c["x"], c["y"]), xytext=(10, 10) if fs >= 9 else (14, -34), textcoords="offset points", fontsize=fs,
                        color="0.1", zorder=8, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.4", lw=0.6, alpha=0.9))


# --------------------------------------------------------------------------------------------------------------- page 1
def page_map(ctx, r, evnum, bits, m, clusters, title_extra):
    E = {int(mi): float(e) for mi, e in zip(r["mi"], r["emain"]) if np.isfinite(e)}
    ev = np.array(list(E.values()))
    fig = plt.figure(figsize=(18, 9.4))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.0, 0.03, 1.0, 0.03], wspace=0.68, left=0.05, right=0.965, top=0.88, bottom=0.085)
    ax, cax, az, cax2 = (fig.add_subplot(gs[j]) for j in range(4))
    hyc = [i for i in range(len(ctx.mods)) if ctx.mods[i].is_hycal()]       # LMS reference and veto modules are not HyCal channels
    vmin = max(ev.min() / 1.3, 1.0) if ev.size else 1.0
    norm = LogNorm(vmin=vmin, vmax=max(ev.max(), vmin * 10) if ev.size else 10)
    draw_modules(ax, ctx, hyc, E, norm)
    h = np.array(hyc)
    ax.set_xlim((ctx.x[h] - ctx.sx[h] / 2).min() - 5, (ctx.x[h] + ctx.sx[h] / 2).max() + 5)
    ax.set_ylim((ctx.y[h] - ctx.sy[h] / 2).min() - 5, (ctx.y[h] + ctx.sy[h] / 2).max() + 5)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.spines[["top", "right"]].set_visible(True)
    draw_clusters(ax, ctx, list(enumerate(clusters)))
    hl = rect(ctx, m)
    hl.set(fill=False, edgecolor="#d62728", linewidth=3.0, zorder=9)
    ax.add_patch(hl)
    ax.annotate(ctx.names[m], (ctx.x[m], ctx.y[m]), xytext=(30, 40) if ctx.x[m] < 0 else (-60, 40), textcoords="offset points", fontsize=13, fontweight="bold",
                color="#d62728", zorder=10, arrowprops=dict(arrowstyle="-|>", color="#d62728", lw=2),
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#d62728", lw=1.5))
    ax.set_title("HyCal, all modules (log colour scale)", loc="left")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=CMAP)
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("module energy E_main (MeV)  [log scale]")
    # zoom
    half = 4.5 * ctx.sx[m]
    x0, y0 = ctx.x[m], ctx.y[m]
    ov = [i for i in hyc if (abs(ctx.x[i] - x0) + 0.05 < half + ctx.sx[i] / 2) and (abs(ctx.y[i] - y0) + 0.05 < half + ctx.sy[i] / 2)]
    zin = np.array([E[i] for i in ov if i in E])
    znorm = Normalize(vmin=0, vmax=max(zin.max(), 1.0) if zin.size else 1.0)
    draw_modules(az, ctx, ov, E, znorm, lw=0.8)
    for i in ov:
        e = E.get(i, np.nan)
        fc = CMAP(znorm(e)) if np.isfinite(e) else (GREY)
        tc = text_color(fc) if np.isfinite(e) else "0.3"
        az.text(ctx.x[i], ctx.y[i] + ctx.sy[i] * 0.22, ctx.names[i], ha="center", va="center", fontsize=7.5, color=tc, clip_on=True, zorder=8)
        az.text(ctx.x[i], ctx.y[i] - ctx.sy[i] * 0.22, (f"{e:.0f}" if e >= 9.5 else f"{e:.1f}") if np.isfinite(e) else "-", ha="center", va="center",
                fontsize=11, fontweight="bold", color=tc, clip_on=True, zorder=9)
    az.set_xlim(x0 - half, x0 + half)
    az.set_ylim(y0 - half, y0 + half)
    az.set_aspect("equal")
    az.set_xlabel("x (mm)")
    az.set_ylabel("y (mm)")
    az.spines[["top", "right"]].set_visible(True)
    hl2 = rect(ctx, m)
    hl2.set(fill=False, edgecolor="#d62728", linewidth=4.0, zorder=9)
    az.add_patch(hl2)
    draw_clusters(az, ctx, [(n, c) for n, c in enumerate(clusters) if abs(c["x"] - x0) < half and abs(c["y"] - y0) < half], labels=True, fs=8.5)
    az.set_title(f"zoom: 9 x 9 modules around {ctx.names[m]} (linear, MeV per cell)", loc="left")
    cb2 = fig.colorbar(plt.cm.ScalarMappable(norm=znorm, cmap=CMAP), cax=cax2)
    cb2.set_label("module energy E_main (MeV)  [linear]")
    fig.suptitle(f"run {ctx.run}  event {evnum}  trigger bits 0x{bits:x}  beam {ctx.beam:g} MeV  highlighted {ctx.names[m]}\n"
                 f"colour = E_main = (largest fitted pulse amplitude, ADC) x gain (MeV/ADC), the pileup_rates chain definition", x=0.5, y=0.99, fontsize=14)
    fig.text(0.05, 0.012, "grey: no kept pulse (channel not read out [sparse readout], bad waveform, A < %g ADC, truncated main pulse or no calibration factor); hatched: dead module; red box: highlighted module;\n"
             "white/black box + star: seed and centroid of a reconstructed cluster (HyCalCluster, legacy, seed >= %g MeV); C<n>: cluster energy (hit energy = largest in-window pulse, 130-200 ns, which can differ from E_main).%s"
             % (pr.KEEP_ADC, MIN_CENTER_E, title_extra), fontsize=9.5, color="0.25")
    return fig


# --------------------------------------------------------------------------------------------------------------- page 2
def draw_wave(ax, ctx, r, i, hl=False, ylim=None, label=True):
    y, ped = r["Y"][i].astype(float), r["ped"][i]
    k = int(r["k"][i])
    ax.axhline(0, color="0.65", lw=0.8)
    ax.plot(TT, y - ped, "k.", ms=3.5, label="data" if label else None)
    if k:
        ax.plot(TT, r["comp"][i].sum(0), "r-", lw=1.5, label="fit (sum)" if label else None)
        for j in range(k):
            ax.plot(TT, r["comp"][i][j], "--", color=PULSE_COLORS[j], lw=1.3, label=f"pulse {j + 1}" if label else None)
    ax.set_xlim(0, NS * DT)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.tick_params(labelsize=9)


def panel_title(ctx, r, i, extra=""):
    m = int(r["mi"][i])
    rms, _, _ = pedestal_stats(r["Y"][i])
    s = f"{ctx.names[m]}  k={int(r['k'][i])}"
    if r["k"][i]:
        s += f"  chi2/dof={r['chi2'][i]:.1f}"
    s += f"\nped rms (clipped) {rms:.2f} ADC (sigma_ref {ctx.sig_ref[m]:.2f}), baseline {r['ped'][i]:.1f}" + (
        f" (30-sample {r['ped30'][i]:.1f})" if abs(r["ped30"][i] - r["ped"][i]) > 1.0 else "")
    e = r["emain"][i]
    s += f"  E_main {e:.0f} MeV" if np.isfinite(e) else ("  bad waveform" if r["bad"][i] else "  no energy")
    return s + extra


def page_neighbours(ctx, r, evnum, m, grid):
    fig, axs = plt.subplots(3, 3, figsize=(15, 11.5))
    for a in range(3):
        for b in range(3):
            ax, mm = axs[a, b], grid[a][b]
            if mm is None or mm not in r["row"]:
                ax.set_axis_off()
                ax.text(0.5, 0.5, "no module" if mm is None else f"{ctx.names[mm]}: not read out", ha="center", va="center", transform=ax.transAxes, color="0.4")
                continue
            i = r["row"][mm]
            draw_wave(ax, ctx, r, i, label=(mm == m))
            ax.set_title(panel_title(ctx, r, i), fontsize=10, color="#d62728" if mm == m else "0.12", loc="left")
            if mm == m:
                for s in ax.spines.values():
                    s.set_visible(True), s.set_color("#d62728"), s.set_linewidth(2.5)
                ax.legend(fontsize=8, loc="upper right", frameon=False)
            if a == 2:
                ax.set_xlabel("t (ns)")
            if b == 0:
                ax.set_ylabel("ADC (pedestal-subtracted)")
    fig.suptitle(f"run {ctx.run} event {evnum}: {ctx.names[m]} (red frame) and its 8 neighbours, geometric arrangement (x right, y up), 0-{NS * DT:.0f} ns",
                 y=0.995)
    fig.text(0.01, 0.005, "black: raw samples minus the fit baseline; red: total pile-up fit; dashed: fitted pulses in time order. Each panel has its own y range.\n"
             "Fit baseline = the chain's start value (30-sample pedestal, replaced by the module's b_ref when farther than 3 sigma_ref) plus the fitted offset (|offset| <= 1 sigma_ref).",
             fontsize=9.5, color="0.25")
    fig.tight_layout(rect=(0, 0.035, 1, 0.96))
    return fig


# --------------------------------------------------------------------------------------------------------------- page 3
def corr(a, b, mask):
    if mask.sum() < 10:
        return np.nan, int(mask.sum())
    a, b = a[mask] - a[mask].mean(), b[mask] - b[mask].mean()
    d = np.sqrt((a * a).sum() * (b * b).sum())
    return (float((a * b).sum() / d) if d > 0 else np.nan), int(mask.sum())


def fit_windows(r, i):
    """Mask [NS] of the fit windows of the fitted pulses of row i (all False without pulses)."""
    k = int(r["k"][i])
    if not k:
        return np.zeros(NS, bool)
    t0 = np.full((1, pr.KMAX), 0.0)
    t0[0, :k] = r["t"][i][:k]
    v = np.zeros((1, pr.KMAX), bool)
    v[0, :k] = True
    return pr.MF.windows(t0, v, NS)[0]


def coherence(ctx, r, m):
    """Same-board channels, common-mode estimates and the correlations of module m's waveform with them."""
    crate, slot, _ = ctx.daq[m]
    board = {c: ctx.by_daq.get((crate, slot, c)) for c in range(16)}
    ih = r["row"][m]
    yh = r["Y"][ih].astype(float) - r["ped"][ih]
    wh = fit_windows(r, ih)
    nop = (r["k"] == 0) & ~r["bad"] & (r["nclip"] == 0)
    nop[ih] = False
    P = r["Y"].astype(float) - pr.pedestal(r["Y"])[:, None]            # 30-sample pedestal for the no-pulse channels
    crates = np.array([ctx.daq[int(mi)][0] for mi in r["mi"]])
    slots = np.array([ctx.daq[int(mi)][1] for mi in r["mi"]])
    cm = dict(event=nop, crate=nop & (crates == crate), board=nop & (crates == crate) & (slots == slot))
    means = {k: (P[v].mean(0) if v.sum() else None) for k, v in cm.items()}
    counts = {k: int(v.sum()) for k, v in cm.items()}
    cors = []
    for c in range(16):
        mm = board[c]
        if mm is None or mm == m or mm not in r["row"]:
            continue
        j = r["row"][mm]
        yj = r["Y"][j].astype(float) - r["ped"][j]
        wj = fit_windows(r, j)
        cors.append(dict(c=c, name=str(ctx.names[mm]), k=int(r["k"][j]), out=corr(yh, yj, ~wh & ~wj), full=corr(yh, yj, np.ones(NS, bool))))
    cmc = {k: dict(out=corr(yh, means[k], ~wh), full=corr(yh, means[k], np.ones(NS, bool))) if means[k] is not None
           else dict(out=(np.nan, 0), full=(np.nan, 0)) for k in means}
    return dict(board=board, crate=crate, slot=slot, ih=ih, yh=yh, wh=wh, means=means, counts=counts, cors=cors, cmc=cmc)


def page_board(ctx, r, evnum, m, co):
    fig = plt.figure(figsize=(16, 13.5))
    gs = fig.add_gridspec(5, 4, height_ratios=[1, 1, 1, 1, 1.35], hspace=0.62, wspace=0.22, left=0.06, right=0.985, top=0.92, bottom=0.06)
    yh = co["yh"]
    lo, hi = yh.min(), yh.max()
    pad = 0.15 * (hi - lo + 1)
    ylim = (lo - pad, hi + pad)
    for c in range(16):
        ax = fig.add_subplot(gs[c // 4, c % 4])
        mm = co["board"][c]
        if mm is None or mm not in r["row"]:
            ax.set_axis_off()
            ax.text(0.5, 0.5, f"ch {c}: " + ("no module" if mm is None else f"{ctx.names[mm]} not read out"), ha="center", va="center", transform=ax.transAxes, color="0.4")
            continue
        i = r["row"][mm]
        draw_wave(ax, ctx, r, i, ylim=ylim, label=False)
        off = max(abs(r["Y"][i].astype(float) - r["ped"][i]).max() - max(abs(ylim[0]), abs(ylim[1])), 0)
        rms, _, _ = pedestal_stats(r["Y"][i])
        ax.set_title(f"ch {c} {ctx.names[mm]}  k={int(r['k'][i])}  rms {rms:.1f}" + ("  (pulse off scale)" if off > 0 else ""), fontsize=9.5,
                     color="#d62728" if mm == m else "0.12", loc="left")
        if mm == m:
            for s in ax.spines.values():
                s.set_visible(True), s.set_color("#d62728"), s.set_linewidth(2.5)
        if c // 4 == 3:
            ax.set_xlabel("t (ns)")
        if c % 4 == 0:
            ax.set_ylabel("ADC")
    ax = fig.add_subplot(gs[4, :2])
    ax.plot(TT, yh, "k.-", lw=0.8, ms=4, label=f"{ctx.names[m]} (pedestal-subtracted)")
    for k, col, lab in (("event", "#d62728", "event-wide"), ("crate", "#1f77b4", f"crate {co['crate']}"), ("board", "#2ca02c", f"board slot {co['slot']}")):
        if co["means"][k] is not None:
            ax.plot(TT, co["means"][k], color=col, lw=2.0, label=f"mean of no-pulse channels, {lab} (N={co['counts'][k]})")
        else:
            ax.plot([], [], color=col, lw=2.0, label=f"{lab}: no channel without a fitted pulse in the event (N=0)")
    ax.axhline(0, color="0.65", lw=0.8)
    ax.set_xlim(0, NS * DT)
    ax.set_xlabel("t (ns)")
    ax.set_ylabel("ADC")
    ax.set_title("common-mode estimate versus the highlighted channel", loc="left")
    ax.legend(fontsize=8.5, frameon=False, loc="best")
    ax = fig.add_subplot(gs[4, 2:])
    labels = [f"ch{x['c']} {x['name']}" for x in co["cors"]] + [f"CM {k}" for k in ("board", "crate", "event")]
    vals = [x["out"][0] for x in co["cors"]] + [co["cmc"][k]["out"][0] for k in ("board", "crate", "event")]
    valsf = [x["full"][0] for x in co["cors"]] + [co["cmc"][k]["full"][0] for k in ("board", "crate", "event")]
    xx = np.arange(len(labels))
    ax.bar(xx - 0.2, np.array(vals, float), 0.4, color="#4c78a8", label="outside the fitted windows of the pair")
    ax.bar(xx + 0.2, np.array(valsf, float), 0.4, color="#9ecae9", label="all 100 samples")
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xticks(xx)
    ax.set_xticklabels(labels, rotation=75, ha="right", fontsize=8)
    ax.set_ylim(-1, 1)
    if not co["cors"]:
        ax.text(0.5, 0.75, "no other channel of this board is read out in the event (sparse readout)", transform=ax.transAxes, ha="center", fontsize=10, color="0.3")
    if all(v is None for v in co["means"].values()):
        ax.text(0.5, 0.6, "no channel without a fitted pulse in the event: no common-mode estimate", transform=ax.transAxes, ha="center", fontsize=10, color="0.3")
    ax.set_ylabel(f"Pearson r with {ctx.names[m]}")
    ax.set_title("waveform correlation with the highlighted channel", loc="left")
    ax.legend(fontsize=8.5, frameon=False, loc="upper right")
    fig.suptitle(f"run {ctx.run} event {evnum}: same-board coherence of {ctx.names[m]} (crate {co['crate']}, slot {co['slot']}); all panels y = {ylim[0]:.0f} to {ylim[1]:.0f} ADC (red frame: {ctx.names[m]}), 0-{NS * DT:.0f} ns",
                 y=0.985, fontsize=13)
    return fig


# -------------------------------------------------------------------------------------------------------- persistence
def persistence(ctx, f, evnum, m, stride, window):
    """Baseline noise of module m and its 8 neighbours over a stride sample of EVIO file f plus a contiguous window around the event, in the
    events where the module was read out (the readout is sparse).  Per module: event numbers, std of samples 0-29 (r30) and the
    30-sample baseline minus the module's b_ref in sigma_ref (off; NaN without a b_ref); n_read / n_stride = readout fraction of the stride sample."""
    grid = neighbours3x3(ctx, m)
    nb = [x for row in grid for x in row if x is not None]
    ch = ctx.open(f)
    rows = {x: ([], []) for x in nb}
    n_phys = n_stride = 0
    n_read = dict.fromkeys(nb, 0)
    for n in physics_blocks(ch):
        for i in range(n):
            ch.select_event(i)
            en = int(ch.info().event_number)
            n_phys += 1
            on_stride = n_phys % stride == 0
            if not on_stride and abs(en - evnum) > window:
                continue
            d = ch.decode_event(i, with_ssp=False)
            if not d["ok"]:
                continue
            n_stride += on_stride
            for mi, smp in ctx.channels(d["event"]):
                if mi in rows:
                    rows[mi][0].append(en)
                    rows[mi][1].append(smp)
                    n_read[mi] += on_stride
    ch.close()
    out = {}
    for x in nb:
        if not rows[x][0]:
            out[x] = dict(ev=np.zeros(0, int), r30=np.zeros(0), off=np.zeros(0))
            continue
        Y = np.stack(rows[x][1]).astype(float)
        ped = pr.pedestal(Y)
        r30 = (Y[:, :pr.PED_NSAMPLES] - ped[:, None]).std(axis=1)
        out[x] = dict(ev=np.array(rows[x][0]), r30=r30, off=(ped - ctx.b_ref[x]) / ctx.sig_ref[x])
    return dict(mods=nb, d=out, n_phys=n_phys, n_stride=n_stride, n_read=n_read, stride=stride, window=window, evnum=evnum, m=m)


def page_persistence(ctx, ps):
    m, d = ps["m"], ps["d"]
    nb = [x for x in ps["mods"] if x != m]
    fig, axs = plt.subplots(2, 2, figsize=(15, 10))
    for ax, key, lab in ((axs[0, 0], "r30", "pedestal rms: std of samples 0-29 (ADC)"), (axs[0, 1], "off", "30-sample baseline minus b_ref, in sigma_ref")):
        fin = lambda v: v[np.isfinite(v)]    # noqa: E731  (off is NaN for modules without a b_ref in the file's noise reference)
        tf = (lambda v: np.sign(v) * np.log10(1.0 + np.abs(v))) if key == "off" else (lambda v: v)   # off: signed log axis (sigma_ref)
        pooled, own = tf(fin(np.concatenate([d[x][key] for x in nb]))), tf(fin(d[m][key]))
        allv = np.concatenate([pooled, own])
        bins = np.linspace(allv.min() - 0.02, allv.max() + 0.02, 70)
        ax.hist(pooled, bins=bins, histtype="step", lw=2, color="#4c78a8", density=True, label="8 neighbours pooled")
        ax.hist(own, bins=bins, histtype="step", lw=2.5, color="#d62728", density=True, label=ctx.names[m])
        t = np.nonzero(d[m]["ev"] == ps["evnum"])[0]
        if t.size:
            ax.axvline(tf(d[m][key][t[0]]), color="k", ls="--", lw=1.5, label=f"this event {d[m][key][t[0]]:.2f}")
        ax.set_yscale("log")
        if key == "off":
            tk = np.array([-3, -1, 0, 1, 3, 10, 30, 100, 300])
            ax.set_xticks(tf(tk.astype(float)))
            ax.set_xticklabels([str(v) for v in tk])
        ax.set_xlabel(lab)
        ax.set_ylabel("probability density (events read out)")
        ax.legend(frameon=False, fontsize=9)
    for ax, key, lab in ((axs[1, 0], "r30", "pedestal rms (ADC)"), (axs[1, 1], "off", "30-sample baseline minus b_ref (sigma_ref)")):
        for x in nb:
            ax.plot(d[x]["ev"], d[x][key], ".", ms=2, color="#4c78a8", alpha=0.5, label="8 neighbours" if x == nb[0] else None)
        ax.plot(d[m]["ev"], d[m][key], ".", ms=3, color="#d62728", label=ctx.names[m])
        ax.axvline(ps["evnum"], color="k", ls="--", lw=1, label="this event")
        if key == "off":
            ax.set_yscale("symlog", linthresh=3.0)
        ax.set_xlabel("event number")
        ax.set_ylabel(lab)
        ax.legend(frameon=False, fontsize=9, markerscale=4)
    fig.suptitle(f"run {ctx.run}: noise persistence of {ctx.names[m]} (events where it was read out): {'every ' + str(ps['stride']) + 'th physics event' if ps['stride'] > 1 else 'every physics event'} of the split "
                 f"\n({ps['n_stride']} decoded; {ps['n_phys']} physics events) + all events within +-{ps['window']} event numbers of {ps['evnum']}", fontsize=11, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fig


def persistence_text(ctx, ps):
    ctx_far = pr.START_NSIGMA
    m, d = ps["m"], ps["d"]
    nb = [x for x in ps["mods"] if x != m]
    ev = d[m]["ev"]
    L = [f"  persistence of {ctx.names[m]} baseline noise in its EVIO split: stride sample = {'every ' + str(ps['stride']) + 'th physics event' if ps['stride'] > 1 else 'every physics event'} ({ps['n_stride']} decoded of "
         f"{ps['n_phys']}), plus all events within +-{ps['window']} event numbers of the event; statistics are over the events where the channel was read out (sparse readout)",
         "  readout fraction in the stride sample (events with the channel present / decoded): "
         + ", ".join(f"{ctx.names[x]}{'*' if x == m else ''} {100 * ps['n_read'][x] / max(ps['n_stride'], 1):.2f} %" for x in ps["mods"])]
    tgt = np.nonzero(ev == ps["evnum"])[0]
    for key, name in (("r30", "pedestal rms (samples 0-29, ADC)"), ("off", "baseline offset (30-sample baseline - b_ref) / sigma_ref")):
        a = d[m][key]
        pooled = np.concatenate([d[x][key] for x in nb])
        pooled = pooled[np.isfinite(pooled)]
        if tgt.size:
            L.append(f"  {name}: this event {a[tgt[0]]:.2f}{' ADC' if key == 'r30' else ' sigma_ref'}; fraction of {ctx.names[m]}'s read-out events below it {100 * (a < a[tgt[0]]).mean():.2f} %, "
                     f"of the pooled neighbours' {100 * (pooled < a[tgt[0]]).mean():.3f} %")
        q = lambda v: ", ".join(f"{p:g}%: {np.percentile(v, p):.2f}" for p in (50, 90, 99, 99.9)) + f", max {v.max():.2f} (N={v.size})"   # noqa: E731
        L.append(f"    {ctx.names[m]}: {q(a)}")
        L.append(f"    neighbours pooled: {q(pooled)}")
        for x in nb:
            v = d[x][key][np.isfinite(d[x][key])]
            L.append(f"      {ctx.names[x]}: " + (f"median {np.median(v):.2f}, 99 % {np.percentile(v, 99):.2f}, max {v.max():.2f} (N={v.size})" if v.size else "no valid b_ref in this file's noise reference"))
        thr = 3.0 * np.median(pooled) if key == "r30" else ctx_far
        far = np.abs(ev - ps["evnum"]) > ps["window"]
        hi = (a > thr) if key == "r30" else (np.abs(a) > thr)
        L.append(f"    {ctx.names[m]} read-out events with " + (f"rms above 3 x the pooled-neighbour median ({thr:.2f} ADC)" if key == "r30" else f"|offset| > {thr:g} sigma_ref (the chain's 'far' criterion, baseline start replaced by b_ref)") + f": {int(hi.sum())} of {a.size} ({100 * hi.mean():.2f} %); "
                 f"inside the +-{ps['window']} window {int((hi & ~far).sum())} of {int((~far).sum())}; outside it {int((hi & far).sum())} of {int(far.sum())} "
                 f"({100 * (hi & far).sum() / max(far.sum(), 1):.2f} %)")
        if hi.any():
            L.append(f"    those events (event number: value): " + ", ".join(f"{e}: {v:.1f}" for e, v in zip(ev[hi], a[hi])))
        for x in nb:
            dd = d[x][key]
            hh = (dd > thr) if key == "r30" else (np.abs(dd) > thr)
            L.append(f"      {ctx.names[x]}: " + (f"{int(hh.sum())} of {int(np.isfinite(dd).sum())} read-out events above the same threshold" if np.isfinite(dd).any() else "no valid b_ref in this file's noise reference") + (f" (event numbers {', '.join(str(e) for e in d[x]['ev'][hh])})" if hh.any() else ""))
    return L


# --------------------------------------------------------------------------------------------------------------- summary
def summary_text(ctx, r, evnum, bits, m, clusters, chain, co, grid, fpath, split, extra):
    i = r["row"][m]
    rms, slope, rms_raw = pedestal_stats(r["Y"][i])
    k = int(r["k"][i])
    L = [f"=== run {ctx.run} event {evnum} module {ctx.names[m]} ===",
         f"EVIO file {Path(fpath).name} (split ordinal {split}); trigger bits 0x{bits:x}; beam {ctx.beam:g} MeV; run type {ctx.rt}; calibration {ctx.cfg['calib']}",
         f"{len(r['mi'])} HyCal channels in the event, {int(r['bad'].sum())} bad waveforms, {int(r['keep'].sum())} kept (A_main >= {pr.KEEP_ADC:g} ADC)",
         f"{ctx.names[m]}: DAQ crate {co['crate']} slot {co['slot']} channel {ctx.daq[m][2]}; gain {ctx.gain[m]:.4f} MeV/ADC; sigma_ref {ctx.sig_ref[m]:.3f} ADC; "
         f"b_ref {ctx.b_ref[m]:.2f}",
         f"  E_main {r['emain'][i]:.1f} MeV (A_main {r['amain'][i]:.1f} ADC x gain)" if np.isfinite(r["emain"][i]) else
         f"  E_main: none (bad waveform {bool(r['bad'][i])}, A_main {r['amain'][i]:.1f} ADC, truncated main {bool(r['trunc'][i])})",
         (f"  baseline: 30-sample pedestal {r['ped30'][i]:.2f} ADC, b_ref " + (
             f"{ctx.b_ref[m]:.2f} ({(r['ped30'][i] - ctx.b_ref[m]) / ctx.sig_ref[m]:+.1f} sigma_ref; "
             + ("farther than" if abs(r['ped30'][i] - ctx.b_ref[m]) > pr.START_NSIGMA * ctx.sig_ref[m] else "within")
             + f" {pr.START_NSIGMA:g} sigma_ref, so the fit starts at "
             + ("b_ref" if abs(r['ped30'][i] - ctx.b_ref[m]) > pr.START_NSIGMA * ctx.sig_ref[m] else "the 30-sample pedestal") + ")"
             if np.isfinite(ctx.b_ref[m]) else "none in this file's noise reference (no baseline override)")
          + f"; fit baseline {r['ped'][i]:.2f} ADC; rms (clipped samples 0-29) {rms:.3f} ADC (raw std {rms_raw:.3f}); slope {slope:+.3f} ADC per 100 ns"),
         f"  fitted pulses k = {k}, chi2/dof {r['chi2'][i]:.2f}, clipped samples {int(r['nclip'][i])}:"]
    for j in range(k):
        L.append(f"    pulse {j + 1}: t0 = {r['t'][i][j]:.1f} ns, A = {r['a'][i][j]:.1f} ADC (E = {r['a'][i][j] * ctx.gain[m]:.1f} MeV)")
    c = next((c for c in clusters if m in c["mi"]), None)
    if c is None:
        L.append(f"  cluster (HyCalCluster, legacy, seed >= {MIN_CENTER_E:g} MeV): not in a cluster; clusters of the event: "
                 + (", ".join(f"[seed {ctx.names[q['seed']]} {q['e']:.0f} MeV size {q['n']}]" for q in clusters) or "none"))
    else:
        L.append(f"  cluster: seed {ctx.names[c['seed']]} (E {c['seed_e']:.1f} MeV), E_cluster {c['e']:.1f} MeV (sum of members {c['esum']:.1f}), size {c['n']} modules, "
                 f"E/E_beam {c['e'] / ctx.beam:.3f}, position ({c['x']:.1f}, {c['y']:.1f}) mm; {ctx.names[m]} has {c['me'][c['mi'].index(m)]:.1f} MeV "
                 f"({100 * c['me'][c['mi'].index(m)] / c['esum']:.0f} % of the sum); seed distance {np.hypot(ctx.x[m] - ctx.x[c['seed']], ctx.y[m] - ctx.y[c['seed']]):.1f} mm")
        L.append("  all clusters: " + ", ".join(f"[seed {ctx.names[q['seed']]} E {q['e']:.0f} MeV size {q['n']} E/Eb {q['e'] / ctx.beam:.2f} at ({q['x']:.0f}, {q['y']:.0f}) mm]" for q in clusters))
    if chain is None:
        L.append("  pileup_rates chain cluster (seed + members within 1.5 x size, in-window pulses): none")
    else:
        L.append(f"  pileup_rates chain cluster: seed {ctx.names[chain[0]]}; {ctx.names[m]} is " + ("a member" if m in chain[1] else "not a member") + f" ({len(chain[1])} members)")
    L.append("  neighbours (row top to bottom = y larger to smaller; col left to right): E_main (MeV) | k | ped rms")
    for a in range(3):
        cells = []
        for b in range(3):
            mm = grid[a][b]
            if mm is None or mm not in r["row"]:
                cells.append("---")
                continue
            j = r["row"][mm]
            e = r["emain"][j]
            cells.append(f"{ctx.names[mm]}{'*' if mm == m else ''} {e:.1f} | k{int(r['k'][j])} | {pedestal_stats(r['Y'][j])[0]:.2f}" if np.isfinite(e)
                         else f"{ctx.names[mm]}{'*' if mm == m else ''} - | k{int(r['k'][j])} | {pedestal_stats(r['Y'][j])[0]:.2f}")
        L.append("    " + "   ".join(f"{c:<26}" for c in cells))
    # pulses of neighbours
    L.append("  neighbour pulses (t0 ns : A ADC):")
    for a in range(3):
        for b in range(3):
            mm = grid[a][b]
            if mm is None or mm == m or mm not in r["row"]:
                continue
            j = r["row"][mm]
            kk = int(r["k"][j])
            L.append(f"    {ctx.names[mm]}: " + ("; ".join(f"{r['t'][j][q]:.1f}:{r['a'][j][q]:.1f}" for q in range(kk)) if kk else "no pulse"))
    L.append(f"  coherence (page 3): fADC board crate {co['crate']} slot {co['slot']}; no-pulse channels in the mean: event-wide {co['counts']['event']}, crate {co['counts']['crate']}, board {co['counts']['board']}")
    L.append(f"  highlighted fit windows cover {int(co['wh'].sum())} of {NS} samples; r = Pearson correlation of pedestal-subtracted waveforms; "
             f"'out' = samples outside both fitted-pulse windows (n), 'all' = all {NS} samples")
    for x in co["cors"]:
        L.append(f"    ch{x['c']:>2} {x['name']:<5} k={x['k']}: r_out {x['out'][0]:+.3f} (n={x['out'][1]}), r_all {x['full'][0]:+.3f}")
    for kname in ("board", "crate", "event"):
        v = co["cmc"][kname]
        L.append(f"    common mode ({kname}, N={co['counts'][kname]}): r_out {v['out'][0]:+.3f} (n={v['out'][1]}), r_all {v['full'][0]:+.3f}; mean-waveform rms " + (f"{np.std(co['means'][kname]):.3f} ADC" if co['means'][kname] is not None else "n/a"))
    L.extend(extra)
    return L


# ------------------------------------------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", type=int, required=True)
    ap.add_argument("--event", type=int, nargs="+", required=True, help="EVIO event number(s) (fe.info.event_number)")
    ap.add_argument("--module", nargs="+", required=True, help="module name to highlight, one per event (or one for all)")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--evio-dir", default="/Users/L00338853/work/PRad/data/evio")
    ap.add_argument("--persistence", nargs="*", default=[], metavar="MODULE", help="modules (among --module) for which page 4 and the noise statistics over the split are made")
    ap.add_argument("--stride", type=int, default=1, help="persistence: decode every Nth physics event of the split")
    ap.add_argument("--window", type=int, default=5000, help="persistence: also decode all events within this many event numbers of the event")
    ap.add_argument("--persistence-cache", default=None, metavar="FILE.pkl", help="pickle of the persistence pass: reused if it exists and matches, else written")
    ap.add_argument("--keep-png", action="store_true", help="also save one PNG per page")
    a = ap.parse_args()
    mods = a.module * len(a.event) if len(a.module) == 1 else a.module
    if len(mods) != len(a.event):
        sys.exit("[ERROR] give one --module per --event (or a single one)")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ctx = Ctx(a.run)
    nbr = pr.neighbours(ctx.names)
    refs = {}
    summ = []
    for evnum, mname in zip(a.event, mods):
        t0 = time.monotonic()
        if mname not in ctx.index:
            sys.exit(f"[ERROR] module {mname} not in the HyCal map")
        m = ctx.index[mname]
        fpath, split, files = locate(ctx, a.evio_dir, evnum)
        print(f"[locate] event {evnum} in {Path(fpath).name} ({time.monotonic() - t0:.1f} s)", flush=True)
        fe = find_event(ctx, fpath, evnum)
        if fe is None:
            sys.exit(f"[ERROR] event {evnum} not found in {fpath}")
        bits = int(fe.info.trigger_bits)
        if fpath not in refs:
            refs[fpath] = noise_reference(ctx, fpath)
            print(f"[noise reference] {Path(fpath).name} ({time.monotonic() - t0:.1f} s)", flush=True)
        _, b_ref, sig_ref = refs[fpath]
        ctx.b_ref, ctx.sig_ref = b_ref, sig_ref
        r = fit_event(ctx, fe, b_ref, sig_ref)
        if m not in r["row"]:
            sys.exit(f"[ERROR] module {mname} has no channel with {NS} samples in event {evnum}")
        clusters, _ = cluster_event(ctx, r)
        chain = chain_cluster(ctx, r, nbr)
        grid = neighbours3x3(ctx, m)
        co = coherence(ctx, r, m)
        pdf_path = out / f"event_display_run{a.run}_ev{evnum}_{mname}.pdf"
        extra = []
        figs = [page_map(ctx, r, evnum, bits, m, clusters, ""), page_neighbours(ctx, r, evnum, m, grid), page_board(ctx, r, evnum, m, co)]
        if mname in a.persistence:
            ps = None
            key = (a.run, evnum, mname, a.stride, a.window)
            if a.persistence_cache and Path(a.persistence_cache).exists():
                with open(a.persistence_cache, "rb") as fh:
                    cached = pickle.load(fh)
                ps = cached["ps"] if cached["key"] == key else None
            if ps is None:
                ps = persistence(ctx, fpath, evnum, m, a.stride, a.window)
                if a.persistence_cache:
                    with open(a.persistence_cache, "wb") as fh:
                        pickle.dump(dict(key=key, ps=ps), fh)
            figs.append(page_persistence(ctx, ps))
            extra = persistence_text(ctx, ps)
            print(f"[persistence] {ps['n_stride']} stride events decoded ({time.monotonic() - t0:.1f} s)", flush=True)
        with PdfPages(pdf_path) as pdf:
            for n, fig in enumerate(figs):
                pdf.savefig(fig)
                if a.keep_png:
                    fig.savefig(out / f"{pdf_path.stem}_p{n + 1}.png", dpi=110)
                plt.close(fig)
        summ += summary_text(ctx, r, evnum, bits, m, clusters, chain, co, grid, fpath, split, extra) + [""]
        print(f"[done] {pdf_path} ({time.monotonic() - t0:.0f} s)", flush=True)
    (out / "summary.txt").write_text("\n".join(summ) + "\n")


if __name__ == "__main__":
    main()
