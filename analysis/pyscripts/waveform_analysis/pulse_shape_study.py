#!/usr/bin/env python3
"""pulse_shape_study.py -- smooth vs bumpy falling tails of good PbWO4 pulses.

    python3 waveform_analysis/pulse_shape_study.py scan <out.npz> <evio file> [<evio file> ...]
    python3 waveform_analysis/pulse_shape_study.py report <scan dir> <out dir>

scan    all EVIO split files of one run -> npz of every good pulse (module, event, amplitude,
        peak, descent length, bump size) plus a random reservoir of example waveforms per
        class (A >= AMP_MIN only).
report  pool the <scan dir>/run_*.npz files -> examples.pdf, hot_examples.pdf,
        bumpy_rate_map.png/.pdf, bumpy_rate_hist.png, modules.tsv, hot_modules.tsv, summary.txt
        in <out dir>.  Everything derived from the pulses keeps A >= AMP_MIN; modules.tsv also
        keeps the all-A columns.

Good pulse = the production `template_extract.py select` cuts (trigger gate of the run
type, LMS bit veto, <= 20 modules/event, clip veto, A >= 50, peak <= 59, first_thr >= 20,
pedestal, tail baseline, no secondary bump >= 10 %, e_rough < 1.2 E_beam) except its smoothness
cut c10 (monotonic rise, no tail sample up by more than 2 rms): the tail shape is what this
study measures.
Tail classes (w = samples - WaveAnalyzer pedestal; tail = the TAIL_NSAMP samples after the peak,
peak+1 .. peak+TAIL_NSAMP, cut earlier at the first sample with w < TAIL_END_ADC (back at pedestal
level)):
  bumpy   = a dip sample k in the tail with w[k] < w[k-1], followed within the tail by a rise
            to a secondary maximum (also in the tail) >= BUMP_MIN_ADC above w[k].  A plateau
            followed by continued descent is not a bump; a flat run followed by a rise >=
            BUMP_MIN_ADC is (the dip is the first sample of the flat bottom).  With several
            qualifying dips the one with the largest rise is recorded.
  plateau = not bumpy, but the strictly descending run after the peak (PASS 1 `falling`)
            breaks within TAIL_NSAMP samples.
  smooth  = otherwise.
Required rise = max(BUMP_MIN_ADC, BUMP_MIN_FRAC * A).  BUMP_MIN_ADC = 3 ADC is about 3 sigma
of the ~0.9 ADC pedestal noise, but alone it is noise dominated: on split 25308/50, 74 % of
the A >= 200 pulses were "bumpy" with the secondary maximum 50 samples after the peak
(tail jitter of 1-5 ADC grows with A, and the rise was a max over a ~60 sample tail, before the
TAIL_NSAMP window).
Requiring 3 % of A as well leaves 0.13 %; 1 %, 2 % are still inside the noise tail.
Rate = bumpy / all good pulses with A >= AMP_MIN; plateau is reported separately.
scan also keeps per module and run up to K_HOT random smooth and bumpy waveforms (A >=
AMP_MIN) for the hot-module examples.
"""
from __future__ import annotations

import glob
import os
import sys
import time
from pathlib import Path

import numpy as np

import template_extract as te

TAIL_NSAMP = 8                         # samples after the peak: the tail window (dip, secondary max, plateau test)
BUMP_MIN_ADC = 3.0                     # ADC: minimum rise above the dip (~3 sigma of the pedestal noise)
BUMP_MIN_FRAC = 0.03                   # ... and at least this fraction of A (tail noise scales with A)
TAIL_END_ADC = 2.0                     # ADC: the tail ends at the first sample below this (pedestal level)
K_HOT = 4                              # per module and class example waveforms kept for hot modules
MIN_PULSES_MAP = 20                    # heat map: PbWO4 modules with fewer pulses at A >= AMP_MIN are grey
HOT_RATE = 0.10                        # hot module: bumpy rate (A >= AMP_MIN) above this ...
HOT_MIN_N = 20                         # ... with at least this many pulses
SMOOTH, PLATEAU, BUMPY = 0, 1, 2
CLS_NAME = {SMOOTH: "smooth", PLATEAU: "plateau", BUMPY: "bumpy"}
N_RESERVOIR = 200                      # example waveforms kept per class and run
MIN_N_HIST = 20                        # modules with >= this many good pulses in the histogram
SEED = 1                               # random seed of the example reservoirs
N_EXAMPLES = 50                        # examples per class in examples.pdf (<= N_RESERVOIR)
AMP_MIN = 200.0                        # ADC: report and example reservoir keep A >= AMP_MIN; scan records keep all
# run -> template_extract run type (data/evio/README.md; best50 pulses/ directories)
RUN_TYPE_OF = {
    24590: "0p73GeV_H2", 24603: "0p73GeV_H2",
    25228: "2p24GeV_H2", 25257: "2p24GeV_H2", 25260: "2p24GeV_H2",
    25299: "2p24GeV_H2", 25320: "2p24GeV_H2",
    25308: "2p24GeV_C",
    25800: "2p24GeV_Ta_X17", 26138: "2p24GeV_Ta_X17", 26230: "2p24GeV_Ta_X17",
    24698: "3p49GeV_H2", 24911: "3p49GeV_H2", 24937: "3p49GeV_H2",
    # previous 3p49GeV_H2 runs (replaced 2026-11-10)
    # 24386: "3p49GeV_H2", 24503: "3p49GeV_H2", 24512: "3p49GeV_H2",
}
REC = dict(midx=np.int16, event_number=np.int64, amp=np.float32, peak=np.int16,
           falling=np.int16, cls=np.int8, dip=np.int16, bump_adc=np.float32,
           bump=np.float32, gap=np.int16)


def classify_tails(W, peak, falling):
    """Tail class, dip sample, bump height (ADC and / A) and peak-to-secondary-maximum gap
    (samples) of each row of W (pedestal-subtracted); dip = -1 and bump = gap = 0 when not
    bumpy.  Dip and secondary maximum lie in the tail window peak+1 .. peak+TAIL_NSAMP."""
    n, ns = W.shape
    rows = np.arange(n)
    idx = np.arange(ns)[None, :]
    below = (idx > peak[:, None]) & (W < TAIL_END_ADC)
    end = np.where(below.any(1), below.argmax(1), ns)          # first sample back at pedestal
    end = np.minimum(end, peak + TAIL_NSAMP + 1)               # tail window: peak+1 .. peak+TAIL_NSAMP
    inside = idx < end[:, None]
    sufmax = np.maximum.accumulate(np.where(inside, W, -np.inf)[:, ::-1], axis=1)[:, ::-1]
    nxt = np.concatenate([sufmax[:, 1:], np.full((n, 1), -np.inf)], axis=1)   # max w[k+1:end]
    rise = nxt - W
    prev = np.concatenate([W[:, :1], W[:, :-1]], axis=1)
    cand = inside & (idx > peak[:, None]) & (W < prev) & (rise >= np.maximum(BUMP_MIN_ADC, BUMP_MIN_FRAC * W[rows, peak])[:, None])
    score = np.where(cand, rise, -1.0)
    dip = score.argmax(1)
    bumpy = cand.any(1)
    h = np.where(bumpy, score[rows, dip], 0.0)
    sec = np.where(inside & (idx > dip[:, None]), W, -np.inf).argmax(1)
    gap = np.where(bumpy, sec - peak, 0)
    cls = np.where(bumpy, BUMPY, np.where(falling < TAIL_NSAMP, PLATEAU, SMOOTH))
    return (cls.astype(np.int8), np.where(bumpy, dip, -1).astype(np.int16), h,
            h / W[rows, peak], gap.astype(np.int16))


class StudySelector(te.Selector):
    """template_extract's Selector without the smoothness cut, with every survivor recorded as a
    per-pulse record and classified by its falling tail."""
    smooth_cut = False

    def __init__(self, *a, seed=0, **k):
        super().__init__(*a, **k)
        self.rec_out = {f: [] for f in REC}
        self.ex = None                  # global reservoir: dict of arrays incl. key and cls
        self.hx = None                  # per-module reservoir (hot-module examples)
        self.rng = np.random.default_rng(seed)
        self.n_good = 0
        self.split = 0                  # split number of the file being scanned

    def kept(self):
        return self.n_good

    def _keep(self, f, fields):
        d = fields(f)
        n = f.size
        if n == 0:
            return
        self.n_good += n
        W = d["samples"].astype(np.float64) - d["ped"][:, None]
        cls, dip, bump_adc, bump, gap = classify_tails(W, d["peak_sample"], d["falling"])
        r = dict(midx=d["midx"], event_number=d["event_number"], amp=d["amp"],
                 peak=d["peak_sample"], falling=d["falling"], cls=cls, dip=dip,
                 bump_adc=bump_adc, bump=bump, gap=gap)
        for k, dt in REC.items():
            self.rec_out[k].append(np.asarray(r[k]).astype(dt))
        new = dict(r, samples=d["samples"], ped=d["ped"], key=self.rng.random(n),
                   split=np.full(n, self.split, np.int16))
        self._reservoir(new)
        self._hot_reservoir(new)

    def _hot_reservoir(self, new):
        """Per module and class (smooth, bumpy): the K_HOT smallest random keys among
        pulses with A >= AMP_MIN."""
        ns = new["samples"].shape[1]
        if self.hx is None:
            nm = self.nm
            self.hx = dict(key=np.full((nm, 2, K_HOT), np.inf),
                           samples=np.zeros((nm, 2, K_HOT, ns), np.uint16),
                           ped=np.zeros((nm, 2, K_HOT), np.float32),
                           event_number=np.zeros((nm, 2, K_HOT), np.int64),
                           amp=np.zeros((nm, 2, K_HOT), np.float32),
                           peak=np.zeros((nm, 2, K_HOT), np.int16),
                           dip=np.zeros((nm, 2, K_HOT), np.int16),
                           gap=np.zeros((nm, 2, K_HOT), np.int16),
                           split=np.zeros((nm, 2, K_HOT), np.int16))
        hx = self.hx
        c = (new["cls"] == BUMPY).astype(int)
        ok = (new["cls"] != PLATEAU) & (new["amp"] >= AMP_MIN)
        ok &= new["key"] < hx["key"][new["midx"], c].max(axis=1)
        for i in np.nonzero(ok)[0]:
            m, cc = int(new["midx"][i]), int(c[i])
            s = int(hx["key"][m, cc].argmax())
            if new["key"][i] >= hx["key"][m, cc, s]:
                continue
            hx["key"][m, cc, s] = new["key"][i]
            hx["samples"][m, cc, s] = new["samples"][i]
            for k in ("ped", "event_number", "amp", "peak", "dip", "gap", "split"):
                hx[k][m, cc, s] = new[k][i]

    def _reservoir(self, new):
        """Keep the N_RESERVOIR smallest random keys per class (smooth, bumpy) among pulses with
        A >= AMP_MIN (a uniform sample of that population)."""
        if self.ex is not None:
            new = {k: np.concatenate([self.ex[k], new[k]]) for k in new}
        keep = np.zeros(new["key"].size, dtype=bool)
        for cls in (SMOOTH, BUMPY):
            ix = np.nonzero((new["cls"] == cls) & (new["amp"] >= AMP_MIN))[0]
            keep[ix[np.argsort(new["key"][ix])[:N_RESERVOIR]]] = True
        self.ex = {k: v[keep] for k, v in new.items()}

    def records(self):
        return {k: (np.concatenate(v) if v else np.zeros(0, REC[k]))
                for k, v in self.rec_out.items()}


def scan(out, evios):
    """`scan`: all split files of one run -> per-pulse records + example reservoir npz."""
    t_start = time.monotonic()
    found = []
    for evio in evios:
        mt = te._SPLIT_PAT.search(os.path.basename(evio))
        if not mt:
            sys.exit(f"[ERROR] cannot read run/split from {evio}")
        found.append((int(mt.group(1)), int(mt.group(2))))
    run = found[0][0]
    if any(r != run for r, _ in found):
        sys.exit(f"[ERROR] the files are not all of run {run}")
    if run not in RUN_TYPE_OF:
        sys.exit(f"[ERROR] run {run} has no run type in RUN_TYPE_OF")
    rt_name = RUN_TYPE_OF[run]

    # the same setup as template_extract.select_candidates
    S = te.select_setup(run, rt_name)
    p, mods, gate = S.p, S.mods, S.gate
    splits = [s for _, s in found]
    print(f"[setup] run {run} splits {splits} run type {rt_name}; trigger {S.trig_rule}",
          flush=True)

    sel = StudySelector(p, S.cal, S.emax_mev, seed=SEED * 1000003 + run * 100000)
    n_seen = n_trig = n_excl = 0
    for evio, split in zip(evios, splits):
        sel.split = split
        counts = te._scan_events(evio, p, sel, gate, t_start)
        n_seen, n_trig, n_excl = (a + b for a, b in zip((n_seen, n_trig, n_excl), counts))
    wall = time.monotonic() - t_start
    r = sel.records()
    ex = sel.ex or {}
    te._save_npz(out, dict(
        run=np.int64(run), splits=np.array(splits, np.int64), run_type=np.array(rt_name),
        n_events_seen=np.int64(n_seen), n_events_trig=np.int64(n_trig),
        n_events_excl=np.int64(n_excl),
        wall_s=np.float64(wall),
        mod_name=np.array([m.name for m in mods]),
        mod_x=np.array([m.x for m in mods]), mod_y=np.array([m.y for m in mods]),
        mod_sx=np.array([m.size_x for m in mods]), mod_sy=np.array([m.size_y for m in mods]),
        mod_pwo4=np.array([bool(m.is_pwo4()) for m in mods]),
        **r, **{f"ex_{k}": v for k, v in ex.items()},
        **{f"hx_{k}": v for k, v in (sel.hx or {}).items()}))
    nb = int((r["cls"] == BUMPY).sum())
    npl = int((r["cls"] == PLATEAU).sum())
    print(f"[done] {out}: {n_seen} events, {r['midx'].size} good pulses, {nb} bumpy "
          f"({nb / max(r['midx'].size, 1):.2%}), {npl} plateau; {wall:.0f} s", flush=True)


# report ======================================================================
def load_scans(scan_dir):
    files = sorted(glob.glob(os.path.join(scan_dir, "run_*.npz")))
    if not files:
        sys.exit(f"[ERROR] no run_*.npz files in {scan_dir}")
    zs = [np.load(f, allow_pickle=False) for f in files]
    geo = {k: zs[0][k] for k in ("mod_name", "mod_x", "mod_y", "mod_sx", "mod_sy", "mod_pwo4")}
    if any(not np.array_equal(z["mod_name"], geo["mod_name"]) for z in zs):
        sys.exit("[ERROR] the scan files have different module lists")
    rec = {k: np.concatenate([z[k] for z in zs]) for k in REC}
    ex = {}
    for z in zs:
        if "ex_key" not in z.files:
            continue
        n = z["ex_key"].size
        part = {k[3:]: z[k] for k in z.files if k.startswith("ex_")}
        part["run"] = np.full(n, int(z["run"]))
        for k, v in part.items():
            ex.setdefault(k, []).append(v)
    ex = {k: np.concatenate(v) for k, v in ex.items()}
    hx = {}
    for z in zs:
        if "hx_key" not in z.files:
            continue
        part = {k[3:]: z[k] for k in z.files if k.startswith("hx_")}
        part["run"] = np.full(part["key"].shape, int(z["run"]))
        for k, v in part.items():
            hx.setdefault(k, []).append(v)
    hx = {k: np.stack(v) for k, v in hx.items()}      # [run, module, class, K(, ns)]
    ex["hx"] = hx
    tot = {k: sum(int(z[k]) for z in zs) for k in
           ("n_events_seen", "n_events_trig", "n_events_excl")}
    tot["wall_s"] = [float(z["wall_s"]) for z in zs]
    tot["runs"] = sorted({int(z["run"]) for z in zs})
    tot["n_splits"] = sum(int(z["splits"].size) for z in zs)
    return files, geo, rec, ex, tot


def mark_pulse(ax, w, pk, dip, gap):
    ax.plot(np.arange(w.size), w, "k.-", lw=0.7, ms=2.5)
    ax.plot(pk, w[pk], "o", color="tab:blue", ms=5)
    if dip >= 0:
        ax.plot(dip, w[dip], "v", color="tab:red", ms=6)
        ax.plot(pk + gap, w[pk + gap], "^", color="tab:orange", ms=6)


def plot_hot_examples(hx, hot, geo, path):
    """One row per hot module: 2 smooth then 2 bumpy pulses (random across runs)."""
    import matplotlib.pyplot as plt
    # plot style: large bold axis labels, large legends
    plt.rcParams.update({
        "axes.labelsize": 16, "axes.labelweight": "bold", "axes.titlesize": 16, "axes.titleweight": "bold",
        "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14, "axes.linewidth": 1.8,
        "xtick.major.width": 1.5, "ytick.major.width": 1.5, "xtick.major.size": 6, "ytick.major.size": 6,
        "lines.linewidth": 2.0, "figure.titlesize": 18, "figure.titleweight": "bold"})
    from matplotlib.backends.backend_pdf import PdfPages
    rng = np.random.default_rng(12345)
    per_page = 6
    with PdfPages(path) as pdf:
        for p0 in range(0, max(len(hot), 1), per_page):
            fig, axs = plt.subplots(per_page, 4, figsize=(22, 22))
            for r, ax_row in enumerate(axs):
                m = hot[p0 + r] if p0 + r < len(hot) else None
                picks = {}
                for cc in (0, 1):
                    if m is None:
                        break
                    s_ix, k_ix = np.nonzero(np.isfinite(hx["key"][:, m, cc]))
                    picks[cc] = (s_ix, k_ix, rng.permutation(s_ix.size)[:2])
                for c, ax in enumerate(ax_row):
                    if m is None:
                        ax.axis("off")
                        continue
                    cc = c // 2
                    s_ix, k_ix, pick = picks[cc]
                    if (c % 2) >= pick.size:
                        ax.axis("off")
                        continue
                    j = pick[c % 2]
                    s, k = s_ix[j], k_ix[j]
                    w = hx["samples"][s, m, cc, k].astype(float) - hx["ped"][s, m, cc, k]
                    pk = int(hx["peak"][s, m, cc, k])
                    mark_pulse(ax, w, pk, int(hx["dip"][s, m, cc, k]), int(hx["gap"][s, m, cc, k]))
                    ax.set_xlim(max(pk - 10, 0), w.size - 1)
                    ax.tick_params(labelsize=11)
                    ax.set_title(f"{geo['mod_name'][m]} {'bumpy' if cc else 'smooth'}  A={hx['amp'][s, m, cc, k]:.0f}\n"
                                 f"run {hx['run'][s, m, cc, k]} split {hx['split'][s, m, cc, k]} "
                                 f"ev {hx['event_number'][s, m, cc, k]}", fontsize=11)
            fig.suptitle("hot modules: 2 smooth + 2 bumpy pulses each (blue = peak, red = dip, "
                         "orange = secondary max)")
            fig.tight_layout(rect=(0, 0, 1, 0.97))
            pdf.savefig(fig)
            plt.close(fig)


def plot_examples(ex, geo, n_ex, path):
    """n_ex smooth then n_ex bumpy examples (the smallest random keys = uniform sample of
    every scanned good pulse of that class with A >= AMP_MIN), 10 per page."""
    import matplotlib.pyplot as plt
    # plot style: large bold axis labels, large legends
    plt.rcParams.update({
        "axes.labelsize": 16, "axes.labelweight": "bold", "axes.titlesize": 16, "axes.titleweight": "bold",
        "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14, "axes.linewidth": 1.8,
        "xtick.major.width": 1.5, "ytick.major.width": 1.5, "xtick.major.size": 6, "ytick.major.size": 6,
        "lines.linewidth": 2.0, "figure.titlesize": 18, "figure.titleweight": "bold"})
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(path) as pdf:
        for cls, label in ((SMOOTH, "smooth"), (BUMPY, "bumpy")):
            keep = (ex["cls"] == cls) & (ex["amp"] >= AMP_MIN) if "cls" in ex else np.zeros(0, bool)
            ix = np.nonzero(keep)[0]
            ix = ix[np.argsort(ex["key"][ix])[:n_ex]] if ix.size else ix
            for p0 in range(0, ix.size, 10):
                fig, axs = plt.subplots(5, 2, figsize=(14, 18), sharex=True)
                for ax, i in zip(axs.flat, list(ix[p0:p0 + 10]) + [None] * 10):
                    if i is None:
                        ax.axis("off")
                        continue
                    w = ex["samples"][i].astype(float) - ex["ped"][i]
                    pk = int(ex["peak"][i])
                    mark_pulse(ax, w, pk, int(ex["dip"][i]), int(ex["gap"][i]))
                    ax.set_xlim(max(pk - 10, 0), w.size - 1)
                    ax.set_title(f"{geo['mod_name'][ex['midx'][i]]}  run {ex['run'][i]} "
                                 f"split {ex['split'][i]}  ev {ex['event_number'][i]}\n"
                                 f"{label}: A={ex['amp'][i]:.0f}, bump={ex['bump_adc'][i]:.1f} ADC "
                                 f"({100 * ex['bump'][i]:.1f}% A), gap={ex['gap'][i]}", fontsize=11)
                for ax in axs[-1]:
                    ax.set_xlabel("sample (4 ns)")
                fig.suptitle(f"{label} tails with A >= {AMP_MIN:.0f} ADC "
                             f"({p0 + 1}-{min(p0 + 10, ix.size)} of {ix.size}); "
                             "blue = peak, red = dip, orange = secondary max")
                fig.tight_layout(rect=(0, 0, 1, 0.97))
                pdf.savefig(fig)
                plt.close(fig)


def plot_map(geo, n_good, rate, out):
    import matplotlib.pyplot as plt
    # plot style: large bold axis labels, large legends
    plt.rcParams.update({
        "axes.labelsize": 16, "axes.labelweight": "bold", "axes.titlesize": 16, "axes.titleweight": "bold",
        "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14, "axes.linewidth": 1.8,
        "xtick.major.width": 1.5, "ytick.major.width": 1.5, "xtick.major.size": 6, "ytick.major.size": 6,
        "lines.linewidth": 2.0, "figure.titlesize": 18, "figure.titleweight": "bold"})
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle
    pw = geo["mod_pwo4"].astype(bool)
    x, y, sx, sy = geo["mod_x"], geo["mod_y"], geo["mod_sx"], geo["mod_sy"]
    ids = np.nonzero(pw)[0]
    rect = {i: Rectangle((x[i] - sx[i] / 2, y[i] - sy[i] / 2), sx[i], sy[i]) for i in ids}
    has = n_good >= MIN_PULSES_MAP
    fig, ax = plt.subplots(figsize=(11, 9.5))
    ax.add_collection(PatchCollection([rect[i] for i in ids if not has[i]],
                                      facecolor="0.8", edgecolor="0.6", lw=0.2))
    sel = [i for i in ids if has[i]]
    if sel:
        pc = PatchCollection([rect[i] for i in sel], cmap="viridis", edgecolor="0.3", lw=0.2)
        pc.set_array(rate[sel])
        ax.add_collection(pc)
        cb = fig.colorbar(pc, ax=ax)
        cb.set_label("bumpy-tail rate (dip + secondary rise)", fontsize=16, fontweight="bold")
        cb.ax.tick_params(labelsize=14)
    ax.set_xlim((x[pw] - sx[pw] / 2).min(), (x[pw] + sx[pw] / 2).max())
    ax.set_ylim((y[pw] - sy[pw] / 2).min(), (y[pw] + sy[pw] / 2).max())
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title(f"Bumpy-tail rate per PbWO4 module, A >= {AMP_MIN:.0f} ADC "
                 f"\n(grey: < {MIN_PULSES_MAP} pulses; {len(sel)} modules)")
    for ext in ("png", "pdf"):
        fig.savefig(out / f"bumpy_rate_map.{ext}", dpi=150)
    plt.close(fig)


def plot_hist(n_good, rate, out):
    import matplotlib.pyplot as plt
    # plot style: large bold axis labels, large legends
    plt.rcParams.update({
        "axes.labelsize": 16, "axes.labelweight": "bold", "axes.titlesize": 16, "axes.titleweight": "bold",
        "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14, "axes.linewidth": 1.8,
        "xtick.major.width": 1.5, "ytick.major.width": 1.5, "xtick.major.size": 6, "ytick.major.size": 6,
        "lines.linewidth": 2.0, "figure.titlesize": 18, "figure.titleweight": "bold"})
    m = n_good >= MIN_N_HIST
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.hist(rate[m], bins=50, range=(0, max(rate[m].max(), 1e-3)) if m.any() else (0, 1),
            color="tab:blue")
    ax.set_xlabel("bumpy-tail rate per module (dip + secondary rise)")
    ax.set_ylabel("modules")
    ax.set_title(f"{int(m.sum())} modules with >= {MIN_N_HIST} good pulses "
                 f"(A >= {AMP_MIN:.0f} ADC)")
    fig.tight_layout()
    fig.savefig(out / "bumpy_rate_hist.png", dpi=150)
    plt.close(fig)


def report(scan_dir, out):
    import matplotlib
    matplotlib.use("Agg")
    files, geo, rec, ex, tot = load_scans(scan_dir)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    nm = geo["mod_name"].size
    bumpy = rec["cls"] == BUMPY
    amp_ok = rec["amp"] >= AMP_MIN          # all products but the all-A modules.tsv columns
    n_good = np.bincount(rec["midx"], minlength=nm)
    n_bumpy = np.bincount(rec["midx"][bumpy], minlength=nm)
    n_good_a = np.bincount(rec["midx"][amp_ok], minlength=nm)
    n_bumpy_a = np.bincount(rec["midx"][bumpy & amp_ok], minlength=nm)
    ns_per_sample = 4.0                     # waveform sample spacing (ns)
    gap_sum_a = np.bincount(rec["midx"][bumpy & amp_ok], weights=rec["gap"][bumpy & amp_ok],
                            minlength=nm)
    with np.errstate(invalid="ignore", divide="ignore"):
        rate = np.where(n_good > 0, n_bumpy / n_good, np.nan)
        rate_a = np.where(n_good_a > 0, n_bumpy_a / n_good_a, np.nan)
        gap_mean_a = np.where(n_bumpy_a > 0, gap_sum_a / n_bumpy_a * ns_per_sample, np.nan)

    with open(out / "modules.tsv", "w") as fh:
        fh.write("module\tx\ty\tn_good\tn_not_bumpy\tn_bumpy\trate_bumpy\t"
                 "n_good_a200\tn_bumpy_a200\trate_bumpy_a200\tgap_mean_ns_a200\n")
        for i in range(nm):
            fh.write(f"{geo['mod_name'][i]}\t{geo['mod_x'][i]:.2f}\t{geo['mod_y'][i]:.2f}\t"
                     f"{n_good[i]}\t{n_good[i] - n_bumpy[i]}\t{n_bumpy[i]}\t{rate[i]:.4f}\t"
                     f"{n_good_a[i]}\t{n_bumpy_a[i]}\t{rate_a[i]:.4f}\t{gap_mean_a[i]:.2f}\n")
    plot_map(geo, n_good_a, rate_a, out)
    plot_hist(n_good_a, rate_a, out)
    plot_examples(ex, geo, N_EXAMPLES, out / "examples.pdf")
    hx = ex.pop("hx")
    pwo = geo["mod_pwo4"].astype(bool)
    hot = [i for i in np.argsort(-np.nan_to_num(rate_a)) if pwo[i] and n_good_a[i] >= HOT_MIN_N
           and rate_a[i] > HOT_RATE]
    with open(out / "hot_modules.tsv", "w") as fh:
        fh.write("module\tx\ty\tn_good_a200\trate_bumpy_a200\n")
        for i in hot:
            fh.write(f"{geo['mod_name'][i]}\t{geo['mod_x'][i]:.2f}\t{geo['mod_y'][i]:.2f}\t"
                     f"{n_good_a[i]}\t{rate_a[i]:.4f}\n")
    print(f"hot modules (rate > {HOT_RATE:.0%}, n >= {HOT_MIN_N}): {len(hot)}")
    for i in hot:
        print(f"  {geo['mod_name'][i]:>6} x={geo['mod_x'][i]:8.2f} y={geo['mod_y'][i]:8.2f} "
              f"n={n_good_a[i]:6d} rate={rate_a[i]:.3f}")
    if hx:
        plot_hot_examples(hx, hot, geo, out / "hot_examples.pdf")

    n_all = int(bumpy.size)
    n, nb = int(amp_ok.sum()), int((bumpy & amp_ok).sum())
    plateau = (rec["cls"] == PLATEAU) & amp_ok
    w = np.array(tot["wall_s"])
    lines = [
        f"scan files: {len(files)} ({scan_dir}); runs {tot['runs']}; splits {tot['n_splits']}",
        f"events: {tot['n_events_seen']} seen, {tot['n_events_trig']} after the trigger "
        f"gate, {tot['n_events_excl']} vetoed by the LMS bit",
        f"good pulses: {n_all} (PbWO4); modules with good pulses: {int((n_good > 0).sum())}",
        f"good pulses with A >= {AMP_MIN:.0f} ADC: {n} ({n / max(n_all, 1):.2%} of all); "
        f"modules with such pulses: {int((n_good_a > 0).sum())}",
        f"smooth  (A >= {AMP_MIN:.0f} ADC): {n - nb - int(plateau.sum())} "
        f"({(n - nb - int(plateau.sum())) / max(n, 1):.2%})",
        f"plateau (A >= {AMP_MIN:.0f} ADC, no bump, descent breaks in < {TAIL_NSAMP} samples; "
        f"not in the rate): {int(plateau.sum())} ({int(plateau.sum()) / max(n, 1):.2%})",
        f"bumpy   (A >= {AMP_MIN:.0f} ADC, dip + rise >= max({BUMP_MIN_ADC:g} ADC, {BUMP_MIN_FRAC:.0%} A)): {nb} "
        f"({nb / max(n, 1):.2%}); bumpy rate = bumpy / all good",
        f"BUMP_MIN_ADC = {BUMP_MIN_ADC:g}, BUMP_MIN_FRAC = {BUMP_MIN_FRAC:g}, TAIL_END_ADC = {TAIL_END_ADC:g}, TAIL_NSAMP = {TAIL_NSAMP}, hot modules "
        f"(rate > {HOT_RATE:.0%}, n >= {HOT_MIN_N}): {len(hot)}",
        f"bumpy-rate per module (>= {MIN_N_HIST} good, A >= {AMP_MIN:.0f} ADC): median "
        f"{np.nanmedian(rate_a[n_good_a >= MIN_N_HIST]) if (n_good_a >= MIN_N_HIST).any() else float('nan'):.4f}",
        f"scan wall time per run: median {np.median(w):.0f} s, max {w.max():.0f} s",
    ]
    bs = bumpy & amp_ok
    lines += ["", f"dip position (dip - peak, samples) of bumpy pulses with A >= {AMP_MIN:.0f} ADC:"]
    for f, c in zip(*np.unique((rec["dip"] - rec["peak"])[bs], return_counts=True)):
        lines.append(f"  {int(f)}: {int(c)}")
    lines += ["", "gap (peak to secondary maximum, samples), same pulses:"]
    for f, c in zip(*np.unique(rec["gap"][bs], return_counts=True)):
        lines.append(f"  {int(f)}: {int(c)}")
    gap_ns = rec["gap"][bs] * ns_per_sample
    g_mean, g_std, g_med = ((gap_ns.mean(), gap_ns.std(), np.median(gap_ns)) if gap_ns.size
                            else (np.nan, np.nan, np.nan))
    gap_line = (f"gap main peak -> secondary peak (bumpy, A >= {AMP_MIN:.0f} ADC): "
                f"mean {g_mean:.2f} ns, sigma {g_std:.2f} ns (std), median {g_med:.2f} ns, "
                f"N = {gap_ns.size}")
    lines += ["", gap_line]
    (out / "summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:10]))
    print(gap_line)
    print(f"[report] wrote {out}/: examples.pdf hot_examples.pdf bumpy_rate_map.png/.pdf "
          "bumpy_rate_hist.png modules.tsv hot_modules.tsv summary.txt")


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) >= 3 and argv[0] == "scan":
        scan(argv[1], argv[2:])
    elif len(argv) == 3 and argv[0] == "report":
        report(argv[1], argv[2])
    else:
        sys.exit(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main())
