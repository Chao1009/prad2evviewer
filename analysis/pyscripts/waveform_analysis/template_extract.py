#!/usr/bin/env python3
"""template_extract.py -- per-module PbWO4 pulse templates (cubic B-splines) from clean single pulses.

    python3 template_extract.py select <out.npz> <run type> <evio file> [<evio file> ...]
    python3 template_extract.py build <run type>

select  one run (all its EVIO split files) -> one npz of every pulse that passes the data-quality
        cuts on its raw samples (Selector.flush, in order), each with its raw samples, module,
        amplitude, constant-fraction onset t0_cfd, rough energy e_rough, run and split.
build   pool OUT_ROOT/<RT>/pulses/run_*.npz.  Per module the fit pool is, in this order, the first whose
        fit succeeds: the pulses with e_rough >= POOL_EFRAC * E_beam (if at least POOL_MIN_EP), the pulses with
        A in [A_LO, A_HI] (if at least POOL_MIN_GEN; `pool_general`), the elastic pulses (at least 3), the single
        pulse of largest amplitude; else the run-type average.  Starting from
        t0_cfd, the pool is aligned iteratively to its own penalised cubic B-spline (lambda chosen on
        held-out pulses).  shape_sigma(t), the pulse-to-pulse spread about the template (what the matched
        filter reads), is measured for every module on its general band pool (A_LO <= A <= A_HI, up to POOL_MAX
        random pulses) aligned to the module's final template, whichever pool the template came from; fewer
        than SIGMA_MIN_POOL band pulses: the run-type median (flag sigma_borrowed).  a_band is the median amplitude of that
        band pool (ADC; run-type median where borrowed): shape_sigma is measured at a_band, the matched filter scales it
        by sqrt(a_band / A).  shape_sigma_elastic is the
        spread of the elastic fit pool itself (NaN unless that pool has >= SIGMA_MIN_POOL pulses).  Per module
        the same band pulses give the per-pulse 10-90 % rise and 90-10 % fall times of the raw samples: rise_jitter_ns /
        fall_jitter_ns (1.4826 MAD) and rise_med_ns / fall_med_ns (NaN below SIGMA_MIN_POOL pulses).
        Writes OUT_ROOT/<RT>/templates_<RT>.npz, .pdf, .root and n_templates_map.png.

Standalone: numpy, scipy, matplotlib, uproot / awkward and the compiled prad2py decoder binding.
Units: ADC, ns, samples of 4 ns.  All parameters are the constants below.
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import textwrap
import time
import warnings
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from prad2py import dec, det          # compiled PRad-II decoder / detector binding

HERE = Path(__file__).resolve().parent

# =============================================================================
# parameters
# =============================================================================
# database directory: $PRAD2_DATABASE_DIR (exported by setup_env.sh), else <repo>/database
DB_DIR = Path(os.environ.get("PRAD2_DATABASE_DIR") or HERE.parents[2] / "database")
OUT_ROOT = HERE / "pulse_templates"                                          # this script's output
HYCAL_MAP = DB_DIR / "hycal_map.json"                                        # module geometry (mm)

# run type: select trigger gate (physics: trigger_bits == 0x100; physics-any: trigger_bits & 0xf00)
# and calibration file of e_rough (None: the file database/runinfo/general.json assigns to the run)
RUN_TYPES = {
    "0p73GeV_H2":     dict(trigger="physics",     calib_file=None),
    "2p24GeV_H2":     dict(trigger="physics",     calib_file=None),
    "2p24GeV_C":      dict(trigger="physics",
                           calib_file=DB_DIR / "calibration/old/calibration_factor_2p2_fadc.json"),
    "2p24GeV_Ta_X17": dict(trigger="physics-any", calib_file=None),
    "3p49GeV_H2":     dict(trigger="physics",     calib_file=None),
}

# --- select
SCHEMA = "spline_cand_v4"       # select-file schema
NS_EXPECT = 100                 # samples per PbWO4 readout
DT_NS = 4.0                     # sample period (ns)
BLOCK = 256                     # events per batched pass
MAX_EVENTS = 0                  # stop after this many decoded events per file (0 = all)
PHYSICS_BITS_EXACT = 0x100
PHYSICS_MASK_ANY = 0xf00
EXCLUDE_TRIGGER_BITS = 0x1000000    # reject events with (trigger_bits & mask) != 0 (LMS bit)
MAX_MODULES_PER_EVENT = 20      # reject events with more channels at A >= MODULE_COUNT_HEIGHT_MIN
MODULE_COUNT_HEIGHT_MIN = 50.0
CLIP_ADC = 4095                 # c1: no raw sample >= this
A_STORE_MIN = 50.0              # c5: A = raw peak - WaveAnalyzer pedestal >= this
PEAK_MAX_SAMPLE = 59            # c4: peak sample <= this (>= 40 samples of tail)
THR_ADC = 2.0                   # c3: first_thr = first sample before the peak above this ...
PRE_SAMPLES = 20                # ... at sample >= this, so a PRE_SAMPLES pre-pulse window exists
PED_RMS_MAX = 1.5               # c8: WaveAnalyzer rms
PED_SLOPE_MAX = 0.05            # c8: |slope| of the pre-pulse window (ADC / sample)
TAIL_OFF0, TAIL_OFF1 = 40, 50   # c9: samples [peak+40, min(peak+50, 99)] ...
TAIL_ABS_ADC, TAIL_RMS_MULT, TAIL_FRAC_A = 2.0, 2.0, 0.04   # ... |mean| <= max(2, 2 rms) + 0.04 A
SEC_FRAC = 0.10                 # c2: no sample after the peak >= this * A above the lowest sample since the peak
TAIL_NSAMP = 8                  # c10: no decrease on the rising edge (the samples above THR_ADC up to the peak); of peak+1 .. peak+8 ...
NOISE_ALLOW = 2.0               # ... none above its predecessor by more than this * rms
CFD_FRAC = 0.5                  # t0_cfd: rising-edge crossing of this * A, linear between samples
E_PRE = 3                       # e_rough sums the samples from this many before floor(t0_cfd / 4 ns) to 99
EMAX_FRAC = 1.2                 # stored: e_rough < EMAX_FRAC * E_beam
CUTFLOW_NAMES = ["nseen", "c1_clip_ok", "c5_amp_ge50", "c4_contained", "c3_first_thr",
                 "c8_pedestal", "c9_tail", "c2_no_secondary", "c10_smooth", "stored"]

# --- build
POOL_EFRAC = 0.5                # fit pool: e_rough >= this * E_beam (the seed crystal holds ~70 % of an elastic cluster)
HOLDOUT_FRAC = 0.2              # held out (lambda choice): the next pulses below the pool in e_rough, this share of pool + held out
LAMBDA_MIN_POOL = 10            # lambda chosen on held-out pulses from this pool size on; smaller pools: LAMBDA_FALLBACK
SIGMA_MIN_POOL = 50             # own shape_sigma, jitter and shape_sigma_elastic from this pool size on; smaller pools: the run-type median / NaN
ALIGN_PASSES = 3                # spline fit -> onset refit of every pulse, this many times, then the final fit
ALIGN_HW = 4.0                  # onset refit within +-this (ns) of the CFD onset
POOL_MIN_EP = 20                # elastic pools of at least this many pulses are tried first; smaller ones only after the amplitude band
A_LO, A_HI = 200.0, 1000.0      # amplitude band (ADC) of that general pool: above the noise floor, below the elastic amplitudes
POOL_MIN_GEN = 20               # fewer band pulses than this: the module keeps its elastic pool, or with no elastic pulse its single largest pulse
POOL_MAX = 1000                 # random subsample cap: sigma(t) and the spline saturate well below this; keeps build minutes
SEED = 1                        # random seed (plotted points, run-type average draw)
SIG_MIN = 0.3                   # ADC floor on the pedestal rms
WIN = (-20.0, 160.0)            # alignment / chi2 window (ns after onset)
POOL_TAU_MIN = -20.0            # pooled points: every sample with tau >= this
BIN_W = 0.5                     # bin width of the pooled points (ns)
SHAPE_SIGMA = 0.012             # (SHAPE_SIGMA * ref(tau))^2 added to the point variance of the fit
SIGMA_T_NS = 0.15               # timing-jitter term of the point variance (ns)
LAM_GRID = (1e2, 1e3, 1e4)      # smoothing lambda: the value with the lowest held-out chi2/dof
LAMBDA_FALLBACK = 1e3           # lambda of pools below LAMBDA_MIN_POOL
SMOOTH_LAMBDA = 1e4             # lambda of the run-type average and the wiggle cap
OSC_PROM = 0.002                # spurious extremum: prominence >= this * peak ...
OSC_WIN = (2.0, 150.0)          # ... in this window (ns) ...
OSC_EXCL = 3.0                  # ... farther than this from the main peak (ns) ...
OSC_PEAK_WIN = (3.0, 30.0)      # ... the largest value in here
OSC_POOL_MAX = 100              # wiggle check only below this pool size: larger pools fix the shape at every sample phase
ONSET_FRAC = 0.005              # onset = last time before the peak at <= 0.5 % of it
T_END = 400.0                   # stored template extent (ns after onset)
T_AXIS_END = 416.0              # spline domain end on the pre-shift axis
T_LEFT = -20.0                  # spline domain start
_NCOL = 52                      # alignment: sample columns kept around the seed
TAIL_BIN_MIN_PTS = 5            # T_data: last bin with >= this many points
TAIL_FIT_FROM = 140.0           # exponential tail fitted on [140, T_data]
TAIL_PSEUDO_DT = 4.0            # spacing of the tail pseudo-points (ns)
TAIL_TAU_DEFAULT = 170.0        # tail time constant if the fit fails (ns)
GRID_DT = 0.1                   # t_grid step of the stored shapes (ns)
T_GRID = np.round(-8.0 + GRID_DT * np.arange(int(round((T_END + 8.0) / GRID_DT)) + 1), 6)
LUT_DT = 0.05                   # step of the sampled-shape lookup table (ns)
NIPH = 256                      # onset phases per 4 ns sample of the energy integral
TS = DT_NS * np.arange(NS_EXPECT)
# shape_sigma(t): spread of y/A - spline over the pool, pedestal noise subtracted in quadrature
SIGMA_BIN_NS = 2.0              # tau bin width (ns)
SIGMA_BIN_MIN = 8               # points needed in a bin
PRE_ONSET_NS = -20.0            # samples at tau < this give the pedestal noise
CHI2_1_MED = 0.4549364          # median of a chi2 with 1 dof
MAP_LIM = 350.0                 # heat map: x and y range [-MAP_LIM, MAP_LIM] mm
EX_RING2, EX_BORROWED, EX_PER_MODULE, EX_COLS = 5, 5, 8, 4   # PDF examples: one page per module (every innermost-ring module, EX_RING2 of the second ring, EX_BORROWED outer modules with a borrowed shape_sigma), EX_PER_MODULE random pool pulses each, EX_COLS panels per row


# shared helpers ==============================================================
def _git_short():
    try:
        return subprocess.check_output(["git", "-C", str(HERE), "rev-parse", "--short", "HEAD"],
                                       stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:                                   # noqa: BLE001
        return "unknown"


def _sha1(path):
    h = hashlib.sha1()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def _save_npz(path, arrays):
    """np.savez_compressed through a temp file + rename."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "wb") as fh:
        np.savez_compressed(fh, **arrays)
    os.replace(tmp, path)


def _constants(names):
    g = globals()
    return {k: (list(g[k]) if isinstance(g[k], tuple) else g[k]) for k in names}


# detector and calibration ===================================================
def dead_modules(run):
    """(from_run, set of names) of the hycal_dead_modules entry of database/runinfo/general.json
    with the largest from_run <= run (independent lookup, no chaining)."""
    ents = json.load(open(DB_DIR / "runinfo" / "general.json"))["hycal_dead_modules"]
    e = max((e for e in ents if e["from_run"] <= run), key=lambda e: e["from_run"])
    return e["from_run"], set(e["dead_modules"])


def setup_detector(run):
    """prad2py setup of `run`: cfg (DAQ config), crate_map (ROC tag -> crate), hycal, wave_ana
    (WaveAnalyzer), geo (geo.Ebeam in MeV), core (keeps them alive), dead (module indices of
    dead_modules(run)) and dead_from (its from_run)."""
    b = det.PipelineBuilder()
    b.set_run_number(int(run))
    core = b.build()
    dead_from, dead = dead_modules(int(run))
    return SimpleNamespace(core=core, cfg=core.daq_cfg, geo=core.run_cfg, hycal=core.hycal,
                           crate_map={int(r.tag): int(r.crate) for r in core.daq_cfg.roc_tags},
                           wave_ana=dec.WaveAnalyzer(core.daq_cfg.wave_cfg), dead_from=dead_from,
                           dead={i for i in range(core.hycal.module_count())
                                 if core.hycal.module(i).name in dead})


def wave_config(wcfg):
    """WaveAnalyzer settings that define the energy integral."""
    return dict(smooth_order=int(wcfg.smooth_order), int_tail_ratio=float(wcfg.int_tail_ratio),
                tail_break_n=int(wcfg.tail_break_n), ped_nsamples=int(wcfg.ped_nsamples))


def resolve_calib(run):
    """Calibration file that database/runinfo/general.json assigns to `run`: the entry with the
    largest from_run <= run, inheriting the file of the earlier entries."""
    cfgs = json.load(open(DB_DIR / "runinfo" / "general.json"))["configurations"]
    cfgs = sorted(((c.get("from_run", c.get("run_number")), c) for c in cfgs
                   if c.get("from_run", c.get("run_number")) is not None), key=lambda t: t[0])
    f = ""
    for r0, c in cfgs:
        if r0 > run:
            break
        f = c.get("calibration", {}).get("file", f)
    if not f:
        sys.exit(f"[ERROR] runinfo has no calibration file for run {run}")
    path = Path(f) if os.path.isabs(f) else DB_DIR / f
    if not path.exists():
        sys.exit(f"[ERROR] calibration file for run {run} not found: {path}")
    return str(path)


def read_calib(path, mods, min_n=100):
    """Per-module calibration factor (NaN where <= 0 or missing)."""
    by_name = {m["name"]: float(m.get("factor", 0.0)) for m in json.load(open(path))}
    cal = np.array([by_name.get(m.name, 0.0) for m in mods])
    if (cal > 0).sum() < min_n:
        sys.exit(f"[ERROR] calibration {path}: only {(cal > 0).sum()} modules with factor > 0")
    return np.where(cal > 0, cal, np.nan)


def smooth_tri(x, order):
    """WaveAnalyzer triangular smoothing along the last axis (edge samples unsmoothed)."""
    x = np.asarray(x, dtype=np.float64)
    if order <= 1:
        return x.copy()
    n = x.shape[-1]
    val = x.copy()
    wsum = np.ones(n)
    for j in range(1, order):
        if 2 * j >= n:
            break
        w = 1.0 - j / float(order + 1)
        val[..., j:n - j] += w * (x[..., :n - 2 * j] + x[..., 2 * j:])
        wsum[j:n - j] += 2.0 * w
    return val / wsum


def walk_integral(v, ratio, nbreak):
    """WaveAnalyzer peak integral of the smoothed rows v[R, N]: walk outward from the maximum; a
    sample is 'below' if v < ratio * peak or v < 0; nbreak consecutive below samples end the walk
    (those pending samples are dropped)."""
    R, N = v.shape
    ar = np.arange(R)
    pk = np.argmax(v, axis=1)
    h = v[ar, pk]
    cut = ratio * h
    integral = h.copy()
    nb = max(1, int(nbreak))
    for side in (-1, 1):
        run = np.zeros(R, dtype=np.int64)
        pend = np.zeros(R)
        active = np.ones(R, dtype=bool)
        for d in range(1, N):
            pos = pk + side * d
            inr = active & (pos >= 0) & (pos < N)
            if not inr.any():
                break
            vv = v[ar, np.clip(pos, 0, N - 1)]
            below = (vv < cut) | (vv < 0.0)
            mb = inr & below
            run = np.where(mb, run + 1, run)
            pend = np.where(mb, pend + vv, pend)
            active &= ~(mb & (run >= nb))
            ma = inr & ~below
            integral += np.where(ma, pend + vv, 0.0)
            pend = np.where(ma, 0.0, pend)
            run = np.where(ma, 0, run)
    return integral


def wa_integral(shape, wave):
    """WaveAnalyzer integral of the unit shape (LutShape) sampled every 4 ns, mean over NIPH onset
    phases (ADC x samples per unit peak)."""
    tau = DT_NS * (np.arange(64) - 12)[None, :] + (DT_NS / NIPH) * np.arange(NIPH)[:, None]
    v = smooth_tri(shape.eval(tau), wave["smooth_order"])
    return float(walk_integral(v, wave["int_tail_ratio"], wave["tail_break_n"]).mean())


# select ======================================================================
_SPLIT_PAT = re.compile(r"prad_0*(\d+)\.evio\.(\d+)", re.IGNORECASE)
# per-pulse arrays of the select npz and their dtypes
_CAND = dict(event_number=np.int64, midx=np.int64, samples=np.uint16, ped=np.float32,
             rms=np.float32, amp=np.float32, peak_sample=np.int16, t0_cfd=np.float32,
             e_rough=np.float32)


def pulse_features(W):
    """Pedestal-subtracted W[n, ns] (A >= A_STORE_MIN) -> peak (first maximum), amp, first_thr
    (first sample <= peak above THR_ADC), falling (strictly falling samples after the peak)."""
    n, ns = W.shape
    peak = np.argmax(W, axis=1)
    above = (W > THR_ADC) & (np.arange(ns)[None, :] <= peak[:, None])
    desc = (W[:, 1:] < W[:, :-1]) | (np.arange(ns - 1)[None, :] < peak[:, None])
    falling = np.where(desc.all(axis=1), ns - 1, np.argmin(desc, axis=1)) - peak
    return dict(peak=peak, amp=W[np.arange(n), peak], first_thr=np.argmax(above, axis=1), falling=falling)


def cfd_onset(W, amp):
    """t0_cfd (ns): where the rising edge first reaches CFD_FRAC * amp, linear between the samples."""
    lev = CFD_FRAC * amp
    i1 = np.argmax(W >= lev[:, None], axis=1)
    r = np.arange(W.shape[0])
    w0, w1 = W[r, i1 - 1], W[r, i1]
    return DT_NS * (i1 - 1 + (lev - w0) / (w1 - w0))


class Selector:
    """Block-batched selection: the rows of a block (raw PbWO4 waveforms, BLOCK events) are cut
    with numpy on their raw samples.  cal: calibration factor per module index (MeV per ADC
    sample, NaN: none)."""
    smooth_cut = True               # c10 applied (pulse_shape_study counts it but keeps the pulses)

    def __init__(self, p, cal, emax_mev):
        self.p, self.cal, self.emax_mev = p, cal, emax_mev
        self.nm = p.hycal.module_count()
        self.cutflow = np.zeros((self.nm, len(CUTFLOW_NAMES)), dtype=np.int64)
        self.rows, self.meta = [], []
        self.out = {k: [] for k in _CAND}
        self.n_kept = self.n_rows_mult_rej = self.n_above_emax = 0

    def add(self, row, ev, evnum, trig, midx):
        self.rows.append(row)
        self.meta.append((ev, evnum, trig, midx))

    def kept(self):
        return self.n_kept

    def _count(self, name, midx):
        np.add.at(self.cutflow[:, CUTFLOW_NAMES.index(name)], midx, 1)

    def flush(self):
        """Cut the stacked block and keep the survivors."""
        if not self.rows:
            return
        Y = np.stack(self.rows)
        meta = np.array(self.meta, dtype=np.int64)               # ev, evnum, trig, midx
        self.rows, self.meta = [], []
        ana = self.p.wave_ana.analyze
        pr = np.array([ana(y)[:2] for y in Y]).reshape(-1, 2)    # WaveAnalyzer pedestal, rms
        ped, rms = pr[:, 0], pr[:, 1]
        # busy events (LMS pulses light up the whole calorimeter)
        _, inv = np.unique(meta[:, 0], return_inverse=True)
        cnt = np.bincount(inv, weights=(Y.max(axis=1) - ped) >= MODULE_COUNT_HEIGHT_MIN)
        k = cnt[inv] <= MAX_MODULES_PER_EVENT
        self.n_rows_mult_rej += int((~k).sum())
        Y, meta, ped, rms = Y[k], meta[k], ped[k], rms[k]
        mid = meta[:, 3]
        self._count("nseen", mid)
        k = Y.max(axis=1) < CLIP_ADC
        self._count("c1_clip_ok", mid[k])
        k &= (Y.max(axis=1) - ped) >= A_STORE_MIN
        self._count("c5_amp_ge50", mid[k])
        k &= np.argmax(Y, axis=1) <= PEAK_MAX_SAMPLE
        self._count("c4_contained", mid[k])
        s = np.nonzero(k)[0]
        if s.size == 0:
            return
        Y, meta, ped, rms, mid = Y[s], meta[s], ped[s], rms[s], mid[s]
        W = Y.astype(np.float64) - ped[:, None]
        P = pulse_features(W)
        pk = P["peak"][:, None]
        k = P["first_thr"] >= PRE_SAMPLES
        self._count("c3_first_thr", mid[k])
        g = np.nonzero(k)[0]
        cols = P["first_thr"][g, None] - PRE_SAMPLES + np.arange(PRE_SAMPLES)
        xs = np.arange(PRE_SAMPLES, dtype=np.float64) - (PRE_SAMPLES - 1) / 2.0
        slope = (W[g[:, None], cols] * xs).sum(axis=1) / (xs * xs).sum()
        k[g] &= (rms[g] <= PED_RMS_MAX) & (np.abs(slope) <= PED_SLOPE_MAX)
        self._count("c8_pedestal", mid[k])
        idx = np.arange(NS_EXPECT)[None, :]
        tm = (idx >= pk + TAIL_OFF0) & (idx <= np.minimum(pk + TAIL_OFF1, NS_EXPECT - 1))
        tail = (W * tm).sum(axis=1) / tm.sum(axis=1)
        k &= np.abs(tail) <= np.maximum(TAIL_ABS_ADC, TAIL_RMS_MULT * rms) + TAIL_FRAC_A * P["amp"]
        self._count("c9_tail", mid[k])
        after = np.where(idx >= pk, W, np.inf)
        bump = np.where(idx > pk, W - np.minimum.accumulate(after, axis=1), 0.0).max(axis=1)
        k &= bump < SEC_FRAC * P["amp"]
        self._count("c2_no_secondary", mid[k])
        d, j = np.diff(W, axis=1), idx[:, :-1]                   # d[:, j] = W[:, j+1] - W[:, j]
        edge0 = NS_EXPECT - np.argmax(((W <= THR_ADC) & (idx < pk))[:, ::-1], axis=1)   # rising edge start
        rise_ok = ((d >= 0) | (j < edge0[:, None]) | (j >= pk)).all(axis=1)
        tail_ok = ((d <= NOISE_ALLOW * rms[:, None]) | (j < pk) | (j >= pk + TAIL_NSAMP)).all(axis=1)
        smooth = k & rise_ok & tail_ok
        self._count("c10_smooth", mid[smooth])
        if self.smooth_cut:
            k = smooth
        s = np.nonzero(k)[0]
        if s.size == 0:
            return
        t0 = cfd_onset(W[s], P["amp"][s])
        k0 = np.floor(t0 / DT_NS).astype(np.int64) - E_PRE
        E = self.cal[mid[s]] * (W[s] * (idx >= k0[:, None])).sum(axis=1)
        ok = E < self.emax_mev
        self.n_above_emax += int((~ok).sum())
        f = s[ok]
        self._keep(f, lambda ix: dict(
            event_number=meta[ix, 1], midx=mid[ix], samples=Y[ix], ped=ped[ix], rms=rms[ix],
            amp=P["amp"][ix], peak_sample=P["peak"][ix], falling=P["falling"][ix], t0_cfd=t0[ok],
            e_rough=E[ok]))

    def _keep(self, f, fields):
        """Append the stored fields of the survivors f (indices into the block's cut rows)."""
        d = fields(f)
        for key in self.out:
            self.out[key].append(d[key])
        self._count("stored", d["midx"])
        self.n_kept += f.size

    def arrays(self):
        return {k: (np.concatenate(self.out[k]).astype(d) if self.out[k] else
                    np.zeros((0, NS_EXPECT) if k == "samples" else (0,), d)) for k, d in _CAND.items()}


def select_setup(run, run_type):
    """Detector setup, calibration factor per module index (cal), E_max and trigger gate of `run`."""
    rt = RUN_TYPES[run_type]
    p = setup_detector(run)
    mods = [p.hycal.module(i) for i in range(p.hycal.module_count())]
    calib_path = str(rt["calib_file"]) if rt["calib_file"] else resolve_calib(run)
    ebeam = float(p.geo.Ebeam)
    if rt["trigger"] == "physics":
        gate = lambda fe: int(fe.info.trigger_bits) == PHYSICS_BITS_EXACT              # noqa: E731
        trig_rule = f"trigger_bits == 0x{PHYSICS_BITS_EXACT:x}"
    else:
        gate = lambda fe: (int(fe.info.trigger_bits) & PHYSICS_MASK_ANY) != 0          # noqa: E731
        trig_rule = f"(trigger_bits & 0x{PHYSICS_MASK_ANY:x}) != 0"
    dead = sorted(mods[i].name for i in p.dead)
    print(f"[setup] run {run}, trigger {rt['trigger']}, calib {calib_path}, E_beam {ebeam:.1f} MeV; "
          f"dead modules excluded ({len(dead)}, general.json hycal_dead_modules from_run {p.dead_from}): "
          f"{' '.join(dead)}", flush=True)
    return SimpleNamespace(
        p=p, mods=mods, wave_cfg=wave_config(p.wave_ana.cfg), gate=gate, trig_rule=trig_rule,
        emax_mev=EMAX_FRAC * ebeam, cal=read_calib(calib_path, mods),
        meta=dict(calib_file=calib_path, calib_sha1=_sha1(calib_path), e_beam_mev=ebeam,
                  e_max_mev=EMAX_FRAC * ebeam, dead_modules=dead,
                  dead_modules_source=f"runinfo/general.json hycal_dead_modules from_run {p.dead_from}"))


def _scan_events(evio, p, sel, gate, t_start):
    """Decode the file, gate the events and feed the PbWO4 channels (dead modules p.dead skipped)
    with NS_EXPECT samples to `sel` (flushed every BLOCK events).  Returns (events seen, after the gate, LMS-rejected)."""
    ch = dec.EvChannel()
    ch.set_config(p.cfg)
    if ch.open_auto(str(evio)) != dec.Status.success:
        sys.exit(f"[ERROR] cannot open {evio}")
    cache = {}                                  # (roc tag, slot, ch) -> midx | -1
    n_seen = n_trig = n_blk = n_excl = 0
    done = False
    while not done and ch.read() == dec.Status.success:
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            d = ch.decode_event(i, with_ssp=False)
            if not d["ok"]:
                continue
            ev_idx = n_seen
            n_seen += 1
            fe = d["event"]
            done = bool(MAX_EVENTS and n_seen >= MAX_EVENTS)
            trig = int(fe.info.trigger_bits)
            if gate(fe):
                n_trig += 1
                if trig & EXCLUDE_TRIGGER_BITS:
                    n_excl += 1
                else:
                    evnum = int(fe.info.event_number)
                    for ri in range(fe.nrocs):
                        roc = fe.roc(ri)
                        crate = p.crate_map.get(roc.tag) if roc.present else None
                        if crate is None:
                            continue
                        for s in roc.present_slots():
                            slot = roc.slot(s)
                            for c in slot.present_channels():
                                mi = cache.get((roc.tag, s, c))
                                if mi is None:
                                    m = p.hycal.module_by_daq(crate, s, c)
                                    mi = cache[(roc.tag, s, c)] = (
                                        m.index if m is not None and m.is_hycal() and m.is_pwo4()
                                        and m.index not in p.dead else -1)
                                cd = slot.channel(c)
                                if mi >= 0 and cd.nsamples == NS_EXPECT:
                                    sel.add(np.asarray(cd.samples, dtype=np.uint16), ev_idx, evnum, trig, mi)
                    n_blk += 1
                    if n_blk >= BLOCK:
                        sel.flush()
                        n_blk = 0
            if n_seen % 20000 == 0:
                print(f"  [progress] {n_seen} events, stored {sel.kept()}, "
                      f"{n_seen / (time.monotonic() - t_start):.0f} ev/s", flush=True)
            if done:
                break
    sel.flush()
    return n_seen, n_trig, n_excl


def select_candidates(out, run_type, evios):
    """`select`: the EVIO split files of one run -> one npz of `run_type` (per-pulse `split`)."""
    t_start = time.monotonic()
    found = []
    for e in evios:
        mt = _SPLIT_PAT.search(os.path.basename(e))
        if not mt:
            sys.exit(f"[ERROR] cannot read run/split from {e} (expected prad_0RRRRR.evio.SSSSS)")
        found.append((int(mt.group(1)), int(mt.group(2))))
    run = found[0][0]
    if any(r != run for r, _ in found):
        sys.exit(f"[ERROR] the EVIO files belong to different runs: {sorted(set(r for r, _ in found))}")
    if not str(out).endswith(".npz"):
        sys.exit("[ERROR] <out.npz> must end in .npz")
    S = select_setup(run, run_type)
    parts, splits, per_split = [], [], []
    cutflow = np.zeros((len(S.mods), len(CUTFLOW_NAMES)), dtype=np.int64)
    n_mult = n_emax = 0
    for evio, (_, split) in zip(evios, found):
        sel = Selector(S.p, S.cal, S.emax_mev)
        n_s, n_t, n_x = _scan_events(evio, S.p, sel, S.gate, time.monotonic())
        if not MAX_EVENTS and n_s == 0:
            sys.exit(f"[ERROR] no events decoded in {evio}")
        arr = sel.arrays()
        parts.append(arr)
        splits.append(np.full(arr["midx"].size, split, dtype=np.int16))
        cutflow += sel.cutflow
        n_mult += sel.n_rows_mult_rej
        n_emax += sel.n_above_emax
        per_split.append(dict(evio=os.path.abspath(evio), split=split, n_events_seen=n_s,
                              n_events_after_trigger=n_t, n_lms_rejected=n_x,
                              n_candidates=int(arr["midx"].size)))
    arr = {k: np.concatenate([a[k] for a in parts]) for k in parts[0]}
    names = np.array([m.name for m in S.mods])
    meta = dict(schema=SCHEMA, run=run, run_type=run_type, splits=per_split,
                cuts=_constants(["CLIP_ADC", "A_STORE_MIN", "PEAK_MAX_SAMPLE", "THR_ADC",
                                 "PRE_SAMPLES", "PED_RMS_MAX", "PED_SLOPE_MAX", "TAIL_OFF0",
                                 "TAIL_OFF1", "TAIL_ABS_ADC", "TAIL_RMS_MULT", "TAIL_FRAC_A",
                                 "SEC_FRAC", "TAIL_NSAMP", "NOISE_ALLOW", "CFD_FRAC", "E_PRE",
                                 "EMAX_FRAC", "MAX_MODULES_PER_EVENT", "MODULE_COUNT_HEIGHT_MIN",
                                 "EXCLUDE_TRIGGER_BITS", "MAX_EVENTS"]),
                trigger=RUN_TYPES[run_type]["trigger"], wave_config=S.wave_cfg, **S.meta,
                git=_git_short(), created=datetime.now().isoformat(timespec="seconds"))
    _save_npz(out, dict(run=np.int64(run), split=np.concatenate(splits), name=names[arr["midx"]],
                        **arr, cutflow_names=np.array(CUTFLOW_NAMES), cutflow=cutflow,
                        meta=np.array(json.dumps(meta))))
    pm = np.array([m.is_pwo4() for m in S.mods])
    nk = np.bincount(arr["midx"], minlength=len(S.mods))[pm]
    print(f"[done] {out}  run {run}, {len(evios)} files, {sum(s['n_events_seen'] for s in per_split)} "
          f"events in {time.monotonic() - t_start:.0f} s; {Path(out).stat().st_size / 1e6:.1f} MB")
    print("[cutflow] " + "  ".join(f"{n} {int(v)}" for n, v in zip(CUTFLOW_NAMES, cutflow[pm].sum(0))))
    print(f"[stored] {arr['midx'].size} pulses, per module median {np.median(nk):.0f}, "
          f"{int((nk == 0).sum())} of {int(pm.sum())} modules empty; busy-event rows dropped {n_mult}, "
          f"e_rough >= E_max or no calibration factor {n_emax}", flush=True)


# build =======================================================================
# Times tau are ns after the template onset unless noted.
class LutShape:
    """Shape sampled on t = 0, 0.05, ..., T_END; linear interpolation, 0 for tau <= 0 or >= T_END."""

    def __init__(self, v):
        self.v = np.asarray(v, dtype=np.float64)
        self.n = self.v.size

    @classmethod
    def from_spline(cls, t, c):
        return cls(spline_values(t, c, LUT_DT * np.arange(int(round(T_END / LUT_DT)) + 1)))

    def eval(self, tau):
        tau = np.asarray(tau, dtype=np.float64)
        x = tau * (1.0 / LUT_DT)
        i = np.clip(np.floor(x), 0, self.n - 2).astype(np.int64)
        f = x - i
        v = self.v[i] * (1.0 - f) + self.v[i + 1] * f
        return np.where((tau > 0.0) & (tau < T_END), v, 0.0)


class PulseSet:
    """Pedestal-subtracted pulses Y[n, 100] of the candidate indices idx, their noise sig, amp and
    CFD onset cfd (ns)."""

    def __init__(self, D, idx):
        self.Y = D["samples"][idx].astype(np.float64) - D["ped"][idx, None]
        self.sig = np.maximum(D["rms"][idx], SIG_MIN)
        self.amp = D["amp"][idx].astype(np.float64)
        self.cfd = D["t0_cfd"][idx].astype(np.float64)
        self.n = idx.size


def _align_columns(t_seed):
    """Sample columns [base, base + _NCOL) around each seed time (ns)."""
    base = np.clip(np.floor((t_seed - 26.0) / DT_NS).astype(np.int64), 0, NS_EXPECT - _NCOL)
    return base[:, None] + np.arange(_NCOL)[None, :]


def _chi2_at(Ys, cols, sig, t0, shape, win):
    """A and chi2/dof of every pulse at fixed t0 for the samples Ys of `cols`."""
    tau = TS[cols] - t0[:, None]
    msk = (tau >= win[0]) & (tau <= win[1])
    m = shape.eval(tau) * msk
    den = (m * m).sum(1)
    A = np.where(den > 1e-12, (Ys * m).sum(1) / np.maximum(den, 1e-12), 0.0)
    r = (Ys - A[:, None] * m) * msk
    nwin = msk.sum(1)
    chi2 = (r * r).sum(1) / (sig * sig)
    ok = (A > 0) & np.isfinite(chi2) & (nwin > 10)
    return dict(t0=t0, A=A, chi2dof=np.where(ok, chi2 / np.maximum(nwin - 2, 1), np.inf), ok=ok)


def align_pulses(Y, sig, t_seed, shape, halfwidth, win=WIN):
    """(A, t0) per pulse against `shape`: closed-form A = sum(y m) / sum(m^2) over the samples in
    `win`; t0 searched within +-halfwidth of the seed in nested 0.25 / 0.05 / 0.01 ns passes."""
    n = Y.shape[0]
    cols = _align_columns(np.asarray(t_seed, dtype=np.float64))
    Ys = np.take_along_axis(Y, cols, axis=1)

    def sweep(center, offs):
        bt, bs = center.copy(), np.full(n, -np.inf)
        for o in offs:
            t0 = center + o
            tau = TS[cols] - t0[:, None]
            m = shape.eval(tau) * ((tau >= win[0]) & (tau <= win[1]))
            num, den = (Ys * m).sum(1), (m * m).sum(1)
            s = np.where(den > 1e-12, num * np.abs(num) / np.maximum(den, 1e-12), -np.inf)
            up = s > bs
            bs, bt = np.where(up, s, bs), np.where(up, t0, bt)
        return bt

    best = sweep(np.asarray(t_seed, dtype=np.float64), np.arange(-halfwidth, halfwidth + 1e-9, 0.25))
    best = sweep(best, np.arange(-0.25, 0.25 + 1e-9, 0.05))
    best = sweep(best, np.arange(-0.05, 0.05 + 1e-9, 0.01))
    return _chi2_at(Ys, cols, sig, np.round(best, 6), shape, win)


def _chi2_full(Y, sig, t0, A, shape, win=WIN):
    """chi2/dof with the variance sigma^2 + (A SIGMA_T_NS f')^2 + (A SHAPE_SIGMA f)^2 over `win`
    (the held-out measure of the lambda selection)."""
    cols = _align_columns(t0)
    Ys, tau = np.take_along_axis(Y, cols, axis=1), TS[cols] - t0[:, None]
    msk = (tau >= win[0]) & (tau <= win[1])
    f = shape.eval(tau)
    slope = (shape.eval(tau + 0.05) - shape.eval(tau - 0.05)) / 0.1
    var = sig[:, None] ** 2 + (A[:, None] * SIGMA_T_NS * slope) ** 2 + (A[:, None] * SHAPE_SIGMA * f) ** 2
    return (((Ys - A[:, None] * f) ** 2 / var) * msk).sum(1) / np.maximum(msk.sum(1) - 2, 1)


def _grouped_median(vals, b, nb):
    """Per-bin median of vals (bins b) and the counts; NaN in empty bins."""
    cnt = np.bincount(b, minlength=nb)
    vs = vals[np.lexsort((vals, b))]
    start = np.cumsum(cnt) - cnt
    ok = cnt > 0
    med = 0.5 * (vs[np.where(ok, start + (cnt - 1) // 2, 0)] + vs[np.where(ok, start + cnt // 2, 0)])
    return np.where(ok, med, np.nan), cnt


def pool_points(Y, sig, t0, A, sel, ref):
    """Flat (tau, y/A, sigma/A, y/A detrended to the BIN_W bin centre by ref, bin) of the pulses in
    `sel` for every sample with tau >= POOL_TAU_MIN."""
    idx = np.flatnonzero(sel)
    tau = TS[None, :] - t0[idx, None]
    keep = tau >= POOL_TAU_MIN
    a = A[idx][:, None]
    y = (Y[idx] / a)[keep]
    sg = np.broadcast_to((sig[idx][:, None] / a), tau.shape)[keep]
    tau = tau[keep]
    b = np.floor((tau - POOL_TAU_MIN) / BIN_W).astype(np.int64)
    ctr = POOL_TAU_MIN + (b + 0.5) * BIN_W
    return dict(tau=tau, y=y, sg=sg, yd=y - (ref.eval(tau) - ref.eval(ctr)), bin=b)


def bin_stats(pp, ref):
    """BIN_W bins: centre, n, median, sigma = max(1.4826 MAD / sqrt(n), floor), floor =
    sqrt(noise^2 + (SIGMA_T_NS ref')^2) / sqrt(n) with noise the mean sigma_pulse / A of the bin."""
    b = pp["bin"]
    nb = int(b.max()) + 1
    med, cnt = _grouped_median(pp["yd"], b, nb)
    mad, _ = _grouped_median(np.abs(pp["yd"] - np.nan_to_num(med)[b]), b, nb)
    ok = cnt > 0
    ctr = POOL_TAU_MIN + (np.arange(nb) + 0.5) * BIN_W
    msg = np.bincount(b, weights=pp["sg"], minlength=nb) / np.maximum(cnt, 1)
    slope = (ref.eval(ctr + 0.05) - ref.eval(ctr - 0.05)) / 0.1
    rn = np.sqrt(np.maximum(cnt, 1))
    sigma = np.maximum(1.4826 * np.nan_to_num(mad) / rn, np.sqrt(msg ** 2 + (SIGMA_T_NS * slope) ** 2) / rn)
    return dict(x=ctr[ok], n=cnt[ok], y=med[ok], s=sigma[ok])


def knot_plan(n_pool, scale=1):
    """Interior knots (ns, pre-shift axis): -10, -1, then 0.5 (pool >= 100) or 1 ns to 10 ns, 1 to
    20, 2 to 40, 4 to 80, 8 to 160, 16 to 400; `scale` multiplies every spacing."""
    s0 = 0.5 if n_pool >= 100 else 1.0
    k = [-10.0, -1.0]
    for lo, hi, st in [(-1.0, 10.0, s0), (10.0, 20.0, 1.0), (20.0, 40.0, 2.0),
                       (40.0, 80.0, 4.0), (80.0, 160.0, 8.0), (160.0, 400.0, 16.0)]:
        k.extend(np.arange(lo + st * scale, hi + 1e-9, st * scale))
    k = np.unique(np.round(np.asarray(k), 6))
    return k[(k > T_LEFT) & (k < T_AXIS_END)]


def tail_extension(bs):
    """T_data (last bin with >= TAIL_BIN_MIN_PTS points), exponential a exp(-t/tau) fitted on
    [TAIL_FIT_FROM, T_data] and its pseudo-points every TAIL_PSEUDO_DT beyond T_data."""
    from scipy.optimize import curve_fit
    big = bs["n"] >= TAIL_BIN_MIN_PTS
    t_data = float(bs["x"][big].max()) if big.any() else float(bs["x"].max())
    sel = (bs["x"] >= TAIL_FIT_FROM) & (bs["x"] <= t_data)
    a_s = tau_s = np.nan
    if sel.sum() >= 5:
        x, y, s = bs["x"][sel], bs["y"][sel], bs["s"][sel]
        try:
            pos = y > 0
            p0 = (0.02, TAIL_TAU_DEFAULT)
            if pos.sum() >= 3:
                c1 = np.polyfit(x[pos] - TAIL_FIT_FROM, np.log(y[pos]), 1)
                if c1[0] < 0:
                    p0 = (math.exp(c1[1]), -1.0 / c1[0])
            popt, _ = curve_fit(lambda t, a, tau: a * np.exp(-(t - TAIL_FIT_FROM) / tau), x, y, p0=p0,
                                sigma=s, absolute_sigma=True, bounds=([0.0, 20.0], [1.0, 2000.0]), maxfev=2000)
            tau_s = float(popt[1])
            a_s = float(popt[0]) * math.exp(TAIL_FIT_FROM / tau_s)
        except Exception:                                           # noqa: BLE001
            a_s = tau_s = np.nan
    if not np.isfinite(tau_s):
        tau_s = TAIL_TAU_DEFAULT
        a_s = max(float(np.median(bs["y"][bs["x"] >= t_data - 20.0])), 0.0) * math.exp(t_data / tau_s)
    tx = np.arange(t_data + TAIL_PSEUDO_DT, T_AXIS_END + 1e-9, TAIL_PSEUDO_DT)
    return dict(T_data=t_data, tx=tx, ty=a_s * np.exp(-tx / tau_s), w_pseudo=float(np.median(1.0 / bs["s"])))


_PEN_CACHE = {}


def _penalty(t):
    """Gram matrix int B_i'' B_j'' dt of the cubic basis on knots t (0.05 ns Riemann sum), cached."""
    from scipy.interpolate import BSpline
    key = t.tobytes()
    if key not in _PEN_CACHE:
        g = np.arange(T_LEFT, T_AXIS_END, 0.05)
        d2 = BSpline(t, np.eye(t.size - 4), 3).derivative(2)(g)
        _PEN_CACHE[key] = (d2.T @ d2) * 0.05
        if len(_PEN_CACHE) > 8:
            _PEN_CACHE.pop(next(iter(_PEN_CACHE)))
    return _PEN_CACHE[key]


def fit_spline(bs, n_pool, pp, ref, lam):
    """Cubic B-spline fitted to the pooled points y/A (tau <= T_data, weight 1/sigma, sigma^2 =
    (sigma_pulse/A)^2 + (SIGMA_T_NS ref')^2 + (SHAPE_SIGMA ref)^2) plus the tail pseudo-points,
    penalised by lam * int f''^2.  The trigger-synchronous readout puts most pulses of a module at
    one sub-sample phase, so an unpenalised spline is free to oscillate between the populated
    phases; the penalty acts where the data leave it undetermined.  A failed solve doubles every
    knot spacing.  Returns dict(t, c, tail) on the pre-shift axis, or None."""
    from scipy.interpolate import BSpline
    tail = tail_extension(bs)
    dm = pp["tau"] <= tail["T_data"]
    xd = pp["tau"][dm]
    slope = (ref.eval(xd + 0.05) - ref.eval(xd - 0.05)) / 0.1
    sgd = np.sqrt(pp["sg"][dm] ** 2 + (SIGMA_T_NS * slope) ** 2 + (SHAPE_SIGMA * ref.eval(xd)) ** 2)
    x = np.concatenate([xd, tail["tx"]])
    y = np.concatenate([pp["y"][dm], tail["ty"]])
    w = np.concatenate([1.0 / sgd, np.full(tail["tx"].size, tail["w_pseudo"])])
    o = np.argsort(x, kind="stable")
    x, y, w = x[o], y[o], w[o]
    for scale in (1, 2, 4):
        t = np.concatenate([[T_LEFT] * 4, knot_plan(n_pool, scale), [T_AXIS_END] * 4])
        try:
            Bx = BSpline.design_matrix(x, t, 3).toarray()
            c = np.linalg.solve(Bx.T @ (Bx * (w * w)[:, None]) + lam * _penalty(t), Bx.T @ (w * w * y))
        except (ValueError, np.linalg.LinAlgError):
            continue
        if np.all(np.isfinite(c)):
            return dict(t=t, c=c, tail=tail)
    return None


def finalize_spline(t, c):
    """Unit peak and onset at 0 (where it last rises through ONSET_FRAC of the peak): dict(t, c)
    on the shifted axis, t_on (the shift) and vmax (the pre-normalisation peak)."""
    g = np.arange(-15.0, 80.0, 0.01)
    v = spline_values(t, c, g)
    ip = int(np.argmax(v))
    vmax = float(v[ip])
    below = np.flatnonzero(v[:ip] <= ONSET_FRAC * vmax)
    t_on = float(g[0])
    if below.size:
        j = int(below[-1])
        v0, v1 = v[j], v[j + 1]
        t_on = float(g[j] + ((ONSET_FRAC * vmax - v0) / (v1 - v0) if v1 > v0 else 0.0) * 0.01)
    return dict(t=t - t_on, c=c / vmax, t_on=t_on, vmax=vmax)


def spline_values(t, c, tt):
    from scipy.interpolate import BSpline
    return BSpline(t, c, 3)(tt)


def peak_times(t, c):
    """(t_pk, t50) of the shape (0 outside (0, T_END)) on the shifted axis (ns)."""
    g = np.arange(0.0, 60.0, 0.01)
    v = np.where((g > 0) & (g < T_END), spline_values(t, c, np.clip(g, t[3], t[-4])), 0.0)
    ip = int(np.argmax(v))
    i1 = int(np.argmax(v[:ip + 1] >= 0.5))
    i0 = max(i1 - 1, 0)
    return float(g[ip]), float(g[i0] + ((0.5 - v[i0]) / (v[i1] - v[i0]) if v[i1] > v[i0] else 0.0) * 0.01)


def rise_fall(f):
    """(10-90 % rise time, 90-10 % fall time), ns, of a shape f on T_GRID."""
    f = np.nan_to_num(np.asarray(f, dtype=np.float64))
    ip = int(np.argmax(f))

    def cross(level, step):
        j = ip
        while 0 <= j + step < f.size and f[j] >= level:
            j += step
        if f[j] >= level:
            return np.nan
        k = j - step
        return float(T_GRID[j] + (level - f[j]) / (f[k] - f[j]) * (T_GRID[k] - T_GRID[j]))
    pk = f[ip]
    return cross(0.9 * pk, -1) - cross(0.1 * pk, -1), cross(0.1 * pk, 1) - cross(0.9 * pk, 1)


def build_shape(ps, lam):
    """Template of one pulse set.  Pass 0 pools the pulses at (t0_cfd, amp).  Each pass fits the
    spline (penalty lam), puts t = 0 at its onset (finalize_spline) and refits the onset of every
    pulse against it within +-ALIGN_HW of its CFD onset (align_pulses); after ALIGN_PASSES passes
    comes the final fit.  Returns None if a fit fails, else dict(keep, pp, bs, sp, fin, ref, shift
    (final axis = CFD axis + shift), rms (per pass: std of t0 - t0_cfd over the pool), final (t0, A,
    chi2/dof, ok of every pulse against the final template))."""
    t0, A, keep = ps.cfd, ps.amp, np.ones(ps.n, dtype=bool)
    ref = LutShape(np.zeros(int(round(T_END / LUT_DT)) + 1))      # no shape yet: no detrending
    shift, rms = 0.0, []
    for it in range(ALIGN_PASSES + 1):
        pp = pool_points(ps.Y, ps.sig, t0, A, keep, ref)
        bs = bin_stats(pp, ref)
        sp = fit_spline(bs, int(keep.sum()), pp, ref, lam)
        if sp is None:
            return None
        fin = finalize_spline(sp["t"], sp["c"])
        shift += fin["t_on"]
        ref = LutShape.from_spline(fin["t"], fin["c"])
        if it == ALIGN_PASSES:
            break
        al = align_pulses(ps.Y, ps.sig, ps.cfd + shift, ref, ALIGN_HW)
        t0, A, keep = al["t0"], al["A"], al["ok"]
        if not keep.any():
            return None
        rms.append(float(np.std((t0 - ps.cfd)[keep])))
    t0 = t0 + fin["t_on"]
    cols = _align_columns(t0)
    final = _chi2_at(np.take_along_axis(ps.Y, cols, axis=1), cols, ps.sig, t0, ref, WIN)
    return dict(keep=keep, pp=pp, bs=bs, sp=sp, fin=fin, ref=ref, shift=shift, rms=rms, final=final)


def _osc_count(r):
    """Spurious extrema of the final shape of r: prominence >= OSC_PROM * peak, in OSC_WIN, farther
    than OSC_EXCL from the main peak (the largest value in OSC_PEAK_WIN)."""
    from scipy.signal import find_peaks
    f = spline_values(r["fin"]["t"], r["fin"]["c"], T_GRID)
    mp = (T_GRID >= OSC_PEAK_WIN[0]) & (T_GRID <= OSC_PEAK_WIN[1])
    ip = int(np.flatnonzero(mp)[np.argmax(f[mp])])
    w = np.flatnonzero((T_GRID >= OSC_WIN[0]) & (T_GRID <= OSC_WIN[1]))
    n = 0
    for sg in (f[w], -f[w]):
        j = w[find_peaks(sg, prominence=OSC_PROM * f[ip])[0]]
        n += int((np.abs(T_GRID[j] - T_GRID[ip]) > OSC_EXCL).sum())
    return n


def fit_module(ps, hs):
    """Module template: for a pool of >= LAMBDA_MIN_POOL pulses with held-out pulses hs,
    build_shape for each LAM_GRID value; hs are aligned to each final template (+-ALIGN_HW around
    their CFD onset) and the lambda with the lowest mean held-out chi2/dof (_chi2_full, over the
    pulses ok for every lambda) wins; otherwise LAMBDA_FALLBACK.  While the shape has spurious
    extrema (pools below OSC_POOL_MAX), lambda x10 up to SMOOTH_LAMBDA and rebuild.  Returns
    (template or None if it failed or is still wiggly, lambda used)."""
    if ps.n >= LAMBDA_MIN_POOL and hs.n:
        rs = {lam: r for lam in LAM_GRID if (r := build_shape(ps, lam)) is not None}
        if not rs:
            return None, np.nan
        als = {lam: align_pulses(hs.Y, hs.sig, hs.cfd + r["shift"], r["ref"], ALIGN_HW) for lam, r in rs.items()}
        ok = np.all([al["ok"] for al in als.values()], axis=0)
        if ok.any():
            hc = {lam: float(np.mean(_chi2_full(hs.Y, hs.sig, al["t0"], al["A"], rs[lam]["ref"])[ok]))
                  for lam, al in als.items()}
            lam = min(hc, key=hc.get)
        else:
            lam = LAMBDA_FALLBACK if LAMBDA_FALLBACK in rs else next(iter(rs))
        r = rs[lam]
    else:
        lam = LAMBDA_FALLBACK
        r = build_shape(ps, lam)
        if r is None:
            return None, np.nan
    n_osc = _osc_count(r) if ps.n < OSC_POOL_MAX else 0
    while n_osc > 0 and lam < SMOOTH_LAMBDA:
        lam = min(lam * 10, SMOOTH_LAMBDA)
        r2 = build_shape(ps, lam)
        if r2 is None:
            break
        r, n_osc = r2, _osc_count(r2)
    return (None if n_osc > 0 else r), lam


def energy_pools(D, e_beam):
    """Per module the candidate fit pools, tried in this order until a fit succeeds, as
    {module: [(pool indices, held-out indices, pool size before the cap, kind)]}: (1) "elastic", the pulses
    with e_rough >= POOL_EFRAC * e_beam, if there are at least POOL_MIN_EP; (2) "general", a random subset
    (rng per module from SEED) of the pulses with amplitude in [A_LO, A_HI], if there are at least
    POOL_MIN_GEN; (3) "elastic" with 3 to POOL_MIN_EP - 1 pulses; (4) "single", the pulse of largest
    amplitude.  Elastic: the next int(n_pool * HOLDOUT_FRAC / (1 - HOLDOUT_FRAC) + 0.5) pulses below the
    pool in e_rough are held out, and a pool above POOL_MAX is cut to a random subset of POOL_MAX (the
    held-out count follows the used pool).  General: HOLDOUT_FRAC of the band is not used for the pool,
    and the next pulses of it, as many as for the elastic rule, are held out.  A module without pulses has
    no entry."""
    E, mi, amp = D["e_rough"].astype(np.float64), D["midx"], D["amp"]
    plans = {}
    order = np.lexsort((-E, mi))                  # by module, e_rough descending
    for g in np.split(order, np.flatnonzero(np.diff(mi[order])) + 1):
        if g.size == 0:
            continue
        m = int(mi[g[0]])
        rng = np.random.default_rng([SEED, m])
        n = int((E[g] >= POOL_EFRAC * e_beam).sum())
        band = g[(amp[g] >= A_LO) & (amp[g] <= A_HI)]
        single = (np.array([g[np.argmax(amp[g])]]), np.array([], dtype=g.dtype), 1, "single")
        used = min(n, POOL_MAX)
        elastic = (np.sort(rng.choice(g[:n], POOL_MAX, replace=False) if n > POOL_MAX else g[:n]),
                   np.sort(g[n:n + int(used * HOLDOUT_FRAC / (1.0 - HOLDOUT_FRAC) + 0.5)]), n, "elastic")
        plan = [elastic] if n >= POOL_MIN_EP else []
        if band.size >= POOL_MIN_GEN:
            band = rng.permutation(band)
            used = min(POOL_MAX, int(band.size * (1.0 - HOLDOUT_FRAC) + 0.5))
            plan += [(np.sort(band[:used]), np.sort(band[used:used + int(used * HOLDOUT_FRAC / (1.0 - HOLDOUT_FRAC) + 0.5)]),
                      int(band.size), "general")]
        plan += [elastic] if 3 <= n < POOL_MIN_EP else []
        plan += [single]
        plans[m] = plan
    return plans


def _gsigma(R, nv, g, ng):
    """Noise-subtracted robust spread per group: the s >= 0 with
    median((R - median R)^2 / (s^2 + nv)) = CHI2_1_MED, i.e. 1.4826 MAD with the per-point noise
    variance nv removed in quadrature; NaN below SIGMA_BIN_MIN points.  Bisection: that median
    exceeds CHI2_1_MED exactly when more than half of the points do."""
    med, cnt = _grouped_median(R, g, ng)
    r2 = (R - med[g]) ** 2
    lo, hi = np.zeros(ng), np.full(ng, r2.max() / CHI2_1_MED)
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        low = np.bincount(g, r2 > CHI2_1_MED * (mid[g] + nv), ng) > 0.5 * cnt   # s^2 = mid too small
        lo, hi = np.where(low, mid, lo), np.where(low, hi, mid)
    return np.where(cnt >= SIGMA_BIN_MIN, np.sqrt(hi), np.nan)


def shape_sigma(Y, t0, A, t, c, t_end, n_ped):
    """Pulse-to-pulse spread of pulses Y (pedestal-subtracted) about a template (spline t, c) on T_GRID (fraction
    of the peak): the pulses at their onset t0 and amplitude A give R = y / A - f(tau) for every sample with
    T_GRID[0] <= tau <= t_end; per SIGMA_BIN_NS bin sigma = _gsigma(R, nv), nv = v_ped (1 + 1 /
    n_ped) / A^2, v_ped = mean over the pulses of the variance of the samples at tau < PRE_ONSET_NS
    (+ 1 / n_ped: the pedestal error).  Linear interpolation between the bin centres; NaN outside."""
    tau = TS[None, :] - t0[:, None]
    pre = tau < PRE_ONSET_NS
    npre = pre.sum(1)
    mean_pre = (Y * pre).sum(1) / npre
    v_ped = float(np.mean(((Y - mean_pre[:, None]) ** 2 * pre).sum(1) / (npre - 1)))
    m = (tau >= T_GRID[0]) & (tau <= t_end)
    ts = tau[m]
    R = (Y / A[:, None])[m] - spline_values(t, c, ts)
    nv = np.broadcast_to((v_ped * (1.0 + 1.0 / n_ped) / A ** 2)[:, None], tau.shape)[m]
    b = np.floor((ts - T_GRID[0]) / SIGMA_BIN_NS).astype(np.int64)
    sg = _gsigma(R, nv, b, int(b.max()) + 1)
    ctr = T_GRID[0] + (np.arange(sg.size) + 0.5) * SIGMA_BIN_NS
    ok = np.isfinite(sg)
    if not ok.any():
        return np.full(T_GRID.size, np.nan, dtype=np.float32)
    return np.interp(T_GRID, ctr[ok], sg[ok], left=np.nan, right=np.nan).astype(np.float32)


def edge_times(Y):
    """Per pulse (rows of the baseline-subtracted samples Y) the 10-90 % rise time and 90-10 % fall time (ns)
    relative to the pulse's own peak sample, linear interpolation between samples; the rising crossings are the
    last samples below the level before the peak.  NaN without a crossing."""
    n, ns = Y.shape
    r, col = np.arange(n), np.arange(ns)[None, :]
    ip = Y.argmax(1)
    pk = Y[r, ip]

    def rise_at(f):
        lev = f * pk
        j = np.clip(np.where((Y < lev[:, None]) & (col < ip[:, None]), col, -1).max(1), 0, ns - 2)
        y0, y1 = Y[r, j], Y[r, j + 1]
        return np.where((Y[r, j] < lev) & (j < ip), TS[j] + (lev - y0) / (y1 - y0) * DT_NS, np.nan)

    def fall_at(f):
        lev = f * pk
        k = np.where((Y < lev[:, None]) & (col > ip[:, None]), col, ns).min(1)
        kk = np.clip(k, 1, ns - 1)
        y0, y1 = Y[r, kk - 1], Y[r, kk]
        return np.where(k < ns, TS[kk - 1] + (y0 - lev) / (y0 - y1) * DT_NS, np.nan)
    return rise_at(0.9) - rise_at(0.1), fall_at(0.1) - fall_at(0.9)


def band_stats(D, midx, t, c, shift, t_end, n_ped):
    """General band pool of a module (pulses with A_LO <= A <= A_HI, at most POOL_MAX random ones, own rng from
    SEED) aligned to the module's final template (spline t, c; shift = final axis - CFD axis, t_end its data
    end): shape_sigma(t) of the pool about the template, and the pulse-to-pulse jitter (1.4826 MAD) and median
    of the 10-90 % rise and 90-10 % fall times of the raw samples; a_band = the median fitted amplitude (ADC) of the pool,
    the amplitude at which shape_sigma is measured.  None below SIGMA_MIN_POOL aligned pulses."""
    sel = np.flatnonzero((D["midx"] == midx) & (D["amp"] >= A_LO) & (D["amp"] <= A_HI))
    if sel.size > POOL_MAX:
        sel = np.sort(np.random.default_rng([SEED, int(midx), 1]).choice(sel, POOL_MAX, replace=False))
    if sel.size < SIGMA_MIN_POOL:
        return None
    ps = PulseSet(D, sel)
    al = align_pulses(ps.Y, ps.sig, ps.cfd + shift, LutShape.from_spline(t, c), ALIGN_HW)
    ok = al["ok"]
    if ok.sum() < SIGMA_MIN_POOL:
        return None
    out = dict(sigma=shape_sigma(ps.Y[ok], al["t0"][ok], al["A"][ok], t, c, t_end, n_ped),
               a_band=float(np.median(al["A"][ok])))
    for key, x in zip(("rise", "fall"), edge_times(ps.Y[ok])):
        x = x[np.isfinite(x)]
        med = float(np.median(x)) if x.size >= SIGMA_MIN_POOL else np.nan
        out[key + "_med"] = med
        out[key + "_jit"] = 1.4826 * float(np.median(np.abs(x - med))) if x.size >= SIGMA_MIN_POOL else np.nan
    return out


def template_rows(fin, bs, sp):
    """shape, shape_err (statistical error of the BIN_W bin median, NaN beyond T_data) and resid
    (bin median - spline) on T_GRID, float32."""
    shape = spline_values(fin["t"], fin["c"], T_GRID)
    d = bs["x"] <= sp["tail"]["T_data"]
    xd = bs["x"][d] - fin["t_on"]
    err = np.interp(T_GRID, xd, bs["s"][d] / fin["vmax"], left=np.nan, right=np.nan)
    res = np.interp(T_GRID, xd, bs["y"][d] / fin["vmax"] - spline_values(fin["t"], fin["c"], xd),
                    left=np.nan, right=np.nan)
    return shape.astype(np.float32), err.astype(np.float32), res.astype(np.float32)


def process_module(ctx, midx):
    """One PbWO4 module: template from the first of its candidate pools (energy_pools) that fits,
    shape_sigma_elastic for elastic pools of >= SIGMA_MIN_POOL; fallback 2 (run-type average) without pulses or when
    every pool gives a failed or wiggly fit (n_pool, single and the other fields then describe the first pool)."""
    D = ctx["D"]
    plan = ctx["plans"].get(midx, [])
    res = dict(midx=midx, fallback=2, n_pool=0, n_pool_total=0, single=0, pool_general=0,
               lam=np.nan, chi2=np.nan, rms=[np.nan] * ALIGN_PASSES, sigma_el=None)
    if plan:
        res.update(n_pool=int(plan[0][0].size), n_pool_total=plan[0][2], single=int(plan[0][3] == "single"))
    for pool, held, total, kind in plan:
        ps = PulseSet(D, pool)
        r, lam = fit_module(ps, PulseSet(D, held))
        if r is not None:
            break
    else:
        return res
    fin = r["fin"]
    shape, err, resid = template_rows(fin, r["bs"], r["sp"])
    keep = r["keep"] & r["final"]["ok"]
    rng = np.random.default_rng([SEED, int(midx)])
    sel = np.sort(rng.choice(r["pp"]["tau"].size, min(3000, r["pp"]["tau"].size), replace=False))
    res.update(fallback=0, n_pool=int(pool.size), n_pool_total=total, single=int(kind == "single"),
               pool_general=int(kind == "general"), lam=lam, chi2=float(np.median(r["final"]["chi2dof"][keep])),
               shape=shape, err=err, resid=resid, t=fin["t"], c=fin["c"], cum_shift=r["shift"], rms=r["rms"],
               t_end=r["sp"]["tail"]["T_data"] - fin["t_on"],
               sigma_el=(shape_sigma(ps.Y[keep], r["final"]["t0"][keep], r["final"]["A"][keep], fin["t"], fin["c"],
                                     r["sp"]["tail"]["T_data"] - fin["t_on"], ctx["n_ped"])
                         if kind == "elastic" and pool.size >= SIGMA_MIN_POOL else None),
               pts=(np.float32(r["pp"]["tau"][sel] - fin["t_on"]), np.float32(r["pp"]["y"][sel] / fin["vmax"])),
               kept_gidx=pool[keep], kept_t0=r["final"]["t0"][keep], kept_A=r["final"]["A"][keep])
    return res


def run_average(D, results):
    """Run-type average: build_shape (lambda SMOOTH_LAMBDA) of <= 20 kept pool pulses of each
    own-template module that is not a general-pool one (elastic shapes only)."""
    rng = np.random.default_rng([SEED, 987654321])
    use = [r for r in results if r["fallback"] == 0 and not r["pool_general"]]
    if not use:
        sys.exit("[ERROR] no module with an own template: cannot build the run-type average")
    g = np.concatenate([rng.choice(r["kept_gidx"], min(20, r["kept_gidx"].size), replace=False) for r in use])
    return build_shape(PulseSet(D, np.sort(g)), SMOOTH_LAMBDA), g.size


def load_candidates(files, fields=("midx", "samples", "ped", "rms", "amp", "t0_cfd", "e_rough", "event_number")):
    """Concatenated select files: the fields plus run, split, cutflow (summed), wave config, inputs,
    dead (modules dead in every file), e_beam (MeV, the same in every file)."""
    out = {k: [] for k in fields}
    runs, splits, inputs, cutflow, wave, dead, ebeam = [], [], [], 0, None, None, set()
    for f in files:
        z = np.load(f, allow_pickle=False)
        meta = json.loads(str(z["meta"]))
        if meta.get("schema") != SCHEMA:
            sys.exit(f"[ERROR] {f}: schema {meta.get('schema')} != {SCHEMA}")
        if wave is not None and meta["wave_config"] != wave:
            sys.exit(f"[ERROR] {f}: wave config differs from the first file")
        wave = meta["wave_config"]
        for k in fields:
            out[k].append(z[k])
        runs.append(np.full(z["midx"].size, int(z["run"]), dtype=np.int64))
        splits.append(z["split"].astype(np.int64))
        cutflow = cutflow + z["cutflow"]
        ebeam.add(round(meta["e_beam_mev"], 1))
        dead = set(meta["dead_modules"]) if dead is None else dead & set(meta["dead_modules"])
        inputs.append(dict(path=os.path.abspath(f), sha1=_sha1(f), run=int(z["run"]), n_pulses=int(z["midx"].size)))
    if len(ebeam) != 1:
        sys.exit(f"[ERROR] the select files have different beam energies: {sorted(ebeam)} MeV")
    D = {k: np.concatenate(v) for k, v in out.items()}
    D.update(run=np.concatenate(runs), split=np.concatenate(splits), cutflow=cutflow, wave=wave,
             inputs=inputs, dead=dead, e_beam=ebeam.pop())
    return D


# writers =====================================================================
def write_root(path, arr):
    """ROOT file (uproot): TTree `templates` (one entry per PbWO4 module: module, module_id, valid,
    fallback, single, pool_general, sigma_borrowed, a_band, n_pulses, lambda, chi2, rise_ns, fall_ns, t_pk, iint_wa, knots[], coeffs[],
    rise_jitter_ns, fall_jitter_ns, rise_med_ns, fall_med_ns, shape / shape_err / shape_sigma / shape_sigma_elastic on t_grid), TH1D shape/<module> (0.1 ns bins, content shape,
    error shape_sigma, 0 where undefined) and TObjString meta."""
    import awkward as ak
    import uproot
    from uproot.writing.identify import to_TAxis, to_TH1x
    nr, nt = arr["shape"].shape
    rf = np.array([rise_fall(f) for f in arr["shape"]])
    cols = dict(module=ak.Array([str(x) for x in arr["names"]]),
                module_id=arr["module_id"].astype(np.int32), valid=arr["valid"].astype(np.int32),
                fallback=arr["fallback"].astype(np.int32), single=arr["single"].astype(np.int32),
                pool_general=arr["pool_general"].astype(np.int32), sigma_borrowed=arr["sigma_borrowed"].astype(np.int32), n_pulses=arr["n_pool"],
                **{"lambda": arr["lam_module"]}, chi2=arr["chi2_ndf"], rise_ns=rf[:, 0], fall_ns=rf[:, 1],
                t_pk=arr["t_pk"], iint_wa=arr["iint_wa"],
                knots=ak.Array([arr["spline_t"][i, :arr["spline_nt"][i]] for i in range(nr)]),
                coeffs=ak.Array([arr["spline_c"][i, :arr["spline_nc"][i]] for i in range(nr)]),
                shape=arr["shape"], shape_err=arr["shape_err"], shape_sigma=arr["shape_sigma"], a_band=arr["a_band"],
                shape_sigma_elastic=arr["shape_sigma_elastic"],
                **{k: arr[k] for k in ("rise_jitter_ns", "fall_jitter_ns", "rise_med_ns", "fall_med_ns")})
    types = {k: (v.type if isinstance(v, ak.Array) else np.dtype((v.dtype, v.shape[1:])) if v.ndim > 1
                 else v.dtype) for k, v in cols.items()}
    lo, hi = T_GRID[0] - GRID_DT / 2, T_GRID[-1] + GRID_DT / 2
    with uproot.recreate(path) as f:
        f.mktree("templates", types)
        f["templates"].extend(cols)
        for i, name in enumerate(arr["names"]):
            y = np.nan_to_num(arr["shape"][i].astype(np.float64))
            e = np.nan_to_num(arr["shape_sigma"][i].astype(np.float64))
            f[f"shape/{name}"] = to_TH1x(
                str(name), f"{name} unit-peak template;t (ns after onset);shape",
                np.concatenate([[0.0], y, [0.0]]), 0.0, y.sum(), 0.0, 0.0, 0.0,
                np.concatenate([[0.0], e * e, [0.0]]), to_TAxis("xaxis", "", nt, lo, hi))
        f["meta"] = str(arr["meta"])


def pool_map_figure(arr, run_type):
    """PbWO4 heat map of n_pool per module (HYCAL_MAP geometry); run-type-average modules hatched
    with x, general-pool modules with /, single-pulse modules outlined in red, dead modules white."""
    import matplotlib.pyplot as plt
    # plot style: large bold axis labels, large legends
    plt.rcParams.update({
        "axes.labelsize": 16, "axes.labelweight": "bold", "axes.titlesize": 16, "axes.titleweight": "bold",
        "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14, "axes.linewidth": 1.8,
        "xtick.major.width": 1.5, "ytick.major.width": 1.5, "xtick.major.size": 6, "ytick.major.size": 6,
        "lines.linewidth": 2.0, "figure.titlesize": 18, "figure.titleweight": "bold"})
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle
    geo = {m["n"]: m["geo"] for m in json.load(open(HYCAL_MAP))}
    rect = [Rectangle((geo[n]["x"] - geo[n]["sx"] / 2, geo[n]["y"] - geo[n]["sy"] / 2), geo[n]["sx"], geo[n]["sy"])
            for n in (str(x) for x in arr["names"])]
    fig, ax = plt.subplots(figsize=(11, 10))
    pc = PatchCollection(rect, cmap="viridis", edgecolor="0.3", lw=0.2)
    pc.set_array(arr["n_pool"].astype(np.float64))
    ax.add_collection(pc)
    dead = arr["dead"] == 1
    fb, sg = (arr["fallback"] == 2) & ~dead, (arr["single"] == 1) & ~dead
    gen = arr["pool_general"] == 1
    ax.add_collection(PatchCollection([r for r, f in zip(rect, dead) if f], facecolor="w", edgecolor="0.5", lw=0.5))
    ax.add_collection(PatchCollection([r for r, f in zip(rect, fb) if f], facecolor="none", edgecolor="k",
                                      hatch="xxx", lw=0.3))
    ax.add_collection(PatchCollection([r for r, f in zip(rect, gen) if f], facecolor="none", edgecolor="k",
                                      hatch="///", lw=0.3))
    ax.add_collection(PatchCollection([r for r, f in zip(rect, sg) if f], facecolor="none", edgecolor="r", lw=1.0))
    ax.set(xlim=(-MAP_LIM, MAP_LIM), ylim=(-MAP_LIM, MAP_LIM), aspect="equal", xlabel="x (mm)", ylabel="y (mm)")
    cb = fig.colorbar(pc, ax=ax)
    cb.set_label("pulses in the fit pool", fontsize=16, fontweight="bold")
    cb.ax.tick_params(labelsize=14)
    ax.set_title(f"{run_type}: pool per module: elastic e_rough >= {POOL_EFRAC:g} E_beam,\nor general {A_LO:g}-{A_HI:g} ADC band where elastic < {POOL_MIN_EP};\nx hatch: run-type average "
                 f"({int(fb.sum())}); / hatch: general pool ({int(gen.sum())});\nred outline: single pulse ({int(sg.sum())}); white: dead ({int(dead.sum())})",
                 fontsize=12)
    fig.tight_layout()
    return fig


def example_pulses(results):
    """PDF examples, one page per module with EX_PER_MODULE random pool pulses (all of them if fewer): all
    own-template modules of the innermost ring around the beam hole (ring = Chebyshev distance max(|x|, |y|)
    from the centre in module pitches, innermost = smallest), the first EX_RING2 (by module id) of the second ring, and the
    EX_BORROWED outermost modules whose shape_sigma is borrowed (sigma_borrowed, i.e. fewer than SIGMA_MIN_POOL
    band pulses).  Dead modules (and fallback-2 ones, which keep no pulses) are skipped.
    Returns [(label, module result, [indices into its kept pulses])]."""
    geo = {m["n"]: m["geo"] for m in json.load(open(HYCAL_MAP))}
    ok = [r for r in results if r["fallback"] == 0 and not r["dead"]]
    cheb = lambda r: int(max(abs(geo[r["name"]]["x"]), abs(geo[r["name"]]["y"])) // 20.77)   # noqa: E731  ring index: floor of the distance in module pitches (20.77 mm; centres jitter by < 0.1 mm)
    rings = sorted({cheb(r) for r in ok})
    by_id = lambda rs: sorted(rs, key=lambda r: r["module_id"])   # noqa: E731
    brw = [r for r in ok if r["sigma_borrowed"]]
    groups = [("innermost ring", by_id([r for r in ok if cheb(r) == rings[0]])),
              ("second ring", by_id([r for r in ok if cheb(r) == rings[1]])[:EX_RING2]),
              ("outer, borrowed shape_sigma", by_id(sorted(brw, key=lambda r: -cheb(r))[:EX_BORROWED]))]
    return [(lab, r, [int(j) for j in np.sort(np.random.default_rng([SEED, r["midx"]]).choice(
        r["kept_gidx"].size, min(EX_PER_MODULE, r["kept_gidx"].size), replace=False))]) for lab, rs in groups for r in rs]


def write_pdf(path, rows, arr, run_type, map_png, D):
    """PDF: a summary page (rise / fall times, pool sizes, array median and 90th percentile of
    shape_sigma(t) of the modules with their own band-pool sigma, solid, and of shape_sigma_elastic, dashed;
    pulse rise / fall jitter), the n_pool heat map (also saved as map_png), then 3 x 3 modules per page: the
    pooled points (y/A, a random 3000), the spline with a +-shape_sigma band and below it the
    binned residual with the same band (dead modules left out).  Fallback-2 modules show the
    run-type average.  Last, the example pulses (example_pulses), one page per listed module, EX_PER_MODULE pulses in EX_COLS columns: samples - pedestal
    and the module's template at the pulse's fitted onset and amplitude."""
    import matplotlib.pyplot as plt
    # plot style: large bold axis labels, large legends
    plt.rcParams.update({
        "axes.labelsize": 16, "axes.labelweight": "bold", "axes.titlesize": 16, "axes.titleweight": "bold",
        "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14, "axes.linewidth": 1.8,
        "xtick.major.width": 1.5, "ytick.major.width": 1.5, "xtick.major.size": 6, "ytick.major.size": 6,
        "lines.linewidth": 2.0, "figure.titlesize": 18, "figure.titleweight": "bold"})
    from matplotlib.backends.backend_pdf import PdfPages
    fb = arr["fallback"]
    own = fb == 0
    rf = np.array([rise_fall(f) for f in arr["shape"]])
    v = T_GRID <= 200.0
    tg = T_GRID[v]
    with PdfPages(path) as pdf:
        fig, ax = plt.subplot_mosaic([["rise", "fall", "pool"], ["sig", "sig", "txt"]], figsize=(17, 14), layout="constrained",
                                       height_ratios=[2, 3])
        for key, x, lab in (("rise", rf[own, 0], "10-90 % rise time (ns)"), ("fall", rf[own, 1], "90-10 % fall time (ns)"),
                            ("pool", np.log10(arr["n_pool"][own]), "log10 pool size (used)")):
            ax[key].hist(x[np.isfinite(x)], bins=40)
            ax[key].set(xlabel=lab, ylabel="modules (own template)")
        sg = 100.0 * arr["shape_sigma"][own & (arr["sigma_borrowed"] == 0)][:, v].astype(np.float64)
        se = 100.0 * arr["shape_sigma_elastic"][own][:, v].astype(np.float64)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)       # all-NaN columns
            for k, (q, lab) in enumerate(((50, "median"), (90, "90th percentile"))):
                ax["sig"].plot(tg, np.nanpercentile(sg, q, axis=0), c=f"C{k}", label=f"band pool (matched filter), {lab}")
                ax["sig"].plot(tg, np.nanpercentile(se, q, axis=0), c=f"C{k}", ls="--", label=f"elastic pool, {lab}")
        ax["sig"].set(xlim=(-10, 200), xlabel="t (ns after onset)", ylabel="shape_sigma (% of peak)")
        ax["sig"].legend()
        ax["txt"].axis("off")
        ax["txt"].text(0, 1, "\n".join([
            f"run type {run_type}", f"{len(fb)} PbWO4 modules",
            f"own template {int(own.sum())}", f"  (general pool {int(arr['pool_general'].sum())})", f"run-type average {int((fb == 2).sum())}",
            f"dead (excluded) {int(arr['dead'].sum())}",
            f"single pulse {int(arr['single'].sum())}", f"shape_sigma borrowed {int(arr['sigma_borrowed'][own].sum())}",
            f"rise median {np.nanmedian(rf[own, 0]):.2f} ns", f"fall median {np.nanmedian(rf[own, 1]):.1f} ns",
            f"shape_sigma_elastic own {int(np.isfinite(se).any(1).sum())}",
            f"median chi2/dof {np.nanmedian(arr['chi2_ndf'][own]):.2f}",
            "pulse-to-pulse jitter (band pool):",
            f"  rise {np.nanmedian(arr['rise_jitter_ns'][own]):.2f} ns (median of modules)",
            f"  fall {np.nanmedian(arr['fall_jitter_ns'][own]):.2f} ns",
            f"  rise time median {np.nanmedian(arr['rise_med_ns'][own]):.2f} ns",
            f"  fall time median {np.nanmedian(arr['fall_med_ns'][own]):.2f} ns"]), va="top", family="monospace", fontsize=14)
        fig.suptitle(f"{run_type}: cubic B-spline templates")
        pdf.savefig(fig)
        plt.close(fig)
        fig = pool_map_figure(arr, run_type)
        fig.savefig(map_png, dpi=150)
        pdf.savefig(fig)
        plt.close(fig)
        live = [i for i in range(len(rows)) if not arr["dead"][i]]
        for p0 in range(0, len(live), 9):
            fig = plt.figure(figsize=(24, 18))
            gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.2, left=0.05, right=0.99, top=0.96, bottom=0.06)
            for k, i in enumerate(live[p0:p0 + 9]):
                r = rows[i]
                sub = gs[k // 3, k % 3].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
                a = fig.add_subplot(sub[0])
                b = fig.add_subplot(sub[1], sharex=a)
                shp = arr["shape"][i][v].astype(float)
                err = np.nan_to_num(arr["shape_sigma"][i][v].astype(float))
                if r.get("pts") is not None:
                    a.scatter(*r["pts"], s=1, c="0.6", lw=0, rasterized=True)
                a.fill_between(tg, shp - err, shp + err, color="C1", alpha=0.4, lw=0, rasterized=True)
                a.plot(tg, shp, c="C0", lw=0.8, rasterized=True)
                flag = " single" * arr["single"][i] + " general" * arr["pool_general"][i] + " sigma borrowed" * arr["sigma_borrowed"][i]
                a.set_title(f"{r['name']}  n {r['n_pool']}/{r['n_pool_total']}  lam {r['lam']:g}  chi2 {r['chi2']:.2f}{flag}" if fb[i] == 0
                            else f"{r['name']}  run-type average", fontsize=11, pad=3)
                a.set_xlim(-10, 200)
                a.tick_params(labelsize=12, labelbottom=False)
                b.fill_between(tg, -err, err, color="C1", alpha=0.4, lw=0, rasterized=True)
                b.plot(tg, r["resid"][v], c="k", lw=0.5, rasterized=True)
                b.axhline(0, c="C0", lw=0.5)
                b.set_ylim(-0.03, 0.03)
                b.tick_params(labelsize=12)
            fig.supxlabel("t (ns after onset); lower panels: bin median - spline", fontsize=16, fontweight="bold")
            pdf.savefig(fig, dpi=100)
            plt.close(fig)
        for label, r, js in example_pulses(rows):
            nr = -(-EX_PER_MODULE // EX_COLS)
            fig, axs = plt.subplots(nr, EX_COLS, figsize=(6 * EX_COLS, 5 * nr + 1), sharey=False, layout="constrained", squeeze=False)
            for a, j in zip(axs.ravel(), js):
                g, t0, A = r["kept_gidx"][j], r["kept_t0"][j], r["kept_A"][j]
                tt = np.arange(0.0, TS[-1] + DT_NS, 0.5)
                a.plot(tt, A * np.where(tt > t0, spline_values(r["t"], r["c"], tt - t0), 0.0), c="C0", lw=0.8)
                a.plot(TS, D["samples"][g] - D["ped"][g], "k.", ms=3)
                a.set_xlim(0.0, 400.0)
                a.set_title(f"run {D['run'][g]} ev {D['event_number'][g]}\nA {A:.0f} ADC (fit), "
                            f"e_rough {D['e_rough'][g]:.0f} MeV", fontsize=11)
                a.tick_params(labelsize=12)
            for a in axs.ravel()[len(js):]:
                a.axis("off")
            fig.suptitle(f"{r['name']} ({label}): {len(js)} random pool pulses of {r['n_pool']}"
                         + (f" (fewer than {EX_PER_MODULE} available)" if len(js) < EX_PER_MODULE else ""))
            fig.supxlabel("t (ns)", fontsize=16, fontweight="bold")
            fig.supylabel("samples - pedestal (ADC)", fontsize=16, fontweight="bold")
            pdf.savefig(fig)
            plt.close(fig)


# driver ======================================================================
def build_templates(run_type, pulse_glob, out_stem):
    """`build`: pool the select files matching pulse_glob and fit the templates of run_type;
    writes <out_stem>.npz / .pdf / .root and n_templates_map.png next to it."""
    t_start = time.monotonic()
    files = sorted(glob.glob(str(pulse_glob)))
    if not files:
        sys.exit(f"[ERROR] no files match {pulse_glob}")
    stem = Path(out_stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    D = load_candidates(files)
    runs = sorted(set(int(x) for x in D["run"]))
    plans = energy_pools(D, D["e_beam"])
    print(f"[load] {len(files)} files, runs {runs}, {D['midx'].size} pulses, E_beam {D['e_beam']:.1f} MeV; "
          f"first-choice fit pools {sum(p[0][0].size for p in plans.values())}, held out {sum(p[0][1].size for p in plans.values())}", flush=True)
    p = setup_detector(runs[0])
    mods = [p.hycal.module(i) for i in range(p.hycal.module_count())]
    names = [m.name for m in mods]
    ctx = dict(D=D, plans=plans, n_ped=int(D["wave"]["ped_nsamples"]))
    pwo4 = [i for i, m in enumerate(mods) if m.is_pwo4()]
    results = []
    for k, mi in enumerate(pwo4):
        results.append(process_module(ctx, mi))
        if (k + 1) % 100 == 0:
            print(f"  [modules] {k + 1}/{len(pwo4)}  {time.monotonic() - t_start:.0f} s", flush=True)

    av, n_avg = run_average(D, results)
    afin = av["fin"]
    ashape, aerr, ares = template_rows(afin, av["bs"], av["sp"])
    for r in results:
        r["name"], r["module_id"] = names[r["midx"]], int(mods[r["midx"]].id)
        r["dead"] = int(r["name"] in D["dead"])
        if r["fallback"] == 2:
            r.update(shape=ashape, err=aerr, resid=ares, t=afin["t"], c=afin["c"], cum_shift=av["shift"], pts=None,
                     t_end=av["sp"]["tail"]["T_data"] - afin["t_on"])
        b = None if r["dead"] else band_stats(D, r["midx"], r["t"], r["c"], r["cum_shift"], r["t_end"], ctx["n_ped"])
        r["band"] = b or dict(sigma=None, a_band=np.nan, rise_med=np.nan, fall_med=np.nan, rise_jit=np.nan, fall_jit=np.nan)
        r["sigma_borrowed"] = int(b is None)
    own_sig = [r["band"]["sigma"] for r in results if r["band"]["sigma"] is not None]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)           # all-NaN columns
        sig_med = (np.nanmedian(np.stack(own_sig), axis=0).astype(np.float32) if own_sig
                   else np.full(T_GRID.size, np.nan, dtype=np.float32))
    own_a = [r["band"]["a_band"] for r in results if r["band"]["sigma"] is not None]
    a_med = float(np.median(own_a)) if own_a else np.nan
    for r in results:
        r["a_band"] = a_med if r["band"]["sigma"] is None else r["band"]["a_band"]
        r["sigma"] = sig_med if r["band"]["sigma"] is None else r["band"]["sigma"]
        r["sigma_elastic"] = np.full(T_GRID.size, np.nan, dtype=np.float32) if r["sigma_el"] is None else r["sigma_el"]
        sh = LutShape.from_spline(r["t"], r["c"])
        r["iint_wa"] = wa_integral(sh, D["wave"])
        r["t_pk"], r["t50"] = peak_times(r["t"], r["c"])

    g = lambda k, dt=np.float64: np.array([r[k] for r in results], dtype=dt)   # noqa: E731
    nk, nc = max(r["t"].size for r in results), max(r["c"].size for r in results)
    pad = lambda key, n: np.array([np.pad(r[key], (0, n - r[key].size), constant_values=np.nan)   # noqa: E731
                                   for r in results])
    meta = dict(
        schema="spline_templates_v3", run_type=run_type, runs=runs, inputs=D["inputs"], e_beam_mev=D["e_beam"],
        wave_config=D["wave"],
        build=_constants(["POOL_EFRAC", "POOL_MIN_EP", "A_LO", "A_HI", "POOL_MIN_GEN", "HOLDOUT_FRAC", "LAMBDA_MIN_POOL", "SIGMA_MIN_POOL", "ALIGN_PASSES",
                          "ALIGN_HW", "SEED", "WIN", "POOL_TAU_MIN", "BIN_W", "SHAPE_SIGMA", "SIGMA_T_NS",
                          "LAM_GRID", "LAMBDA_FALLBACK", "SMOOTH_LAMBDA", "OSC_PROM", "OSC_WIN", "OSC_EXCL",
                          "OSC_PEAK_WIN", "OSC_POOL_MAX", "POOL_MAX", "ONSET_FRAC", "SIGMA_BIN_NS", "SIGMA_BIN_MIN", "PRE_ONSET_NS",
                          "EX_RING2", "EX_BORROWED", "EX_PER_MODULE", "EX_COLS"]),
        notes=dict(
            shape_sigma="pulse-to-pulse spread about the template (fraction of the peak) of the module's general band pool "
                        "(A_LO <= A <= A_HI, <= POOL_MAX random pulses) aligned to the final template, for every module whichever "
                        "pool the template came from; modules with < SIGMA_MIN_POOL band pulses (sigma_borrowed = 1) carry the "
                        "run-type median of the others",
            a_band="median fitted amplitude (ADC) of the same band pool (the run-type median of the own-sigma modules where "
                   "sigma_borrowed = 1): shape_sigma is measured at a_band, the matched filter scales it by sqrt(a_band / A)",
            shape_sigma_elastic="the same spread for the elastic fit pool itself; NaN unless the template is from an elastic "
                                "pool of >= SIGMA_MIN_POOL pulses",
            jitter="rise_jitter_ns / fall_jitter_ns: 1.4826 x median absolute deviation, rise_med_ns / fall_med_ns: median, of "
                   "the per-pulse 10-90 % rise and 90-10 % fall times (ns) of the same band pulses (baseline-subtracted raw samples, "
                   "linear interpolation between samples, relative to the pulse's own peak sample); NaN below SIGMA_MIN_POOL pulses"),
        dead_modules=sorted(D["dead"]), average_pulses=int(n_avg), git=_git_short(),
        created=datetime.now().isoformat(timespec="seconds"))
    rms = np.array([r["rms"] for r in results], dtype=np.float64)        # [module, pass]
    arr = dict(
        names=np.array(names)[g("midx", np.int64)], module_id=g("module_id", np.int64),
        valid=1 - g("dead", np.int8), dead=g("dead", np.int8), fallback=g("fallback", np.int8),
        single=g("single", np.int8), pool_general=g("pool_general", np.int8), sigma_borrowed=g("sigma_borrowed", np.int8), n_pool=g("n_pool", np.int32),
        n_pool_total=g("n_pool_total", np.int32),
        lam_module=g("lam"), chi2_ndf=g("chi2"), align_rms_ns=rms[:, -1],
        t_grid=T_GRID, shape=np.stack([r["shape"] for r in results]),
        shape_err=np.stack([r["err"] for r in results]), shape_sigma=np.stack([r["sigma"] for r in results]), a_band=g("a_band"),
        shape_sigma_elastic=np.stack([r["sigma_elastic"] for r in results]),
        **{f"{k}_ns": np.array([r["band"][v] for r in results]) for k, v in (
            ("rise_jitter", "rise_jit"), ("fall_jitter", "fall_jit"), ("rise_med", "rise_med"), ("fall_med", "fall_med"))},
        spline_k=np.int64(3), spline_t=pad("t", nk), spline_c=pad("c", nc),
        spline_nt=np.array([r["t"].size for r in results], dtype=np.int32),
        spline_nc=np.array([r["c"].size for r in results], dtype=np.int32),
        avg_spline_t=afin["t"], avg_spline_c=afin["c"], T_end=np.float64(T_END),
        onset_shift_ns=g("cum_shift"), t_pk=g("t_pk"), t50=g("t50"), iint_wa=g("iint_wa"),
        meta=np.array(json.dumps(meta)))
    out_npz = stem.with_suffix(".npz")
    _save_npz(out_npz, arr)
    write_root(str(stem.with_suffix(".root")), arr)
    write_pdf(str(stem.with_suffix(".pdf")), results, arr, run_type, str(stem.parent / "n_templates_map.png"), D)

    fb = arr["fallback"] == 0
    cf = D["cutflow"][[m.is_pwo4() for m in mods]].sum(0)
    own_s = fb & (arr["sigma_borrowed"] == 0)
    rr = np.arange(int(own_s.sum()))
    sg = arr["shape_sigma"][own_s].astype(np.float64)
    at = lambda t: np.rint((arr[t][own_s] - T_GRID[0]) / GRID_DT).astype(np.int64)   # noqa: E731
    el = fb & np.isfinite(arr["shape_sigma_elastic"]).any(1)
    ie = np.rint((arr["t50"][el] - T_GRID[0]) / GRID_DT).astype(np.int64)
    nu = arr["n_pool"][fb]
    print("[cutflow] (select, PbWO4) " + "  ".join(f"{n} {int(v)}" for n, v in zip(CUTFLOW_NAMES, cf)))
    print(f"[summary] own template {int(fb.sum())} (general pool {int(arr['pool_general'].sum())}), run-type average {int((~fb).sum() - arr['dead'].sum())} "
          f"({n_avg} pulses), dead {int(arr['dead'].sum())} ({' '.join(sorted(D['dead']))}); single pulse "
          f"{int(arr['single'].sum())}, shape_sigma borrowed {int(arr['sigma_borrowed'][fb].sum())}")
    print(f"[summary] pool size min {nu.min()} p10 {np.percentile(nu, 10):.0f} median {np.median(nu):.0f} "
          f"max {nu.max()}, >= {SIGMA_MIN_POOL}: {int((nu >= SIGMA_MIN_POOL).sum())}; total before cap median "
          f"{np.median(arr['n_pool_total'][fb]):.0f} max {arr['n_pool_total'][fb].max()}, capped {int((arr['n_pool'][fb] == POOL_MAX).sum())}; lambda "
          + " ".join(f"{x:g}:{int((arr['lam_module'] == x).sum())}" for x in LAM_GRID))
    print("[summary] alignment rms vs CFD (median over own modules with pool >= 2), per pass: "
          + " ".join(f"{np.nanmedian(rms[fb & (arr['n_pool'] >= 2), k]):.3f}" for k in range(ALIGN_PASSES)) + " ns")
    print(f"[summary] own shape_sigma median over modules: peak {100 * np.nanmedian(sg[rr, at('t_pk')]):.2f} %, "
          f"leading edge (t50) {100 * np.nanmedian(sg[rr, at('t50')]):.2f} % of the peak; elastic ({int(el.sum())} modules) "
          f"t50 {100 * np.nanmedian(arr['shape_sigma_elastic'][el][np.arange(int(el.sum())), ie]):.2f} %")
    print(f"[summary] pulse jitter median over own-sigma modules: rise {np.nanmedian(arr['rise_jitter_ns'][own_s]):.2f} ns, "
          f"fall {np.nanmedian(arr['fall_jitter_ns'][own_s]):.2f} ns")
    print(f"[done] {out_npz} ({out_npz.stat().st_size / 1e6:.1f} MB), .pdf, .root, n_templates_map.png; "
          f"{time.monotonic() - t_start:.0f} s", flush=True)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select", help="EVIO split files of one run -> pulse npz")
    s.add_argument("out", help="output .npz")
    s.add_argument("run_type", choices=list(RUN_TYPES))
    s.add_argument("evio", nargs="+", help="EVIO split files of one run, e.g. prad_025320.evio.00050 ...")
    b = sub.add_parser("build", help="OUT_ROOT/<RT>/pulses/run_*.npz -> templates npz / pdf / root")
    b.add_argument("run_type", choices=list(RUN_TYPES))
    a = ap.parse_args(argv)
    if a.cmd == "select":
        select_candidates(a.out, a.run_type, a.evio)
        return 0
    rt_dir = OUT_ROOT / a.run_type
    return build_templates(a.run_type, rt_dir / "pulses" / "run_*.npz", rt_dir / f"templates_{a.run_type}")


if __name__ == "__main__":
    sys.exit(main())
