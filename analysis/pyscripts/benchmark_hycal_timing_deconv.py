#!/usr/bin/env python3
"""benchmark_hycal_timing_deconv.py — 2×2 HyCal timing-gate × neighbour-selection
benchmark using the C++ template-fit (deconvolution) path for pulse amplitudes.

Phase 1 scope: template-fit amplitudes from WaveAnalyzer.deconvolve + CFD timing
from WaveResult.peaks[k].time; comparison across four clustering configurations:

  legacy_largestE  — seed_time_window ≤ 0, one largest-integral peak per module
  legacy_smallestDt — seed_time_window ≤ 0, one closest-to-seed-time peak per module
  gated_largestE   — seed_time_window = W ns, all in-window peaks; C++ BFS picks
                     largest-energy neighbour
  gated_smallestDt — seed_time_window = W ns, Python pre-filter to one
                     closest-to-seed-time peak per module (effectively smallest-Δt)

Headline question: does using template-fit amplitudes + CFD times + smallest-Δt
neighbour selection restore the expected σ_E/E(gated) ≤ σ_E/E(legacy) ordering?

Phase 1 non-goals
-----------------
* Matched-filter in the detection loop (Phase 2 — needs binding changes).
* Isolated CFD call (baked into WaveAnalyzer::findPeaks, not callable from Python).
* C++ BFS policy change to smallest-Δt (Phase 3 — deferred until Phase 1 data).
* No C++ changes, no pybind11 changes.

Dependencies
------------
Hard:   numpy, matplotlib, prad2py (dec, det modules)
Soft:   scipy.optimize.curve_fit  — Gaussian/CB fits; script runs without it
        (NaN fit parameters, stderr warning).

Output artifacts
----------------
  <out>_dt.tsv          — 14-column (seed, neighbour) timing rows (§8.1)
  <out>_spectrum.png    — four-config overlaid single-cluster energy spectra (§8.2)
  <out>_summary.csv     — four rows + #-prefixed preamble (§8.3)
  <out>_calibrated.png  — four-config overlay with CB fits (--calibrate only, §8.4)
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import _common as C
from prad2py import dec, det  # noqa: E402

try:
    from scipy.optimize import curve_fit
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CONFIG_NAMES = (
    "legacy_largestE",
    "legacy_smallestDt",
    "gated_largestE",
    "gated_smallestDt",
)

TSV_COLUMNS = [
    "event_num", "trigger_bits",
    "seed_module", "seed_id", "seed_time", "seed_energy",
    "neighbor_module", "neighbor_id", "neighbor_time", "neighbor_energy",
    "dt", "dx_q", "dy_q", "dist_q",
]

# Default template file path
_SCRIPT_DIR = Path(__file__).resolve().parent
_DEFAULT_TEMPLATE = _SCRIPT_DIR / "output" / \
    "pulse_templates_025308_h50_RawSum_max20modules.json"

# Colour map for four configs: cool tones for legacy, warm for gated
_CONFIG_COLOURS = {
    "legacy_largestE":  "#7f7f7f",   # grey
    "legacy_smallestDt": "#1f77b4",  # blue
    "gated_largestE":   "#ff7f0e",   # orange
    "gated_smallestDt": "#d62728",   # red
}

# ---------------------------------------------------------------------------
# Math helpers (mirrored from benchmark_hycal_timing.py)
# ---------------------------------------------------------------------------

def gauss(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def crystal_ball(x, A, mu, sigma, alpha, n):
    """Crystal-Ball PDF (left-side power-law tail).  α > 0, n > 1."""
    z = (np.asarray(x, dtype=float) - mu) / sigma
    abs_a = abs(alpha)
    A_coef = (n / abs_a) ** n * np.exp(-0.5 * abs_a ** 2)
    B_coef = n / abs_a - abs_a
    gauss_part = np.exp(-0.5 * z * z)
    arg = B_coef - z
    arg = np.where(arg > 0, arg, 1e-30)
    tail_part = A_coef * arg ** (-n)
    return A * np.where(z > -alpha, gauss_part, tail_part)


def _nan() -> float:
    return float("nan")


def _ratio(s: float, m: float) -> float:
    if m and m > 0 and not math.isnan(m) and not math.isnan(s):
        return s / m
    return _nan()

# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class PeakRecord:
    """One candidate peak from a single HyCal channel, post-pre-window filter."""
    module_index: int
    module_id:    int
    module_name:  str
    module_type:  str       # mod.type.name — "PbWO4" or "PbGlass"
    energy_MeV:   float     # from out.integral[k] via mod.energize (deconv success)
                            # or pk.integral via mod.energize (fallback)
    time_ns:      float     # wres.peaks[k].time — CFD-timed (baked into findPeaks)
    deconv_state: int       # out.state (one of dec.Q_DECONV_*)
    chi2_per_dof: float     # out.chi2_per_dof (NaN on fallback)
    used_fallback: bool     # True iff pk.integral was used instead of out.integral


@dataclass
class ConfigAccumulator:
    """Per-config state accumulated across all events."""
    name:       str
    hist:       np.ndarray      # int64 count array, length = N energy bins
    n_events:   int = 0         # events contributing exactly one cluster > signal_min
    # Lists populated only with --calibrate (memory cost is opt-in):
    energies:   List[float] = field(default_factory=list)
    cal_events: List[Tuple] = field(default_factory=list)
    # Fit results:
    fit_mu:        float = field(default=float("nan"))
    fit_sigma:     float = field(default=float("nan"))
    fit_mu_err:    float = field(default=float("nan"))
    fit_sigma_err: float = field(default=float("nan"))
    fit_n_in_window: int = 0
    # CB results (populated only with --calibrate):
    cb_mu:    float = field(default=float("nan"))
    cb_sigma: float = field(default=float("nan"))
    cb_alpha: float = field(default=float("nan"))
    cb_n:     float = field(default=float("nan"))
    # Post-calibration Gaussian:
    cal_mu:    float = field(default=float("nan"))
    cal_sigma: float = field(default=float("nan"))


@dataclass
class RunContext:
    """Immutable run configuration + mutable run-level tallies."""
    # Grid / cuts
    e_lo:   float
    e_hi:   float
    bin_width: float
    window: float           # seed_time_window for gated configs (ns)
    pre_window: Tuple[float, float]  # (lo, hi) ns
    require_trigger: int    # 0 = accept all
    max_fired_modules: int
    max_qdist: float
    signal_min: float
    # Calibration
    calibrate: bool
    target_peak: float
    inner_radius: float
    min_seed_energy: float
    min_per_module: int
    exclude_rowcol: Optional[Tuple[int, int, int, int]]
    cal_method: str
    # Meta
    scipy_available: bool
    template_file: Path
    # Derived
    e_edges: np.ndarray     # np.arange(e_lo, e_hi + bin_width, bin_width)
    # ---- Run-level mutable tallies (populated during event loop) ----
    n_events_total:              int   = 0
    n_events_fired_module_cut:   int   = 0
    n_peaks_deconv_per_channel:  int   = 0   # Q_DECONV_APPLIED
    n_peaks_deconv_global:       int   = 0   # Q_DECONV_FALLBACK_GLOBAL
    n_peaks_deconv_failed:       int   = 0   # anything else or n mismatch
    chi2_sum:                    float = 0.0
    chi2_count:                  int   = 0

# ---------------------------------------------------------------------------
# Template loading
# ---------------------------------------------------------------------------

def load_templates(path: Path, wave_cfg) -> Tuple[object, Dict[str, object]]:
    """Build PulseTemplateStore and pre-fetch per-type templates.

    Returns (store, dict) where dict maps type-name string ("PbWO4", "PbGlass")
    to a PulseTemplate object (or None if missing for that type).

    Missing-file / parse-error policy: load_from_file returns False on any
    failure; this function catches that, logs one WARN line, and returns
    ({}, {}).  Does NOT raise.  Downstream, store_by_type.get(mod.type.name)
    returns None for every channel → pk.integral fallback branch in extract_peaks.
    """
    store = dec.PulseTemplateStore()
    wcfg = (dec.WaveConfig(wave_cfg)
            if not isinstance(wave_cfg, dec.WaveConfig)
            else wave_cfg)
    if not store.load_from_file(str(path), wcfg):
        sys.stderr.write(
            f"[WARN] load_templates: PulseTemplateStore.load_from_file({path!r}) "
            f"returned False — template missing / parse error.  All channels "
            f"will use pk.integral fallback.\n"
        )
        return store, {}

    store_by_type: Dict[str, object] = {}
    for type_name in ("PbWO4", "PbGlass"):
        tmpl = store.type_template(type_name)
        if tmpl is None:
            sys.stderr.write(
                f"[WARN] load_templates: no template for type {type_name!r} in "
                f"{path} — channels of this type will use pk.integral fallback.\n"
            )
        store_by_type[type_name] = tmpl  # may be None; caller checks

    return store, store_by_type

# ---------------------------------------------------------------------------
# Per-channel peak extraction (deconv pipeline)
# ---------------------------------------------------------------------------

def extract_peaks(wa, store_by_type: Dict[str, object],
                  cd, mod, ctx: RunContext) -> List[PeakRecord]:
    """Run the deconv pipeline for one HyCal channel-data object.

    Pipeline:
      wa.analyze_result(samples) → wres (CFD timing, peak candidates)
      wa.deconvolve(samples, wres, tmpl) → out (template-fit amplitudes)

    Both Q_DECONV_APPLIED and Q_DECONV_FALLBACK_GLOBAL are treated as success.
    Deconv-health counters on ctx are updated for the entire channel (before
    the per-peak pre-window cut).

    Returns zero or more PeakRecord entries, one per peak passing pre-window.
    """
    wres = wa.analyze_result(cd.samples)
    if wres.npeaks == 0:
        return []

    tmpl = store_by_type.get(mod.type.name)   # None if template missing
    out = (wa.deconvolve(cd.samples, wres, tmpl)
           if tmpl is not None else None)

    # Determine deconv success
    used_deconv = (
        out is not None
        and out.state in (dec.Q_DECONV_APPLIED, dec.Q_DECONV_FALLBACK_GLOBAL)
        and out.n == wres.npeaks
    )

    # Update run-level deconv-health counters (before pre-window cut)
    if out is not None and out.n == wres.npeaks:
        if out.state == dec.Q_DECONV_APPLIED:
            ctx.n_peaks_deconv_per_channel += out.n
            ctx.chi2_sum   += out.chi2_per_dof * out.n
            ctx.chi2_count += out.n
        elif out.state == dec.Q_DECONV_FALLBACK_GLOBAL:
            ctx.n_peaks_deconv_global += out.n
            ctx.chi2_sum   += out.chi2_per_dof * out.n
            ctx.chi2_count += out.n
        else:
            ctx.n_peaks_deconv_failed += out.n
    else:
        # tmpl was None, or out.n != wres.npeaks — every raw peak is a failure
        ctx.n_peaks_deconv_failed += wres.npeaks

    records: List[PeakRecord] = []
    pre_lo, pre_hi = ctx.pre_window
    for k in range(wres.npeaks):
        pk = wres.peaks[k]
        if pk.time <= pre_lo or pk.time >= pre_hi:
            continue
        if used_deconv:
            adc = out.integral[k]
        else:
            adc = pk.integral   # fallback: same field used by sibling scripts
        records.append(PeakRecord(
            module_index  = mod.index,
            module_id     = mod.id,
            module_name   = mod.name,
            module_type   = mod.type.name,
            energy_MeV    = mod.energize(adc),
            time_ns       = float(pk.time),
            deconv_state  = (out.state if out is not None
                             else dec.Q_DECONV_NOT_RUN),
            chi2_per_dof  = (out.chi2_per_dof if used_deconv else _nan()),
            used_fallback = not used_deconv,
        ))
    return records

# ---------------------------------------------------------------------------
# Pulse selection helpers
# ---------------------------------------------------------------------------

def max_energy_per_module(peaks: List[PeakRecord]) -> List[PeakRecord]:
    """Keep the single largest-energy peak per module (legacy_largestE input)."""
    best: Dict[int, PeakRecord] = {}
    for pk in peaks:
        mi = pk.module_index
        if mi not in best or pk.energy_MeV > best[mi].energy_MeV:
            best[mi] = pk
    return list(best.values())


def pre_filter_smallest_dt(peaks: List[PeakRecord],
                            seed_time: float) -> List[PeakRecord]:
    """Group peaks by module_index; keep the single peak with minimum
    |time_ns - seed_time| per module.  Tie-break: larger energy wins.

    Used for legacy_smallestDt and gated_smallestDt input streams.

    Edge cases: zero peaks → empty list; one peak per module → unchanged.
    """
    best: Dict[int, PeakRecord] = {}
    for pk in peaks:
        mi = pk.module_index
        dt_abs = abs(pk.time_ns - seed_time)
        if mi not in best:
            best[mi] = pk
        else:
            prev = best[mi]
            prev_dt = abs(prev.time_ns - seed_time)
            if dt_abs < prev_dt or (dt_abs == prev_dt
                                    and pk.energy_MeV > prev.energy_MeV):
                best[mi] = pk
    return list(best.values())


def _global_seed_time(peaks: List[PeakRecord]) -> float:
    """Largest-amplitude peak across all modules, deterministic tie-break.

    Primary:   energy_MeV (descending)
    Secondary: -module_index (lower index preferred on exact tie)
    Tertiary:  -time_ns (earlier peak preferred on double tie)

    Mirrors HyCalCluster::FindSeeds' global energy maximum.
    """
    seed = max(peaks, key=lambda p: (p.energy_MeV,
                                      -p.module_index,
                                      -p.time_ns))
    return seed.time_ns

# ---------------------------------------------------------------------------
# Four-config clustering
# ---------------------------------------------------------------------------

def build_four_clusterers(p, ctx: RunContext):
    """Construct four HyCalCluster instances differing only in seed_time_window.

    Uses the HyCalClusterConfig → set_config(cfg) pattern (benchmark_hycal_timing.py:309-322).
    Returns (cl_legacy_bigE, cl_legacy_smDt, cl_gated_bigE, cl_gated_smDt).
    """
    base_cfg = p.hc_clusterer.get_config()

    def _make_cluster(seed_window: float):
        cfg = det.HyCalClusterConfig()
        for f in ("min_module_energy", "min_center_energy",
                  "min_cluster_energy", "min_cluster_size",
                  "corner_conn", "split_iter", "least_split",
                  "log_weight_thres"):
            setattr(cfg, f, getattr(base_cfg, f))
        cfg.seed_time_window = seed_window
        cl = det.HyCalCluster(p.hycal)
        cl.set_config(cfg)
        return cl

    cl_legacy_bigE = _make_cluster(-1.0)
    cl_legacy_smDt = _make_cluster(-1.0)
    cl_gated_bigE  = _make_cluster(ctx.window)
    cl_gated_smDt  = _make_cluster(ctx.window)
    return cl_legacy_bigE, cl_legacy_smDt, cl_gated_bigE, cl_gated_smDt


def run_configs(peaks: List[PeakRecord],
                cl_legacy_bigE, cl_legacy_smDt,
                cl_gated_bigE, cl_gated_smDt,
                ctx: RunContext) -> Optional[Dict[str, list]]:
    """Feed peaks into the four clusterers and return reconstructed clusters.

    Returns dict keyed by CONFIG_NAMES, each value a list of reconstructed
    cluster hits.  Returns None if peaks is empty.
    """
    if not peaks:
        return None

    seed_time = _global_seed_time(peaks)
    smDt_peaks = pre_filter_smallest_dt(peaks, seed_time)
    bigE_peaks  = max_energy_per_module(peaks)

    # Feed each clusterer the appropriate stream
    for cl, stream in [
        (cl_legacy_bigE, bigE_peaks),
        (cl_legacy_smDt, smDt_peaks),
        (cl_gated_bigE,  peaks),        # all surviving pulses — C++ BFS gates
        (cl_gated_smDt,  smDt_peaks),   # one per module → effectively smallest-Δt
    ]:
        cl.clear()
        for pk in stream:
            cl.add_hit(pk.module_index, pk.energy_MeV, pk.time_ns)
        cl.form_clusters()

    result = {}
    for name, cl in zip(CONFIG_NAMES, [cl_legacy_bigE, cl_legacy_smDt,
                                        cl_gated_bigE,  cl_gated_smDt]):
        result[name] = cl.reconstruct_matched()
    return result

# ---------------------------------------------------------------------------
# Spectrum accumulation
# ---------------------------------------------------------------------------

def accumulate_spectrum(clusters_by_cfg: Dict[str, list],
                        accs: Dict[str, ConfigAccumulator],
                        hycal,
                        ctx: RunContext) -> None:
    """For each config: if exactly one cluster passes --signal-min, accumulate."""
    for cfg_name, recs in clusters_by_cfg.items():
        big = [r for r in recs if r.hit.energy > ctx.signal_min]
        if len(big) != 1:
            continue
        r = big[0]
        E = r.hit.energy
        acc = accs[cfg_name]
        acc.n_events += 1
        # Histogram
        idx = int((E - ctx.e_lo) / ctx.bin_width)
        if 0 <= idx < len(acc.hist):
            acc.hist[idx] += 1
        # Calibration data (opt-in)
        if ctx.calibrate:
            acc.energies.append(E)
            sm = hycal.module_by_id(r.hit.center_id)
            acc.cal_events.append((
                r.hit.center_id,
                sm.x, sm.y,
                sm.row, sm.column,
                r.cluster.center.energy,
                E,
            ))

# ---------------------------------------------------------------------------
# Fit helpers (mirrors benchmark_hycal_timing.py)
# ---------------------------------------------------------------------------

def _fit_peak_gauss(values: List[float], e_lo: float, e_hi: float,
                    bin_width: float):
    """Gaussian fit around histogram argmax.  Returns (mu, sigma, A, n_in_fit)."""
    nan = _nan()
    edges  = np.arange(e_lo, e_hi + bin_width, bin_width)
    counts, _ = np.histogram(values, bins=edges)
    centres = 0.5 * (edges[:-1] + edges[1:])
    if not HAVE_SCIPY or counts.sum() == 0:
        return nan, nan, nan, 0, edges, counts

    mu0  = float(centres[int(np.argmax(counts))])
    win  = max(150.0, 4.0 * bin_width)
    mask = (centres > mu0 - win) & (centres < mu0 + win)
    x, y = centres[mask], counts[mask].astype(float)
    if y.sum() < 20:
        return nan, nan, nan, int(y.sum()), edges, counts
    try:
        popt, _ = curve_fit(
            gauss, x, y,
            p0     = [float(y.max()), mu0, max(bin_width, 60.0)],
            bounds = ([0, e_lo, 1.0], [np.inf, e_hi, win]))
        A, mu, sig = popt
        return float(mu), abs(float(sig)), float(A), int(y.sum()), edges, counts
    except Exception as exc:
        sys.stderr.write(f"[WARN] gauss fit failed ({exc})\n")
        return nan, nan, nan, int(y.sum()), edges, counts


def _fit_peak_cb(values: List[float], e_lo: float, e_hi: float,
                 bin_width: float, mu0_hint: Optional[float] = None):
    """Crystal-Ball fit in [e_lo, e_hi].  Returns (mu, sigma, alpha, n, A, n_in_fit)."""
    nan = _nan()
    edges  = np.arange(e_lo, e_hi + bin_width, bin_width)
    counts, _ = np.histogram(values, bins=edges)
    centres = 0.5 * (edges[:-1] + edges[1:])
    if not HAVE_SCIPY or counts.sum() == 0:
        return nan, nan, nan, nan, nan, 0

    imax = int(np.argmax(counts))
    mu0  = float(mu0_hint) if mu0_hint is not None else float(centres[imax])
    win  = max(400.0, 12.0 * bin_width)
    mask = (centres > mu0 - win) & (centres < mu0 + 0.5 * win)
    x, y = centres[mask], counts[mask].astype(float)
    if y.sum() < 50:
        return nan, nan, nan, nan, nan, int(y.sum())
    try:
        p0     = [float(y.max()), mu0, max(60.0, bin_width * 3), 1.5, 2.0]
        bounds = ([0,        e_lo - 200.0, 1.0,  0.1,  1.01],
                  [np.inf,   e_hi + 200.0, win,  20.0, 50.0])
        popt, _ = curve_fit(crystal_ball, x, y, p0=p0, bounds=bounds,
                            maxfev=10000)
        A, mu, sig, alpha, n = popt
        return float(mu), abs(float(sig)), float(alpha), float(n), float(A), int(y.sum())
    except Exception as exc:
        sys.stderr.write(f"[WARN] crystal-ball fit failed: {exc}\n")
        return nan, nan, nan, nan, nan, int(y.sum())


def _per_module_peak(values, bin_width: float) -> float:
    """Per-module peak position for gain refinement (mirrors benchmark_hycal_timing.py)."""
    vals = np.asarray(values, dtype=float)
    if len(vals) == 0:
        return _nan()
    med = float(np.median(vals))
    if len(vals) < 10 or not HAVE_SCIPY:
        return med

    e_lo_m, e_hi_m = med - 600.0, med + 600.0
    edges  = np.arange(e_lo_m, e_hi_m + bin_width, bin_width)
    counts, _ = np.histogram(vals, bins=edges)
    if counts.sum() < 10:
        return med
    centres = 0.5 * (edges[:-1] + edges[1:])
    mu0 = float(centres[int(np.argmax(counts))])

    win = 200.0
    mask = (centres > mu0 - win) & (centres < mu0 + 0.5 * win)
    x, y = centres[mask], counts[mask].astype(float)
    if y.sum() < 10:
        return mu0
    try:
        popt, _ = curve_fit(
            crystal_ball, x, y,
            p0     = [float(y.max()), mu0, 50.0, 1.5, 2.0],
            bounds = ([0, e_lo_m, 5.0, 0.1, 1.01],
                      [np.inf, e_hi_m, win, 20.0, 50.0]),
            maxfev = 5000)
        return float(popt[1])
    except Exception:
        pass
    try:
        popt, _ = curve_fit(
            gauss, x, y,
            p0     = [float(y.max()), mu0, 50.0],
            bounds = ([0, e_lo_m, 5.0], [np.inf, e_hi_m, win]),
            maxfev = 2000)
        return float(popt[1])
    except Exception:
        pass
    return mu0


def fit_gaussian_per_config(acc: ConfigAccumulator, ctx: RunContext) -> None:
    """Gaussian fit on acc.hist with --signal-min as a left cut-off.
    Populates acc.fit_mu, acc.fit_sigma, acc.fit_n_in_window.
    """
    if not acc.energies and acc.n_events == 0:
        return
    # Reconstruct values list from histogram for fit (or use energies if available)
    if acc.energies:
        values = acc.energies
    else:
        # Reconstruct values from histogram counts (approximate, used when
        # --calibrate is off and energies list is empty)
        centres = 0.5 * (ctx.e_edges[:-1] + ctx.e_edges[1:])
        values = []
        for c, cnt in zip(centres, acc.hist):
            values.extend([c] * int(cnt))

    mu, sig, A, n_in, edges, counts = _fit_peak_gauss(
        values, ctx.e_lo, ctx.e_hi, ctx.bin_width)
    acc.fit_mu           = mu
    acc.fit_sigma        = sig
    acc.fit_n_in_window  = n_in


def fit_crystal_ball_per_config(acc: ConfigAccumulator, ctx: RunContext,
                                 corrected_values: List[float]) -> None:
    """Crystal-Ball fit on corrected (post-calibration) values.
    Populates acc.cb_mu, acc.cb_sigma, acc.cb_alpha, acc.cb_n.
    Also refits Gaussian to get acc.cal_mu, acc.cal_sigma.
    """
    cal_lo = ctx.target_peak - 800.0
    cal_hi = ctx.target_peak + 800.0
    mu, sig, A, n_in, edges, counts = _fit_peak_gauss(
        corrected_values, cal_lo, cal_hi, ctx.bin_width)
    acc.cal_mu    = mu
    acc.cal_sigma = sig

    cb_mu, cb_sig, cb_alpha, cb_n, cb_A, _ = _fit_peak_cb(
        corrected_values, cal_lo, cal_hi, ctx.bin_width, mu0_hint=mu)
    acc.cb_mu    = cb_mu
    acc.cb_sigma = cb_sig
    acc.cb_alpha = cb_alpha
    acc.cb_n     = cb_n

# ---------------------------------------------------------------------------
# Per-module gain calibration (mirrors benchmark_hycal_timing.fit_per_module_gains)
# ---------------------------------------------------------------------------

def calibrate_per_module(cal_events: List[Tuple], ctx: RunContext
                          ) -> Tuple[Dict[int, float], List[float], List[float]]:
    """Per-seed-module gain refinement on the inner-ring high-statistics sample.

    cal_events: list of (seed_id, x, y, row, col, seed_energy_MeV, cluster_energy_MeV)

    Returns (gains, raw_energies, corrected_energies) where raw/corrected are
    parallel lists for events whose seed module had ≥ min_per_module entries.
    """
    inner_r2 = ctx.inner_radius * ctx.inner_radius
    excl = ctx.exclude_rowcol

    def keep(e):
        _sid, x, y, row, col, se, _ce = e
        if (x * x + y * y) >= inner_r2 or se < ctx.min_seed_energy:
            return False
        if excl is not None and (excl[0] <= row <= excl[1]
                                  and excl[2] <= col <= excl[3]):
            return False
        return True

    pool = [e for e in cal_events if keep(e)]
    if not pool:
        return {}, [], []

    seed_ids = sorted({e[0] for e in pool})
    gains: Dict[int, float] = {sid: 1.0 for sid in seed_ids}

    n_iter = 2
    for _ in range(n_iter):
        per_mod: Dict[int, List[float]] = {sid: [] for sid in seed_ids}
        for sid, *_, ce in pool:
            per_mod[sid].append(ce * gains[sid])
        for sid, vals in per_mod.items():
            if len(vals) < ctx.min_per_module:
                continue
            if ctx.cal_method == "cb-peak":
                ref = _per_module_peak(vals, ctx.bin_width)
            else:
                ref = float(np.median(vals))
            if ref and ref > 0 and not math.isnan(ref):
                gains[sid] *= ctx.target_peak / ref

    counts_m = {sid: 0 for sid in seed_ids}
    for sid, *_ in pool:
        counts_m[sid] += 1
    converged = {sid: g for sid, g in gains.items()
                 if counts_m[sid] >= ctx.min_per_module}

    raw  = [e[6] for e in pool if e[0] in converged]
    corr = [e[6] * converged[e[0]] for e in pool if e[0] in converged]
    return converged, raw, corr

# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def write_dt_tsv(all_rows: list, path: Path) -> None:
    """Write the 14-column (seed, neighbour) timing TSV (§8.1)."""
    with open(path, "w", newline="", encoding="utf-8") as fh:
        fh.write("\t".join(TSV_COLUMNS) + "\n")
        for row in all_rows:
            fh.write("\t".join(str(v) for v in row) + "\n")


def _sha1(path: Path) -> str:
    try:
        h = hashlib.sha1()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(65536), b""):
                h.update(chunk)
        return h.hexdigest()[:12]
    except Exception:
        return "unknown"


def write_summary_csv(accs: Dict[str, ConfigAccumulator],
                      ctx: RunContext,
                      evio_path: str,
                      path: Path,
                      run_num: int = -1) -> None:
    """Write per-config summary rows + #-prefixed preamble (§8.3)."""
    total_deconv = (ctx.n_peaks_deconv_per_channel
                    + ctx.n_peaks_deconv_global
                    + ctx.n_peaks_deconv_failed)
    failed_frac = (ctx.n_peaks_deconv_failed / total_deconv
                   if total_deconv > 0 else 0.0)
    mean_chi2 = (ctx.chi2_sum / ctx.chi2_count
                 if ctx.chi2_count > 0 else _nan())
    fired_frac = (ctx.n_events_fired_module_cut / max(ctx.n_events_total, 1))

    tmpl_sha = _sha1(ctx.template_file)
    ts = time.strftime("%Y-%m-%dT%H:%M:%S")

    preamble_lines = [
        f"# benchmark_hycal_timing_deconv — Phase 1 summary",
        f"# generated={ts}",
        f"# evio_path={evio_path}",
        f"# run_num={run_num}",
        f"# template_file={ctx.template_file}  sha1={tmpl_sha}",
        f"# window={ctx.window}",
        f"# pre_window={ctx.pre_window}",
        f"# require_trigger={ctx.require_trigger:#x}",
        f"# signal_min={ctx.signal_min}",
        f"# e_range=[{ctx.e_lo},{ctx.e_hi}]",
        f"# bin_width={ctx.bin_width}",
        f"# max_fired_modules={ctx.max_fired_modules}",
        f"# max_qdist={ctx.max_qdist}",
        f"# calibrate={ctx.calibrate}",
        f"# cal_method={ctx.cal_method}",
        f"# target_peak={ctx.target_peak}",
        f"# scipy_available={ctx.scipy_available}",
        # Deconv health (§14.3)
        f"# n_peaks_deconv_per_channel={ctx.n_peaks_deconv_per_channel}",
        f"# n_peaks_deconv_global={ctx.n_peaks_deconv_global}",
        f"# n_peaks_deconv_failed={ctx.n_peaks_deconv_failed}",
        f"# deconv_failed_fraction={failed_frac:.2f}",
        f"# mean_chi2_per_dof={mean_chi2:.4f}" if not math.isnan(mean_chi2)
            else "# mean_chi2_per_dof=nan",
        # Event counts (§13)
        f"# n_events_total={ctx.n_events_total}",
        f"# n_events_fired_module_cut={ctx.n_events_fired_module_cut}",
        f"# fired_module_cut_fraction={fired_frac:.2f}",
    ]

    columns = [
        "config_name", "n_events", "n_in_fit_window",
        "mu_MeV", "sigma_MeV", "sigma_over_E_pct",
        "mu_calibrated", "sigma_calibrated",
        "sigma_CB", "sigma_over_E_CB",
    ]

    with open(path, "w", newline="", encoding="utf-8") as fh:
        for line in preamble_lines:
            fh.write(line + "\n")
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(columns)
        for cfg_name in CONFIG_NAMES:
            acc = accs[cfg_name]
            mu   = acc.fit_mu
            sig  = acc.fit_sigma
            soe  = (_ratio(sig, mu) * 100.0
                    if not (math.isnan(sig) or math.isnan(mu)) else _nan())
            cal_mu  = acc.cal_mu
            cal_sig = acc.cal_sigma
            cb_sig  = acc.cb_sigma
            cb_soe  = (_ratio(cb_sig, acc.cb_mu) * 100.0
                       if not (math.isnan(cb_sig) or math.isnan(acc.cb_mu))
                       else _nan())
            w.writerow([
                cfg_name,
                acc.n_events,
                acc.fit_n_in_window,
                f"{mu:.3f}" if not math.isnan(mu) else "nan",
                f"{sig:.3f}" if not math.isnan(sig) else "nan",
                f"{soe:.4f}" if not math.isnan(soe) else "nan",
                f"{cal_mu:.3f}" if not math.isnan(cal_mu) else "nan",
                f"{cal_sig:.3f}" if not math.isnan(cal_sig) else "nan",
                f"{cb_sig:.3f}" if not math.isnan(cb_sig) else "nan",
                f"{cb_soe:.4f}" if not math.isnan(cb_soe) else "nan",
            ])


def plot_spectrum(accs: Dict[str, ConfigAccumulator],
                  ctx: RunContext, out_path: Path) -> None:
    """Four overlaid histograms (one per config) with Gaussian fits (§8.2)."""
    centres = 0.5 * (ctx.e_edges[:-1] + ctx.e_edges[1:])
    fig, ax = plt.subplots(figsize=(10, 6))

    for cfg_name in CONFIG_NAMES:
        acc  = accs[cfg_name]
        col  = _CONFIG_COLOURS[cfg_name]
        hist = acc.hist.astype(float)
        ax.step(centres, hist, where="mid", color=col, lw=1.5,
                label=f"{cfg_name}  N={acc.n_events}")
        mu, sig = acc.fit_mu, acc.fit_sigma
        if HAVE_SCIPY and not math.isnan(mu) and not math.isnan(sig) and sig > 0:
            x_fit = np.linspace(mu - 4 * sig, mu + 4 * sig, 200)
            # Reconstruct A from histogram to scale the curve
            idx = int((mu - ctx.e_lo) / ctx.bin_width)
            A_est = float(hist[idx]) if 0 <= idx < len(hist) else float(hist.max())
            soe_pct = _ratio(sig, mu) * 100.0
            ax.plot(x_fit, gauss(x_fit, A_est, mu, sig), color=col,
                    ls="--", lw=1.0,
                    label=f"   μ={mu:.0f}  σ={sig:.0f} MeV  ({soe_pct:.2f}%)")

    ax.axvline(ctx.signal_min, color="#444", lw=0.8, ls=":",
               label=f"cut E > {ctx.signal_min:.0f} MeV")
    ax.set_xlabel("cluster energy (MeV)")
    ax.set_ylabel("events / bin")
    ax.set_title("Single-cluster energy — deconv template-fit amplitudes\n"
                 f"pre-window [{ctx.pre_window[0]:.0f}, {ctx.pre_window[1]:.0f}] ns  "
                 f"gate W={ctx.window:.1f} ns")
    ax.set_xlim(ctx.e_lo, ctx.e_hi)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def plot_calibrated(accs: Dict[str, ConfigAccumulator],
                    ctx: RunContext, out_path: Path,
                    corr_data: Dict[str, Tuple[list, list]]) -> None:
    """Four-config overlay with post-calibration CB fits (§8.4, --calibrate only)."""
    cal_lo = ctx.target_peak - 800.0
    cal_hi = ctx.target_peak + 800.0
    edges  = np.arange(cal_lo, cal_hi + ctx.bin_width, ctx.bin_width)
    centres_cal = 0.5 * (edges[:-1] + edges[1:])

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharey=False)
    ax_map = dict(zip(CONFIG_NAMES, axes.flat))

    for cfg_name in CONFIG_NAMES:
        acc = accs[cfg_name]
        col = _CONFIG_COLOURS[cfg_name]
        ax  = ax_map[cfg_name]
        raw, corr = corr_data.get(cfg_name, ([], []))

        if raw:
            cnt_raw, _ = np.histogram(raw, bins=edges)
            ax.step(centres_cal, cnt_raw, where="mid", color=col,
                    alpha=0.4, lw=1.3, label="pre-cal")
        if corr:
            cnt_corr, _ = np.histogram(corr, bins=edges)
            ax.step(centres_cal, cnt_corr, where="mid", color=col,
                    lw=1.8, label=f"post-cal  μ={acc.cal_mu:.0f} σ={acc.cal_sigma:.0f}")
            # CB overlay
            cb_mu, cb_sig, cb_alpha, cb_n = (acc.cb_mu, acc.cb_sigma,
                                              acc.cb_alpha, acc.cb_n)
            if not math.isnan(cb_mu):
                xfine = np.linspace(cb_mu - 6 * cb_sig, cb_mu + 4 * cb_sig, 400)
                A_est = float(cnt_corr.max()) if len(cnt_corr) > 0 else 1.0
                ax.plot(xfine,
                        crystal_ball(xfine, A_est, cb_mu, cb_sig, cb_alpha, cb_n),
                        color="#d62728", lw=1.4,
                        label=f"CB  μ={cb_mu:.0f} σ={cb_sig:.0f}")

        ax.axvline(ctx.target_peak, color="#444", ls=":", lw=0.8)
        soe_cb = _ratio(acc.cb_sigma, acc.cb_mu) * 100.0
        ax.set_title(f"{cfg_name}\nσ_CB/E={soe_cb:.2f}%")
        ax.set_xlabel("cluster energy (MeV)")
        ax.set_ylabel("events / bin")
        ax.grid(alpha=0.3)
        ax.legend(loc="upper right", fontsize=8)

    fig.suptitle(
        f"Per-module gain refinement  (r < {ctx.inner_radius:.0f} mm, "
        f"seed E ≥ {ctx.min_seed_energy:.0f} MeV, "
        f"≥ {ctx.min_per_module} ev/mod,  target={ctx.target_peak:.1f} MeV)",
        fontsize=10.5)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)

# ---------------------------------------------------------------------------
# Argument validation
# ---------------------------------------------------------------------------

def _validate_args(ns: argparse.Namespace) -> None:
    errs = []
    if ns.window <= 0.0:
        errs.append(f"--window must be > 0 (got {ns.window})")
    if ns.bin_width <= 0.0:
        errs.append(f"--bin-width must be > 0 (got {ns.bin_width})")
    lo, hi = ns.pre_window
    if lo >= hi:
        errs.append(f"--pre-window LO must be < HI (got {lo} >= {hi})")
    e_lo, e_hi = ns.e_range
    if e_lo >= e_hi:
        errs.append(f"--e-range LO must be < HI (got {e_lo} >= {e_hi})")
    if errs:
        sys.stderr.write("[ERROR] argument validation:\n")
        for e in errs:
            sys.stderr.write(f"  {e}\n")
        sys.exit(2)

# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=(
            "benchmark_hycal_timing_deconv — 2×2 HyCal timing-gate × "
            "neighbour-selection benchmark using WaveAnalyzer.deconvolve "
            "template-fit amplitudes and CFD timing."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Output artifacts (all share the OUT_PATH stem):
  <stem>_dt.tsv           14-column (seed, neighbour) timing rows
  <stem>_spectrum.png     four-config overlaid single-cluster energy spectra
  <stem>_summary.csv      four rows + #-prefixed preamble (deconv health + counters)
  <stem>_calibrated.png   four-config CB-fit overlay (--calibrate only)
""")
    C.add_common_args(ap)

    # Input / pipeline
    ap.add_argument("--template-file", type=Path,
                    default=_DEFAULT_TEMPLATE,
                    metavar="PATH",
                    help=f"Per-type pulse-template JSON consumed by "
                         f"PulseTemplateStore.load_from_file. "
                         f"Default: {_DEFAULT_TEMPLATE}")

    # Timing-gate configs
    ap.add_argument("--window", type=float, default=4.0,
                    metavar="W",
                    help="seed_time_window in ns for the two gated configs "
                         "(default 4.0).")

    # Pulse pre-filter
    ap.add_argument("--pre-window", type=float, nargs=2,
                    default=(100.0, 200.0), metavar=("LO", "HI"),
                    help="Keep peaks with LO < pk.time < HI in ns "
                         "(default 100 200).")

    # Event-level cuts (script-local; NOT in _common.add_common_args)
    ap.add_argument("--require-trigger", type=lambda s: int(s, 0),
                    default=0, metavar="MASK",
                    help="Trigger-bits gate (hex 0x... or decimal). "
                         "Default 0 = accept any physics trigger (unbiased). "
                         "Pass e.g. 0x100 for SSP0 strict equality.")
    ap.add_argument("--max-fired-modules", type=int, default=50,
                    metavar="N",
                    help="Skip events where more than N HyCal channels produced "
                         "≥1 peak passing --pre-window (default 50).")

    # CollectNeighborTiming
    ap.add_argument("--max-qdist", type=float, default=5.0,
                    metavar="Q",
                    help="Max quantized module-distance for CollectNeighborTiming "
                         "(default 5.0).")

    # Spectrum selection
    ap.add_argument("--signal-min", type=float, default=500.0,
                    metavar="E_MeV",
                    help="One-cluster energy cut in MeV (default 500, tuned for "
                         "2.2 GeV e-p).")
    ap.add_argument("--e-range", type=float, nargs=2,
                    default=(500.0, 2500.0), metavar=("LO", "HI"),
                    help="Energy histogram range in MeV (default 500 2500).")
    ap.add_argument("--bin-width", type=float, default=20.0,
                    metavar="BW",
                    help="Histogram bin width in MeV (default 20).")

    # Calibration
    ap.add_argument("--calibrate", action="store_true",
                    help="Per-seed-module gain refinement + Crystal-Ball re-fit.")
    ap.add_argument("--target-peak", type=float, default=2239.0,
                    metavar="E_MeV",
                    help="Calibration target peak position in MeV "
                         "(default 2239.0, e-p elastic at 2.2 GeV).")
    ap.add_argument("--inner-radius", type=float, default=200.0,
                    metavar="R_mm",
                    help="Radius cut for per-module gain fit in mm (default 200).")
    ap.add_argument("--min-seed-energy", type=float, default=200.0,
                    metavar="E_MeV",
                    help="Seed-energy cut for per-module gain fit in MeV "
                         "(default 200).")
    ap.add_argument("--min-per-module", type=int, default=30,
                    metavar="N",
                    help="Minimum events per module for per-module gain fit "
                         "(default 30).")
    ap.add_argument("--exclude-rowcol", type=int, nargs=4, default=None,
                    metavar=("R_LO", "R_HI", "C_LO", "C_HI"),
                    help="Exclude seed modules with row in [R_LO, R_HI] AND "
                         "col in [C_LO, C_HI] from per-module calibration sample.")
    ap.add_argument("--cal-method", choices=["median", "cb-peak"],
                    default="median",
                    help="Per-module reference for gain refinement: 'median' or "
                         "'cb-peak' (Crystal-Ball fitted peak; default: median).")

    args = ap.parse_args(argv)
    if not HAVE_SCIPY:
        print("[WARN] scipy not available — Gaussian/CB fits will be NaN", file=sys.stderr)
    _validate_args(args)

    # Build RunContext
    e_lo, e_hi = args.e_range
    bw = args.bin_width
    e_edges = np.arange(e_lo, e_hi + bw, bw)
    n_bins  = len(e_edges) - 1

    ctx = RunContext(
        e_lo               = e_lo,
        e_hi               = e_hi,
        bin_width          = bw,
        window             = args.window,
        pre_window         = tuple(args.pre_window),
        require_trigger    = args.require_trigger,
        max_fired_modules  = args.max_fired_modules,
        max_qdist          = args.max_qdist,
        signal_min         = args.signal_min,
        calibrate          = args.calibrate,
        target_peak        = args.target_peak,
        inner_radius       = args.inner_radius,
        min_seed_energy    = args.min_seed_energy,
        min_per_module     = args.min_per_module,
        exclude_rowcol     = (tuple(args.exclude_rowcol)
                              if args.exclude_rowcol else None),
        cal_method         = args.cal_method,
        scipy_available    = HAVE_SCIPY,
        template_file      = args.template_file,
        e_edges            = e_edges,
    )

    # Setup pipeline
    p = C.setup_pipeline(
        evio_path     = args.evio_path,
        max_events    = args.max_events,
        run_num       = args.run_num,
        gem_ped_file  = args.gem_ped_file,
        gem_cm_file   = args.gem_cm_file,
        hc_calib_file = args.hc_calib_file,
        daq_config    = args.daq_config,
        gem_map_file  = args.gem_map_file,
        hc_map_file   = args.hc_map_file,
    )
    print(f"[setup] template   : {args.template_file}", flush=True)
    print(f"[setup] pre-window : [{ctx.pre_window[0]}, {ctx.pre_window[1]}] ns",
          flush=True)
    print(f"[setup] gate window: ±{ctx.window} ns (gated configs)", flush=True)
    print(f"[setup] signal cut : single cluster E > {ctx.signal_min:.0f} MeV",
          flush=True)
    print(f"[setup] trigger    : {ctx.require_trigger:#x} "
          f"({'unbiased' if ctx.require_trigger == 0 else 'gated'})", flush=True)
    print(f"[setup] max-fired  : {ctx.max_fired_modules} modules", flush=True)

    # Load templates
    _store, store_by_type = load_templates(args.template_file, p.cfg.wave_cfg)

    # Build four clusterers
    cl_legacy_bigE, cl_legacy_smDt, cl_gated_bigE, cl_gated_smDt = \
        build_four_clusterers(p, ctx)

    # Build a scratch clusterer for CollectNeighborTiming
    base_cfg_cnt = p.hc_clusterer.get_config()
    cl_cnt = det.HyCalCluster(p.hycal)
    cl_cnt.set_config(base_cfg_cnt)

    # Accumulators
    accs: Dict[str, ConfigAccumulator] = {
        name: ConfigAccumulator(name=name, hist=np.zeros(n_bins, dtype=np.int64))
        for name in CONFIG_NAMES
    }

    # Event loop
    dt_rows: list = []
    ch = dec.EvChannel()
    ch.set_config(p.cfg)

    t0 = time.monotonic()
    n_read = n_phys = n_kept = 0
    n_files_open = 0
    unknown_type_warned: set = set()

    try:
        for fpath in p.evio_files:
            if ch.open_auto(fpath) != dec.Status.success:
                print(f"[WARN] skip (cannot open): {fpath}", flush=True)
                continue
            n_files_open += 1
            print(f"[file {n_files_open}/{len(p.evio_files)}] {fpath}",
                  flush=True)

            done = False
            while ch.read() == dec.Status.success:
                n_read += 1
                if not ch.scan():
                    continue
                if ch.get_event_type() != dec.EventType.Physics:
                    continue

                for i in range(ch.get_n_events()):
                    decoded = ch.decode_event(i, with_ssp=False)
                    if not decoded["ok"]:
                        continue
                    n_phys += 1
                    fadc_evt = decoded["event"]

                    # --- Step 1: trigger filter ---
                    trigger_bits = int(fadc_evt.info.trigger_bits)
                    if ctx.require_trigger and trigger_bits != ctx.require_trigger:
                        if args.max_events > 0 and n_phys >= args.max_events:
                            done = True; break
                        continue
                    n_kept += 1
                    ctx.n_events_total += 1
                    event_num = int(fadc_evt.info.event_number)

                    # --- Step 2: collect peaks per channel ---
                    peaks: List[PeakRecord] = []
                    fired_module_set: set = set()

                    for ri in range(fadc_evt.nrocs):
                        roc = fadc_evt.roc(ri)
                        if not roc.present:
                            continue
                        crate = p.crate_map.get(roc.tag)
                        if crate is None:
                            continue
                        for s in roc.present_slots():
                            slot = roc.slot(s)
                            for c in slot.present_channels():
                                mod = p.hycal.module_by_daq(crate, s, c)
                                if mod is None or not mod.is_hycal():
                                    continue
                                cd = slot.channel(c)
                                if cd.nsamples <= 0:
                                    continue

                                # Warn once per unknown module type
                                tname = mod.type.name
                                if (tname not in ("PbWO4", "PbGlass")
                                        and tname not in unknown_type_warned):
                                    sys.stderr.write(
                                        f"[WARN] unexpected module type {tname!r} "
                                        f"on module {mod.name} — using pk.integral "
                                        f"fallback\n")
                                    unknown_type_warned.add(tname)

                                ch_peaks = extract_peaks(
                                    p.wave_ana, store_by_type, cd, mod, ctx)
                                for pk in ch_peaks:
                                    peaks.append(pk)
                                    fired_module_set.add(mod.index)

                    # --- Step 3: fired-module count cut ---
                    if len(fired_module_set) > ctx.max_fired_modules:
                        ctx.n_events_fired_module_cut += 1
                        if args.max_events > 0 and n_phys >= args.max_events:
                            done = True; break
                        continue

                    # --- Step 4: CollectNeighborTiming on raw peak stream ---
                    if peaks:
                        cl_cnt.clear()
                        for pk in peaks:
                            cl_cnt.add_hit(pk.module_index, pk.energy_MeV,
                                           pk.time_ns)
                        cnt_rows = cl_cnt.collect_neighbor_timing(ctx.max_qdist)
                        for r in cnt_rows:
                            seed_mod = p.hycal.module(r.seed_module)
                            nbr_mod  = p.hycal.module(r.neighbor_module)
                            dist_q   = math.sqrt(r.dx_q * r.dx_q + r.dy_q * r.dy_q)
                            dt_rows.append([
                                event_num, trigger_bits,
                                r.seed_module, seed_mod.id,
                                f"{r.seed_time:.3f}", f"{r.seed_energy:.4f}",
                                r.neighbor_module, nbr_mod.id,
                                f"{r.neighbor_time:.3f}", f"{r.neighbor_energy:.4f}",
                                f"{r.dt:.3f}",
                                f"{r.dx_q:.4f}", f"{r.dy_q:.4f}", f"{dist_q:.4f}",
                            ])

                    # --- Step 5: build and fill the four configs ---
                    if peaks:
                        clusters_by_cfg = run_configs(
                            peaks,
                            cl_legacy_bigE, cl_legacy_smDt,
                            cl_gated_bigE,  cl_gated_smDt,
                            ctx)
                        if clusters_by_cfg is not None:
                            accumulate_spectrum(clusters_by_cfg, accs,
                                                p.hycal, ctx)

                    if args.max_events > 0 and n_phys >= args.max_events:
                        done = True; break

                if done:
                    break

                # Progress line every 5000 physics events
                if n_phys > 0 and n_phys % 5000 == 0:
                    elapsed_p = time.monotonic() - t0
                    rate = n_phys / max(elapsed_p, 1e-3)
                    evts = {name: accs[name].n_events for name in CONFIG_NAMES}
                    print(f"[progress] {n_phys} physics  kept={n_kept}  "
                          f"dt_rows={len(dt_rows)}  "
                          f"legBigE={evts['legacy_largestE']}  "
                          f"legSmDt={evts['legacy_smallestDt']}  "
                          f"gatBigE={evts['gated_largestE']}  "
                          f"gatSmDt={evts['gated_smallestDt']}  "
                          f"({rate:.1f} ev/s)", flush=True)

            ch.close()
            if done:
                break

    finally:
        elapsed = time.monotonic() - t0

    # ---- Post-loop: Gaussian fits per config --------------------------------
    for cfg_name in CONFIG_NAMES:
        fit_gaussian_per_config(accs[cfg_name], ctx)

    # ---- Optional calibration -----------------------------------------------
    corr_data: Dict[str, Tuple[list, list]] = {}
    if args.calibrate:
        print("\n[calibrate] running per-module gain refinement ...", flush=True)
        for cfg_name in CONFIG_NAMES:
            acc = accs[cfg_name]
            if not acc.cal_events:
                print(f"  {cfg_name}: no cal events — skipping", flush=True)
                corr_data[cfg_name] = ([], [])
                continue
            _gains, raw, corr = calibrate_per_module(acc.cal_events, ctx)
            corr_data[cfg_name] = (raw, corr)
            if corr:
                fit_crystal_ball_per_config(acc, ctx, corr)
                soe_cb = _ratio(acc.cb_sigma, acc.cb_mu) * 100.0
                print(f"  {cfg_name}: {len(corr)} inner events  "
                      f"μ={acc.cal_mu:.1f}  σ_CB={acc.cb_sigma:.1f}  "
                      f"σ_CB/E={soe_cb:.2f}%", flush=True)

    # ---- Write output artifacts (stable order) --------------------------------
    stem = args.out_path
    dt_tsv_path     = Path(f"{stem}_dt.tsv")
    spectrum_path   = Path(f"{stem}_spectrum.png")
    summary_path    = Path(f"{stem}_summary.csv")
    calibrated_path = Path(f"{stem}_calibrated.png")

    # Ensure output directory exists
    dt_tsv_path.parent.mkdir(parents=True, exist_ok=True)

    write_dt_tsv(dt_rows, dt_tsv_path)
    print(f"\n  wrote: {dt_tsv_path}")

    plot_spectrum(accs, ctx, spectrum_path)
    print(f"  wrote: {spectrum_path}")

    write_summary_csv(accs, ctx, args.evio_path, summary_path, run_num=args.run_num)
    print(f"  wrote: {summary_path}")

    if args.calibrate:
        plot_calibrated(accs, ctx, calibrated_path, corr_data)
        print(f"  wrote: {calibrated_path}")

    # ---- Deconv health summary to stderr -------------------------------------
    total_deconv = (ctx.n_peaks_deconv_per_channel
                    + ctx.n_peaks_deconv_global
                    + ctx.n_peaks_deconv_failed)
    failed_frac = (ctx.n_peaks_deconv_failed / total_deconv
                   if total_deconv > 0 else 0.0)
    mean_chi2 = (ctx.chi2_sum / ctx.chi2_count
                 if ctx.chi2_count > 0 else _nan())
    chi2_str = f"{mean_chi2:.4f}" if not math.isnan(mean_chi2) else "nan"
    sys.stderr.write(
        f"[deconv health] per_channel={ctx.n_peaks_deconv_per_channel}  "
        f"global={ctx.n_peaks_deconv_global}  "
        f"failed={ctx.n_peaks_deconv_failed}\n"
        f"                mean_chi2/dof={chi2_str}  "
        f"failed_frac={failed_frac * 100:.1f}%\n"
    )

    # ---- One-line summary ---------------------------------------------------
    rate_final = n_phys / max(elapsed, 1e-3)
    print()
    print("=" * 72)
    print("benchmark_hycal_timing_deconv — Phase 1 done")
    print("=" * 72)
    print(f"  files {n_files_open}/{len(p.evio_files)},  "
          f"records {n_read},  physics {n_phys},  triggered {n_kept}")
    print(f"  fired-module cut: {ctx.n_events_fired_module_cut} / "
          f"{ctx.n_events_total} events skipped")
    print(f"  elapsed {elapsed:.1f} s  ({rate_final:.1f} ev/s)")
    print()
    print(f"  {'config':<20} {'N_events':>10}  {'μ (MeV)':>10}  "
          f"{'σ (MeV)':>10}  {'σ/E (%)':>10}")
    print(f"  {'-'*64}")
    for cfg_name in CONFIG_NAMES:
        acc = accs[cfg_name]
        mu_s  = f"{acc.fit_mu:.1f}" if not math.isnan(acc.fit_mu) else "nan"
        sig_s = f"{acc.fit_sigma:.1f}" if not math.isnan(acc.fit_sigma) else "nan"
        soe   = _ratio(acc.fit_sigma, acc.fit_mu) * 100.0
        soe_s = f"{soe:.2f}" if not math.isnan(soe) else "nan"
        print(f"  {cfg_name:<20} {acc.n_events:>10}  {mu_s:>10}  "
              f"{sig_s:>10}  {soe_s:>10}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
