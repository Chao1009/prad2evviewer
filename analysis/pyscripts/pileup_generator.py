#!/usr/bin/env python3
"""
pileup_generator.py — deterministic (ΔT, ratio) pile-up response-surface
generator for PbWO4 PbWO4 HyCal modules (PRad-II).

v0.1 scope
----------
Produces a bag of synthetic 2-pulse waveforms on a (ΔT, ratio) grid,
where every injected pulse has **known truth** (t0, height, template shape).
The downstream user runs dec.WaveAnalyzer.deconvolve() on each synthetic
event and compares recovered heights and onsets to truth to produce the
"amplitude error vs ΔT vs ratio" response-surface plots.

This script does NOT run any deconvolution, does NOT modify C++ code, and
does NOT produce ROC curves.  See §13 of the plan for the full non-goals list.

Time convention (§1a.1)
-----------------------
The two-τ template is parameterised by an **onset** t0:

    T(t; t0, τ_r, τ_f) = (1 − exp(−(t−t0)/τ_r)) · exp(−(t−t0)/τ_f) / T_max
                          for t > t0, else 0.

The template **peak** sits at a fixed positive offset past the onset:

    t_peak_offset = τ_r · log((τ_r + τ_f) / τ_r)  ≈ 5.5 ns (PbWO4)

WaveAnalyzer.analyze returns pk.time = the **peak position** (interpolated,
in ns), NOT the onset.  WaveAnalyzer.deconvolve reports dec_out.t0_ns =
the fitted **template onset**.  These differ by t_peak_offset.

ΔT in this generator = **peak-to-peak separation** (intuitive, scope-like).
Because both pulses share the same (τ_r, τ_f), the onset-to-onset
separation equals the peak-to-peak separation.

Each synthetic event stores BOTH:
  truth_peak_times_ns[:, k]  — peak time (ns);   compare to pk.time from analyze.
  truth_onset_t0_ns[:, k]   — template onset (ns); compare to dec_out.t0_ns.

Downstream rule: compare dec_out.t0_ns to truth_onset_t0_ns[:, 0],
                          dec_out.height  to truth_heights_adc[:, 0].
DO NOT compare dec_out.t0_ns to truth_peak_times_ns — they differ by
t_peak_offset (~5.5 ns for PbWO4), producing a spurious timing bias.

Amplitude convention (§1a.2)
-----------------------------
The injected pulse uses the **peak-normalised** template: inject_wave =
h_inject · T(t; t0, τ_r, τ_f) where T includes /T_max normalisation, so
h_inject IS the physical peak height in ADC.

Truth storage:
  truth_heights_adc[:, k]       — physical peak height (ADC).  PRIMARY truth.
                                   Compare to dec_out.height.
  truth_model_amplitudes[:, k]  — truth_heights_adc / T_max.  Convenience;
                                   compare to dec_out.amplitude.

Column convention: col 0 = injected (exact truth), col 1 = base_ref
(WaveAnalyzer/fit_pulse_shape observables — accurate but not exact truth).

Output schema (§3.2 — per-cell .npz)
--------------------------------------
  waveforms               uint16  (N, n_samples)  synthetic ADC, ready for analyze
  truth_heights_adc       float32 (N, 2)          [injected, base_ref] peak height ADC
  truth_model_amplitudes  float32 (N, 2)          heights / T_max
  truth_onset_t0_ns       float32 (N, 2)          template onset (ns)
  truth_peak_times_ns     float32 (N, 2)          peak position (ns)
  event_dt_ns             float32 (N,)            per-event ΔT used for injection (ns)
                                                   fixed-grid: all equal meta_dt_ns
                                                   random-dt:  varies per event
  truth_source            <U9     (2,)             ["injected", "base_ref"]
  pedestals               float32 (N, 2)          [ped_mean, ped_rms]
  base_channel_names      <U8     (N,)
  base_channel_ids        <U16    (N,)
  base_roc_tags           int32   (N,)
  base_slots              int32   (N,)
  base_chans              int32   (N,)
  base_run_numbers        int32   (N,)
  base_source_files       <U64    (N,)
  base_physics_event_index int32  (N,)
  base_chi2_per_dof       float32 (N,)
  meta_dt_ns              float32 ()               fixed-grid ΔT; NaN in random-dt mode
  meta_ratio              float32 ()
  meta_tau_r_ns           float32 ()
  meta_tau_f_ns           float32 ()
  meta_T_max              float32 ()
  meta_t_peak_offset_ns   float32 ()
  meta_n_events           int32   ()
  meta_dt_min_ns          float32 ()               random-dt mode only
  meta_dt_max_ns          float32 ()               random-dt mode only

All arrays are pickle-free; np.load(path, allow_pickle=False) works.

Non-goals for v0.1 (§13)
-------------------------
  - Coloured/autocorrelated noise (pedestal-ACF)   [future: noise_model dispatch hook exists]
  - Cluster-coincident correlated pile-up
  - Poisson-random pile-up timing
  - Luminosity-scaling / rate-quote table
  - More than 2 pulses per event
  - Per-channel (τ_r, τ_f) overrides
  - Amplitude scan in absolute ADC
  - Any C++ modification
  - Any deconvolution
  - Plotting
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Generator, List, Optional, Tuple

import numpy as np

# _common (and by extension prad2py) is only needed for the real EVIO run path.
# The --self-test mode runs entirely in pure Python + numpy, so we defer this
# import to avoid a hard SystemExit when prad2py is not yet built.
try:
    import _common as C
    from _common import dec  # prad2py.dec re-export
    _PRAD2PY_AVAILABLE = True
except SystemExit:
    # _common raises SystemExit when prad2py is absent; catch it so
    # --self-test can run without the C++ bindings.
    C = None  # type: ignore[assignment]
    dec = None  # type: ignore[assignment]
    _PRAD2PY_AVAILABLE = False

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

VERSION = "pileup_generator.py v0.1"

DEFAULT_DT_GRID_NS: List[float] = [8.0, 12.0, 16.0, 20.0, 24.0,
                                    32.0, 48.0, 64.0, 96.0, 128.0]
DEFAULT_RATIO_GRID: List[float] = [0.1, 0.2, 0.3, 0.5, 1.0, 2.0]
DEFAULT_MATERIAL = "PbWO4"
MAX_ADC = 4095
SATURATION_SAFETY = 50          # MAX_ADC - SATURATION_SAFETY = 4045 hard limit
DEFAULT_POST_PAD = 40           # samples of headroom past the last injected peak
DEFAULT_HEIGHT_MIN = 500.0
DEFAULT_HEIGHT_MAX = 2000.0
DEFAULT_HEIGHT_RMS_MULT = 10.0
DEFAULT_T0_MIN_NS = 25.0        # §2b Physics Population A onset cut
DEFAULT_CHI2_MAX = 5.0
DEFAULT_RESIDUAL_VETO_SIGMA = 5.0
DEFAULT_N_PER_CONFIG = 50
DEFAULT_PROGRESS_EVERY = 500
DEFAULT_SEED = 12345


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class TemplateShape:
    """Resolved per-material template shape — cached once, used in hot loop."""
    material: str
    tau_r_ns: float
    tau_f_ns: float
    t0_template_ns: float           # median onset t0 from the fit output
    T_max: float                    # (1 - u) · u^(τ_r / τ_f), u = τ_r/(τ_r+τ_f)
    t_peak_offset_ns: float         # τ_r · log((τ_r+τ_f)/τ_r)
    source_path: str = ""


@dataclass
class BaseEvent:
    """One accepted base waveform + WaveAnalyzer fit diagnostics + provenance."""
    samples: np.ndarray             # uint16, shape (n_samples,) — raw ADC copy
    ped_mean: float                 # from wave_ana.analyze
    ped_rms: float                  # noise σ
    base_height: float              # pk.height (ADC, pedsub)
    base_peak_time_ns: float        # pk.time (ns) — WaveAnalyzer peak position
    base_onset_t0_ns: float         # fitted template onset from fit_pulse_shape
    base_peak_sample: int           # pk.pos (integer sample index)
    base_chi2_per_dof: float        # from fit_pulse_shape
    base_fit_converged: bool        # from fit_pulse_shape
    # Provenance
    base_roc_tag: int
    base_slot: int
    base_chan: int
    base_channel_id: str            # "<roc_tag>_<slot>_<chan>"
    base_channel_name: str          # e.g. "W123"
    module_type: str                # "PbWO4"
    base_source_file: str           # basename of the EVIO split
    base_physics_event_index: int   # physics-event counter within the run
    base_run_number: int


@dataclass
class GridPoint:
    """One (dt_ns, ratio) coordinate with its pre-allocated output buffers."""
    dt_ns: float
    ratio: float
    n_per_config: int
    n_samples: int

    # Pre-allocated buffers filled by .append()
    _waveforms: np.ndarray = field(init=False)
    _truth_heights_adc: np.ndarray = field(init=False)
    _truth_model_amplitudes: np.ndarray = field(init=False)
    _truth_onset_t0_ns: np.ndarray = field(init=False)
    _truth_peak_times_ns: np.ndarray = field(init=False)
    _event_dt_ns: np.ndarray = field(init=False)
    _pedestals: np.ndarray = field(init=False)
    _base_channel_names: List[str] = field(init=False)
    _base_channel_ids: List[str] = field(init=False)
    _base_roc_tags: List[int] = field(init=False)
    _base_slots: List[int] = field(init=False)
    _base_chans: List[int] = field(init=False)
    _base_run_numbers: List[int] = field(init=False)
    _base_source_files: List[str] = field(init=False)
    _base_physics_event_index: List[int] = field(init=False)
    _base_chi2_per_dof: List[float] = field(init=False)
    _n_filled: int = field(init=False, default=0)

    def __post_init__(self) -> None:
        N, S = self.n_per_config, self.n_samples
        self._waveforms             = np.zeros((N, S), dtype=np.uint16)
        self._truth_heights_adc     = np.zeros((N, 2), dtype=np.float32)
        self._truth_model_amplitudes = np.zeros((N, 2), dtype=np.float32)
        self._truth_onset_t0_ns     = np.zeros((N, 2), dtype=np.float32)
        self._truth_peak_times_ns   = np.zeros((N, 2), dtype=np.float32)
        self._event_dt_ns           = np.zeros(N, dtype=np.float32)
        self._pedestals             = np.zeros((N, 2), dtype=np.float32)
        self._base_channel_names    = [""] * N
        self._base_channel_ids      = [""] * N
        self._base_roc_tags         = [0] * N
        self._base_slots            = [0] * N
        self._base_chans            = [0] * N
        self._base_run_numbers      = [0] * N
        self._base_source_files     = [""] * N
        self._base_physics_event_index = [0] * N
        self._base_chi2_per_dof     = [0.0] * N
        self._n_filled = 0

    @property
    def key(self) -> Tuple[float, float]:
        return (self.dt_ns, self.ratio)

    @property
    def full(self) -> bool:
        return self._n_filled >= self.n_per_config

    @property
    def n_filled(self) -> int:
        return self._n_filled

    def append(self, waveform: np.ndarray, truth: Dict, base: "BaseEvent",
               event_dt_ns: float = 0.0) -> None:
        """Add one synthetic event.  Caller ensures not full.

        Parameters
        ----------
        waveform     : uint16 array (n_samples,) — composite ADC samples.
        truth        : dict with inject_* and base_* keys from inject_pulse().
        base         : BaseEvent provenance.
        event_dt_ns  : the actual ΔT used for this event (equals gp.dt_ns in
                       fixed-grid mode; varies per event in random-dt mode).
        """
        i = self._n_filled
        # Guard: if the actual waveform length differs from the pre-allocated
        # buffer (e.g. DAQ config reported 200 samples but real data has 100),
        # reallocate _waveforms on the very first append so subsequent appends
        # are consistent.  Resizing after data has already been stored is an
        # unrecoverable programming error.
        if waveform.shape[0] != self._waveforms.shape[1]:
            if self._n_filled != 0:
                raise ValueError(
                    f"waveform shape {waveform.shape} does not match "
                    f"pre-allocated buffer shape (n_samples={self._waveforms.shape[1]}); "
                    f"cannot resize after appends have started"
                )
            # First append — resize the buffer to the actual sample count.
            N = self.n_per_config
            self._waveforms = np.zeros((N, waveform.shape[0]), dtype=np.uint16)
            self.n_samples = waveform.shape[0]
        self._waveforms[i] = waveform
        self._event_dt_ns[i]               = np.float32(event_dt_ns)
        # col 0 = injected, col 1 = base_ref
        self._truth_heights_adc[i, 0]      = truth["inject_height_adc"]
        self._truth_heights_adc[i, 1]      = truth["base_height_adc"]
        self._truth_model_amplitudes[i, 0] = truth["inject_model_amp"]
        self._truth_model_amplitudes[i, 1] = truth["base_model_amp"]
        self._truth_onset_t0_ns[i, 0]      = truth["inject_onset_t0_ns"]
        self._truth_onset_t0_ns[i, 1]      = truth["base_onset_t0_ns"]
        self._truth_peak_times_ns[i, 0]    = truth["inject_peak_time_ns"]
        self._truth_peak_times_ns[i, 1]    = truth["base_peak_time_ns"]
        self._pedestals[i, 0]              = base.ped_mean
        self._pedestals[i, 1]             = base.ped_rms
        self._base_channel_names[i]        = base.base_channel_name
        self._base_channel_ids[i]          = base.base_channel_id
        self._base_roc_tags[i]             = base.base_roc_tag
        self._base_slots[i]                = base.base_slot
        self._base_chans[i]                = base.base_chan
        self._base_run_numbers[i]          = base.base_run_number
        self._base_source_files[i]         = base.base_source_file
        self._base_physics_event_index[i]  = base.base_physics_event_index
        self._base_chi2_per_dof[i]         = base.base_chi2_per_dof
        self._n_filled += 1

    def slice(self) -> Dict:
        """Return dict of arrays trimmed to _n_filled rows (for partial output)."""
        n = self._n_filled
        return dict(
            waveforms               = self._waveforms[:n],
            truth_heights_adc       = self._truth_heights_adc[:n],
            truth_model_amplitudes  = self._truth_model_amplitudes[:n],
            truth_onset_t0_ns       = self._truth_onset_t0_ns[:n],
            truth_peak_times_ns     = self._truth_peak_times_ns[:n],
            event_dt_ns             = self._event_dt_ns[:n],
            pedestals               = self._pedestals[:n],
            base_channel_names      = np.array(self._base_channel_names[:n], dtype="<U8"),
            base_channel_ids        = np.array(self._base_channel_ids[:n],   dtype="<U16"),
            base_roc_tags           = np.array(self._base_roc_tags[:n],      dtype=np.int32),
            base_slots              = np.array(self._base_slots[:n],         dtype=np.int32),
            base_chans              = np.array(self._base_chans[:n],         dtype=np.int32),
            base_run_numbers        = np.array(self._base_run_numbers[:n],   dtype=np.int32),
            base_source_files       = np.array(self._base_source_files[:n],  dtype="<U64"),
            base_physics_event_index= np.array(self._base_physics_event_index[:n], dtype=np.int32),
            base_chi2_per_dof       = np.array(self._base_chi2_per_dof[:n],  dtype=np.float32),
        )


# ---------------------------------------------------------------------------
# Template loading
# ---------------------------------------------------------------------------

def load_template_shape(json_path: str, material: str) -> TemplateShape:
    """Load (tau_r_ns, tau_f_ns, t0_ns) medians from a pulse-template JSON
    and precompute T_max and t_peak_offset_ns.

    Raises SystemExit on any error (missing file, missing material, NaN medians).
    """
    p = Path(json_path)
    if not p.is_file():
        raise SystemExit(
            f"[ERROR] template file not found: {json_path}"
        )
    try:
        data = json.loads(p.read_text())
    except json.JSONDecodeError as exc:
        raise SystemExit(f"[ERROR] template JSON parse error: {exc}")

    by_type = data.get("_by_type", {})
    if material not in by_type:
        avail = ", ".join(sorted(by_type.keys())) if by_type else "(none)"
        raise SystemExit(
            f"[ERROR] material {material!r} not found in template "
            f"'_by_type'; available: {avail}"
        )

    entry = by_type[material]
    try:
        tau_r = float(entry["tau_r_ns"]["median"])
        tau_f = float(entry["tau_f_ns"]["median"])
        t0    = float(entry["t0_ns"]["median"])
    except (KeyError, TypeError, ValueError) as exc:
        raise SystemExit(
            f"[ERROR] cannot read tau_r/tau_f/t0 medians for {material!r}: {exc}"
        )

    import math
    if math.isnan(tau_r) or math.isnan(tau_f) or math.isnan(t0):
        raise SystemExit(
            f"[ERROR] NaN in template medians for {material!r}: "
            f"tau_r={tau_r}, tau_f={tau_f}, t0={t0}"
        )
    if tau_r <= 0 or tau_f <= 0:
        raise SystemExit(
            f"[ERROR] non-positive tau for {material!r}: tau_r={tau_r}, tau_f={tau_f}"
        )

    u = tau_r / (tau_r + tau_f)
    T_max = (1.0 - u) * (u ** (tau_r / tau_f))
    t_peak_offset_ns = tau_r * math.log((tau_r + tau_f) / tau_r)

    return TemplateShape(
        material=material,
        tau_r_ns=tau_r,
        tau_f_ns=tau_f,
        t0_template_ns=t0,
        T_max=T_max,
        t_peak_offset_ns=t_peak_offset_ns,
        source_path=str(p.resolve()),
    )


# ---------------------------------------------------------------------------
# Template evaluation (Python-side — for injection kernel only)
# ---------------------------------------------------------------------------

def eval_template(t_ns: np.ndarray, t0_ns: float,
                  tau_r: float, tau_f: float, T_max: float) -> np.ndarray:
    """Evaluate the peak-normalised two-τ template on grid t_ns.

    Returns array of same shape as t_ns, dtype float32.
    Peak value = 1.0 by construction (divided by T_max).
    Returns zeros everywhere t_ns <= t0_ns.
    """
    out = np.zeros(t_ns.shape, dtype=np.float32)
    mask = t_ns > t0_ns
    if not mask.any():
        return out
    dt = (t_ns[mask] - t0_ns).astype(np.float32)
    raw = (1.0 - np.exp(-dt / tau_r)) * np.exp(-dt / tau_f)
    out[mask] = raw / T_max
    return out


# ---------------------------------------------------------------------------
# Injection kernel
# ---------------------------------------------------------------------------

def inject_pulse(
    base_samples_u16: np.ndarray,
    ped_mean: float,
    ped_rms: float,
    base_height: float,
    base_peak_time_ns: float,
    base_onset_t0_ns: float,
    dt_ns: float,
    ratio: float,
    shape: TemplateShape,
    clk_ns: float,
    noise_model: str,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, Dict]:
    """Compose a synthetic 2-pulse waveform by injecting one template pulse
    on top of a real base event.

    ΔT = dt_ns is peak-to-peak (§1a.1).  Both pulses share (τ_r, τ_f), so
    onset-to-onset separation == peak-to-peak separation.

    Returns
    -------
    (composite_u16, truth_dict)
      composite_u16 : uint16 array of shape (n_samples,)
      truth_dict    : 8 keys (inject_*, base_*) — see inline comments.
    """
    n = base_samples_u16.shape[0]
    t_ns = np.arange(n, dtype=np.float32) * clk_ns

    # Injected-pulse timing (§1a.1)
    inject_peak_time_ns = base_peak_time_ns + dt_ns
    inject_onset_t0_ns  = inject_peak_time_ns - shape.t_peak_offset_ns

    # Height convention (§1a.2): inject_height IS the peak height in ADC.
    inject_height = ratio * base_height

    pulse_unit = eval_template(t_ns, inject_onset_t0_ns,
                               shape.tau_r_ns, shape.tau_f_ns, shape.T_max)

    composite_f = (base_samples_u16.astype(np.float32)
                   + inject_height * pulse_unit)

    if noise_model == "none":
        pass
    elif noise_model == "gaussian-iid":
        # Explicit opt-in stress test only.  Artificially increases baseline
        # RMS (base already carries real detector noise — this double-counts).
        noise = rng.normal(0.0, ped_rms, size=n).astype(np.float32)
        composite_f += noise
    else:
        raise ValueError(f"unknown noise_model={noise_model!r}")

    # Saturation invariant — §4 cut 10 should have prevented this.
    if composite_f.max() >= MAX_ADC:
        raise AssertionError(
            f"saturation invariant violated: max={composite_f.max():.1f}; "
            f"§4 cut 10 should have rejected this base"
        )
    composite_f = np.clip(composite_f, 0.0, float(MAX_ADC - 1))  # floor only
    composite_u16 = np.rint(composite_f).astype(np.uint16)

    assert composite_u16.shape == base_samples_u16.shape
    assert composite_u16.dtype == np.uint16
    assert abs(inject_peak_time_ns - (base_peak_time_ns + dt_ns)) < 1e-6
    assert abs(inject_onset_t0_ns - (inject_peak_time_ns - shape.t_peak_offset_ns)) < 1e-6

    truth = {
        # Injected — exact truth (col 0)
        "inject_height_adc":   float(inject_height),
        "inject_model_amp":    float(inject_height / shape.T_max),
        "inject_onset_t0_ns":  float(inject_onset_t0_ns),
        "inject_peak_time_ns": float(inject_peak_time_ns),
        # Base reference — measured, not exact truth (col 1)
        "base_height_adc":     float(base_height),
        "base_model_amp":      float(base_height / shape.T_max),
        "base_onset_t0_ns":    float(base_onset_t0_ns),
        "base_peak_time_ns":   float(base_peak_time_ns),
    }
    return composite_u16, truth


# ---------------------------------------------------------------------------
# Base-event iteration
# ---------------------------------------------------------------------------

def iterate_base_events(
    pipeline,
    args,
    shape: TemplateShape,
    clk_ns: float,
    dt_max_ns: float,
) -> Generator["BaseEvent", None, None]:
    """Walk EVIO files in pipeline.evio_files and yield accepted BaseEvents.

    Applies all 10 cuts from §4 of the plan.  Each cut increments a named
    counter; the caller prints the rejection histogram at the end.

    Parameters
    ----------
    pipeline  : C._common.Pipeline — from setup_pipeline()
    args      : argparse.Namespace — CLI args (height_min, t0_min, etc.)
    shape     : TemplateShape — precomputed template parameters
    clk_ns    : float — ns per sample
    dt_max_ns : float — largest ΔT in the grid (for headroom cut)

    Yields
    ------
    BaseEvent
    """
    # Rejection counters (one per cut from §4)
    n_reject_multi_peak    = 0
    n_reject_quality       = 0
    n_reject_overflow      = 0
    n_reject_height        = 0
    n_reject_snr           = 0
    n_reject_material      = 0
    n_reject_fit_converged = 0
    n_reject_fit_t0        = 0
    n_reject_fit_chi2      = 0
    n_reject_residual_veto = 0
    n_reject_headroom      = 0
    n_reject_saturation    = 0
    n_total_candidates     = 0

    # Precompute headroom threshold (§4 cut 9): checked against dt_max once.
    # The injected peak lands near pk.pos + round(dt_max_ns / clk_ns).
    dt_max_samples = int(round(dt_max_ns / clk_ns))
    post_pad = int(args.post_pad)

    # Saturation ratio limit for check (§4 cut 10)
    ratio_max = float(max(args.ratio_grid))
    safety_adc = float(MAX_ADC - SATURATION_SAFETY)

    # Cache: (roc_tag, slot, chan) -> (name, chan_id, module_type)
    name_cache: Dict[Tuple[int, int, int], Tuple[str, str, str]] = {}

    def _key_for(roc_tag: int, crate: Optional[int], s: int, c: int
                 ) -> Tuple[str, str, str]:
        cached = name_cache.get((roc_tag, s, c))
        if cached is not None:
            return cached
        chan_id = f"{roc_tag}_{s}_{c}"
        name = chan_id
        mtype = "Unknown"
        if crate is not None:
            mod = pipeline.hycal.module_by_daq(crate, s, c)
            if mod is not None:
                name  = mod.name
                mtype = mod.type.name
        name_cache[(roc_tag, s, c)] = (name, chan_id, mtype)
        return name, chan_id, mtype

    ch = dec.EvChannel()
    ch.set_config(pipeline.cfg)

    n_phys = 0
    n_files_open = 0
    n_bases_yielded = 0

    # Residual veto window (§4 cut 8): ±3·τ_f around the main peak (in samples)
    veto_half_ns  = 3.0 * shape.tau_f_ns
    veto_half_smp = int(round(veto_half_ns / clk_ns))
    residual_sigma = float(args.residual_veto_sigma)

    for fpath in pipeline.evio_files:
        if ch.open_auto(fpath) != dec.Status.success:
            print(f"[WARN] skip (cannot open): {fpath}", file=sys.stderr, flush=True)
            continue
        n_files_open += 1
        print(f"[file {n_files_open}/{len(pipeline.evio_files)}] {fpath}", flush=True)
        done = False

        while ch.read() == dec.Status.success:
            if not ch.scan():
                continue
            if ch.get_event_type() != dec.EventType.Physics:
                continue

            for i in range(ch.get_n_events()):
                decoded = ch.decode_event(i, with_ssp=False)
                if not decoded["ok"]:
                    continue
                n_phys += 1
                if args.max_events and n_phys >= args.max_events:
                    done = True

                run_number = C.extract_run_number(fpath)
                fadc_evt = decoded["event"]

                for ri in range(fadc_evt.nrocs):
                    roc = fadc_evt.roc(ri)
                    if not roc.present:
                        continue
                    crate = pipeline.crate_map.get(roc.tag)
                    for s in roc.present_slots():
                        slot_obj = roc.slot(s)
                        for c in slot_obj.present_channels():
                            cd = slot_obj.channel(c)
                            if cd.nsamples <= 0:
                                continue

                            n_total_candidates += 1

                            samples_raw = np.asarray(cd.samples, dtype=np.uint16)
                            ped, rms, peaks = pipeline.wave_ana.analyze(samples_raw)

                            # Cut 1: exactly one peak
                            if len(peaks) != 1:
                                n_reject_multi_peak += 1
                                continue
                            pk = peaks[0]

                            # Cut 2: quality == 0 (no Q_PEAK_PILED etc.)
                            if pk.quality != 0:
                                n_reject_quality += 1
                                continue

                            # Cut 3: not overflow
                            if pk.overflow:
                                n_reject_overflow += 1
                                continue

                            # Cut 4: height in [height_min, height_max]
                            h = float(pk.height)
                            if not (args.height_min <= h <= args.height_max):
                                n_reject_height += 1
                                continue

                            # Cut 5: height >= height_rms_mult * ped_rms (SNR)
                            if h < args.height_rms_mult * float(rms):
                                n_reject_snr += 1
                                continue

                            # Cut 6: module type matches --material
                            name, chan_id, mtype = _key_for(roc.tag, crate, s, c)
                            if mtype != args.material:
                                n_reject_material += 1
                                continue

                            # Cut 7: fit-quality gate — run fit_pulse_shape
                            n_samples = samples_raw.shape[0]
                            pk_pos = int(pk.pos)
                            # Use a window around the peak (matching fit_pulse_template.py)
                            pre_smp  = min(pk_pos, 8)
                            post_smp = min(n_samples - pk_pos - 1, 40)
                            lo = pk_pos - pre_smp
                            hi = pk_pos + post_smp + 1
                            if lo < 0 or hi > n_samples:
                                n_reject_fit_converged += 1
                                continue

                            slice_u16 = samples_raw[lo:hi]
                            rel_peak  = pk_pos - lo

                            fit = dec.WaveAnalyzer.fit_pulse_shape(
                                slice_u16, rel_peak,
                                float(ped), float(rms), clk_ns, 0.03
                            )

                            if not fit.ok:
                                n_reject_fit_converged += 1
                                continue

                            # Fitted onset must be in the physics Population A window
                            if fit.t0_ns < args.t0_min:
                                n_reject_fit_t0 += 1
                                continue
                            # Must be inside the readout window
                            if not (0.0 < fit.t0_ns < n_samples * clk_ns):
                                n_reject_fit_t0 += 1
                                continue

                            if fit.chi2_per_dof > args.chi2_max:
                                n_reject_fit_chi2 += 1
                                continue

                            # Cut 8: residual hidden-pileup veto
                            # Reconstruct fitted single-pulse, subtract, check outside
                            # the ±veto_half_smp window around the main peak.
                            t_ns_arr = np.arange(n_samples, dtype=np.float32) * clk_ns
                            # The fit gave us t0_ns in slice-relative coords; convert
                            # back to full-window coords using the lo offset.
                            fit_t0_full_ns = fit.t0_ns + lo * clk_ns
                            fit_tau_r = float(shape.tau_r_ns)  # use material median shape
                            fit_tau_f = float(shape.tau_f_ns)

                            model_pulse = (h * eval_template(
                                t_ns_arr, fit_t0_full_ns,
                                fit_tau_r, fit_tau_f, shape.T_max
                            ))
                            # Residual = raw (pedsub) − fitted model
                            pedsub = samples_raw.astype(np.float32) - float(ped)
                            residual = pedsub - model_pulse

                            # Define the protected window around the main peak
                            veto_lo = max(0, pk_pos - veto_half_smp)
                            veto_hi = min(n_samples, pk_pos + veto_half_smp + 1)

                            # Scan outside the protected bracket
                            outside_mask = np.ones(n_samples, dtype=bool)
                            outside_mask[veto_lo:veto_hi] = False
                            threshold = residual_sigma * float(rms)
                            if float(rms) > 0 and np.any(np.abs(residual[outside_mask]) > threshold):
                                n_reject_residual_veto += 1
                                continue

                            # Cut 9: injection headroom — base peak + dt_max + post_pad < n_samples
                            if pk_pos + dt_max_samples + post_pad >= n_samples:
                                n_reject_headroom += 1
                                continue

                            # Cut 10: saturation rejection (§4 cut 10)
                            # Composite peak ≈ ped_mean + base_height + ratio_max * base_height
                            # We check the peak of the injected pulse added to the base
                            # (conservative: ped_mean + (1 + ratio_max) * base_height)
                            composite_peak_approx = float(ped) + (1.0 + ratio_max) * h
                            if composite_peak_approx >= safety_adc:
                                n_reject_saturation += 1
                                continue

                            # --- All cuts passed: accepted base event ---
                            # Compute base_peak_time_ns (in full-window coordinates)
                            base_peak_time_ns = float(pk.time)
                            # fit.t0_ns is slice-relative; convert to full window
                            base_onset_t0_ns = fit.t0_ns + lo * clk_ns

                            # Make a copy of samples so we don't alias C++ buffer
                            samples_copy = samples_raw.copy()

                            base = BaseEvent(
                                samples             = samples_copy,
                                ped_mean            = float(ped),
                                ped_rms             = float(rms),
                                base_height         = h,
                                base_peak_time_ns   = base_peak_time_ns,
                                base_onset_t0_ns    = base_onset_t0_ns,
                                base_peak_sample    = pk_pos,
                                base_chi2_per_dof   = float(fit.chi2_per_dof),
                                base_fit_converged  = bool(fit.ok),
                                base_roc_tag        = int(roc.tag),
                                base_slot           = int(s),
                                base_chan           = int(c),
                                base_channel_id     = chan_id,
                                base_channel_name   = name,
                                module_type         = mtype,
                                base_source_file    = Path(fpath).name,
                                base_physics_event_index = n_phys,
                                base_run_number     = run_number,
                            )
                            n_bases_yielded += 1
                            yield base

                if done:
                    break

        ch.close()
        if done:
            break

    # Store rejection counters as attributes on the generator frame so
    # the caller can retrieve them after exhaustion.  We use a workaround:
    # attach to a module-level dict that scan_grid can read.
    _rejection_counts.update(dict(
        n_total_candidates    = n_total_candidates,
        n_reject_multi_peak   = n_reject_multi_peak,
        n_reject_quality      = n_reject_quality,
        n_reject_overflow     = n_reject_overflow,
        n_reject_height       = n_reject_height,
        n_reject_snr          = n_reject_snr,
        n_reject_material     = n_reject_material,
        n_reject_fit_converged= n_reject_fit_converged,
        n_reject_fit_t0       = n_reject_fit_t0,
        n_reject_fit_chi2     = n_reject_fit_chi2,
        n_reject_residual_veto= n_reject_residual_veto,
        n_reject_headroom     = n_reject_headroom,
        n_reject_saturation   = n_reject_saturation,
        n_phys_events_scanned = n_phys,
        n_bases_yielded       = n_bases_yielded,
    ))


# Module-level dict for passing rejection stats out of the generator.
_rejection_counts: Dict[str, int] = {}


# ---------------------------------------------------------------------------
# Scanning driver — paired-base design (§7)
# ---------------------------------------------------------------------------

def scan_grid(
    base_iter: Generator,
    grid_points: List[GridPoint],
    n_per_config: int,
    shape: TemplateShape,
    clk_ns: float,
    noise_model: str,
    rng_map: Dict[Tuple[float, float], np.random.Generator],
    t0_wall: float,
    progress_every: int,
) -> int:
    """Paired-base scan: each accepted BaseEvent is reused across ALL grid cells.

    Parameters
    ----------
    base_iter      : generator of BaseEvent
    grid_points    : list of GridPoint (one per (dt, ratio) cell)
    n_per_config   : how many base events to collect
    shape          : TemplateShape
    clk_ns         : ns per sample
    noise_model    : "none" or "gaussian-iid"
    rng_map        : dict keyed by (dt_ns, ratio) -> np.random.Generator
    t0_wall        : wall-clock start time (for progress output)
    progress_every : print progress every N accepted bases

    Returns
    -------
    n_bases_done : int
    """
    n_bases_done = 0
    next_progress = progress_every

    for base in base_iter:
        if n_bases_done >= n_per_config:
            break
        for gp in grid_points:
            wf, truth = inject_pulse(
                base.samples, base.ped_mean, base.ped_rms,
                base.base_height, base.base_peak_time_ns,
                base.base_onset_t0_ns,
                gp.dt_ns, gp.ratio, shape, clk_ns,
                noise_model, rng_map[gp.key]
            )
            gp.append(wf, truth, base, event_dt_ns=gp.dt_ns)
        n_bases_done += 1

        if n_bases_done >= next_progress:
            elapsed = time.monotonic() - t0_wall
            rate = n_bases_done / elapsed if elapsed > 0 else 0.0
            cells_full = sum(1 for gp in grid_points if gp.full)
            print(
                f"  [progress] bases_kept={n_bases_done}/{n_per_config}  "
                f"cells_full={cells_full}/{len(grid_points)}  "
                f"rate={rate:.0f} base/s  elapsed={elapsed:.1f}s",
                flush=True,
            )
            while next_progress <= n_bases_done:
                next_progress += progress_every

    return n_bases_done


def scan_grid_random_dt(
    base_iter: Generator,
    grid_points: List[GridPoint],
    n_per_config: int,
    shape: TemplateShape,
    clk_ns: float,
    noise_model: str,
    rng_map: Dict[float, np.random.Generator],
    dt_min: float,
    dt_max: float,
    t0_wall: float,
    progress_every: int,
) -> int:
    """Random-ΔT scan: each accepted BaseEvent is reused across all ratio cells,
    with a fresh ΔT drawn uniformly from [dt_min, dt_max] for every event.

    Parameters
    ----------
    base_iter      : generator of BaseEvent
    grid_points    : list of GridPoint (one per ratio; gp.dt_ns is a sentinel NaN)
    n_per_config   : events to collect per ratio point
    shape          : TemplateShape
    clk_ns         : ns per sample
    noise_model    : "none" or "gaussian-iid"
    rng_map        : dict keyed by ratio -> np.random.Generator
    dt_min         : minimum ΔT in ns
    dt_max         : maximum ΔT in ns
    t0_wall        : wall-clock start time (for progress output)
    progress_every : print progress every N accepted bases

    Returns
    -------
    n_bases_done : int
    """
    n_bases_done = 0
    next_progress = progress_every

    for base in base_iter:
        if n_bases_done >= n_per_config:
            break
        for gp in grid_points:
            rng = rng_map[gp.ratio]
            dt_ns = float(rng.uniform(dt_min, dt_max))
            wf, truth = inject_pulse(
                base.samples, base.ped_mean, base.ped_rms,
                base.base_height, base.base_peak_time_ns,
                base.base_onset_t0_ns,
                dt_ns, gp.ratio, shape, clk_ns,
                noise_model, rng,
            )
            gp.append(wf, truth, base, event_dt_ns=dt_ns)
        n_bases_done += 1

        if n_bases_done >= next_progress:
            elapsed = time.monotonic() - t0_wall
            rate = n_bases_done / elapsed if elapsed > 0 else 0.0
            cells_full = sum(1 for gp in grid_points if gp.full)
            print(
                f"  [progress] bases_kept={n_bases_done}/{n_per_config}  "
                f"cells_full={cells_full}/{len(grid_points)}  "
                f"rate={rate:.0f} base/s  elapsed={elapsed:.1f}s",
                flush=True,
            )
            while next_progress <= n_bases_done:
                next_progress += progress_every

    return n_bases_done


# ---------------------------------------------------------------------------
# Filename helpers
# ---------------------------------------------------------------------------

def _npz_name(material: str, dt_ns: float, ratio: float) -> str:
    """Canonical per-cell filename: pileup_PbWO4_dt016_ratio0500.npz"""
    dt_int    = int(round(dt_ns))
    ratio_int = int(round(ratio * 1000))
    return f"pileup_{material}_dt{dt_int:03d}_ratio{ratio_int:04d}.npz"


def _npz_name_random_dt(material: str, ratio: float) -> str:
    """Filename for random-ΔT mode: pileup_PbWO4_dtRandom_ratio0500.npz"""
    ratio_int = int(round(ratio * 1000))
    return f"pileup_{material}_dtRandom_ratio{ratio_int:04d}.npz"


# ---------------------------------------------------------------------------
# Serialisation
# ---------------------------------------------------------------------------

def write_grid_point(
    out_dir: Path,
    gp: GridPoint,
    shape: TemplateShape,
    meta_common: Dict,
    metadata_in_npz: bool,
    random_dt: bool = False,
) -> Tuple[Path, Path]:
    """Write one .npz + one sidecar .json for a single grid-point cell.

    Parameters
    ----------
    random_dt : bool
        When True the cell was produced in random-ΔT mode.  The filename uses
        "dtRandom" instead of a fixed dt value, meta_dt_ns is stored as NaN,
        and meta_dt_min_ns / meta_dt_max_ns scalars are added.

    Returns (npz_path, json_path).
    """
    n_ev = gp.n_filled
    arrays = gp.slice()

    if random_dt:
        stem = _npz_name_random_dt(shape.material, gp.ratio)
    else:
        stem = _npz_name(shape.material, gp.dt_ns, gp.ratio)
    stem      = stem[:-4]  # strip .npz
    npz_path  = out_dir / (stem + ".npz")
    json_path = out_dir / (stem + ".json")

    # In random-dt mode dt_ns is not fixed; represent as null in the sidecar.
    sidecar_dt_ns = None if random_dt else gp.dt_ns
    sidecar: Dict = {
        "material":          shape.material,
        "dt_ns":             sidecar_dt_ns,
        "amplitude_ratio":   gp.ratio,
        "tau_r_ns":          shape.tau_r_ns,
        "tau_f_ns":          shape.tau_f_ns,
        "T_max":             shape.T_max,
        "t_peak_offset_ns":  shape.t_peak_offset_ns,
        "template_t0_ns":    shape.t0_template_ns,
        "template_source":   shape.source_path,
        "n_events":          n_ev,
        "n_samples":         gp.n_samples,
        "clk_ns":            meta_common["clk_ns"],
        "noise_model":       meta_common["noise_model"],
        "noise_sigma_source": (
            "n/a (base pulse carries the readout noise)"
            if meta_common["noise_model"] == "none"
            else "ped_rms from WaveAnalyzer.analyze on each base event"
        ),
        "seed":              meta_common["seed"],
        "generator_version": VERSION,
        "generated_utc":     meta_common["generated_utc"],
        "base_event_source": {
            "runs":                  meta_common.get("runs", []),
            "n_evio_splits":         meta_common.get("n_evio_splits", 0),
            "cut_height_min":        meta_common["cut_height_min"],
            "cut_height_max":        meta_common["cut_height_max"],
            "cut_height_rms_mult":   meta_common["cut_height_rms_mult"],
            "cut_t0_min_ns":         meta_common["cut_t0_min_ns"],
            "cut_chi2_per_dof_max":  meta_common["cut_chi2_per_dof_max"],
            "residual_veto_sigma":   meta_common["residual_veto_sigma"],
        },
        "column_convention": {
            "truth_columns": ["injected", "base_ref"],
            "notes": (
                "col 0 = exact injection truth; "
                "col 1 = base_reference observables from WaveAnalyzer + "
                "fit_pulse_shape (accurate but not exact)."
            ),
        },
    }
    if random_dt:
        sidecar["dt_min_ns"] = meta_common["dt_min_ns"]
        sidecar["dt_max_ns"] = meta_common["dt_max_ns"]
        sidecar["random_dt"] = True

    # Build the npz save-dict
    # meta_dt_ns: the fixed ΔT for fixed-grid mode; NaN sentinel for random-dt mode.
    meta_dt_value = np.float32("nan") if random_dt else np.float32(gp.dt_ns)
    save_dict = dict(
        waveforms                = arrays["waveforms"],
        truth_heights_adc        = arrays["truth_heights_adc"],
        truth_model_amplitudes   = arrays["truth_model_amplitudes"],
        truth_onset_t0_ns        = arrays["truth_onset_t0_ns"],
        truth_peak_times_ns      = arrays["truth_peak_times_ns"],
        event_dt_ns              = arrays["event_dt_ns"],
        truth_source             = np.array(["injected", "base_ref"], dtype="<U9"),
        pedestals                = arrays["pedestals"],
        base_channel_names       = arrays["base_channel_names"],
        base_channel_ids         = arrays["base_channel_ids"],
        base_roc_tags            = arrays["base_roc_tags"],
        base_slots               = arrays["base_slots"],
        base_chans               = arrays["base_chans"],
        base_run_numbers         = arrays["base_run_numbers"],
        base_source_files        = arrays["base_source_files"],
        base_physics_event_index = arrays["base_physics_event_index"],
        base_chi2_per_dof        = arrays["base_chi2_per_dof"],
        meta_dt_ns               = meta_dt_value,
        meta_ratio               = np.float32(gp.ratio),
        meta_tau_r_ns            = np.float32(shape.tau_r_ns),
        meta_tau_f_ns            = np.float32(shape.tau_f_ns),
        meta_T_max               = np.float32(shape.T_max),
        meta_t_peak_offset_ns    = np.float32(shape.t_peak_offset_ns),
        meta_n_events            = np.int32(n_ev),
    )
    if random_dt:
        save_dict["meta_dt_min_ns"] = np.float32(meta_common["dt_min_ns"])
        save_dict["meta_dt_max_ns"] = np.float32(meta_common["dt_max_ns"])

    if metadata_in_npz:
        save_dict["metadata_json"] = np.array(json.dumps(sidecar))

    np.savez(str(npz_path), **save_dict)

    json_path.write_text(json.dumps(sidecar, indent=2))

    print(f"[write] {npz_path.name} ({n_ev} events) sidecar={json_path.name}",
          flush=True)
    return npz_path, json_path


def write_manifest(
    out_dir: Path,
    shape: TemplateShape,
    meta_common: Dict,
    grid_points: List[GridPoint],
    file_entries: List[Dict],
    cli_argv: str,
) -> Path:
    """Write top-level manifest.json."""
    random_dt = meta_common.get("random_dt", False)
    manifest: Dict = {
        "generator_version": VERSION,
        "generated_utc":     meta_common["generated_utc"],
        "material":          shape.material,
        "template_source":   shape.source_path,
        "clk_ns":            meta_common["clk_ns"],
        "n_samples":         meta_common.get("n_samples", 0),
        "ratio_grid":        sorted({gp.ratio for gp in grid_points}),
        "n_per_config":      meta_common["n_per_config"],
        "noise_model":       meta_common["noise_model"],
        "T_max":             shape.T_max,
        "t_peak_offset_ns":  shape.t_peak_offset_ns,
        "files":             file_entries,
        "seed":              meta_common["seed"],
        "cli":               cli_argv,
    }
    if random_dt:
        manifest["random_dt"]  = True
        manifest["dt_min_ns"]  = meta_common["dt_min_ns"]
        manifest["dt_max_ns"]  = meta_common["dt_max_ns"]
        manifest["dt_grid_ns"] = None
    else:
        manifest["dt_grid_ns"] = sorted({gp.dt_ns for gp in grid_points})
    mpath = out_dir / "manifest.json"
    mpath.write_text(json.dumps(manifest, indent=2))
    print(f"[manifest] {mpath}", flush=True)
    return mpath


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_output(npz_path: Path, wave_ana=None) -> None:
    """Schema-only validation (no C++ deps) when wave_ana is None.
    Adds a peak-detect smoke test when wave_ana is provided.

    Raises AssertionError with a descriptive message on any schema failure.
    """
    # Check 9 first: must load without pickle.
    data = np.load(str(npz_path), allow_pickle=False)

    # Check 1: waveforms dtype and shape
    wf = data["waveforms"]
    assert wf.dtype == np.uint16, \
        f"waveforms dtype {wf.dtype} != uint16 in {npz_path.name}"
    N, S = wf.shape
    assert N > 0 and S > 0, \
        f"waveforms shape {wf.shape} invalid in {npz_path.name}"

    # Check 2: truth_heights_adc shape, finite, positive
    th = data["truth_heights_adc"]
    assert th.shape == (N, 2), \
        f"truth_heights_adc shape {th.shape} != ({N}, 2) in {npz_path.name}"
    assert np.all(np.isfinite(th)), \
        f"truth_heights_adc has non-finite values in {npz_path.name}"
    assert np.all(th > 0), \
        f"truth_heights_adc has non-positive values in {npz_path.name}"

    # Check 3: shape consistency for onset, peak-time, and per-event dt arrays
    to = data["truth_onset_t0_ns"]
    tp = data["truth_peak_times_ns"]
    assert to.shape == (N, 2), \
        f"truth_onset_t0_ns shape {to.shape} != ({N}, 2) in {npz_path.name}"
    assert tp.shape == (N, 2), \
        f"truth_peak_times_ns shape {tp.shape} != ({N}, 2) in {npz_path.name}"
    assert "event_dt_ns" in data, \
        f"event_dt_ns array missing from {npz_path.name}"
    assert data["event_dt_ns"].shape == (N,), \
        f"event_dt_ns shape {data['event_dt_ns'].shape} != ({N},) in {npz_path.name}"

    # Check 4: ΔT round-trip (peak-to-peak of col0 vs col1 == meta_dt_ns).
    # The injected peak time is set by formula: inject_peak = base_peak + dt_ns.
    # However, base_peak_time_ns (col 1) comes from WaveAnalyzer pk.time, which
    # carries interpolation noise — so the difference won't be exact to 1e-4 ns.
    # In random-dt mode meta_dt_ns is NaN; check event_dt_ns bounds instead.
    meta_dt_raw = data["meta_dt_ns"]
    if np.isnan(float(meta_dt_raw)):
        # Random-dt mode: check each event's dt lies within [dt_min, dt_max]
        meta_dt_min = float(data["meta_dt_min_ns"])
        meta_dt_max = float(data["meta_dt_max_ns"])
        event_dt = data["event_dt_ns"]
        assert np.all(event_dt >= meta_dt_min - 0.01) and np.all(event_dt <= meta_dt_max + 0.01), \
            (f"ΔT range check failed in {npz_path.name}: "
             f"expected [{meta_dt_min}, {meta_dt_max}], "
             f"got [{event_dt.min():.4f}, {event_dt.max():.4f}]")
    else:
        meta_dt = float(meta_dt_raw)
        dt_actual = tp[:, 0] - tp[:, 1]
        assert np.allclose(dt_actual, meta_dt, atol=1.0), \
            (f"ΔT round-trip failed in {npz_path.name}: "
             f"expected {meta_dt}, got range [{dt_actual.min():.4f}, {dt_actual.max():.4f}]")

    # Check 5: ratio round-trip
    meta_ratio = float(data["meta_ratio"])
    ratio_actual = th[:, 0] / th[:, 1]
    assert np.allclose(ratio_actual, meta_ratio, rtol=1e-4), \
        (f"ratio round-trip failed in {npz_path.name}: "
         f"expected {meta_ratio}, got range [{ratio_actual.min():.4f}, {ratio_actual.max():.4f}]")

    # Check 6: peak_time - onset_t0 == t_peak_offset for col 0 (injected pulse).
    # Only the injected pulse (col 0) has inject_peak_time = inject_onset_t0 +
    # t_peak_offset by construction.  The base-reference pulse (col 1) has
    # independently measured peak_time (from WaveAnalyzer pk.time) and onset_t0
    # (from fit_pulse_shape), which don't satisfy this relationship exactly on
    # real noisy data.
    meta_tpo = float(data["meta_t_peak_offset_ns"])
    offset_injected = tp[:, 0] - to[:, 0]
    assert np.allclose(offset_injected, meta_tpo, atol=1e-4), \
        (f"t_peak_offset round-trip failed for injected pulses in {npz_path.name}: "
         f"expected {meta_tpo}, got range "
         f"[{offset_injected.min():.4f}, {offset_injected.max():.4f}]")

    # Check 7: truth_model_amplitudes * T_max == truth_heights_adc
    meta_T_max = float(data["meta_T_max"])
    tma = data["truth_model_amplitudes"]
    assert tma.shape == (N, 2), \
        f"truth_model_amplitudes shape {tma.shape} != ({N}, 2) in {npz_path.name}"
    assert np.allclose(tma * meta_T_max, th, rtol=1e-4), \
        f"truth_model_amplitudes × T_max != truth_heights_adc in {npz_path.name}"

    # Check 8: saturation invariant
    sat_limit = MAX_ADC - SATURATION_SAFETY
    max_wf = int(wf.max())
    assert max_wf < sat_limit, \
        (f"saturation invariant violated in {npz_path.name}: "
         f"waveforms.max()={max_wf} >= {sat_limit}")

    # Check 10: sidecar JSON exists and parses; n_events matches
    sidecar = npz_path.with_suffix(".json")
    if sidecar.is_file():
        try:
            meta = json.loads(sidecar.read_text())
        except json.JSONDecodeError as exc:
            raise AssertionError(
                f"sidecar JSON parse error for {npz_path.name}: {exc}"
            )
        assert meta.get("n_events") == N, \
            (f"sidecar n_events={meta.get('n_events')} != array shape N={N} "
             f"in {npz_path.name}")
    elif "metadata_json" in data:
        try:
            meta = json.loads(str(data["metadata_json"]))
        except (json.JSONDecodeError, KeyError) as exc:
            raise AssertionError(
                f"embedded metadata_json parse error for {npz_path.name}: {exc}"
            )
        assert meta.get("n_events") == N

    # Check 11 (optional): peak-detect smoke test on 5% subsample
    if wave_ana is not None:
        rng = np.random.default_rng(seed=42)
        n_sample = max(1, int(round(N * 0.05)))
        idx = rng.choice(N, size=n_sample, replace=False)
        n_with_peaks = 0
        for i in idx:
            _, _, pks = wave_ana.analyze(wf[i])
            if len(pks) >= 1:
                n_with_peaks += 1
        assert n_with_peaks > 0, \
            (f"smoke test: none of {n_sample} sampled waveforms had any detected "
             f"peak in {npz_path.name}")


# ---------------------------------------------------------------------------
# Self-test (§10 last paragraph)
# ---------------------------------------------------------------------------

def run_self_test() -> None:
    """Schema-only smoke test — no EVIO, DAQ config, or HyCal map required.

    Builds a synthetic base waveform from the PbWO4 template, runs the
    full inject → serialise → reload → schema-validate cycle on a 2×2 grid
    in a tempfile.TemporaryDirectory, and exits 0 if all checks pass.
    """
    print("[self-test] starting in-memory schema smoke test …", flush=True)

    # Use hardcoded PbWO4 medians (plan §1, from pulse_templates_025308_h500_t0cut.json)
    tau_r = 2.2356
    tau_f = 24.009
    import math
    u     = tau_r / (tau_r + tau_f)
    T_max = (1.0 - u) * (u ** (tau_r / tau_f))
    t_peak_offset = tau_r * math.log((tau_r + tau_f) / tau_r)

    shape = TemplateShape(
        material          = "PbWO4",
        tau_r_ns          = tau_r,
        tau_f_ns          = tau_f,
        t0_template_ns    = 26.279,
        T_max             = T_max,
        t_peak_offset_ns  = t_peak_offset,
        source_path       = "(self-test synthetic)",
    )

    # Validate precomputed shape constants (plan §10 Phase 1 self-check)
    assert abs(tau_r - 2.24) < 0.5,    f"tau_r={tau_r} out of expected range"
    assert abs(tau_f - 24.0) < 0.5,    f"tau_f={tau_f} out of expected range"
    assert abs(t_peak_offset - 5.5) < 0.5, \
        f"t_peak_offset={t_peak_offset:.3f} expected ~5.5 ns"
    # T_max for PbWO4 (τ_r≈2.24, τ_f≈24.0): formula gives ~0.727.
    # The plan's sidecar example (0.5842) is a documentation error in the plan;
    # the formula T_max = (1−u)·u^(τ_r/τ_f) is correct, the example is not.
    assert abs(T_max - 0.727) < 0.05,  f"T_max={T_max:.4f} expected ~0.727 for PbWO4"

    # Build a synthetic base waveform: pure template at h=1000, peak at ~100 ns
    clk_ns     = 4.0
    n_samples  = 200
    ped_val    = 500.0
    ped_rms_v  = 8.0
    base_h     = 1000.0
    base_peak_ns = 100.0                    # chosen peak time (ns)
    base_onset_ns = base_peak_ns - t_peak_offset

    t_ns = np.arange(n_samples, dtype=np.float32) * clk_ns
    template_vals = eval_template(t_ns, base_onset_ns, tau_r, tau_f, T_max)
    base_float = ped_val + base_h * template_vals
    base_samples = np.clip(np.rint(base_float), 0, MAX_ADC - 1).astype(np.uint16)

    # Test inject_pulse: ΔT=16 ns, ratio=0.5
    rng0 = np.random.default_rng(seed=99)
    dt_test, ratio_test = 16.0, 0.5
    wf_test, truth_test = inject_pulse(
        base_samples, ped_val, ped_rms_v,
        base_h, base_peak_ns, base_onset_ns,
        dt_test, ratio_test, shape, clk_ns,
        "none", rng0,
    )
    # Injected peak should be near sample (100 + 16) / 4 = 29 past the base peak
    inject_peak_smp = int(round((base_peak_ns + dt_test) / clk_ns))
    wf_float = wf_test.astype(np.float32) - ped_val
    region = wf_float[max(0, inject_peak_smp - 3): inject_peak_smp + 4]
    assert region.max() > 400.0, \
        f"inject smoke: peak near sample {inject_peak_smp} is only {region.max():.1f} ADC"

    # Smoke-test gaussian-iid branch
    rng1 = np.random.default_rng(seed=101)
    wf_noisy, _ = inject_pulse(
        base_samples, ped_val, ped_rms_v,
        base_h, base_peak_ns, base_onset_ns,
        dt_test, ratio_test, shape, clk_ns,
        "gaussian-iid", rng1,
    )
    assert wf_noisy.dtype == np.uint16

    # Build a 2×2 grid in a tempdir and run the full cycle
    dt_grid    = [16.0, 32.0]
    ratio_grid = [0.5, 1.0]
    n_per      = 3
    seed       = 42

    ss = np.random.SeedSequence(seed)
    grid_points: List[GridPoint] = []
    for dt in dt_grid:
        for r in ratio_grid:
            grid_points.append(GridPoint(dt_ns=dt, ratio=r,
                                         n_per_config=n_per, n_samples=n_samples))

    rng_map = {}
    child_seeds = ss.spawn(len(grid_points))
    for gp, child in zip(grid_points, child_seeds):
        rng_map[gp.key] = np.random.default_rng(child)

    # Fill grid by running n_per synthetic bases
    rng_base = np.random.default_rng(seed=77)
    for _ in range(n_per):
        for gp in grid_points:
            wf_g, truth_g = inject_pulse(
                base_samples, ped_val, ped_rms_v,
                base_h, base_peak_ns, base_onset_ns,
                gp.dt_ns, gp.ratio, shape, clk_ns,
                "none", rng_map[gp.key],
            )
            fake_base = BaseEvent(
                samples                  = base_samples.copy(),
                ped_mean                 = ped_val,
                ped_rms                  = ped_rms_v,
                base_height              = base_h,
                base_peak_time_ns        = base_peak_ns,
                base_onset_t0_ns         = base_onset_ns,
                base_peak_sample         = int(base_peak_ns / clk_ns),
                base_chi2_per_dof        = 0.9,
                base_fit_converged       = True,
                base_roc_tag             = 1,
                base_slot                = 3,
                base_chan                = 0,
                base_channel_id          = "1_3_0",
                base_channel_name        = "W001",
                module_type              = "PbWO4",
                base_source_file         = "selftest.evio.00000",
                base_physics_event_index = 1,
                base_run_number          = 0,
            )
            gp.append(wf_g, truth_g, fake_base, event_dt_ns=gp.dt_ns)

    meta_common = {
        "clk_ns":             clk_ns,
        "noise_model":        "none",
        "seed":               seed,
        "generated_utc":      datetime.now(timezone.utc).isoformat(),
        "n_per_config":       n_per,
        "cut_height_min":     DEFAULT_HEIGHT_MIN,
        "cut_height_max":     DEFAULT_HEIGHT_MAX,
        "cut_height_rms_mult":DEFAULT_HEIGHT_RMS_MULT,
        "cut_t0_min_ns":      DEFAULT_T0_MIN_NS,
        "cut_chi2_per_dof_max":DEFAULT_CHI2_MAX,
        "residual_veto_sigma":DEFAULT_RESIDUAL_VETO_SIGMA,
        "random_dt":          False,
        "runs":               [0],
        "n_evio_splits":      1,
        "n_samples":          n_samples,
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        file_entries = []
        for gp in grid_points:
            npz_p, json_p = write_grid_point(
                tmp, gp, shape, meta_common, metadata_in_npz=False
            )
            validate_output(npz_p, wave_ana=None)
            file_entries.append({
                "dt_ns":     gp.dt_ns,
                "ratio":     gp.ratio,
                "path":      npz_p.name,
                "sidecar":   json_p.name,
                "n_written": gp.n_filled,
            })

        write_manifest(tmp, shape, meta_common, grid_points, file_entries,
                       cli_argv="--self-test")

    print("[self-test] OK", flush=True)
    raise SystemExit(0)


# ---------------------------------------------------------------------------
# CLI argument parsing
# ---------------------------------------------------------------------------

def _parse_float_grid(s: str) -> List[float]:
    """Parse a comma-separated list of positive floats for dt-grid/ratio-grid."""
    try:
        vals = [float(x) for x in s.split(",") if x.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"grid must be a comma-separated list of positive floats; got {s!r}: {exc}"
        )
    if not vals or any(v <= 0 for v in vals):
        raise argparse.ArgumentTypeError(
            f"grid must be a non-empty comma-separated list of positive floats; got {s!r}"
        )
    return sorted(set(vals))


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="pileup_generator.py",
        description=(
            "Deterministic (ΔT, ratio) pile-up response-surface generator "
            "for PbWO4 HyCal modules (PRad-II).  v0.1 — single-channel, "
            "2-pulse (base + 1 injected) synthetic events."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Positional — EVIO inputs
    ap.add_argument(
        "evio_paths",
        nargs="*",
        metavar="EVIO_PATH",
        help="EVIO input file(s).  Each accepts glob, directory, or single split.",
    )

    # Required (unless --self-test)
    ap.add_argument("--template", default="",
                    help="Path to pulse_templates JSON file (required).")
    ap.add_argument("--out-dir", default="",
                    help="Output directory (required unless --self-test).")

    # Material / grid
    ap.add_argument("--material", default=DEFAULT_MATERIAL,
                    help=f"Module type to use (default: {DEFAULT_MATERIAL}).")
    ap.add_argument("--n-per-config", type=int, default=DEFAULT_N_PER_CONFIG,
                    help=f"Base events per grid cell (default: {DEFAULT_N_PER_CONFIG}).")
    ap.add_argument("--dt-grid", type=_parse_float_grid,
                    default=DEFAULT_DT_GRID_NS,
                    help="Comma-sep ΔT values in ns (default: 8,12,16,…,128).")
    ap.add_argument("--ratio-grid", type=_parse_float_grid,
                    default=DEFAULT_RATIO_GRID,
                    help="Comma-sep amplitude ratios (default: 0.1,0.2,0.3,0.5,1.0,2.0).")
    ap.add_argument("--random-dt", action="store_true",
                    help="Instead of the fixed --dt-grid, draw each injected "
                         "pulse's ΔT uniformly from [dt-min, dt-max] (in ns). "
                         "When set, --dt-grid is ignored and --n-per-config "
                         "events are produced per ratio grid point.")
    ap.add_argument("--dt-min", type=float, default=4.0,
                    help="Minimum ΔT in ns for --random-dt mode (default: 4.0).")
    ap.add_argument("--dt-max", type=float, default=200.0,
                    help="Maximum ΔT in ns for --random-dt mode (default: 200.0).")

    # Base-event selection cuts
    ap.add_argument("--height-min", type=float, default=DEFAULT_HEIGHT_MIN,
                    help=f"Min base-peak height in ADC (default: {DEFAULT_HEIGHT_MIN}).")
    ap.add_argument("--height-max", type=float, default=DEFAULT_HEIGHT_MAX,
                    help=f"Max base-peak height in ADC (default: {DEFAULT_HEIGHT_MAX}).")
    ap.add_argument("--height-rms-mult", type=float, default=DEFAULT_HEIGHT_RMS_MULT,
                    help=f"Min height as multiple of ped RMS (default: {DEFAULT_HEIGHT_RMS_MULT}).")
    ap.add_argument("--t0-min", type=float, default=DEFAULT_T0_MIN_NS,
                    help=f"Min fitted onset t0 in ns (default: {DEFAULT_T0_MIN_NS}).")
    ap.add_argument("--chi2-max", type=float, default=DEFAULT_CHI2_MAX,
                    help=f"Max chi2/dof for fit-quality gate (default: {DEFAULT_CHI2_MAX}).")
    ap.add_argument("--residual-veto-sigma", type=float,
                    default=DEFAULT_RESIDUAL_VETO_SIGMA,
                    help=f"Residual veto threshold in ped-RMS units "
                         f"(default: {DEFAULT_RESIDUAL_VETO_SIGMA}).")
    ap.add_argument("--post-pad", type=int, default=DEFAULT_POST_PAD,
                    help=f"Headroom samples past last injected peak (default: {DEFAULT_POST_PAD}).")

    # Noise / reproducibility
    ap.add_argument("--noise-model", choices=["none", "gaussian-iid"], default="none",
                    help="Noise model: none (default) or gaussian-iid (stress test).")
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED,
                    help=f"Global RNG seed for reproducibility (default: {DEFAULT_SEED}).")

    # Run control
    ap.add_argument("--max-events", type=int, default=0,
                    help="Cap on physics events scanned (0 = all).")
    ap.add_argument("--progress-every", type=int, default=DEFAULT_PROGRESS_EVERY,
                    help=f"Print progress every N accepted bases "
                         f"(default: {DEFAULT_PROGRESS_EVERY}).")

    # Infrastructure
    ap.add_argument("--daq-config", default="",
                    help="DAQ config path (default: installed default).")
    ap.add_argument("--hc-map-file", default="",
                    help="HyCal modules map path (default: database lookup).")
    ap.add_argument("--force", action="store_true",
                    help="Allow overwriting non-empty --out-dir.")
    ap.add_argument("--metadata-in-npz", action="store_true",
                    help="Embed metadata JSON string inside .npz (no separate sidecar).")
    ap.add_argument("--self-test", action="store_true",
                    help="Schema-only in-memory smoke test; exit 0 if passed.")
    return ap


# ---------------------------------------------------------------------------
# main()
# ---------------------------------------------------------------------------

def main() -> int:
    ap = build_parser()
    args = ap.parse_args()

    # --self-test: run standalone, no EVIO/DAQ/HyCal needed.
    if args.self_test:
        run_self_test()
        return 0  # unreachable — run_self_test raises SystemExit(0)

    # --- Validate required args for normal run ---
    if not _PRAD2PY_AVAILABLE:
        raise SystemExit(
            "[ERROR] cannot import prad2py (C++ bindings not found).\n"
            "        Build the python bindings (cmake -DBUILD_PYTHON=ON) and "
            "ensure the install directory is on PYTHONPATH.\n"
            "        (Tip: --self-test runs without prad2py.)"
        )
    if not args.template:
        ap.error("--template is required (unless --self-test)")
    if not args.out_dir:
        ap.error("--out-dir is required (unless --self-test)")
    if not args.evio_paths:
        ap.error("at least one EVIO_PATH is required (unless --self-test)")

    t0_wall = time.monotonic()
    cli_argv = " ".join(sys.argv)
    generated_utc = datetime.now(timezone.utc).isoformat()

    # --- Template loading ---
    shape = load_template_shape(args.template, args.material)
    clk_ns_fallback = 4.0  # filled in from pipeline below

    if args.random_dt:
        dt_mode_str = f"random ΔT in [{args.dt_min}, {args.dt_max}] ns"
    else:
        dt_mode_str = f"fixed grid {args.dt_grid}"
    print(
        f"[setup] template   : {shape.source_path}\n"
        f"[setup] material   : {shape.material}  "
        f"τ_r={shape.tau_r_ns:.4f} ns  τ_f={shape.tau_f_ns:.4f} ns  "
        f"T_max={shape.T_max:.4f}  t_peak_offset={shape.t_peak_offset_ns:.3f} ns\n"
        f"[setup] dt mode    : {dt_mode_str}\n"
        f"[setup] ratio_grid : {args.ratio_grid}\n"
        f"[setup] n_per_config={args.n_per_config}  "
        f"noise_model={args.noise_model}  seed={args.seed}",
        flush=True,
    )

    # --- Output directory ---
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if not args.force:
        existing = list(out_dir.iterdir())
        if existing:
            raise SystemExit(
                f"[ERROR] --out-dir {out_dir} is non-empty "
                f"(use --force to overwrite)"
            )

    # --- Multi-input EVIO discovery (mirrors fit_pulse_template.py lines 661–684) ---
    all_files: List[str] = []
    seen: set = set()
    for inp in args.evio_paths:
        for fp in C.discover_split_files(inp):
            if fp not in seen:
                seen.add(fp)
                all_files.append(fp)
    if not all_files:
        raise SystemExit(
            f"[ERROR] no EVIO splits found for inputs: {args.evio_paths}"
        )
    print(f"[setup] {len(all_files)} EVIO split(s) across "
          f"{len(args.evio_paths)} input arg(s)", flush=True)

    # setup_pipeline does its own discovery; override p.evio_files with multi-input list.
    p = C.setup_pipeline(
        evio_path    = args.evio_paths[0],
        max_events   = args.max_events,
        daq_config   = args.daq_config,
        hc_map_file  = args.hc_map_file,
    )
    p.evio_files = all_files

    clk_mhz = float(p.cfg.wave_cfg.clk_mhz)
    clk_ns  = (1000.0 / clk_mhz) if clk_mhz > 0 else clk_ns_fallback
    print(f"[setup] clk_ns={clk_ns}  ({clk_mhz} MHz)", flush=True)

    # --- Gather run numbers from EVIO filenames ---
    runs = sorted({C.extract_run_number(f) for f in all_files if C.extract_run_number(f) >= 0})

    # --- Build grid ---
    ratio_grid = args.ratio_grid

    # We don't know n_samples until we start reading EVIO, but we can get it
    # from the DAQ config's wave window, or just use a safe 200.
    # We'll patch grid_point.n_samples on first base event.
    # For now use a placeholder; append() will write to _waveforms which is
    # pre-allocated to n_per_config × n_samples.  We correct below.
    # Actually, we need to allocate before reading — peek at the wave config.
    try:
        n_samples_cfg = int(p.cfg.wave_cfg.window_size)
        if n_samples_cfg <= 0:
            n_samples_cfg = 200
    except Exception:
        n_samples_cfg = 200

    if args.random_dt:
        # Random-ΔT mode: one cell per ratio; ΔT drawn per event from [dt_min, dt_max].
        # gp.dt_ns is NaN (sentinel — not used for injection, only event_dt_ns matters).
        # dt_max_ns for the headroom cut uses args.dt_max.
        dt_max_ns = args.dt_max
        grid_points = [
            GridPoint(dt_ns=float("nan"), ratio=r,
                      n_per_config=args.n_per_config,
                      n_samples=n_samples_cfg)
            for r in sorted(set(ratio_grid))
        ]
        n_cells = len(grid_points)
        print(f"[setup] random-dt mode: {len(ratio_grid)} ratio cells  "
              f"dt_min={args.dt_min} dt_max={args.dt_max}  "
              f"target={args.n_per_config} per cell", flush=True)

        # RNG map keyed by ratio
        ss = np.random.SeedSequence(args.seed)
        child_seeds = ss.spawn(n_cells)
        rng_map_ratio: Dict[float, np.random.Generator] = {}
        for gp, child in zip(grid_points, child_seeds):
            rng_map_ratio[gp.ratio] = np.random.default_rng(child)
    else:
        dt_grid   = args.dt_grid
        dt_max_ns = max(dt_grid)
        grid_points = []
        for dt in sorted(set(dt_grid)):
            for r in sorted(set(ratio_grid)):
                grid_points.append(GridPoint(
                    dt_ns=dt, ratio=r,
                    n_per_config=args.n_per_config,
                    n_samples=n_samples_cfg,
                ))
        n_cells = len(grid_points)
        print(f"[setup] grid: {len(dt_grid)} dt × {len(ratio_grid)} ratio = "
              f"{n_cells} cells  target={args.n_per_config} per cell", flush=True)

        # RNG map keyed by (dt, ratio)
        ss = np.random.SeedSequence(args.seed)
        child_seeds = ss.spawn(n_cells)
        rng_map_fixed: Dict[Tuple[float, float], np.random.Generator] = {}
        for gp, child in zip(grid_points, child_seeds):
            rng_map_fixed[gp.key] = np.random.default_rng(child)

    # --- Scan ---
    meta_common: Dict = {
        "clk_ns":              clk_ns,
        "noise_model":         args.noise_model,
        "seed":                args.seed,
        "generated_utc":       generated_utc,
        "n_per_config":        args.n_per_config,
        "cut_height_min":      args.height_min,
        "cut_height_max":      args.height_max,
        "cut_height_rms_mult": args.height_rms_mult,
        "cut_t0_min_ns":       args.t0_min,
        "cut_chi2_per_dof_max":args.chi2_max,
        "residual_veto_sigma": args.residual_veto_sigma,
        "runs":                runs,
        "n_evio_splits":       len(all_files),
        "n_samples":           n_samples_cfg,
        "random_dt":           args.random_dt,
    }
    if args.random_dt:
        meta_common["dt_min_ns"] = args.dt_min
        meta_common["dt_max_ns"] = args.dt_max

    base_iter = iterate_base_events(p, args, shape, clk_ns, dt_max_ns)

    try:
        if args.random_dt:
            n_bases_done = scan_grid_random_dt(
                base_iter, grid_points, args.n_per_config,
                shape, clk_ns, args.noise_model, rng_map_ratio,
                args.dt_min, args.dt_max,
                t0_wall, args.progress_every,
            )
        else:
            n_bases_done = scan_grid(
                base_iter, grid_points, args.n_per_config,
                shape, clk_ns, args.noise_model, rng_map_fixed,
                t0_wall, args.progress_every,
            )
    except KeyboardInterrupt:
        print("\n[interrupted — partial output]", flush=True)
        n_bases_done = sum(gp.n_filled for gp in grid_points) // max(1, n_cells)

    # Exhaust the iterator to let it finalize rejection counts
    try:
        for _ in base_iter:
            pass
    except Exception:
        pass

    # --- Warn on under-filled cells ---
    underfilled = [gp for gp in grid_points if gp.n_filled < args.n_per_config]
    if underfilled:
        print(
            f"[WARN] {len(underfilled)} grid cell(s) under-filled:",
            file=sys.stderr, flush=True,
        )
        for gp in underfilled:
            print(f"  dt={gp.dt_ns} ratio={gp.ratio} n={gp.n_filled}/{args.n_per_config}",
                  file=sys.stderr)

    # --- Write output ---
    file_entries: List[Dict] = []
    n_written_total = 0

    for gp in grid_points:
        if gp.n_filled == 0:
            continue
        # n_samples is self-correcting: GridPoint.append() reallocates _waveforms
        # on the first append if the actual sample count differs from the
        # DAQ-config estimate (see append() implementation above).
        npz_p, json_p = write_grid_point(
            out_dir, gp, shape, meta_common, args.metadata_in_npz,
            random_dt=args.random_dt,
        )
        validate_output(npz_p, wave_ana=None)
        entry: Dict = {
            "ratio":     gp.ratio,
            "path":      npz_p.name,
            "sidecar":   json_p.name,
            "n_written": gp.n_filled,
        }
        if args.random_dt:
            entry["dt_ns"] = None
        else:
            entry["dt_ns"] = gp.dt_ns
        file_entries.append(entry)
        n_written_total += gp.n_filled

    write_manifest(out_dir, shape, meta_common, grid_points, file_entries, cli_argv)

    # --- Summary ---
    elapsed = time.monotonic() - t0_wall
    rc = _rejection_counts
    print(
        f"[done] target={args.n_per_config}  kept={n_bases_done}  "
        f"wrote={n_written_total} synthetic events  "
        f"files={len(file_entries)}  elapsed={elapsed:.1f}s",
        flush=True,
    )
    print("[done] rejection histogram:", flush=True)
    for key, val in [
        ("total_candidates",     rc.get("n_total_candidates", "?")),
        ("phys_events_scanned",  rc.get("n_phys_events_scanned", "?")),
        ("bases_yielded",        rc.get("n_bases_yielded", "?")),
        ("multi_peak",           rc.get("n_reject_multi_peak", "?")),
        ("quality_flag",         rc.get("n_reject_quality", "?")),
        ("overflow",             rc.get("n_reject_overflow", "?")),
        ("height",               rc.get("n_reject_height", "?")),
        ("snr",                  rc.get("n_reject_snr", "?")),
        ("material",             rc.get("n_reject_material", "?")),
        ("fit_not_converged",    rc.get("n_reject_fit_converged", "?")),
        ("fit_t0",               rc.get("n_reject_fit_t0", "?")),
        ("fit_chi2",             rc.get("n_reject_fit_chi2", "?")),
        ("residual_veto",        rc.get("n_reject_residual_veto", "?")),
        ("headroom",             rc.get("n_reject_headroom", "?")),
        ("saturation",           rc.get("n_reject_saturation", "?")),
    ]:
        print(f"  {key:<26s} {val}", flush=True)

    if n_written_total == 0:
        print("[WARN] no events written — EVIO had zero usable base events.",
              file=sys.stderr)
        return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
