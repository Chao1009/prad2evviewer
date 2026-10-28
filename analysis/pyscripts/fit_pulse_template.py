#!/usr/bin/env python3
"""
fit_pulse_template.py — per-channel pulse-shape characterisation for the
PRad-II FADC250 readout.

For every channel that carries waveform samples in the input EVIO file(s)
we collect "clean" isolated pulses (single peak, no pile-up flag, no
overflow, peak well inside the buffer) and fit each one with a two-tau
pulse model:

    p(t) = A · (1 − exp(−(t−t0)/τ_r)) · exp(−(t−t0)/τ_f)    for t ≥ t0
           0                                                 otherwise

Per-pulse fits give a population of (τ_r, τ_f, A, t0, χ²/dof) values per
channel.  We summarise them with median + MAD (robust against the noisy
tail of the population) and emit one entry per channel into a JSON file.

The output feeds the pile-up deconvolution (prad2dec PulseTemplateStore);
see docs/technical_notes/waveform_analysis/wave_analysis.md.

Plotting
--------
With `--plot-dir <dir>` set, we:
  * keep a small cache of raw pulses per channel during fitting,
  * auto-pick the N best-fitting and N worst-fitting channels (criterion:
    median χ²/dof; "good" = passed the good_fit gate, "bad" = enough
    pulses but failed the gate),
  * write one diagnostic PNG per picked channel (raw pulses + median-fit
    overlay) under <plot-dir>/{good,bad,explicit}/<name>.png,
  * write summary PNGs with the global τ_r / τ_f / χ² distributions
    (per-channel medians and pooled per-pulse fits).

Use `--plot-channels W001,W002,...` to force-plot specific channels in
addition to the auto-picked ones.

Module-multiplicity filter
--------------------------
With `--min-modules-per-event N` and/or `--max-modules-per-event N` set, we do
a preliminary pass on every channel to count how many have a qualifying pulse
(defined by `--module-count-height-min` and `--module-count-require-clean`),
then reject the entire event if the count is outside the requested range.
Motivation: isolating candidate multi-module discharge events, where a
discharge between two adjacent HyCal modules produces exactly two firing
channels on an otherwise unrelated trigger. Typical use:
`--min-modules-per-event 2 --max-modules-per-event 2` to select the
two-channel topology.

Usage
-----
    python fit_pulse_template.py <evio_path> [<evio_path> ...] -o out.json
                                 [--max-events N]
                                 [--max-pulses-per-channel N]
                                 [--min-pulses N]
                                 [--chi2-max X]
                                 [--height-min ADC] [--height-rms-mult X]
                                 [--model two_tau|two_tau_p]
                                 [--model-err-floor F]
                                 [--pre-samples N] [--post-samples N]
                                 [--channels W001,W002,…]
                                 [--plot-dir plots/]
                                 [--n-plot-good N] [--n-plot-bad N]
                                 [--plot-channels W001,W002,…]
                                 [--no-summary-plot]
                                 [--min-modules-per-event N]
                                 [--max-modules-per-event N]
                                 [--module-count-height-min ADC]
                                 [--module-count-require-clean | --module-count-allow-dirty]
                                 [--progress-every N]
                                 [--daq-config FILE] [--hc-map-file FILE]

Each `<evio_path>` accepts the same shapes as the other analysis scripts —
glob, directory, or single split file.  Multiple paths are concatenated.
Run `--help` for option details.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

import _common as C
from _common import dec  # prad2py.dec re-export


# Number of raw pulses we keep per channel for the diagnostic plots.
# 20 × 50 samples × 8 bytes × ~2000 channels ≈ 16 MB — negligible.
PULSE_CACHE = 20
PULSE_CACHE_PER_BIN = 30  # per-amplitude-bin waveform cache (Test A diagnostic)

# Amplitude bins for Test A: (label, lo, hi).  Bins are inclusive on lo,
# exclusive on hi, using -inf / +inf for open ends.
AMP_BINS = [
    ("lt_300",     float("-inf"), 300.0),
    ("300_to_500", 300.0,         500.0),
    ("gt_500",     500.0,         float("inf")),
]

# Stable colour assignments for the per-type summary panels.  Types that
# don't appear in this map fall through to matplotlib's default cycle.
TYPE_COLORS = {"PbGlass": "C0", "PbWO4": "C1",
               "LMS":     "C2", "Veto":  "C3", "Unknown": "0.5"}


def _resolve_trigger_bits(spec: str, bits_file: Path) -> int:
    """Resolve a comma-separated trigger spec ('SSP0,LMS' or '8,24') to a
    bitmask (int). Empty spec returns 0. Names are looked up in
    trigger_bits.json's 'trigger_bits' list."""
    if not spec.strip():
        return 0
    # Build name -> bit lookup
    with open(bits_file, "r", encoding="utf-8") as f:
        bits_json = json.load(f)
    name_to_bit = {}
    for entry in bits_json.get("trigger_bits", []):
        name_to_bit[entry["name"]] = entry["bit"]
    mask = 0
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if token in name_to_bit:
            mask |= (1 << name_to_bit[token])
        else:
            try:
                bit_num = int(token, 0)  # supports "8", "0x8", etc.
                mask |= (1 << bit_num)
            except ValueError:
                raise SystemExit(
                    f"[ERROR] unknown trigger name or bit: '{token}'\n"
                    f"        valid names: {sorted(name_to_bit.keys())}\n"
                    f"        or use bit numbers (e.g., '8', '0x10')")
    return mask


def _resolve_trigger_event_types(spec: str, bits_file: Path) -> set:
    """Resolve a comma-separated trigger event-type spec ('SSP_RawSum,LMS')
    to a set of integer type codes. Empty spec returns empty set (= no filter).
    Names are looked up in trigger_bits.json's 'trigger_type' list."""
    if not spec.strip():
        return set()
    with open(bits_file, "r", encoding="utf-8") as f:
        bits_json = json.load(f)
    name_to_type = {}
    for entry in bits_json.get("trigger_type", []):
        name_to_type[entry["name"]] = int(entry["type"], 0)
    result = set()
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if token in name_to_type:
            result.add(name_to_type[token])
        else:
            try:
                result.add(int(token, 0))
            except ValueError:
                raise SystemExit(
                    f"[ERROR] unknown trigger event-type name: '{token}'\n"
                    f"        valid names: {sorted(name_to_type.keys())}\n"
                    f"        or use hex/int codes (e.g., '0x29', '41')")
    return result


def _replace_nan_recursive(obj):
    """Replace any NaN/Inf floats in a nested dict/list with None.

    JSON standard forbids NaN and Inf literals; nlohmann/json (used by the
    C++ PulseTemplateStore) rejects them with a parse error.  Channels with
    no converged fits are skipped before reaching this point (Option A), but
    this sanitiser is kept as a belt-and-braces fallback (Option B) so that
    any pathological NaN that slips through causes a clear ValueError from
    json.dump(allow_nan=False) rather than silently producing broken JSON.
    """
    import math
    if isinstance(obj, dict):
        return {k: _replace_nan_recursive(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_replace_nan_recursive(v) for v in obj]
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    return obj


# ---- Pulse model ----
# Fitting itself runs in C++ via dec.WaveAnalyzer.fit_pulse_shape() —
# three-parameter Levenberg-Marquardt on the unit-amplitude two-tau model
# (T(t; t0, τ_r, τ_f) / T_max), with the per-pulse peak-height
# normalisation that decouples shape from amplitude.  See
# prad2dec/src/WaveAnalyzer.cpp::FitPulseShape for the implementation.

def two_tau_p_unit(t: np.ndarray, t0: float, tau_r: float, tau_f: float,
                   p: float) -> np.ndarray:
    """[1 - exp(-(t-t0)/τ_r)]^p · exp(-(t-t0)/τ_f), peak normalised to 1;
    p = 1 is the plain two-tau model.  For plot overlays only; the C++
    FitPulseShape / FitPulseShapeTwoTauP are the real fitters."""
    out = np.zeros_like(t, dtype=np.float64)
    mask = t > t0
    if not mask.any():
        return out
    dt = t[mask] - t0
    u = 1.0 - np.exp(-dt / tau_r)
    u = np.clip(u, 0.0, None)
    raw = (u ** p) * np.exp(-dt / tau_f)
    denom = tau_r + p * tau_f
    u_pk = p * tau_f / denom
    t_max = (u_pk ** p) * ((tau_r / denom) ** (tau_r / tau_f))
    out[mask] = raw / t_max
    return out


# ---- Per-channel accumulator ----

@dataclass
class ChannelStats:
    name: str
    channel_id: str            # roc_<tag>_<slot>_<channel>
    module_type: str = "Unknown"   # mod.type.name from HyCalSystem; "Unknown"
                                   # for non-HyCal channels (tagger TDC, etc).
    n_total: int = 0 # total in the evio pool
    n_trigger: int = 0 # triggered by total sum or LMS
    n_good: int = 0 # good pulse
    n_after_cut: int = 0 # after pulse height
    n_good_fit: int = 0 #good fits 
    n_good_fit_time_cut: int = 0 #good fits after time cut
    n_attempted: int = 0
    n_used: int = 0
    tau_r:    List[float] = field(default_factory=list)
    tau_f:    List[float] = field(default_factory=list)
    t0:       List[float] = field(default_factory=list)
    p_list:   List[float] = field(default_factory=list)   # populated only when --model=two_tau_p
    chi2:     List[float] = field(default_factory=list)
    peak_amp: List[float] = field(default_factory=list)   # original ADC peak height per pulse
    ped_mean_sum: float = 0.0
    ped_rms_sum:  float = 0.0
    ped_n:        int   = 0
    # First few pedsub pulses, kept for the diagnostic plot only.  Stored
    # *un-normalised* — plotting code divides by the matching peak_amp on
    # the fly so users see the raw shapes alongside the normalised stack.
    sample_pulses: List[np.ndarray] = field(default_factory=list)
    sample_peak_amps: List[float]   = field(default_factory=list)
    # Per-amplitude-bin caches for Test A (background-vs-signal shape check).
    # Keys are amp-bin labels ("lt_300", "300_to_500", "gt_500"); values are
    # lists of pedsub pulse arrays up to PULSE_CACHE_PER_BIN per bin.
    # Populated only when plotting is enabled AND the amplitude bin isn't full.
    sample_pulses_by_amp: Dict[str, List[np.ndarray]] = field(default_factory=dict)
    sample_peak_amps_by_amp: Dict[str, List[float]]   = field(default_factory=dict)


def _median_mad(values: List[float]) -> Tuple[float, float]:
    if not values:
        return (float("nan"), float("nan"))
    arr = np.asarray(values, dtype=np.float64)
    med = float(np.median(arr))
    mad = float(np.median(np.abs(arr - med)))
    return med, mad


def aggregate_by_type(summaries: List[Dict], min_pulses: int) -> Dict[str, Dict]:
    """For each module type present in `summaries`, take the channels with
    enough pulses and compute median ± MAD across their per-channel
    medians for every fit parameter.  These per-type aggregates are the
    templates the C++ PulseTemplateStore actually feeds to the deconv —
    one shape per category (PbGlass / PbWO4 / LMS / Veto), shared by
    every channel of that type.  Type strings come from
    HyCalSystem.Module.type (sourced from hycal_map.json), so
    Veto/LMS/PbWO4/PbGlass are whatever the JSON spelled them."""
    by_type: Dict[str, Dict] = {}
    types_present = sorted({s["module_type"] for s in summaries})
    for t in types_present:
        chans = [s for s in summaries
                 if s["module_type"] == t
                 and s["n_pulses_used"] >= min_pulses]
        if not chans:
            continue
        rec: Dict = {
            "n_channels":        len(chans),
            "n_channels_good":   int(sum(s["good_fit"] for s in chans)),
            "n_pulses_total":    int(sum(s["n_pulses_used"] for s in chans)),
        }
        # Aggregate every numeric "{median, mad}" parameter that is present
        # on at least one channel.  Channel medians become the per-type
        # population — we report median+MAD over them.
        for key in ("t0_ns", "tau_r_ns", "tau_f_ns", "p", "peak_amp_adc"):
            vals = [s[key]["median"] for s in chans if key in s]
            if not vals:
                continue
            med, mad = _median_mad(vals)
            rec[key] = {"median": med, "mad": mad}
        chi2_vals = [s["chi2_per_dof"]["median"] for s in chans]
        rec["chi2_per_dof"] = {
            "median": float(np.median(chi2_vals)),
            "mean":   float(np.mean(chi2_vals)),
        }
        by_type[t] = rec
    return by_type


def finalize_channel(s: ChannelStats, min_pulses: int, chi2_max: float
                     ) -> Dict:
    tr_med, tr_mad = _median_mad(s.tau_r)
    tf_med, tf_mad = _median_mad(s.tau_f)
    t0_med, t0_mad = _median_mad(s.t0)
    pk_med, pk_mad = _median_mad(s.peak_amp)
    chi2_med = float(np.median(s.chi2)) if s.chi2 else float("nan")
    chi2_mean = float(np.mean(s.chi2)) if s.chi2 else float("nan")
    good = bool(s.n_used >= min_pulses and chi2_med < chi2_max)
    ped_mean = s.ped_mean_sum / s.ped_n if s.ped_n else float("nan")
    ped_rms  = s.ped_rms_sum  / s.ped_n if s.ped_n else float("nan")
    rec = {
        "channel_id": s.channel_id,
        "module_type": s.module_type,
        "n_pulses_attempted": s.n_attempted,
        "n_pulses_used":      s.n_used,
        "n_good":             s.n_good,
        "n_after_cut":        s.n_after_cut,
        "n_good_fit":           s.n_good_fit,
        "n_good_fit_time_cut":  s.n_good_fit_time_cut,
        "tau_r_ns": {"median": tr_med, "mad": tr_mad},
        "tau_f_ns": {"median": tf_med, "mad": tf_mad},
        "t0_ns":    {"median": t0_med, "mad": t0_mad},
        "peak_amp_adc": {"median": pk_med, "mad": pk_mad},
        "chi2_per_dof": {"median": chi2_med, "mean": chi2_mean},
        "ped": {"mean": float(ped_mean), "rms": float(ped_rms)},
        "good_fit": good,
    }
    if s.p_list:                     # only present for --model=two_tau_p
        p_med, p_mad = _median_mad(s.p_list)
        rec["p"] = {"median": p_med, "mad": p_mad}
    return rec


# ---- Plotting ----

def plot_channel(plt, st: ChannelStats, clk_ns: float, pre: int, post: int,
                 out_png: Path, label: str = "") -> None:
    """Two-panel diagnostic — top: raw pulses + raw-scale median fit;
    bottom: normalised stack + unit-amplitude median-shape fit (the view
    the LM actually sees)."""

    fig, (ax_raw, ax_norm) = plt.subplots(
        2, 1, figsize=(8, 8),
        sharex=True,
        gridspec_kw=dict(height_ratios=[1, 1], hspace=0.08))

    has_fit = bool(st.tau_r)

    if has_fit:
        med_t0    = float(np.median(st.t0))
        med_tau_r = float(np.median(st.tau_r))
        med_tau_f = float(np.median(st.tau_f))
        med_pk    = float(np.median(st.peak_amp))
        t_dense   = np.linspace(0, (pre + post) * clk_ns, 400)
        med_p     = float(np.median(st.p_list)) if st.p_list else None
        unit_curve = two_tau_p_unit(t_dense, med_t0, med_tau_r, med_tau_f,
                                    1.0 if med_p is None else med_p)

        # Top — raw scale.
        for s in st.sample_pulses:
            t = np.arange(len(s)) * clk_ns
            ax_raw.plot(t, s, color="C0", alpha=0.25, lw=0.7)
        ax_raw.plot(t_dense, unit_curve * med_pk,
                    color="C3", lw=2.0, label="median shape × median peak")

        # Bottom — normalised stack.
        for s, pk in zip(st.sample_pulses, st.sample_peak_amps):
            if pk <= 0:
                continue
            t = np.arange(len(s)) * clk_ns
            ax_norm.plot(t, s / pk, color="C0", alpha=0.25, lw=0.7)
        ax_norm.plot(t_dense, unit_curve, color="C3", lw=2.0,
                     label="unit-amplitude fit")
        ax_norm.axhline(1.0, color="0.6", lw=0.6, ls="--")

        chi2_med = float(np.median(st.chi2)) if st.chi2 else float("nan")
        p_str = f"  p={med_p:.2f}" if med_p is not None else ""
        title = (f"{label + ' ' if label else ''}{st.name}  "
                 f"({st.channel_id})\n"
                 f"n_used={st.n_used}  "
                 f"τ_r={med_tau_r:.2f} ns  τ_f={med_tau_f:.2f} ns{p_str}  "
                 f"peak={med_pk:.0f} ADC  χ²/dof med={chi2_med:.2f}")
    else:
        for s in st.sample_pulses:
            t = np.arange(len(s)) * clk_ns
            ax_raw.plot(t, s, color="C0", alpha=0.25, lw=0.7)
        title = f"{st.name}: no successful fits"

    ax_raw.set_ylabel("pedsub ADC")
    ax_raw.grid(True, alpha=0.3)
    if has_fit:
        ax_raw.legend(loc="best", fontsize=8)
    ax_raw.set_title(title)

    ax_norm.set_xlabel("time within window (ns)")
    ax_norm.set_ylabel("pedsub / peak  (unit)")
    ax_norm.grid(True, alpha=0.3)
    if has_fit:
        ax_norm.legend(loc="best", fontsize=8)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)


def _split_hist(ax, arrays: Dict[str, np.ndarray], bins, xlabel: str,
                title: str, ylabel: str, count_sym: str,
                log: bool = False) -> None:
    """One histogram per module type, overlaid; legend shows per-type counts."""
    if log:
        ax.set_xscale("log")
    for t, arr in arrays.items():
        if arr.size == 0:
            continue
        ax.hist(arr, bins=bins,
                color=TYPE_COLORS.get(t),
                alpha=0.65, label=f"{t} ({count_sym}={arr.size})")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)


def _wide_bottom_axis(fig, axes):
    """Merge the bottom row of a 3×3 subplot grid into one full-width axis."""
    gs = axes[2, 0].get_gridspec()
    for a in axes[2, :]:
        a.remove()
    return fig.add_subplot(gs[2, :])


def plot_summary(plt, summaries: List[Dict], min_pulses: int, chi2_max: float,
                 out_png: Path) -> None:
    """Global histograms of every fit parameter, split by module type
    (PbGlass / PbWO4 / LMS / Veto — whatever HyCalSystem reported per
    channel).  Bottom row: τ_r vs τ_f scatter coloured by χ²/dof."""
    have = [s for s in summaries if s["n_pulses_used"] >= min_pulses]
    if not have:
        print("[plot] no channels with enough pulses for summary", flush=True)
        return

    types_arr = np.array([s["module_type"] for s in have])
    # Order types by population (largest first) for legend stability.
    types_present = sorted(set(types_arr.tolist()),
                           key=lambda t: -int((types_arr == t).sum()))
    type_masks = {t: (types_arr == t) for t in types_present}

    t0       = np.array([s["t0_ns"]["median"]        for s in have])
    tau_r    = np.array([s["tau_r_ns"]["median"]     for s in have])
    tau_f    = np.array([s["tau_f_ns"]["median"]     for s in have])
    peak_amp = np.array([s["peak_amp_adc"]["median"] for s in have])
    chi2     = np.array([s["chi2_per_dof"]["median"] for s in have])
    good     = np.array([s["good_fit"]               for s in have])
    has_p    = all("p" in s for s in have)
    p_arr    = (np.array([s["p"]["median"] for s in have])
                if has_p else None)

    def _hist(ax, x, bins, xlabel, title, log=False):
        _split_hist(ax, {t: x[type_masks[t]] for t in types_present},
                    bins, xlabel, title, "channels", "N", log)

    fig, axes = plt.subplots(3, 3, figsize=(15, 12))

    _hist(axes[0, 0], t0,    60, "t0 (ns) — channel median",      "pulse onset t0")
    _hist(axes[0, 1], tau_r, 60, "τ_r (ns) — channel median",     "rise time τ_r")
    _hist(axes[0, 2], tau_f, 60, "τ_f (ns) — channel median",     "fall time τ_f")

    if has_p:
        _hist(axes[1, 0], p_arr, 60, "p — channel median",        "rise-edge exponent p")
    else:
        axes[1, 0].set_visible(False)

    pk_bins = np.logspace(np.log10(max(peak_amp.min(), 1.0)),
                          np.log10(max(peak_amp.max(), 10.0)), 60)
    _hist(axes[1, 1], peak_amp, pk_bins, "peak amplitude (ADC) — channel median",
          "peak amplitude", log=True)

    chi2_bins = np.logspace(np.log10(max(chi2.min(), 1e-2)),
                            np.log10(max(chi2.max(), 10.0)), 60)
    ax = axes[1, 2]
    ax.set_xscale("log")
    ax.hist(chi2[good],  bins=chi2_bins, color="C2", alpha=0.7, label="good_fit")
    ax.hist(chi2[~good], bins=chi2_bins, color="C3", alpha=0.7, label="bad")
    ax.axvline(chi2_max, color="k", ls="--", lw=1, label=f"χ²_max={chi2_max}")
    ax.set_xlabel("χ²/dof — channel median")
    ax.set_ylabel("channels")
    ax.set_title(f"fit quality   (good={int(good.sum())}, "
                 f"bad={int((~good).sum())})")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)

    # Bottom row: τ_r vs τ_f scatter, full width.
    ax_sc = _wide_bottom_axis(fig, axes)
    sc = ax_sc.scatter(tau_r, tau_f, c=np.log10(np.clip(chi2, 1e-2, None)),
                       s=10, cmap="viridis")
    fig.colorbar(sc, ax=ax_sc, label="log₁₀(χ²/dof)")
    ax_sc.set_xlabel("τ_r (ns)")
    ax_sc.set_ylabel("τ_f (ns)")
    ax_sc.set_title("τ_r vs τ_f, coloured by χ²/dof")
    ax_sc.grid(True, alpha=0.3)

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def plot_summary_per_pulse(plt, stats: Dict[str, "ChannelStats"],
                           min_pulses: int, out_png: Path) -> None:
    """Same parameter grid as plot_summary, but every per-pulse fit is one
    histogram entry — pooled by module type.  Tells us whether shape really
    is type-uniform (tight per-pulse distributions) and resolves whether
    bound-hits in the per-channel-median view are genuinely active per-pulse
    or just dragging the median from the tail."""
    chans = [s for s in stats.values() if s.n_used >= min_pulses]
    if not chans:
        print("[plot] no channels for per-pulse summary", flush=True)
        return

    # Pool per-pulse arrays by type.  All lists are populated only for
    # converged fits, so they're all the same length per channel.
    pool: Dict[str, Dict[str, List[float]]] = {}
    for s in chans:
        d = pool.setdefault(s.module_type, {
            "t0": [], "tau_r": [], "tau_f": [], "p": [],
            "peak_amp": [], "chi2": [],
        })
        d["t0"].extend(s.t0)
        d["tau_r"].extend(s.tau_r)
        d["tau_f"].extend(s.tau_f)
        d["p"].extend(s.p_list)
        d["peak_amp"].extend(s.peak_amp)
        d["chi2"].extend(s.chi2)

    types_present = sorted(pool.keys(),
                           key=lambda t: -len(pool[t]["tau_r"]))
    has_p = any(len(pool[t]["p"]) > 0 for t in types_present)

    def _hist(ax, key, bins, xlabel, title, log=False):
        _split_hist(ax, {t: np.asarray(pool[t][key], dtype=np.float64)
                         for t in types_present},
                    bins, xlabel, title, "pulses", "n", log)

    # Auto-range bins from the pooled data so type-tails are visible.
    def _range_bins(key, n=80, log=False):
        flat = np.concatenate([np.asarray(pool[t][key]) for t in types_present
                               if pool[t][key]])
        lo, hi = float(np.nanmin(flat)), float(np.nanmax(flat))
        if log:
            lo = max(lo, 1e-2)
            hi = max(hi, lo * 10.0)
            return np.logspace(np.log10(lo), np.log10(hi), n)
        return np.linspace(lo, hi, n)

    fig, axes = plt.subplots(3, 3, figsize=(15, 12))

    _hist(axes[0, 0], "t0",    _range_bins("t0"),
          "t0 (ns) — per pulse",       "pulse onset t0")
    _hist(axes[0, 1], "tau_r", _range_bins("tau_r"),
          "τ_r (ns) — per pulse",      "rise time τ_r")
    _hist(axes[0, 2], "tau_f", _range_bins("tau_f"),
          "τ_f (ns) — per pulse",      "fall time τ_f")

    if has_p:
        _hist(axes[1, 0], "p", _range_bins("p"),
              "p — per pulse",         "rise-edge exponent p")
    else:
        axes[1, 0].set_visible(False)

    _hist(axes[1, 1], "peak_amp", _range_bins("peak_amp", log=True),
          "peak amplitude (ADC) — per pulse",
          "peak amplitude", log=True)

    _hist(axes[1, 2], "chi2", _range_bins("chi2", log=True),
          "χ²/dof — per pulse", "fit quality", log=True)

    # Bottom row: per-pulse τ_r vs τ_f scatter, full width, by type.
    ax_sc = _wide_bottom_axis(fig, axes)
    for t in types_present:
        tr = np.asarray(pool[t]["tau_r"])
        tf = np.asarray(pool[t]["tau_f"])
        if tr.size == 0:
            continue
        ax_sc.scatter(tr, tf, s=2, alpha=0.25,
                      color=TYPE_COLORS.get(t), label=f"{t} (n={tr.size})")
    ax_sc.set_xlabel("τ_r (ns) — per pulse")
    ax_sc.set_ylabel("τ_f (ns) — per pulse")
    ax_sc.set_title("τ_r vs τ_f, per pulse, coloured by type")
    ax_sc.legend(loc="best", fontsize=8, markerscale=4)
    ax_sc.grid(True, alpha=0.3)

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=120)
    plt.close(fig)


def select_plot_targets(stats: Dict[str, ChannelStats], min_pulses: int,
                        chi2_max: float, n_good: int, n_bad: int
                        ) -> Tuple[List[str], List[str]]:
    """Auto-pick the N best- and N worst-fitting channels by median χ²/dof.
    "Good" = enough pulses AND median χ² < chi2_max (then sort ascending).
    "Bad" = enough pulses but median χ² >= chi2_max (then sort descending)."""
    good, bad = [], []
    for name, s in stats.items():
        if s.n_used < min_pulses or not s.chi2:
            continue
        med = float(np.median(s.chi2))
        (good if med < chi2_max else bad).append((med, name))
    good.sort()
    bad.sort(reverse=True)
    return [n for _, n in good[:n_good]], [n for _, n in bad[:n_bad]]


# ---- Main ----

def main() -> int:
    # Handle --list-triggers early, before argparse validates required args.
    if "--list-triggers" in sys.argv:
        # Resolve trigger_bits.json path
        bits_file_arg = None
        for i, arg in enumerate(sys.argv):
            if arg == "--trigger-bits-file" and i + 1 < len(sys.argv):
                bits_file_arg = sys.argv[i + 1]
                break
        if bits_file_arg:
            bits_file = Path(bits_file_arg)
        elif os.environ.get("PRAD2_DATABASE_DIR"):
            bits_file = Path(os.environ["PRAD2_DATABASE_DIR"]) / "trigger_bits.json"
        else:
            bits_file = Path(__file__).resolve().parents[2] / "database" / "trigger_bits.json"

        if not bits_file.is_file():
            print(f"[ERROR] trigger_bits.json not found at {bits_file}", file=sys.stderr)
            return 1

        with open(bits_file, "r", encoding="utf-8") as f:
            bits_json = json.load(f)

        print(f"Trigger bits (from {bits_file}):")
        print(f"  {'name':10s}  {'bit':>4s}  label")
        print(f"  {'-'*10}  {'-'*4}  {'-'*30}")
        for entry in bits_json.get("trigger_bits", []):
            print(f"  {entry['name']:10s}  {entry['bit']:>4d}  {entry.get('label', '')}")

        trigger_types = bits_json.get("trigger_type", [])
        if trigger_types:
            print()
            print("Trigger event types:")
            print(f"  {'name':16s}  {'type':>6s}  primary_bit  label")
            print(f"  {'-'*16}  {'-'*6}  {'-'*11}  {'-'*30}")
            for entry in trigger_types:
                print(f"  {entry['name']:16s}  {entry['type']:>6s}  "
                      f"{entry.get('primary_bit', ''):>11}  {entry.get('label', '')}")
        return 0

    # Resolve trigger bits file for help-text listing.  Same resolution order
    # as the runtime code below.
    _default_bits_file = None
    if os.environ.get("PRAD2_DATABASE_DIR"):
        _cand = Path(os.environ["PRAD2_DATABASE_DIR"]) / "trigger_bits.json"
        if _cand.is_file():
            _default_bits_file = _cand
    if _default_bits_file is None:
        _cand = Path(__file__).resolve().parents[2] / "database" / "trigger_bits.json"
        if _cand.is_file():
            _default_bits_file = _cand

    _known_trigger_names = "(unavailable — trigger_bits.json not found)"
    if _default_bits_file is not None:
        try:
            with open(_default_bits_file, "r", encoding="utf-8") as f:
                _known_trigger_names = ", ".join(
                    e["name"] for e in json.load(f).get("trigger_bits", []))
        except Exception:
            pass

    _known_event_type_names = "(unavailable)"
    if _default_bits_file is not None:
        try:
            with open(_default_bits_file, "r", encoding="utf-8") as f:
                _known_event_type_names = ", ".join(
                    e["name"] for e in json.load(f).get("trigger_type", []))
        except Exception:
            pass

    ap = argparse.ArgumentParser(
        description="Per-channel pulse-shape fit on FADC250 waveforms.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("evio_paths", nargs="+",
                    help="One or more EVIO inputs (each may be a glob, "
                         "directory, or single split file).  Splits found "
                         "for every input are concatenated.")
    ap.add_argument("-o", "--out", required=True,
                    help="Output JSON path.")
    ap.add_argument("--max-events", type=int, default=0,
                    help="Stop after N raw physics sub-events across all "
                         "files (0 = process everything).")
    ap.add_argument("--max-pulses-per-channel", type=int, default=500,
                    help="Stop fitting a channel once it has this many "
                         "accepted pulses (0 = unlimited, no cap).")
    ap.add_argument("--min-pulses", type=int, default=50,
                    help="Channel needs at least this many converged fits "
                         "to be marked good_fit=true (also gates summary "
                         "plot inclusion).")
    ap.add_argument("--chi2-max", type=float, default=3.0,
                    help="Median chi2/dof threshold for good_fit=true.")
    ap.add_argument("--height-min", type=float, default=20.0,
                    help="Minimum peak height (ADC, pedsub) for fitting.")
    ap.add_argument("--height-rms-mult", type=float, default=10.0,
                    help="Also require peak height >= this × ped.rms.")
    ap.add_argument("--t0-min", type=float, default=float("-inf"),
                    help="Post-fit filter: reject pulses with fit t0_ns < this "
                         "value (ns).  Default -inf (no lower cut).  Useful for "
                         "excluding early-arriving background pulses that fit "
                         "with atypically small t0.")
    ap.add_argument("--t0-max", type=float, default=float("inf"),
                    help="Post-fit filter: reject pulses with fit t0_ns > this "
                         "value (ns).  Default +inf (no upper cut).")
    ap.add_argument("--model", choices=["two_tau", "two_tau_p"],
                    default="two_tau",
                    help="Pulse-shape model.  'two_tau': "
                         "(1 - exp(-(t-t0)/τ_r)) · exp(-(t-t0)/τ_f).  "
                         "'two_tau_p': adds a rise-edge exponent so the "
                         "rise *shape* is free, not just the timescale: "
                         "[1 - exp(-(t-t0)/τ_r)]^p · exp(-(t-t0)/τ_f).  "
                         "Use two_tau_p when the standard model's onset "
                         "is too sharp on the data (high-amp PMT pulses "
                         "in particular).")
    ap.add_argument("--model-err-floor", type=float, default=0.01,
                    help="Floor on per-sample sigma in the normalised fit, as "
                         "a fraction of unit peak.  Caps how much a small "
                         "model-vs-data mismatch (sub-sample alignment, slow "
                         "second-tail component) can inflate chi2/dof on "
                         "high-amplitude pulses.  0 disables the floor.")
    ap.add_argument("--pre-samples",  type=int, default=8,
                    help="Samples kept before the peak for the fit window.")
    ap.add_argument("--post-samples", type=int, default=40,
                    help="Samples kept after the peak for the fit window.")
    ap.add_argument("--channels", default="",
                    help="Comma-separated channel name filter "
                         "(e.g. W001,W002).  Empty = all channels.")
    ap.add_argument("--plot-dir", default="",
                    help="Directory for diagnostic PNGs.  Empty = no plots.")
    ap.add_argument("--n-plot-good", type=int, default=5,
                    help="Auto-pick this many best-fitting channels for "
                         "diagnostic plots (sorted ascending χ²).")
    ap.add_argument("--n-plot-bad", type=int, default=5,
                    help="Auto-pick this many worst-fitting channels for "
                         "diagnostic plots (sorted descending χ²).")
    ap.add_argument("--plot-channels", default="",
                    help="Comma-separated extra channels to plot regardless "
                         "of fit quality (e.g. W001,W002).")
    ap.add_argument("--no-summary-plot", action="store_true",
                    help="Skip the global τ_r/τ_f/χ² summary PNG.")
    ap.add_argument("--progress-every", type=int, default=2000,
                    help="Print a progress line every N physics sub-events "
                         "(file index, events, fits, rate, elapsed).")
    ap.add_argument("--daq-config", default="",
                    help="DAQ config (default: installed default).")
    ap.add_argument("--hc-map-file", default="",
                    help="HyCal modules map (default: database lookup).")
    ap.add_argument("--trigger-accept", default="",
                    help="Comma-separated trigger names/bits to accept. "
                         f"Known names: {_known_trigger_names}. "
                         "Numeric bits (decimal or 0x) also accepted. "
                         "Empty (default) = accept all events. "
                         "See --list-triggers for full details.")
    ap.add_argument("--trigger-reject", default="",
                    help="Comma-separated trigger names/bits to reject "
                         "(applied after --trigger-accept). "
                         f"Known names: {_known_trigger_names}. "
                         "See --list-triggers for full details.")
    ap.add_argument("--trigger-bits-file", default="",
                    help="Path to trigger_bits.json (default: try "
                         "$PRAD2_DATABASE_DIR/trigger_bits.json, then "
                         "<repo>/database/trigger_bits.json).")
    ap.add_argument("--list-triggers", action="store_true",
                    help="Print all known trigger names and bits from "
                         "trigger_bits.json, then exit.")
    ap.add_argument("--trigger-event-type", default="",
                    help="Comma-separated trigger event-type names to accept "
                         "(e.g. 'SSP_RawSum,LMS'). Names are resolved from "
                         "trigger_bits.json's trigger_type list. "
                         f"Known names: {_known_event_type_names}. "
                         "Filters on fadc_evt.info.trigger_type (a single "
                         "per-event value), NOT on the trigger bitmask. "
                         "Empty (default) = accept all event types. "
                         "See --list-triggers for full details.")
    ap.add_argument("--min-modules-per-event", type=int, default=0,
                    metavar="N",
                    help="Reject events where fewer than N channels have a "
                         "qualifying pulse. Default 0 = no lower bound. A "
                         "'qualifying pulse' is defined by "
                         "--module-count-height-min and "
                         "--module-count-require-clean. Useful for isolating "
                         "multi-channel events (e.g. --min-modules-per-event 2 "
                         "--max-modules-per-event 2 selects candidate "
                         "two-channel discharge events).")
    ap.add_argument("--max-modules-per-event", type=int, default=0,
                    metavar="N",
                    help="Reject events where more than N channels have a "
                         "qualifying pulse. Default 0 = no upper bound.")
    ap.add_argument("--module-count-height-min", type=float, default=-1.0,
                    metavar="FLOAT",
                    help="Minimum peak height (ADC, pedsub) for a channel to "
                         "count as 'fired' toward the module-multiplicity "
                         "filter. Default -1.0 = reuse --height-min. Set this "
                         "independently when you want the multiplicity gate to "
                         "use a different (typically lower) threshold than the "
                         "fit-quality gate.")
    _mc_grp = ap.add_mutually_exclusive_group()
    _mc_grp.add_argument("--module-count-require-clean",
                         dest="module_count_require_clean",
                         action="store_true", default=True,
                         help="Require pk.quality == 0 for a channel to count "
                              "as 'fired' (default). Excludes piled-up or "
                              "overflow channels from the multiplicity count.")
    _mc_grp.add_argument("--module-count-allow-dirty",
                         dest="module_count_require_clean",
                         action="store_false",
                         help="Include piled-up and overflow channels in the "
                              "multiplicity count.")
    args = ap.parse_args()

    if args.trigger_bits_file:
        bits_file = Path(args.trigger_bits_file)
    elif os.environ.get("PRAD2_DATABASE_DIR"):
        bits_file = Path(os.environ["PRAD2_DATABASE_DIR"]) / "trigger_bits.json"
    else:
        # Fall back to repo location relative to this script
        bits_file = Path(__file__).resolve().parents[2] / "database" / "trigger_bits.json"

    accept_mask = _resolve_trigger_bits(args.trigger_accept, bits_file)
    reject_mask = _resolve_trigger_bits(args.trigger_reject, bits_file)

    if accept_mask or reject_mask:
        print(f"[setup] trigger accept mask = 0x{accept_mask:x}", flush=True)
        print(f"[setup] trigger reject mask = 0x{reject_mask:x}", flush=True)

    event_type_accept = _resolve_trigger_event_types(args.trigger_event_type, bits_file)
    if event_type_accept:
        print(f"[setup] trigger event-type accept = "
              f"{{{', '.join(f'0x{t:02x}' for t in sorted(event_type_accept))}}} "
              f"({args.trigger_event_type})", flush=True)

    chan_filter   = {s.strip() for s in args.channels.split(",") if s.strip()}
    extra_plot    = {s.strip() for s in args.plot_channels.split(",") if s.strip()}
    plotting      = bool(args.plot_dir)
    plot_dir      = Path(args.plot_dir) if plotting else None

    # Discover EVIO files for every input arg, concatenate, dedupe order-
    # preserving so users see the same iteration order as the CLI list.
    all_files: List[str] = []
    seen = set()
    for inp in args.evio_paths:
        for fp in C.discover_split_files(inp):
            if fp not in seen:
                seen.add(fp)
                all_files.append(fp)
    if not all_files:
        raise SystemExit(f"[ERROR] no EVIO splits found for inputs: "
                         f"{args.evio_paths}")
    print(f"[setup] {len(all_files)} EVIO split(s) across "
          f"{len(args.evio_paths)} input arg(s)", flush=True)

    # setup_pipeline still does its own discovery, but we override the
    # resulting evio_files with the concatenated multi-input list below.
    p = C.setup_pipeline(
        evio_path=args.evio_paths[0],
        daq_config=args.daq_config,
        hc_map_file=args.hc_map_file,
    )
    p.evio_files = all_files

    clk_ns = float(p.cfg.wave_cfg.clk_ns)
    pre, post = args.pre_samples, args.post_samples
    with_p = args.model == "two_tau_p"
    fit_pulse = (dec.WaveAnalyzer.fit_pulse_shape_two_tau_p if with_p
                 else dec.WaveAnalyzer.fit_pulse_shape)

    stats: Dict[str, ChannelStats] = {}
    # cache: (roc_tag, slot, ch) -> (name, chan_id, module_type)
    name_cache: Dict[Tuple[int, int, int], Tuple[str, str, str]] = {}

    def _key_for(roc_tag: int, s: int, c: int) -> Tuple[str, str, str]:
        cached = name_cache.get((roc_tag, s, c))
        if cached is not None:
            return cached
        chan_id = f"{roc_tag}_{s}_{c}"
        name = chan_id
        mtype = "Unknown"
        crate = p.crate_map.get(roc_tag)
        if crate is not None:
            mod = p.hycal.module_by_daq(crate, s, c)
            if mod is not None:
                name  = mod.name
                # mod.type is the C++ HyCalSystem ModuleType enum, populated
                # from hycal_map.json's "t" field — single source of truth.
                mtype = mod.type.name
        name_cache[(roc_tag, s, c)] = (name, chan_id, mtype)
        return name, chan_id, mtype

    n_pulses_attempted = n_pulses_used = 0
    loop_stats = C.LoopStats()

    n_events_trigger_filtered = 0
    n_events_multiplicity_filtered = 0
    n_trigger = n_good = n_after_cut = n_good_fit = n_good_fit_time_cut = 0
    _mult_active = (args.min_modules_per_event > 0
                    or args.max_modules_per_event > 0)
    _eff_height_min = (args.module_count_height_min
                       if args.module_count_height_min >= 0
                       else args.height_min)

    def _emit_progress(st: C.LoopStats) -> None:
        elapsed = st.elapsed
        rate = st.n_phys / elapsed if elapsed > 0 else 0.0
        print(f"  [progress] file {st.n_files_open}/{len(p.evio_files)} "
              f"{Path(st.cur_file).name}  "
              f"phys={st.n_phys}  fit_attempted={n_pulses_attempted}  "
              f"fit_used={n_pulses_used}  "
              f"rate={rate:.0f} ev/s  elapsed={elapsed:.1f}s",
              flush=True)

    try:
        for fadc_evt, _ in C.iter_physics_events(
                p, loop_stats, with_ssp=False, max_events=args.max_events,
                accept=None, progress=_emit_progress,
                progress_every=args.progress_every):
            # Trigger filter (optional; default no-op)
            if accept_mask or reject_mask:
                trigger_bits = int(fadc_evt.info.trigger_bits)
                if accept_mask and not (trigger_bits & accept_mask):
                    n_events_trigger_filtered += 1
                    continue
                if reject_mask and (trigger_bits & reject_mask):
                    n_events_trigger_filtered += 1
                    continue

            # Trigger event-type filter (optional; default no-op)
            if event_type_accept:
                if int(fadc_evt.info.trigger_type) not in event_type_accept:
                    n_events_trigger_filtered += 1
                    continue

            # Module-multiplicity filter (optional; default no-op): count
            # channels with a qualifying pulse before the full fit loop.
            if _mult_active:
                n_fired = 0
                for _, _, _, _cd in C.iter_fadc_channels(fadc_evt):
                    _ped, _rms, _peaks = p.wave_ana.analyze(
                        np.asarray(_cd.samples, dtype=np.uint16))
                    if len(_peaks) < 1:
                        continue
                    _pk = _peaks[0]
                    if args.module_count_require_clean and (
                            _pk.quality != 0 or _pk.overflow):
                        continue
                    if _pk.height < _eff_height_min:
                        continue
                    n_fired += 1
                    if (args.max_modules_per_event > 0
                            and n_fired > args.max_modules_per_event):
                        break
                if (args.min_modules_per_event > 0
                        and n_fired < args.min_modules_per_event):
                    n_events_multiplicity_filtered += 1
                    continue
                if (args.max_modules_per_event > 0
                        and n_fired > args.max_modules_per_event):
                    n_events_multiplicity_filtered += 1
                    continue

            n_trigger += 1

            for roc_tag, s, c, cd in C.iter_fadc_channels(fadc_evt):
                name, chan_id, mtype = _key_for(roc_tag, s, c)
                if chan_filter and name not in chan_filter:
                    continue

                samples = np.asarray(cd.samples, dtype=np.uint16)
                ped, rms, peaks = p.wave_ana.analyze(samples)

                # Good pulse selection
                if len(peaks) != 1:
                    continue
                pk = peaks[0]
                if pk.quality != 0:
                    continue
                if pk.overflow:
                    continue
                lo, hi = pk.pos - pre, pk.pos + post + 1
                if lo < 0 or hi > samples.shape[0]:
                    continue
                n_good += 1
                # Additional cuts
                if pk.height < args.height_min:
                    continue
                if pk.height < args.height_rms_mult * rms:
                    continue
                n_after_cut += 1   # after height

                # The C++ fitter takes the raw uint16 slice and does the
                # pedestal subtraction itself.
                slice_u16 = samples[lo:hi]
                rel_peak  = pk.pos - lo

                st = stats.get(name)
                if st is None:
                    st = ChannelStats(name=name, channel_id=chan_id,
                                      module_type=mtype)
                    stats[name] = st
                if args.max_pulses_per_channel and st.n_used >= args.max_pulses_per_channel:
                    continue
                st.n_good += 1  # per-channel clean-pulse count
                st.n_after_cut += 1  # per-channel: after height cuts

                st.n_attempted += 1
                n_pulses_attempted += 1
                st.ped_mean_sum += float(ped)
                st.ped_rms_sum  += float(rms)
                st.ped_n        += 1

                fit = fit_pulse(slice_u16, rel_peak, float(ped), float(rms),
                                clk_ns, args.model_err_floor)
                if not fit.ok:
                    continue
                n_good_fit += 1
                st.n_good_fit += 1
                if fit.t0_ns < args.t0_min or fit.t0_ns > args.t0_max:
                    continue
                n_good_fit_time_cut += 1
                st.n_good_fit_time_cut += 1
                st.tau_r.append(fit.tau_r_ns)
                st.tau_f.append(fit.tau_f_ns)
                st.t0.append(fit.t0_ns)
                if with_p:
                    st.p_list.append(fit.p)
                st.peak_amp.append(fit.peak_amp)
                st.chi2.append(fit.chi2_per_dof)
                st.n_used += 1
                n_pulses_used += 1

                # Cache raw pulses + their peak amps for the diagnostic
                # plots (normalised stack vs median fit).  Pedsub on the
                # fly here so the cache lives in absolute ADC like the
                # fit's input.
                if plotting and len(st.sample_pulses) < PULSE_CACHE:
                    st.sample_pulses.append(
                        slice_u16.astype(np.float64) - ped)
                    st.sample_peak_amps.append(fit.peak_amp)

                # Test A caching: also stratify a copy of this pulse into
                # an amplitude bin so per-channel plots can show separate
                # low/mid/high-amp stacks.
                if plotting:
                    ph = float(pk.height)
                    for label, blo, bhi in AMP_BINS:
                        if blo <= ph < bhi:
                            bin_pulses = st.sample_pulses_by_amp.setdefault(label, [])
                            bin_amps   = st.sample_peak_amps_by_amp.setdefault(label, [])
                            if len(bin_pulses) < PULSE_CACHE_PER_BIN:
                                bin_pulses.append(
                                    slice_u16.astype(np.float64) - ped)
                                bin_amps.append(fit.peak_amp)
                            break
    except KeyboardInterrupt:
        print("\n[interrupted — writing partial results]", flush=True)

    elapsed = loop_stats.elapsed
    print(f"[done] {loop_stats.n_phys} phys events  /  {n_pulses_attempted} fits attempted"
          f"  /  {n_pulses_used} converged  /  {len(stats)} channels"
          f"  /  trigger_filtered={n_events_trigger_filtered}"
          f"  /  mult_filtered={n_events_multiplicity_filtered}"
          f"  /  {elapsed:.1f}s", flush=True)

    out: Dict = {
        "_meta": {
            "inputs": list(args.evio_paths),
            "n_evio_splits": len(p.evio_files),
            "n_phys_events": loop_stats.n_phys,
            "n_total_events": loop_stats.n_phys,
            "n_events_after_trigger": n_trigger,
            "n_pulses_attempted": n_pulses_attempted,
            "n_pulses_used": n_pulses_used,
            "n_pulses_good": n_good,
            "n_pulses_after_cut": n_after_cut,
            "n_pulses_good_fit":          n_good_fit,
            "n_pulses_good_fit_time_cut": n_good_fit_time_cut,
            "n_channels": len(stats),
            "model": args.model,
            "model_formula": (
                "((1 - exp(-(t-t0)/tau_r)) * exp(-(t-t0)/tau_f)) / T_max(tau_r, tau_f),  t > t0"
                if args.model == "two_tau" else
                "([1 - exp(-(t-t0)/tau_r)]^p * exp(-(t-t0)/tau_f)) / T_max(tau_r, tau_f, p),  t > t0   "
                "[two-tau with rise-edge exponent]"
            ),
            "fit_params": (["t0_ns", "tau_r_ns", "tau_f_ns"]
                           if args.model == "two_tau"
                           else ["t0_ns", "tau_r_ns", "tau_f_ns", "p"]),
            "amplitude_handling": "each pulse divided by its own pedsub peak height before the fit; peak_amp_adc records the un-normalised peak amplitude per channel for downstream diagnostics",
            "sigma_per_sample": "ped.rms / pulse_peak_amp (i.e. relative noise on the normalised pulse) — keeps chi2/dof comparable across channels with very different signal amplitudes",
            "clk_ns": clk_ns,
            "pre_samples": pre,
            "post_samples": post,
            "min_pulses": args.min_pulses,
            "chi2_max": args.chi2_max,
            "height_min": args.height_min,
            "height_rms_mult": args.height_rms_mult,
            "t0_min": None if args.t0_min == float("-inf") else args.t0_min,
            "t0_max": None if args.t0_max == float("inf") else args.t0_max,
            "model_err_floor": args.model_err_floor,
            "trigger_accept_spec": args.trigger_accept or None,
            "trigger_reject_spec": args.trigger_reject or None,
            "trigger_accept_mask": accept_mask or None,
            "trigger_reject_mask": reject_mask or None,
            "trigger_event_type_spec": args.trigger_event_type or None,
            "trigger_event_type_accept": sorted(event_type_accept) if event_type_accept else None,
            "n_events_trigger_filtered": n_events_trigger_filtered,
            "min_modules_per_event": args.min_modules_per_event,
            "max_modules_per_event": args.max_modules_per_event,
            "module_count_height_min": (
                args.module_count_height_min
                if args.module_count_height_min >= 0
                else args.height_min
            ),
            "module_count_require_clean": args.module_count_require_clean,
            "n_events_multiplicity_filtered": n_events_multiplicity_filtered,
        }
    }
    summaries: List[Dict] = []
    n_skipped_empty = 0
    for name in sorted(stats):
        if stats[name].n_used == 0:
            n_skipped_empty += 1
            continue
        rec = finalize_channel(stats[name], args.min_pulses, args.chi2_max)
        out[name] = rec
        summaries.append({"name": name, **rec})

    if n_skipped_empty:
        print(f"[write] skipped {n_skipped_empty} channels with n_pulses_used=0 "
              f"(all fits filtered by t0 cut or LM failure)", flush=True)

    # The per-channel entries above are kept for diagnostics and so the
    # store can build its (roc, slot, channel) → module_type lookup at
    # load time.
    out["_by_type"] = aggregate_by_type(summaries, args.min_pulses)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_sanitized = _replace_nan_recursive(out)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out_sanitized, f, indent=2, sort_keys=False, allow_nan=False)
    print(f"[write] {out_path}", flush=True)

    # Per-pulse (peak_amp, chi2, module_type, t0, tau_r, tau_f, name) dump for
    # χ²-vs-amplitude and τ-vs-amplitude diagnostics.
    # (sigma_per_sample = ped_rms / peak_amp, so low-amplitude pulses carry
    # a larger relative noise weighting; this dump lets plot_chi2_vs_amp.py
    # test whether the chi2/dof threshold preferentially rejects them, and
    # plot_tau_vs_amp.py test whether pulse shape parameters are amplitude-
    # dependent — a signature of PMT saturation.)
    if plotting:
        all_amps  = []
        all_chi2  = []
        all_types = []
        all_t0    = []
        all_tau_r = []
        all_tau_f = []
        all_p     = []
        all_names = []
        has_p     = False
        for st in stats.values():
            # st.peak_amp and st.chi2 are always the same length (both
            # appended together only on fit convergence).
            n = len(st.peak_amp)
            if n == 0:
                continue
            all_amps.extend(st.peak_amp)
            all_chi2.extend(st.chi2)
            all_types.extend([st.module_type] * n)
            all_t0.extend(st.t0)
            all_tau_r.extend(st.tau_r)
            all_tau_f.extend(st.tau_f)
            all_names.extend([st.name] * n)
            # p_list is only populated for --model two_tau_p, else empty
            if st.p_list:
                has_p = True
                all_p.extend(st.p_list)
        if all_amps:
            npz_path = plot_dir / "per_pulse_amp_chi2.npz"
            npz_path.parent.mkdir(parents=True, exist_ok=True)
            save_kwargs = dict(
                amp=np.asarray(all_amps,  dtype=np.float32),
                chi2=np.asarray(all_chi2, dtype=np.float32),
                mtype=np.asarray(all_types),   # numpy will pick a string dtype
                t0=np.asarray(all_t0,     dtype=np.float32),
                tau_r=np.asarray(all_tau_r, dtype=np.float32),
                tau_f=np.asarray(all_tau_f, dtype=np.float32),
                name=np.asarray(all_names),
            )
            if has_p:
                save_kwargs["p"] = np.asarray(all_p, dtype=np.float32)
            np.savez_compressed(npz_path, **save_kwargs)
            print(f"[write] {npz_path}  ({len(all_amps)} pulses)", flush=True)

    # Test A dump: per-channel raw waveforms stratified by amplitude bin.
    # One .npz per channel that has at least one bin populated; downstream
    # plot_raw_pulses_by_amp.py consumes these.
    if plotting:
        wf_dir = plot_dir / "waveforms_by_amp"
        n_written = 0
        for name, st in stats.items():
            if not st.sample_pulses_by_amp:
                continue
            # Only write channels with at least 5 pulses in at least two bins,
            # otherwise the overlay comparison is uninformative.
            populated = [lbl for lbl, pulses in st.sample_pulses_by_amp.items()
                         if len(pulses) >= 5]
            if len(populated) < 2:
                continue
            wf_dir.mkdir(parents=True, exist_ok=True)
            # One npz per channel; keys are the bin labels; each stores a
            # 2D array of shape (n_pulses, n_samples) plus an "amps" 1D array.
            save_kwargs = {"module_type": np.asarray(st.module_type)}
            for lbl in populated:
                pulses = st.sample_pulses_by_amp[lbl]
                amps   = st.sample_peak_amps_by_amp[lbl]
                # Pulses on this channel all have the same slice length by
                # construction (pre+post samples, checked in the event loop),
                # so stacking is safe.
                save_kwargs[f"pulses_{lbl}"] = np.stack(pulses).astype(np.float32)
                save_kwargs[f"amps_{lbl}"]   = np.asarray(amps, dtype=np.float32)
            np.savez_compressed(wf_dir / f"{name}.npz", **save_kwargs)
            n_written += 1
        if n_written:
            print(f"[write] {wf_dir}/  ({n_written} channels)", flush=True)

    # ---- plotting ----------------------------------------------------------
    if not plotting:
        return 0
    plt = C.import_pyplot()
    if plt is None:
        print("[plot] matplotlib unavailable — skipping", flush=True)
        return 0

    good_names, bad_names = select_plot_targets(
        stats, args.min_pulses, args.chi2_max,
        args.n_plot_good, args.n_plot_bad)

    print(f"[plot] dir={plot_dir}  "
          f"good={len(good_names)}  bad={len(bad_names)}  "
          f"explicit={len(extra_plot)}", flush=True)

    targets = ([(n, "good", "[good]") for n in good_names]
               + [(n, "bad", "[bad]") for n in bad_names]
               + [(n, "explicit", "") for n in sorted(extra_plot)])
    for name, sub, label in targets:
        st = stats.get(name)
        if st is None or not st.sample_pulses:
            print(f"[plot] no pulses for {sub} channel {name}", flush=True)
            continue
        plot_channel(plt, st, clk_ns, pre, post,
                     plot_dir / sub / f"{name}.png", label=label)

    if not args.no_summary_plot:
        plot_summary(plt, summaries, args.min_pulses, args.chi2_max,
                     plot_dir / "summary.png")
        print(f"[plot] wrote {plot_dir / 'summary.png'}", flush=True)
        plot_summary_per_pulse(plt, stats, args.min_pulses,
                               plot_dir / "summary_per_pulse.png")
        print(f"[plot] wrote {plot_dir / 'summary_per_pulse.png'}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
