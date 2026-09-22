#!/usr/bin/env python3
"""
benchmark_matched_filter.py — benchmark a matched-filter (MF) peak finder on
synthetic pile-up events from ``pileup_generator.py``.

Overview
--------
The matched filter cross-correlates a pedestal-subtracted waveform with an
analytic two-τ template kernel.  Peaks in the MF output correspond to pulse
arrivals; the MF peak height is proportional to the pulse amplitude.  This
algorithm is compared against the existing ``WaveAnalyzer.findPeaks`` routine
on the same synthetic 2-pulse waveforms to quantify improvement in pile-up
detection at small ΔT and low amplitude ratios — the regime where the current
peak finder fails.

Matched-filter coordinate convention
--------------------------------------
The kernel has ``kernel_samples`` entries, evaluated at sample offsets 0, 1,
…, kernel_samples-1 (i.e. the template onset is pinned to sample 0 of the
kernel).  ``np.correlate(waveform, kernel, mode='valid')`` therefore produces
an array of length ``n_samples - kernel_samples + 1``.  MF output index ``j``
corresponds to waveform sample ``j`` (the correlation window covers waveform
samples j … j+kernel_samples-1).  In other words:

    mf_output[j]  ≈  amplitude of a pulse whose template-onset is at
                     waveform sample j.

Because the truth is stored as *peak* times (not onset times), and the template
peak sits at sample offset ``t_peak_offset_ns / clk_ns`` from the onset, the
truth peak position in MF-output coordinates is:

    truth_mf_idx = truth_peak_times_ns / clk_ns   (already in peak coords)

The truth-matching window (default ±2 samples = ±8 ns) is applied directly
in MF-output / sample-index space.  Because mf_output[j] captures a pulse
onset at sample j, and the template peak is ~1–2 samples later, the ±2
window is deliberately generous to cover this sub-sample offset.

Usage
-----
    python benchmark_matched_filter.py <synth_dir> \\
        [--template pulse_templates.json] \\
        [--out-dir results/] \\
        [--daq-config database/daq_config.json] \\
        [--kernel-samples 30] \\
        [--mf-nsigma 5.0] \\
        [--min-separation 3] \\
        [--plot] \\
        [--plot-examples N]

Output files (written to <out-dir>, defaulting to <synth_dir>)
--------------------------------------------------------------
    mf_benchmark.json          — per-cell stats + global summary
    mf_vs_wa_detection.png     — side-by-side WA vs MF efficiency heatmap
                                 (only when --plot)
    mf_example_dt*.png         — waveform + MF overlay examples
                                 (only when --plot-examples > 0)
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# prad2py import — deferred so argparse --help works even without bindings.
# ---------------------------------------------------------------------------

try:
    from prad2py import dec as _dec
    _PRAD2PY_AVAILABLE = True
except ImportError as _exc:
    _dec = None  # type: ignore[assignment]
    _PRAD2PY_AVAILABLE = False
    _PRAD2PY_IMPORT_ERR = str(_exc)


def _require_prad2py() -> None:
    """Raise SystemExit if prad2py is not available."""
    if not _PRAD2PY_AVAILABLE:
        raise SystemExit(
            f"[ERROR] cannot import prad2py: {_PRAD2PY_IMPORT_ERR}\n"
            "        Build the python bindings (cmake -DBUILD_PYTHON=ON) and "
            "ensure the install directory is on PYTHONPATH."
        )


# ---------------------------------------------------------------------------
# Path helper (mirrors benchmark_deconv.py)
# ---------------------------------------------------------------------------

def _resolve_db_path(p: str) -> str:
    """Resolve a possibly-relative path against PRAD2_DATABASE_DIR."""
    if not p or os.path.isabs(p):
        return p
    db = os.environ.get("PRAD2_DATABASE_DIR")
    return os.path.join(db, p) if db else p


# ---------------------------------------------------------------------------
# Manifest loading
# ---------------------------------------------------------------------------

def _load_manifest(synth_dir: Path) -> Dict:
    """Load and return manifest.json from synth_dir.  Raises SystemExit on error."""
    mpath = synth_dir / "manifest.json"
    if not mpath.is_file():
        raise SystemExit(f"[ERROR] manifest.json not found in {synth_dir}")
    try:
        return json.loads(mpath.read_text())
    except json.JSONDecodeError as exc:
        raise SystemExit(f"[ERROR] manifest.json parse error: {exc}")


# ---------------------------------------------------------------------------
# Template parameter loading
# ---------------------------------------------------------------------------

def _load_template_params(
    manifest: Dict,
    template_override: str,
) -> Tuple[float, float]:
    """Return (tau_r_ns, tau_f_ns) from the template JSON's _by_type entry.

    Prefers ``template_override`` CLI argument; falls back to the manifest's
    ``template_source`` field.  Raises SystemExit on any error.
    """
    tmpl_rel = template_override or manifest.get("template_source", "")
    tmpl_path = _resolve_db_path(tmpl_rel)
    if not tmpl_path or not Path(tmpl_path).is_file():
        raise SystemExit(
            f"[ERROR] pulse template file not found: {tmpl_path!r}\n"
            "        Use --template to specify the path explicitly."
        )

    try:
        data = json.loads(Path(tmpl_path).read_text())
    except json.JSONDecodeError as exc:
        raise SystemExit(f"[ERROR] template JSON parse error: {exc}")

    material = manifest.get("material", "PbWO4")
    by_type = data.get("_by_type", {})
    if material not in by_type:
        avail = ", ".join(sorted(by_type.keys())) if by_type else "(none)"
        raise SystemExit(
            f"[ERROR] material {material!r} not found in template '_by_type'; "
            f"available: {avail}"
        )

    entry = by_type[material]
    try:
        tau_r = float(entry["tau_r_ns"]["median"])
        tau_f = float(entry["tau_f_ns"]["median"])
    except (KeyError, TypeError, ValueError) as exc:
        raise SystemExit(
            f"[ERROR] cannot read tau_r_ns/tau_f_ns medians for "
            f"{material!r}: {exc}"
        )

    if math.isnan(tau_r) or math.isnan(tau_f):
        raise SystemExit(
            f"[ERROR] NaN in template medians for {material!r}: "
            f"tau_r={tau_r}, tau_f={tau_f}"
        )
    if tau_r <= 0 or tau_f <= 0:
        raise SystemExit(
            f"[ERROR] non-positive tau for {material!r}: "
            f"tau_r={tau_r}, tau_f={tau_f}"
        )

    print(f"[setup] template   : {tmpl_path}", flush=True)
    print(f"[setup] material   : {material}  "
          f"tau_r={tau_r:.3f} ns  tau_f={tau_f:.3f} ns", flush=True)
    return tau_r, tau_f


# ---------------------------------------------------------------------------
# WaveAnalyzer setup (mirrors benchmark_deconv.py)
# ---------------------------------------------------------------------------

def _setup_wave_analyzer(daq_config: str):
    """Load DAQ config and return (WaveAnalyzer, clk_ns)."""
    daq_cfg_path = daq_config or _resolve_db_path("daq_config.json")
    cfg = _dec.load_daq_config(daq_cfg_path) if daq_cfg_path else _dec.load_daq_config()
    print(f"[setup] DAQ config : {daq_cfg_path or '(default)'}", flush=True)
    wcfg = _dec.WaveConfig(cfg.wave_cfg)
    wa = _dec.WaveAnalyzer(wcfg)
    clk_mhz = float(wcfg.clk_mhz) if wcfg.clk_mhz > 0 else 250.0
    clk_ns = 1000.0 / clk_mhz
    print(f"[setup] clk_ns     : {clk_ns:.4f}  ({clk_mhz} MHz)", flush=True)
    return wa, clk_ns


# ---------------------------------------------------------------------------
# Matched-filter kernel construction
# ---------------------------------------------------------------------------

def build_kernel(
    tau_r: float,
    tau_f: float,
    clk_ns: float,
    kernel_samples: int = 30,
) -> Tuple[np.ndarray, float]:
    """Build the matched-filter kernel from the analytic two-τ template.

    The kernel is the peak-normalised two-τ pulse shape sampled at
    ``kernel_samples`` points starting from t=0 (template onset at t=0).
    Dividing the correlation by ``kernel_norm`` = ||kernel||² converts the
    raw correlation into an amplitude estimate at each waveform position.

    Parameters
    ----------
    tau_r          : rise time constant (ns)
    tau_f          : fall time constant (ns)
    clk_ns         : nanoseconds per ADC sample
    kernel_samples : number of kernel samples (default 30 = 120 ns at 4 ns/sample)

    Returns
    -------
    kernel      : float64 array of shape (kernel_samples,), peak-normalised
    kernel_norm : float, sum of squared kernel values (for amplitude scaling)
    """
    t = np.arange(kernel_samples) * clk_ns
    u = tau_r / (tau_r + tau_f)
    T_max = (1.0 - u) * (u ** (tau_r / tau_f))
    kernel = np.zeros(kernel_samples, dtype=np.float64)
    mask = t > 0
    kernel[mask] = (
        (1.0 - np.exp(-t[mask] / tau_r)) * np.exp(-t[mask] / tau_f) / T_max
    )
    kernel_norm = float(np.sum(kernel ** 2))
    return kernel, kernel_norm


# ---------------------------------------------------------------------------
# Matched-filter application
# ---------------------------------------------------------------------------

def matched_filter(
    waveform_pedsub: np.ndarray,
    kernel: np.ndarray,
    kernel_norm: float,
) -> np.ndarray:
    """Cross-correlate pedestal-subtracted waveform with the template kernel.

    Uses ``np.correlate(..., mode='valid')`` so the output length is
    ``n_samples - kernel_samples + 1``.  MF output index ``j`` corresponds to
    waveform sample ``j`` (pulse onset at sample j).  See the module docstring
    for the full coordinate convention.

    Parameters
    ----------
    waveform_pedsub : 1-D float array, length n_samples
    kernel          : 1-D float array, length kernel_samples
    kernel_norm     : sum(kernel**2) — used to convert to amplitude units

    Returns
    -------
    mf_output : 1-D float array, length n_samples - kernel_samples + 1
    """
    raw = np.correlate(waveform_pedsub.astype(np.float64), kernel, mode="valid")
    return raw / kernel_norm


# ---------------------------------------------------------------------------
# MF peak finding
# ---------------------------------------------------------------------------

def find_mf_peaks(
    mf_output: np.ndarray,
    threshold: float,
    min_separation: int,
    max_peaks: int = 6,
    min_prominence_frac: float = 0.3,
) -> List[Tuple[int, float]]:
    """Find peaks in matched-filter output above threshold.

    Algorithm
    ---------
    1. Scan for all local maxima strictly above ``threshold``.
    2. Reject peaks whose prominence is below ``min_prominence_frac`` of their
       own amplitude.  Prominence is peak_value minus the higher of the two
       local minima found by walking left and right until a higher sample or
       a distance of ``min_separation * 3`` is reached.  This suppresses
       shallow tail-crossing artefacts that appear just above threshold on the
       falling edge of the MF response.
    3. Sort by amplitude descending.
    4. Greedily keep peaks whose position is at least ``min_separation``
       samples from any already-kept peak.
    5. Return up to ``max_peaks`` peaks.

    Parameters
    ----------
    mf_output           : 1-D float array (MF output from ``matched_filter()``)
    threshold           : minimum amplitude to consider (mf_nsigma × ped_rms)
    min_separation      : minimum inter-peak separation in MF-output samples
    max_peaks           : maximum number of peaks to return
    min_prominence_frac : minimum prominence as a fraction of the detection
                          threshold (default 0.3).  Set to 0 to disable.

    Returns
    -------
    List of (sample_index, mf_amplitude) tuples sorted by amplitude descending.
    ``sample_index`` is in MF-output coordinates, which equals waveform sample
    coordinates (see module docstring for the coordinate mapping).
    """
    peaks: List[Tuple[int, float]] = []
    n = len(mf_output)
    for i in range(1, n - 1):
        if (mf_output[i] > mf_output[i - 1]
                and mf_output[i] >= mf_output[i + 1]
                and mf_output[i] > threshold):
            peaks.append((i, float(mf_output[i])))

    # Prominence filter: reject peaks whose rise above the surrounding baseline
    # is less than min_prominence_frac × peak_amplitude.  Real pulse peaks have
    # high prominence; tail-crossing artefacts have near-zero prominence.
    if min_prominence_frac > 0:
        filtered: List[Tuple[int, float]] = []
        for pos, amp in peaks:
            # Walk left to find the minimum before hitting a higher sample.
            left_min = amp
            for j in range(pos - 1, max(0, pos - min_separation * 3) - 1, -1):
                if mf_output[j] < left_min:
                    left_min = float(mf_output[j])
                if mf_output[j] > amp:
                    break  # hit a higher peak — stop

            # Walk right to find the minimum before hitting a higher sample.
            right_min = amp
            for j in range(pos + 1, min(n, pos + min_separation * 3 + 1)):
                if mf_output[j] < right_min:
                    right_min = float(mf_output[j])
                if mf_output[j] > amp:
                    break  # hit a higher peak — stop

            base = max(left_min, right_min)
            prominence = amp - base

            if prominence >= min_prominence_frac * threshold:
                filtered.append((pos, amp))
        peaks = filtered

    # Sort by amplitude descending.
    peaks.sort(key=lambda p: -p[1])

    # Greedy non-maximum suppression by proximity.
    kept: List[Tuple[int, float]] = []
    for pos, amp in peaks:
        too_close = any(abs(pos - kpos) < min_separation for kpos, _ in kept)
        if not too_close:
            kept.append((pos, amp))
        if len(kept) >= max_peaks:
            break

    return kept


# ---------------------------------------------------------------------------
# Per-cell accumulator
# ---------------------------------------------------------------------------

class MFCellStats:
    """Accumulate per-event results for one (ΔT, ratio) grid cell."""

    def __init__(self, dt_ns: float, ratio: float) -> None:
        self.dt_ns = dt_ns
        self.ratio = ratio
        self.n_events: int = 0
        self.n_wa_detected: int = 0   # WaveAnalyzer found ≥ 2 peaks
        self.n_mf_detected: int = 0   # MF found peak near truth injected position

    def record(self, wa_detected: bool, mf_detected: bool) -> None:
        self.n_events += 1
        if wa_detected:
            self.n_wa_detected += 1
        if mf_detected:
            self.n_mf_detected += 1

    @property
    def wa_eff(self) -> float:
        return (self.n_wa_detected / self.n_events * 100.0
                if self.n_events > 0 else float("nan"))

    @property
    def mf_eff(self) -> float:
        return (self.n_mf_detected / self.n_events * 100.0
                if self.n_events > 0 else float("nan"))

    @property
    def improvement(self) -> float:
        wa = self.wa_eff
        mf = self.mf_eff
        if math.isnan(wa) or math.isnan(mf):
            return float("nan")
        return mf - wa

    def to_dict(self) -> Dict:
        return {
            "dt_ns":        self.dt_ns,
            "ratio":        self.ratio,
            "n_events":     self.n_events,
            "n_wa_detected": self.n_wa_detected,
            "n_mf_detected": self.n_mf_detected,
            "wa_eff_pct":   self.wa_eff,
            "mf_eff_pct":   self.mf_eff,
            "improvement_pct": self.improvement,
        }


# ---------------------------------------------------------------------------
# ΔT-bin helpers for random-dt mode (mirrors benchmark_deconv.py)
# ---------------------------------------------------------------------------

def _make_dt_bins(dt_min: float, dt_max: float, n_bins: int) -> np.ndarray:
    return np.linspace(dt_min, dt_max, n_bins + 1)


def _dt_bin_index(dt: float, edges: np.ndarray) -> int:
    idx = int(np.searchsorted(edges, dt, side="right")) - 1
    return max(0, min(len(edges) - 2, idx))


def _bin_center(edges: np.ndarray, k: int) -> float:
    return float(0.5 * (edges[k] + edges[k + 1]))


# ---------------------------------------------------------------------------
# Example-plot collection (populated during processing)
# ---------------------------------------------------------------------------

# Stores (pedsub, mf_output, mf_threshold, mf_peaks, truth_inj_sample,
#         truth_base_sample, dt_ns, ratio, ev_index) for each example.
_example_buffer: List[Dict] = []
_example_cells_seen: set = set()  # (dt_key, ratio_key) cells already sampled


# ---------------------------------------------------------------------------
# Core per-.npz processing loop
# ---------------------------------------------------------------------------

def _process_npz(
    npz_path: Path,
    wa: object,
    kernel: np.ndarray,
    kernel_norm: float,
    clk_ns: float,
    mf_nsigma: float,
    min_separation: int,
    match_window: int,
    manifest: Dict,
    cell_map: Dict[Tuple[float, float], MFCellStats],
    random_dt: bool,
    dt_edges: Optional[np.ndarray],
    collect_examples: int,
    min_prominence_frac: float = 0.3,
    n_warn_limit: int = 5,
) -> int:
    """Process one .npz file.  Returns the number of events processed."""

    try:
        data = np.load(str(npz_path), allow_pickle=False)
    except Exception as exc:
        print(f"[WARN] cannot load {npz_path.name}: {exc}", file=sys.stderr)
        return 0

    waveforms           = data["waveforms"]           # uint16 (N, n_samples)
    truth_peak_times_ns = data["truth_peak_times_ns"] # float32 (N, 2)
    event_dt_ns_arr     = data["event_dt_ns"]         # float32 (N,)
    pedestals           = data["pedestals"]           # float32 (N, 2): [ped_mean, ped_rms]

    meta_ratio = round(float(data["meta_ratio"]), 4)

    N = waveforms.shape[0]
    n_warn = 0

    for i in range(N):
        samples = waveforms[i]       # uint16 (n_samples,)
        dt_ev = float(event_dt_ns_arr[i])

        # Determine cell bucket.
        if random_dt:
            if dt_edges is None:
                continue
            bin_k = _dt_bin_index(dt_ev, dt_edges)
            dt_key = round(_bin_center(dt_edges, bin_k), 2)
        else:
            dt_key = round(float(dt_ev), 2)

        cell_key = (dt_key, round(meta_ratio, 4))
        cell = cell_map.get(cell_key)
        if cell is None:
            if n_warn < n_warn_limit:
                print(f"[WARN] no cell for key {cell_key}; skipping event {i}",
                      file=sys.stderr)
                n_warn += 1
            continue

        # ---- Pedestal from stored truth (pileup_generator stores ped_mean,
        #      ped_rms in pedestals[:, 0] and pedestals[:, 1]).
        ped_mean = float(pedestals[i, 0])
        ped_rms  = float(pedestals[i, 1])

        # Also run WaveAnalyzer.analyze_result to get its peak count.
        try:
            wres = wa.analyze_result(samples)
            wa_n_peaks = len(list(wres.peaks))
        except Exception as exc:
            if n_warn < n_warn_limit:
                print(f"[WARN] analyze_result failed for event {i}: {exc}",
                      file=sys.stderr)
                n_warn += 1
            wa_n_peaks = 0

        wa_detected = wa_n_peaks >= 2

        # ---- Pedestal-subtracted waveform ----
        pedsub = samples.astype(np.float64) - ped_mean

        # ---- Matched filter ----
        mf_output = matched_filter(pedsub, kernel, kernel_norm)

        # Threshold = mf_nsigma × ped_rms, expressed in amplitude units.
        # Since kernel_norm = sum(kernel^2) and a noise sample has std = ped_rms,
        # the noise std of the MF output is ped_rms * sqrt(sum(kernel^2)) / kernel_norm
        # = ped_rms / sqrt(kernel_norm).  We threshold on this effective noise.
        if ped_rms > 0.0:
            mf_threshold = mf_nsigma * ped_rms / math.sqrt(kernel_norm)
        else:
            mf_threshold = mf_nsigma  # fallback if ped_rms unknown

        mf_peaks = find_mf_peaks(mf_output, mf_threshold, min_separation,
                                  max_peaks=6, min_prominence_frac=min_prominence_frac)

        # ---- Truth matching ----
        # Truth injected pulse peak position in sample (= MF output) coordinates.
        # truth_peak_times_ns[:, 0] is the *peak* time of the injected pulse.
        # MF output index j corresponds to waveform sample j (onset pinned there).
        # The MF peak for a pulse with onset at sample j_onset will appear near
        # j_onset (kernel captures the full rising+falling shape from onset).
        # We match against the peak time directly: truth_peak_sample ≈ peak_time/clk_ns.
        truth_inj_peak_sample = float(truth_peak_times_ns[i, 0]) / clk_ns

        mf_detected = any(
            abs(pos - truth_inj_peak_sample) <= match_window
            for pos, _ in mf_peaks
        )

        cell.record(wa_detected=wa_detected, mf_detected=mf_detected)

        # ---- Collect example plots (cells where MF detects but WA doesn't) ----
        if (collect_examples > 0
                and len(_example_buffer) < collect_examples
                and mf_detected and not wa_detected
                and cell_key not in _example_cells_seen):
            _example_cells_seen.add(cell_key)
            truth_base_sample = float(truth_peak_times_ns[i, 1]) / clk_ns
            _example_buffer.append({
                "pedsub":             pedsub.copy(),
                "mf_output":          mf_output.copy(),
                "mf_threshold":       mf_threshold,
                "mf_peaks":           list(mf_peaks),
                "truth_inj_sample":   truth_inj_peak_sample,
                "truth_base_sample":  truth_base_sample,
                "dt_ns":              dt_ev,
                "ratio":              meta_ratio,
                "ev_index":           i,
                "clk_ns":             clk_ns,
                "kernel_samples":     len(kernel),
            })

    return N


# ---------------------------------------------------------------------------
# Summary table printer
# ---------------------------------------------------------------------------

def _print_table(cells: List[MFCellStats]) -> None:
    """Print fixed-width ASCII comparison table to stdout."""
    print()
    header = (f"{'ΔT':>6}  {'ratio':>6}  {'N':>5}  "
              f"{'WA_det':>6}  {'MF_det':>6}  "
              f"{'WA_eff%':>7}  {'MF_eff%':>7}  {'improvement':>11}")
    sep    = (f"{'-----':>6}  {'------':>6}  {'---':>5}  "
              f"{'------':>6}  {'------':>6}  "
              f"{'-------':>7}  {'-------':>7}  {'-----------':>11}")
    print(header)
    print(sep)
    for c in cells:
        dt_str    = f"{c.dt_ns:6.1f}" if math.isfinite(c.dt_ns) else "   rnd"
        ratio_str = f"{c.ratio:6.3f}"
        n_str     = f"{c.n_events:5d}"
        wa_d_str  = f"{c.n_wa_detected:6d}"
        mf_d_str  = f"{c.n_mf_detected:6d}"
        wa_e = c.wa_eff
        mf_e = c.mf_eff
        imp  = c.improvement
        wa_eff_str  = f"{wa_e:6.1f}%" if math.isfinite(wa_e) else "    NaN"
        mf_eff_str  = f"{mf_e:6.1f}%" if math.isfinite(mf_e) else "    NaN"
        imp_str     = (f"{imp:+.1f}%" if math.isfinite(imp) else "    NaN")
        print(f"{dt_str}  {ratio_str}  {n_str}  "
              f"{wa_d_str}  {mf_d_str}  "
              f"{wa_eff_str:>7}  {mf_eff_str:>7}  {imp_str:>11}")
    print()


# ---------------------------------------------------------------------------
# Heatmap plotting
# ---------------------------------------------------------------------------

def _make_heatmaps(cells: List[MFCellStats], out_dir: Path) -> None:
    """Produce side-by-side WA vs MF detection efficiency heatmaps."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dt_vals    = sorted({c.dt_ns for c in cells if math.isfinite(c.dt_ns)})
    ratio_vals = sorted({c.ratio for c in cells})

    if not dt_vals or not ratio_vals:
        print("[plot] insufficient data for heatmaps (no finite ΔT values).",
              file=sys.stderr)
        return

    n_dt    = len(dt_vals)
    n_ratio = len(ratio_vals)
    dt_idx  = {v: k for k, v in enumerate(dt_vals)}
    r_idx   = {v: k for k, v in enumerate(ratio_vals)}

    wa_eff = np.full((n_ratio, n_dt), float("nan"))
    mf_eff = np.full((n_ratio, n_dt), float("nan"))

    for c in cells:
        if not math.isfinite(c.dt_ns):
            continue
        ki = dt_idx.get(c.dt_ns)
        ri = r_idx.get(c.ratio)
        if ki is None or ri is None:
            continue
        wa_eff[ri, ki] = c.wa_eff / 100.0  # store as fraction [0,1]
        mf_eff[ri, ki] = c.mf_eff / 100.0

    # Shared colour limits for direct visual comparison.
    vmin, vmax = 0.0, 1.0

    fig, axes = plt.subplots(
        1, 2,
        figsize=(max(10, n_dt * 1.6 + 2), max(4, n_ratio * 0.7 + 1.5)),
        sharey=True,
    )

    xtick_labels = [f"{v:.0f}" for v in dt_vals]
    ytick_labels = [f"{v:.2g}" for v in ratio_vals]

    common_kw = dict(
        aspect="auto", origin="lower",
        vmin=vmin, vmax=vmax, cmap="viridis",
        extent=[-0.5, n_dt - 0.5, -0.5, n_ratio - 0.5],
    )

    for ax, data, title_label in (
        (axes[0], wa_eff, "WaveAnalyzer (≥2 peaks)"),
        (axes[1], mf_eff, "Matched Filter"),
    ):
        img = ax.imshow(data, **common_kw)
        ax.set_xticks(range(n_dt))
        ax.set_xticklabels(xtick_labels, rotation=45, ha="right")
        ax.set_yticks(range(n_ratio))
        ax.set_yticklabels(ytick_labels)
        ax.set_xlabel("ΔT (ns)")
        ax.set_title(f"Detection efficiency\n{title_label}")

    axes[0].set_ylabel("Amplitude ratio (inj / base)")

    # Single shared colorbar on the right.
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    sm = plt.cm.ScalarMappable(cmap="viridis",
                               norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label("Detection efficiency")

    p = out_dir / "mf_vs_wa_detection.png"
    fig.savefig(str(p), dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] wrote {p}", flush=True)


# ---------------------------------------------------------------------------
# Example waveform + MF overlay plots
# ---------------------------------------------------------------------------

def _plot_examples(out_dir: Path) -> None:
    """Produce overlay plots for the collected example events."""
    if not _example_buffer:
        return

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for ex in _example_buffer:
        pedsub     = ex["pedsub"]
        mf_output  = ex["mf_output"]
        threshold  = ex["mf_threshold"]
        mf_peaks   = ex["mf_peaks"]
        t_inj      = ex["truth_inj_sample"]
        t_base     = ex["truth_base_sample"]
        dt_ns      = ex["dt_ns"]
        ratio      = ex["ratio"]
        ev_idx     = ex["ev_index"]
        clk_ns     = ex["clk_ns"]
        k_len      = ex["kernel_samples"]

        n_samp     = len(pedsub)
        mf_len     = len(mf_output)
        t_wave_ns  = np.arange(n_samp) * clk_ns
        # MF output index j → waveform sample j (see module docstring).
        t_mf_ns    = np.arange(mf_len) * clk_ns

        fig, (ax_top, ax_bot) = plt.subplots(
            2, 1, figsize=(9, 5), sharex=False
        )

        # ---- Top panel: pedsub waveform ----
        ax_top.plot(t_wave_ns, pedsub, color="0.3", lw=1.2,
                    label="pedsub waveform")
        ax_top.axvline(t_inj * clk_ns, color="C1", ls="--", lw=1.2,
                       label=f"truth inj peak ({t_inj * clk_ns:.1f} ns)")
        ax_top.axvline(t_base * clk_ns, color="C2", ls="--", lw=1.2,
                       label=f"truth base peak ({t_base * clk_ns:.1f} ns)")
        # Truth markers on the waveform
        inj_sample  = int(round(t_inj))
        base_sample = int(round(t_base))
        if 0 <= inj_sample < len(pedsub):
            ax_top.plot([t_inj * clk_ns], [pedsub[inj_sample]], "v",
                        color="C1", ms=12, mec="k", mew=0.8,
                        label="injected peak (truth)")
        if 0 <= base_sample < len(pedsub):
            ax_top.plot([t_base * clk_ns], [pedsub[base_sample]], "v",
                        color="C2", ms=12, mec="k", mew=0.8,
                        label="base peak (truth)")
        ax_top.set_ylabel("ADC − ped")
        ax_top.set_title(
            f"MF example: ΔT={dt_ns:.1f} ns  ratio={ratio:.3f}  ev={ev_idx}"
        )
        ax_top.legend(fontsize=8, loc="upper right")
        ax_top.grid(True, alpha=0.25)

        # ---- Bottom panel: MF output ----
        ax_bot.plot(t_mf_ns, mf_output, color="C0", lw=1.2, label="MF output")
        ax_bot.axhline(threshold, color="C3", ls=":", lw=1.0,
                       label=f"MF threshold ({threshold:.1f} ADC)")
        for idx_pk, (pos, amp) in enumerate(mf_peaks):
            lbl = f"detected MF peaks (n={len(mf_peaks)})" if idx_pk == 0 else ""
            ax_bot.plot([pos * clk_ns], [amp], "D", color="C3",
                        ms=10, mfc="C3", mec="k", mew=1.0, label=lbl)
        ax_bot.axvline(t_inj * clk_ns, color="C1", ls="--", lw=2.0, alpha=0.5)
        ax_bot.axvline(t_base * clk_ns, color="C2", ls="--", lw=2.0, alpha=0.5)
        # Truth markers on the MF curve
        for t_truth, color, lbl in [(t_inj, "C1", "truth inj"),
                                     (t_base, "C2", "truth base")]:
            mf_idx = int(round(t_truth))
            if 0 <= mf_idx < len(mf_output):
                ax_bot.plot([t_truth * clk_ns], [mf_output[mf_idx]], "^",
                            color=color, ms=10, mec="k", mew=0.8, label=lbl)
        ax_bot.set_xlabel("time (ns)")
        ax_bot.set_ylabel("MF amplitude (ADC)")
        n_mf = len(mf_peaks)
        ax_bot.legend(
            title=f"MF peaks detected: {n_mf}", fontsize=8, loc="upper right"
        )
        ax_bot.grid(True, alpha=0.25)

        fig.tight_layout()

        dt_int    = int(round(dt_ns))
        ratio_int = int(round(ratio * 10000))
        fname = (f"mf_example_dt{dt_int:03d}_ratio{ratio_int:05d}"
                 f"_ev{ev_idx:04d}.png")
        p = out_dir / fname
        fig.savefig(str(p), dpi=130)
        plt.close(fig)
        print(f"[plot] wrote {p}", flush=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:  # noqa: C901
    ap = argparse.ArgumentParser(
        description="Benchmark the matched-filter peak finder on synthetic pile-up.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("synth_dir",
                    help="Directory with pileup_generator.py output "
                         "(manifest.json + .npz).")
    ap.add_argument("--template", default="",
                    help="Pulse template JSON (default: from manifest).")
    ap.add_argument("--out-dir", default="",
                    help="Output directory (default: synth_dir).")
    ap.add_argument("--daq-config", default="",
                    help="DAQ config path.")
    ap.add_argument("--kernel-samples", type=int, default=30,
                    help="Number of samples in the MF kernel (default: 30 = 120 ns).")
    ap.add_argument("--mf-nsigma", type=float, default=5.0,
                    help="MF peak threshold in units of pedestal RMS.")
    ap.add_argument("--min-separation", type=int, default=3,
                    help="Minimum peak separation in samples (default: 3 = 12 ns).")
    ap.add_argument("--min-prominence", type=float, default=0.3,
                    help="Minimum peak prominence as a fraction of the MF threshold. "
                         "Rejects shallow tail-crossing artifacts where the MF output "
                         "barely rises above the local baseline (default: 0.3).")
    ap.add_argument("--match-window", type=int, default=2,
                    help="Truth-matching half-window in samples (default: 2 = ±8 ns).")
    ap.add_argument("--dt-bins", type=int, default=10,
                    help="Number of ΔT bins for random-dt mode summary table.")
    ap.add_argument("--plot", action="store_true",
                    help="Produce comparison heatmaps (mf_vs_wa_detection.png).")
    ap.add_argument("--plot-examples", type=int, default=0,
                    help="Number of example waveform+MF overlay plots to produce "
                         "(0 = none); picks cells where MF detects but WA doesn't.")
    args = ap.parse_args()

    _require_prad2py()

    synth_dir = Path(args.synth_dir).resolve()
    if not synth_dir.is_dir():
        raise SystemExit(f"[ERROR] synth_dir not found: {synth_dir}")

    out_dir = Path(args.out_dir).resolve() if args.out_dir else synth_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- Load manifest ----
    manifest = _load_manifest(synth_dir)
    material   = manifest.get("material", "PbWO4")
    random_dt  = bool(manifest.get("random_dt", False))
    file_entries = manifest.get("files", [])

    print(f"[setup] synth_dir  : {synth_dir}", flush=True)
    print(f"[setup] material   : {material}", flush=True)
    print(f"[setup] random_dt  : {random_dt}", flush=True)
    print(f"[setup] n_files    : {len(file_entries)}", flush=True)

    # ---- Template params (tau_r, tau_f) for MF kernel ----
    tau_r, tau_f = _load_template_params(manifest, args.template)

    # ---- Build MF kernel ----
    kernel, kernel_norm = build_kernel(
        tau_r, tau_f,
        clk_ns=manifest.get("clk_ns", 4.0),
        kernel_samples=args.kernel_samples,
    )
    print(f"[setup] MF kernel  : {args.kernel_samples} samples  "
          f"kernel_norm={kernel_norm:.4f}", flush=True)

    # ---- WaveAnalyzer + DAQ config ----
    wa, clk_ns = _setup_wave_analyzer(args.daq_config)

    # Re-build kernel with DAQ-derived clk_ns (may differ from manifest value
    # if the DAQ config was updated; the WaveAnalyzer clock is authoritative).
    manifest_clk = float(manifest.get("clk_ns", clk_ns))
    if abs(manifest_clk - clk_ns) > 0.01:
        print(f"[setup] NOTE: manifest clk_ns={manifest_clk:.4f} differs from "
              f"DAQ config clk_ns={clk_ns:.4f}; rebuilding kernel.", flush=True)
        kernel, kernel_norm = build_kernel(
            tau_r, tau_f, clk_ns=clk_ns, kernel_samples=args.kernel_samples
        )

    print(f"[setup] MF params  : mf_nsigma={args.mf_nsigma}  "
          f"min_sep={args.min_separation} smp  "
          f"match_window=±{args.match_window} smp", flush=True)

    # ---- Build cell_map skeleton ----
    dt_edges: Optional[np.ndarray] = None
    cell_map: Dict[Tuple[float, float], MFCellStats] = {}

    if random_dt:
        dt_min_ns = float(manifest.get("dt_min_ns", 4.0))
        dt_max_ns = float(manifest.get("dt_max_ns", 200.0))
        n_bins = max(1, args.dt_bins)
        dt_edges = _make_dt_bins(dt_min_ns, dt_max_ns, n_bins)
        ratio_grid: List[float] = sorted(set(
            float(e.get("ratio", 0.0)) for e in file_entries
        ))
        for bin_k in range(n_bins):
            dt_center = round(_bin_center(dt_edges, bin_k), 2)
            for ratio in ratio_grid:
                cell_map[(dt_center, ratio)] = MFCellStats(dt_center, ratio)
        print(f"[setup] random-dt bins: {n_bins}  "
              f"dt_min={dt_min_ns:.1f}  dt_max={dt_max_ns:.1f} ns", flush=True)
    else:
        for entry in file_entries:
            dt_ns  = round(float(entry.get("dt_ns") or 0.0), 2)
            ratio  = round(float(entry.get("ratio", 0.0)), 4)
            cell_map[(dt_ns, ratio)] = MFCellStats(dt_ns, ratio)
        print(f"[setup] grid cells  : {len(cell_map)}", flush=True)

    # ---- Process .npz files ----
    n_total_events = 0
    for entry in file_entries:
        fname    = entry.get("path", "")
        npz_path = synth_dir / fname
        if not npz_path.is_file():
            print(f"[WARN] .npz not found: {npz_path}", file=sys.stderr)
            continue
        print(f"[file] {fname} ...", end="", flush=True)
        n_ev = _process_npz(
            npz_path             = npz_path,
            wa                   = wa,
            kernel               = kernel,
            kernel_norm          = kernel_norm,
            clk_ns               = clk_ns,
            mf_nsigma            = args.mf_nsigma,
            min_separation       = args.min_separation,
            match_window         = args.match_window,
            manifest             = manifest,
            cell_map             = cell_map,
            random_dt            = random_dt,
            dt_edges             = dt_edges,
            collect_examples     = args.plot_examples,
            min_prominence_frac  = args.min_prominence,
        )
        n_total_events += n_ev
        print(f" {n_ev} events", flush=True)

    print(f"\n[summary] total events processed: {n_total_events}", flush=True)

    # ---- Sort cells for output ----
    sorted_cells = sorted(cell_map.values(), key=lambda c: (c.dt_ns, c.ratio))

    # ---- Print comparison table ----
    _print_table(sorted_cells)

    # ---- Global summary ----
    n_all     = sum(c.n_events      for c in sorted_cells)
    n_wa_det  = sum(c.n_wa_detected for c in sorted_cells)
    n_mf_det  = sum(c.n_mf_detected for c in sorted_cells)
    global_summary = {
        "n_total_events":    n_all,
        "n_wa_detected":     n_wa_det,
        "n_mf_detected":     n_mf_det,
        "global_wa_eff_pct": n_wa_det / n_all * 100.0 if n_all > 0 else float("nan"),
        "global_mf_eff_pct": n_mf_det / n_all * 100.0 if n_all > 0 else float("nan"),
        "global_improvement_pct": (
            (n_mf_det - n_wa_det) / n_all * 100.0 if n_all > 0 else float("nan")
        ),
    }

    print(f"[summary] global WA eff: "
          f"{global_summary['global_wa_eff_pct']:.1f}%  "
          f"MF eff: {global_summary['global_mf_eff_pct']:.1f}%  "
          f"improvement: {global_summary['global_improvement_pct']:+.1f}%",
          flush=True)

    # ---- Save results JSON ----
    mf_params = {
        "tau_r_ns":        tau_r,
        "tau_f_ns":        tau_f,
        "kernel_samples":  args.kernel_samples,
        "kernel_norm":     kernel_norm,
        "clk_ns":          clk_ns,
        "mf_nsigma":       args.mf_nsigma,
        "min_separation":  args.min_separation,
        "match_window":    args.match_window,
    }
    result = {
        "synth_dir":  str(synth_dir),
        "material":   material,
        "random_dt":  random_dt,
        "mf_params":  mf_params,
        "cells":      [c.to_dict() for c in sorted_cells],
        "global":     global_summary,
    }
    json_out = out_dir / "mf_benchmark.json"
    json_out.write_text(json.dumps(result, indent=2))
    print(f"[output] {json_out}", flush=True)

    # ---- Optional heatmap plots ----
    if args.plot:
        _make_heatmaps(sorted_cells, out_dir)

    # ---- Optional example waveform plots ----
    if args.plot_examples > 0:
        _plot_examples(out_dir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
