#!/usr/bin/env python3
"""
benchmark_deconv.py — benchmark the LM deconvolver on synthetic pile-up events.

Reads synthetic 2-pulse waveforms produced by ``pileup_generator.py``
(manifest.json + .npz files), runs ``WaveAnalyzer.deconvolve()`` on
every event, compares the recovered peak amplitudes to the stored truth
values, and writes a performance summary table + JSON + optional 2D
heatmap PNGs.

The script covers both fixed-grid (ΔT, ratio) datasets and random-ΔT
datasets.  In random-ΔT mode the per-event ΔT values are binned into
``--dt-bins`` equally-spaced intervals covering the dataset range, and
all statistics are reported per bin.

Usage
-----
    python benchmark_deconv.py <synth_dir> \\
        [--template pulse_templates.json] \\
        [--out-dir results/] \\
        [--daq-config database/daq_config.json] \\
        [--plot] \\
        [--dt-bins 10]

Output files (written to <out-dir>, defaulting to <synth_dir>)
--------------------------------------------------------------
    deconv_benchmark.json       — per-cell stats + global summary (JSON)
    deconv_detection_efficiency.png   — 2D heatmap: WA ≥2 peak fraction
    deconv_amplitude_error.png        — 2D heatmap: median |inj amp err| %

Amplitude-comparison convention (from pileup_generator.py §1a.2)
-----------------------------------------------------------------
Compare ``dec_out.height[k]`` to ``truth_heights_adc[i, col]``.
  col 0  = injected pulse (exact truth)
  col 1  = base-reference pulse (WaveAnalyzer observable — accurate,
            not exact)

Peak-matching convention
------------------------
See the ``_match_peaks`` helper below.  Recovered peak times come from
``wres.peaks[k].time``; recovered heights from ``dec_out.height[k]``.

Edge cases counted per event
----------------------------
  WA found < 2 peaks ("wa_1peak")   → deconvolution skipped entirely.
  dec_out.n == 0                     → "deconv_zero".
  dec_out.n == 1                     → "deconv_merged".
  dec_out.n >= 2                     → "deconv_ok" (2 used, extras discarded).
  dec_out.n > 2                      → additionally counted as "deconv_spurious".
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
# Path helper (mirrors deconv_pileup_demo.py)
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
# Template setup
# ---------------------------------------------------------------------------

def _setup_template(manifest: Dict, template_override: str, wave_cfg) -> object:
    """Load PulseTemplateStore and return the per-type template for the material.

    Parameters
    ----------
    manifest        : parsed manifest.json dict
    template_override : CLI --template override (may be empty string)
    wave_cfg        : dec.WaveConfig (already constructed from daq config)

    Returns
    -------
    tmpl : per-type PulseTemplate object suitable for WaveAnalyzer.deconvolve()
    """
    tmpl_rel = template_override or manifest.get("template_source", "")
    tmpl_path = _resolve_db_path(tmpl_rel)
    if not tmpl_path or not Path(tmpl_path).is_file():
        raise SystemExit(
            f"[ERROR] pulse template file not found: {tmpl_path!r}\n"
            "        Use --template to specify the path explicitly."
        )
    store = _dec.PulseTemplateStore()
    wcfg = _dec.WaveConfig(wave_cfg) if not isinstance(wave_cfg, _dec.WaveConfig) else wave_cfg
    if not store.load_from_file(tmpl_path, wcfg):
        raise SystemExit(
            f"[ERROR] PulseTemplateStore.load_from_file({tmpl_path!r}) failed"
        )
    material = manifest.get("material", "PbWO4")
    tmpl = store.type_template(material)
    if tmpl is None:
        raise SystemExit(
            f"[ERROR] no template for material {material!r} in {tmpl_path}"
        )
    print(f"[setup] templates  : {tmpl_path}  "
          f"({store.n_types_loaded} types,  material={material})", flush=True)
    return tmpl


# ---------------------------------------------------------------------------
# WaveAnalyzer setup
# ---------------------------------------------------------------------------

def _setup_wave_analyzer(daq_config: str):
    """Load DAQ config and return (WaveAnalyzer, clk_ns, daq_cfg)."""
    daq_cfg_path = daq_config or _resolve_db_path("daq_config.json")
    cfg = _dec.load_daq_config(daq_cfg_path) if daq_cfg_path else _dec.load_daq_config()
    print(f"[setup] DAQ config : {daq_cfg_path or '(default)'}", flush=True)
    wcfg = _dec.WaveConfig(cfg.wave_cfg)
    wa = _dec.WaveAnalyzer(wcfg)
    clk_mhz = float(wcfg.clk_mhz) if wcfg.clk_mhz > 0 else 250.0
    clk_ns = 1000.0 / clk_mhz
    print(f"[setup] clk_ns     : {clk_ns:.4f}  ({clk_mhz} MHz)", flush=True)
    return wa, clk_ns, cfg.wave_cfg


# ---------------------------------------------------------------------------
# Peak-matching logic
# ---------------------------------------------------------------------------

def _match_peaks(
    rec_times: List[float],
    rec_heights: List[float],
    truth_inj_t: float,
    truth_base_t: float,
    truth_h_inj: float,
    truth_h_base: float,
) -> Tuple[float, float]:
    """Match recovered peaks to injected / base truth by time proximity.

    Parameters
    ----------
    rec_times    : recovered peak times from wres.peaks (ns)
    rec_heights  : recovered peak heights from dec_out.height
    truth_inj_t  : truth_peak_times_ns[i, 0]  (injected pulse peak time)
    truth_base_t : truth_peak_times_ns[i, 1]  (base pulse peak time)
    truth_h_inj  : truth_heights_adc[i, 0]    (injected pulse height)
    truth_h_base : truth_heights_adc[i, 1]    (base pulse height)

    Returns
    -------
    (inj_err_pct, base_err_pct) — signed percentage errors vs. truth heights.
    Both are NaN if assignment is not possible (e.g. < 2 rec peaks provided).
    """
    nan = float("nan")
    if len(rec_times) < 2 or len(rec_heights) < 2:
        return nan, nan

    # Use first two recovered times for the nearest-truth assignment.
    d0_inj = abs(rec_times[0] - truth_inj_t)
    d1_inj = abs(rec_times[1] - truth_inj_t)
    if d0_inj <= d1_inj:
        inj_idx, base_idx = 0, 1
    else:
        inj_idx, base_idx = 1, 0

    inj_err_pct = (100.0 * (rec_heights[inj_idx] - truth_h_inj) / truth_h_inj
                   if truth_h_inj != 0.0 else nan)
    base_err_pct = (100.0 * (rec_heights[base_idx] - truth_h_base) / truth_h_base
                    if truth_h_base != 0.0 else nan)
    return inj_err_pct, base_err_pct


# ---------------------------------------------------------------------------
# Per-cell accumulator
# ---------------------------------------------------------------------------

class CellStats:
    """Accumulate per-event results for one (ΔT, ratio) grid cell."""

    def __init__(self, dt_ns: float, ratio: float) -> None:
        self.dt_ns = dt_ns
        self.ratio = ratio
        self.n_events: int = 0
        self.n_wa_1peak: int = 0    # WaveAnalyzer found < 2 peaks
        self.n_wa_2peak: int = 0    # WaveAnalyzer found ≥ 2 peaks
        self.n_deconv_zero: int = 0
        self.n_deconv_merged: int = 0   # dec_out.n == 1
        self.n_deconv_spurious: int = 0 # dec_out.n >= 3
        self.n_deconv_ok: int = 0       # dec_out.n >= 2 (may also be spurious)
        self._inj_errs: List[float] = []
        self._base_errs: List[float] = []

    def record_wa_1peak(self) -> None:
        self.n_events += 1
        self.n_wa_1peak += 1

    def record_deconv_zero(self) -> None:
        self.n_events += 1
        self.n_wa_2peak += 1
        self.n_deconv_zero += 1

    def record_deconv_merged(self) -> None:
        self.n_events += 1
        self.n_wa_2peak += 1
        self.n_deconv_merged += 1

    def record_deconv_ok(self, inj_err: float, base_err: float,
                         spurious: bool = False) -> None:
        self.n_events += 1
        self.n_wa_2peak += 1
        self.n_deconv_ok += 1
        if spurious:
            self.n_deconv_spurious += 1
        if math.isfinite(inj_err):
            self._inj_errs.append(inj_err)
        if math.isfinite(base_err):
            self._base_errs.append(base_err)

    # ---- summary accessors -------------------------------------------------

    @property
    def detection_eff(self) -> float:
        """Fraction of events where WaveAnalyzer found ≥ 2 peaks."""
        return self.n_wa_2peak / self.n_events if self.n_events > 0 else float("nan")

    @property
    def deconv_ok_eff(self) -> float:
        """Fraction of events where deconv produced ≥ 2 peaks."""
        return self.n_deconv_ok / self.n_events if self.n_events > 0 else float("nan")

    @property
    def median_inj_err(self) -> float:
        return float(np.median(self._inj_errs)) if self._inj_errs else float("nan")

    @property
    def mean_inj_err(self) -> float:
        return float(np.mean(self._inj_errs)) if self._inj_errs else float("nan")

    @property
    def median_abs_inj_err(self) -> float:
        if not self._inj_errs:
            return float("nan")
        return float(np.median(np.abs(self._inj_errs)))

    @property
    def median_base_err(self) -> float:
        return float(np.median(self._base_errs)) if self._base_errs else float("nan")

    @property
    def mean_base_err(self) -> float:
        return float(np.mean(self._base_errs)) if self._base_errs else float("nan")

    @property
    def median_abs_base_err(self) -> float:
        if not self._base_errs:
            return float("nan")
        return float(np.median(np.abs(self._base_errs)))

    def to_dict(self) -> Dict:
        return {
            "dt_ns":              self.dt_ns,
            "ratio":              self.ratio,
            "n_events":           self.n_events,
            "n_wa_1peak":         self.n_wa_1peak,
            "n_wa_2peak":         self.n_wa_2peak,
            "n_deconv_zero":      self.n_deconv_zero,
            "n_deconv_merged":    self.n_deconv_merged,
            "n_deconv_spurious":  self.n_deconv_spurious,
            "n_deconv_ok":        self.n_deconv_ok,
            "detection_eff":      self.detection_eff,
            "deconv_ok_eff":      self.deconv_ok_eff,
            "median_inj_err_pct": self.median_inj_err,
            "mean_inj_err_pct":   self.mean_inj_err,
            "median_abs_inj_err_pct": self.median_abs_inj_err,
            "median_base_err_pct": self.median_base_err,
            "mean_base_err_pct":  self.mean_base_err,
            "median_abs_base_err_pct": self.median_abs_base_err,
        }


# ---------------------------------------------------------------------------
# ΔT-bin helpers for random-dt mode
# ---------------------------------------------------------------------------

def _make_dt_bins(dt_min: float, dt_max: float, n_bins: int) -> np.ndarray:
    """Return edges of n_bins equal-width ΔT bins spanning [dt_min, dt_max]."""
    return np.linspace(dt_min, dt_max, n_bins + 1)


def _dt_bin_index(dt: float, edges: np.ndarray) -> int:
    """Return 0-based bin index for dt; clips to [0, n_bins-1]."""
    idx = int(np.searchsorted(edges, dt, side="right")) - 1
    return max(0, min(len(edges) - 2, idx))


def _bin_center(edges: np.ndarray, k: int) -> float:
    return float(0.5 * (edges[k] + edges[k + 1]))


# ---------------------------------------------------------------------------
# Core benchmark loop
# ---------------------------------------------------------------------------

def _process_npz(
    npz_path: Path,
    wa: object,
    tmpl: object,
    manifest: Dict,
    cell_map: Dict[Tuple[float, float], CellStats],
    random_dt: bool,
    dt_edges: Optional[np.ndarray],
    ratios: List[float],
    n_warn_limit: int = 5,
) -> int:
    """Process one .npz file.  Returns the number of events processed."""

    try:
        data = np.load(str(npz_path), allow_pickle=False)
    except Exception as exc:
        print(f"[WARN] cannot load {npz_path.name}: {exc}", file=sys.stderr)
        return 0

    waveforms = data["waveforms"]                        # uint16 (N, n_samples)
    truth_h = data["truth_heights_adc"]                  # float32 (N, 2)
    truth_peak_times_ns = data["truth_peak_times_ns"]    # float32 (N, 2)
    event_dt_ns_arr = data["event_dt_ns"]                # float32 (N,)

    meta_ratio = round(float(data["meta_ratio"]), 4)

    N = waveforms.shape[0]
    n_warn = 0

    for i in range(N):
        samples = waveforms[i]          # uint16 (n_samples,)
        dt_ev = float(event_dt_ns_arr[i])

        # Determine which CellStats bucket this event belongs to.
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
            # Shouldn't happen — cell_map is pre-populated before this call.
            if n_warn < n_warn_limit:
                print(f"[WARN] no cell for key {cell_key}; skipping event {i}",
                      file=sys.stderr)
                n_warn += 1
            continue

        truth_inj_t = float(truth_peak_times_ns[i, 0])
        truth_base_t = float(truth_peak_times_ns[i, 1])
        truth_h_inj = float(truth_h[i, 0])
        truth_h_base = float(truth_h[i, 1])

        # ---- Step 1: WaveAnalyzer analyze_result ----
        try:
            wres = wa.analyze_result(samples)
        except Exception as exc:
            if n_warn < n_warn_limit:
                print(f"[WARN] analyze_result failed for event {i}: {exc}",
                      file=sys.stderr)
                n_warn += 1
            cell.record_wa_1peak()
            continue

        n_peaks_found = len(list(wres.peaks))

        if n_peaks_found < 2:
            cell.record_wa_1peak()
            continue

        # ---- Step 2: deconvolve ----
        try:
            out = wa.deconvolve(samples, wres, tmpl)
        except Exception as exc:
            if n_warn < n_warn_limit:
                print(f"[WARN] deconvolve failed for event {i}: {exc}",
                      file=sys.stderr)
                n_warn += 1
            cell.record_deconv_zero()
            continue

        n_deconv = out.n

        if n_deconv == 0:
            cell.record_deconv_zero()
            continue

        if n_deconv == 1:
            cell.record_deconv_merged()
            continue

        # n_deconv >= 2: extract recovered peak times + heights.
        # Guard: deconv should never produce more peaks than WA found,
        # but clamp to be safe.
        n_use = min(n_deconv, n_peaks_found, len(list(wres.peaks)))
        rec_times = [list(wres.peaks)[k].time for k in range(n_use)]
        rec_heights = [out.height[k] for k in range(n_deconv)]

        inj_err, base_err = _match_peaks(
            rec_times, rec_heights,
            truth_inj_t, truth_base_t,
            truth_h_inj, truth_h_base,
        )

        spurious = n_deconv >= 3
        cell.record_deconv_ok(inj_err, base_err, spurious=spurious)

    return N


# ---------------------------------------------------------------------------
# Summary table printer
# ---------------------------------------------------------------------------

def _print_table(cells: List[CellStats]) -> None:
    """Print a fixed-width ASCII summary table to stdout."""
    print()
    header = (f"{'ΔT':>6}  {'ratio':>6}  {'N':>5}  "
              f"{'WA_1pk':>6}  {'WA_2pk':>6}  {'dc_ok':>5}  "
              f"{'inj_err%':>8}  {'base_err%':>9}")
    sep = (f"{'-----':>6}  {'------':>6}  {'---':>5}  "
           f"{'------':>6}  {'------':>6}  {'-----':>5}  "
           f"{'--------':>8}  {'---------':>9}")
    print(header)
    print(sep)
    for c in cells:
        dt_str = f"{c.dt_ns:6.1f}" if math.isfinite(c.dt_ns) else "   rnd"
        ratio_str = f"{c.ratio:6.3f}"
        n_str = f"{c.n_events:5d}"
        wa1_str = f"{c.n_wa_1peak:6d}"
        wa2_str = f"{c.n_wa_2peak:6d}"
        dc_str = f"{c.n_deconv_ok:5d}"
        inj = c.median_inj_err
        base = c.median_base_err
        inj_str = f"{inj:+.1f}%" if math.isfinite(inj) else "     NaN"
        base_str = f"{base:+.1f}%" if math.isfinite(base) else "      NaN"
        print(f"{dt_str}  {ratio_str}  {n_str}  "
              f"{wa1_str}  {wa2_str}  {dc_str}  "
              f"{inj_str:>8}  {base_str:>9}")
    print()


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _make_heatmaps(cells: List[CellStats], out_dir: Path) -> None:
    """Produce and save two 2D heatmap PNGs to out_dir."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors

    # Gather unique sorted axes.
    dt_vals = sorted({c.dt_ns for c in cells if math.isfinite(c.dt_ns)})
    ratio_vals = sorted({c.ratio for c in cells})

    if not dt_vals or not ratio_vals:
        print("[plot] insufficient data for heatmaps (no finite ΔT values).",
              file=sys.stderr)
        return

    n_dt = len(dt_vals)
    n_ratio = len(ratio_vals)

    # Build lookup for fast indexing.
    dt_idx = {v: k for k, v in enumerate(dt_vals)}
    r_idx = {v: k for k, v in enumerate(ratio_vals)}

    det_eff = np.full((n_ratio, n_dt), float("nan"))
    amp_err = np.full((n_ratio, n_dt), float("nan"))

    for c in cells:
        if not math.isfinite(c.dt_ns):
            continue
        ki = dt_idx.get(c.dt_ns)
        ri = r_idx.get(c.ratio)
        if ki is None or ri is None:
            continue
        det_eff[ri, ki] = c.detection_eff
        amp_err[ri, ki] = c.median_abs_inj_err

    # ---- Plot 1: detection efficiency ----
    fig, ax = plt.subplots(figsize=(max(6, n_dt * 0.9), max(4, n_ratio * 0.7 + 1.2)))
    img = ax.imshow(det_eff, aspect="auto", origin="lower",
                    vmin=0.0, vmax=1.0, cmap="viridis",
                    extent=[-0.5, n_dt - 0.5, -0.5, n_ratio - 0.5])
    ax.set_xticks(range(n_dt))
    ax.set_xticklabels([f"{v:.0f}" for v in dt_vals], rotation=45, ha="right")
    ax.set_yticks(range(n_ratio))
    ax.set_yticklabels([f"{v:.2g}" for v in ratio_vals])
    ax.set_xlabel("ΔT (ns)")
    ax.set_ylabel("Amplitude ratio (inj / base)")
    ax.set_title("Deconvolution detection efficiency\n"
                 "(fraction of events: WaveAnalyzer ≥ 2 peaks)")
    cbar = fig.colorbar(img, ax=ax, shrink=0.85)
    cbar.set_label("Detection efficiency")
    fig.tight_layout()
    p1 = out_dir / "deconv_detection_efficiency.png"
    fig.savefig(str(p1), dpi=130)
    plt.close(fig)
    print(f"[plot] wrote {p1}", flush=True)

    # ---- Plot 2: median |injected amplitude error| % ----
    # Use log colour scale if the error range spans more than one decade.
    finite_err = amp_err[np.isfinite(amp_err)]
    use_log = (finite_err.size > 0
               and finite_err.max() > 0
               and finite_err.min() > 0
               and finite_err.max() / finite_err.min() > 10.0)

    fig, ax = plt.subplots(figsize=(max(6, n_dt * 0.9), max(4, n_ratio * 0.7 + 1.2)))
    if use_log:
        norm = mcolors.LogNorm(vmin=max(finite_err.min(), 1e-2),
                               vmax=finite_err.max())
    else:
        norm = mcolors.Normalize(vmin=0.0,
                                 vmax=finite_err.max() if finite_err.size > 0 else 1.0)

    img2 = ax.imshow(amp_err, aspect="auto", origin="lower",
                     norm=norm, cmap="plasma",
                     extent=[-0.5, n_dt - 0.5, -0.5, n_ratio - 0.5])
    ax.set_xticks(range(n_dt))
    ax.set_xticklabels([f"{v:.0f}" for v in dt_vals], rotation=45, ha="right")
    ax.set_yticks(range(n_ratio))
    ax.set_yticklabels([f"{v:.2g}" for v in ratio_vals])
    ax.set_xlabel("ΔT (ns)")
    ax.set_ylabel("Amplitude ratio (inj / base)")
    scale_note = " (log scale)" if use_log else ""
    ax.set_title(f"Median |injected amplitude error|{scale_note}\n"
                 "(events with deconvolution ≥ 2 peaks, injected pulse)")
    cbar2 = fig.colorbar(img2, ax=ax, shrink=0.85)
    cbar2.set_label("Median |amplitude error| (%)")
    fig.tight_layout()
    p2 = out_dir / "deconv_amplitude_error.png"
    fig.savefig(str(p2), dpi=130)
    plt.close(fig)
    print(f"[plot] wrote {p2}", flush=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:  # noqa: C901
    ap = argparse.ArgumentParser(
        description="Benchmark the LM deconvolver on synthetic pile-up events.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("synth_dir",
                    help="Directory containing pileup_generator.py output "
                         "(manifest.json + .npz files).")
    ap.add_argument("--template",
                    help="Pulse template JSON (default: read from manifest).")
    ap.add_argument("--out-dir", default="",
                    help="Output directory for plots/tables (default: synth_dir).")
    ap.add_argument("--daq-config", default="",
                    help="DAQ config path.")
    ap.add_argument("--plot", action="store_true",
                    help="Produce diagnostic plots (2D heatmaps).")
    ap.add_argument("--dt-bins", type=int, default=10,
                    help="Number of ΔT bins for random-dt mode summary table.")
    args = ap.parse_args()

    _require_prad2py()

    synth_dir = Path(args.synth_dir).resolve()
    if not synth_dir.is_dir():
        raise SystemExit(f"[ERROR] synth_dir not found: {synth_dir}")

    out_dir = Path(args.out_dir).resolve() if args.out_dir else synth_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- Load manifest ----
    manifest = _load_manifest(synth_dir)
    material = manifest.get("material", "PbWO4")
    random_dt = bool(manifest.get("random_dt", False))
    file_entries = manifest.get("files", [])

    print(f"[setup] synth_dir  : {synth_dir}", flush=True)
    print(f"[setup] material   : {material}", flush=True)
    print(f"[setup] random_dt  : {random_dt}", flush=True)
    print(f"[setup] n_files    : {len(file_entries)}", flush=True)

    # ---- WaveAnalyzer + DAQ config ----
    wa, clk_ns, wave_cfg = _setup_wave_analyzer(args.daq_config)

    # ---- Template ----
    tmpl = _setup_template(manifest, args.template or "", wave_cfg)

    # ---- Build cell_map skeleton ----
    # For fixed-grid mode: one cell per (dt_ns, ratio) entry in the manifest.
    # For random-dt mode:  bin edges span [dt_min_ns, dt_max_ns]; one row per
    #                      (bin_center, ratio).
    dt_edges: Optional[np.ndarray] = None
    cell_map: Dict[Tuple[float, float], CellStats] = {}

    if random_dt:
        dt_min_ns = float(manifest.get("dt_min_ns", 4.0))
        dt_max_ns = float(manifest.get("dt_max_ns", 200.0))
        n_bins = max(1, args.dt_bins)
        dt_edges = _make_dt_bins(dt_min_ns, dt_max_ns, n_bins)
        ratio_grid: List[float] = sorted(set(float(e.get("ratio", 0.0))
                                              for e in file_entries))
        for bin_k in range(n_bins):
            dt_center = round(_bin_center(dt_edges, bin_k), 2)
            for ratio in ratio_grid:
                cell_map[(dt_center, ratio)] = CellStats(dt_center, ratio)
        print(f"[setup] random-dt bins: {n_bins}  "
              f"dt_min={dt_min_ns:.1f}  dt_max={dt_max_ns:.1f} ns", flush=True)
    else:
        for entry in file_entries:
            dt_ns = round(float(entry.get("dt_ns") or 0.0), 2)
            ratio = round(float(entry.get("ratio", 0.0)), 4)
            cell_map[(dt_ns, ratio)] = CellStats(dt_ns, ratio)
        print(f"[setup] grid cells  : {len(cell_map)}", flush=True)

    # ---- Process .npz files ----
    n_total_events = 0
    for entry in file_entries:
        fname = entry.get("path", "")
        npz_path = synth_dir / fname
        if not npz_path.is_file():
            print(f"[WARN] .npz not found: {npz_path}", file=sys.stderr)
            continue
        print(f"[file] {fname} ...", end="", flush=True)
        n_ev = _process_npz(
            npz_path, wa, tmpl, manifest,
            cell_map, random_dt, dt_edges,
            list({float(e.get("ratio", 0.0)) for e in file_entries}),
        )
        n_total_events += n_ev
        print(f" {n_ev} events", flush=True)

    print(f"\n[summary] total events processed: {n_total_events}", flush=True)

    # ---- Sort cells for table output ----
    sorted_cells = sorted(cell_map.values(), key=lambda c: (c.dt_ns, c.ratio))

    # ---- Print summary table ----
    _print_table(sorted_cells)

    # ---- Global summary statistics ----
    n_all = sum(c.n_events for c in sorted_cells)
    n_ok = sum(c.n_deconv_ok for c in sorted_cells)
    n_wa2 = sum(c.n_wa_2peak for c in sorted_cells)
    global_summary = {
        "n_total_events": n_all,
        "n_wa_2peak_total": n_wa2,
        "n_deconv_ok_total": n_ok,
        "global_detection_eff": n_wa2 / n_all if n_all > 0 else float("nan"),
        "global_deconv_ok_eff": n_ok / n_all if n_all > 0 else float("nan"),
    }

    # ---- Save JSON ----
    result = {
        "synth_dir": str(synth_dir),
        "material":  material,
        "random_dt": random_dt,
        "cells": [c.to_dict() for c in sorted_cells],
        "global": global_summary,
    }
    json_out = out_dir / "deconv_benchmark.json"
    json_out.write_text(json.dumps(result, indent=2))
    print(f"[output] {json_out}", flush=True)

    # ---- Optional heatmap plots ----
    if args.plot:
        _make_heatmaps(sorted_cells, out_dir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
