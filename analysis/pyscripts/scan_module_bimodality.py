#!/usr/bin/env python3
"""scan_module_bimodality.py — Detect bimodal pre-calibration energy spectra
across PbWO4 modules in a clusters TSV.

Reads a clusters TSV produced by benchmark_hycal_timing.py --write-cluster-tsv,
filters to config_name == 'legacy' and PbWO4 modules (seed_module_id starts
with 'W'), then for each module with sufficient statistics detects whether its
raw energy spectrum is bimodal using scipy.signal.find_peaks with Gaussian
peak fitting.

Outputs:
  <output-dir>/module_bimodality_<run-label>.tsv       — per-module summary
  <output-dir>/module_bimodality_<run-label>_map.png   — detector face map
  <output-dir>/module_bimodality_<run-label>_spatial.png — spatial analysis
  <output-dir>/module_bimodality_<run-label>_top20.pdf — worst offenders
  <output-dir>/module_bimodality_<run-label>_findings.md — findings summary
"""

from __future__ import annotations

import argparse
import math
import sys
from datetime import date
from pathlib import Path
from typing import Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

try:
    import pandas as pd
    HAVE_PANDAS = True
except ImportError:
    HAVE_PANDAS = False

try:
    from scipy.signal import find_peaks
    from scipy.optimize import curve_fit
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

# ---------------------------------------------------------------------------
# Script-level constants
# ---------------------------------------------------------------------------

_SCRIPT_DIR = Path(__file__).resolve().parent
_NAN = float("nan")


# ---------------------------------------------------------------------------
# Gaussian helper
# ---------------------------------------------------------------------------

def _gauss(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


# ---------------------------------------------------------------------------
# Per-module peak detection and Gaussian fitting
# ---------------------------------------------------------------------------

def _fit_gaussian_window(centres, counts, peak_bin: int,
                          bin_width: float, window_mev: float = 100.0):
    """Fit a Gaussian in ±window_mev around peak_bin.

    Returns (mu, sigma, A) or (NaN, NaN, NaN) on failure.
    """
    mu0 = float(centres[peak_bin])
    mask = np.abs(centres - mu0) <= window_mev
    x = centres[mask]
    y = counts[mask].astype(float)
    if y.sum() < 5 or len(x) < 3:
        return _NAN, _NAN, _NAN
    try:
        p0 = [float(y.max()), mu0, max(bin_width * 2, 40.0)]
        lo = centres[0] - bin_width
        hi = centres[-1] + bin_width
        bounds = ([0.0, lo, 1.0], [np.inf, hi, window_mev * 2])
        popt, _ = curve_fit(_gauss, x, y, p0=p0, bounds=bounds, maxfev=5000)
        return float(popt[1]), abs(float(popt[2])), float(popt[0])
    except Exception:
        return _NAN, _NAN, _NAN


def _estimate_sigma_fwhm(centres, counts, peak_bin: int):
    """Estimate sigma from FWHM: sigma = FWHM / 2.355."""
    peak_amp = float(counts[peak_bin])
    half_max = peak_amp / 2.0
    # Walk left
    left = peak_bin
    while left > 0 and counts[left] > half_max:
        left -= 1
    # Walk right
    right = peak_bin
    while right < len(counts) - 1 and counts[right] > half_max:
        right += 1
    fwhm = float(centres[right] - centres[left])
    if fwhm <= 0:
        fwhm = float(centres[1] - centres[0]) * 3  # rough fallback
    return fwhm / 2.355


def analyze_module(raw_vals: np.ndarray, e_lo: float, e_hi: float,
                   bin_width: float, prominence_frac: float,
                   bimodal_ratio: float) -> dict:
    """Detect peaks and classify bimodality for one module.

    Returns a dict with keys:
        n_peaks_detected, peak_1_mu, peak_1_sigma, peak_1_amp,
        peak_2_mu, peak_2_sigma, peak_2_amp, peak_2_over_peak_1_amp, is_bimodal
    """
    edges = np.arange(e_lo, e_hi + bin_width, bin_width)
    raw_counts, _ = np.histogram(raw_vals, bins=edges)
    centres = 0.5 * (edges[:-1] + edges[1:])

    # Light smoothing: 3-bin running average
    kernel = np.ones(3) / 3.0
    smooth = np.convolve(raw_counts.astype(float), kernel, mode="same")

    max_smooth = float(smooth.max()) if smooth.max() > 0 else 1.0
    prominence = max(5.0, prominence_frac * max_smooth)
    min_distance = max(1, int(100.0 / bin_width))  # 100 MeV separation

    result = dict(
        n_peaks_detected=0,
        peak_1_mu_MeV=_NAN, peak_1_sigma_MeV=_NAN, peak_1_amp=_NAN,
        peak_2_mu_MeV=_NAN, peak_2_sigma_MeV=_NAN, peak_2_amp=_NAN,
        peak_2_over_peak_1_amp=_NAN,
        is_bimodal="N",
    )

    if not HAVE_SCIPY:
        return result

    peak_bins, properties = find_peaks(smooth,
                                        prominence=prominence,
                                        distance=min_distance)

    if len(peak_bins) == 0:
        return result

    # For each detected peak, fit a Gaussian in ±100 MeV
    fitted_peaks = []  # list of (mu, sigma, amp)
    for pb in peak_bins:
        mu, sigma, amp = _fit_gaussian_window(centres, smooth, pb, bin_width,
                                               window_mev=100.0)
        if math.isnan(mu):
            # Fallback: use histogram argmax position + FWHM sigma
            mu = float(centres[pb])
            sigma = _estimate_sigma_fwhm(centres, smooth, pb)
            amp = float(smooth[pb])
        fitted_peaks.append((mu, sigma, amp))

    # Sort by amplitude descending
    fitted_peaks.sort(key=lambda t: t[2], reverse=True)

    result["n_peaks_detected"] = len(fitted_peaks)

    # Peak 1 (dominant)
    mu1, sig1, amp1 = fitted_peaks[0]
    result["peak_1_mu_MeV"] = mu1
    result["peak_1_sigma_MeV"] = sig1
    result["peak_1_amp"] = amp1

    if len(fitted_peaks) >= 2:
        mu2, sig2, amp2 = fitted_peaks[1]
        result["peak_2_mu_MeV"] = mu2
        result["peak_2_sigma_MeV"] = sig2
        result["peak_2_amp"] = amp2
        ratio = amp2 / amp1 if amp1 > 0 else 0.0
        result["peak_2_over_peak_1_amp"] = ratio
        if ratio >= bimodal_ratio:
            result["is_bimodal"] = "Y"

    return result


# ---------------------------------------------------------------------------
# CLI argument parser
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--input-tsv",
        default="output/per_module_025320_split00010_clusters.tsv",
        help="Path to clusters TSV (default: %(default)s)",
    )
    ap.add_argument(
        "--output-dir",
        default="output/",
        help="Output directory for all deliverables (default: %(default)s)",
    )
    ap.add_argument(
        "--run-label",
        default="025320_split00010",
        help="Short label for output filenames (default: %(default)s)",
    )
    ap.add_argument(
        "--min-events", type=int, default=100,
        help="Minimum events per module (default: %(default)s)",
    )
    ap.add_argument(
        "--bimodal-ratio", type=float, default=0.30,
        help="Minimum peak_2/peak_1 amplitude ratio to flag as bimodal "
             "(default: %(default)s)",
    )
    ap.add_argument(
        "--prominence-frac", type=float, default=0.03,
        help="Prominence as fraction of max counts (default: %(default)s, minimum 5)",
    )
    ap.add_argument(
        "--top-n-pdf", type=int, default=20,
        help="Number of worst offenders in the detailed PDF (default: %(default)s)",
    )
    ap.add_argument(
        "--e-range", type=float, nargs=2, default=[500.0, 3500.0],
        metavar=("LO", "HI"),
        help="Energy range in MeV (default: 500 3500)",
    )
    ap.add_argument(
        "--bin-width", type=float, default=20.0,
        help="Bin width in MeV (default: %(default)s)",
    )
    ap.add_argument(
        "--target-peak", type=float, default=2239.5,
        help="Target peak position in MeV (default: %(default)s)",
    )
    return ap


# ---------------------------------------------------------------------------
# Path resolution (mirrors pattern in plot_per_module_pre_post.py)
# ---------------------------------------------------------------------------

def _resolve_input(raw: str) -> Path:
    p = Path(raw)
    if p.is_absolute():
        return p
    # Try CWD first, then script dir
    for base in (Path.cwd(), _SCRIPT_DIR):
        candidate = base / p
        if candidate.exists():
            return candidate
    return Path.cwd() / p  # return so error message is clear


def _resolve_output_dir(raw: str) -> Path:
    p = Path(raw)
    if p.is_absolute():
        return p
    # Prefer script-dir/output if that matches the default pattern
    candidate = _SCRIPT_DIR / p
    return candidate


# ---------------------------------------------------------------------------
# TSV reading
# ---------------------------------------------------------------------------

def read_tsv(input_path: Path, min_events: int) -> tuple:
    """Read and filter the clusters TSV.

    Returns (df_full, module_data_dict, skipped_modules_positions)
      where module_data_dict maps module_name -> dict with fields:
        raw_vals, seed_row, seed_col, seed_x_mm, seed_y_mm, n_events
      and skipped_modules_positions is a list of (x_mm, y_mm) for low-stat modules.
    """
    print("Reading TSV...", flush=True)

    if not HAVE_PANDAS:
        sys.stderr.write("[ERROR] pandas is required. Install with: pip install pandas\n")
        sys.exit(1)

    try:
        df = pd.read_csv(str(input_path), sep="\t", dtype=str, low_memory=False)
    except Exception as exc:
        sys.stderr.write(f"[ERROR] Failed to read TSV: {exc}\n")
        sys.exit(1)

    if df.empty:
        sys.stderr.write("[ERROR] TSV is empty.\n")
        sys.exit(1)

    # Filter to legacy config and PbWO4
    if "config_name" in df.columns:
        df = df[df["config_name"] == "legacy"].copy()
    if "seed_module_id" in df.columns:
        df = df[df["seed_module_id"].str.startswith("W", na=False)].copy()

    if df.empty:
        sys.stderr.write("[ERROR] No PbWO4 (W-prefix) legacy rows found in TSV.\n")
        sys.exit(1)

    # Convert numeric columns
    for col in ("cluster_E_raw_MeV", "cluster_E_cal_MeV",
                "seed_row", "seed_col", "seed_x_mm", "seed_y_mm",
                "seed_energy_MeV"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    unique_modules = df["seed_module_id"].unique()
    n_all = len(unique_modules)

    module_data: dict[str, dict] = {}
    skipped_positions: list[tuple[float, float]] = []

    for mid in unique_modules:
        sub = df[df["seed_module_id"] == mid]
        raw_vals = sub["cluster_E_raw_MeV"].dropna().values.astype(float)
        n = len(raw_vals)

        # Position (mode — all events for a given module should have same x,y)
        x_mm = float(sub["seed_x_mm"].dropna().mode().iloc[0]) \
               if not sub["seed_x_mm"].dropna().empty else _NAN
        y_mm = float(sub["seed_y_mm"].dropna().mode().iloc[0]) \
               if not sub["seed_y_mm"].dropna().empty else _NAN
        row = int(sub["seed_row"].dropna().mode().iloc[0]) \
              if not sub["seed_row"].dropna().empty else -1
        col = int(sub["seed_col"].dropna().mode().iloc[0]) \
              if not sub["seed_col"].dropna().empty else -1

        if n < min_events:
            skipped_positions.append((x_mm, y_mm))
            continue

        module_data[mid] = dict(
            raw_vals=raw_vals,
            seed_row=row,
            seed_col=col,
            seed_x_mm=x_mm,
            seed_y_mm=y_mm,
            n_events=n,
        )

    n_analyzed = len(module_data)
    print(f"Found {n_all} PbWO4 modules, {n_analyzed} with ≥{min_events} events",
          flush=True)

    if n_analyzed == 0:
        sys.stderr.write(
            f"[ERROR] No PbWO4 modules with ≥{min_events} events found.\n"
        )
        sys.exit(1)

    return df, module_data, skipped_positions


# ---------------------------------------------------------------------------
# Main analysis loop
# ---------------------------------------------------------------------------

def run_analysis(module_data: dict, e_lo: float, e_hi: float,
                 bin_width: float, prominence_frac: float,
                 bimodal_ratio: float) -> list[dict]:
    """Analyze all modules. Returns list of result dicts, one per module."""
    print("Analyzing...", flush=True)

    results = []
    for mid, mdata in module_data.items():
        raw_vals = mdata["raw_vals"]
        x_mm = mdata["seed_x_mm"]
        y_mm = mdata["seed_y_mm"]
        d_bh = math.sqrt(x_mm**2 + y_mm**2) if not (math.isnan(x_mm) or math.isnan(y_mm)) else _NAN

        res = analyze_module(raw_vals, e_lo, e_hi, bin_width,
                             prominence_frac, bimodal_ratio)

        row_dict = dict(
            module_name=mid,
            seed_row=mdata["seed_row"],
            seed_col=mdata["seed_col"],
            seed_x_mm=x_mm,
            seed_y_mm=y_mm,
            n_events=mdata["n_events"],
            distance_from_beam_hole_mm=d_bh,
            n_peaks_detected=res["n_peaks_detected"],
            peak_1_mu_MeV=res["peak_1_mu_MeV"],
            peak_1_sigma_MeV=res["peak_1_sigma_MeV"],
            peak_1_amp=res["peak_1_amp"],
            peak_2_mu_MeV=res["peak_2_mu_MeV"],
            peak_2_sigma_MeV=res["peak_2_sigma_MeV"],
            peak_2_amp=res["peak_2_amp"],
            peak_2_over_peak_1_amp=res["peak_2_over_peak_1_amp"],
            is_bimodal=res["is_bimodal"],
            # Keep raw_vals reference for plotting
            _raw_vals=mdata["raw_vals"],
        )
        results.append(row_dict)

        if res["is_bimodal"] == "Y":
            mu1_str = f"{res['peak_1_mu_MeV']:.0f}" if not math.isnan(res["peak_1_mu_MeV"]) else "N/A"
            mu2_str = f"{res['peak_2_mu_MeV']:.0f}" if not math.isnan(res["peak_2_mu_MeV"]) else "N/A"
            ratio_str = f"{res['peak_2_over_peak_1_amp']:.2f}" if not math.isnan(res["peak_2_over_peak_1_amp"]) else "N/A"
            d_str = f"{d_bh:.1f}" if not math.isnan(d_bh) else "N/A"
            print(f"  bimodal: {mid} (row={mdata['seed_row']}, col={mdata['seed_col']}, "
                  f"d_bh={d_str} mm) ratio={ratio_str} "
                  f"μ1={mu1_str} μ2={mu2_str}", flush=True)

    # Sort: bimodal first (Y > N), then by ratio descending
    results.sort(
        key=lambda r: (
            0 if r["is_bimodal"] == "Y" else 1,
            -(r["peak_2_over_peak_1_amp"] if not math.isnan(r["peak_2_over_peak_1_amp"]) else 0.0),
        )
    )

    n_bimodal = sum(1 for r in results if r["is_bimodal"] == "Y")
    n_total = len(results)
    pct = 100.0 * n_bimodal / n_total if n_total else 0.0
    print(f"Done: {n_bimodal} bimodal modules out of {n_total} analyzed "
          f"({n_bimodal}/{n_total} = {pct:.1f}%)", flush=True)

    return results


# ---------------------------------------------------------------------------
# Output 1: Summary TSV
# ---------------------------------------------------------------------------

_TSV_COLS = [
    "module_name", "seed_row", "seed_col", "seed_x_mm", "seed_y_mm",
    "n_events", "distance_from_beam_hole_mm", "n_peaks_detected",
    "peak_1_mu_MeV", "peak_1_sigma_MeV", "peak_1_amp",
    "peak_2_mu_MeV", "peak_2_sigma_MeV", "peak_2_amp",
    "peak_2_over_peak_1_amp", "is_bimodal",
]


def write_summary_tsv(results: list[dict], out_path: Path) -> None:
    with open(str(out_path), "w", encoding="utf-8") as fh:
        fh.write("\t".join(_TSV_COLS) + "\n")
        for r in results:
            vals = []
            for col in _TSV_COLS:
                v = r.get(col, "")
                if isinstance(v, float):
                    if math.isnan(v):
                        vals.append("")
                    else:
                        # Format floats cleanly
                        vals.append(f"{v:.4f}" if col not in ("seed_row", "seed_col", "n_events", "n_peaks_detected") else str(int(v)))
                else:
                    vals.append(str(v) if v is not None else "")
            fh.write("\t".join(vals) + "\n")
    print(f"Wrote: {out_path}", flush=True)


# ---------------------------------------------------------------------------
# Beam-hole exclusion geometry helpers
# ---------------------------------------------------------------------------

def _compute_exclusion_box_mm(results: list[dict],
                               excl_rows=(15, 20), excl_cols=(15, 20)):
    """Estimate the exclusion rectangle in mm from row/col boundaries.

    Finds modules at the boundary rows/cols and uses their x/y positions
    to determine the box extent.  Falls back to a half-bin offset if the
    exact boundary modules are not in the dataset.
    """
    # Build a lookup: (row, col) -> (x_mm, y_mm)
    pos_by_rowcol: dict[tuple, tuple] = {}
    for r in results:
        key = (r["seed_row"], r["seed_col"])
        pos_by_rowcol[key] = (r["seed_x_mm"], r["seed_y_mm"])

    # We want the box from row=excl_rows[0] to row=excl_rows[1]
    # and col=excl_cols[0] to col=excl_cols[1].
    # Modules just outside the hole: rows 14 and 21, cols 14 and 21
    # Find the half-width of one module from adjacent modules
    # Approximate: find any pair of adjacent (in col) modules
    x_vals_at_cols: list[float] = []
    for row in range(33):
        for col in range(excl_cols[0], excl_cols[1] + 1):
            if (row, col) in pos_by_rowcol:
                x_vals_at_cols.append(pos_by_rowcol[(row, col)][0])
    y_vals_at_rows: list[float] = []
    for col in range(33):
        for row in range(excl_rows[0], excl_rows[1] + 1):
            if (row, col) in pos_by_rowcol:
                y_vals_at_rows.append(pos_by_rowcol[(row, col)][1])

    # Estimate module half-size from the step between adjacent cols
    module_half = 10.25  # 20.5 mm / 2

    # x range: find modules at col=excl_cols[0] and col=excl_cols[1]
    x_lo_candidates = [pos_by_rowcol[(r, excl_cols[0])][0]
                        for r in range(33) if (r, excl_cols[0]) in pos_by_rowcol]
    x_hi_candidates = [pos_by_rowcol[(r, excl_cols[1])][0]
                        for r in range(33) if (r, excl_cols[1]) in pos_by_rowcol]
    y_lo_candidates = [pos_by_rowcol[(excl_rows[0], c)][1]
                        for c in range(33) if (excl_rows[0], c) in pos_by_rowcol]
    y_hi_candidates = [pos_by_rowcol[(excl_rows[1], c)][1]
                        for c in range(33) if (excl_rows[1], c) in pos_by_rowcol]

    if x_lo_candidates:
        x_lo = float(np.median(x_lo_candidates)) - module_half
    else:
        x_lo = excl_cols[0] * 20.5 - module_half
    if x_hi_candidates:
        x_hi = float(np.median(x_hi_candidates)) + module_half
    else:
        x_hi = excl_cols[1] * 20.5 + module_half
    if y_lo_candidates:
        # rows go top-to-bottom or bottom-to-top; take the min/max
        y_lo = float(np.median(y_lo_candidates)) - module_half
    else:
        y_lo = excl_rows[0] * 20.5 - module_half
    if y_hi_candidates:
        y_hi = float(np.median(y_hi_candidates)) + module_half
    else:
        y_hi = excl_rows[1] * 20.5 + module_half

    return x_lo, x_hi, y_lo, y_hi


# ---------------------------------------------------------------------------
# Output 2: Detector face map PNG
# ---------------------------------------------------------------------------

def write_map_png(results: list[dict], skipped_positions: list[tuple],
                   out_path: Path, run_label: str,
                   bimodal_ratio: float) -> None:
    fig, ax = plt.subplots(figsize=(10, 10))

    # --- Colour coding ---
    x_grey, y_grey = [], []
    x_green, y_green = [], []
    x_yellow, y_yellow = [], []
    x_red, y_red = [], []

    for x_mm, y_mm in skipped_positions:
        if not (math.isnan(x_mm) or math.isnan(y_mm)):
            x_grey.append(x_mm)
            y_grey.append(y_mm)

    for r in results:
        x_mm = r["seed_x_mm"]
        y_mm = r["seed_y_mm"]
        if math.isnan(x_mm) or math.isnan(y_mm):
            continue
        ratio = r["peak_2_over_peak_1_amp"]
        if r["is_bimodal"] == "N":
            x_green.append(x_mm)
            y_green.append(y_mm)
        elif not math.isnan(ratio) and ratio >= 0.50:
            x_red.append(x_mm)
            y_red.append(y_mm)
        else:
            x_yellow.append(x_mm)
            y_yellow.append(y_mm)

    ms = 150  # marker size
    ax.scatter(x_grey, y_grey, s=30, c="grey", alpha=0.3, label="Low-stat (<100 ev)")
    ax.scatter(x_green, y_green, s=ms, c="green", alpha=0.7, label="Unimodal")
    ax.scatter(x_yellow, y_yellow, s=ms, c="gold", alpha=0.85,
               edgecolors="darkorange", linewidths=0.8, label="Mild bimodal (0.30–0.50)")
    ax.scatter(x_red, y_red, s=ms, c="red", alpha=0.85,
               edgecolors="darkred", linewidths=0.8, label="Strong bimodal (≥0.50)")

    # --- Beam-hole exclusion box ---
    x_lo, x_hi, y_lo, y_hi = _compute_exclusion_box_mm(results)
    box_x = [x_lo, x_hi, x_hi, x_lo, x_lo]
    box_y = [y_lo, y_lo, y_hi, y_hi, y_lo]
    ax.plot(box_x, box_y, "k--", lw=1.5, label="Beam-hole excl. (rows/cols 15–20)")

    # --- Inner-ring circle r=200 mm ---
    theta = np.linspace(0, 2 * np.pi, 300)
    r_inner = 200.0
    ax.plot(r_inner * np.cos(theta), r_inner * np.sin(theta),
            "r--", lw=1.2, label="Inner-ring r=200 mm")

    ax.set_xlabel("x (mm)", fontsize=11)
    ax.set_ylabel("y (mm)", fontsize=11)
    ax.set_title(
        f"PbWO4 module bimodality map — {run_label}\n"
        "(Red=strong bimodal, Yellow=mild, Green=unimodal, Grey=low-stat)",
        fontsize=11,
    )
    ax.grid(alpha=0.3)
    ax.set_aspect("equal")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(str(out_path), dpi=150)
    plt.close(fig)
    print(f"Wrote: {out_path}", flush=True)


# ---------------------------------------------------------------------------
# Output 3: Spatial summary PNG (2x2)
# ---------------------------------------------------------------------------

def write_spatial_png(results: list[dict], out_path: Path,
                       run_label: str, bimodal_ratio: float) -> None:
    bimodal = [r for r in results if r["is_bimodal"] == "Y"]
    unimodal = [r for r in results if r["is_bimodal"] == "N"]

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle(f"Bimodality spatial analysis — {run_label}", fontsize=13)

    # ---- Panel 1 (top-left): Distance histograms ----------------------------
    ax = axes[0, 0]
    bins_dist = np.arange(0, 420, 20)
    d_bi = [r["distance_from_beam_hole_mm"] for r in bimodal
            if not math.isnan(r["distance_from_beam_hole_mm"])]
    d_uni = [r["distance_from_beam_hole_mm"] for r in unimodal
             if not math.isnan(r["distance_from_beam_hole_mm"])]
    ax.hist(d_uni, bins=bins_dist, color="steelblue", alpha=0.6,
            label=f"Unimodal (N={len(d_uni)})")
    ax.hist(d_bi, bins=bins_dist, color="red", alpha=0.6, hatch="//",
            label=f"Bimodal (N={len(d_bi)})")
    ax.set_xlabel("Distance from beam hole (mm)", fontsize=10)
    ax.set_ylabel("Module count", fontsize=10)
    ax.set_title("Distance from beam hole", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # ---- Panel 2 (top-right): Distance vs ratio scatter ----------------------
    ax = axes[0, 1]
    two_peak = [r for r in results
                if not math.isnan(r.get("peak_2_over_peak_1_amp", _NAN))]
    d_vals = [r["distance_from_beam_hole_mm"] for r in two_peak]
    ratio_vals = [r["peak_2_over_peak_1_amp"] for r in two_peak]
    colors = ["red" if r["is_bimodal"] == "Y" else "steelblue" for r in two_peak]
    ax.scatter(d_vals, ratio_vals, c=colors, alpha=0.7, s=60)
    ax.axhline(bimodal_ratio, color="k", ls="--", lw=1.2,
               label=f"Bimodal threshold ({bimodal_ratio:.2f})")
    ax.set_xlabel("Distance from beam hole (mm)", fontsize=10)
    ax.set_ylabel("Peak 2 / Peak 1 amplitude", fontsize=10)
    ax.set_title("Peak ratio vs distance from beam hole", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # ---- Panel 3 (bottom-left): Row/col scatter of bimodal modules ----------
    ax = axes[1, 0]
    if bimodal:
        cols_bm = [r["seed_col"] for r in bimodal]
        rows_bm = [r["seed_row"] for r in bimodal]
        ratios_bm = [r["peak_2_over_peak_1_amp"]
                     if not math.isnan(r["peak_2_over_peak_1_amp"]) else 0.0
                     for r in bimodal]
        sc = ax.scatter(cols_bm, rows_bm, c=ratios_bm, cmap="Reds",
                        s=120, alpha=0.85, edgecolors="darkred", linewidths=0.5)
        plt.colorbar(sc, ax=ax, label="Peak 2 / Peak 1 amplitude")
    # Exclusion rectangle in row/col space (rows 15-20, cols 15-20)
    excl_box_x = [15, 20, 20, 15, 15]
    excl_box_y = [15, 15, 20, 20, 15]
    ax.plot(excl_box_x, excl_box_y, "k--", lw=1.5, label="Excl. rows/cols 15–20")
    ax.set_xlabel("seed_col", fontsize=10)
    ax.set_ylabel("seed_row", fontsize=10)
    ax.set_title("Spatial location of bimodal modules (row/col)", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    ax.invert_yaxis()

    # ---- Panel 4 (bottom-right): Fraction bimodal vs distance bin ----------
    ax = axes[1, 1]
    dist_bins = [0, 50, 100, 150, 200, 300, 400]
    bin_labels = ["0–50", "50–100", "100–150", "150–200", "200–300", "300–400"]
    fracs = []
    for i in range(len(dist_bins) - 1):
        lo_b, hi_b = dist_bins[i], dist_bins[i + 1]
        mods_in_bin = [r for r in results
                       if not math.isnan(r["distance_from_beam_hole_mm"])
                       and lo_b <= r["distance_from_beam_hole_mm"] < hi_b]
        if mods_in_bin:
            n_bm = sum(1 for r in mods_in_bin if r["is_bimodal"] == "Y")
            fracs.append(n_bm / len(mods_in_bin))
        else:
            fracs.append(0.0)
    x_pos = np.arange(len(bin_labels))
    bar_colors = ["steelblue" if f < bimodal_ratio else "salmon" for f in fracs]
    ax.bar(x_pos, fracs, color=bar_colors, edgecolor="black", linewidth=0.5)
    ax.axhline(bimodal_ratio, color="k", ls="--", lw=1.0,
               label=f"Bimodal threshold ({bimodal_ratio:.2f})")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(bin_labels, rotation=15, ha="right", fontsize=8)
    ax.set_xlabel("Distance from beam hole (mm)", fontsize=10)
    ax.set_ylabel("Fraction bimodal", fontsize=10)
    ax.set_title("Fraction bimodal vs distance from beam hole", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, axis="y")
    ax.set_ylim(0, 1.05)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(str(out_path), dpi=150)
    plt.close(fig)
    print(f"Wrote: {out_path}", flush=True)


# ---------------------------------------------------------------------------
# Output 4: Focused PDF of worst N offenders
# ---------------------------------------------------------------------------

def write_top_pdf(results: list[dict], module_data: dict,
                   out_path: Path, run_label: str,
                   top_n: int, e_lo: float, e_hi: float,
                   bin_width: float, target_peak: float) -> None:
    # Get bimodal modules sorted by ratio descending
    bimodal = [r for r in results if r["is_bimodal"] == "Y"]
    bimodal_sorted = sorted(
        bimodal,
        key=lambda r: r["peak_2_over_peak_1_amp"] if not math.isnan(r["peak_2_over_peak_1_amp"]) else 0.0,
        reverse=True,
    )
    top_modules = bimodal_sorted[:top_n]

    n_pages = len(top_modules)
    if n_pages == 0:
        # Write empty PDF with note
        with PdfPages(str(out_path)) as pdf:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.text(0.5, 0.5, "No bimodal modules detected.",
                    ha="center", va="center", fontsize=14, transform=ax.transAxes)
            ax.axis("off")
            fig.suptitle(f"Bimodal modules — {run_label}", fontsize=12)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
        print(f"Wrote: {out_path} (0 pages — no bimodal modules)", flush=True)
        return

    edges = np.arange(e_lo, e_hi + bin_width, bin_width)
    centres = 0.5 * (edges[:-1] + edges[1:])
    x_fit = np.linspace(e_lo, e_hi, 600)

    with PdfPages(str(out_path)) as pdf:
        for r in top_modules:
            mid = r["module_name"]
            # Get raw_vals from results dict (stored as _raw_vals)
            raw_vals = r.get("_raw_vals")
            if raw_vals is None:
                # Fallback: try module_data
                if mid in module_data:
                    raw_vals = module_data[mid]["raw_vals"]
                else:
                    continue

            counts, _ = np.histogram(raw_vals, bins=edges)
            ratio = r["peak_2_over_peak_1_amp"]
            ratio_str = f"{ratio:.2f}" if not math.isnan(ratio) else "N/A"

            fig, ax = plt.subplots(figsize=(10, 6))
            ax.step(centres, counts, where="mid", color="#555555", lw=1.4,
                    label=f"N = {r['n_events']}")

            # Overlay Gaussian peaks if fitted
            mu1 = r["peak_1_mu_MeV"]
            sig1 = r["peak_1_sigma_MeV"]
            amp1 = r["peak_1_amp"]
            mu2 = r["peak_2_mu_MeV"]
            sig2 = r["peak_2_sigma_MeV"]
            amp2 = r["peak_2_amp"]

            if not math.isnan(mu1) and not math.isnan(sig1) and not math.isnan(amp1):
                y1 = _gauss(x_fit, amp1, mu1, sig1)
                ax.plot(x_fit, y1, "b-", lw=1.8,
                        label=f"Peak 1: μ={mu1:.0f} MeV, σ={sig1:.0f} MeV")
            if not math.isnan(mu2) and not math.isnan(sig2) and not math.isnan(amp2):
                y2 = _gauss(x_fit, amp2, mu2, sig2)
                ax.plot(x_fit, y2, "r-", lw=1.8,
                        label=f"Peak 2: μ={mu2:.0f} MeV, σ={sig2:.0f} MeV")

            ax.axvline(target_peak, color="#444444", ls="--", lw=1.2,
                       label=f"Target {target_peak:.1f} MeV")

            # Text box
            mu1_str = f"{mu1:.0f}" if not math.isnan(mu1) else "N/A"
            mu2_str = f"{mu2:.0f}" if not math.isnan(mu2) else "N/A"
            info_txt = (
                f"Module: {mid}\n"
                f"Row={r['seed_row']}, Col={r['seed_col']}\n"
                f"N events: {r['n_events']}\n"
                f"Peak 1 μ: {mu1_str} MeV\n"
                f"Peak 2 μ: {mu2_str} MeV\n"
                f"Ratio: {ratio_str}"
            )
            props = dict(boxstyle="round", facecolor="lightyellow", alpha=0.85)
            ax.text(0.97, 0.97, info_txt, transform=ax.transAxes,
                    fontsize=8, va="top", ha="right", bbox=props)

            ax.set_xlabel("cluster_E_raw_MeV (MeV)", fontsize=10)
            ax.set_ylabel("Events / bin", fontsize=10)
            ax.set_xlim(0, 3500)
            ax.set_title(
                f"Module {mid} — bimodal spectrum, ratio={ratio_str}, "
                f"row={r['seed_row']}, col={r['seed_col']}",
                fontsize=11,
            )
            ax.legend(loc="upper left", fontsize=9)
            ax.grid(alpha=0.3)

            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"Wrote: {out_path} ({n_pages} pages)", flush=True)


# ---------------------------------------------------------------------------
# Output 5: Findings Markdown
# ---------------------------------------------------------------------------

def write_findings_md(results: list[dict], input_path: Path,
                       out_path: Path, run_label: str,
                       min_events: int) -> None:
    bimodal = [r for r in results if r["is_bimodal"] == "Y"]
    unimodal = [r for r in results if r["is_bimodal"] == "N"]

    n_total = len(results)
    n_bimodal = len(bimodal)
    pct_bi = 100.0 * n_bimodal / n_total if n_total else 0.0

    strong = [r for r in bimodal
              if not math.isnan(r["peak_2_over_peak_1_amp"])
              and r["peak_2_over_peak_1_amp"] >= 0.50]
    mild = [r for r in bimodal
            if not math.isnan(r["peak_2_over_peak_1_amp"])
            and 0.30 <= r["peak_2_over_peak_1_amp"] < 0.50]

    def _mean_med(lst: list[float]) -> tuple[str, str]:
        if not lst:
            return "N/A", "N/A"
        arr = np.array(lst)
        return f"{arr.mean():.1f}", f"{np.median(arr):.1f}"

    d_bi = [r["distance_from_beam_hole_mm"] for r in bimodal
            if not math.isnan(r["distance_from_beam_hole_mm"])]
    d_uni = [r["distance_from_beam_hole_mm"] for r in unimodal
             if not math.isnan(r["distance_from_beam_hole_mm"])]
    bi_mean, bi_med = _mean_med(d_bi)
    uni_mean, uni_med = _mean_med(d_uni)

    # Top-10 table (sorted by ratio descending)
    top10 = sorted(bimodal,
                   key=lambda r: r["peak_2_over_peak_1_amp"]
                   if not math.isnan(r["peak_2_over_peak_1_amp"]) else 0.0,
                   reverse=True)[:10]

    # Spatial interpretation
    if d_bi:
        mean_d = float(np.mean(d_bi))
        if mean_d < 150:
            spatial_interp = (
                "Bimodal modules are concentrated at small distances from the beam hole, "
                "suggesting proximity to the high-rate central region may be a factor."
            )
        elif mean_d > 250:
            spatial_interp = (
                "Bimodal modules tend to cluster at larger radii from the beam hole, "
                "suggesting an outer-ring or edge-region phenomenon."
            )
        else:
            spatial_interp = (
                "Bimodal modules are distributed across a range of distances from "
                "the beam hole with no strong radial concentration. "
                "The effect appears spatially scattered rather than clustered."
            )
    else:
        spatial_interp = "No bimodal modules detected; spatial interpretation not applicable."

    today = date.today().isoformat()

    lines = [
        f"# Bimodality Findings — run {run_label}",
        "",
        f"**Date generated:** {today}",
        f"**Input TSV:** `{input_path}`",
        "",
        "## Summary",
        "",
        f"| Item | Value |",
        f"|---|---|",
        f"| Total PbWO4 modules in TSV | (all with ≥1 event) |",
        f"| Modules with ≥{min_events} events (analyzed) | {n_total} |",
        f"| Bimodal modules | {n_bimodal} ({pct_bi:.1f}%) |",
        f"| Strongly bimodal (ratio ≥ 0.50) | {len(strong)} |",
        f"| Mildly bimodal (0.30 ≤ ratio < 0.50) | {len(mild)} |",
        "",
        "## Distance from Beam Hole",
        "",
        f"| Group | Mean (mm) | Median (mm) |",
        f"|---|---|---|",
        f"| Bimodal | {bi_mean} | {bi_med} |",
        f"| Unimodal | {uni_mean} | {uni_med} |",
        "",
        "## Top 10 Most Severe Bimodal Modules",
        "",
        "| Module | Row | Col | N events | d_BH (mm) | Peak 1 μ (MeV) "
        "| Peak 2 μ (MeV) | Ratio |",
        "|---|---|---|---|---|---|---|---|",
    ]

    for r in top10:
        mu1_s = f"{r['peak_1_mu_MeV']:.0f}" if not math.isnan(r["peak_1_mu_MeV"]) else "N/A"
        mu2_s = f"{r['peak_2_mu_MeV']:.0f}" if not math.isnan(r["peak_2_mu_MeV"]) else "N/A"
        rat_s = f"{r['peak_2_over_peak_1_amp']:.3f}" if not math.isnan(r["peak_2_over_peak_1_amp"]) else "N/A"
        d_s = f"{r['distance_from_beam_hole_mm']:.1f}" if not math.isnan(r["distance_from_beam_hole_mm"]) else "N/A"
        lines.append(
            f"| {r['module_name']} | {r['seed_row']} | {r['seed_col']} "
            f"| {r['n_events']} | {d_s} | {mu1_s} | {mu2_s} | {rat_s} |"
        )

    lines += [
        "",
        "## Spatial Interpretation",
        "",
        spatial_interp,
        "",
        "---",
        f"*Generated by `scan_module_bimodality.py` on {today}*",
    ]

    with open(str(out_path), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"Wrote: {out_path}", flush=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = build_parser()
    args = ap.parse_args(argv)

    if not HAVE_SCIPY:
        sys.stderr.write(
            "[ERROR] scipy is required. Install with: pip install scipy\n"
        )
        return 1
    if not HAVE_PANDAS:
        sys.stderr.write(
            "[ERROR] pandas is required. Install with: pip install pandas\n"
        )
        return 1

    e_lo, e_hi = args.e_range
    bin_width = args.bin_width
    target_peak = args.target_peak
    run_label = args.run_label
    min_events = args.min_events
    bimodal_ratio = args.bimodal_ratio
    prominence_frac = args.prominence_frac
    top_n = args.top_n_pdf

    # Resolve input
    input_path = _resolve_input(args.input_tsv)
    if not input_path.exists():
        sys.stderr.write(
            f"[ERROR] Input TSV not found: {input_path}\n"
            f"        Generate it with benchmark_hycal_timing.py --write-cluster-tsv\n"
        )
        return 1

    # Resolve output dir
    output_dir = _resolve_output_dir(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Output paths
    tsv_out = output_dir / f"module_bimodality_{run_label}.tsv"
    map_out = output_dir / f"module_bimodality_{run_label}_map.png"
    spatial_out = output_dir / f"module_bimodality_{run_label}_spatial.png"
    pdf_out = output_dir / f"module_bimodality_{run_label}_top{top_n}.pdf"
    md_out = output_dir / f"module_bimodality_{run_label}_findings.md"

    # Read data
    df_full, module_data, skipped_positions = read_tsv(input_path, min_events)

    # Run analysis
    results = run_analysis(module_data, e_lo, e_hi, bin_width,
                           prominence_frac, bimodal_ratio)

    # Output 1: Summary TSV
    write_summary_tsv(results, tsv_out)

    # Output 2: Detector face map
    write_map_png(results, skipped_positions, map_out, run_label, bimodal_ratio)

    # Output 3: Spatial summary
    write_spatial_png(results, spatial_out, run_label, bimodal_ratio)

    # Output 4: Focused PDF
    write_top_pdf(results, module_data, pdf_out, run_label,
                  top_n, e_lo, e_hi, bin_width, target_peak)

    # Output 5: Findings markdown
    write_findings_md(results, input_path, md_out, run_label, min_events)

    return 0


if __name__ == "__main__":
    sys.exit(main())
