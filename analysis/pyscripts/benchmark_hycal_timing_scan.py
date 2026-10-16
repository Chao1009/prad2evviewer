#!/usr/bin/env python3
"""benchmark_hycal_timing_scan.py — Multi-run EVIO scan for secondary peaks.

Scans all prad_*.evio.* files in an EVIO directory, invokes
benchmark_hycal_timing.py (legacy path + --calibrate) on each split file,
detects any secondary peaks in the energy spectrum (prominence threshold +
min distance from primary), and aggregates σ_E/E trends across files / runs.

Secondary-peak detection flags any peak detected by find_peaks() that is NOT
the primary and is more than min_dist_mev away from it, regardless of whether
it appears above or below the primary peak centre.

By default the full EVIO file is processed (no --max-events cap). Pass an
explicit --max-events value to impose a limit.

Outputs
-------
  output/scan_${RUN}/split_${SPLITNUM}/         per-file benchmark outputs
  output/benchmark_hycal_timing_scan_summary.tsv  aggregate (one row per file)
  output/benchmark_hycal_timing_scan_sigma_E_vs_file.png  trend plot
  output/benchmark_hycal_timing_scan_all_calibrated_plots.pdf  combined PDF
  output/benchmark_hycal_timing_scan_findings.md  human-readable summary
  benchmark_hycal_timing_scan.log                 persistent log

Usage
-----
  python analysis/pyscripts/benchmark_hycal_timing_scan.py [options]
  python analysis/pyscripts/benchmark_hycal_timing_scan.py --skip-benchmark
  python analysis/pyscripts/benchmark_hycal_timing_scan.py --parallel 2

Reference: user specification §1–10 (June 2025 inline request).
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime
import json
import logging
import math
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

try:
    from scipy.signal import find_peaks
    from scipy.optimize import curve_fit
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

# ---------------------------------------------------------------------------
# Repository / path constants
# ---------------------------------------------------------------------------
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_BENCHMARK_SCRIPT = str(
    _REPO_ROOT / "analysis" / "pyscripts" / "benchmark_hycal_timing.py"
)

_DEFAULT_EVIO_DIR  = "/Users/L00338853/work/PRad/data/evio"
_DEFAULT_OUTPUT_DIR = str(_REPO_ROOT / "analysis" / "pyscripts" / "output")
_DEFAULT_MAX_EVENTS: Optional[int] = None   # None → process full EVIO file
_DEFAULT_PROMINENCE = 40
_DEFAULT_PEAK_DIST_MEV = 100.0
_DEFAULT_PARALLEL   = 1

# ---------------------------------------------------------------------------
# Run-number → target peak (MeV) mapping  (spec §2)
# ---------------------------------------------------------------------------
_TARGET_PEAK_MAP: dict[int, float] = {
    24386: 3488.5,
    24503: 3488.5,
    24512: 3488.5,
    24590:  728.9,
    25308: 2239.5,
    25320: 2239.5,
    26138: 2239.5,
}
_DEFAULT_TARGET_PEAK = 3488.5  # fallback for unknown runs

# ---------------------------------------------------------------------------
# Run-number → plot colour  (spec §6)
# ---------------------------------------------------------------------------
_RUN_COLORS: dict[int, str] = {
    24386: "blue",
    24503: "green",
    24512: "orange",
    24590: "red",
    25308: "purple",
    25320: "brown",
    26138: "pink",
}
_DEFAULT_COLOR = "gray"

# ---------------------------------------------------------------------------
# Summary column schema  (spec §5)
# ---------------------------------------------------------------------------
_SUMMARY_COLS = [
    "run_number", "split_number", "filename",
    "n_events", "n_modules_used", "sigma_over_E_CB",
    "primary_peak_mu_MeV", "primary_peak_sigma_MeV",
    "n_peaks_detected", "secondary_peak_present",
    "secondary_peak_mu_MeV", "secondary_peak_sigma_MeV",
    "peak_prominence_values",
]

# Runinfo pre-window [lo, hi] ns used for each run family
_PRE_WINDOW_MAP: dict[int, tuple[float, float]] = {
    24386: (130.0, 200.0),
    24503: (130.0, 200.0),
    24512: (130.0, 200.0),
    24590: (130.0, 200.0),
    25308: (100.0, 200.0),
    25320: (100.0, 200.0),
    26138: (100.0, 200.0),
}
_DEFAULT_PRE_WINDOW = (100.0, 200.0)

# signal-min cut — we use a broad cut so the headline captures secondary peaks
# below the elastic line; 500 MeV is compatible with all runs
_SIGNAL_MIN = 500.0

# Energy histogram range for find_peaks input
_E_RANGE_MAP: dict[int, tuple[float, float]] = {
    24386: (500.0, 4500.0),
    24503: (500.0, 4500.0),
    24512: (500.0, 4500.0),
    24590: (500.0, 2000.0),
    25308: (500.0, 3500.0),
    25320: (500.0, 3500.0),
    26138: (500.0, 3500.0),
}
_DEFAULT_E_RANGE = (500.0, 4500.0)

# ---------------------------------------------------------------------------
# Logging setup — module-level; callers add file handler after args are parsed
# ---------------------------------------------------------------------------
log = logging.getLogger("scan")


# ===========================================================================
# §1 — File discovery
# ===========================================================================

_EVIO_PAT = re.compile(r"^prad_0*(\d+)\.evio\.(\d+)$", re.IGNORECASE)


def discover_files(evio_dir: Path) -> dict[int, list[tuple[int, Path]]]:
    """Return {run_number: [(split_int, Path), ...]} sorted by split."""
    runs: dict[int, list[tuple[int, Path]]] = {}
    for p in sorted(evio_dir.iterdir()):
        if not p.is_file():
            continue
        m = _EVIO_PAT.match(p.name)
        if not m:
            continue
        run   = int(m.group(1))
        split = int(m.group(2))
        runs.setdefault(run, []).append((split, p))
    for run in runs:
        runs[run].sort(key=lambda t: t[0])
    return runs


# ===========================================================================
# §2 — Per-file benchmark invocation
# ===========================================================================

def benchmark_out_dir(output_dir: Path, run: int, split: int) -> Path:
    return output_dir / f"scan_{run:06d}" / f"split_{split:05d}"


def benchmark_tsv_path(output_dir: Path, run: int, split: int) -> Path:
    d = benchmark_out_dir(output_dir, run, split)
    return d / "result.tsv"


def run_benchmark(
    *,
    evio_file: Path,
    run: int,
    split: int,
    output_dir: Path,
    max_events: Optional[int],
    skip_benchmark: bool,
) -> Path:
    """Invoke benchmark_hycal_timing.py; return path to the TSV.

    If skip_benchmark is True and the TSV already exists, skip the invocation.
    If max_events is None, --max-events is NOT passed to the subprocess so the
    full EVIO file is processed.
    Raises RuntimeError on subprocess failure.
    """
    out_d   = benchmark_out_dir(output_dir, run, split)
    tsv_p   = benchmark_tsv_path(output_dir, run, split)
    log_p   = out_d / "benchmark.log"
    out_stem = str(out_d / "result")

    if skip_benchmark and tsv_p.exists():
        log.info("skip (TSV exists): %s", tsv_p)
        return tsv_p

    out_d.mkdir(parents=True, exist_ok=True)

    target_peak = _TARGET_PEAK_MAP.get(run, _DEFAULT_TARGET_PEAK)
    pre_lo, pre_hi = _PRE_WINDOW_MAP.get(run, _DEFAULT_PRE_WINDOW)
    e_lo, e_hi = _E_RANGE_MAP.get(run, _DEFAULT_E_RANGE)

    cmd = [
        sys.executable,
        _BENCHMARK_SCRIPT,
        str(evio_file),
        out_stem,
        "--calibrate",
        "--target-peak", str(target_peak),
        "--exclude-rowcol", "15", "20", "15", "20",
        "--signal-min", str(_SIGNAL_MIN),
        "--pre-window", str(pre_lo), str(pre_hi),
        "--e-range", str(e_lo), str(e_hi),
        "--bin-width", "20",
        "--inner-radius", "200",
        "--min-seed-energy", "200",
        "--min-per-module", "30",
        "--cal-method", "cb-peak",
    ]
    # Only pass --max-events when an explicit limit is requested.
    # When max_events is None the subprocess processes the full EVIO file.
    if max_events is not None:
        cmd.extend(["--max-events", str(max_events)])

    max_ev_label = str(max_events) if max_events is not None else "full file"
    log.info("run %06d split %05d → %s  (max_events=%s)",
             run, split, out_stem, max_ev_label)
    print(f"Processing run {run} split {split:05d}  (max_events={max_ev_label})...",
          flush=True)

    t0 = time.monotonic()
    try:
        with open(log_p, "w") as lf:
            result = subprocess.run(
                cmd,
                stdout=lf,
                stderr=subprocess.STDOUT,
                check=False,
            )
        elapsed = time.monotonic() - t0
        if result.returncode != 0:
            msg = (
                f"benchmark failed (rc={result.returncode}) "
                f"run={run} split={split}. Log: {log_p}"
            )
            log.error(msg)
            raise RuntimeError(msg)
        log.info("  done in %.1f s → %s", elapsed, tsv_p)
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"benchmark script not found: {_BENCHMARK_SCRIPT}"
        ) from exc

    if not tsv_p.exists():
        raise RuntimeError(
            f"Expected TSV not created: {tsv_p} (log: {log_p})"
        )
    return tsv_p


# ===========================================================================
# §2b — TSV parsing (reusing sweep_pre_window_24386 conventions)
# ===========================================================================

def _parse_float(s: str) -> float:
    s = s.strip()
    if not s:
        return float("nan")
    try:
        return float(s)
    except ValueError:
        return float("nan")


def parse_tsv(tsv_path: Path) -> dict:
    """Parse benchmark_hycal_timing TSV.

    Returns dict with keys:
      'headline'     : {'legacy': {...}, 'new': {...}}  (bulk path metrics)
      'calibration'  : {'legacy': {...}, 'new': {...}}  (inner-ring cal metrics)
      'histogram'    : {'centres': [...], 'legacy': [...], 'new': [...]}

    Raises RuntimeError loudly if headline block is missing (spec §10).
    calibration block absence only raises if --calibrate was expected.
    """
    text  = tsv_path.read_text(encoding="utf-8")
    lines = text.splitlines()

    headline_hdr_lineno    = None
    calibration_hdr_lineno = None
    hist_hdr_lineno        = None

    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith("path\tn_events\t") and headline_hdr_lineno is None:
            headline_hdr_lineno = i
        if s.startswith("# --calibrate inner ring"):
            calibration_hdr_lineno = i
        if s.startswith("e_centre\tcount_legacy\t") and hist_hdr_lineno is None:
            hist_hdr_lineno = i

    if headline_hdr_lineno is None:
        raise RuntimeError(
            f"Headline block missing in {tsv_path} "
            "(expected line starting 'path\\tn_events\\t...')"
        )

    # Parse headline block
    headline_cols = lines[headline_hdr_lineno].strip().split("\t")
    headline: dict[str, dict] = {}
    for line in lines[headline_hdr_lineno + 1:]:
        s = line.strip()
        if not s or s.startswith("#"):
            break
        parts = s.split("\t")
        if len(parts) < len(headline_cols):
            continue
        row  = dict(zip(headline_cols, parts))
        name = row.pop("path")
        headline[name] = {k: _parse_float(v) for k, v in row.items()}

    if not headline:
        raise RuntimeError(
            f"Headline block present in {tsv_path} but contained no data rows."
        )

    # Parse histogram block
    hist: dict = {}
    if hist_hdr_lineno is not None:
        centres_list, legacy_list, new_list = [], [], []
        for line in lines[hist_hdr_lineno + 1:]:
            s = line.strip()
            if not s or s.startswith("#"):
                break
            parts = s.split("\t")
            if len(parts) < 3:
                continue
            try:
                centres_list.append(float(parts[0]))
                legacy_list.append(int(parts[1]))
                new_list.append(int(parts[2]))
            except ValueError:
                continue
        hist = {
            "centres": np.array(centres_list, dtype=float),
            "legacy":  np.array(legacy_list,  dtype=int),
            "new":     np.array(new_list,      dtype=int),
        }

    # Parse calibration block (optional — only if --calibrate was used)
    calibration: dict[str, dict] = {}
    if calibration_hdr_lineno is not None:
        cal_data_hdr = None
        for i in range(calibration_hdr_lineno + 1, len(lines)):
            s = lines[i].strip()
            if s.startswith("path\tn_modules\t"):
                cal_data_hdr = i
                break
        if cal_data_hdr is not None:
            cal_cols = lines[cal_data_hdr].strip().split("\t")
            for line in lines[cal_data_hdr + 1:]:
                s = line.strip()
                if not s or s.startswith("#"):
                    break
                parts = s.split("\t")
                if len(parts) < len(cal_cols):
                    continue
                row  = dict(zip(cal_cols, parts))
                name = row.pop("path")
                calibration[name] = {k: _parse_float(v) for k, v in row.items()}

    return {
        "headline":    headline,
        "calibration": calibration,
        "histogram":   hist,
    }


# ===========================================================================
# §3 — Peak finding & Gaussian fitting
# ===========================================================================

def _gauss(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def fit_gaussian_peak(centres: np.ndarray, counts: np.ndarray,
                      peak_idx: int, window_mev: float = 150.0
                      ) -> tuple[float, float, float]:
    """Fit a Gaussian around centres[peak_idx] ± window_mev.

    Returns (mu, sigma, amplitude); NaN on failure.
    """
    nan = float("nan")
    if not HAVE_SCIPY or len(centres) == 0:
        return nan, nan, nan
    peak_e = centres[peak_idx]
    mask   = (centres >= peak_e - window_mev) & (centres <= peak_e + window_mev)
    x, y   = centres[mask], counts[mask].astype(float)
    if y.sum() < 20:
        return nan, nan, nan
    try:
        popt, _ = curve_fit(
            _gauss, x, y,
            p0     = [float(y.max()), peak_e, window_mev / 4.0],
            bounds = ([0, peak_e - window_mev, 1.0],
                      [np.inf, peak_e + window_mev, window_mev]),
            maxfev = 5000,
        )
        return float(popt[1]), abs(float(popt[2])), float(popt[0])
    except Exception:
        return nan, nan, nan


def detect_peaks(
    histogram: dict,
    target_peak: float,
    prominence_threshold: float,
    min_dist_mev: float,
) -> dict:
    """Detect primary and secondary peaks in the legacy energy spectrum.

    Secondary peaks are flagged for any peak distinct from the primary
    (>min_dist_mev away), regardless of whether the peak appears above or
    below the primary peak centre.

    Returns dict with keys:
      n_peaks            : int
      primary_mu         : float  (MeV)
      primary_sigma      : float  (MeV)
      secondary_present  : bool   (True for any secondary peak found)
      secondary_mu       : float  (MeV, NaN if absent)
      secondary_sigma    : float  (MeV, NaN if absent)
      prominences        : list[float]
    """
    nan = float("nan")
    result = {
        "n_peaks":           0,
        "primary_mu":        nan,
        "primary_sigma":     nan,
        "secondary_present": False,
        "secondary_mu":      nan,
        "secondary_sigma":   nan,
        "prominences":       [],
    }

    if not HAVE_SCIPY or not histogram or "centres" not in histogram:
        return result
    if histogram["centres"].size == 0:
        return result

    centres = histogram["centres"]
    counts  = histogram["legacy"].astype(float)

    # Bin width (MeV) → convert distance threshold to bins
    if len(centres) < 2:
        return result
    bin_w    = float(centres[1] - centres[0])
    dist_bins = max(1, int(min_dist_mev / bin_w))

    peaks, props = find_peaks(
        counts,
        prominence=prominence_threshold,
        distance=dist_bins,
    )
    if len(peaks) == 0:
        return result

    prominences = list(props["prominences"])
    result["prominences"] = [float(p) for p in prominences]
    result["n_peaks"] = len(peaks)

    # Identify the primary peak as the one closest to target_peak
    dists_to_target = np.abs(centres[peaks] - target_peak)
    primary_idx_in_peaks = int(np.argmin(dists_to_target))
    primary_bin   = peaks[primary_idx_in_peaks]
    primary_mu, primary_sigma, _ = fit_gaussian_peak(centres, counts, primary_bin)
    result["primary_mu"]    = primary_mu
    result["primary_sigma"] = primary_sigma

    # Secondary peak: any peak that is > min_dist_mev away from the primary peak
    # centre. Position relative to the primary is not considered.
    primary_e = centres[primary_bin]
    # Use the Gaussian-fitted centre when available; fall back to histogram bin centre.
    primary_peak_center = primary_mu if not math.isnan(primary_mu) else primary_e
    secondary_bins = [
        pk for pk in peaks
        if abs(centres[pk] - primary_e) > min_dist_mev
    ]

    if secondary_bins:
        result["secondary_present"] = True
        # Pick the secondary with the highest prominence
        secondary_proms = [
            props["prominences"][np.where(peaks == sb)[0][0]]
            for sb in secondary_bins
        ]
        best_sec_idx = int(np.argmax(secondary_proms))
        sec_bin = secondary_bins[best_sec_idx]
        sec_mu, sec_sigma, _ = fit_gaussian_peak(centres, counts, sec_bin)
        result["secondary_mu"]    = sec_mu
        result["secondary_sigma"] = sec_sigma

        log.info(
            "  Found secondary peak at %.0f MeV "
            "(σ=%.0f MeV, prominence=%.1f, primary=%.0f MeV)",
            sec_mu if not math.isnan(sec_mu) else centres[sec_bin],
            sec_sigma if not math.isnan(sec_sigma) else float("nan"),
            secondary_proms[best_sec_idx],
            primary_peak_center,
        )
        print(
            f"  Found secondary peak at "
            f"{sec_mu if not math.isnan(sec_mu) else centres[sec_bin]:.0f} MeV"
            f" (primary at {primary_peak_center:.0f} MeV)",
            flush=True,
        )

    return result


# ===========================================================================
# §4 — Aggregate per-file metrics
# ===========================================================================

def extract_row(
    *,
    run: int,
    split: int,
    evio_file: Path,
    tsv_path: Path,
    prominence_threshold: float,
    min_dist_mev: float,
) -> dict:
    """Parse TSV and peak-detect; return one aggregate row dict."""
    nan = float("nan")
    row: dict = {
        "run_number":                  run,
        "split_number":                split,
        "filename":                    evio_file.name,
        "n_events":                    nan,
        "n_modules_used":              nan,
        "sigma_over_E_CB":             nan,
        "primary_peak_mu_MeV":         nan,
        "primary_peak_sigma_MeV":      nan,
        "n_peaks_detected":        0,
        "secondary_peak_present":  "N",
        "secondary_peak_mu_MeV":   nan,
        "secondary_peak_sigma_MeV":    nan,
        "peak_prominence_values":      "[]",
    }

    parsed = parse_tsv(tsv_path)
    hl  = parsed["headline"].get("legacy", {})
    cal = parsed["calibration"].get("legacy", {})
    hist = parsed["histogram"]

    row["n_events"]       = hl.get("n_events", nan)
    row["n_modules_used"] = cal.get("n_modules", nan)
    row["sigma_over_E_CB"] = cal.get("sigma_over_E_cb", nan)

    # Peak detection on full headline spectrum
    target_peak = _TARGET_PEAK_MAP.get(run, _DEFAULT_TARGET_PEAK)
    peaks_info  = detect_peaks(hist, target_peak, prominence_threshold, min_dist_mev)

    row["n_peaks_detected"]       = peaks_info["n_peaks"]
    row["primary_peak_mu_MeV"]    = peaks_info["primary_mu"]
    row["primary_peak_sigma_MeV"] = peaks_info["primary_sigma"]
    row["secondary_peak_present"] = "Y" if peaks_info["secondary_present"] else "N"
    row["secondary_peak_mu_MeV"]  = peaks_info["secondary_mu"]
    row["secondary_peak_sigma_MeV"] = peaks_info["secondary_sigma"]
    row["peak_prominence_values"]    = json.dumps(
        [round(p, 2) for p in peaks_info["prominences"]]
    )

    return row


# ===========================================================================
# §5 — Write aggregate TSV
# ===========================================================================

def _fmt(val, col: str) -> str:
    if isinstance(val, float):
        if math.isnan(val):
            return "nan"
        if col in ("run_number", "split_number", "n_events",
                   "n_modules_used", "n_peaks_detected"):
            return f"{val:.0f}"
        if col in ("sigma_over_E_CB",):
            return f"{val:.4f}"
        return f"{val:.2f}"
    return str(val)


# ---------------------------------------------------------------------------
# §5b — Calibrated PNG path helper
# ---------------------------------------------------------------------------

def calibrated_png_path(output_dir: Path, run: int, split: int) -> Path:
    """Return the expected path for a per-file calibrated result PNG."""
    return (
        output_dir / f"scan_{run:06d}" / f"split_{split:05d}"
        / "result_calibrated.png"
    )


def write_summary_tsv(rows: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(_SUMMARY_COLS) + "\n")
        for row in rows:
            fh.write(
                "\t".join(_fmt(row.get(c, float("nan")), c) for c in _SUMMARY_COLS)
                + "\n"
            )


# ===========================================================================
# §6 — Trend plot
# ===========================================================================

def _label(run: int, split: int) -> str:
    return f"run_{run}.split_{split:05d}"


def make_trend_plot(rows: list[dict], out_png: Path) -> None:
    """σ_E/E (CB) vs. EVIO file, colour by run, star for secondary peaks."""
    out_png.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(max(12, len(rows) * 0.18 + 4), 6))

    # Group rows by run for run-level legends
    run_groups: dict[int, list[dict]] = {}
    for row in rows:
        run_groups.setdefault(int(row["run_number"]), []).append(row)

    # Draw per-run lines + per-point markers
    all_x_labels: list[str] = []
    all_x_positions: list[int] = []
    x_pos = 0

    x_by_run: dict[int, list[int]]   = {}
    y_by_run: dict[int, list[float]] = {}
    sec_x:    list[int]   = []
    sec_y:    list[float] = []

    for run in sorted(run_groups.keys()):
        color = _RUN_COLORS.get(run, _DEFAULT_COLOR)
        run_rows = run_groups[run]
        run_x: list[int]   = []
        run_y: list[float] = []

        for row in run_rows:
            lbl = _label(int(row["run_number"]), int(row["split_number"]))
            all_x_labels.append(lbl)
            all_x_positions.append(x_pos)

            sigma_e = row.get("sigma_over_E_CB", float("nan"))
            if isinstance(sigma_e, float) and not math.isnan(sigma_e):
                sigma_e_pct = sigma_e * 100.0
                run_x.append(x_pos)
                run_y.append(sigma_e_pct)
                if row.get("secondary_peak_present") == "Y":
                    sec_x.append(x_pos)
                    sec_y.append(sigma_e_pct)

            x_pos += 1

        x_by_run[run] = run_x
        y_by_run[run] = run_y

        if run_x:
            ax.plot(run_x, run_y, "-o", color=color, lw=1.4, ms=5,
                    label=f"run {run}", zorder=3)

    # Overlay star markers for files with secondary peaks
    if sec_x:
        ax.scatter(sec_x, sec_y, marker="*", s=180, color="crimson",
                   zorder=5, label="secondary peak detected")

    ax.set_xticks(all_x_positions)
    ax.set_xticklabels(all_x_labels, rotation=45, ha="right", fontsize=7)
    ax.set_xlabel("EVIO file")
    ax.set_ylabel("σ_E/E CB (%)")
    ax.set_title(
        "σ_E/E (Crystal Ball) across EVIO files — all runs\n"
        "(stars = secondary peak detected)"
    )
    ax.legend(loc="upper right", fontsize=8, ncol=2)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"\nWrote trend plot: {out_png}", flush=True)


# ===========================================================================
# §7 — Summary markdown
# ===========================================================================

def write_findings_md(rows: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Group by run
    run_groups: dict[int, list[dict]] = {}
    for row in rows:
        run_groups.setdefault(int(row["run_number"]), []).append(row)

    lines: list[str] = []
    lines.append("# benchmark_hycal_timing_scan — Findings\n")
    lines.append(
        f"_Generated {datetime.datetime.now().isoformat(timespec='seconds')}_\n"
    )
    lines.append("")
    lines.append(
        "> **Note**: secondary peaks are flagged when a distinct peak is detected "
        "(prominence > threshold) at least `min_dist_mev` away from the primary "
        "(elastic) peak, regardless of position.\n"
    )
    lines.append("")
    lines.append("## Per-Run Summary\n")

    for run in sorted(run_groups.keys()):
        rrows = run_groups[run]
        sigmas = [
            r["sigma_over_E_CB"] * 100.0
            for r in rrows
            if isinstance(r["sigma_over_E_CB"], float)
            and not math.isnan(r["sigma_over_E_CB"])
        ]
        sec_rows = [r for r in rrows if r.get("secondary_peak_present") == "Y"]
        n_files  = len(rrows)
        mean_s   = float(np.mean(sigmas)) if sigmas else float("nan")
        std_s    = float(np.std(sigmas, ddof=1)) if len(sigmas) > 1 else float("nan")
        beam_e   = _TARGET_PEAK_MAP.get(run, _DEFAULT_TARGET_PEAK)

        lines.append(f"### Run {run}  (beam energy: {beam_e:.1f} MeV)\n")
        lines.append(f"- Files scanned: **{n_files}**")
        lines.append(f"- σ_E/E (CB): mean = **{mean_s:.2f}%**, stdev = {std_s:.2f}%"
                     if not math.isnan(mean_s) else "- σ_E/E (CB): insufficient data")
        lines.append(
            f"- Secondary peaks detected: **{len(sec_rows)}** / {n_files} files"
        )

        if sec_rows:
            sec_splits = [str(int(r["split_number"])) for r in sec_rows]
            lines.append(
                f"- Splits with secondary peaks: {', '.join(sec_splits)}"
            )

        # Anomaly: >10% relative variation in σ_E/E within run
        if len(sigmas) > 1 and not math.isnan(mean_s) and mean_s > 0:
            variation_pct = std_s / mean_s * 100.0
            if variation_pct > 10.0:
                lines.append(
                    f"- ⚠️  **Anomaly**: σ_E/E varies by {variation_pct:.1f}% "
                    f"(>10%) within this run."
                )

        lines.append("")

    # Overall summary table
    lines.append("## Summary Table: Secondary Peaks by Run\n")
    lines.append(
        "| Run | Beam E (MeV) | Files | Secondary peaks | Mean σ_E/E CB (%) | Assessment |"
    )
    lines.append(
        "|-----|-------------|-------|-----------------|-------------------|------------|"
    )

    for run in sorted(run_groups.keys()):
        rrows    = run_groups[run]
        sec_rows = [r for r in rrows if r.get("secondary_peak_present") == "Y"]
        sigmas   = [
            r["sigma_over_E_CB"] * 100.0
            for r in rrows
            if isinstance(r["sigma_over_E_CB"], float)
            and not math.isnan(r["sigma_over_E_CB"])
        ]
        mean_s  = float(np.mean(sigmas)) if sigmas else float("nan")
        has_sec = len(sec_rows) > 0
        beam_e  = _TARGET_PEAK_MAP.get(run, _DEFAULT_TARGET_PEAK)

        if has_sec:
            assessment = "⚠️ secondary peak(s) detected — review PDF for interpretation"
        else:
            assessment = "✅ no secondary peaks detected"

        mean_str = f"{mean_s:.2f}" if not math.isnan(mean_s) else "N/A"
        lines.append(
            f"| {run} | {beam_e:.1f} | {len(rrows)} | {len(sec_rows)} "
            f"| {mean_str} | {assessment} |"
        )

    lines.append("")

    # Data quality assessment
    lines.append("## Data Quality Assessment\n")
    lines.append("### Recommended for precision physics\n")
    clean_runs = [
        run for run in sorted(run_groups.keys())
        if not any(r.get("secondary_peak_present") == "Y"
                   for r in run_groups[run])
    ]
    if clean_runs:
        lines.append(
            "Runs with no secondary peaks detected in any split:"
        )
        for run in clean_runs:
            lines.append(f"- Run {run}")
    else:
        lines.append("_(No runs are fully free of secondary peaks.)_")

    lines.append("")
    lines.append("### Flagged for contamination / double-peak issues\n")
    contaminated = [
        run for run in sorted(run_groups.keys())
        if any(r.get("secondary_peak_present") == "Y"
               for r in run_groups[run])
    ]
    if contaminated:
        for run in contaminated:
            sec_splits = [
                str(int(r["split_number"]))
                for r in run_groups[run]
                if r.get("secondary_peak_present") == "Y"
            ]
            lines.append(
                f"- **Run {run}**: secondary peaks in splits {', '.join(sec_splits)}"
            )
    else:
        lines.append("_(No contamination detected.)_")

    lines.append("")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nWrote findings: {out_path}", flush=True)


# ===========================================================================
# §7b — Combined calibrated-plot PDF
# ===========================================================================

def generate_combined_pdf(
    rows: list[dict],
    output_dir: Path,
    pdf_path: Path,
) -> None:
    """Generate a multi-page PDF of all per-file calibrated result PNGs.

    Each page contains one calibrated plot from::

        output/scan_RRRRRR/split_SSSSS/result_calibrated.png

    Pages are ordered by run number then split number.  If a PNG is missing
    (the benchmark failed for that file) a placeholder page is emitted instead.

    The page title encodes run/split, σ_E/E, and whether a secondary peak was
    detected (Y/N) and at what energy.

    Parameters
    ----------
    rows:
        Aggregate rows from extract_row(), already sorted by (run, split).
    output_dir:
        Base output directory (parent of scan_RRRRRR subdirectories).
    pdf_path:
        Destination path for the combined PDF.
    """
    pdf_path.parent.mkdir(parents=True, exist_ok=True)

    # Sort defensively (rows should already be sorted, but be safe)
    sorted_rows = sorted(
        rows,
        key=lambda r: (int(r["run_number"]), int(r["split_number"])),
    )

    with PdfPages(pdf_path) as pdf:
        for row in sorted_rows:
            run   = int(row["run_number"])
            split = int(row["split_number"])
            png_p = calibrated_png_path(output_dir, run, split)

            sigma_e = row.get("sigma_over_E_CB", float("nan"))
            try:
                sigma_e_f = float(sigma_e)
            except (TypeError, ValueError):
                sigma_e_f = float("nan")

            sigma_str = (
                f"{sigma_e_f * 100.0:.2f}%" if not math.isnan(sigma_e_f) else "N/A"
            )

            sec_present = row.get("secondary_peak_present", "N")
            sec_mu = row.get("secondary_peak_mu_MeV", float("nan"))
            try:
                sec_mu_f = float(sec_mu)
            except (TypeError, ValueError):
                sec_mu_f = float("nan")

            if sec_present == "Y" and not math.isnan(sec_mu_f):
                sec_str = f"Y (at {sec_mu_f:.0f} MeV)"
            elif sec_present == "Y":
                sec_str = "Y"
            else:
                sec_str = "N"

            title = (
                f"Run {run:06d}  Split {split:05d}"
                f"   \u2014   \u03c3_E/E = {sigma_str}"
                f"   \u2014   Secondary peak: {sec_str}"
            )

            fig = plt.figure(figsize=(11, 8.5))

            if png_p.exists():
                img = plt.imread(str(png_p))
                ax = fig.add_axes([0.0, 0.0, 1.0, 0.93])   # leave room for title
                ax.imshow(img)
                ax.axis("off")
            else:
                ax = fig.add_axes([0.0, 0.0, 1.0, 0.93])
                ax.set_facecolor("#f0f0f0")
                ax.text(
                    0.5, 0.5,
                    f"Missing:\nscan_{run:06d}/split_{split:05d}\n"
                    f"result_calibrated.png",
                    ha="center", va="center",
                    fontsize=18, color="#880000",
                    transform=ax.transAxes,
                )
                ax.axis("off")

            fig.suptitle(title, fontsize=9, y=0.98)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"\nWrote combined calibrated PDF: {pdf_path}", flush=True)


# ===========================================================================
# §8 / §9 — CLI + main
# ===========================================================================

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--evio-dir", default=_DEFAULT_EVIO_DIR,
        help=f"Directory containing prad_*.evio.* files "
             f"(default: {_DEFAULT_EVIO_DIR})",
    )
    ap.add_argument(
        "--output-dir", default=_DEFAULT_OUTPUT_DIR,
        help=f"Output base directory (default: {_DEFAULT_OUTPUT_DIR})",
    )
    ap.add_argument(
        "--max-events", type=int, default=_DEFAULT_MAX_EVENTS,
        help=(
            "Maximum events to process per EVIO file. "
            "Default: None (process full file — no cap). "
            "Set an explicit integer to impose a limit."
        ),
    )
    ap.add_argument(
        "--prominence-threshold", type=float, default=_DEFAULT_PROMINENCE,
        help=f"scipy.signal.find_peaks prominence threshold "
             f"(default: {_DEFAULT_PROMINENCE})",
    )
    ap.add_argument(
        "--peak-distance", type=float, default=_DEFAULT_PEAK_DIST_MEV,
        help=f"Minimum distance between detected peaks in MeV "
             f"(default: {_DEFAULT_PEAK_DIST_MEV})",
    )
    ap.add_argument(
        "--parallel", type=int, default=_DEFAULT_PARALLEL,
        help=f"Number of parallel benchmark invocations "
             f"(default: {_DEFAULT_PARALLEL}; recommend ≤4)",
    )
    ap.add_argument(
        "--skip-benchmark", action="store_true",
        help="Reuse existing TSV files from a previous run (re-plot / re-analyse only).",
    )
    return ap


def main(argv=None) -> int:
    ap    = build_parser()
    args  = ap.parse_args(argv)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- Logging setup -------------------------------------------------------
    log_path = output_dir / "benchmark_hycal_timing_scan.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-7s  %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
        handlers=[
            logging.FileHandler(log_path, mode="a", encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )

    start_ts = datetime.datetime.now().isoformat(timespec="seconds")
    log.info("=" * 70)
    log.info("benchmark_hycal_timing_scan  started %s", start_ts)
    log.info("  evio_dir    : %s", args.evio_dir)
    log.info("  output_dir  : %s", output_dir)
    log.info("  max_events  : %s",
             str(args.max_events) if args.max_events is not None else "full file")
    log.info("  prominence  : %.1f", args.prominence_threshold)
    log.info("  peak_dist   : %.1f MeV", args.peak_distance)
    log.info("  parallel    : %d", args.parallel)
    log.info("  skip_bench  : %s", args.skip_benchmark)
    log.info("=" * 70)

    if not HAVE_SCIPY:
        log.warning(
            "scipy not available — peak detection disabled (install scipy)."
        )

    # ---- File discovery -----------------------------------------------------
    evio_dir = Path(args.evio_dir)
    if not evio_dir.is_dir():
        log.error("EVIO directory not found: %s", evio_dir)
        return 1

    run_map = discover_files(evio_dir)
    if not run_map:
        log.error("No prad_*.evio.* files found in %s", evio_dir)
        return 1

    total_files = sum(len(splits) for splits in run_map.values())
    log.info(
        "Discovered %d run(s), %d total files: %s",
        len(run_map), total_files,
        {r: len(s) for r, s in sorted(run_map.items())},
    )

    # ---- Build work list (run, split, path) ---------------------------------
    work_items: list[tuple[int, int, Path]] = []
    for run in sorted(run_map.keys()):
        for split, fpath in run_map[run]:
            work_items.append((run, split, fpath))

    # ---- Run benchmark (serial or parallel) ---------------------------------
    def _run_one_benchmark(item: tuple[int, int, Path]) -> tuple[int, int, Path, Optional[Path], Optional[str]]:
        run, split, fpath = item
        try:
            tsv = run_benchmark(
                evio_file=fpath,
                run=run,
                split=split,
                output_dir=output_dir,
                max_events=args.max_events,
                skip_benchmark=args.skip_benchmark,
            )
            return run, split, fpath, tsv, None
        except Exception as exc:
            err_msg = str(exc)
            log.error(
                "SKIP run=%d split=%05d: %s", run, split, err_msg
            )
            return run, split, fpath, None, err_msg

    benchmark_results: list[tuple[int, int, Path, Optional[Path], Optional[str]]] = []
    if args.parallel <= 1:
        for item in work_items:
            benchmark_results.append(_run_one_benchmark(item))
    else:
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=args.parallel
        ) as pool:
            futures = {pool.submit(_run_one_benchmark, item): item
                       for item in work_items}
            for fut in concurrent.futures.as_completed(futures):
                benchmark_results.append(fut.result())
        # Sort back into stable order (run, split)
        benchmark_results.sort(key=lambda t: (t[0], t[1]))

    # ---- Extract per-file aggregate rows ------------------------------------
    print("\nAggregating results...", flush=True)
    all_rows: list[dict] = []
    errors:   list[str]  = []

    for run, split, fpath, tsv_path, err_msg in benchmark_results:
        if tsv_path is None or err_msg is not None:
            errors.append(f"run={run} split={split}: {err_msg}")
            continue
        try:
            row = extract_row(
                run=run,
                split=split,
                evio_file=fpath,
                tsv_path=tsv_path,
                prominence_threshold=args.prominence_threshold,
                min_dist_mev=args.peak_distance,
            )
            all_rows.append(row)
            n_ev = row.get("n_events", float("nan"))
            n_ev_str = (
                f"{int(n_ev):,}" if not math.isnan(float(n_ev)) else "unknown"
            )
            log.info(
                "  run=%d split=%05d  events=%s  σ_E/E=%.3f%%  secondary-peak=%s",
                run, split,
                n_ev_str,
                row["sigma_over_E_CB"] * 100.0
                if not math.isnan(float(row["sigma_over_E_CB"])) else float("nan"),
                row["secondary_peak_present"],
            )
        except Exception as exc:
            msg = f"run={run} split={split} TSV parse failed: {exc}"
            log.error(msg)
            errors.append(msg)

    if not all_rows:
        log.error("No rows collected — all benchmark runs failed or TSVs missing.")
        return 1

    # Sort by run then split
    all_rows.sort(key=lambda r: (int(r["run_number"]), int(r["split_number"])))

    # ---- Write aggregate TSV ------------------------------------------------
    summary_tsv = output_dir / "benchmark_hycal_timing_scan_summary.tsv"
    write_summary_tsv(all_rows, summary_tsv)
    print(f"\nWrote summary TSV: {summary_tsv}", flush=True)
    log.info("Wrote summary TSV: %s", summary_tsv)

    # ---- Trend plot ---------------------------------------------------------
    trend_png = output_dir / "benchmark_hycal_timing_scan_sigma_E_vs_file.png"
    try:
        make_trend_plot(all_rows, trend_png)
        if not trend_png.exists():
            log.error("Trend plot not generated: %s", trend_png)
            return 1
        log.info("Wrote trend plot: %s", trend_png)
    except Exception as exc:
        log.error("Trend plot generation failed: %s", exc)
        return 1

    # ---- Combined calibrated-plot PDF ----------------------------------------
    combined_pdf = (
        output_dir / "benchmark_hycal_timing_scan_all_calibrated_plots.pdf"
    )
    try:
        generate_combined_pdf(all_rows, output_dir, combined_pdf)
        log.info("Wrote combined calibrated PDF: %s", combined_pdf)
    except Exception as exc:
        log.error("Combined PDF generation failed: %s", exc)
        # Non-fatal: log and continue

    # ---- Findings markdown --------------------------------------------------
    findings_md = output_dir / "benchmark_hycal_timing_scan_findings.md"
    write_findings_md(all_rows, findings_md)
    log.info("Wrote findings: %s", findings_md)

    # ---- Final summary ------------------------------------------------------
    n_secondary = sum(
        1 for r in all_rows if r.get("secondary_peak_present") == "Y"
    )
    n_runs      = len(run_map)
    end_ts      = datetime.datetime.now().isoformat(timespec="seconds")
    summary_msg = (
        f"\n[DONE]  {end_ts}\n"
        f"  {n_runs} run(s) × {len(work_items)} files = "
        f"{len(all_rows)} files processed, "
        f"{n_secondary} files with secondary peaks.\n"
    )
    if errors:
        summary_msg += f"  {len(errors)} error(s) (see log for details).\n"
    print(summary_msg, flush=True)
    log.info(summary_msg.strip())

    return 0 if not errors else 2


if __name__ == "__main__":
    sys.exit(main())
