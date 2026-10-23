#!/usr/bin/env python3
"""batch_bimodality_scan.py — Per-split bimodality analysis for 2.2 GeV H2 runs.

Runs bimodality analysis for every split of the 2.2 GeV H2 runs 25260, 25299,
and 25320, then consolidates all per-split bimodality maps into a single
comparison PDF.

Workflow per (run, split) pair
-------------------------------
Step 1: Generate clusters TSV (benchmark_hycal_timing.py --write-cluster-tsv)
Step 2: Run bimodality scan (scan_module_bimodality.py)

Output layout
-------------
  output/bimodality_scan/{run:06d}/split_{split:05d}/
      benchmark_result.tsv
      benchmark_result.png
      benchmark_result_calibrated.png
      benchmark_result_clusters.tsv          ← input for Step 2
      module_bimodality_{run}_split{split:05d}.tsv
      module_bimodality_{run}_split{split:05d}_map.png
      module_bimodality_{run}_split{split:05d}_spatial.png
      module_bimodality_{run}_split{split:05d}_top20.pdf
      module_bimodality_{run}_split{split:05d}_findings.md

  output/bimodality_scan/combined_maps_2p2_h2.pdf  ← multi-page comparison
  output/bimodality_scan/batch_bimodality_scan.log  ← persistent log

Usage
-----
  python analysis/pyscripts/batch_bimodality_scan.py [options]
  python analysis/pyscripts/batch_bimodality_scan.py --parallel 4
  python analysis/pyscripts/batch_bimodality_scan.py --skip-existing

Reference: user specification (inline, Oct 2026).
"""

from __future__ import annotations

import argparse
import concurrent.futures
import datetime
import glob
import logging
import math
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

try:
    import numpy as np
    HAVE_NUMPY = True
except ImportError:
    HAVE_NUMPY = False

try:
    import pandas as pd
    HAVE_PANDAS = True
except ImportError:
    HAVE_PANDAS = False

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RUNS_2P2_H2 = [25260, 25299, 25320]
EVIO_DIR = Path("/Users/L00338853/work/PRad/data/evio")

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_SCRIPT_DIR = Path(__file__).resolve().parent

_BENCHMARK_SCRIPT = str(_SCRIPT_DIR / "benchmark_hycal_timing.py")
_BIMODALITY_SCRIPT = str(_SCRIPT_DIR / "scan_module_bimodality.py")

_DEFAULT_OUTPUT_DIR = str(_SCRIPT_DIR / "output" / "bimodality_scan")
_DEFAULT_PARALLEL = 8
_DEFAULT_TIMEOUT = 1200  # seconds (20 min)

# DAQ config / HyCal map from environment (set by setup_env.sh)
_DAQ_CONFIG = os.environ.get("DAQ_CONFIG", "")
_HC_MAP = os.environ.get("HC_MAP", "")

# Benchmark parameters for 2.2 GeV H2 runs
_TARGET_PEAK = 2239.5
_PRE_WINDOW = (100.0, 200.0)
_WINDOW = 16
_SIGNAL_MIN = 500
_E_RANGE = (500, 3500)
_BIN_WIDTH = 20
_EXCLUDE_ROWCOL = (15, 20, 15, 20)
_INNER_RADIUS = 200
_MIN_SEED_ENERGY = 200
_MIN_PER_MODULE = 30

# Bimodality scan parameters
_BIMODAL_MIN_EVENTS = 100
_BIMODAL_RATIO = 0.30

# EVIO filename pattern
_EVIO_PAT = re.compile(r"^prad_0*(\d+)\.evio\.(\d+)$", re.IGNORECASE)

# Per-run colors for the bar chart
_RUN_COLORS = {
    25260: "#1f77b4",  # blue
    25299: "#ff7f0e",  # orange
    25320: "#2ca02c",  # green
}

log = logging.getLogger("batch_bimodality")


# ---------------------------------------------------------------------------
# Split discovery
# ---------------------------------------------------------------------------

def discover_splits(evio_dir: Path, runs: list[int]) -> dict[int, list[tuple[int, Path]]]:
    """Return {run: [(split_int, Path), ...]} sorted by split, for the given runs."""
    result: dict[int, list[tuple[int, Path]]] = {r: [] for r in runs}
    pattern = str(evio_dir / "prad_*.evio.*")
    for fpath in sorted(glob.glob(pattern)):
        p = Path(fpath)
        m = _EVIO_PAT.match(p.name)
        if not m:
            continue
        run = int(m.group(1))
        split = int(m.group(2))
        if run in result:
            result[run].append((split, p))
    for run in result:
        result[run].sort(key=lambda t: t[0])
    return result


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def split_output_dir(output_dir: Path, run: int, split: int) -> Path:
    return output_dir / f"{run:06d}" / f"split_{split:05d}"


def benchmark_stem(output_dir: Path, run: int, split: int) -> str:
    return str(split_output_dir(output_dir, run, split) / "benchmark_result")


def clusters_tsv_path(output_dir: Path, run: int, split: int) -> Path:
    return split_output_dir(output_dir, run, split) / "benchmark_result_clusters.tsv"


def bimodality_tsv_path(output_dir: Path, run: int, split: int) -> Path:
    label = f"{run}_split{split:05d}"
    return split_output_dir(output_dir, run, split) / f"module_bimodality_{label}.tsv"


def bimodality_map_path(output_dir: Path, run: int, split: int) -> Path:
    label = f"{run}_split{split:05d}"
    return split_output_dir(output_dir, run, split) / f"module_bimodality_{label}_map.png"


# ---------------------------------------------------------------------------
# Step 1: Run benchmark_hycal_timing.py --write-cluster-tsv
# ---------------------------------------------------------------------------

def run_benchmark(
    *,
    evio_file: Path,
    run: int,
    split: int,
    output_dir: Path,
    skip_existing: bool,
    timeout_seconds: int,
) -> tuple[bool, float, str]:
    """Run benchmark_hycal_timing.py for one split.

    Returns (success, elapsed_seconds, error_message).
    """
    out_stem = benchmark_stem(output_dir, run, split)
    clusters_tsv = clusters_tsv_path(output_dir, run, split)
    log_path = split_output_dir(output_dir, run, split) / "benchmark.log"

    if skip_existing and clusters_tsv.exists():
        log.info("skip benchmark (clusters TSV exists): %s", clusters_tsv)
        return True, 0.0, ""

    split_output_dir(output_dir, run, split).mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        _BENCHMARK_SCRIPT,
        str(evio_file),
        out_stem,
        "--pre-window", str(_PRE_WINDOW[0]), str(_PRE_WINDOW[1]),
        "--window", str(_WINDOW),
        "--signal-min", str(_SIGNAL_MIN),
        "--e-range", str(_E_RANGE[0]), str(_E_RANGE[1]),
        "--bin-width", str(_BIN_WIDTH),
        "--calibrate",
        "--target-peak", str(_TARGET_PEAK),
        "--cal-method", "cb-peak",
        "--exclude-rowcol",
            str(_EXCLUDE_ROWCOL[0]), str(_EXCLUDE_ROWCOL[1]),
            str(_EXCLUDE_ROWCOL[2]), str(_EXCLUDE_ROWCOL[3]),
        "--inner-radius", str(_INNER_RADIUS),
        "--min-seed-energy", str(_MIN_SEED_ENERGY),
        "--min-per-module", str(_MIN_PER_MODULE),
        "--write-cluster-tsv",
    ]
    if _DAQ_CONFIG:
        cmd.extend(["--daq-config", _DAQ_CONFIG])
    if _HC_MAP:
        cmd.extend(["--hc-map-file", _HC_MAP])

    t0 = time.monotonic()
    try:
        with open(log_path, "w") as lf:
            try:
                result = subprocess.run(
                    cmd,
                    stdout=lf,
                    stderr=subprocess.STDOUT,
                    check=False,
                    timeout=timeout_seconds if timeout_seconds > 0 else None,
                )
            except subprocess.TimeoutExpired:
                elapsed = time.monotonic() - t0
                msg = (
                    f"benchmark TIMEOUT after {elapsed:.1f}s "
                    f"(limit: {timeout_seconds}s) run={run} split={split}"
                )
                log.error(msg)
                with open(log_path, "a") as lf2:
                    lf2.write(f"\n\n[ORCHESTRATOR] TIMEOUT after {elapsed:.1f}s\n")
                return False, elapsed, msg
        elapsed = time.monotonic() - t0
        if result.returncode != 0:
            msg = (
                f"benchmark failed (rc={result.returncode}) "
                f"run={run} split={split}. Log: {log_path}"
            )
            log.error(msg)
            return False, elapsed, msg
    except FileNotFoundError as exc:
        elapsed = time.monotonic() - t0
        msg = f"benchmark script not found: {_BENCHMARK_SCRIPT}: {exc}"
        log.error(msg)
        return False, elapsed, msg

    elapsed = time.monotonic() - t0

    if not clusters_tsv.exists():
        msg = (
            f"Expected clusters TSV not created: {clusters_tsv} "
            f"(log: {log_path})"
        )
        log.error(msg)
        return False, elapsed, msg

    log.info("  benchmark done in %.1f s → %s", elapsed, clusters_tsv)
    return True, elapsed, ""


# ---------------------------------------------------------------------------
# Step 2: Run scan_module_bimodality.py
# ---------------------------------------------------------------------------

def run_bimodality_scan(
    *,
    run: int,
    split: int,
    output_dir: Path,
    skip_existing: bool,
    timeout_seconds: int,
) -> tuple[bool, str, int, int]:
    """Run scan_module_bimodality.py for one split.

    Returns (success, error_message, n_analyzed, n_bimodal).
    """
    label = f"{run}_split{split:05d}"
    input_tsv = clusters_tsv_path(output_dir, run, split)
    out_dir = split_output_dir(output_dir, run, split)
    result_tsv = bimodality_tsv_path(output_dir, run, split)
    log_path = out_dir / "bimodality.log"

    if skip_existing and result_tsv.exists():
        log.info("skip bimodality scan (TSV exists): %s", result_tsv)
        n_analyzed, n_bimodal = _read_bimodality_counts(result_tsv)
        return True, "", n_analyzed, n_bimodal

    if not input_tsv.exists():
        msg = f"clusters TSV missing, cannot run bimodality: {input_tsv}"
        log.error(msg)
        return False, msg, 0, 0

    cmd = [
        sys.executable,
        _BIMODALITY_SCRIPT,
        "--input-tsv", str(input_tsv),
        "--output-dir", str(out_dir),
        "--run-label", label,
        "--target-peak", str(_TARGET_PEAK),
        "--e-range", str(_E_RANGE[0]), str(_E_RANGE[1]),
        "--min-events", str(_BIMODAL_MIN_EVENTS),
        "--bimodal-ratio", str(_BIMODAL_RATIO),
    ]

    try:
        with open(log_path, "w") as lf:
            try:
                proc = subprocess.run(
                    cmd,
                    stdout=lf,
                    stderr=subprocess.STDOUT,
                    check=False,
                    timeout=timeout_seconds if timeout_seconds > 0 else None,
                )
            except subprocess.TimeoutExpired:
                msg = (
                    f"bimodality scan TIMEOUT after {timeout_seconds}s "
                    f"run={run} split={split}"
                )
                log.error(msg)
                return False, msg, 0, 0
        if proc.returncode != 0:
            msg = (
                f"bimodality scan failed (rc={proc.returncode}) "
                f"run={run} split={split}. Log: {log_path}"
            )
            log.error(msg)
            return False, msg, 0, 0
    except FileNotFoundError as exc:
        msg = f"bimodality script not found: {_BIMODALITY_SCRIPT}: {exc}"
        log.error(msg)
        return False, msg, 0, 0

    n_analyzed, n_bimodal = _read_bimodality_counts(result_tsv)
    return True, "", n_analyzed, n_bimodal


def _read_bimodality_counts(tsv_path: Path) -> tuple[int, int]:
    """Read n_analyzed and n_bimodal from a module_bimodality TSV.

    Returns (0, 0) if file missing or unreadable.
    """
    if not tsv_path.exists():
        return 0, 0
    if not HAVE_PANDAS:
        # Fallback: count lines manually
        try:
            lines = tsv_path.read_text(encoding="utf-8").splitlines()
            data = [l for l in lines[1:] if l.strip()]
            n_bimodal = sum(1 for l in data if l.split("\t")[-1].strip() == "Y")
            return len(data), n_bimodal
        except Exception:
            return 0, 0
    try:
        df = pd.read_csv(str(tsv_path), sep="\t", dtype=str, low_memory=False)
        n_analyzed = len(df)
        n_bimodal = int((df.get("is_bimodal", pd.Series(dtype=str)) == "Y").sum())
        return n_analyzed, n_bimodal
    except Exception:
        return 0, 0


# ---------------------------------------------------------------------------
# Per-split processing (benchmark + bimodality)
# ---------------------------------------------------------------------------

class SplitResult:
    """Holds the outcome of processing one (run, split) pair."""

    def __init__(self, run: int, split: int, evio_file: Path):
        self.run = run
        self.split = split
        self.evio_file = evio_file
        self.benchmark_ok = False
        self.benchmark_elapsed = 0.0
        self.benchmark_error = ""
        self.bimodality_ok = False
        self.bimodality_error = ""
        self.n_analyzed = 0
        self.n_bimodal = 0

    @property
    def bimodality_pct(self) -> float:
        if self.n_analyzed > 0:
            return 100.0 * self.n_bimodal / self.n_analyzed
        return float("nan")

    @property
    def completed(self) -> bool:
        return self.benchmark_ok and self.bimodality_ok

    @property
    def label(self) -> str:
        return f"run {self.run} split {self.split:05d}"


def process_split(
    run: int,
    split: int,
    evio_file: Path,
    output_dir: Path,
    skip_existing: bool,
    timeout_seconds: int,
) -> SplitResult:
    """Run benchmark + bimodality for one split. Returns SplitResult."""
    res = SplitResult(run, split, evio_file)
    msg_prefix = f"Processing run {run} split {split:05d}"
    print(f"{msg_prefix}...", flush=True)
    log.info("%s: EVIO=%s", msg_prefix, evio_file)

    # Step 1: benchmark
    ok, elapsed, err = run_benchmark(
        evio_file=evio_file,
        run=run,
        split=split,
        output_dir=output_dir,
        skip_existing=skip_existing,
        timeout_seconds=timeout_seconds,
    )
    res.benchmark_ok = ok
    res.benchmark_elapsed = elapsed
    res.benchmark_error = err
    if ok:
        skipped = skip_existing and elapsed == 0.0
        if skipped:
            print(f"  benchmark skipped (output exists)", flush=True)
        else:
            print(f"  benchmark completed in {elapsed:.0f}s", flush=True)
        log.info("  %s: benchmark ok in %.1fs", msg_prefix, elapsed)
    else:
        print(f"  benchmark FAILED: {err}", flush=True)
        log.error("  %s: benchmark failed: %s", msg_prefix, err)
        return res

    # Step 2: bimodality scan
    ok2, err2, n_analyzed, n_bimodal = run_bimodality_scan(
        run=run,
        split=split,
        output_dir=output_dir,
        skip_existing=skip_existing,
        timeout_seconds=timeout_seconds,
    )
    res.bimodality_ok = ok2
    res.bimodality_error = err2
    res.n_analyzed = n_analyzed
    res.n_bimodal = n_bimodal
    if ok2:
        pct = res.bimodality_pct
        pct_str = f"{pct:.0f}%" if not math.isnan(pct) else "N/A"
        print(
            f"  bimodality: {n_analyzed} modules analyzed, "
            f"{n_bimodal} bimodal ({pct_str})",
            flush=True,
        )
        log.info(
            "  %s: bimodality ok — %d analyzed, %d bimodal (%.1f%%)",
            msg_prefix, n_analyzed, n_bimodal,
            pct if not math.isnan(pct) else 0.0,
        )
    else:
        print(f"  bimodality scan FAILED: {err2}", flush=True)
        log.error("  %s: bimodality scan failed: %s", msg_prefix, err2)

    return res


# ---------------------------------------------------------------------------
# Combined comparison PDF
# ---------------------------------------------------------------------------

def build_combined_pdf(
    results: list[SplitResult],
    output_dir: Path,
    pdf_path: Path,
) -> None:
    """Build a multi-page PDF with one page per split's bimodality map.

    First page: summary with bar chart.
    Subsequent pages: full-size map image (14x10 inches each).
    Pages ordered by run then split.
    """
    pdf_path.parent.mkdir(parents=True, exist_ok=True)

    # Sort defensively
    sorted_results = sorted(results, key=lambda r: (r.run, r.split))
    completed = [r for r in sorted_results if r.bimodality_ok]
    n_attempted = len(sorted_results)
    n_completed = len(completed)
    n_failed = n_attempted - n_completed

    # Gather per-run stats for bar chart
    run_groups: dict[int, list[SplitResult]] = {}
    for r in sorted_results:
        run_groups.setdefault(r.run, []).append(r)

    with PdfPages(str(pdf_path)) as pdf:
        # ---- Page 1: Summary -------------------------------------------------
        _write_summary_page(pdf, sorted_results, run_groups,
                             n_attempted, n_completed, n_failed)

        # ---- Pages 2+: one map per split ------------------------------------
        for r in sorted_results:
            map_path = bimodality_map_path(output_dir, r.run, r.split)
            _write_map_page(pdf, r, map_path)

    log.info("Combined PDF written: %s", pdf_path)
    print(f"\nCombined PDF: {pdf_path}", flush=True)


def _write_summary_page(
    pdf: PdfPages,
    sorted_results: list[SplitResult],
    run_groups: dict[int, list[SplitResult]],
    n_attempted: int,
    n_completed: int,
    n_failed: int,
) -> None:
    """Write the first summary page with a bar chart."""
    fig = plt.figure(figsize=(14, 10))

    # Text summary (top half)
    text_ax = fig.add_axes([0.05, 0.55, 0.90, 0.38])
    text_ax.axis("off")

    lines = [
        f"Batch Bimodality Scan — 2.2 GeV H2 Runs (25260, 25299, 25320)",
        f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"",
        f"Total splits attempted: {n_attempted}   |   Completed: {n_completed}   "
        f"|   Failed/skipped: {n_failed}",
        f"",
    ]

    # Per-run breakdown
    for run in RUNS_2P2_H2:
        rlist = run_groups.get(run, [])
        done = [r for r in rlist if r.bimodality_ok]
        with_bimodal = [r for r in done if r.n_bimodal > 0]
        pcts = [r.bimodality_pct for r in done if not math.isnan(r.bimodality_pct)]
        if pcts and HAVE_NUMPY:
            median_pct = float(np.median(pcts))
            med_str = f"{median_pct:.1f}%"
        else:
            med_str = "N/A"
        lines.append(
            f"  Run {run}: {len(rlist)} splits attempted, {len(done)} completed, "
            f"{len(with_bimodal)} with any bimodality, median %bimodal = {med_str}"
        )

    text_ax.text(
        0.0, 1.0,
        "\n".join(lines),
        transform=text_ax.transAxes,
        fontsize=11,
        va="top", ha="left",
        family="monospace",
    )

    # Bar chart (bottom half): bimodality % per split, grouped by run
    bar_ax = fig.add_axes([0.08, 0.08, 0.88, 0.42])
    _draw_bar_chart(bar_ax, sorted_results, run_groups)

    fig.suptitle(
        "Per-Split Bimodality Summary — 2.2 GeV H2 Runs",
        fontsize=14, fontweight="bold", y=0.97,
    )
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def _draw_bar_chart(
    ax: plt.Axes,
    sorted_results: list[SplitResult],
    run_groups: dict[int, list[SplitResult]],
) -> None:
    """Draw a bar chart: bimodality % per split, color-coded by run."""
    if not HAVE_NUMPY:
        ax.text(0.5, 0.5, "numpy unavailable — bar chart skipped",
                ha="center", va="center", transform=ax.transAxes, fontsize=10)
        return

    x_positions = []
    y_values = []
    bar_colors = []
    x_labels = []

    gap_between_runs = 0.5
    x = 0.0
    bar_width = 0.7

    for run in RUNS_2P2_H2:
        rlist = run_groups.get(run, [])
        color = _RUN_COLORS.get(run, "gray")
        first_in_run = True
        for r in rlist:
            pct = r.bimodality_pct if r.bimodality_ok else float("nan")
            x_positions.append(x)
            y_values.append(pct if not math.isnan(pct) else 0.0)
            bar_colors.append(color if not math.isnan(pct) else "#cccccc")
            lbl = f"{r.split:05d}"
            if first_in_run:
                lbl = f"{run}\n{lbl}"
                first_in_run = False
            x_labels.append(lbl)
            x += 1.0
        x += gap_between_runs

    if not x_positions:
        ax.text(0.5, 0.5, "No splits processed.", ha="center", va="center",
                transform=ax.transAxes, fontsize=10)
        return

    bars = ax.bar(x_positions, y_values, width=bar_width, color=bar_colors,
                  edgecolor="black", linewidth=0.4, alpha=0.85)

    # Value labels on bars
    for bar, y in zip(bars, y_values):
        if y > 0:
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                bar.get_height() + 0.3,
                f"{y:.0f}%",
                ha="center", va="bottom", fontsize=6, rotation=90,
            )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("% modules bimodal", fontsize=10)
    ax.set_title(
        "Bimodality % per split  (grey = failed/missing)",
        fontsize=10,
    )
    ax.grid(alpha=0.3, axis="y")
    ax.set_ylim(0, max(max(y_values) * 1.15, 5.0))

    # Legend
    import matplotlib.patches as mpatches
    legend_patches = [
        mpatches.Patch(facecolor=_RUN_COLORS.get(r, "gray"),
                       edgecolor="black", label=f"Run {r}")
        for r in RUNS_2P2_H2
    ]
    legend_patches.append(
        mpatches.Patch(facecolor="#cccccc", edgecolor="black",
                       label="Failed/missing")
    )
    ax.legend(handles=legend_patches, loc="upper right", fontsize=9)


def _write_map_page(
    pdf: PdfPages,
    r: SplitResult,
    map_path: Path,
) -> None:
    """Write one page: the bimodality map for a split, or a placeholder."""
    fig = plt.figure(figsize=(14, 10))

    # Build page title
    if r.bimodality_ok and not math.isnan(r.bimodality_pct):
        title = (
            f"Run {r.run}  Split {r.split:05d}  \u2014  "
            f"{r.n_analyzed} modules analyzed, "
            f"{r.n_bimodal} bimodal ({r.bimodality_pct:.0f}%)"
        )
    else:
        title = f"Run {r.run}  Split {r.split:05d}  \u2014  (failed / no data)"

    if map_path.exists():
        try:
            img = plt.imread(str(map_path))
            ax = fig.add_axes([0.0, 0.0, 1.0, 0.93])
            ax.imshow(img)
            ax.axis("off")
        except Exception as exc:
            ax = fig.add_axes([0.0, 0.0, 1.0, 0.93])
            ax.set_facecolor("#f8f8f8")
            ax.text(
                0.5, 0.5,
                f"Error loading image:\n{map_path.name}\n{exc}",
                ha="center", va="center", fontsize=14, color="#880000",
                transform=ax.transAxes,
            )
            ax.axis("off")
    else:
        ax = fig.add_axes([0.0, 0.0, 1.0, 0.93])
        ax.set_facecolor("#f0f0f0")
        if r.benchmark_error:
            detail = f"Benchmark failed:\n{r.benchmark_error}"
        elif r.bimodality_error:
            detail = f"Bimodality scan failed:\n{r.bimodality_error}"
        else:
            detail = "Map PNG not found."
        ax.text(
            0.5, 0.5,
            f"Missing: run {r.run} split {r.split:05d}\n{detail}",
            ha="center", va="center",
            fontsize=14, color="#880000",
            transform=ax.transAxes,
        )
        ax.axis("off")

    fig.suptitle(title, fontsize=11, y=0.98)
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------

def setup_logging(log_path: Path) -> None:
    """Configure root logger with a file handler and a stdout handler."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fmt = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(
        level=logging.DEBUG,
        format=fmt,
        handlers=[
            logging.FileHandler(str(log_path), mode="a", encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    # Quieten the stdout handler so INFO+ goes to file but only WARNING+ goes
    # to stdout (prints already handle the progress output).
    for h in logging.root.handlers:
        if isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler):
            h.setLevel(logging.WARNING)


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Runs are hardcoded to 25260, 25299, 25320 (2.2 GeV H2 dataset).\n"
            "Splits are discovered from EVIO_DIR = "
            f"{EVIO_DIR}\n"
        ),
    )
    ap.add_argument(
        "--parallel",
        type=int,
        default=_DEFAULT_PARALLEL,
        metavar="N",
        help=(
            "Maximum number of benchmark subprocesses to run concurrently "
            f"(default: {_DEFAULT_PARALLEL}).  The bimodality scan step is "
            "lightweight and always runs sequentially after its benchmark."
        ),
    )
    ap.add_argument(
        "--skip-existing",
        action="store_true",
        default=False,
        help=(
            "Skip steps whose output files already exist. "
            "Default: False (overwrite / re-run everything)."
        ),
    )
    ap.add_argument(
        "--timeout-seconds",
        type=int,
        default=_DEFAULT_TIMEOUT,
        metavar="T",
        help=(
            "Per-benchmark subprocess timeout in seconds "
            f"(default: {_DEFAULT_TIMEOUT}).  "
            "Set 0 to disable."
        ),
    )
    ap.add_argument(
        "--output-dir",
        default=_DEFAULT_OUTPUT_DIR,
        metavar="DIR",
        help=(
            "Base output directory for all bimodality scan outputs "
            f"(default: {_DEFAULT_OUTPUT_DIR})."
        ),
    )
    return ap


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = build_parser()
    args = ap.parse_args(argv)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log_path = output_dir / "batch_bimodality_scan.log"
    setup_logging(log_path)

    log.info(
        "batch_bimodality_scan starting — runs=%s  parallel=%d  "
        "skip_existing=%s  timeout=%ds  output_dir=%s",
        RUNS_2P2_H2,
        args.parallel,
        args.skip_existing,
        args.timeout_seconds,
        output_dir,
    )
    print(
        f"\n{'='*70}\n"
        f"batch_bimodality_scan.py\n"
        f"Runs: {RUNS_2P2_H2}\n"
        f"EVIO dir: {EVIO_DIR}\n"
        f"Output dir: {output_dir}\n"
        f"Parallel: {args.parallel}  |  Skip existing: {args.skip_existing}  "
        f"|  Timeout: {args.timeout_seconds}s\n"
        f"{'='*70}\n",
        flush=True,
    )

    # Discover splits
    splits_by_run = discover_splits(EVIO_DIR, RUNS_2P2_H2)
    all_work: list[tuple[int, int, Path]] = []  # (run, split, evio_path)
    for run in RUNS_2P2_H2:
        splits = splits_by_run.get(run, [])
        if not splits:
            log.warning("No EVIO splits found for run %d in %s", run, EVIO_DIR)
            print(
                f"[WARN] No EVIO splits found for run {run} in {EVIO_DIR}",
                flush=True,
            )
        for split, path in splits:
            all_work.append((run, split, path))

    n_total = len(all_work)
    if n_total == 0:
        print("[ERROR] No splits found for any run. Exiting.", flush=True)
        log.error("No splits found for any run. Exiting.")
        return 1

    print(f"Found {n_total} total splits to process.\n", flush=True)

    # Process splits — benchmark runs in parallel, bimodality scan follows
    # each benchmark sequentially within the same worker thread.
    all_results: list[SplitResult] = []
    errors: list[str] = []

    def _worker(work: tuple[int, int, Path]) -> SplitResult:
        run, split, evio_file = work
        return process_split(
            run=run,
            split=split,
            evio_file=evio_file,
            output_dir=output_dir,
            skip_existing=args.skip_existing,
            timeout_seconds=args.timeout_seconds,
        )

    n_workers = max(1, min(args.parallel, n_total))
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as pool:
        future_to_work = {pool.submit(_worker, w): w for w in all_work}
        for future in concurrent.futures.as_completed(future_to_work):
            work = future_to_work[future]
            run, split, _ = work
            try:
                res = future.result()
            except Exception as exc:
                msg = f"Unexpected exception processing run {run} split {split}: {exc}"
                log.error(msg)
                print(f"[ERROR] {msg}", flush=True)
                errors.append(msg)
                # Create a failed placeholder result
                res = SplitResult(run, split, work[2])
                res.benchmark_error = str(exc)
            all_results.append(res)
            if not res.benchmark_ok and res.benchmark_error:
                errors.append(
                    f"benchmark run={run} split={split}: {res.benchmark_error}"
                )
            if res.benchmark_ok and not res.bimodality_ok and res.bimodality_error:
                errors.append(
                    f"bimodality run={run} split={split}: {res.bimodality_error}"
                )

    # Sort by (run, split)
    all_results.sort(key=lambda r: (r.run, r.split))

    # Build combined PDF
    pdf_path = output_dir / "combined_maps_2p2_h2.pdf"
    try:
        build_combined_pdf(all_results, output_dir, pdf_path)
    except Exception as exc:
        msg = f"Failed to build combined PDF: {exc}"
        log.error(msg)
        print(f"[ERROR] {msg}", flush=True)
        errors.append(msg)

    # Final summary
    n_completed = sum(1 for r in all_results if r.completed)
    n_bench_ok = sum(1 for r in all_results if r.benchmark_ok)
    n_bm_failed = sum(
        1 for r in all_results if r.benchmark_ok and not r.bimodality_ok
    )
    n_bench_failed = sum(1 for r in all_results if not r.benchmark_ok)
    n_skipped = sum(
        1 for r in all_results
        if (r.benchmark_ok and r.benchmark_elapsed == 0.0)
    )

    summary = (
        f"\n{'='*70}\n"
        f"[DONE] {n_total} splits attempted\n"
        f"       {n_completed} fully completed (benchmark + bimodality)\n"
        f"       {n_bench_failed} benchmark failed/timed-out\n"
        f"       {n_bm_failed} bimodality scan failed after benchmark\n"
        f"       {n_skipped} skipped (output existed)\n"
    )
    if errors:
        summary += f"\n{len(errors)} error(s) collected:\n"
        for i, e in enumerate(errors, 1):
            summary += f"  [{i}] {e}\n"
    summary += f"\nCombined PDF: {pdf_path}\n"
    summary += f"Log: {log_path}\n"
    summary += f"{'='*70}\n"

    print(summary, flush=True)
    log.info(summary.strip())

    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())
