#!/usr/bin/env python3
"""benchmark_hycal_timing_scan.py — Multi-run EVIO scan for secondary peaks.

Scans all prad_*.evio.* files in an EVIO directory, invokes
benchmark_hycal_timing_test_20261021.py (LMS-corrected, elastic-window
--calibrate) on each split file,
detects any secondary peaks in the energy spectrum (prominence threshold +
min distance from primary), and aggregates σ_E/E trends across files / runs.

Secondary-peak detection uses the **post-calibration inner-ring histogram**
(header ``e_centre  legacy_pre  legacy_post  new_pre  new_post``), not the
bulk/headline spectrum.  This histogram covers only inner-ring events
(r < 200 mm) after per-module gain refinement, which means:

- The inner-ring cut (r < 200 mm) removes most edge/leakage events that
  inflate the bulk spectrum's inelastic shoulder.
- Per-module gain refinement sharpens the elastic peak, making it easier to
  resolve any secondary structure nearby.
- Below-primary inelastic contamination (e.g., the ~1951 MeV shoulder visible
  in full-detector spectra for run 25320) is largely suppressed in this view.
- Any secondary peak that *does* remain visible here — e.g., the ~2600 MeV
  pile-up feature in run 25320 — is a genuine physics contamination (pile-up
  or a real double-population) that is interpretable and actionable.

Secondary-peak detection flags only peaks that are BOTH more than
min_dist_mev away from the primary AND located above the primary peak centre
(higher energy).  Below-primary bumps (inelastic tail, Δ-excitation, Møller)
are expected physics and are not flagged.

By default the full EVIO file is processed (no --max-events cap). Pass an
explicit --max-events value to impose a limit.

Outputs
-------
  output/scan_${RUN}/split_${SPLITNUM}/         per-file benchmark outputs
  output/benchmark_hycal_timing_scan_summary.tsv  aggregate (one row per file)
  output/benchmark_hycal_timing_scan_sigma_E_vs_file_pre.png
      σ_E/E (Gaussian, PRE-calibration) trend plot
  output/benchmark_hycal_timing_scan_sigma_E_vs_file_post.png
      σ_E/E (Gaussian, POST-calibration) trend plot
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
import os
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
    _REPO_ROOT / "analysis" / "pyscripts" / "benchmark_hycal_timing_test_20261021.py"
)
_DATABASE_DIR = Path(os.environ.get("PRAD2_DATABASE_DIR",
                                    str(_REPO_ROOT / "database")))

_DEFAULT_EVIO_DIR  = "/Users/L00338853/work/PRad/data/evio"
_DEFAULT_OUTPUT_DIR = str(_REPO_ROOT / "analysis" / "pyscripts" / "output")
_DEFAULT_MAX_EVENTS: Optional[int] = None   # None → process full EVIO file
_DEFAULT_PROMINENCE = 40
_DEFAULT_PEAK_DIST_MEV = 100.0
_DEFAULT_PARALLEL   = 1
# Pulse-method tagging of outputs: rawsum keeps the legacy names/titles unchanged.
_PULSE_LABEL = {"rawsum": "rawsum", "mf": "mf (template)",
                "rawsum_pc": "rawsum + per-crystal", "mf_pc": "mf + per-crystal"}
_STEM_BASE = "benchmark_hycal_timing_scan"
# --gain-file-map auto: pooled per-crystal gain table per run (percrystal_calib.py fit)
_PERCRYSTAL_DIR = Path(__file__).resolve().parent / "output" / "percrystal"


def _stem(pulse_method: str, spl: bool = False) -> str:
    """Output-file stem; mf gets an ``_mf`` suffix, rawsum is unchanged; the
    per-crystal label keys (``*_pc``) append ``_pc``; spl (--mf-templates-map,
    sampled templates) appends ``_spl``."""
    return (_STEM_BASE + {"rawsum": "", "mf": "_mf", "rawsum_pc": "_pc",
                          "mf_pc": "_mf_pc"}[pulse_method]
            + ("_spl" if spl else ""))


def gain_file_path(run: int, gain_dir: Path = _PERCRYSTAL_DIR) -> Path:
    """Pooled per-crystal gain table of a run (--gain-file-map auto)."""
    return Path(gain_dir) / f"{run:06d}" / f"gains_{run:06d}_all.tsv"


def mf_templates_path(run: int, pattern: str) -> Path:
    """Template file of a run (--mf-templates-map): ``{runtype}`` in pattern ->
    gate_scan_20261022.runtype_slug(run), e.g. '2p24GeV_H2'."""
    import gate_scan_20261022            # lazy: it imports this module
    return Path(pattern.replace("{runtype}", gate_scan_20261022.runtype_slug(run)))


# Runs left out of the scan entirely (user decision); value = reason shown in SKIP line.
_EXCLUDED_RUNS: dict[int, str] = {24578: "non-standard trigger"}

# ---------------------------------------------------------------------------
# Run-number → target peak (MeV) mapping  (spec §2)
# ---------------------------------------------------------------------------
_TARGET_PEAK_MAP: dict[int, float] = {
    24386: 3488.5,
    24503: 3488.5,
    24512: 3488.5,
    24578:  728.9,  # H2, 729 MeV (same as 24590)
    24590:  728.9,
    24603:  728.9,  # H2, 729 MeV (same as 24590)
    25228: 2239.5,  # H2, 2.24 GeV
    25257: 2239.5,  # H2, 2.24 GeV
    25260: 2239.5,  # H2, 2.24 GeV
    25299: 2239.5,  # H2, 2.24 GeV
    25308: 2239.5,
    25320: 2239.5,
    25800: 2239.5,  # X17 (Ta target), 2.24 GeV
    26138: 2239.5,
    26230: 2239.5,  # X17 (Ta target), 2.24 GeV
}
_DEFAULT_TARGET_PEAK = 3488.5  # fallback for unknown runs

# ---------------------------------------------------------------------------
# Run-number → plot colour  (spec §6)
# ---------------------------------------------------------------------------
_RUN_COLORS: dict[int, str] = {
    24386: "blue",
    24503: "green",
    24512: "orange",
    24578: "lightcoral",  # new 729 MeV run (distinct from 24590's "red")
    24590: "red",
    25260: "cyan",      # new H2 run
    24603: "darkred",   # 729 MeV H2
    25228: "navy",      # 2.24 GeV H2
    25257: "goldenrod", # 2.24 GeV H2
    25299: "magenta",   # new H2 run
    25308: "purple",
    25320: "brown",
    25800: "olive",     # new X17 run
    26138: "pink",
    26230: "teal",      # new X17 run
}
_DEFAULT_COLOR = "gray"

# Run type (beam energy + target) → marker shape for the σ_E/E trend plots.
# Targets: H2 unless listed (25308 carbon; X17 runs Ta).
_RUN_TARGET_LABEL: dict[int, str] = {
    25308: "C", 25800: "Ta (X17)", 26138: "Ta (X17)", 26230: "Ta (X17)",
}
_RUN_TYPE_MARKERS: dict[str, str] = {
    "0.73 GeV H2": "o",
    "2.24 GeV H2": "s",
    "2.24 GeV C": "D",
    "2.24 GeV Ta (X17)": "^",
    "3.49 GeV H2": "v",
}
_DEFAULT_MARKER = "X"


def _run_type(run: int) -> str:
    """Run-type label '<E> GeV <target>' from the beam energy and target."""
    e_mev = _TARGET_PEAK_MAP.get(run)
    if e_mev is None:
        return "unknown run type"
    e_gev = e_mev / 1000.0
    return f"{e_gev:.2f} GeV {_RUN_TARGET_LABEL.get(run, 'H2')}"

# ---------------------------------------------------------------------------
# Summary column schema  (spec §5)
# ---------------------------------------------------------------------------
_SUMMARY_COLS = [
    "run_number", "split_number", "filename",
    "n_events", "n_modules_used",
    "sigma_over_E_pre_G", "sigma_over_E_post_G",
    "sigma_over_E_cb_pre", "sigma_over_E_cb_post",   # renamed from sigma_over_E_CB; added _pre
    "primary_peak_mu_pre_MeV", "primary_peak_sigma_pre_MeV",
    "secondary_peak_pre_present", "secondary_peak_mu_pre_MeV", "secondary_peak_sigma_pre_MeV",
    "primary_peak_mu_post_MeV", "primary_peak_sigma_post_MeV",
    "secondary_peak_post_present", "secondary_peak_mu_post_MeV", "secondary_peak_sigma_post_MeV",
    "n_peaks_detected", "secondary_peak_present",
    "peak_prominence_values",
    "sigma_over_E_CB",   # backward-compat alias for sigma_over_E_cb_post
]

# Runinfo pre-window [lo, hi] ns used for each run family
_PRE_WINDOW_MAP: dict[int, tuple[float, float]] = {
    24386: (130.0, 200.0),
    24503: (130.0, 200.0),
    24512: (130.0, 200.0),
    24578: (130.0, 200.0),  # same as 24590
    24590: (130.0, 200.0),
    25260: (100.0, 200.0),  # same as other 25xxx runs
    24603: (130.0, 200.0),  # same as 24590
    25228: (100.0, 200.0),
    25257: (100.0, 200.0),
    25299: (100.0, 200.0),
    25308: (100.0, 200.0),
    25320: (100.0, 200.0),
    25800: (100.0, 200.0),
    26138: (100.0, 200.0),
    26230: (100.0, 200.0),
}
_DEFAULT_PRE_WINDOW = (100.0, 200.0)

# NOTE: _SIGNAL_MIN_MAP and _E_RANGE_MAP are no longer passed to the
# benchmark; benchmark_hycal_timing_test_20261021.py derives --signal-min,
# --e-range, --min-seed-energy, --bin-width and --target-peak from E_beam.
# The tables are kept for reference.
# signal-min cut (MeV) per run — the lowest single-cluster energy accepted
# before histogramming. For low-beam-energy runs (e.g. 24590 at 729 MeV),
# a 500 MeV cut would eat into the elastic peak, so a smaller value is used.
_SIGNAL_MIN_MAP: dict[int, float] = {
    24386: 500.0,
    24503: 500.0,
    24512: 500.0,
    24578: 200.0,   # 729 MeV beam — need lower cut to keep elastic peak
    24590: 200.0,   # 729 MeV beam — need lower cut to keep the elastic peak
    24603: 200.0,   # same as 24590
    25228: 500.0,   # same as 25260
    25257: 500.0,   # same as 25260
    25260: 500.0,
    25299: 500.0,
    25308: 500.0,
    25320: 500.0,
    25800: 500.0,
    26138: 500.0,
    26230: 500.0,
}
_DEFAULT_SIGNAL_MIN = 500.0

# Energy histogram range for find_peaks input
_E_RANGE_MAP: dict[int, tuple[float, float]] = {
    24386: (500.0, 4500.0),
    24503: (500.0, 4500.0),
    24512: (500.0, 4500.0),
    24578: (500.0, 2000.0),  # same as 24590 (low-energy run)
    24590: (500.0, 2000.0),
    24603: (500.0, 2000.0),  # same as 24590
    25308: (500.0, 3500.0),
    25320: (500.0, 3500.0),
    26138: (500.0, 3500.0),
}
_DEFAULT_E_RANGE = (500.0, 4500.0)

# HyCal calibration file per run (FADC-integral calibrations), relative to
# $PRAD2_DATABASE_DIR/calibration.  X17 runs use the _fadc variant of their
# runinfo period file.  Runs not listed fall back to the runinfo default,
# which the benchmark checks against E_beam (exit 2 on mismatch).
_HC_CALIB_MAP: dict[int, str] = {
    **{r: "calibration_factor_fadc_3p5.json" for r in (24386, 24503, 24512)},
    **{r: "calibration_factor_fadc_0p7.json" for r in (24578, 24590, 24603)},
    **{r: "calibration_factor_2p2_fadc.json"
       for r in (25228, 25257, 25260, 25299, 25308, 25320)},
    25800: "x17_25825_fadc.json",
    26138: "x17_26094_fadc.json",
    26230: "x17_26214_fadc.json",
}


def gain_corr_path(run: int) -> Path:
    """LMS gain_corr file the benchmark's --lms-corr auto would use."""
    return (_DATABASE_DIR / "gain_factor" / "gain_correction"
            / f"prad_{run:06d}_gain_corr.root")


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


def filter_excluded_runs(run_map: dict[int, list[tuple[int, Path]]]
                         ) -> dict[int, list[tuple[int, Path]]]:
    """Drop runs listed in _EXCLUDED_RUNS, printing one SKIP line per dropped run."""
    kept: dict[int, list[tuple[int, Path]]] = {}
    for run, splits in run_map.items():
        if run in _EXCLUDED_RUNS:
            print(f"SKIP run {run:06d}: {_EXCLUDED_RUNS[run]}", flush=True)
            continue
        kept[run] = splits
    return kept


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
    timeout_seconds: Optional[int] = None,
    pulse_method: str = "rawsum",
    window: Optional[float] = None,
    pre_window: Optional[tuple[float, float]] = None,
    gain_file: Optional[Path] = None,
    write_cluster_tsv: bool = False,
    mf_templates: Optional[Path] = None,
    mf_scale: Optional[float] = None,
) -> Path:
    """Invoke benchmark_hycal_timing_test_20261021.py; return path to the TSV.

    gain_file, when given, is passed as ``--gain-file`` (per-crystal gains at
    the member level, the seed-indexed elastic-window refinement on top); a
    missing file raises RuntimeError (the split is skipped and logged).
    write_cluster_tsv appends ``--write-cluster-tsv`` (result_clusters.tsv).
    mf_templates, when given, is passed as ``--mf-templates`` (a missing file
    raises RuntimeError like gain_file); mf_scale as ``--mf-scale``.

    pulse_method "mf" appends ``--pulse-method mf``; window / pre_window, when
    given, override the benchmark default window and the per-run pre-window.

    If skip_benchmark is True and the TSV already exists, skip the invocation.
    If max_events is None, --max-events is NOT passed to the subprocess so the
    full EVIO file is processed.
    If timeout_seconds is a positive integer, the subprocess is killed after
    that many wall-clock seconds; a RuntimeError is raised so the caller can
    record the file as a timeout failure and continue.
    Raises RuntimeError on subprocess failure.
    """
    out_d   = benchmark_out_dir(output_dir, run, split)
    tsv_p   = benchmark_tsv_path(output_dir, run, split)
    log_p   = out_d / "benchmark.log"
    out_stem = str(out_d / "result")

    if skip_benchmark and tsv_p.exists():
        log.info("skip (TSV exists): %s", tsv_p)
        return tsv_p

    if gain_file is not None and not Path(gain_file).is_file():
        raise RuntimeError(f"gain file missing: {gain_file}")
    if mf_templates is not None and not Path(mf_templates).is_file():
        raise RuntimeError(f"mf templates missing: {mf_templates}")
    out_d.mkdir(parents=True, exist_ok=True)

    pre_lo, pre_hi = (pre_window if pre_window is not None
                      else _PRE_WINDOW_MAP.get(run, _DEFAULT_PRE_WINDOW))

    cmd = [
        sys.executable,
        _BENCHMARK_SCRIPT,
        str(evio_file),
        out_stem,
        "--calibrate",
        "--exclude-rowcol", "15", "20", "15", "20",
        "--pre-window", str(pre_lo), str(pre_hi),
        "--inner-radius", "200",
        "--min-per-module", "30",
        # elastic-window (default); target from the benchmark's run map
        "--lms-corr", "auto" if gain_corr_path(run).exists() else "off",
    ]
    if run in _HC_CALIB_MAP:
        cmd.extend(["--hc-calib-file",
                    str(_DATABASE_DIR / "calibration" / _HC_CALIB_MAP[run])])
    # Only pass --max-events when an explicit limit is requested.
    # When max_events is None the subprocess processes the full EVIO file.
    if max_events is not None:
        cmd.extend(["--max-events", str(max_events)])
    if window is not None:
        cmd.extend(["--window", str(window)])
    if pulse_method == "mf":
        cmd.extend(["--pulse-method", "mf"])
    if mf_templates is not None:
        cmd.extend(["--mf-templates", str(mf_templates)])
    if mf_scale is not None:
        cmd.extend(["--mf-scale", f"{mf_scale:g}"])
    if gain_file is not None:
        cmd.extend(["--gain-file", str(gain_file)])
    if write_cluster_tsv:
        cmd.append("--write-cluster-tsv")

    max_ev_label = str(max_events) if max_events is not None else "full file"
    log.info("run %06d split %05d → %s  (max_events=%s)",
             run, split, out_stem, max_ev_label)
    print(f"Processing run {run} split {split:05d}  (max_events={max_ev_label})...",
          flush=True)

    t0 = time.monotonic()
    try:
        with open(log_p, "w") as lf:
            timeout_arg = timeout_seconds if (timeout_seconds and timeout_seconds > 0) else None
            try:
                result = subprocess.run(
                    cmd,
                    stdout=lf,
                    stderr=subprocess.STDOUT,
                    check=False,
                    timeout=timeout_arg,
                )
            except subprocess.TimeoutExpired as te:
                elapsed = time.monotonic() - t0
                msg = (
                    f"benchmark TIMEOUT after {elapsed:.1f}s "
                    f"(limit: {timeout_arg}s) run={run} split={split}. "
                    f"File likely triggers the known hang bug. Log: {log_p}"
                )
                log.error(msg)
                # Also log to the per-file log so there's a record
                with open(log_p, "a") as lf2:
                    lf2.write(f"\n\n[ORCHESTRATOR] TIMEOUT after {elapsed:.1f}s (limit {timeout_arg}s)\n")
                raise RuntimeError(msg)
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
      'headline'          : {'legacy': {...}, 'new': {...}}  (bulk path metrics)
      'calibration'       : {'legacy': {...}, 'new': {...}}  (inner-ring cal metrics)
      'histogram'         : {'centres': [...], 'legacy': [...], 'new': [...]}
                            Bulk/headline histogram (all single-cluster events,
                            full detector).  Kept for reference; NOT used for
                            peak detection.
      'postcal_histogram' : {'centres': [...], 'legacy_pre': [...],
                             'legacy_post': [...], 'new_pre': [...],
                             'new_post': [...]}
                            Inner-ring post-calibration histogram (r < 200 mm,
                            after per-module gain refinement).  Used for
                            secondary peak detection.

    Raises RuntimeError loudly if headline block is missing (spec §10).
    calibration block absence only raises if --calibrate was expected.
    """
    text  = tsv_path.read_text(encoding="utf-8")
    lines = text.splitlines()

    headline_hdr_lineno    = None
    calibration_hdr_lineno = None
    hist_hdr_lineno        = None
    postcal_hdr_lineno     = None

    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith("path\tn_events\t") and headline_hdr_lineno is None:
            headline_hdr_lineno = i
        if s.startswith("# --calibrate inner ring"):
            calibration_hdr_lineno = i
        if s.startswith("e_centre\tcount_legacy\t") and hist_hdr_lineno is None:
            hist_hdr_lineno = i
        if s.startswith("e_centre\tlegacy_pre\t") and postcal_hdr_lineno is None:
            postcal_hdr_lineno = i

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

    # Parse bulk histogram block (kept for reference; not used for peak detection)
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

    # Parse post-calibration inner-ring histogram block
    # Header: e_centre  legacy_pre  legacy_post  new_pre  new_post
    # Only present when --calibrate is active (which the scan always uses).
    postcal: dict = {}
    if postcal_hdr_lineno is not None:
        pc_centres, pc_leg_pre, pc_leg_post, pc_new_pre, pc_new_post = (
            [], [], [], [], []
        )
        for line in lines[postcal_hdr_lineno + 1:]:
            s = line.strip()
            if not s or s.startswith("#"):
                break
            parts = s.split("\t")
            if len(parts) < 5:
                continue
            try:
                pc_centres.append(float(parts[0]))
                pc_leg_pre.append(int(parts[1]))
                pc_leg_post.append(int(parts[2]))
                pc_new_pre.append(int(parts[3]))
                pc_new_post.append(int(parts[4]))
            except ValueError:
                continue
        if pc_centres:
            postcal = {
                "centres":    np.array(pc_centres,  dtype=float),
                "legacy_pre": np.array(pc_leg_pre,  dtype=int),
                "legacy_post":np.array(pc_leg_post, dtype=int),
                "new_pre":    np.array(pc_new_pre,  dtype=int),
                "new_post":   np.array(pc_new_post, dtype=int),
            }
    if not postcal:
        log.warning(
            "WARNING: post-calibration histogram missing in %s; "
            "peak detection skipped",
            tsv_path,
        )

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
        "headline":          headline,
        "calibration":       calibration,
        "histogram":         hist,       # bulk, kept for reference
        "postcal_histogram": postcal,    # inner-ring post-cal, used for peak detection
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
    postcal_histogram: dict,
    target_peak: float,
    prominence_threshold: float,
    min_dist_mev: float,
    min_prominence_ratio: float = 0.10,  # NEW: minimum secondary/primary prominence ratio
    counts_key: str = "legacy_post",
    label: str = "post-cal",
) -> dict:
    """Detect primary and secondary peaks in a post-calibration inner-ring spectrum.

    Uses the inner-ring post-calibration histogram (r < 200 mm, per-module
    gain-refined legacy counts) rather than the bulk/headline spectrum.  This
    view suppresses most inelastic contamination and sharpens the elastic peak,
    making secondary peaks more interpretable:

    - The inner-ring cut removes most edge/leakage events.
    - Per-module gain refinement sharpens the elastic peak.
    - Any secondary peak visible here (e.g. above-primary pile-up) is a
      genuine contamination feature rather than an artefact of the bulk
      detector acceptance.

    A secondary peak is identified only when it is BOTH more than
    min_dist_mev away from the primary peak AND located above the primary
    peak centre (higher energy).  Below-primary bumps (inelastic tail) are
    expected physics and are not flagged.

    Parameters
    ----------
    postcal_histogram:
        Dict with ``centres`` array and per-histogram count arrays (e.g.
        ``legacy_pre``, ``legacy_post``).
    target_peak:
        Expected energy (MeV) of the primary (elastic) peak.
    prominence_threshold:
        scipy find_peaks prominence threshold.
    min_dist_mev:
        Minimum MeV separation between primary and any secondary peak.
    min_prominence_ratio:
        Minimum ratio (secondary_prominence / primary_prominence) for a
        detected peak to qualify as a real secondary. Default 0.10 (secondary
        must be at least 10% as prominent as the primary). Rejects false
        positives from noise where scipy finds a weak bump that satisfies the
        absolute prominence threshold but is tiny relative to the elastic peak.
    counts_key:
        Key in ``postcal_histogram`` to use for counts (default: ``legacy_post``).
    label:
        Short label included in log/print messages (default: ``post-cal``).

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

    if not HAVE_SCIPY or not postcal_histogram or "centres" not in postcal_histogram:
        return result
    if postcal_histogram["centres"].size == 0:
        return result
    if counts_key not in postcal_histogram:
        log.warning("  [%s] counts_key '%s' not found in postcal_histogram", label, counts_key)
        return result

    centres = postcal_histogram["centres"]
    counts  = postcal_histogram[counts_key].astype(float)

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

    # Always report the primary peak position.
    primary_e = centres[primary_bin]
    # Use the Gaussian-fitted centre when available; fall back to histogram bin centre.
    primary_peak_center = primary_mu if not math.isnan(primary_mu) else primary_e
    log.info("  [%s] Primary peak at %.0f MeV", label, primary_peak_center)
    print(f"  [{label}] Primary peak at {primary_peak_center:.0f} MeV", flush=True)

    # Secondary peak: a peak that is BOTH:
    #   (a) more than min_dist_mev away from the primary peak centre, AND
    #   (b) located ABOVE the primary peak centre (i.e. higher energy), AND
    #   (c) at least min_prominence_ratio × primary prominence (ratio check).
    # Peaks below the primary are expected physics (inelastic / Δ-excitation /
    # Møller) in the post-cal inner-ring histogram's low-energy tail and are
    # NOT flagged. Only above-primary peaks (pile-up signatures) are flagged.

    # Primary peak prominence (for ratio check)
    primary_prominence = float(props["prominences"][primary_idx_in_peaks])

    secondary_bins = []
    for pk in peaks:
        if pk == primary_bin:
            continue  # skip the primary itself
        # Spatial criteria
        if abs(centres[pk] - primary_e) <= min_dist_mev:
            continue
        if centres[pk] <= primary_peak_center:
            continue  # must be above the primary
        # Prominence ratio criterion (NEW)
        pk_idx_in_peaks = int(np.where(peaks == pk)[0][0])
        pk_prominence = float(props["prominences"][pk_idx_in_peaks])
        if primary_prominence > 0 and pk_prominence / primary_prominence < min_prominence_ratio:
            continue  # too small relative to primary
        secondary_bins.append(pk)

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
            "  [%s] Found above-primary secondary peak at %.0f MeV "
            "(σ=%.0f MeV, prominence=%.1f, primary=%.0f MeV)",
            label,
            sec_mu if not math.isnan(sec_mu) else centres[sec_bin],
            sec_sigma if not math.isnan(sec_sigma) else float("nan"),
            secondary_proms[best_sec_idx],
            primary_peak_center,
        )
        print(
            f"  [{label}] Found above-primary secondary peak at "
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
    min_prominence_ratio: float = 0.10,  # NEW
) -> dict:
    """Parse TSV and peak-detect; return one aggregate row dict."""
    nan = float("nan")
    row: dict = {
        "run_number":                    run,
        "split_number":                  split,
        "filename":                      evio_file.name,
        "n_events":                      nan,
        "n_modules_used":                nan,
        "sigma_over_E_pre_G":            nan,
        "sigma_over_E_post_G":           nan,
        "sigma_over_E_cb_pre":           nan,
        "sigma_over_E_cb_post":          nan,
        "sigma_over_E_CB":               nan,   # backward-compat alias
        # Pre-calibration peak columns
        "primary_peak_mu_pre_MeV":       nan,
        "primary_peak_sigma_pre_MeV":    nan,
        "secondary_peak_pre_present":    "N",
        "secondary_peak_mu_pre_MeV":     nan,
        "secondary_peak_sigma_pre_MeV":  nan,
        # Post-calibration peak columns
        "primary_peak_mu_post_MeV":      nan,
        "primary_peak_sigma_post_MeV":   nan,
        "secondary_peak_post_present":   "N",
        "secondary_peak_mu_post_MeV":    nan,
        "secondary_peak_sigma_post_MeV": nan,
        # Aggregate peak info (post-cal aliased for backward compatibility)
        "n_peaks_detected":              0,
        "secondary_peak_present":        "N",
        "peak_prominence_values":        "[]",
    }

    parsed = parse_tsv(tsv_path)
    hl  = parsed["headline"].get("legacy", {})
    cal = parsed["calibration"].get("legacy", {})

    row["n_events"]             = hl.get("n_events", nan)
    row["n_modules_used"]       = cal.get("n_modules", nan)
    row["sigma_over_E_cb_post"] = cal.get("sigma_over_E_cb", nan)
    row["sigma_over_E_cb_pre"]  = cal.get("sigma_over_E_cb_pre", nan)
    row["sigma_over_E_CB"]      = row["sigma_over_E_cb_post"]   # backward-compat alias

    # Pre-calibration Gaussian σ_E/E = sigma_pre_MeV / mu_pre_MeV
    mu_pre    = cal.get("mu_pre_MeV",    nan)
    sigma_pre = cal.get("sigma_pre_MeV", nan)
    if (not math.isnan(mu_pre) and not math.isnan(sigma_pre)
            and mu_pre != 0.0):
        row["sigma_over_E_pre_G"] = sigma_pre / mu_pre
    else:
        row["sigma_over_E_pre_G"] = nan

    # Post-calibration Gaussian σ_E/E — use sigma_over_E_post from TSV directly;
    # fall back to sigma_post_MeV / mu_post_MeV if the direct field is absent.
    sigma_over_E_post = cal.get("sigma_over_E_post", nan)
    if not math.isnan(sigma_over_E_post):
        row["sigma_over_E_post_G"] = sigma_over_E_post
    else:
        mu_post    = cal.get("mu_post_MeV",    nan)
        sigma_post = cal.get("sigma_post_MeV", nan)
        if (not math.isnan(mu_post) and not math.isnan(sigma_post)
                and mu_post != 0.0):
            row["sigma_over_E_post_G"] = sigma_post / mu_post
        else:
            row["sigma_over_E_post_G"] = nan

    # Peak detection on the inner-ring histogram — run on BOTH pre and post columns.
    target_peak = _TARGET_PEAK_MAP.get(run, _DEFAULT_TARGET_PEAK)

    # Pre-calibration inner-ring histogram
    pre_peaks = detect_peaks(
        postcal_histogram=parsed["postcal_histogram"],
        target_peak=target_peak,
        prominence_threshold=prominence_threshold,
        min_dist_mev=min_dist_mev,
        min_prominence_ratio=min_prominence_ratio,
        counts_key="legacy_pre",
        label="pre-cal",
    )

    # Post-calibration inner-ring histogram
    post_peaks = detect_peaks(
        postcal_histogram=parsed["postcal_histogram"],
        target_peak=target_peak,
        prominence_threshold=prominence_threshold,
        min_dist_mev=min_dist_mev,
        min_prominence_ratio=min_prominence_ratio,
        counts_key="legacy_post",
        label="post-cal",
    )

    # Populate pre-calibration peak columns
    row["primary_peak_mu_pre_MeV"]       = pre_peaks["primary_mu"]
    row["primary_peak_sigma_pre_MeV"]    = pre_peaks["primary_sigma"]
    row["secondary_peak_pre_present"]    = "Y" if pre_peaks["secondary_present"] else "N"
    row["secondary_peak_mu_pre_MeV"]     = pre_peaks["secondary_mu"]
    row["secondary_peak_sigma_pre_MeV"]  = pre_peaks["secondary_sigma"]

    # Populate post-calibration peak columns
    row["primary_peak_mu_post_MeV"]      = post_peaks["primary_mu"]
    row["primary_peak_sigma_post_MeV"]   = post_peaks["primary_sigma"]
    row["secondary_peak_post_present"]   = "Y" if post_peaks["secondary_present"] else "N"
    row["secondary_peak_mu_post_MeV"]    = post_peaks["secondary_mu"]
    row["secondary_peak_sigma_post_MeV"] = post_peaks["secondary_sigma"]

    # Aggregate / backward-compat columns (keyed off post-cal)
    row["n_peaks_detected"]    = post_peaks["n_peaks"]
    row["secondary_peak_present"] = "Y" if post_peaks["secondary_present"] else "N"
    row["peak_prominence_values"] = json.dumps(
        [round(p, 2) for p in post_peaks["prominences"]]
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
        if col in ("sigma_over_E_pre_G", "sigma_over_E_post_G",
                   "sigma_over_E_cb_pre", "sigma_over_E_cb_post", "sigma_over_E_CB"):
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


def make_trend_plot(
    rows: list[dict],
    output_path: Path,
    value_col_g: str = "sigma_over_E_post_G",
    value_col_cb: str = "sigma_over_E_cb_post",
    title_suffix: str = "POST-calibration",
    secondary_flag_col: str = "secondary_peak_post_present",
    # Legacy positional-keyword compat shim: ignored if callers pass value_col_g.
    pulse_method: str = "rawsum",
) -> None:
    """σ_E/E vs. EVIO file, overlaying Gaussian and Crystal Ball fits.

    Plots TWO lines per run on the same axes:
    - Solid line + open markers: Gaussian σ_E/E (``value_col_g``)
    - Dashed line + filled markers: Crystal Ball σ_E/E (``value_col_cb``)
    Marker shape encodes the run type (beam energy + target).

    Both lines for a given run share the same colour.  Star markers (crimson)
    indicate files where a secondary peak was detected; they are placed on the
    CB values to reduce clutter.  NaN values are skipped naturally by
    matplotlib (line gaps).

    Parameters
    ----------
    rows:
        Aggregate row dicts (one per EVIO file), already sorted by
        (run_number, split_number).
    output_path:
        Destination PNG file path.
    value_col_g:
        Column name for the Gaussian σ_E/E (a raw ratio).  Multiplied by 100
        for display.
    value_col_cb:
        Column name for the Crystal Ball σ_E/E (a raw ratio).  Multiplied by
        100 for display.  Will be NaN for TSVs produced before the CB-pre fit
        was added — those files simply show a gap in the CB line.
    title_suffix:
        String appended to the plot title, e.g. ``"PRE-calibration"`` or
        ``"POST-calibration"``.
    secondary_flag_col:
        Column used to decide where to draw star markers.
        Pre-cal: ``"secondary_peak_pre_present"``; post-cal (default):
        ``"secondary_peak_post_present"``.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(max(12, len(rows) * 0.18 + 4), 6))

    # Group rows by run for run-level legends
    run_groups: dict[int, list[dict]] = {}
    for row in rows:
        run_groups.setdefault(int(row["run_number"]), []).append(row)

    # Draw per-run lines + per-point markers (two series per run)
    all_x_labels: list[str] = []
    all_x_positions: list[int] = []
    x_pos = 0

    # Stars placed on CB values (only one star per file to reduce clutter)
    sec_x: list[int]   = []
    sec_y: list[float] = []

    # Collect all valid y-values to determine y-axis range
    all_y_values: list[float] = []

    types_seen: list[str] = []
    for run in sorted(run_groups.keys()):
        color    = _RUN_COLORS.get(run, _DEFAULT_COLOR)
        rtype    = _run_type(run)
        marker   = _RUN_TYPE_MARKERS.get(rtype, _DEFAULT_MARKER)
        if rtype not in types_seen:
            types_seen.append(rtype)
        run_rows = run_groups[run]

        run_x_g:  list[int]   = []
        run_y_g:  list[float] = []
        run_x_cb: list[int]   = []
        run_y_cb: list[float] = []

        for row in run_rows:
            lbl = _label(int(row["run_number"]), int(row["split_number"]))
            all_x_labels.append(lbl)
            all_x_positions.append(x_pos)

            # Gaussian value
            sigma_g = row.get(value_col_g, float("nan"))
            if isinstance(sigma_g, float) and not math.isnan(sigma_g):
                sigma_g_pct = sigma_g * 100.0
                run_x_g.append(x_pos)
                run_y_g.append(sigma_g_pct)
                all_y_values.append(sigma_g_pct)

            # Crystal Ball value
            sigma_cb = row.get(value_col_cb, float("nan"))
            if isinstance(sigma_cb, float) and not math.isnan(sigma_cb):
                sigma_cb_pct = sigma_cb * 100.0
                run_x_cb.append(x_pos)
                run_y_cb.append(sigma_cb_pct)
                all_y_values.append(sigma_cb_pct)
                # Star on CB line when secondary peak present
                if row.get(secondary_flag_col) == "Y":
                    sec_x.append(x_pos)
                    sec_y.append(sigma_cb_pct)
            elif row.get(secondary_flag_col) == "Y" and run_x_g:
                # CB value absent — fall back to placing star on the G value
                last_g_idx = run_x_g.index(x_pos) if x_pos in run_x_g else -1
                if last_g_idx >= 0:
                    sec_x.append(x_pos)
                    sec_y.append(run_y_g[last_g_idx])

            x_pos += 1

        # Solid line + open markers: Gaussian
        if run_x_g:
            ax.plot(run_x_g, run_y_g, "-", marker=marker, color=color, lw=1.4,
                    ms=6, mfc="none", mec=color, zorder=3)

        # Dashed line + filled markers: Crystal Ball
        if run_x_cb:
            ax.plot(run_x_cb, run_y_cb, "--", marker=marker, color=color,
                    lw=1.2, ms=5, zorder=3)

        # One colour entry per run (shape/fill explained separately below)
        ax.plot([], [], "-", color=color, lw=2, label=f"run {run}")

    # Overlay star markers for files with secondary peaks (on CB when available)
    if sec_x:
        ax.scatter(sec_x, sec_y, marker="*", s=180, color="crimson",
                   zorder=5, label="secondary peak detected")

    # Y-axis range: accommodate both G and CB values
    if all_y_values:
        max_y = max(all_y_values)
        ax.set_ylim(0, max_y * 1.1)

    ax.set_xticks(all_x_positions)
    ax.set_xticklabels(all_x_labels, rotation=45, ha="right", fontsize=7)
    ax.set_xlabel("EVIO file")
    ax.set_ylabel("σ_E/E (%)")
    if pulse_method != "rawsum":
        title_suffix = f"{title_suffix}, {_PULSE_LABEL[pulse_method]}"
    ax.set_title(
        f"σ_E/E ({title_suffix}) across EVIO files — all runs\n"
        "(open marker + solid = Gaussian, filled marker + dashed = CB, "
        "marker shape = run type, stars = secondary peak detected)"
    )
    # Proxy handles: run-type shapes, then the fill convention
    from matplotlib.lines import Line2D
    run_handles, run_labels = ax.get_legend_handles_labels()
    types_seen = ([t for t in _RUN_TYPE_MARKERS if t in types_seen]
                  + [t for t in types_seen if t not in _RUN_TYPE_MARKERS])
    type_handles = [Line2D([], [], ls="none", marker=_RUN_TYPE_MARKERS.get(t, _DEFAULT_MARKER),
                           color="k", mfc="k", ms=6, label=t) for t in types_seen]
    fill_handles = [
        Line2D([], [], ls="--", marker="o", color="k", mfc="k", ms=6,
               label="filled + dashed: CB"),
        Line2D([], [], ls="-", marker="o", color="k", mfc="none", ms=6,
               label="open + solid: Gaussian"),
    ]
    if pulse_method != "rawsum":
        fill_handles.append(Line2D([], [], ls="none", label=
                                   f"pulse method: {_PULSE_LABEL[pulse_method]}"))
    ax.legend(handles=run_handles + type_handles + fill_handles,
              loc="upper right", fontsize=7, ncol=3)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


# ===========================================================================
# §7 — Summary markdown
# ===========================================================================

def write_findings_md(rows: list[dict], out_path: Path,
                      pulse_method: str = "rawsum",
                      window: Optional[float] = None,
                      pre_window: Optional[tuple[float, float]] = None,
                      stem: str = _STEM_BASE) -> None:
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
    if pulse_method != "rawsum" or window is not None or pre_window is not None:
        wp = [f"pulse method: **{_PULSE_LABEL[pulse_method]}**"]
        if window is not None:
            wp.append(f"window = {window:g} ns")
        if pre_window is not None:
            wp.append(f"pre-window = {pre_window[0]:g}\u2013{pre_window[1]:g} ns")
        lines.append("_" + "; ".join(wp) + "_\n")
    lines.append("")
    lines.append("## Trend Plots\n")
    lines.append(
        "Two σ_E/E trend plots are generated (one row per EVIO split file,\n"
        "colour by run, star markers for files with secondary peaks detected).\n"
        "Each plot overlays **both** the Gaussian fit (solid line, open markers)\n"
        "and the Crystal Ball fit (dashed line, filled markers) on the same axes;\n"
        "marker shape = run type (beam energy + target).\n"
        "CB values will be absent (line gap) for TSVs produced before the CB-pre\n"
        "fit was added; they will populate once the benchmark is re-run.\n"
    )
    lines.append(
        "- **Pre-calibration (Gaussian + CB)**: "
        f"`{stem}_sigma_E_vs_file_pre.png`  \n"
        "  Gaussian: σ_E/E = σ_pre / μ_pre from the raw inner-ring Gaussian fit.  \n"
        "  Crystal Ball: σ_E/E = σ_cb_pre from the CB fit to the pre-calibration spectrum."
    )
    lines.append(
        "- **Post-calibration (Gaussian + CB)**: "
        f"`{stem}_sigma_E_vs_file_post.png`  \n"
        "  Gaussian: σ_E/E = σ_over_E_post from the inner-ring Gaussian fit after "
        "per-module gain refinement.  \n"
        "  Crystal Ball: σ_E/E = σ_cb_post from the CB fit to the post-calibration spectrum."
    )
    lines.append("")
    lines.append(
        "> **Note**: Secondary peaks are flagged only when a distinct peak is "
        "detected ABOVE the primary (elastic) peak centre in the post-calibration "
        "inner-ring histogram, with prominence > threshold and at least "
        "`min_dist_mev` away from the primary.  Peaks below the primary "
        "(inelastic tail, Δ-excitation, Møller) are expected physics and are "
        "not flagged.  Peak detection is run on **both** the pre-calibration "
        "(`legacy_pre`) and post-calibration (`legacy_post`) inner-ring "
        "histograms; both primary and secondary peak energies are recorded in "
        "the summary TSV for comparison.\n"
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
    pulse_method: str = "rawsum",
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
            sec_mu = row.get("secondary_peak_mu_post_MeV", float("nan"))
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
            if pulse_method != "rawsum":
                title = f"[{_PULSE_LABEL[pulse_method]}]  " + title

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
        "--min-prominence-ratio",
        type=float,
        default=0.10,
        help=f"Minimum ratio (secondary_prominence / primary_prominence) for a "
             f"secondary peak to be flagged. Default: 0.10 (secondary must be "
             f"at least 10%% as prominent as the primary). Rejects weak noise "
             f"bumps that pass the absolute prominence threshold.",
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
    ap.add_argument(
        "--timeout-seconds",
        type=int,
        default=1200,  # 20 minutes
        help="Maximum wall-clock seconds for each benchmark subprocess. "
             "Workers that exceed this are killed, their EVIO file is logged "
             "as a timeout, and the scan continues. Set to 0 to disable timeout. "
             "Default: 1200 (20 minutes).",
    )
    ap.add_argument(
        "--pulse-method", choices=["rawsum", "mf"], default="rawsum",
        help="Pulse reconstruction passed to the benchmark (default rawsum; "
             "mf appends '--pulse-method mf' and tags output names/titles "
             "with 'mf').",
    )
    ap.add_argument(
        "--window", type=float, default=None, metavar="W",
        help="seed_time_window (ns) passed to the benchmark for every run. "
             "Default: not passed (benchmark default).",
    )
    ap.add_argument(
        "--pre-window", type=float, nargs=2, default=None, metavar=("LO", "HI"),
        help="Pre-window [LO, HI] ns for every run, overriding the per-run "
             "map. Default: per-run map.",
    )
    ap.add_argument(
        "--gain-file-map", choices=["off", "auto"], default="off",
        help="auto: pass --gain-file output/percrystal/<run>/gains_<run>_all.tsv "
             "(pooled per-crystal gains) for every run, with the seed-indexed "
             "elastic-window refinement on top; output names/titles get the "
             "'per-crystal' tag.  The tables are fitted on mf dumps with "
             "W = 16 ns and pre-window 100-200 ns: use matching flags.  "
             "Default off (commands unchanged).",
    )
    ap.add_argument(
        "--gain-dir", type=Path, default=_PERCRYSTAL_DIR, metavar="DIR",
        help="Directory of the --gain-file-map auto tables "
             "(<DIR>/<run>/gains_<run>_all.tsv).  Default output/percrystal.",
    )
    ap.add_argument(
        "--mf-templates-map", default=None, metavar="PATTERN",
        help="mf only: per-run template file, passed as --mf-templates; "
             "{runtype} in PATTERN -> run-type slug (2p24GeV_H2, 2p24GeV_C, "
             "3p49GeV_H2, 0p73GeV_H2, 2p24GeV_Ta_X17), e.g. "
             "output/templates_spline/{runtype}/templates_{runtype}.npz.  A "
             "missing file skips the split like a missing gain file.  Output "
             "names get the '_spl' stem label.  Default off (commands "
             "unchanged).",
    )
    ap.add_argument(
        "--mf-scale", type=float, default=None, metavar="F",
        help="mf only: pass --mf-scale F to the benchmark (energy scale kappa).  "
             "Default: not passed (benchmark default).",
    )
    ap.add_argument(
        "--write-cluster-tsv", action="store_true",
        help="Pass --write-cluster-tsv to the benchmark (per-event "
             "result_clusters.tsv per split, input of "
             "plot_module_resolution.py).  Default off (commands unchanged).",
    )
    return ap


def main(argv=None) -> int:
    ap    = build_parser()
    args  = ap.parse_args(argv)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- Logging setup -------------------------------------------------------
    # label key for names/titles; the benchmark still gets args.pulse_method
    lab = args.pulse_method + ("_pc" if args.gain_file_map == "auto" else "")
    if args.pulse_method != "mf" and (args.mf_templates_map is not None
                                      or args.mf_scale is not None):
        ap.error("--mf-templates-map and --mf-scale need --pulse-method mf")
    stem = _stem(lab, spl=args.mf_templates_map is not None)
    pre_window = tuple(args.pre_window) if args.pre_window is not None else None
    log_path = output_dir / f"{stem}.log"
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
    log.info("  min_prom_ratio: %.2f", args.min_prominence_ratio)
    log.info("  parallel    : %d", args.parallel)
    log.info("  skip_bench  : %s", args.skip_benchmark)
    log.info("  timeout     : %s",
             f"{args.timeout_seconds}s" if args.timeout_seconds > 0 else "disabled")
    if lab != "rawsum" or args.window is not None or pre_window is not None:
        log.info("  pulse_method: %s", _PULSE_LABEL[lab])
        log.info("  window      : %s",
                 f"{args.window:g} ns" if args.window is not None else "benchmark default")
        log.info("  pre_window  : %s",
                 f"{pre_window[0]:g}-{pre_window[1]:g} ns" if pre_window is not None
                 else "per-run map")
    if args.mf_templates_map is not None:
        log.info("  mf templates: %s", args.mf_templates_map)
    if args.mf_scale is not None:
        log.info("  mf scale    : %g", args.mf_scale)
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
    run_map = filter_excluded_runs(run_map)   # drop excluded runs before scheduling

    total_files = sum(len(splits) for splits in run_map.values())
    log.info(
        "Discovered %d run(s), %d total files: %s",
        len(run_map), total_files,
        {r: len(s) for r, s in sorted(run_map.items())},
    )

    no_lms = sorted(r for r in run_map if not gain_corr_path(r).exists())
    if no_lms:
        log.warning("No gain_corr file, running with --lms-corr off: %s", no_lms)

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
                timeout_seconds=args.timeout_seconds,   # NEW
                pulse_method=args.pulse_method,
                window=args.window,
                pre_window=pre_window,
                gain_file=(gain_file_path(run, args.gain_dir)
                           if args.gain_file_map == "auto"
                           else None),
                write_cluster_tsv=args.write_cluster_tsv,
                mf_templates=(mf_templates_path(run, args.mf_templates_map)
                              if args.mf_templates_map is not None else None),
                mf_scale=args.mf_scale,
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
                min_prominence_ratio=args.min_prominence_ratio,  # NEW
            )
            all_rows.append(row)
            n_ev = row.get("n_events", float("nan"))
            n_ev_str = (
                f"{int(n_ev):,}" if not math.isnan(float(n_ev)) else "unknown"
            )
            # Format primary peak energies for both pre and post
            mu_pre_f  = float(row.get("primary_peak_mu_pre_MeV",  float("nan")))
            mu_post_f = float(row.get("primary_peak_mu_post_MeV", float("nan")))
            mu_pre_str  = f"{mu_pre_f:.0f} MeV"  if not math.isnan(mu_pre_f)  else "N/A"
            mu_post_str = f"{mu_post_f:.0f} MeV" if not math.isnan(mu_post_f) else "N/A"
            sec_pre  = row.get("secondary_peak_pre_present",  "N")
            sec_post = row.get("secondary_peak_post_present", "N")
            sigma_cb = float(row["sigma_over_E_CB"])
            sigma_str = f"{sigma_cb * 100.0:.3f}%" if not math.isnan(sigma_cb) else "N/A"
            log.info(
                "  run=%d split=%05d  events=%s"
                "  pre: μ=%s sec=%s  post: μ=%s sec=%s  σ_E/E=%s",
                run, split, n_ev_str,
                mu_pre_str, sec_pre,
                mu_post_str, sec_post,
                sigma_str,
            )
        except Exception as exc:
            msg = f"run={run} split={split} TSV parse failed: {exc}"
            log.error(msg)
            errors.append(msg)

    # Excluded runs stay out of summary, plots and PDF even if stale scan_* folders exist.
    all_rows = [r for r in all_rows if int(r["run_number"]) not in _EXCLUDED_RUNS]

    if not all_rows:
        log.error("No rows collected — all benchmark runs failed or TSVs missing.")
        return 1

    # Sort by run then split
    all_rows.sort(key=lambda r: (int(r["run_number"]), int(r["split_number"])))

    # ---- Write aggregate TSV ------------------------------------------------
    summary_tsv = output_dir / f"{stem}_summary.tsv"
    write_summary_tsv(all_rows, summary_tsv)
    print(f"\nWrote summary TSV: {summary_tsv}", flush=True)
    log.info("Wrote summary TSV: %s", summary_tsv)

    # ---- Trend plots (pre-calibration and post-calibration) -----------------
    trend_pre_png  = (
        output_dir / f"{stem}_sigma_E_vs_file_pre.png"
    )
    trend_post_png = (
        output_dir / f"{stem}_sigma_E_vs_file_post.png"
    )
    try:
        # Pre-calibration plot (Gaussian + Crystal Ball overlaid)
        make_trend_plot(
            rows=all_rows,
            output_path=trend_pre_png,
            value_col_g="sigma_over_E_pre_G",
            value_col_cb="sigma_over_E_cb_pre",
            title_suffix="PRE-calibration",
            secondary_flag_col="secondary_peak_pre_present",
            pulse_method=lab,
        )
        print(f"\nWrote pre-calibration trend plot: {trend_pre_png}", flush=True)
        log.info("Wrote pre-calibration trend plot: %s", trend_pre_png)

        # Post-calibration plot (Gaussian + Crystal Ball overlaid)
        make_trend_plot(
            rows=all_rows,
            output_path=trend_post_png,
            value_col_g="sigma_over_E_post_G",
            value_col_cb="sigma_over_E_cb_post",
            title_suffix="POST-calibration",
            secondary_flag_col="secondary_peak_post_present",
            pulse_method=lab,
        )
        print(f"\nWrote post-calibration trend plot: {trend_post_png}", flush=True)
        log.info("Wrote post-calibration trend plot: %s", trend_post_png)

        if not trend_pre_png.exists() or not trend_post_png.exists():
            log.error(
                "One or both trend plots not generated: %s, %s",
                trend_pre_png, trend_post_png,
            )
            return 1
    except Exception as exc:
        log.error("Trend plot generation failed: %s", exc)
        return 1

    # ---- Combined calibrated-plot PDF ----------------------------------------
    combined_pdf = (
        output_dir / f"{stem}_all_calibrated_plots.pdf"
    )
    try:
        generate_combined_pdf(all_rows, output_dir, combined_pdf,
                              pulse_method=lab)
        log.info("Wrote combined calibrated PDF: %s", combined_pdf)
    except Exception as exc:
        log.error("Combined PDF generation failed: %s", exc)
        # Non-fatal: log and continue

    # ---- Findings markdown --------------------------------------------------
    findings_md = output_dir / f"{stem}_findings.md"
    write_findings_md(all_rows, findings_md, pulse_method=lab,
                      window=args.window, pre_window=pre_window, stem=stem)
    log.info("Wrote findings: %s", findings_md)

    # Separate counting for timeouts vs. other errors
    timeout_errors = [e for e in errors if "TIMEOUT" in e]
    other_errors   = [e for e in errors if "TIMEOUT" not in e]

    if timeout_errors:
        log.warning("%d file(s) timed out (hang bug — moved to error list):",
                    len(timeout_errors))
        for e in timeout_errors:
            log.warning("  %s", e)
        log.warning(
            "Consider moving these files to a quarantine directory to avoid "
            "re-processing them in future scans."
        )

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
        summary_msg += (
            f"  {len(errors)} error(s) total: "
            f"{len(timeout_errors)} timeout(s), {len(other_errors)} other.\n"
        )
    print(summary_msg, flush=True)
    log.info(summary_msg.strip())

    return 0 if not errors else 2


if __name__ == "__main__":
    sys.exit(main())
