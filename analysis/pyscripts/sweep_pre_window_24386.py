#!/usr/bin/env python3
"""sweep_pre_window_24386.py — Pre-window sweep driver for run 24386.

Invokes `benchmark_hycal_timing.py` over a grid of pre-window lower edges
(LO ns, fixed HI=200 ns) and aggregates per-window metrics into a single
summary TSV and three trend PNGs.

Reference design: docs/plans/plan-timing-gate-benchmark-24386.md §4 (sweep
grid), §5 (TSV schema), §6 (trend plots).

Usage example (run from repo root, with setup_env.sh sourced):
    python analysis/pyscripts/sweep_pre_window_24386.py

Override a single point for a cheap smoke-test:
    python analysis/pyscripts/sweep_pre_window_24386.py --pre-window-los 100

Skip already-completed LOs:
    python analysis/pyscripts/sweep_pre_window_24386.py --skip-existing

Inputs (defaults, all overridable via CLI):
  EVIO  : /Users/L00338853/work/PRad/data/evio/prad_024386.evio.00050
  calib : /Users/L00338853/work/PRad/data/physics_calib_out/Physics_calib/
            run24386/calib_factor_iter1.json
  grid  : LO ∈ {50, 75, 100, 125, 130, 150, 175} ns   HI = 200 ns

Outputs (under analysis/pyscripts/output/):
  timing_gate_benchmark_24386/          per-window .tsv / .png / .log files
  timing_gate_benchmark_24386.log       top-level sweep log
  timing_gate_benchmark_24386_summary.tsv  aggregate (14 rows, 2 paths × 7 LOs)
  timing_gate_benchmark_24386_curves_sigma_E.png
  timing_gate_benchmark_24386_curves_acceptance.png
  timing_gate_benchmark_24386_curves_mu.png
"""

from __future__ import annotations

import argparse
import csv
import datetime
import math
import subprocess
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ---------------------------------------------------------------------------
# Repository root is two levels above this file (analysis/pyscripts/).
# All relative paths are anchored here so the script is location-independent.
# ---------------------------------------------------------------------------
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Fixed parameters from hycal_clustering.md reference recipe
# (plan §4.2 — do NOT change without updating the plan).
# ---------------------------------------------------------------------------
_FIXED_ARGS = [
    "--signal-min", "500",
    "--calibrate",
    "--cal-method", "cb-peak",
    "--exclude-rowcol", "15", "20", "15", "20",
    "--target-peak", "3488.5",
    "--inner-radius", "200",
    "--min-seed-energy", "200",
    "--min-per-module", "30",
    "--bin-width", "20",
    "--e-range", "2500", "4000",
]

# Approved sweep grid (plan §4.1).
_DEFAULT_PRE_LOS = [50, 75, 100, 125, 130, 150, 175]
_DEFAULT_PRE_HI  = 200
_DEFAULT_WINDOW  = 16

_DEFAULT_EVIO = (
    "/Users/L00338853/work/PRad/data/evio/prad_024386.evio.00050"
)
_DEFAULT_CALIB = (
    "/Users/L00338853/work/PRad/data/physics_calib_out/"
    "Physics_calib/run24386/calib_factor_iter1.json"
)
_DEFAULT_MAX_EVENTS = 50000

_BENCHMARK_SCRIPT = str(
    _REPO_ROOT / "analysis" / "pyscripts" / "benchmark_hycal_timing.py"
)

# Output defaults (relative to repo root).
_DEFAULT_OUT_DIR = str(
    _REPO_ROOT / "analysis" / "pyscripts" / "output" / "timing_gate_benchmark_24386"
)
_DEFAULT_SUMMARY_TSV = str(
    _REPO_ROOT / "analysis" / "pyscripts" / "output" /
    "timing_gate_benchmark_24386_summary.tsv"
)
_DEFAULT_CURVES_PNG_STEM = str(
    _REPO_ROOT / "analysis" / "pyscripts" / "output" /
    "timing_gate_benchmark_24386_curves"
)
_DEFAULT_TOP_LOG = str(
    _REPO_ROOT / "analysis" / "pyscripts" / "output" /
    "timing_gate_benchmark_24386.log"
)

# Runinfo LO for reference marking in plots (plan §4.1 note).
_RUNINFO_LO_NS = 130


# ===========================================================================
# Step 1 — Run one benchmark invocation
# ===========================================================================

def run_one(lo: int, hi: int, cfg: argparse.Namespace, top_log_fh) -> Path:
    """Invoke benchmark_hycal_timing.py for a single (lo, hi) pair.

    Returns the path to the emitted .tsv file.  Raises RuntimeError if the
    subprocess exits non-zero (after capturing the log).
    """
    out_dir  = Path(cfg.out_dir)
    basename = f"pre_{lo}_{hi}"
    out_stem = out_dir / basename
    tsv_path = out_dir / f"{basename}.tsv"
    log_path = out_dir / f"{basename}.log"

    # Honour --skip-existing: reuse prior outputs if TSV already exists.
    if cfg.skip_existing and tsv_path.exists():
        msg = f"  [skip] {tsv_path} already exists — reusing.\n"
        print(msg, end="")
        top_log_fh.write(msg)
        top_log_fh.flush()
        return tsv_path

    cmd = [
        sys.executable,
        _BENCHMARK_SCRIPT,
        cfg.evio_path,
        str(out_stem),
        "--max-events", str(cfg.max_events),
        "--pre-window", str(lo), str(hi),
        "--window", str(cfg.window),
        "--hc-calib-file", cfg.hc_calib_file,
    ] + _FIXED_ARGS

    msg = (
        f"\n[run] pre-window [{lo}, {hi}] ns  →  {out_stem}\n"
        f"  cmd: {' '.join(cmd)}\n"
    )
    print(msg, end="")
    top_log_fh.write(msg)
    top_log_fh.flush()

    t0 = time.monotonic()
    try:
        with open(log_path, "w") as lf:
            result = subprocess.run(
                cmd,
                stdout=lf,
                stderr=subprocess.STDOUT,
                check=False,          # we inspect returncode ourselves
            )
        elapsed = time.monotonic() - t0
        status_msg = (
            f"  returncode={result.returncode}  elapsed={elapsed:.1f}s\n"
        )
        if result.returncode != 0:
            err = (
                f"[ERROR] benchmark failed for pre_window=[{lo},{hi}]  "
                f"(rc={result.returncode}).  "
                f"Log: {log_path}\n"
            )
            print(err, end="", file=sys.stderr)
            top_log_fh.write(err)
            top_log_fh.flush()
            raise RuntimeError(
                f"benchmark_hycal_timing.py exited with code "
                f"{result.returncode} for LO={lo}. See {log_path}"
            )
        print(status_msg, end="")
        top_log_fh.write(status_msg)
        top_log_fh.flush()
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"benchmark script not found: {_BENCHMARK_SCRIPT}"
        ) from exc

    if not tsv_path.exists():
        raise RuntimeError(
            f"Expected TSV was not created: {tsv_path}  "
            f"(check log: {log_path})"
        )
    return tsv_path


# ===========================================================================
# Step 2 — Parse one per-window TSV
# ===========================================================================

def _parse_float(s: str) -> float:
    """Convert string to float, treating 'nan', 'inf', empty → float('nan')."""
    s = s.strip()
    if not s:
        return float("nan")
    try:
        return float(s)
    except ValueError:
        return float("nan")


def parse_tsv(tsv_path: Path) -> dict:
    """Parse the structured TSV emitted by benchmark_hycal_timing.py.

    Returns a dict with keys ``headline`` and ``calibration``, each mapping
    path name ("legacy", "new") → dict of metric fields.

    Raises RuntimeError loudly if either block is absent (plan §4.5 contract).

    Headline block columns (after the header comment line):
        path  n_events  gauss_mu_MeV  gauss_sigma_MeV  n_in_fit

    Calibration block columns (after the # --calibrate header comment line):
        path  n_modules  n_inner_events  mu_pre_MeV  sigma_pre_MeV
        mu_post_MeV  sigma_post_MeV  sigma_over_E_post
        mu_cb_MeV  sigma_cb_MeV  alpha_cb  n_cb  sigma_over_E_cb
    """
    text = tsv_path.read_text(encoding="utf-8")
    lines = text.splitlines()

    # ---- locate the two data sections --------------------------------
    # Headline section: the line "path\tn_events\t..." that follows the
    # opening comment block.  It is always the first tab-separated data line.
    # Calibration section: lines following "# --calibrate inner ring ...".

    headline_hdr_lineno   = None
    calibration_hdr_lineno = None

    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("path\tn_events\t"):
            if headline_hdr_lineno is None:
                headline_hdr_lineno = i
        if stripped.startswith("# --calibrate inner ring"):
            calibration_hdr_lineno = i

    if headline_hdr_lineno is None:
        raise RuntimeError(
            f"Headline block missing in {tsv_path} "
            "(expected line starting 'path\\tn_events\\t...')"
        )
    if calibration_hdr_lineno is None:
        raise RuntimeError(
            f"Calibration block missing in {tsv_path} "
            "(expected comment line '# --calibrate inner ring ...'). "
            "Was --calibrate passed to the benchmark script?"
        )

    # ---- parse headline block ----------------------------------------
    headline_cols = lines[headline_hdr_lineno].strip().split("\t")
    headline: dict[str, dict] = {}
    for line in lines[headline_hdr_lineno + 1 :]:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            break
        parts = stripped.split("\t")
        if len(parts) < len(headline_cols):
            continue
        row = dict(zip(headline_cols, parts))
        path_name = row.pop("path")
        headline[path_name] = {k: _parse_float(v) for k, v in row.items()}
        # n_events and n_in_fit are integers — keep as float for uniformity.

    if not headline:
        raise RuntimeError(
            f"Headline block present in {tsv_path} but contained no data rows."
        )

    # ---- find the calibration data header (line after the comment) ----
    cal_data_hdr_lineno = None
    for i in range(calibration_hdr_lineno + 1, len(lines)):
        stripped = lines[i].strip()
        if stripped.startswith("path\tn_modules\t"):
            cal_data_hdr_lineno = i
            break

    if cal_data_hdr_lineno is None:
        raise RuntimeError(
            f"Calibration data header ('path\\tn_modules\\t...') not found "
            f"after the '# --calibrate' comment in {tsv_path}."
        )

    calibration_cols = lines[cal_data_hdr_lineno].strip().split("\t")
    calibration: dict[str, dict] = {}
    for line in lines[cal_data_hdr_lineno + 1 :]:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            break
        parts = stripped.split("\t")
        if len(parts) < len(calibration_cols):
            continue
        row = dict(zip(calibration_cols, parts))
        path_name = row.pop("path")
        calibration[path_name] = {k: _parse_float(v) for k, v in row.items()}

    if not calibration:
        raise RuntimeError(
            f"Calibration block present in {tsv_path} but contained no data rows."
        )

    return {"headline": headline, "calibration": calibration}


# ===========================================================================
# Step 3 — Write aggregate summary TSV
# ===========================================================================

# Aggregate TSV column order as per plan §5.
_SUMMARY_COLS = [
    "pre_lo_ns", "pre_hi_ns", "width_ns", "path",
    "n_events", "n_inner_events", "n_modules_used",
    "mu_pre_MeV", "sigma_pre_MeV", "mu_post_MeV", "sigma_post_MeV",
    "sigma_over_E_G",
    "mu_cb_MeV", "sigma_cb_MeV", "alpha_cb", "n_cb", "sigma_over_E_CB",
    "acceptance_vs_widest_pct",
    "gauss_mu_headline_MeV", "gauss_sigma_headline_MeV", "n_in_fit",
]


def write_summary(rows: list[dict], out_tsv: Path) -> None:
    """Write aggregate TSV with one row per (path, pre-window).

    ``rows`` is a list of dicts keyed by _SUMMARY_COLS.
    acceptance_vs_widest_pct must already be computed (plan §5).
    """
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=_SUMMARY_COLS, delimiter="\t",
            extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            # Format floats to a readable precision; keep NaN as "nan".
            out_row = {}
            for col in _SUMMARY_COLS:
                val = row.get(col, float("nan"))
                if isinstance(val, float):
                    if math.isnan(val):
                        out_row[col] = "nan"
                    else:
                        # Integer-like columns: no decimal.
                        if col in ("pre_lo_ns", "pre_hi_ns", "width_ns",
                                   "n_events", "n_inner_events",
                                   "n_modules_used", "n_in_fit"):
                            out_row[col] = f"{val:.0f}"
                        elif col in ("acceptance_vs_widest_pct",):
                            out_row[col] = f"{val:.2f}"
                        elif col in ("sigma_over_E_G", "sigma_over_E_CB"):
                            out_row[col] = f"{val:.4f}"
                        else:
                            out_row[col] = f"{val:.2f}"
                else:
                    out_row[col] = str(val)
            writer.writerow(out_row)


# ===========================================================================
# Step 4 — Build aggregate row list from parsed TSV dict
# ===========================================================================

def build_rows(
    lo: int, hi: int, parsed: dict
) -> list[dict]:
    """Convert a parsed TSV dict + (lo, hi) into two aggregate row dicts
    (one for "legacy", one for "new").  acceptance_vs_widest_pct is left
    as NaN here and filled in post-hoc once we have the LO=50 baseline.
    """
    width = hi - lo
    rows = []

    headline    = parsed["headline"]
    calibration = parsed["calibration"]

    for path_name in ("legacy", "new"):
        hl  = headline.get(path_name, {})
        cal = calibration.get(path_name, {})

        mu_post = cal.get("mu_post_MeV", float("nan"))
        sig_g   = cal.get("sigma_post_MeV", float("nan"))
        mu_cb   = cal.get("mu_cb_MeV", float("nan"))
        sig_cb  = cal.get("sigma_cb_MeV", float("nan"))

        # sigma_over_E_G = sigma_post / mu_post (Gaussian on calibrated ring)
        if mu_post > 0 and not math.isnan(mu_post) and not math.isnan(sig_g):
            sigma_over_E_G = sig_g / mu_post
        else:
            sigma_over_E_G = float("nan")

        # sigma_over_E_CB: benchmark script emits sigma_over_E_cb directly.
        sigma_over_E_CB = cal.get("sigma_over_E_cb", float("nan"))

        rows.append({
            "pre_lo_ns":   float(lo),
            "pre_hi_ns":   float(hi),
            "width_ns":    float(width),
            "path":        path_name,
            "n_events":    hl.get("n_events", float("nan")),
            "n_inner_events":  cal.get("n_inner_events", float("nan")),
            "n_modules_used":  cal.get("n_modules", float("nan")),
            "mu_pre_MeV":      cal.get("mu_pre_MeV", float("nan")),
            "sigma_pre_MeV":   cal.get("sigma_pre_MeV", float("nan")),
            "mu_post_MeV":     mu_post,
            "sigma_post_MeV":  sig_g,
            "sigma_over_E_G":  sigma_over_E_G,
            "mu_cb_MeV":       mu_cb,
            "sigma_cb_MeV":    sig_cb,
            "alpha_cb":        cal.get("alpha_cb", float("nan")),
            "n_cb":            cal.get("n_cb", float("nan")),
            "sigma_over_E_CB": sigma_over_E_CB,
            "acceptance_vs_widest_pct": float("nan"),  # filled later
            "gauss_mu_headline_MeV":    hl.get("gauss_mu_MeV", float("nan")),
            "gauss_sigma_headline_MeV": hl.get("gauss_sigma_MeV", float("nan")),
            "n_in_fit":    hl.get("n_in_fit", float("nan")),
        })
    return rows


def compute_acceptance(all_rows: list[dict]) -> None:
    """Fill acceptance_vs_widest_pct in-place.

    Baseline is n_events at LO=50 (the widest window), per path.
    If LO=50 is absent, acceptance stays NaN for all rows of that path.
    """
    baseline: dict[str, float] = {}
    for row in all_rows:
        if int(row["pre_lo_ns"]) == 50:
            path = row["path"]
            baseline[path] = row["n_events"]

    for row in all_rows:
        path = row["path"]
        n_base = baseline.get(path, float("nan"))
        n_this = row["n_events"]
        if (
            not math.isnan(n_base) and not math.isnan(n_this)
            and n_base > 0
        ):
            row["acceptance_vs_widest_pct"] = 100.0 * n_this / n_base
        else:
            row["acceptance_vs_widest_pct"] = float("nan")


# ===========================================================================
# Step 5 — Trend plots
# ===========================================================================

def _split_by_path(summary_tsv: Path):
    """Read summary TSV; return two dicts (legacy, new) each mapping
    width_ns → metric dict, sorted by width ascending (wide → narrow).
    """
    legacy_rows: dict[float, dict] = {}
    new_rows:    dict[float, dict] = {}

    with open(summary_tsv, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            w = _parse_float(row.get("width_ns", "nan"))
            pf = {k: _parse_float(v) for k, v in row.items()
                  if k != "path"}
            if row.get("path") == "legacy":
                legacy_rows[w] = pf
            elif row.get("path") == "new":
                new_rows[w] = pf

    return legacy_rows, new_rows


def _widths_and_metric(rows: dict[float, dict], metric: str):
    """Return (widths, values) arrays for a metric, excluding NaN pairs."""
    pairs = []
    for w, d in sorted(rows.items()):
        v = d.get(metric, float("nan"))
        if not math.isnan(w) and not math.isnan(v):
            pairs.append((w, v))
    if not pairs:
        return [], []
    ws, vs = zip(*pairs)
    return list(ws), list(vs)


def _runinfo_width() -> float:
    """Width corresponding to the runinfo LO=130, HI=200 baseline."""
    return float(_DEFAULT_PRE_HI - _RUNINFO_LO_NS)  # = 70 ns


def plot_curves(summary_tsv: Path, curves_stem: str) -> list[Path]:
    """Generate three trend PNGs from the aggregate TSV.

    Returns list of written PNG paths.
    """
    legacy_rows, new_rows = _split_by_path(summary_tsv)

    run_width = _runinfo_width()  # 70 ns
    target_peak = 3488.5          # MeV

    written = []

    # ----------------------------------------------------------------
    # Plot 1 — σ_E/E (CB) vs. width_ns
    # ----------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))

    wl, vl = _widths_and_metric(legacy_rows, "sigma_over_E_CB")
    wn, vn = _widths_and_metric(new_rows,    "sigma_over_E_CB")

    if wl:
        ax.plot(wl, [v * 100 for v in vl], "o-",
                color="#555555", lw=2, label="legacy (σ_E/E CB)")
    if wn:
        ax.plot(wn, [v * 100 for v in vn], "s--",
                color="#1f77b4", lw=2, label=f"gated W={_DEFAULT_WINDOW} ns (σ_E/E CB)")

    ax.axvline(run_width, color="red", ls="--", lw=1.2,
               label=f"runinfo baseline (LO={_RUNINFO_LO_NS} ns, width={run_width:.0f} ns)")
    ax.set_xlabel("Pre-window width (ns)  [HI − LO,  HI = 200 ns]")
    ax.set_ylabel("σ_E/E CB (%)")
    ax.set_title(
        "Run 24386 — σ_E/E (Crystal Ball) vs. pre-window width\n"
        "(HI=200 ns, W=16 ns, iter1 calib, split 00050, 50k events)"
    )
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    ax.invert_xaxis()  # wide window on left → narrow on right
    fig.tight_layout()
    p1 = Path(f"{curves_stem}_sigma_E.png")
    fig.savefig(p1, dpi=130)
    plt.close(fig)
    written.append(p1)

    # ----------------------------------------------------------------
    # Plot 2 — Event acceptance vs. width_ns
    # ----------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))

    wl, vl = _widths_and_metric(legacy_rows, "acceptance_vs_widest_pct")
    wn, vn = _widths_and_metric(new_rows,    "acceptance_vs_widest_pct")

    if wl:
        ax.plot(wl, vl, "o-", color="#555555", lw=2,
                label="legacy (n_events / n_widest × 100 %)")
    if wn:
        ax.plot(wn, vn, "s--", color="#1f77b4", lw=2,
                label=f"gated W={_DEFAULT_WINDOW} ns")

    ax.axhline(90.0, color="darkorange", ls=":", lw=1.3,
               label="90 % acceptance threshold")
    ax.axvline(run_width, color="red", ls="--", lw=1.2,
               label=f"runinfo baseline (LO={_RUNINFO_LO_NS} ns, width={run_width:.0f} ns)")
    ax.set_xlabel("Pre-window width (ns)  [HI − LO,  HI = 200 ns]")
    ax.set_ylabel("Event acceptance (% of widest window)")
    ax.set_title(
        "Run 24386 — Event acceptance vs. pre-window width\n"
        "(HI=200 ns, W=16 ns, iter1 calib, split 00050, 50k events)"
    )
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    ax.invert_xaxis()
    fig.tight_layout()
    p2 = Path(f"{curves_stem}_acceptance.png")
    fig.savefig(p2, dpi=130)
    plt.close(fig)
    written.append(p2)

    # ----------------------------------------------------------------
    # Plot 3 — Peak position (mu_post_MeV) vs. width_ns
    # ----------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))

    wl_pre, vl_pre = _widths_and_metric(legacy_rows, "mu_pre_MeV")
    wl_pst, vl_pst = _widths_and_metric(legacy_rows, "mu_post_MeV")
    wn_pre, vn_pre = _widths_and_metric(new_rows,    "mu_pre_MeV")
    wn_pst, vn_pst = _widths_and_metric(new_rows,    "mu_post_MeV")

    if wl_pre:
        ax.plot(wl_pre, vl_pre, "o:", color="#999999", lw=1.5,
                label="legacy μ_pre (MeV)")
    if wl_pst:
        ax.plot(wl_pst, vl_pst, "o-", color="#555555", lw=2,
                label="legacy μ_post (MeV)")
    if wn_pre:
        ax.plot(wn_pre, vn_pre, "s:", color="#7fbfff", lw=1.5,
                label="gated μ_pre (MeV)")
    if wn_pst:
        ax.plot(wn_pst, vn_pst, "s--", color="#1f77b4", lw=2,
                label=f"gated μ_post (MeV)  W={_DEFAULT_WINDOW} ns")

    ax.axhline(target_peak, color="green", ls=":", lw=1.2,
               label=f"target peak {target_peak:.1f} MeV")
    ax.axvline(run_width, color="red", ls="--", lw=1.2,
               label=f"runinfo baseline (LO={_RUNINFO_LO_NS} ns, width={run_width:.0f} ns)")

    # Zoom y-axis to ±200 MeV around target.
    ax.set_ylim(target_peak - 200, target_peak + 200)
    ax.set_xlabel("Pre-window width (ns)  [HI − LO,  HI = 200 ns]")
    ax.set_ylabel("Peak position (MeV)")
    ax.set_title(
        "Run 24386 — Peak position vs. pre-window width\n"
        "(HI=200 ns, W=16 ns, iter1 calib, split 00050, 50k events)"
    )
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    ax.invert_xaxis()
    fig.tight_layout()
    p3 = Path(f"{curves_stem}_mu.png")
    fig.savefig(p3, dpi=130)
    plt.close(fig)
    written.append(p3)

    return written


# ===========================================================================
# CLI
# ===========================================================================

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--pre-window-los", type=int, nargs="+",
        default=_DEFAULT_PRE_LOS,
        metavar="LO",
        help=(
            f"Pre-window lower edges to sweep (ns). "
            f"Default: {_DEFAULT_PRE_LOS}"
        ),
    )
    ap.add_argument(
        "--pre-hi", type=int, default=_DEFAULT_PRE_HI,
        metavar="HI",
        help=f"Pre-window upper edge (ns). Default: {_DEFAULT_PRE_HI}",
    )
    ap.add_argument(
        "--window", type=int, default=_DEFAULT_WINDOW,
        help=f"seed_time_window for the gated clusterer (ns). Default: {_DEFAULT_WINDOW}",
    )
    ap.add_argument(
        "--evio-path", default=_DEFAULT_EVIO,
        help=f"EVIO split to process. Default: {_DEFAULT_EVIO}",
    )
    ap.add_argument(
        "--max-events", type=int, default=_DEFAULT_MAX_EVENTS,
        help=f"--max-events forwarded to benchmark script. Default: {_DEFAULT_MAX_EVENTS}",
    )
    ap.add_argument(
        "--hc-calib-file", default=_DEFAULT_CALIB,
        help=f"HyCal calibration JSON. Default: {_DEFAULT_CALIB}",
    )
    ap.add_argument(
        "--out-dir", default=_DEFAULT_OUT_DIR,
        help=(
            "Directory for per-window outputs (.tsv/.png/.log). "
            f"Default: {_DEFAULT_OUT_DIR}"
        ),
    )
    ap.add_argument(
        "--summary-tsv", default=_DEFAULT_SUMMARY_TSV,
        help=f"Aggregate summary TSV path. Default: {_DEFAULT_SUMMARY_TSV}",
    )
    ap.add_argument(
        "--curves-png", default=_DEFAULT_CURVES_PNG_STEM,
        help=(
            "Stem for the three trend PNG files "
            "(suffixes _sigma_E / _acceptance / _mu added). "
            f"Default: {_DEFAULT_CURVES_PNG_STEM}"
        ),
    )
    ap.add_argument(
        "--skip-existing", action="store_true",
        help="Skip re-running a benchmark if the per-window .tsv already exists.",
    )
    return ap


# ===========================================================================
# Main
# ===========================================================================

def main(argv=None) -> int:
    ap = build_parser()
    cfg = ap.parse_args(argv)

    # ---- Pre-flight checks --------------------------------------------------
    calib_path = Path(cfg.hc_calib_file)
    if not calib_path.exists():
        print(
            f"[ERROR] Calibration file not found: {calib_path}\n"
            "        Check --hc-calib-file or verify the iter1 JSON path.",
            file=sys.stderr,
        )
        return 1

    evio_path = Path(cfg.evio_path)
    if not evio_path.exists():
        print(
            f"[ERROR] EVIO file not found: {evio_path}\n"
            "        Check --evio-path.",
            file=sys.stderr,
        )
        return 1

    benchmark_path = Path(_BENCHMARK_SCRIPT)
    if not benchmark_path.exists():
        print(
            f"[ERROR] Benchmark script not found: {benchmark_path}",
            file=sys.stderr,
        )
        return 1

    # ---- Create output directory -------------------------------------------
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- Top-level log file ------------------------------------------------
    top_log_path = Path(_DEFAULT_TOP_LOG)
    top_log_path.parent.mkdir(parents=True, exist_ok=True)
    top_log_fh = open(top_log_path, "a", encoding="utf-8")

    start_ts = datetime.datetime.now().isoformat(timespec="seconds")
    header = (
        f"\n{'='*72}\n"
        f"sweep_pre_window_24386.py  started {start_ts}\n"
        f"  evio        : {cfg.evio_path}\n"
        f"  calib       : {cfg.hc_calib_file}\n"
        f"  pre_los     : {cfg.pre_window_los}\n"
        f"  pre_hi      : {cfg.pre_hi}\n"
        f"  window      : {cfg.window}\n"
        f"  max_events  : {cfg.max_events}\n"
        f"  out_dir     : {cfg.out_dir}\n"
        f"  summary_tsv : {cfg.summary_tsv}\n"
        f"  curves_stem : {cfg.curves_png}\n"
        f"  skip_existing: {cfg.skip_existing}\n"
        f"{'='*72}\n"
    )
    print(header, end="")
    top_log_fh.write(header)
    top_log_fh.flush()

    # ---- Sweep loop --------------------------------------------------------
    all_rows: list[dict] = []
    errors: list[str]    = []

    for lo in cfg.pre_window_los:
        hi = cfg.pre_hi
        print(f"\nRunning pre-window [{lo}, {hi}]...", flush=True)
        top_log_fh.write(f"\nRunning pre-window [{lo}, {hi}]...\n")
        top_log_fh.flush()

        try:
            tsv_path = run_one(lo, hi, cfg, top_log_fh)
            parsed   = parse_tsv(tsv_path)
            rows     = build_rows(lo, hi, parsed)
            all_rows.extend(rows)

            # Flag modules-used drop >10 % vs LO=50 reference (plan §10 risk).
            # We do this after we have the LO=50 row in place if it's the
            # first LO; otherwise we compare incrementally below.
            msg = (
                f"  parsed {tsv_path.name}: "
                f"legacy n_events={rows[0]['n_events']:.0f}  "
                f"new n_events={rows[1]['n_events']:.0f}\n"
            )
            print(msg, end="", flush=True)
            top_log_fh.write(msg)
            top_log_fh.flush()

        except RuntimeError as exc:
            err = f"[ERROR] {exc}\n"
            print(err, end="", file=sys.stderr)
            top_log_fh.write(err)
            top_log_fh.flush()
            errors.append(str(exc))
            # Continue to next LO — do not abort the sweep.

    if not all_rows:
        print(
            "[ERROR] No rows collected — all benchmark runs failed. "
            "Check per-window logs.",
            file=sys.stderr,
        )
        top_log_fh.close()
        return 1

    # ---- Post-hoc: flag modules-used drops > 10 % from widest window --------
    # Find LO=50 reference n_modules_used per path.
    ref_mods: dict[str, float] = {}
    for row in all_rows:
        if int(row["pre_lo_ns"]) == 50:
            ref_mods[row["path"]] = row["n_modules_used"]

    for row in all_rows:
        p = row["path"]
        ref = ref_mods.get(p, float("nan"))
        nm  = row["n_modules_used"]
        if (
            not math.isnan(ref) and not math.isnan(nm)
            and ref > 0 and nm < 0.9 * ref
        ):
            warn = (
                f"[WARN] n_modules_used dropped >10 %: "
                f"path={p}  LO={row['pre_lo_ns']:.0f}  "
                f"n_modules={nm:.0f}  ref(LO=50)={ref:.0f}\n"
            )
            print(warn, end="", file=sys.stderr)
            top_log_fh.write(warn)
            top_log_fh.flush()

    # ---- Compute acceptance -------------------------------------------------
    compute_acceptance(all_rows)

    # ---- Write aggregate TSV ------------------------------------------------
    print("\nAggregating...", flush=True)
    summary_tsv = Path(cfg.summary_tsv)
    summary_tsv.parent.mkdir(parents=True, exist_ok=True)
    write_summary(all_rows, summary_tsv)
    print(f"  wrote: {summary_tsv}", flush=True)
    top_log_fh.write(f"  wrote summary TSV: {summary_tsv}\n")

    # ---- Trend plots --------------------------------------------------------
    print("\nPlotting...", flush=True)
    try:
        pngs = plot_curves(summary_tsv, cfg.curves_png)
        for p in pngs:
            print(f"  wrote: {p}", flush=True)
            top_log_fh.write(f"  wrote trend PNG: {p}\n")
    except Exception as exc:
        warn = f"[WARN] Plot generation failed: {exc}\n"
        print(warn, end="", file=sys.stderr)
        top_log_fh.write(warn)

    # ---- Final status -------------------------------------------------------
    end_ts = datetime.datetime.now().isoformat(timespec="seconds")
    if errors:
        footer = (
            f"\n[DONE with errors]  {end_ts}\n"
            f"  {len(errors)} LO(s) failed:\n"
            + "".join(f"    {e}\n" for e in errors)
            + f"  {len(all_rows)} rows written to {summary_tsv}\n"
        )
        print(footer)
        top_log_fh.write(footer)
        top_log_fh.close()
        return 2  # partial success — outputs exist but some LOs failed
    else:
        footer = (
            f"\n[DONE]  {end_ts}\n"
            f"  {len(all_rows)} rows written to {summary_tsv}\n"
        )
        print(footer)
        top_log_fh.write(footer)
        top_log_fh.close()
        return 0


if __name__ == "__main__":
    sys.exit(main())
