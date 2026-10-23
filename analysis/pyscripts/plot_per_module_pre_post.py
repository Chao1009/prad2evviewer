#!/usr/bin/env python3
"""plot_per_module_pre_post.py — Per-module pre/post calibration energy spectra.

Diagnoses which PbWO4 modules contribute to anomalous features (e.g., the
~2600 MeV secondary bump in run 25320) by plotting each module's pre- and
post-calibration energy distributions side-by-side with Crystal Ball fits.

Input: a clusters TSV produced by benchmark_hycal_timing.py --write-cluster-tsv
  Columns: event_num  config_name  seed_module_id  seed_x_mm  seed_y_mm
           seed_row   seed_col     seed_energy_MeV  cluster_E_raw_MeV
           cluster_E_cal_MeV

Output: a multi-page PDF with one summary page followed by one page per
PbWO4 module (W-prefix IDs only) that passes the minimum-event threshold.

Usage:
  python3 plot_per_module_pre_post.py [options]
  python3 plot_per_module_pre_post.py --help
"""

from __future__ import annotations

import argparse
import math
import sys
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
    from scipy.optimize import curve_fit
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

# ---------------------------------------------------------------------------
# Repository root — for locating benchmark_hycal_timing on the import path
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT   = _SCRIPT_DIR.parent.parent

# ---------------------------------------------------------------------------
# Try importing CB helpers from the sibling script to avoid duplication
# ---------------------------------------------------------------------------
def _import_cb_from_benchmark():
    """Attempt to import crystal_ball and fit_peak_cb from the sibling script."""
    import importlib.util
    bench_path = _SCRIPT_DIR / "benchmark_hycal_timing.py"
    if not bench_path.exists():
        return None, None
    try:
        spec = importlib.util.spec_from_file_location("benchmark_hycal_timing",
                                                       str(bench_path))
        mod  = importlib.util.module_from_spec(spec)
        # Benchmark imports prad2py which may not be available in a pure
        # post-processing context.  Catch ALL exceptions including SystemExit
        # (raised by _common.py when prad2py is absent) so we fall back to
        # the local copies below.
        spec.loader.exec_module(mod)
        return mod.crystal_ball, mod.fit_peak_cb
    except BaseException:
        return None, None


_cb_from_bench, _fit_peak_cb_from_bench = _import_cb_from_benchmark()


# ---------------------------------------------------------------------------
# Local fallback implementations (identical algorithm to benchmark_hycal_timing)
# ---------------------------------------------------------------------------

def crystal_ball(x, A, mu, sigma, alpha, n):
    """Crystal-Ball PDF (left-side power-law tail).  α > 0, n > 1."""
    z      = (np.asarray(x, dtype=float) - mu) / sigma
    abs_a  = abs(alpha)
    A_coef = (n / abs_a) ** n * np.exp(-0.5 * abs_a ** 2)
    B_coef = n / abs_a - abs_a
    gauss_part = np.exp(-0.5 * z * z)
    arg = B_coef - z
    arg = np.where(arg > 0, arg, 1e-30)
    tail_part = A_coef * arg ** (-n)
    return A * np.where(z > -alpha, gauss_part, tail_part)


def gauss(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def fit_peak_cb(values, e_lo, e_hi, bin_width, mu0_hint=None):
    """Crystal-Ball fit of values in [e_lo, e_hi].

    Returns (mu, sigma, alpha, n, A, n_in_fit, edges, counts).
    mu seeded from mu0_hint (else histogram argmax).
    """
    edges  = np.arange(e_lo, e_hi + bin_width, bin_width)
    counts, _ = np.histogram(values, bins=edges)
    centres = 0.5 * (edges[:-1] + edges[1:])
    nan = float("nan")
    if not HAVE_SCIPY or counts.sum() == 0:
        return nan, nan, nan, nan, nan, 0, edges, counts

    imax = int(np.argmax(counts))
    mu0  = float(mu0_hint) if mu0_hint is not None else float(centres[imax])
    win  = max(400.0, 12.0 * bin_width)
    mask = (centres > mu0 - win) & (centres < mu0 + 0.5 * win)
    x = centres[mask]
    y = counts[mask].astype(float)
    if y.sum() < 50:
        return nan, nan, nan, nan, nan, int(y.sum()), edges, counts
    try:
        p0     = [float(y.max()), mu0, max(60.0, bin_width * 3), 1.5, 2.0]
        bounds = ([0,         e_lo - 200.0, 1.0,   0.1,   1.01],
                  [np.inf,   e_hi + 200.0, win,  20.0,  50.0])
        popt, _ = curve_fit(crystal_ball, x, y, p0=p0, bounds=bounds,
                             maxfev=10000)
        A, mu, sig, alpha, n = popt
        return (float(mu), abs(float(sig)), float(alpha), float(n),
                float(A), int(y.sum()), edges, counts)
    except Exception as exc:
        sys.stderr.write(f"[WARN] crystal-ball fit failed: {exc}\n")
        return nan, nan, nan, nan, nan, int(y.sum()), edges, counts


def fit_peak_gauss(values, e_lo, e_hi, bin_width, mu0_hint=None):
    """Fallback Gaussian fit.

    Returns (mu, sigma, A, n_in_fit, edges, counts).
    """
    edges  = np.arange(e_lo, e_hi + bin_width, bin_width)
    counts, _ = np.histogram(values, bins=edges)
    centres = 0.5 * (edges[:-1] + edges[1:])
    nan = float("nan")
    if not HAVE_SCIPY or counts.sum() == 0:
        return nan, nan, nan, 0, edges, counts

    imax = int(np.argmax(counts))
    mu0  = float(mu0_hint) if mu0_hint is not None else float(centres[imax])
    win  = max(150.0, 4.0 * bin_width)
    mask = (centres > mu0 - win) & (centres < mu0 + win)
    x = centres[mask]
    y = counts[mask].astype(float)
    if y.sum() < 20:
        return nan, nan, nan, int(y.sum()), edges, counts
    try:
        p0     = [float(y.max()), mu0, max(bin_width, 60.0)]
        bounds = ([0, e_lo, 1.0], [np.inf, e_hi, win])
        popt, _ = curve_fit(gauss, x, y, p0=p0, bounds=bounds,
                             maxfev=5000)
        A, mu, sig = popt
        return float(mu), abs(float(sig)), float(A), int(y.sum()), edges, counts
    except Exception as exc:
        sys.stderr.write(f"[WARN] gauss fit failed: {exc}\n")
        return nan, nan, nan, int(y.sum()), edges, counts


# ---------------------------------------------------------------------------
# If the benchmark import succeeded, prefer those implementations
# ---------------------------------------------------------------------------
if _cb_from_bench is not None:
    crystal_ball  = _cb_from_bench           # noqa: F811
if _fit_peak_cb_from_bench is not None:
    fit_peak_cb   = _fit_peak_cb_from_bench  # noqa: F811


# ---------------------------------------------------------------------------
# Module decoding helpers
# ---------------------------------------------------------------------------

def _module_sort_key(name: str) -> tuple:
    """Natural sort key: 'W10' < 'W100'."""
    prefix = name[0] if name else "?"
    try:
        number = int(name[1:])
    except ValueError:
        number = 0
    return (prefix, number)


def _secondary_bump_present(values, target_peak, lo_offset=200.0,
                             hi_mev=2800.0):
    """Return True if the POST spectrum has non-trivial content in the
    bump region (target_peak+lo_offset … hi_mev)."""
    if len(values) < 5:
        return False
    bump_lo = target_peak + lo_offset
    bump_hi = hi_mev
    n_bump  = int(np.sum((np.asarray(values) >= bump_lo) &
                         (np.asarray(values) <= bump_hi)))
    n_total = len(values)
    # Flag if bump region contains >3 % of the events and >5 events
    return n_bump > 5 and (n_bump / max(n_total, 1)) > 0.03


# ---------------------------------------------------------------------------
# Per-panel plotting
# ---------------------------------------------------------------------------

def _draw_panel(ax, values, label, e_lo, e_hi, bin_width, target_peak,
                col_hist="#1f77b4", col_fit="#d62728",
                mu_warn_distance=None, diagnose_secondary=False,
                secondary_region=(2440.0, 2800.0)):
    """Draw one histogram + CB fit + reference line panel.

    Returns dict with fit results: mu, sigma, alpha, n, fit_ok, fit_type.
    """
    nan = float("nan")

    # Build histogram
    edges  = np.arange(e_lo, e_hi + bin_width, bin_width)
    counts, _ = np.histogram(values, bins=edges)
    centres = 0.5 * (edges[:-1] + edges[1:])
    n_events = len(values)

    # Draw histogram as step plot
    ax.step(centres, counts, where="mid", color=col_hist, lw=1.4,
            label=f"N = {n_events}")

    # Try CB fit
    fit_ok   = False
    fit_type = "none"
    mu, sigma, alpha_cb, n_cb, A_cb = nan, nan, nan, nan, nan

    if n_events >= 50:
        # mu0 hint: histogram argmax
        mu0_hint = float(centres[int(np.argmax(counts))]) if counts.sum() > 0 else None
        result = fit_peak_cb(values, e_lo, e_hi, bin_width, mu0_hint=mu0_hint)
        mu_cb, sig_cb, a_cb, nc_cb, Acb, _nfit, _edges, _counts = result
        if not math.isnan(mu_cb) and Acb > 0:
            mu, sigma, alpha_cb, n_cb, A_cb = mu_cb, sig_cb, a_cb, nc_cb, Acb
            fit_ok   = True
            fit_type = "CB"
        else:
            # Fallback: Gaussian
            gresult = fit_peak_gauss(values, e_lo, e_hi, bin_width,
                                     mu0_hint=mu0_hint)
            mu_g, sig_g, Ag, _nf, _e, _c = gresult
            if not math.isnan(mu_g) and Ag > 0:
                mu, sigma, A_cb = mu_g, sig_g, Ag
                fit_ok   = True
                fit_type = "Gauss"

    # Draw fit curve
    if fit_ok:
        x_fit = np.linspace(e_lo, e_hi, 500)
        if fit_type == "CB":
            y_fit = crystal_ball(x_fit, A_cb, mu, sigma, alpha_cb, n_cb)
        else:
            y_fit = gauss(x_fit, A_cb, mu, sigma)
        ax.plot(x_fit, y_fit, color=col_fit, lw=1.5,
                label=f"{fit_type} fit")

    # Reference line at target peak
    ax.axvline(target_peak, color="#444", lw=0.9, ls="--",
               label=f"target {target_peak:.0f} MeV")

    # Optional: shade secondary-bump region
    if diagnose_secondary:
        ax.axvspan(secondary_region[0], secondary_region[1],
                   alpha=0.10, color="orange", zorder=0)

    # Fit parameters text box
    if fit_ok and not math.isnan(mu):
        ratio = sigma / mu * 100.0 if mu > 0 else nan
        if fit_type == "CB":
            txt = (f"μ = {mu:.1f} MeV\n"
                   f"σ = {sigma:.1f} MeV\n"
                   f"α = {alpha_cb:.2f}\n"
                   f"n = {n_cb:.2f}\n"
                   f"σ/μ = {ratio:.2f} %")
        else:
            txt = (f"μ = {mu:.1f} MeV\n"
                   f"σ = {sigma:.1f} MeV\n"
                   f"σ/μ = {ratio:.2f} %\n"
                   f"[Gaussian fallback]")
        # Box background color based on diagnostic criteria
        box_color = "lightyellow"
        if mu_warn_distance is not None:
            dist = abs(mu - target_peak)
            if dist > mu_warn_distance:
                box_color = "lightsalmon"  # >threshold from target
        props = dict(boxstyle="round", facecolor=box_color, alpha=0.85)
        ax.text(0.97, 0.97, txt, transform=ax.transAxes,
                fontsize=7.5, va="top", ha="right",
                bbox=props)
    elif n_events > 0:
        ax.text(0.97, 0.97, "fit failed", transform=ax.transAxes,
                fontsize=8, va="top", ha="right",
                bbox=dict(boxstyle="round", facecolor="#ffcccc", alpha=0.8))

    ax.set_xlabel("cluster energy (MeV)", fontsize=8)
    ax.set_ylabel("events / bin", fontsize=8)
    ax.set_xlim(e_lo, e_hi)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_title(label, fontsize=9)

    return dict(mu=mu, sigma=sigma, alpha=alpha_cb, n=n_cb, A=A_cb,
                fit_ok=fit_ok, fit_type=fit_type, n_events=n_events)


# ---------------------------------------------------------------------------
# Summary page
# ---------------------------------------------------------------------------

def _draw_summary_page(pdf, module_results, target_peak):
    """Emit a summary page with aggregate distributions of fit parameters."""
    mu_pre_list  = [r["mu_pre"]  for r in module_results
                    if not math.isnan(r.get("mu_pre",  float("nan")))]
    mu_post_list = [r["mu_post"] for r in module_results
                    if not math.isnan(r.get("mu_post", float("nan")))]
    sig_pre_list  = [r["sigma_pre"]  for r in module_results
                     if not math.isnan(r.get("sigma_pre",  float("nan")))]
    sig_post_list = [r["sigma_post"] for r in module_results
                     if not math.isnan(r.get("sigma_post", float("nan")))]
    # Gain factor: raw μ / target (how much raw differs from calibrated target)
    gain_list = [r["mu_pre"] / target_peak
                 for r in module_results
                 if not math.isnan(r.get("mu_pre", float("nan")))
                 and target_peak > 0]

    fig = plt.figure(figsize=(12, 10))
    fig.suptitle(
        f"Per-Module Calibration Summary  —  {len(module_results)} modules plotted\n"
        f"Target peak: {target_peak:.1f} MeV",
        fontsize=12, y=0.98
    )

    gs  = fig.add_gridspec(3, 2, hspace=0.45, wspace=0.35,
                           top=0.92, bottom=0.07, left=0.09, right=0.97)
    axes = [fig.add_subplot(gs[i, j]) for i in range(3) for j in range(2)]

    # ---- Panel 0: Distribution of pre-cal μ ----------------------------------
    ax = axes[0]
    if mu_pre_list:
        ax.hist(mu_pre_list, bins=30, color="#888", edgecolor="white",
                linewidth=0.4)
        ax.axvline(target_peak, color="#d62728", ls="--", lw=1.0,
                   label=f"target {target_peak:.0f}")
    ax.set_xlabel("pre-cal μ (MeV)", fontsize=8)
    ax.set_ylabel("modules", fontsize=8)
    ax.set_title("Pre-calibration μ distribution", fontsize=8)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.25)

    # ---- Panel 1: Distribution of post-cal μ ---------------------------------
    ax = axes[1]
    if mu_post_list:
        ax.hist(mu_post_list, bins=30, color="#1f77b4", edgecolor="white",
                linewidth=0.4)
        ax.axvline(target_peak, color="#d62728", ls="--", lw=1.0,
                   label=f"target {target_peak:.0f}")
    ax.set_xlabel("post-cal μ (MeV)", fontsize=8)
    ax.set_ylabel("modules", fontsize=8)
    ax.set_title("Post-calibration μ distribution  (should peak at target)", fontsize=8)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.25)

    # ---- Panel 2: Gain factors (μ_pre / target_peak) ------------------------
    ax = axes[2]
    if gain_list:
        ax.hist(gain_list, bins=30, color="#2ca02c", edgecolor="white",
                linewidth=0.4)
        ax.axvline(1.0, color="#d62728", ls="--", lw=1.0, label="gain = 1")
    ax.set_xlabel("pre-cal μ / target peak", fontsize=8)
    ax.set_ylabel("modules", fontsize=8)
    ax.set_title("Gain factor distribution  (μ_pre / target)", fontsize=8)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.25)

    # ---- Panel 3: pre-cal σ distribution ------------------------------------
    ax = axes[3]
    if sig_pre_list:
        ax.hist(sig_pre_list, bins=30, color="#ff7f0e", edgecolor="white",
                linewidth=0.4)
    if sig_post_list:
        ax.hist(sig_post_list, bins=30, color="#1f77b4", edgecolor="white",
                linewidth=0.4, alpha=0.6)
    ax.set_xlabel("σ (MeV)", fontsize=8)
    ax.set_ylabel("modules", fontsize=8)
    ax.set_title("σ distribution: pre (orange) vs post (blue)", fontsize=8)
    ax.grid(alpha=0.25)

    # ---- Panel 4: Scatter pre-cal μ vs post-cal μ ---------------------------
    ax = axes[4]
    mu_pre_arr  = np.array([r["mu_pre"]  for r in module_results
                            if not math.isnan(r.get("mu_pre",  float("nan")))
                            and not math.isnan(r.get("mu_post", float("nan")))])
    mu_post_arr = np.array([r["mu_post"] for r in module_results
                            if not math.isnan(r.get("mu_pre",  float("nan")))
                            and not math.isnan(r.get("mu_post", float("nan")))])
    if len(mu_pre_arr) > 0:
        ax.scatter(mu_pre_arr, mu_post_arr, s=12, alpha=0.6, color="#1f77b4")
        ax.axhline(target_peak, color="#d62728", ls="--", lw=0.9,
                   label=f"y = target {target_peak:.0f}")
    ax.set_xlabel("pre-cal μ (MeV)", fontsize=8)
    ax.set_ylabel("post-cal μ (MeV)", fontsize=8)
    ax.set_title("Scatter: pre-cal μ vs post-cal μ\n(post should be flat at target)",
                 fontsize=8)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.25)

    # ---- Panel 5: Scatter pre-cal μ vs post-cal σ ---------------------------
    ax = axes[5]
    sig_post_arr = np.array([r["sigma_post"] for r in module_results
                              if not math.isnan(r.get("mu_pre",   float("nan")))
                              and not math.isnan(r.get("sigma_post", float("nan")))])
    mu_pre_arr2  = np.array([r["mu_pre"]  for r in module_results
                              if not math.isnan(r.get("mu_pre",   float("nan")))
                              and not math.isnan(r.get("sigma_post", float("nan")))])
    if len(mu_pre_arr2) > 0:
        ax.scatter(mu_pre_arr2, sig_post_arr, s=12, alpha=0.6, color="#9467bd")
    ax.set_xlabel("pre-cal μ (MeV)", fontsize=8)
    ax.set_ylabel("post-cal σ (MeV)", fontsize=8)
    ax.set_title("Scatter: pre-cal μ vs post-cal σ\n(large pre-offset → inflated post-σ?)",
                 fontsize=8)
    ax.grid(alpha=0.25)

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
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
        "--output-pdf",
        default=None,
        help="Output PDF path.  Default: derived from input filename with "
             "'_pre_post_spectra.pdf' suffix.",
    )
    ap.add_argument(
        "--min-events", type=int, default=50,
        help="Minimum events per module to generate a page (default: 50)",
    )
    ap.add_argument(
        "--max-modules", type=int, default=0,
        help="Maximum number of modules to plot (0 = unlimited)",
    )
    ap.add_argument(
        "--sort-by",
        choices=["name", "n_events", "mu_pre_distance"],
        default="name",
        help="Page ordering: 'name' (W1, W2, ...), 'n_events' (descending), "
             "'mu_pre_distance' (furthest pre-cal μ from target first).  "
             "Default: name",
    )
    ap.add_argument(
        "--e-range", type=float, nargs=2, default=(500.0, 3500.0),
        metavar=("LO", "HI"),
        help="Energy histogram range in MeV (default: 500 3500)",
    )
    ap.add_argument(
        "--bin-width", type=float, default=20.0,
        help="Histogram bin width in MeV (default: 20)",
    )
    ap.add_argument(
        "--target-peak", type=float, default=2239.5,
        help="Target peak energy for reference line and diagnostics "
             "(default: 2239.5 MeV)",
    )
    ap.add_argument(
        "--config", choices=["legacy", "new"], default="legacy",
        help="Which config_name to use: 'legacy' or 'new' (default: legacy)",
    )
    ap.add_argument(
        "--diagnose-secondary", action="store_true",
        help="Highlight modules whose POST spectrum shows a visible bump above "
             "target_peak+200 MeV (secondary-bump diagnostic; marks pages with "
             "a red border and asterisk in the title).",
    )
    return ap


def main(argv=None) -> int:
    ap   = build_parser()
    args = ap.parse_args(argv)

    # ---- Resolve paths ------------------------------------------------------
    input_path = Path(args.input_tsv)
    # Support paths relative to the script directory (common usage pattern)
    if not input_path.is_absolute():
        candidate = _SCRIPT_DIR / input_path
        if candidate.exists():
            input_path = candidate
        else:
            # Also try relative to CWD (user may run from repo root)
            cwd_candidate = Path.cwd() / input_path
            if cwd_candidate.exists():
                input_path = cwd_candidate
            # else leave as-is so the error message is clear

    if not input_path.exists():
        sys.stderr.write(
            f"[ERROR] Input TSV not found: {input_path}\n"
            f"\n"
            f"        Generate it with:\n"
            f"          python3 analysis/pyscripts/benchmark_hycal_timing.py \\\n"
            f"              <evio_path> output/per_module_025320_split00010 \\\n"
            f"              --calibrate --write-cluster-tsv \\\n"
            f"              --target-peak 2239.5 --cal-method cb-peak \\\n"
            f"              --signal-min 500 --e-range 500 3500 --bin-width 20\n"
        )
        return 1

    if args.output_pdf is None:
        stem = input_path.stem
        if stem.endswith("_clusters"):
            stem = stem[: -len("_clusters")]
        out_path = input_path.parent / (stem + "_pre_post_spectra.pdf")
    else:
        out_path = Path(args.output_pdf)
        if not out_path.is_absolute():
            # Same logic: try script dir, then CWD
            candidate = _SCRIPT_DIR / out_path
            out_path  = candidate

    e_lo, e_hi  = args.e_range
    bin_width   = args.bin_width
    target_peak = args.target_peak
    cfg_filter  = args.config

    # ---- Read TSV -----------------------------------------------------------
    print("Reading TSV...", flush=True)

    if HAVE_PANDAS:
        try:
            df = pd.read_csv(str(input_path), sep="\t", dtype=str)
        except Exception as exc:
            sys.stderr.write(f"[ERROR] Failed to read TSV: {exc}\n")
            return 1
        # Filter to chosen config
        if "config_name" in df.columns:
            df = df[df["config_name"] == cfg_filter]
        # Filter to W-prefix (PbWO4) only
        if "seed_module_id" in df.columns:
            df = df[df["seed_module_id"].str.startswith("W", na=False)]
        # Convert numeric columns
        for col in ("cluster_E_raw_MeV", "cluster_E_cal_MeV",
                    "seed_row", "seed_col"):
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
    else:
        # Pure-stdlib fallback
        sys.stderr.write("[WARN] pandas not available — using stdlib CSV reader\n")
        import csv
        rows = []
        with open(str(input_path), newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row.get("config_name") != cfg_filter:
                    continue
                mid = row.get("seed_module_id", "")
                if not mid.startswith("W"):
                    continue
                rows.append(row)
        if not rows:
            sys.stderr.write(
                f"[ERROR] No rows found for config='{cfg_filter}' / W-prefix modules\n"
            )
            return 1
        # Build a minimal dict-of-lists structure
        df = None  # signal pandas path not available; use `rows` list below

    # ---- Group by module ----------------------------------------------------
    if HAVE_PANDAS and df is not None:
        module_ids = sorted(df["seed_module_id"].unique(), key=_module_sort_key)
        # Build per-module data
        module_data: dict[str, dict] = {}
        for mid in module_ids:
            sub = df[df["seed_module_id"] == mid]
            raw_vals = sub["cluster_E_raw_MeV"].dropna().values.astype(float)
            cal_vals = sub["cluster_E_cal_MeV"].dropna().values.astype(float)
            if len(raw_vals) < args.min_events:
                continue
            # Mode of row/col
            row_mode = int(sub["seed_row"].dropna().mode().iloc[0]) \
                       if not sub["seed_row"].dropna().empty else -1
            col_mode = int(sub["seed_col"].dropna().mode().iloc[0]) \
                       if not sub["seed_col"].dropna().empty else -1
            module_data[mid] = dict(
                raw_vals=raw_vals,
                cal_vals=cal_vals,
                n_events=len(raw_vals),
                seed_row=row_mode,
                seed_col=col_mode,
            )
    else:
        # stdlib path
        from collections import defaultdict
        accum: dict[str, dict] = defaultdict(
            lambda: {"raw": [], "cal": [], "rows": [], "cols": []}
        )
        for row in rows:
            mid = row["seed_module_id"]
            try:
                raw = float(row["cluster_E_raw_MeV"])
                cal = float(row["cluster_E_cal_MeV"])
            except (ValueError, KeyError):
                continue
            accum[mid]["raw"].append(raw)
            accum[mid]["cal"].append(cal)
            try:
                accum[mid]["rows"].append(int(float(row["seed_row"])))
                accum[mid]["cols"].append(int(float(row["seed_col"])))
            except (ValueError, KeyError):
                pass

        module_data = {}
        for mid, d in accum.items():
            if len(d["raw"]) < args.min_events:
                continue
            rows_arr = d["rows"]
            cols_arr = d["cols"]

            def _mode(lst):
                if not lst:
                    return -1
                from collections import Counter
                return Counter(lst).most_common(1)[0][0]

            module_data[mid] = dict(
                raw_vals=np.array(d["raw"]),
                cal_vals=np.array(d["cal"]),
                n_events=len(d["raw"]),
                seed_row=_mode(rows_arr),
                seed_col=_mode(cols_arr),
            )

    n_total = len(module_data)
    print(f"Found {n_total} modules with ≥ {args.min_events} events; "
          f"writing PDF...", flush=True)

    if n_total == 0:
        sys.stderr.write(
            f"[ERROR] No W-prefix modules found with ≥ {args.min_events} events "
            f"for config='{cfg_filter}'.  Check the TSV and --min-events.\n"
        )
        return 1

    # ---- Sort modules -------------------------------------------------------
    sorted_mids = list(module_data.keys())

    if args.sort_by == "name":
        sorted_mids.sort(key=_module_sort_key)
    elif args.sort_by == "n_events":
        sorted_mids.sort(key=lambda m: module_data[m]["n_events"], reverse=True)
    elif args.sort_by == "mu_pre_distance":
        # Need a rough mu estimate first (just use median for sorting)
        def _pre_dist(mid):
            vals = module_data[mid]["raw_vals"]
            if len(vals) == 0:
                return 0.0
            return abs(float(np.median(vals)) - target_peak)
        sorted_mids.sort(key=_pre_dist, reverse=True)

    # Apply max_modules cap
    if args.max_modules and args.max_modules > 0:
        sorted_mids = sorted_mids[: args.max_modules]

    n_plot = len(sorted_mids)

    # ---- Pre-compute CB fits for all modules (for summary page) -------------
    # We do a first pass so summary stats are ready; the per-page pass
    # re-uses these results.
    module_results: list[dict] = []

    for mid in sorted_mids:
        d = module_data[mid]
        raw_vals = d["raw_vals"]
        cal_vals = d["cal_vals"]

        # CB fit on pre-cal values
        mu0_pre = float(np.median(raw_vals)) if len(raw_vals) > 0 else None
        rp = fit_peak_cb(raw_vals, e_lo, e_hi, bin_width, mu0_hint=mu0_pre)
        mu_pre, sig_pre = rp[0], rp[1]

        # CB fit on post-cal values
        mu0_post = float(np.median(cal_vals)) if len(cal_vals) > 0 else None
        cp = fit_peak_cb(cal_vals, e_lo, e_hi, bin_width, mu0_hint=mu0_post)
        mu_post, sig_post = cp[0], cp[1]

        module_results.append(dict(
            mid=mid,
            n_events=d["n_events"],
            seed_row=d["seed_row"],
            seed_col=d["seed_col"],
            mu_pre=mu_pre,
            sigma_pre=sig_pre,
            mu_post=mu_post,
            sigma_post=sig_post,
        ))

        mu_pre_str  = f"{mu_pre:.0f}"  if not math.isnan(mu_pre)  else "N/A"
        mu_post_str = f"{mu_post:.0f}" if not math.isnan(mu_post) else "N/A"
        print(f"  module {mid}: {d['n_events']} events, "
              f"pre μ={mu_pre_str} MeV, post μ={mu_post_str} MeV",
              flush=True)

    # ---- Open PDF and write pages -------------------------------------------
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(str(out_path)) as pdf:

        # ---- Summary page (first) -------------------------------------------
        _draw_summary_page(pdf, module_results, target_peak)

        # ---- Per-module pages -----------------------------------------------
        for res, mid in zip(module_results, sorted_mids):
            d = module_data[mid]
            raw_vals = d["raw_vals"]
            cal_vals = d["cal_vals"]
            n_ev     = d["n_events"]

            mu_pre   = res["mu_pre"]
            sig_pre  = res["sigma_pre"]
            mu_post  = res["mu_post"]
            sig_post = res["sigma_post"]

            # Detect secondary bump in POST spectrum
            has_bump = False
            if args.diagnose_secondary:
                has_bump = _secondary_bump_present(
                    cal_vals, target_peak,
                    lo_offset=200.0,
                    hi_mev=min(e_hi, 2800.0),
                )

            # Page title
            bump_star = " *** SECONDARY BUMP ***" if has_bump else ""
            title = (
                f"Module {mid}{bump_star}  —  {n_ev} events  "
                f"—  seed_row={d['seed_row']} seed_col={d['seed_col']}"
            )

            fig, (ax_pre, ax_post) = plt.subplots(1, 2, figsize=(14, 6))
            fig.suptitle(title, fontsize=10, fontweight="bold",
                         color="#cc0000" if has_bump else "black")

            # Diagnostic thresholds for text-box color
            # Pre panel: warn if >200 MeV from target
            pre_dist  = abs(mu_pre  - target_peak) if not math.isnan(mu_pre)  else 0.0
            post_dist = abs(mu_post - target_peak) if not math.isnan(mu_post) else 0.0

            pre_warn_dist  = 200.0
            post_warn_dist = 50.0

            # Draw left (pre-cal) panel
            pre_result = _draw_panel(
                ax_pre, raw_vals,
                label="Pre-calibration (cluster_E_raw_MeV)",
                e_lo=e_lo, e_hi=e_hi, bin_width=bin_width,
                target_peak=target_peak,
                col_hist="#888888", col_fit="#d62728",
                mu_warn_distance=pre_warn_dist,
                diagnose_secondary=args.diagnose_secondary,
            )

            # Draw right (post-cal) panel
            post_result = _draw_panel(
                ax_post, cal_vals,
                label="Post-calibration (cluster_E_cal_MeV)",
                e_lo=e_lo, e_hi=e_hi, bin_width=bin_width,
                target_peak=target_peak,
                col_hist="#1f77b4", col_fit="#d62728",
                mu_warn_distance=post_warn_dist,
                diagnose_secondary=args.diagnose_secondary,
            )

            # σ-worsened warning below post panel
            if (not math.isnan(sig_pre) and not math.isnan(sig_post)
                    and sig_post > sig_pre):
                delta_pct = (sig_post - sig_pre) / max(sig_pre, 1.0) * 100.0
                ax_post.text(
                    0.97, 0.55,
                    f"⚠ σ worsened by {delta_pct:.1f}%",
                    transform=ax_post.transAxes,
                    fontsize=7, va="top", ha="right", color="darkred",
                )

            # Red border on page if secondary bump detected
            if has_bump:
                for spine in ax_pre.spines.values():
                    spine.set_edgecolor("#cc0000")
                    spine.set_linewidth(2.0)
                for spine in ax_post.spines.values():
                    spine.set_edgecolor("#cc0000")
                    spine.set_linewidth(2.0)

            fig.tight_layout(rect=[0, 0, 1, 0.95])
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    n_pages = n_plot + 1  # 1 summary + n_plot module pages
    print(f"Done: PDF with {n_pages} pages (1 summary + {n_plot} modules)")
    print(f"  wrote: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
