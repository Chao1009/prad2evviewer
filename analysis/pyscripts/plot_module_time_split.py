#!/usr/bin/env python3
"""plot_module_time_split.py — PMT gain-drift diagnostic via time-bin spectra.

Diagnoses whether two prominent peaks in a module's pre-calibration spectrum
are caused by PMT gain drift (changing over time) or by static detector
pathology.

For each requested module the script splits its events into N equal-count time
bins (ordered by event_num) and plots all per-bin pre-calibration energy
histograms overlaid on a single panel.  A per-module drift metric (spread of
per-bin histogram-argmax peaks) is computed and flagged when > 2 %.

Input: a clusters TSV produced by
    benchmark_hycal_timing.py --write-cluster-tsv
  Columns:
    event_num  config_name  seed_module_id  seed_x_mm  seed_y_mm
    seed_row   seed_col     seed_energy_MeV  cluster_E_raw_MeV
    cluster_E_cal_MeV

Output: a single multi-page PDF (one page per requested module).

Usage:
  python3 plot_module_time_split.py --modules W559 W593 [options]
  python3 plot_module_time_split.py --help
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages

try:
    import pandas as pd
    HAVE_PANDAS = True
except ImportError:
    HAVE_PANDAS = False

# ---------------------------------------------------------------------------
# Repository layout helpers
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT   = _SCRIPT_DIR.parent.parent

# Minimum events per module to generate a page.
_MIN_EVENTS = 50

# Drift threshold (percent) above which the text-box is flagged pale-red.
_DRIFT_FLAG_PCT = 2.0


# ---------------------------------------------------------------------------
# CLI
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
        help="Output PDF path.  Default: <input-stem>_time_split.pdf in same "
             "directory as input.",
    )
    ap.add_argument(
        "--modules",
        nargs="+",
        required=True,
        metavar="MODULE",
        help="Space-separated list of module names to plot (e.g., W559 W593). "
             "Required.",
    )
    ap.add_argument(
        "--n-bins",
        type=int,
        default=3,
        help="Number of equal-count time bins (default: 3)",
    )
    ap.add_argument(
        "--e-range",
        type=float,
        nargs=2,
        default=(500.0, 3500.0),
        metavar=("LO", "HI"),
        help="Energy histogram range in MeV (default: 500 3500)",
    )
    ap.add_argument(
        "--bin-width",
        type=float,
        default=20.0,
        help="Histogram bin width in MeV (default: 20)",
    )
    ap.add_argument(
        "--target-peak",
        type=float,
        default=2239.5,
        help="Target peak for vertical reference dashed line (default: 2239.5 MeV)",
    )
    ap.add_argument(
        "--config",
        choices=["legacy", "new"],
        default="legacy",
        help="Which config_name to filter to: 'legacy' or 'new' (default: legacy)",
    )
    return ap


# ---------------------------------------------------------------------------
# Path resolution (mirrors plot_per_module_pre_post.py convention)
# ---------------------------------------------------------------------------

def _resolve_input(raw: str) -> Path:
    """Return the best-match absolute path for a possibly-relative TSV path."""
    p = Path(raw)
    if p.is_absolute():
        return p
    # Try relative to CWD first (user running from repo root), then script dir
    for candidate in (Path.cwd() / p, _SCRIPT_DIR / p):
        if candidate.exists():
            return candidate
    # Fall back to CWD-relative for a clean error message
    return Path.cwd() / p


def _default_output(input_path: Path) -> Path:
    stem = input_path.stem
    if stem.endswith("_clusters"):
        stem = stem[: -len("_clusters")]
    return input_path.parent / (stem + "_time_split.pdf")


# ---------------------------------------------------------------------------
# Color selection
# ---------------------------------------------------------------------------

def _bin_colors(n: int) -> list:
    """Return N distinct colors for time bins."""
    if n == 3:
        return ["tab:blue", "tab:green", "tab:red"]
    if n == 2:
        return ["tab:blue", "tab:red"]
    if n == 1:
        return ["tab:blue"]
    # General case: viridis colormap
    return [cm.viridis(v) for v in np.linspace(0, 1, n)]


# ---------------------------------------------------------------------------
# Per-module processing
# ---------------------------------------------------------------------------

def _bin_edges_equal_count(event_nums: np.ndarray, n_bins: int) -> List[int]:
    """Return (n_bins+1) indices into sorted event_nums that split the array
    into equal-count groups.  Index 0 is 0, index n_bins is len(event_nums).
    """
    n = len(event_nums)
    edges = [int(round(i * n / n_bins)) for i in range(n_bins + 1)]
    # Ensure strictly increasing and bounded
    edges[0]  = 0
    edges[-1] = n
    return edges


def _histogram_argmax_energy(counts: np.ndarray, edges: np.ndarray) -> float:
    """Return the bin-centre energy at which counts is maximum."""
    if counts.sum() == 0:
        return float("nan")
    centres = 0.5 * (edges[:-1] + edges[1:])
    return float(centres[int(np.argmax(counts))])


def process_module(
    mod_id: str,
    df_mod,           # pandas DataFrame subset, already sorted by event_num
    n_bins: int,
    e_lo: float,
    e_hi: float,
    bin_width: float,
) -> Optional[dict]:
    """
    Split the module's events into time bins and compute per-bin histograms.

    Returns a dict with keys:
      mod_id, n_events, bin_data (list of dicts per bin), drift_pct,
      mu_list, ev_edges
    or None if < _MIN_EVENTS events.
    """
    n_events = len(df_mod)
    if n_events < _MIN_EVENTS:
        return None

    # Sort by event_num (should already be sorted, but be safe)
    df_mod = df_mod.sort_values("event_num")
    event_nums  = df_mod["event_num"].values.astype(int)
    raw_energies = df_mod["cluster_E_raw_MeV"].values.astype(float)

    # Build histogram edges for the energy axis (shared across bins)
    energy_edges = np.arange(e_lo, e_hi + bin_width, bin_width)

    # Equal-count split indices
    idx_edges = _bin_edges_equal_count(event_nums, n_bins)

    bin_data = []
    for k in range(n_bins):
        i0 = idx_edges[k]
        i1 = idx_edges[k + 1]
        ev_lo = int(event_nums[i0])
        ev_hi = int(event_nums[i1 - 1])  # last event_num in this bin
        vals  = raw_energies[i0:i1]
        counts, _ = np.histogram(vals, bins=energy_edges)
        mu = _histogram_argmax_energy(counts, energy_edges)
        bin_data.append(dict(
            k=k + 1,              # 1-based
            i0=i0, i1=i1,
            ev_lo=ev_lo, ev_hi=ev_hi,
            n=i1 - i0,
            counts=counts,
            mu=mu,
        ))

    # Drift metric
    mu_list = [b["mu"] for b in bin_data if not np.isnan(b["mu"])]
    if len(mu_list) >= 2:
        mu_arr  = np.array(mu_list)
        mu_mean = float(np.mean(mu_arr))
        mu_max  = float(np.max(mu_arr))
        mu_min  = float(np.min(mu_arr))
        drift_pct = (mu_max - mu_min) / mu_mean * 100.0 if mu_mean > 0 else 0.0
    else:
        mu_mean = float(mu_list[0]) if mu_list else float("nan")
        mu_max  = mu_mean
        mu_min  = mu_mean
        drift_pct = 0.0

    return dict(
        mod_id=mod_id,
        n_events=n_events,
        bin_data=bin_data,
        energy_edges=energy_edges,
        drift_pct=drift_pct,
        mu_list=mu_list,
        mu_max=mu_max,
        mu_min=mu_min,
        mu_mean=mu_mean,
    )


# ---------------------------------------------------------------------------
# Progress reporting
# ---------------------------------------------------------------------------

def _print_module_summary(result: dict) -> None:
    mod_id    = result["mod_id"]
    n_events  = result["n_events"]
    n_bins    = len(result["bin_data"])
    drift_pct = result["drift_pct"]
    mu_max    = result["mu_max"]
    mu_min    = result["mu_min"]

    print(
        f"Processing module {mod_id}: {n_events} events, {n_bins} time bins",
        flush=True,
    )
    for b in result["bin_data"]:
        mu_str = f"{b['mu']:.0f}" if not np.isnan(b["mu"]) else "N/A"
        print(
            f"  bin {b['k']}: ev [{b['ev_lo']}, {b['ev_hi']}]"
            f"   N={b['n']}   μ={mu_str} MeV",
            flush=True,
        )
    if mu_min > 0:
        ratio = mu_max / mu_min
        sign  = "+" if drift_pct >= 0 else ""
        print(
            f"  drift: {sign}{drift_pct:.1f}% (max/min = {ratio:.3f})",
            flush=True,
        )
    else:
        print(f"  drift: N/A", flush=True)


# ---------------------------------------------------------------------------
# Per-module page rendering
# ---------------------------------------------------------------------------

def _draw_module_page(pdf: "PdfPages", result: dict, target_peak: float) -> None:
    mod_id     = result["mod_id"]
    n_events   = result["n_events"]
    bin_data   = result["bin_data"]
    energy_edges = result["energy_edges"]
    drift_pct  = result["drift_pct"]
    mu_list    = result["mu_list"]
    mu_max     = result["mu_max"]
    mu_min     = result["mu_min"]
    mu_mean    = result["mu_mean"]
    n_bins     = len(bin_data)

    colors = _bin_colors(n_bins)
    centres = 0.5 * (energy_edges[:-1] + energy_edges[1:])

    fig, ax = plt.subplots(figsize=(11, 7))
    fig.suptitle(
        f"Module {mod_id} — pre-calibration spectrum, split by event_num",
        fontsize=11, fontweight="bold",
    )

    # Draw one step histogram per time bin
    for b, color in zip(bin_data, colors):
        label = (
            f"Bin {b['k']}: events {b['ev_lo']}–{b['ev_hi']}"
            f"  N={b['n']}"
        )
        ax.step(centres, b["counts"], where="mid",
                color=color, lw=1.5, label=label, alpha=0.85)

    # Reference line
    ax.axvline(target_peak, color="#444", lw=1.0, ls="--",
               label=f"target {target_peak:.1f} MeV")

    ax.set_xlabel("cluster_E_raw_MeV (MeV)", fontsize=9)
    ax.set_ylabel("counts / bin", fontsize=9)
    ax.set_xlim(energy_edges[0], energy_edges[-1])
    ax.grid(alpha=0.25)
    ax.legend(loc="upper left", fontsize=8)

    # --- Text summary box (upper right) ---
    lines = [f"Module {mod_id}", f"Total events: {n_events}"]

    # Per-bin μ
    for b in bin_data:
        mu_str = f"{b['mu']:.0f} MeV" if not np.isnan(b["mu"]) else "N/A"
        lines.append(f"μ_bin{b['k']} = {mu_str}")

    # Overall drift
    if not np.isnan(mu_mean) and mu_mean > 0:
        sign = "+" if drift_pct >= 0 else ""
        if len(mu_list) >= 2 and mu_min > 0:
            ratio = mu_max / mu_min
            lines.append(f"drift_pct = {sign}{drift_pct:.1f}%")
            lines.append(f"μ_max/μ_min = {ratio:.3f}")
        else:
            lines.append("drift_pct = N/A")
    else:
        lines.append("drift_pct = N/A")

    txt = "\n".join(lines)
    box_color = "mistyrose" if drift_pct > _DRIFT_FLAG_PCT else "lightyellow"
    props = dict(boxstyle="round", facecolor=box_color, alpha=0.90)
    ax.text(
        0.97, 0.97, txt,
        transform=ax.transAxes,
        fontsize=8, va="top", ha="right",
        bbox=props,
        fontfamily="monospace",
    )

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    ap   = build_parser()
    args = ap.parse_args(argv)

    # ---- Resolve paths -------------------------------------------------------
    input_path = _resolve_input(args.input_tsv)

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
        out_path = _default_output(input_path)
    else:
        out_path = Path(args.output_pdf)
        if not out_path.is_absolute():
            out_path = Path.cwd() / out_path

    e_lo, e_hi  = args.e_range
    bin_width   = args.bin_width
    target_peak = args.target_peak
    n_bins      = args.n_bins
    cfg_filter  = args.config
    requested   = list(args.modules)

    if n_bins < 1:
        sys.stderr.write("[ERROR] --n-bins must be >= 1\n")
        return 1

    # ---- Read TSV ------------------------------------------------------------
    print("Reading TSV...", flush=True)

    if not HAVE_PANDAS:
        sys.stderr.write(
            "[ERROR] pandas is required for this script.  "
            "Install with: pip install pandas\n"
        )
        return 1

    try:
        df = pd.read_csv(str(input_path), sep="\t", dtype=str)
    except Exception as exc:
        sys.stderr.write(f"[ERROR] Failed to read TSV: {exc}\n")
        return 1

    if df.empty:
        sys.stderr.write(
            f"[ERROR] TSV is empty: {input_path}\n"
        )
        return 1

    # Validate required columns
    required_cols = {
        "event_num", "config_name", "seed_module_id",
        "cluster_E_raw_MeV",
    }
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        sys.stderr.write(
            f"[ERROR] TSV is missing required columns: "
            f"{', '.join(sorted(missing_cols))}\n"
            f"        Found columns: {', '.join(df.columns)}\n"
        )
        return 1

    # Filter to chosen config
    df = df[df["config_name"] == cfg_filter].copy()
    if df.empty:
        sys.stderr.write(
            f"[ERROR] No rows found for config_name='{cfg_filter}'.  "
            f"Check --config.\n"
        )
        return 1

    # Convert numeric columns
    df["event_num"]        = pd.to_numeric(df["event_num"],        errors="coerce")
    df["cluster_E_raw_MeV"] = pd.to_numeric(df["cluster_E_raw_MeV"], errors="coerce")
    df = df.dropna(subset=["event_num", "cluster_E_raw_MeV"])

    # ---- Check which requested modules exist in TSV --------------------------
    available = set(df["seed_module_id"].unique())
    not_in_tsv = [m for m in requested if m not in available]
    if not_in_tsv:
        for m in not_in_tsv:
            sys.stderr.write(
                f"[WARN] Module '{m}' not found in TSV (config='{cfg_filter}') "
                f"— skipping.\n"
            )

    to_process = [m for m in requested if m in available]
    if not to_process:
        sys.stderr.write(
            "[ERROR] None of the requested modules were found in the TSV.  "
            "Exiting.\n"
        )
        return 1

    # ---- Process each module -------------------------------------------------
    results = []
    for mod_id in to_process:
        df_mod = df[df["seed_module_id"] == mod_id].copy()
        n_ev   = len(df_mod)

        if n_ev < _MIN_EVENTS:
            sys.stderr.write(
                f"[WARN] Module '{mod_id}' has only {n_ev} events "
                f"(minimum is {_MIN_EVENTS}) — skipping.\n"
            )
            continue

        result = process_module(
            mod_id=mod_id,
            df_mod=df_mod,
            n_bins=n_bins,
            e_lo=e_lo,
            e_hi=e_hi,
            bin_width=bin_width,
        )
        if result is None:
            # Shouldn't happen given the guard above, but be safe
            sys.stderr.write(
                f"[WARN] Module '{mod_id}' produced no result — skipping.\n"
            )
            continue

        _print_module_summary(result)
        results.append(result)

    if not results:
        sys.stderr.write(
            "[ERROR] No modules with sufficient events to plot.  "
            "Check --modules and the TSV content.\n"
        )
        return 1

    # ---- Write PDF -----------------------------------------------------------
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(str(out_path)) as pdf:
        for result in results:
            _draw_module_page(pdf, result, target_peak=target_peak)

    n_pages = len(results)
    print(f"Done: PDF with {n_pages} page{'s' if n_pages != 1 else ''}", flush=True)
    print(f"  wrote: {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
