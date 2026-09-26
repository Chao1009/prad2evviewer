#!/usr/bin/env python3
"""
plot_match_summary.py — analyze / plot output of gem_hycal_matching.py.

Reads the per-match TSV (or CSV) produced by gem_hycal_matching.py and
emits four PNG plots:

  1. {prefix}_local_hits.png   2x2 grid of 2D heatmaps — gem_x_local vs
                                gem_y_local, one panel per GEM detector.
                                Hot regions reveal acceptance / dead
                                regions on each GEM in its own frame.

  2. {prefix}_lab_scatter.png  Single scatter — gem_x vs gem_y (lab),
                                color-coded by det_id.  All four GEMs
                                overlaid in the target-centered frame.

  3. {prefix}_peak_adc.png     2x2 grid of histograms — gem_x_peak and
                                gem_y_peak (max-strip ADC per cluster),
                                X and Y overlaid per panel.

  4. {prefix}_timing.png       2x2 grid of histograms — timing of the
                                max-ADC strip in each X / Y cluster, in
                                ns (gem_*_max_tb · ts_period).

Required input columns (emitted by gem_hycal_matching.py):

  det_id, gem_x, gem_y,                             # plot 2
  gem_x_local, gem_y_local,                         # plot 1
  gem_x_peak,   gem_y_peak,                         # plot 3
  gem_x_max_tb, gem_y_max_tb                        # plot 4

Usage
-----
  # default — saves four PNGs next to the input, then pops a GUI window:
  python analysis/pyscripts/plot_match_summary.py match_023867.tsv

  # explicit out-dir, custom binning, no GUI (good for headless / CI):
  python analysis/pyscripts/plot_match_summary.py match.tsv \\
      --out-dir plots/ --bins 200 --no-show

  # other options: --csv (force CSV input), --prefix NAME (output stem),
  # --ts-period NS (time-sample period, default 25)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError as e:
    raise SystemExit(f"[ERROR] pandas required: pip install pandas ({e})")

REQUIRED_COLS = [
    "det_id", "gem_x", "gem_y",
    "gem_x_local", "gem_y_local",
    "gem_x_peak",   "gem_y_peak",
    "gem_x_max_tb", "gem_y_max_tb",
]

# Distinct color per GEM detector (0..3) — kept consistent across plots
# so panels can be cross-referenced at a glance.
DET_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]


def read_table(path: Path, force_csv: bool) -> pd.DataFrame:
    """Auto-detect delimiter from extension unless --csv overrides."""
    sep = "," if (force_csv or path.suffix.lower() == ".csv") else "\t"
    df = pd.read_csv(path, sep=sep)
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise SystemExit(
            f"[ERROR] {path}: missing required columns {missing}.\n"
            "        Re-run gem_hycal_matching.py — the local/peak/timebin "
            "columns were added 2026-04."
        )
    return df


def plot_local_hits(plt, df: pd.DataFrame, bins: int, out_path: Path) -> None:
    """2D heatmap of (gem_x_local, gem_y_local) per detector — 2×2 grid."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 9), constrained_layout=True)
    fig.suptitle("GEM matching hits — local (detector) frame")
    for d, ax in enumerate(axes.flat):
        sub = df[df.det_id == d]
        ax.set_title(f"GEM {d} — {len(sub)} hits", color=DET_COLORS[d])
        ax.set_xlabel("local x [mm]")
        ax.set_ylabel("local y [mm]")
        ax.set_aspect("equal", adjustable="box")
        if sub.empty:
            ax.text(0.5, 0.5, "no hits", ha="center", va="center",
                    transform=ax.transAxes)
            continue
        h = ax.hist2d(sub.gem_x_local, sub.gem_y_local, bins=bins, cmin=1)
        fig.colorbar(h[3], ax=ax, label="hits / bin")
    fig.savefig(out_path, dpi=150)
    print(f"  wrote {out_path}", flush=True)


def plot_lab_scatter(plt, df: pd.DataFrame, out_path: Path) -> None:
    """All matching hits in lab/target-centered frame, colored by det_id."""
    fig, ax = plt.subplots(figsize=(8, 8), constrained_layout=True)
    ax.set_title(f"GEM matching hits — lab frame (n={len(df)})")
    ax.set_xlabel("x_lab [mm]")
    ax.set_ylabel("y_lab [mm]")
    ax.set_aspect("equal", adjustable="box")
    for d in range(4):
        sub = df[df.det_id == d]
        if sub.empty:
            continue
        ax.scatter(sub.gem_x, sub.gem_y, s=2, alpha=0.4,
                   color=DET_COLORS[d], label=f"GEM {d}  ({len(sub)})")
    if len(df):
        ax.legend(markerscale=4, loc="best")
    ax.grid(True, alpha=0.3)
    fig.savefig(out_path, dpi=150)
    print(f"  wrote {out_path}", flush=True)


def _plot_xy_hists(plt, df: pd.DataFrame, xs: pd.Series, ys: pd.Series,
                   bins: int, rng: tuple[float, float], xlabel: str,
                   title: str, out_path: Path) -> None:
    """2×2 grid, one panel per detector, with the X- and Y-cluster
    quantities `xs` / `ys` (aligned with df) overlaid as step histograms."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    fig.suptitle(title)
    for d, ax in enumerate(axes.flat):
        sel = df.det_id == d
        ax.set_title(f"GEM {d}", color=DET_COLORS[d])
        ax.set_xlabel(xlabel)
        ax.set_ylabel("counts")
        if not sel.any():
            ax.text(0.5, 0.5, "no hits", ha="center", va="center",
                    transform=ax.transAxes)
            continue
        ax.hist(xs[sel], bins=bins, range=rng,
                histtype="step", linewidth=1.4, color="C0", label="X cluster")
        ax.hist(ys[sel], bins=bins, range=rng,
                histtype="step", linewidth=1.4, color="C3", label="Y cluster")
        ax.legend(loc="upper right")
    fig.savefig(out_path, dpi=150)
    print(f"  wrote {out_path}", flush=True)


def plot_peak_adc(plt, df: pd.DataFrame, bins: int, out_path: Path) -> None:
    """Histogram of x_peak / y_peak per detector — X and Y overlaid."""
    pmax = float(max(df.gem_x_peak.max(), df.gem_y_peak.max()) or 1.0)
    _plot_xy_hists(plt, df, df.gem_x_peak, df.gem_y_peak, bins,
                   (0.0, pmax * 1.05), "max strip ADC",
                   "Max-strip ADC of matching X / Y clusters", out_path)


def plot_timing(plt, df: pd.DataFrame, bins: int, ts_period: float,
                out_path: Path) -> None:
    """Histogram of max-ADC strip timing (ns) per detector."""
    tx = df.gem_x_max_tb * ts_period
    ty = df.gem_y_max_tb * ts_period
    tmin = float(min(tx.min(), ty.min()))
    tmax = float(max(tx.max(), ty.max()))
    if tmin == tmax:
        tmax = tmin + 1.0  # avoid degenerate range with a single timebin
    _plot_xy_hists(plt, df, tx, ty, bins, (tmin, tmax),
                   "max-ADC strip time [ns]",
                   f"Timing of max-ADC strip — time sample × {ts_period:.1f} ns",
                   out_path)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("input",
                    help="Path to gem_hycal_matching.py output (TSV or CSV).")
    ap.add_argument("--csv", action="store_true",
                    help="Force CSV input (default: detect by extension).")
    ap.add_argument("--out-dir", default=None,
                    help="Where to write PNGs (default: input file's dir).")
    ap.add_argument("--prefix", default=None,
                    help="Output filename prefix (default: input file stem).")
    ap.add_argument("--bins", type=int, default=120,
                    help="Bins per axis for hist2d / 1D histos (default 120).")
    ap.add_argument("--ts-period", type=float, default=25.0,
                    help="Time-sample period in ns (default 25.0).")
    ap.add_argument("--no-show", action="store_true",
                    help="Don't pop a GUI window; just save PNGs.")
    args = ap.parse_args(argv)

    in_path = Path(args.input)
    if not in_path.is_file():
        raise SystemExit(f"[ERROR] not a file: {in_path}")

    out_dir = Path(args.out_dir) if args.out_dir else in_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or in_path.stem

    # Pick backend BEFORE importing pyplot — Agg for headless, default
    # otherwise.  Importing pyplot first locks in whatever it found.
    import matplotlib
    if args.no_show:
        matplotlib.use("Agg")
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise SystemExit(
            f"[ERROR] matplotlib required: pip install matplotlib ({e})")

    print(f"[load] {in_path}", flush=True)
    df = read_table(in_path, args.csv)
    if df.empty:
        raise SystemExit("[ERROR] table has no rows — nothing to plot.")
    dets = sorted(df.det_id.unique().tolist())
    print(f"[load] {len(df)} matched rows; detectors present: {dets}",
          flush=True)

    plot_local_hits (plt, df, args.bins,
                     out_dir / f"{prefix}_local_hits.png")
    plot_lab_scatter(plt, df,
                     out_dir / f"{prefix}_lab_scatter.png")
    plot_peak_adc   (plt, df, args.bins,
                     out_dir / f"{prefix}_peak_adc.png")
    plot_timing     (plt, df, args.bins, args.ts_period,
                     out_dir / f"{prefix}_timing.png")

    if not args.no_show:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
