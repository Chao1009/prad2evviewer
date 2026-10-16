#!/usr/bin/env python3
"""plot_cluster_energy_vs_position.py — diagnostic 4-panel figure from cluster TSV.

Reads the per-event cluster TSV written by benchmark_hycal_timing.py
--write-cluster-tsv and produces a 2×2 figure:

  Panel A  Cluster energy vs seed x position (2D hist)
  Panel B  Cluster energy vs seed y position (2D hist)
  Panel C  Cluster energy vs seed module ID   (2D hist with PbGlass/PbWO4 boundary)
  Panel D  Mean cluster energy mapped onto (x, y) seed position (plasma)

Useful for localising double-peak structure or per-module gain offsets in
the recalibrated spectrum.  No prad2py, no pandas required.

Optional --find-multiple-peaks flag runs a per-module peak-detection pass
and emits a CSV + PDF identifying modules with multi-peak energy distributions
(diagnostic for per-module calibration residuals).

Usage
-----
  python analysis/pyscripts/plot_cluster_energy_vs_position.py \\
      output/run24386_clusters.tsv output/run24386_cluster_diag.png

  python analysis/pyscripts/plot_cluster_energy_vs_position.py \\
      output/run24386_clusters.tsv output/run24386_cluster_diag.png \\
      --config legacy --energy calibrated --e-range 500 3500 --e-bins 120

  python analysis/pyscripts/plot_cluster_energy_vs_position.py \\
      output/run24386_clusters.tsv output/run24386_cluster_diag.png \\
      --find-multiple-peaks --min-events-per-module 50 \\
      --min-peak-separation 100 --peak-prominence 5
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Soft-require scipy
try:
    from scipy.signal import find_peaks as _scipy_find_peaks
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False


# ---------------------------------------------------------------------------
# TSV reader (no pandas)
# ---------------------------------------------------------------------------

def _load_tsv(path: Path, config: str, energy_col: str
              ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Parse the cluster TSV.  Returns (x_mm, y_mm, module_id_int, seed_E, clust_E)
    as float64 arrays filtered to `config` and with finite `energy_col` values.

    `config`     one of 'legacy', 'gated'
    `energy_col` one of 'raw', 'calibrated'
    """
    col_e = "cluster_E_cal_MeV" if energy_col == "calibrated" else "cluster_E_raw_MeV"

    x_list, y_list, sid_list, se_list, ce_list = [], [], [], [], []
    with open(path, encoding="utf-8") as fh:
        header = None
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                try:
                    i_cfg  = header.index("config_name")
                    i_sid  = header.index("seed_module_id")
                    i_x    = header.index("seed_x_mm")
                    i_y    = header.index("seed_y_mm")
                    i_se   = header.index("seed_energy_MeV")
                    i_ce   = header.index(col_e)
                except ValueError as exc:
                    raise SystemExit(
                        f"[ERROR] {path}: missing required column — {exc}\n"
                        "        Was the TSV written by --write-cluster-tsv?"
                    ) from exc
                continue
            if len(fields) <= max(i_cfg, i_sid, i_x, i_y, i_se, i_ce):
                continue
            if fields[i_cfg] != config:
                continue
            try:
                ce = float(fields[i_ce])
            except ValueError:
                continue
            if not (ce == ce):  # NaN check without math.isnan
                continue
            # Decode seed_module_id ("G<n>" → n, "W<n>" → n+1000)
            sid_str = fields[i_sid]
            try:
                if sid_str.startswith("G"):
                    sid_int = int(sid_str[1:])
                elif sid_str.startswith("W"):
                    sid_int = int(sid_str[1:]) + 1000
                else:
                    sid_int = int(sid_str)
            except ValueError:
                sid_int = 0
            x_list.append(float(fields[i_x]))
            y_list.append(float(fields[i_y]))
            sid_list.append(sid_int)
            se_list.append(float(fields[i_se]))
            ce_list.append(ce)

    if not x_list:
        raise SystemExit(
            f"[ERROR] no rows found for config='{config}' with finite "
            f"'{col_e}' in {path}."
        )
    return (np.array(x_list, dtype=float),
            np.array(y_list, dtype=float),
            np.array(sid_list, dtype=float),
            np.array(se_list, dtype=float),
            np.array(ce_list, dtype=float))


# ---------------------------------------------------------------------------
# Figure builder
# ---------------------------------------------------------------------------

def make_figure(x_mm, y_mm, sid, se, ce,
                config: str, energy_label: str,
                e_range: tuple[float, float], e_bins: int,
                xy_bins: int, mod_bins: int, dpi: int) -> plt.Figure:
    """Build and return the 2×2 diagnostic figure."""
    e_lo, e_hi = e_range

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        f"Cluster energy diagnostics  (config={config}, energy={energy_label})",
        fontsize=12
    )

    # ---- Panel A: E vs seed x -----------------------------------------------
    ax = axes[0, 0]
    h2 = ax.hist2d(x_mm, ce,
                   bins=[xy_bins, e_bins],
                   range=[[x_mm.min(), x_mm.max()], [e_lo, e_hi]],
                   cmap="viridis")
    fig.colorbar(h2[3], ax=ax, label="events")
    ax.set_xlabel("seed x (mm)")
    ax.set_ylabel("cluster energy (MeV)")
    ax.set_title("Cluster energy vs seed x position")
    ax.grid(alpha=0.25)

    # ---- Panel B: E vs seed y -----------------------------------------------
    ax = axes[0, 1]
    h2 = ax.hist2d(y_mm, ce,
                   bins=[xy_bins, e_bins],
                   range=[[y_mm.min(), y_mm.max()], [e_lo, e_hi]],
                   cmap="viridis")
    fig.colorbar(h2[3], ax=ax, label="events")
    ax.set_xlabel("seed y (mm)")
    ax.set_ylabel("cluster energy (MeV)")
    ax.set_title("Cluster energy vs seed y position")
    ax.grid(alpha=0.25)

    # ---- Panel C: E vs seed module ID ---------------------------------------
    ax = axes[1, 0]
    sid_lo = float(sid.min())
    sid_hi = float(sid.max())
    h2 = ax.hist2d(sid, ce,
                   bins=[mod_bins, e_bins],
                   range=[[sid_lo, sid_hi], [e_lo, e_hi]],
                   cmap="viridis")
    fig.colorbar(h2[3], ax=ax, label="events")
    # Boundary lines at G/W module ID boundary
    for vx, vlabel in [(900, "G900"), (1001, "W1")]:
        if sid_lo <= vx <= sid_hi:
            ax.axvline(vx, color="#d62728", lw=0.9, ls="--", alpha=0.8,
                       label=f"{vlabel} boundary")
    ax.set_xlabel("seed module ID (PrimEx integer)")
    ax.set_ylabel("cluster energy (MeV)")
    ax.set_title("Cluster energy vs seed module ID")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.25, axis="y")

    # ---- Panel D: mean-E spatial map ----------------------------------------
    ax = axes[1, 1]
    x_lo, x_hi = float(x_mm.min()), float(x_mm.max())
    y_lo, y_hi = float(y_mm.min()), float(y_mm.max())

    # Use numpy to build sum_E and count on a (xy_bins × xy_bins) grid.
    x_edges = np.linspace(x_lo, x_hi, xy_bins + 1)
    y_edges = np.linspace(y_lo, y_hi, xy_bins + 1)

    # np.histogram2d convention: first axis = x (rows), second = y (cols).
    sum_E,   _, _ = np.histogram2d(x_mm, y_mm, bins=[x_edges, y_edges],
                                   weights=ce)
    count_E, _, _ = np.histogram2d(x_mm, y_mm, bins=[x_edges, y_edges])

    with np.errstate(invalid="ignore", divide="ignore"):
        mean_E = np.where(count_E > 0, sum_E / count_E, np.nan)

    # pcolormesh expects (ny, nx) with x along columns → transpose
    mean_E_T = mean_E.T

    pm = ax.pcolormesh(x_edges, y_edges, mean_E_T,
                       cmap="plasma", vmin=e_lo, vmax=e_hi)
    fig.colorbar(pm, ax=ax, label="mean cluster energy (MeV)")

    # Semi-transparent mask where fewer than 5 events fall in the cell
    low_stats = (count_E.T < 5)
    mask_rgba = np.zeros((*low_stats.shape, 4), dtype=float)
    mask_rgba[..., 0] = 0.5   # gray
    mask_rgba[..., 1] = 0.5
    mask_rgba[..., 2] = 0.5
    mask_rgba[..., 3] = np.where(low_stats, 0.55, 0.0)
    ax.pcolormesh(x_edges, y_edges, mask_rgba,
                  shading="flat")

    ax.set_xlabel("seed x (mm)")
    ax.set_ylabel("seed y (mm)")
    ax.set_title("Mean cluster energy by seed position")
    ax.set_aspect("equal", adjustable="box")

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Multi-peak detection
# ---------------------------------------------------------------------------

@dataclass
class ModulePeakInfo:
    module_id: int
    module_name: str
    is_pbwo4: bool
    seed_x_mm: float
    seed_y_mm: float
    n_events: int
    n_peaks: int
    peak_positions_MeV: List[float]
    peak_heights_counts: List[int]
    peak_prominences: List[float]
    secondary_peak_ratio: float


def _mod_name(mod_id: int) -> str:
    """Convert integer module ID to 'G<n>' or 'W<n>' name."""
    if mod_id <= 900:
        return f"G{mod_id}"
    return f"W{mod_id - 1000}"


@dataclass
class _PeakCtx:
    min_events: int
    min_peak_separation: float   # MeV
    peak_prominence: float       # counts floor
    e_range: tuple[float, float]
    e_bins: int


def find_multi_peak_modules(x_mm: np.ndarray, y_mm: np.ndarray,
                            sid: np.ndarray, ce: np.ndarray,
                            ctx: _PeakCtx) -> List[ModulePeakInfo]:
    """Identify HyCal modules whose cluster-energy distribution has ≥2 peaks.

    Parameters
    ----------
    x_mm, y_mm : seed position arrays (float, per-event)
    sid        : seed module ID integer array (per-event)
    ce         : cluster energy MeV array (per-event)
    ctx        : detection thresholds

    Returns list of ModulePeakInfo, one entry per flagged module.
    """
    e_lo, e_hi = ctx.e_range
    bin_width = (e_hi - e_lo) / ctx.e_bins
    bin_edges = np.linspace(e_lo, e_hi, ctx.e_bins + 1)
    bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    min_sep_bins = max(1, int(ctx.min_peak_separation / bin_width))

    # Outer-20% boundary indices (first 10% + last 10%)
    n10 = max(1, ctx.e_bins // 10)

    results: List[ModulePeakInfo] = []

    unique_ids = np.unique(sid.astype(int))
    for mod_id in unique_ids:
        mask = (sid.astype(int) == mod_id)
        e_mod = ce[mask]
        n = len(e_mod)
        if n < ctx.min_events:
            continue

        counts, _ = np.histogram(e_mod, bins=bin_edges)

        # Local background: linear fit to outer 20% bins
        outer_idx = np.concatenate([np.arange(n10), np.arange(ctx.e_bins - n10, ctx.e_bins)])
        outer_x = bin_centres[outer_idx].astype(float)
        outer_y = counts[outer_idx].astype(float)
        try:
            coeffs = np.polyfit(outer_x, outer_y, 1)
            bg_counts = np.polyval(coeffs, bin_centres)
            local_B = float(np.median(np.clip(bg_counts, 0, None)))
        except (np.linalg.LinAlgError, ValueError):
            local_B = float(np.median(counts))
        local_B = max(local_B, 1.0)  # avoid sqrt(0)

        prominence_thresh = max(ctx.peak_prominence, 3.0 * np.sqrt(local_B))

        peak_indices, props = _scipy_find_peaks(
            counts,
            prominence=prominence_thresh,
            distance=min_sep_bins,
        )

        if len(peak_indices) < 2:
            continue

        # Build per-peak info
        peak_positions = [float(bin_centres[i]) for i in peak_indices]
        peak_heights   = [int(counts[i]) for i in peak_indices]
        peak_proms     = [float(p) for p in props["prominences"]]

        # secondary_peak_ratio: position ratio of second-tallest / tallest
        sorted_by_height = sorted(
            zip(peak_heights, peak_positions, peak_proms),
            key=lambda t: t[0], reverse=True
        )
        primary_pos   = sorted_by_height[0][1]
        secondary_pos = sorted_by_height[1][1]
        ratio = secondary_pos / primary_pos if primary_pos != 0.0 else 0.0

        name = _mod_name(mod_id)
        results.append(ModulePeakInfo(
            module_id=mod_id,
            module_name=name,
            is_pbwo4=name.startswith("W"),
            seed_x_mm=float(np.mean(x_mm[mask])),
            seed_y_mm=float(np.mean(y_mm[mask])),
            n_events=n,
            n_peaks=len(peak_indices),
            peak_positions_MeV=peak_positions,
            peak_heights_counts=peak_heights,
            peak_prominences=peak_proms,
            secondary_peak_ratio=ratio,
        ))

    results.sort(key=lambda r: r.secondary_peak_ratio, reverse=True)
    return results


# ---------------------------------------------------------------------------
# Output helpers for multi-peak results
# ---------------------------------------------------------------------------

def _write_multi_peak_csv(results: List[ModulePeakInfo], csv_path: Path) -> None:
    """Write multi-peak module list to CSV."""
    fieldnames = [
        "module_id", "module_name", "is_pbwo4",
        "seed_x_mm", "seed_y_mm",
        "n_events", "n_peaks",
        "peak_positions_MeV", "peak_heights_counts", "peak_prominences",
        "secondary_peak_ratio",
    ]
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow({
                "module_id": r.module_id,
                "module_name": r.module_name,
                "is_pbwo4": r.is_pbwo4,
                "seed_x_mm": f"{r.seed_x_mm:.3f}",
                "seed_y_mm": f"{r.seed_y_mm:.3f}",
                "n_events": r.n_events,
                "n_peaks": r.n_peaks,
                "peak_positions_MeV": ";".join(f"{p:.1f}" for p in r.peak_positions_MeV),
                "peak_heights_counts": ";".join(str(h) for h in r.peak_heights_counts),
                "peak_prominences": ";".join(f"{p:.2f}" for p in r.peak_prominences),
                "secondary_peak_ratio": f"{r.secondary_peak_ratio:.4f}",
            })


def _write_multi_peak_pdf(results: List[ModulePeakInfo], pdf_path: Path,
                          ce_all: np.ndarray, sid_all: np.ndarray,
                          ctx: _PeakCtx) -> None:
    """Write one-page-per-module PDF for flagged modules."""
    e_lo, e_hi = ctx.e_range
    bin_edges = np.linspace(e_lo, e_hi, ctx.e_bins + 1)
    bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # Colour scheme: primary=red, secondary=green, others=grey
    _colours = ["#d62728", "#2ca02c", "#7f7f7f", "#9467bd", "#8c564b"]

    with PdfPages(pdf_path) as pdf:
        for r in results:
            mask = (sid_all.astype(int) == r.module_id)
            e_mod = ce_all[mask]
            counts, _ = np.histogram(e_mod, bins=bin_edges)

            fig, ax = plt.subplots(figsize=(9, 5))
            ax.bar(bin_centres, counts,
                   width=(e_hi - e_lo) / ctx.e_bins * 0.9,
                   color="#1f77b4", alpha=0.7, label="cluster E")

            # Sort peaks by height to assign primary/secondary colours
            order = sorted(range(r.n_peaks),
                           key=lambda i: r.peak_heights_counts[i],
                           reverse=True)
            for rank, idx in enumerate(order):
                colour = _colours[min(rank, len(_colours) - 1)]
                label = ("primary" if rank == 0
                         else "secondary" if rank == 1
                         else f"peak {rank+1}")
                ax.axvline(r.peak_positions_MeV[idx], color=colour,
                           ls="--", lw=1.5, label=f"{label} {r.peak_positions_MeV[idx]:.1f} MeV")

            # Annotation box
            ann_lines = [f"peaks: {r.n_peaks}"]
            for i in range(r.n_peaks):
                ann_lines.append(
                    f"  pos={r.peak_positions_MeV[i]:.1f} MeV  "
                    f"h={r.peak_heights_counts[i]}  "
                    f"prom={r.peak_prominences[i]:.1f}"
                )
            ann_lines.append(f"ratio={r.secondary_peak_ratio:.4f}")
            ax.text(0.97, 0.97, "\n".join(ann_lines),
                    transform=ax.transAxes,
                    va="top", ha="right", fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.4", fc="wheat", alpha=0.8))

            ax.set_xlabel("cluster energy (MeV)")
            ax.set_ylabel("counts / bin")
            ax.set_xlim(e_lo, e_hi)
            ax.set_title(
                f"Module {r.module_name} (id {r.module_id})  "
                f"n_events={r.n_events}  "
                f"position=({r.seed_x_mm:.1f}, {r.seed_y_mm:.1f}) mm"
            )
            ax.legend(fontsize=8, loc="upper left")
            ax.grid(alpha=0.25)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def _print_multi_peak_summary(results: List[ModulePeakInfo],
                              n_tested: int, ctx: _PeakCtx,
                              config: str, energy_label: str,
                              csv_path: Path, pdf_path: Path | None) -> None:
    """Print multi-peak detection summary to stdout."""
    n_flagged = len(results)
    n_pbwo4   = sum(1 for r in results if r.is_pbwo4)
    n_pbglass = n_flagged - n_pbwo4

    sep = "=" * 72
    print(sep)
    print(f"Multi-peak module detection ({config}, {energy_label})")
    print(sep)
    print(f"  modules tested         : {n_tested}  (with >= {ctx.min_events} events)")
    print(f"  modules with ≥2 peaks  : {n_flagged}")
    print(f"  of which PbWO4         : {n_pbwo4}")
    print(f"  of which PbGlass       : {n_pbglass}")
    print()

    if n_flagged == 0:
        print("  No modules with multiple peaks detected.")
    else:
        top_n = min(10, n_flagged)
        print(f"Top {top_n} by severity (secondary / primary peak ratio):")
        hdr = f"  {'module':<8}  {'n_events':>8}  {'n_peaks':>7}  " \
              f"{'primary (MeV)':>13}  {'secondary (MeV)':>15}  {'ratio':>6}"
        print(hdr)
        for r in results[:top_n]:
            # Primary = tallest peak
            order = sorted(range(r.n_peaks),
                           key=lambda i: r.peak_heights_counts[i],
                           reverse=True)
            primary_pos   = r.peak_positions_MeV[order[0]]
            secondary_pos = r.peak_positions_MeV[order[1]]
            print(f"  {r.module_name:<8}  {r.n_events:>8}  {r.n_peaks:>7}  "
                  f"{primary_pos:>13.1f}  {secondary_pos:>15.1f}  "
                  f"{r.secondary_peak_ratio:>6.3f}")
    print()
    print(f"  wrote: {csv_path}")
    if pdf_path is not None:
        print(f"  wrote: {pdf_path}")


def _run_multi_peak(x_mm, y_mm, sid, ce,
                    ctx: _PeakCtx,
                    config: str, energy_label: str,
                    out_dir: Path, fig_stem: str, suffix: str) -> None:
    """Run detection + emit CSV / PDF + print summary for one (config, energy) combo."""
    csv_path = out_dir / f"{fig_stem}_multi_peak_modules{suffix}.csv"
    pdf_path = out_dir / f"{fig_stem}_multi_peak_modules{suffix}.pdf"

    n_tested_ids = np.unique(sid.astype(int))
    # Count how many pass the min-events filter (for summary)
    n_tested = sum(
        1 for mid in n_tested_ids
        if np.sum(sid.astype(int) == mid) >= ctx.min_events
    )

    print(f"[peaks] scanning {len(n_tested_ids)} modules "
          f"({n_tested} with >= {ctx.min_events} events) …", flush=True)

    results = find_multi_peak_modules(x_mm, y_mm, sid, ce, ctx)

    _write_multi_peak_csv(results, csv_path)

    if results:
        _write_multi_peak_pdf(results, pdf_path, ce, sid, ctx)
        pdf_out = pdf_path
    else:
        pdf_out = None

    _print_multi_peak_summary(results, n_tested, ctx,
                              config, energy_label, csv_path, pdf_out)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("input_tsv",
                    help="Cluster TSV from benchmark_hycal_timing.py "
                         "--write-cluster-tsv.")
    ap.add_argument("output_png",
                    help="Output figure path (PNG); a .pdf sibling is also "
                         "written automatically.")
    ap.add_argument("--config", choices=["legacy", "gated", "both"],
                    default="legacy",
                    help="Which config to plot (default: legacy).  "
                         "'both' writes two output files with "
                         "_legacy / _gated suffixes.")
    ap.add_argument("--energy", choices=["raw", "calibrated", "both"],
                    default="calibrated",
                    help="Which energy column to plot (default: calibrated).  "
                         "'both' adds _raw / _calibrated suffixes.")
    ap.add_argument("--e-range", type=float, nargs=2, default=(500.0, 3500.0),
                    metavar=("LO", "HI"),
                    help="MeV range for energy axis (default 500 3500).")
    ap.add_argument("--e-bins", type=int, default=120,
                    help="Number of energy bins (default 120).")
    ap.add_argument("--xy-bins", type=int, default=60,
                    help="Number of x and y bins (default 60).")
    ap.add_argument("--module-bins", type=int, default=100,
                    help="Number of seed-module bins for panel C (default 100).")
    ap.add_argument("--dpi", type=int, default=150,
                    help="Output DPI (default 150).")
    # --- multi-peak detection flags ---
    ap.add_argument("--find-multiple-peaks", action="store_true", default=False,
                    help="Run per-module multi-peak detection after the main "
                         "figure; emit <stem>_multi_peak_modules.csv/.pdf.")
    ap.add_argument("--min-events-per-module", type=int, default=50,
                    metavar="N",
                    help="Skip modules with fewer than N events during peak "
                         "detection (default 50).")
    ap.add_argument("--min-peak-separation", type=float, default=100.0,
                    metavar="MEV",
                    help="Minimum separation between detected peaks in MeV "
                         "(default 100.0).")
    ap.add_argument("--peak-prominence", type=float, default=5.0,
                    metavar="COUNTS",
                    help="Floor on peak prominence threshold; actual threshold "
                         "is max(peak-prominence, 3*sqrt(local_background)) "
                         "(default 5.0).")
    args = ap.parse_args(argv)

    # Warn early if scipy is absent and the feature was requested
    if args.find_multiple_peaks and not HAVE_SCIPY:
        print(
            "[WARNING] --find-multiple-peaks requested but scipy is not "
            "importable.  Peak detection skipped.  Install scipy to enable.",
            file=sys.stderr,
        )

    in_path = Path(args.input_tsv)
    if not in_path.is_file():
        raise SystemExit(f"[ERROR] not a file: {in_path}")

    out_base = Path(args.output_png)
    # Strip .png suffix for stem construction; re-add per output.
    stem = out_base.stem if out_base.suffix.lower() == ".png" else out_base.name
    out_dir = out_base.parent

    configs  = ["legacy", "gated"] if args.config  == "both" else [args.config]
    energies = ["raw", "calibrated"] if args.energy == "both" else [args.energy]

    e_range = tuple(args.e_range)

    peak_ctx = _PeakCtx(
        min_events=args.min_events_per_module,
        min_peak_separation=args.min_peak_separation,
        peak_prominence=args.peak_prominence,
        e_range=e_range,
        e_bins=args.e_bins,
    )

    for cfg in configs:
        for ecol in energies:
            # Build output suffix when iterating over 'both'
            suffix_parts = []
            if args.config == "both":
                suffix_parts.append(cfg)
            if args.energy == "both":
                suffix_parts.append(ecol)
            suffix = ("_" + "_".join(suffix_parts)) if suffix_parts else ""

            print(f"[load] {in_path}  config={cfg}  energy={ecol}", flush=True)
            x_mm, y_mm, sid, se, ce = _load_tsv(in_path, cfg, ecol)
            print(f"[load] {len(ce)} rows loaded", flush=True)

            fig = make_figure(
                x_mm, y_mm, sid, se, ce,
                config=cfg,
                energy_label=ecol,
                e_range=e_range,
                e_bins=args.e_bins,
                xy_bins=args.xy_bins,
                mod_bins=args.module_bins,
                dpi=args.dpi,
            )

            png_path = out_dir / f"{stem}{suffix}.png"
            pdf_path = out_dir / f"{stem}{suffix}.pdf"
            fig.savefig(png_path, dpi=args.dpi)
            print(f"  wrote: {png_path}")
            fig.savefig(pdf_path)
            print(f"  wrote: {pdf_path}")
            plt.close(fig)

            if args.find_multiple_peaks and HAVE_SCIPY:
                _run_multi_peak(x_mm, y_mm, sid, ce,
                                peak_ctx, cfg, ecol,
                                out_dir, stem, suffix)

    return 0


if __name__ == "__main__":
    sys.exit(main())
