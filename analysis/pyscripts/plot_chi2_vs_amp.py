#!/usr/bin/env python3
"""
plot_chi2_vs_amp.py — diagnostic for amplitude bias in the chi2/dof gate
=========================================================================

Background
----------
In WaveAnalyzer::FitPulseShape the per-sample weight used for chi2 is

    sigma_i = ped_rms / peak_amp

i.e. the noise is expressed *relative* to the pulse peak.  This means that
for low-amplitude pulses the absolute noise is small but the relative weight
sigma_i is large, and vice-versa.  The net effect on chi2/dof can create an
**amplitude bias**: pulses near the height-min threshold may systematically
yield lower *or* higher chi2/dof than bright pulses with the same functional
shape, purely because of how sigma scales with amplitude.

This script tests that hypothesis directly by plotting chi2/dof as a
function of peak amplitude using the per-pulse dump produced by
fit_pulse_template.py when --plot-dir is set
(``<plot-dir>/per_pulse_amp_chi2.npz``).

Outputs (written to --out-dir, default = directory containing the .npz)
-------
  chi2_vs_amp_<mtype>.png        — 1×2 per-module-type figure:
      left  : hexbin density (log-log axes, log colour scale)
      right : running median ± 16th/84th percentile in log-spaced amp bins
  chi2_vs_amp_all_types.png      — combined overlay of running-median curves

Console summary
---------------
For each module type the script prints the median chi2/dof in the lowest-
amplitude decile vs. the highest-amplitude decile.  A ratio appreciably
different from 1 confirms an amplitude bias.

Usage
-----
    python plot_chi2_vs_amp.py per_pulse_amp_chi2.npz [--out-dir plots/]
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Match the colour convention used in fit_pulse_template.py
TYPE_COLORS = {
    "PbGlass": "C0",
    "PbWO4":   "C1",
    "LMS":     "C2",
    "Veto":    "C3",
    "Unknown": "0.5",
}

# Horizontal reference line: the default --chi2-max gate in fit_pulse_template
CHI2_REF = 3.0

# Number of log-spaced amplitude bins for the running-percentile panel
N_AMP_BINS = 30

# Font sizes
FS_TITLE  = 11
FS_LABEL  = 10
FS_TICK   = 9
FS_LEGEND = 9


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _running_percentiles(amp: np.ndarray, chi2: np.ndarray,
                          n_bins: int = N_AMP_BINS):
    """
    Return (bin_centres, p16, p50, p84) arrays computed in log-spaced
    amplitude bins.  Bins with fewer than 5 entries are dropped.
    """
    amp_min = amp[amp > 0].min() if (amp > 0).any() else 1.0
    amp_max = amp.max()
    if amp_max <= amp_min:
        amp_max = amp_min * 10.0
    edges = np.logspace(np.log10(amp_min), np.log10(amp_max), n_bins + 1)
    centres, p16_out, p50_out, p84_out = [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (amp >= lo) & (amp < hi)
        if mask.sum() < 5:
            continue
        c2 = chi2[mask]
        centres.append(np.sqrt(lo * hi))   # geometric centre
        p16_out.append(float(np.percentile(c2, 16)))
        p50_out.append(float(np.percentile(c2, 50)))
        p84_out.append(float(np.percentile(c2, 84)))
    return (np.asarray(centres), np.asarray(p16_out),
            np.asarray(p50_out), np.asarray(p84_out))


def _decile_summary(amp: np.ndarray, chi2: np.ndarray, mtype: str):
    """
    Print the median chi2 in the lowest- and highest-amplitude deciles and
    return (low_amp_thresh, high_amp_thresh, low_med, high_med).
    """
    lo_thresh = float(np.percentile(amp, 10))
    hi_thresh = float(np.percentile(amp, 90))
    low_mask  = amp <= lo_thresh
    high_mask = amp >= hi_thresh
    low_med   = float(np.median(chi2[low_mask]))  if low_mask.any()  else float("nan")
    high_med  = float(np.median(chi2[high_mask])) if high_mask.any() else float("nan")
    ratio = high_med / low_med if (low_med > 0 and not np.isnan(low_med)) else float("nan")
    print(
        f"[{mtype}] low-amp decile (< {lo_thresh:.1f} ADC): median chi2 = {low_med:.3f};  "
        f"high-amp decile (> {hi_thresh:.1f} ADC): median chi2 = {high_med:.3f}  "
        f"(ratio high/low = {ratio:.3f})",
        flush=True,
    )
    return lo_thresh, hi_thresh, low_med, high_med


# ---------------------------------------------------------------------------
# Per-type figure  (1×2 layout: hexbin | running percentiles)
# ---------------------------------------------------------------------------

def _plot_single_type(amp: np.ndarray, chi2: np.ndarray,
                      mtype: str, n_pulses: int,
                      out_path: Path) -> None:
    color = TYPE_COLORS.get(mtype, "C4")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(
        f"χ²/dof vs peak amplitude — {mtype}  ({n_pulses:,} accepted pulses)",
        fontsize=FS_TITLE,
    )

    # --- left: hexbin density -----------------------------------------------
    ax = axes[0]
    amp_pos  = amp[amp > 0]
    chi2_pos = chi2[amp > 0]
    if len(amp_pos) > 0:
        hb = ax.hexbin(
            amp_pos, chi2_pos,
            xscale="log", yscale="log",
            bins="log", mincnt=1,
            gridsize=60,
            cmap="viridis",
        )
        cb = fig.colorbar(hb, ax=ax)
        cb.set_label("log₁₀(count)", fontsize=FS_LABEL)
    ax.axhline(CHI2_REF, color="red", lw=1.0, ls="--",
               label=f"chi2 gate = {CHI2_REF}")
    ax.set_xlabel("Peak amplitude (ADC)", fontsize=FS_LABEL)
    ax.set_ylabel("χ²/dof", fontsize=FS_LABEL)
    ax.set_title("Density (hexbin, log–log)", fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.6)
    ax.legend(fontsize=FS_LEGEND)

    # --- right: running percentile bands ------------------------------------
    ax = axes[1]
    centres, p16, p50, p84 = _running_percentiles(amp_pos, chi2_pos)
    if len(centres) > 0:
        ax.fill_between(centres, p16, p84, alpha=0.25, color=color,
                        label="16th–84th pct")
        ax.plot(centres, p50, color=color, lw=1.8, label="median")
    ax.axhline(CHI2_REF, color="red", lw=1.0, ls="--",
               label=f"chi2 gate = {CHI2_REF}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Peak amplitude (ADC)", fontsize=FS_LABEL)
    ax.set_ylabel("χ²/dof", fontsize=FS_LABEL)
    ax.set_title("Running median ± 16/84th percentile", fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.6)
    ax.legend(fontsize=FS_LEGEND)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"[plot] wrote {out_path}", flush=True)


# ---------------------------------------------------------------------------
# Combined overlay figure  (running medians for all types, one axes)
# ---------------------------------------------------------------------------

def _plot_all_types(data_by_type: dict, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_title(
        "χ²/dof vs peak amplitude — all module types (running median)",
        fontsize=FS_TITLE,
    )
    for mtype, (amp, chi2) in sorted(data_by_type.items()):
        color = TYPE_COLORS.get(mtype, "C4")
        amp_pos  = amp[amp > 0]
        chi2_pos = chi2[amp > 0]
        if len(amp_pos) == 0:
            continue
        centres, _, p50, _ = _running_percentiles(amp_pos, chi2_pos)
        if len(centres) == 0:
            continue
        n = len(amp)
        ax.plot(centres, p50, color=color, lw=1.8,
                label=f"{mtype} (n={n:,})")

    ax.axhline(CHI2_REF, color="red", lw=1.0, ls="--",
               label=f"chi2 gate = {CHI2_REF}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Peak amplitude (ADC)", fontsize=FS_LABEL)
    ax.set_ylabel("χ²/dof (median per bin)", fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.6)
    ax.legend(fontsize=FS_LEGEND, loc="best")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"[plot] wrote {out_path}", flush=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot chi2/dof vs peak amplitude from per_pulse_amp_chi2.npz."
    )
    parser.add_argument(
        "npz",
        help="Path to per_pulse_amp_chi2.npz produced by fit_pulse_template.py.",
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Output directory for PNG files (default: same directory as the npz).",
    )
    args = parser.parse_args()

    npz_path = Path(args.npz).resolve()
    if not npz_path.exists():
        sys.exit(f"[error] file not found: {npz_path}")

    out_dir = Path(args.out_dir).resolve() if args.out_dir else npz_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[load] {npz_path}", flush=True)
    data = np.load(npz_path, allow_pickle=False)
    amp   = data["amp"].astype(np.float64)
    chi2  = data["chi2"].astype(np.float64)
    mtype = data["mtype"]          # string array (numpy U* or bytes dtype)

    # Normalise string dtype (numpy may load as bytes on some versions)
    if mtype.dtype.kind == "S":
        mtype = mtype.astype(str)

    print(f"[load] {len(amp):,} pulses total", flush=True)

    # Split by module type
    unique_types = sorted(np.unique(mtype))
    data_by_type: dict = {}
    for mt in unique_types:
        mask = mtype == mt
        data_by_type[mt] = (amp[mask], chi2[mask])
        print(f"  {mt}: {mask.sum():,} pulses", flush=True)

    # Per-type figures
    for mt, (a, c) in data_by_type.items():
        safe_name = mt.replace(" ", "_")
        out_path  = out_dir / f"chi2_vs_amp_{safe_name}.png"
        _plot_single_type(a, c, mt, len(a), out_path)

    # Combined overlay figure
    _plot_all_types(data_by_type, out_dir / "chi2_vs_amp_all_types.png")

    # Console summary: amplitude-bias test
    print("\n[summary] Amplitude-bias test (median chi2 in lowest vs highest decile):",
          flush=True)
    for mt, (a, c) in sorted(data_by_type.items()):
        if len(a) < 20:
            print(f"[{mt}] too few pulses ({len(a)}) — skipped", flush=True)
            continue
        _decile_summary(a, c, mt)


if __name__ == "__main__":
    main()
