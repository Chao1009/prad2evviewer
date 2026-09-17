#!/usr/bin/env python3
"""
plot_tau_vs_amp.py — diagnostic for amplitude-dependent pulse shape parameters
===============================================================================

Background
----------
In the two-tau (and two-tau-p) pulse model the waveform shape is governed by
three time constants:

    t0    — arrival time of the pulse leading edge (ns)
    tau_r — rise time constant (ns)
    tau_f — fall / decay time constant (ns)

If the detector readout chain is linear, these parameters should be
*amplitude-independent*: a dim pulse and a bright pulse from the same module
should have the same shape, just scaled.  Departures from this expectation can
signal two distinct effects:

  1. **PMT saturation** — at large pulse heights the dynode chain is driven into
     a nonlinear regime; tau_r or tau_f shift with amplitude *within* a single
     channel.

  2. **Per-module intrinsic variation** — different physical modules (or
     different HV settings) yield different template parameters, but within any
     one module the shape is amplitude-stable.

This script tests these hypotheses by reading the per-pulse
``per_pulse_amp_chi2.npz`` produced by fit_pulse_template.py (--plot-dir) and
generating:

  tau_r_vs_amp_<material>.png
      Hexbin density (log amplitude axis, linear tau axis) + running
      median ± 16/84th percentile bands in log-spaced amplitude bins,
      for each requested material type.

  tau_f_vs_amp_<material>.png
      Same layout but for tau_f.

  tau_vs_amp_channel_<name>.png
      Per-channel scatter plots of tau_r and tau_f vs amplitude for each
      explicitly requested channel (--channels).  Per-channel counts are
      typically hundreds to a few thousand, so a plain scatter + running
      median is used instead of hexbin.

A console summary is printed at the end:
  - For each material: median tau_r/tau_f in the lowest-amplitude decile
    vs the highest-amplitude decile, with ratio.
  - For each explicitly plotted channel: same low-/high-amp decile summary.

A ratio significantly different from 1 confirms amplitude-dependent pulse
shape (hypothesis 1).  If per-material ratios are near 1 but individual
channels cluster at different tau values, that points to hypothesis 2.

Usage
-----
    python plot_tau_vs_amp.py per_pulse_amp_chi2.npz \\
        [--out-dir plots/] \\
        [--materials PbWO4,PbGlass] \\
        [--channels W042,W100,G100] \\
        [--min-pulses-per-channel 100]

Dependencies: numpy, matplotlib (headless Agg backend), pathlib, argparse.
No prad2py or other project-specific imports are used.
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Constants — match fit_pulse_template.py conventions
# ---------------------------------------------------------------------------

TYPE_COLORS = {
    "PbGlass": "C0",
    "PbWO4":   "C1",
    "LMS":     "C2",
    "Veto":    "C3",
    "Unknown": "0.5",
}

# Number of log-spaced amplitude bins for the running-percentile panel
N_AMP_BINS = 30

# Minimum per-bin count to include in the running-percentile curve
MIN_BIN_COUNT = 5

# Font sizes
FS_TITLE  = 11
FS_LABEL  = 10
FS_TICK   = 9
FS_LEGEND = 9


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _running_percentiles(amp: np.ndarray, tau: np.ndarray,
                          n_bins: int = N_AMP_BINS):
    """Return (bin_centres, p16, p50, p84) in log-spaced amplitude bins.

    Bins with fewer than MIN_BIN_COUNT entries are omitted.  Returns empty
    arrays if there are no valid bins.
    """
    pos = amp > 0
    if not pos.any():
        return (np.array([]), np.array([]), np.array([]), np.array([]))
    amp_min = float(amp[pos].min())
    amp_max = float(amp.max())
    if amp_max <= amp_min:
        amp_max = amp_min * 10.0
    edges = np.logspace(np.log10(amp_min), np.log10(amp_max), n_bins + 1)
    centres, p16_out, p50_out, p84_out = [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (amp >= lo) & (amp < hi)
        if mask.sum() < MIN_BIN_COUNT:
            continue
        t = tau[mask]
        centres.append(float(np.sqrt(lo * hi)))   # geometric centre
        p16_out.append(float(np.percentile(t, 16)))
        p50_out.append(float(np.percentile(t, 50)))
        p84_out.append(float(np.percentile(t, 84)))
    return (np.asarray(centres), np.asarray(p16_out),
            np.asarray(p50_out), np.asarray(p84_out))


def _decile_summary_tau(amp: np.ndarray, tau: np.ndarray,
                         label: str, param_name: str) -> dict:
    """Compute and print low-/high-amplitude decile median for *tau*.

    Returns a dict with keys: lo_thresh, hi_thresh, low_med, high_med, ratio.
    """
    lo_thresh = float(np.percentile(amp, 10))
    hi_thresh = float(np.percentile(amp, 90))
    low_mask  = amp <= lo_thresh
    high_mask = amp >= hi_thresh
    low_med  = float(np.median(tau[low_mask]))  if low_mask.any()  else float("nan")
    high_med = float(np.median(tau[high_mask])) if high_mask.any() else float("nan")
    if low_med > 0 and not np.isnan(low_med):
        ratio = high_med / low_med
    else:
        ratio = float("nan")
    print(
        f"  [{label}] low-amp decile (< {lo_thresh:.1f} ADC): "
        f"median {param_name} = {low_med:.2f} ns",
        flush=True,
    )
    print(
        f"  [{label}] high-amp decile (> {hi_thresh:.1f} ADC): "
        f"median {param_name} = {high_med:.2f} ns",
        flush=True,
    )
    print(
        f"  [{label}] ratio (high/low) = {ratio:.3f}"
        "  ← if != 1, amplitude-dependent shape",
        flush=True,
    )
    return dict(lo_thresh=lo_thresh, hi_thresh=hi_thresh,
                low_med=low_med, high_med=high_med, ratio=ratio)


# ---------------------------------------------------------------------------
# Per-material figures  (hexbin | running percentile, 1×2 layout)
# ---------------------------------------------------------------------------

def _plot_material_tau(amp: np.ndarray, tau: np.ndarray,
                        material: str, param: str,
                        n_pulses: int, out_path: Path) -> None:
    """One 1×2 figure: hexbin density (left) + running percentile (right).

    Parameters
    ----------
    amp      : peak amplitude array (ADC), same length as tau.
    tau      : tau_r or tau_f array (ns).
    material : material label string (e.g. "PbWO4").
    param    : "tau_r" or "tau_f" — used in axis labels and title.
    n_pulses : total number of pulses (for the title).
    out_path : destination .png path.
    """
    color = TYPE_COLORS.get(material, "C4")
    param_label = r"$\tau_r$" if param == "tau_r" else r"$\tau_f$"

    # Interpretation hint shown in the title
    interp_hint = ("flat = amplitude-independent shape; "
                   "slope = amplitude-dependent shape (e.g. PMT saturation)")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        f"{param_label} vs peak amplitude — {material}  "
        f"({n_pulses:,} pulses)\n{interp_hint}",
        fontsize=FS_TITLE,
    )

    amp_pos = amp[amp > 0]
    tau_pos = tau[amp > 0]

    # --- left: hexbin density -----------------------------------------------
    ax = axes[0]
    if len(amp_pos) > 0:
        hb = ax.hexbin(
            amp_pos, tau_pos,
            xscale="log", yscale="linear",
            bins="log", mincnt=1,
            gridsize=60,
            cmap="viridis",
        )
        cb = fig.colorbar(hb, ax=ax)
        cb.set_label("log₁₀(count)", fontsize=FS_LABEL)
    ax.set_xlabel("Peak amplitude (ADC)", fontsize=FS_LABEL)
    ax.set_ylabel(f"{param_label} (ns)", fontsize=FS_LABEL)
    ax.set_title("Density (hexbin, log amp axis)", fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.6)

    # --- right: running percentile bands ------------------------------------
    ax = axes[1]
    centres, p16, p50, p84 = _running_percentiles(amp_pos, tau_pos)
    if len(centres) > 0:
        ax.fill_between(centres, p16, p84, alpha=0.25, color=color,
                        label="16th–84th pct")
        ax.plot(centres, p50, color=color, lw=1.8, label="median")
    ax.set_xscale("log")
    ax.set_xlabel("Peak amplitude (ADC)", fontsize=FS_LABEL)
    ax.set_ylabel(f"{param_label} (ns)", fontsize=FS_LABEL)
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
# Per-channel figure  (scatter + running median, 1×2 layout)
# ---------------------------------------------------------------------------

def _plot_channel_tau(amp: np.ndarray, tau_r: np.ndarray, tau_f: np.ndarray,
                       channel: str, mtype: str, n_pulses: int,
                       out_path: Path) -> None:
    """Two-panel scatter for a single channel: tau_r (left), tau_f (right).

    Plain scatter is used instead of hexbin because per-channel counts are
    typically hundreds to a few thousand.  A contrasting running-median line
    is overlaid.
    """
    color_r = TYPE_COLORS.get(mtype, "C4")
    color_f = "C5"   # contrasting colour for tau_f panel

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        f"τ vs amplitude — channel {channel}  [{mtype}]  ({n_pulses:,} pulses)",
        fontsize=FS_TITLE,
    )

    for ax, tau, param, color in [
        (axes[0], tau_r, r"$\tau_r$", color_r),
        (axes[1], tau_f, r"$\tau_f$", color_f),
    ]:
        amp_pos = amp[amp > 0]
        tau_pos = tau[amp > 0]

        # Scatter of individual pulses
        ax.scatter(
            amp_pos, tau_pos,
            s=4, alpha=0.35, color=color, rasterized=True, label="pulse",
        )

        # Running median overlay
        centres, _, p50, _ = _running_percentiles(amp_pos, tau_pos)
        if len(centres) > 0:
            ax.plot(centres, p50,
                    color="black", lw=1.6, zorder=5, label="running median")

        ax.set_xscale("log")
        ax.set_xlabel("Peak amplitude (ADC)", fontsize=FS_LABEL)
        ax.set_ylabel(f"{param} (ns)", fontsize=FS_LABEL)
        ax.set_title(f"{param} vs amplitude", fontsize=FS_LABEL)
        ax.tick_params(labelsize=FS_TICK)
        ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.6)
        ax.legend(fontsize=FS_LEGEND)

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
        description=(
            "Plot τ_r and τ_f vs peak amplitude from per_pulse_amp_chi2.npz "
            "to diagnose amplitude-dependent pulse shape (e.g. PMT saturation)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "npz_path",
        help="Path to per_pulse_amp_chi2.npz produced by fit_pulse_template.py.",
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Output directory for PNG files (default: same directory as the npz).",
    )
    parser.add_argument(
        "--materials",
        default="PbWO4,PbGlass",
        help="Comma-separated list of material labels to produce per-material "
             "plots for (default: 'PbWO4,PbGlass').",
    )
    parser.add_argument(
        "--channels",
        default="",
        help="Comma-separated channel names for per-channel plots "
             "(e.g. 'W042,W100,G100'; default: none).",
    )
    parser.add_argument(
        "--min-pulses-per-channel",
        type=int,
        default=100,
        metavar="N",
        help="Minimum number of pulses a channel must have to be plotted "
             "individually (default: 100).",
    )
    args = parser.parse_args()

    # --- resolve paths -------------------------------------------------------
    npz_path = Path(args.npz_path).resolve()
    if not npz_path.exists():
        sys.exit(f"[error] file not found: {npz_path}")

    out_dir = Path(args.out_dir).resolve() if args.out_dir else npz_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    requested_materials = [m.strip() for m in args.materials.split(",") if m.strip()]
    requested_channels  = [c.strip() for c in args.channels.split(",")
                           if c.strip()] if args.channels else []

    # --- load ----------------------------------------------------------------
    print(f"[load] {npz_path}", flush=True)
    data = np.load(npz_path, allow_pickle=False)

    amp   = data["amp"].astype(np.float64)
    tau_r = data["tau_r"].astype(np.float64)
    tau_f = data["tau_f"].astype(np.float64)
    mtype = data["mtype"]
    name  = data["name"]

    # Normalise string dtypes (numpy may load as bytes S* on some versions)
    if mtype.dtype.kind == "S":
        mtype = mtype.astype(str)
    if name.dtype.kind == "S":
        name = name.astype(str)

    print(f"[load] {len(amp):,} pulses total", flush=True)

    # --- per-material plots --------------------------------------------------
    print("\n=== per-material τ_r vs amplitude summary ===", flush=True)
    for material in requested_materials:
        mask = mtype == material
        n = int(mask.sum())
        if n == 0:
            print(f"[{material}] no pulses found in npz — skipping", flush=True)
            continue
        a = amp[mask]
        tr = tau_r[mask]
        tf = tau_f[mask]
        print(f"[{material}] N pulses = {n:,}", flush=True)

        # tau_r figure
        safe = material.replace(" ", "_")
        _plot_material_tau(a, tr, material, "tau_r", n,
                           out_dir / f"tau_r_vs_amp_{safe}.png")

        # tau_f figure
        _plot_material_tau(a, tf, material, "tau_f", n,
                           out_dir / f"tau_f_vs_amp_{safe}.png")

        # Console summary (tau_r then tau_f)
        if n >= 20:
            _decile_summary_tau(a, tr, material, "τ_r")
            _decile_summary_tau(a, tf, material, "τ_f")
        else:
            print(f"  [{material}] too few pulses ({n}) — decile summary skipped",
                  flush=True)

    # --- per-channel plots ---------------------------------------------------
    if not requested_channels:
        return

    print("\n=== per-channel summary ===", flush=True)
    for ch in requested_channels:
        mask = name == ch
        n = int(mask.sum())
        if n < args.min_pulses_per_channel:
            print(
                f"[warning] channel {ch}: only {n} pulses "
                f"(< --min-pulses-per-channel {args.min_pulses_per_channel}) "
                f"— skipping",
                flush=True,
            )
            continue

        a  = amp[mask]
        tr = tau_r[mask]
        tf = tau_f[mask]

        # Determine the module type for this channel (majority vote)
        mt_vals = mtype[mask]
        if len(mt_vals) > 0:
            unique_mt, counts = np.unique(mt_vals, return_counts=True)
            ch_mtype = str(unique_mt[counts.argmax()])
        else:
            ch_mtype = "Unknown"

        safe_ch = ch.replace(" ", "_")
        _plot_channel_tau(a, tr, tf, ch, ch_mtype, n,
                          out_dir / f"tau_vs_amp_channel_{safe_ch}.png")

        # Console decile summary
        print(f"[{ch} {ch_mtype} N={n:,}]", flush=True)
        if n >= 20:
            _decile_summary_tau(a, tr, ch, "τ_r")
            _decile_summary_tau(a, tf, ch, "τ_f")
        else:
            print(f"  too few pulses ({n}) — decile summary skipped", flush=True)


if __name__ == "__main__":
    main()
