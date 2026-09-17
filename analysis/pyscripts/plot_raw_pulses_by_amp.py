"""
plot_raw_pulses_by_amp.py — Test A: visualize raw waveforms per channel split
by amplitude range.

PURPOSE
-------
Tests whether low-amplitude and high-amplitude pulses on the same HyCal channel
have distinct waveform shape families, supporting or refuting the hypothesis
that low-amp pulses (< 300 ADC) are backgrounds while high-amp pulses (> 500
ADC) are true physics signals.

INPUT
-----
The ``waveforms_by_amp/`` directory produced by fit_pulse_template.py when run
with --plot (or equivalent plotting-enabled flag).  Each ``<channel>.npz``
contains:

    module_type        : scalar string (e.g. "PbWO4")
    pulses_<label>     : float32 array, shape (N, n_samples) — pedsub ADC
    amps_<label>       : float32 array, shape (N,)           — fitted peak amp

where ``<label>`` is one of ``lt_300``, ``300_to_500``, ``gt_500``.

OUTPUT
------
One PNG per channel: ``<out_dir>/pulses_by_amp_<channel>.png``

LAYOUT
------
2-row × k-column grid (k = number of populated amplitude bins):
  Row 0: raw pedsub ADC waveforms — faint overlaid lines + bold median.
  Row 1: normalized waveforms (each divided by its own peak value) — same style.

The x-axis is in nanoseconds using ``--clk-ns`` (default 4.0 ns/sample).

USAGE
-----
    python plot_raw_pulses_by_amp.py  waveforms_by_amp/  \\
        --out-dir  plots/waveforms_by_amp/_plots          \\
        --channels W042,W043                              \\
        --clk-ns   4.0
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# Bin metadata: display label and axis-title fragment
# ---------------------------------------------------------------------------

_BIN_DISPLAY = {
    "lt_300":     "peak amp < 300 ADC",
    "300_to_500": "300 ≤ peak amp < 500 ADC",
    "gt_500":     "peak amp ≥ 500 ADC",
}

# Preferred order when multiple bins are present
_BIN_ORDER = ["lt_300", "300_to_500", "gt_500"]


# ---------------------------------------------------------------------------
# Core plotting routine
# ---------------------------------------------------------------------------

def plot_channel(channel: str, npz_path: Path, out_dir: Path, clk_ns: float) -> None:
    """Load one channel npz and produce a 2-row overlay figure."""

    data = np.load(npz_path, allow_pickle=False)

    module_type = str(data["module_type"])

    # Discover which bins are actually present in this file
    present_bins = []
    for lbl in _BIN_ORDER:
        if f"pulses_{lbl}" in data:
            present_bins.append(lbl)

    if not present_bins:
        print(f"  [{channel}] no pulse arrays found — skipping", flush=True)
        return

    n_bins = len(present_bins)
    fig, axes = plt.subplots(
        2, n_bins,
        figsize=(4.5 * n_bins, 7),
        squeeze=False,
        constrained_layout=True,
    )

    # Determine a common y-range across all bins (raw row) so panels are
    # directly comparable; same idea for the normalised row.
    all_raw_min, all_raw_max = np.inf, -np.inf
    bin_data: dict = {}

    for lbl in present_bins:
        pulses = data[f"pulses_{lbl}"].astype(np.float64)  # (N, n_samples)
        amps   = data[f"amps_{lbl}"].astype(np.float64)    # (N,)
        bin_data[lbl] = (pulses, amps)
        all_raw_min = min(all_raw_min, float(pulses.min()))
        all_raw_max = max(all_raw_max, float(pulses.max()))

    raw_margin = (all_raw_max - all_raw_min) * 0.05 or 10.0

    # Print summary header
    print(f"\n  [{channel} {module_type}]", flush=True)

    for col, lbl in enumerate(present_bins):
        pulses, amps = bin_data[lbl]
        n_pulses, n_samples = pulses.shape
        t_ns = clk_ns * np.arange(n_samples)

        median_amp = float(np.median(amps))
        display_label = _BIN_DISPLAY.get(lbl, lbl)

        # --- Row 0: raw pedsub ADC ---
        ax_raw = axes[0, col]
        for i in range(n_pulses):
            ax_raw.plot(t_ns, pulses[i], color="C0", alpha=0.3, lw=0.7)
        median_wf = np.median(pulses, axis=0)
        ax_raw.plot(t_ns, median_wf, color="C1", lw=2.0, label="median")
        ax_raw.set_title(f"{display_label}\n(N={n_pulses})", fontsize=9)
        ax_raw.set_ylabel("pedsub ADC" if col == 0 else "")
        ax_raw.set_xlabel("time (ns)")
        ax_raw.set_ylim(
            all_raw_min - raw_margin,
            all_raw_max + raw_margin,
        )
        ax_raw.legend(fontsize=7, loc="upper right")

        # --- Row 1: normalized (peak = 1 per pulse) ---
        ax_norm = axes[1, col]
        norm_all_min, norm_all_max = np.inf, -np.inf
        for i in range(n_pulses):
            peak_val = float(np.max(np.abs(pulses[i])))
            if peak_val == 0.0:
                continue
            norm_wf = pulses[i] / peak_val
            ax_norm.plot(t_ns, norm_wf, color="C2", alpha=0.3, lw=0.7)
            norm_all_min = min(norm_all_min, float(norm_wf.min()))
            norm_all_max = max(norm_all_max, float(norm_wf.max()))

        # Median of the normalized pulses
        peak_vals = np.max(np.abs(pulses), axis=1)
        safe_mask = peak_vals > 0
        if safe_mask.any():
            normed = pulses[safe_mask] / peak_vals[safe_mask, np.newaxis]
            median_norm = np.median(normed, axis=0)
            ax_norm.plot(t_ns, median_norm, color="C3", lw=2.0, label="median")

        ax_norm.set_ylabel("normalized (peak = 1)" if col == 0 else "")
        ax_norm.set_xlabel("time (ns)")
        if np.isfinite(norm_all_min) and np.isfinite(norm_all_max):
            norm_margin = (norm_all_max - norm_all_min) * 0.05 or 0.05
            ax_norm.set_ylim(norm_all_min - norm_margin, norm_all_max + norm_margin)
        ax_norm.legend(fontsize=7, loc="upper right")

        # Console summary line
        print(
            f"    bin={lbl:<12s} N={n_pulses:<4d}  "
            f"median amp = {median_amp:.0f} ADC",
            flush=True,
        )

    # Overall figure title
    total_n = sum(bin_data[lbl][0].shape[0] for lbl in present_bins)
    fig.suptitle(
        f"{channel}  [{module_type}]  —  {total_n} total cached pulses\n"
        "Row 0: raw pedsub ADC    Row 1: normalized (peak = 1 per pulse)",
        fontsize=10,
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"pulses_by_amp_{channel}.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    print(f"    wrote {out_path.name}", flush=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=(
            "Test A diagnostic: overlay raw waveforms per channel split by "
            "amplitude bin.  Reads the waveforms_by_amp/ directory produced "
            "by fit_pulse_template.py --plot."
        )
    )
    p.add_argument(
        "waveform_dir",
        type=Path,
        help="Path to the waveforms_by_amp/ directory produced by fit_pulse_template.py.",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory for PNG files.  Defaults to <waveform_dir>/_plots.",
    )
    p.add_argument(
        "--channels",
        type=str,
        default="",
        help=(
            "Comma-separated list of channel names to plot (e.g. W042,W043). "
            "If omitted, all channels found in the directory are plotted."
        ),
    )
    p.add_argument(
        "--clk-ns",
        type=float,
        default=4.0,
        help="ADC sample clock period in nanoseconds (default: 4.0).",
    )
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    wf_dir: Path = args.waveform_dir
    if not wf_dir.is_dir():
        print(f"ERROR: waveform_dir not found: {wf_dir}", file=sys.stderr)
        sys.exit(1)

    out_dir: Path = args.out_dir if args.out_dir is not None else wf_dir / "_plots"

    # Resolve which channels to process
    if args.channels.strip():
        requested = [c.strip() for c in args.channels.split(",") if c.strip()]
    else:
        # All .npz files in the directory
        requested = sorted(p.stem for p in wf_dir.glob("*.npz"))

    if not requested:
        print("No channels found to plot.", file=sys.stderr)
        sys.exit(1)

    print(
        f"[plot_raw_pulses_by_amp] waveform_dir={wf_dir}  "
        f"out_dir={out_dir}  channels={len(requested)}  "
        f"clk_ns={args.clk_ns}",
        flush=True,
    )

    n_ok = 0
    for channel in requested:
        npz_path = wf_dir / f"{channel}.npz"
        if not npz_path.exists():
            print(f"  [{channel}] npz not found — skipping", flush=True)
            continue
        plot_channel(channel, npz_path, out_dir, args.clk_ns)
        n_ok += 1

    print(f"\n[done] {n_ok}/{len(requested)} channels plotted → {out_dir}", flush=True)


if __name__ == "__main__":
    main()
