#!/usr/bin/env python3
"""
plot_template_by_crate.py — diagnostic: do PbWO4 pulse-shape sub-populations
split by FADC crate (roc_tag)?
=============================================================================

Background
----------
When inspecting the per-channel pulse-template fit parameters produced by
``fit_pulse_template.py``, a clear **bi-modal** distribution was observed in
PbWO4 channels:

  * τ_r shows two peaks (~2–3 ns and ~7–8 ns)
  * τ_f shows two peaks (~20–30 ns and ~35–40 ns)
  * The two populations are strongly correlated with each other

The leading hypothesis is that these sub-populations correspond to distinct
**FADC front-end crates** (equivalently, distinct ``roc_tag`` values derived
from the ``channel_id`` field in the JSON).  If crystal or PMT properties
drove the variation one would expect the two groups to be spatially scrambled;
if electronics drive it one would expect clean crate-level separation.

This script tests that hypothesis by:

  1. Parsing the pulse-template JSON, extracting τ_r, τ_f, and t₀ per channel.
  2. Grouping channels by ``roc_tag`` (the first component of ``channel_id``).
  3. Producing overlaid-histogram and scatter plots coloured by crate.
  4. Printing a concise per-crate summary to the console.

The script is **read-only** with respect to physics data — it only loads the
JSON written by ``fit_pulse_template.py`` and writes PNG diagnostic plots.
No prad2py or decoder imports are needed.

Usage
-----
    python plot_template_by_crate.py pulse_templates_<run>.json \\
        [--out-dir plots/] \\
        [--min-pulses 50] \\
        [--materials PbWO4,PbGlass]

Outputs (per material, in --out-dir)
-------------------------------------
  template_by_crate_<mat>_tau_r.png
  template_by_crate_<mat>_tau_f.png
  template_by_crate_<mat>_t0.png
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# Colour helpers
# ---------------------------------------------------------------------------

def _build_crate_colors(roc_tags: list[str]) -> dict[str, tuple]:
    """Return a dict mapping each roc_tag (sorted) to an RGBA colour tuple.

    Uses tab10 for ≤10 crates, tab20 for 11–20, and viridis for >20.
    """
    tags_sorted = sorted(roc_tags)
    n = len(tags_sorted)
    if n <= 10:
        cmap = plt.cm.tab10
        colors = [cmap(i / 10) for i in range(n)]
    elif n <= 20:
        cmap = plt.cm.tab20
        colors = [cmap(i / 20) for i in range(n)]
    else:
        cmap = plt.cm.viridis
        colors = [cmap(i / (n - 1)) for i in range(n)]
    return {tag: col for tag, col in zip(tags_sorted, colors)}


# ---------------------------------------------------------------------------
# JSON parsing
# ---------------------------------------------------------------------------

def _parse_channel_id(channel_id: str) -> tuple[str, str, str]:
    """Split ``roc_tag_slot_channel`` into (roc_tag, slot, channel).

    The format produced by fit_pulse_template.py is ``<roc>_<slot>_<channel>``.
    Returns ('unknown', '', '') if parsing fails.
    """
    parts = channel_id.split("_")
    if len(parts) >= 3:
        return parts[0], parts[1], parts[2]
    elif len(parts) == 2:
        return parts[0], parts[1], ""
    elif len(parts) == 1:
        return parts[0], "", ""
    return "unknown", "", ""


def load_channels(json_path: Path, min_pulses: int) -> dict[str, list[dict]]:
    """Load the pulse-template JSON and return channels grouped by material.

    Returns
    -------
    dict mapping material name → list of dicts with keys:
        name, roc_tag, slot, channel,
        tau_r, tau_f, t0   (all in ns, median values)
    """
    with json_path.open() as fh:
        data = json.load(fh)

    by_material: dict[str, list[dict]] = defaultdict(list)

    for key, entry in data.items():
        # Skip metadata blocks (keys starting with '_')
        if key.startswith("_"):
            continue
        # Skip entries that are not dicts (safety)
        if not isinstance(entry, dict):
            continue

        n_used = entry.get("n_pulses_used", 0)
        if n_used < min_pulses:
            continue

        material = entry.get("module_type", "Unknown")
        channel_id = entry.get("channel_id", "")

        roc_tag, slot, ch = _parse_channel_id(channel_id)

        # Extract median values; skip channel if any required field is absent
        try:
            tau_r = entry["tau_r_ns"]["median"]
            tau_f = entry["tau_f_ns"]["median"]
            t0    = entry["t0_ns"]["median"]
        except (KeyError, TypeError):
            continue

        by_material[material].append({
            "name":    key,
            "roc_tag": roc_tag,
            "slot":    slot,
            "channel": ch,
            "tau_r":   float(tau_r),
            "tau_f":   float(tau_f),
            "t0":      float(t0),
        })

    return dict(by_material)


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _make_figure(
    material: str,
    param_label: str,
    param_values_by_crate: dict[str, np.ndarray],
    tau_r_by_crate: dict[str, np.ndarray],
    tau_f_by_crate: dict[str, np.ndarray],
    crate_colors: dict[str, tuple],
    n_total: int,
    min_pulses: int,
    x_label: str,
) -> plt.Figure:
    """Build a 1×2 figure: overlaid histogram (left) + τ_r vs τ_f scatter (right).

    Parameters
    ----------
    param_label : str
        Short label for the quantity shown in the left panel (e.g. 'τ_r (ns)').
    param_values_by_crate : dict
        roc_tag → array of the parameter being histogrammed.
    tau_r_by_crate / tau_f_by_crate : dict
        roc_tag → arrays used for the scatter plot.
    x_label : str
        Full axis label for the left panel.
    """
    fig, (ax_hist, ax_scat) = plt.subplots(1, 2, figsize=(14, 5))

    # ---- Histogram panel ---------------------------------------------------
    all_vals = np.concatenate(list(param_values_by_crate.values()))
    v_min, v_max = np.nanmin(all_vals), np.nanmax(all_vals)
    bins = np.linspace(v_min, v_max, 31)

    for roc_tag in sorted(param_values_by_crate):
        vals  = param_values_by_crate[roc_tag]
        color = crate_colors[roc_tag]
        ax_hist.hist(
            vals,
            bins=bins,
            alpha=0.55,
            color=color,
            label=f"ROC={roc_tag}  (N={len(vals)})",
            histtype="stepfilled",
            edgecolor=color,
            linewidth=0.8,
        )

    ax_hist.set_xlabel(x_label, fontsize=11)
    ax_hist.set_ylabel("Channels / bin", fontsize=11)
    ax_hist.legend(fontsize=8, framealpha=0.7)
    ax_hist.grid(True, alpha=0.3)

    # ---- Scatter panel: τ_r vs τ_f coloured by crate ----------------------
    for roc_tag in sorted(tau_r_by_crate):
        color = crate_colors[roc_tag]
        ax_scat.scatter(
            tau_r_by_crate[roc_tag],
            tau_f_by_crate[roc_tag],
            s=12,
            alpha=0.55,
            color=color,
            label=f"ROC={roc_tag}",
            linewidths=0,
        )

    ax_scat.set_xlabel("τ_r (ns)", fontsize=11)
    ax_scat.set_ylabel("τ_f (ns)", fontsize=11)
    ax_scat.legend(fontsize=8, framealpha=0.7)
    ax_scat.grid(True, alpha=0.3)
    ax_scat.set_title("τ_r vs τ_f coloured by crate", fontsize=10)

    # ---- Overall title -----------------------------------------------------
    n_crates = len(param_values_by_crate)
    fig.suptitle(
        f"{material} — {param_label} by FADC crate   "
        f"[{n_total} channels, {n_crates} crates, min_pulses≥{min_pulses}]",
        fontsize=12,
        y=1.01,
    )

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Per-material processing
# ---------------------------------------------------------------------------

def process_material(
    material: str,
    channels: list[dict],
    out_dir: Path,
    min_pulses: int,
) -> None:
    """Produce all three diagnostic PNGs for one material and print summary."""

    # Group by roc_tag
    by_crate: dict[str, list[dict]] = defaultdict(list)
    for ch in channels:
        by_crate[ch["roc_tag"]].append(ch)

    roc_tags = sorted(by_crate)
    crate_colors = _build_crate_colors(roc_tags)
    n_total = len(channels)
    n_crates = len(roc_tags)

    # Build per-crate arrays
    tau_r_by_crate: dict[str, np.ndarray] = {}
    tau_f_by_crate: dict[str, np.ndarray] = {}
    t0_by_crate:    dict[str, np.ndarray] = {}

    for roc_tag in roc_tags:
        chs = by_crate[roc_tag]
        tau_r_by_crate[roc_tag] = np.array([c["tau_r"] for c in chs])
        tau_f_by_crate[roc_tag] = np.array([c["tau_f"] for c in chs])
        t0_by_crate[roc_tag]    = np.array([c["t0"]    for c in chs])

    # ---- Plot 1: τ_r histogram + scatter ----------------------------------
    fig1 = _make_figure(
        material=material,
        param_label="τ_r",
        param_values_by_crate=tau_r_by_crate,
        tau_r_by_crate=tau_r_by_crate,
        tau_f_by_crate=tau_f_by_crate,
        crate_colors=crate_colors,
        n_total=n_total,
        min_pulses=min_pulses,
        x_label="τ_r (ns)",
    )
    p1 = out_dir / f"template_by_crate_{material}_tau_r.png"
    fig1.savefig(p1, dpi=150, bbox_inches="tight")
    plt.close(fig1)
    print(f"  Wrote {p1}")

    # ---- Plot 2: τ_f histogram + scatter ----------------------------------
    fig2 = _make_figure(
        material=material,
        param_label="τ_f",
        param_values_by_crate=tau_f_by_crate,
        tau_r_by_crate=tau_r_by_crate,
        tau_f_by_crate=tau_f_by_crate,
        crate_colors=crate_colors,
        n_total=n_total,
        min_pulses=min_pulses,
        x_label="τ_f (ns)",
    )
    p2 = out_dir / f"template_by_crate_{material}_tau_f.png"
    fig2.savefig(p2, dpi=150, bbox_inches="tight")
    plt.close(fig2)
    print(f"  Wrote {p2}")

    # ---- Plot 3: t₀ histogram + scatter -----------------------------------
    fig3 = _make_figure(
        material=material,
        param_label="t₀",
        param_values_by_crate=t0_by_crate,
        tau_r_by_crate=tau_r_by_crate,
        tau_f_by_crate=tau_f_by_crate,
        crate_colors=crate_colors,
        n_total=n_total,
        min_pulses=min_pulses,
        x_label="t₀ (ns)",
    )
    p3 = out_dir / f"template_by_crate_{material}_t0.png"
    fig3.savefig(p3, dpi=150, bbox_inches="tight")
    plt.close(fig3)
    print(f"  Wrote {p3}")

    # ---- Console summary --------------------------------------------------
    print(
        f"\n=== {material}: {n_total} channels across {n_crates} crates"
        f" (min_pulses={min_pulses}) ==="
    )
    for roc_tag in roc_tags:
        arr_r = tau_r_by_crate[roc_tag]
        arr_f = tau_f_by_crate[roc_tag]
        med_r = float(np.median(arr_r))
        med_f = float(np.median(arr_f))
        print(
            f"  ROC {roc_tag:>6s}:  {len(arr_r):5d} channels"
            f"    τ_r median={med_r:5.1f} ns"
            f"    τ_f median={med_f:5.1f} ns"
        )
    print()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Plot pulse-template parameters (τ_r, τ_f, t₀) split by FADC "
            "crate/ROC to test whether observed multi-modal PbWO4 distributions "
            "correlate with front-end electronics grouping."
        )
    )
    p.add_argument(
        "json_path",
        type=Path,
        help="Path to the pulse_templates_*.json produced by fit_pulse_template.py",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Directory for output PNG files (default: same directory as json_path).",
    )
    p.add_argument(
        "--min-pulses",
        type=int,
        default=50,
        metavar="N",
        help="Skip channels with n_pulses_used < N (default: 50).",
    )
    p.add_argument(
        "--materials",
        type=str,
        default="PbWO4,PbGlass",
        metavar="MAT1,MAT2,...",
        help=(
            "Comma-separated list of module_type values to plot "
            "(default: 'PbWO4,PbGlass')."
        ),
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    json_path: Path = args.json_path.resolve()
    if not json_path.exists():
        print(f"ERROR: JSON file not found: {json_path}", file=sys.stderr)
        return 1

    out_dir: Path = args.out_dir if args.out_dir is not None else json_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    min_pulses: int = args.min_pulses
    materials: list[str] = [m.strip() for m in args.materials.split(",") if m.strip()]

    print(f"Loading {json_path} …")
    by_material = load_channels(json_path, min_pulses)

    for material in materials:
        if material not in by_material:
            print(
                f"WARNING: material '{material}' not found in JSON "
                f"(or no channels passed min_pulses={min_pulses} filter). Skipping."
            )
            continue
        channels = by_material[material]
        print(f"\nProcessing {material}: {len(channels)} channels …")
        process_material(material, channels, out_dir, min_pulses)

    return 0


if __name__ == "__main__":
    sys.exit(main())
