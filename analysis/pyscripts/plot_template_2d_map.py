#!/usr/bin/env python3
"""
plot_template_2d_map.py — 2D spatial maps of HyCal pulse-template parameters
=============================================================================

Background
----------
``fit_pulse_template.py`` produces per-channel medians of the two-tau pulse
shape parameters (τ_r, τ_f, t₀) for every HyCal module that accumulated
enough clean isolated pulses.  Histogram inspection of the PbWO4 population
reveals a clear **bi-modal** structure:

  * τ_r shows two peaks (~2–3 ns and ~7–8 ns)
  * τ_f shows two peaks (~20–30 ns and ~35–40 ns)
  * The two populations are strongly correlated

Two competing hypotheses explain this bi-modality:

  A. **Radiation-damage or crystal-batch clustering** — modules that have
     accumulated more dose (or came from a different production batch) have
     measurably slower rise/fall times.  If true, the two populations should
     form *contiguous spatial patches* on the detector face, because dose
     accumulation is determined by the local flux profile (highest at the
     beam hole edge, falling with radius) and batch assignments were made
     by crystal position at installation time.

  B. **Per-PMT or per-electronics variation** — the shape is set by the
     PMT gain, base, or the FADC front-end, not by the crystal.  If true,
     the two populations should appear *spatially scrambled* (no visible
     clustering on the face map).

This script tests hypothesis A by projecting τ_r, τ_f, and t₀ onto the
physical (x, y) positions of modules as recorded in the HyCal geometry map
(``hycal_map.json``).  Each module is drawn as a filled rectangle at its
actual size — PbWO4 modules are ~20.75 × 20.75 mm, PbGlass modules are
~38.15 × 38.15 mm — and colour-coded by its parameter value.

If spatial clustering is present, the 2D map will show a spatially organised
colour gradient (e.g. a ring of slow-rise modules near the beam hole).
If only per-PMT/electronics variation drives the bi-modality, the colour
pattern will look random or crate-striped (the latter being the signature
tested by ``plot_template_by_crate.py``).

A fourth map shows ``n_pulses_used`` per channel on a log colour scale,
making it easy to see where the high-statistics physics signal lives (central
modules) versus where peripheral or background-dominated channels sit.

Usage
-----
    python plot_template_2d_map.py pulse_templates_<run>.json \\
        [--hc-map-file /path/to/hycal_map.json] \\
        [--out-dir plots/] \\
        [--min-pulses 50] \\
        [--materials PbWO4,PbGlass] \\
        [--vmin 2.0] [--vmax 8.0] \\
        [--n-vmin 10] [--n-vmax 5000]

Outputs (per material, per parameter, in --out-dir)
----------------------------------------------------
    template_2d_<material>_tau_r.png
    template_2d_<material>_tau_f.png
    template_2d_<material>_t0.png
    template_2d_<material>_n_pulses.png

No dependencies beyond numpy, matplotlib, argparse, json, pathlib, and os.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
import numpy as np


# ---------------------------------------------------------------------------
# HyCal map loading
# ---------------------------------------------------------------------------

def load_hycal_map(map_path: Path) -> dict[str, tuple[float, float, float, float]]:
    """Load hycal_map.json and return a lookup of module name → (x, y, sx, sy).

    Parameters
    ----------
    map_path : Path
        Absolute path to hycal_map.json.

    Returns
    -------
    dict mapping module name (e.g. 'W042', 'G127') to
    (x_center_mm, y_center_mm, width_mm, height_mm).
    """
    with map_path.open(encoding="utf-8") as fh:
        records = json.load(fh)

    geo_map: dict[str, tuple[float, float, float, float]] = {}
    for rec in records:
        name = rec.get("n", "")
        if not name:
            continue
        geo = rec.get("geo", {})
        try:
            x  = float(geo["x"])
            y  = float(geo["y"])
            sx = float(geo["sx"])
            sy = float(geo["sy"])
        except (KeyError, TypeError, ValueError):
            print(f"  [WARN] hycal_map: skipping module '{name}' — missing or "
                  f"invalid geo fields", file=sys.stderr)
            continue
        geo_map[name] = (x, y, sx, sy)

    return geo_map


# ---------------------------------------------------------------------------
# Template JSON loading
# ---------------------------------------------------------------------------

def load_template_channels(
    json_path: Path,
    min_pulses: int,
    materials: list[str],
) -> dict[str, list[dict]]:
    """Load pulse-template JSON and return channels grouped by material.

    Only channels whose ``module_type`` is in *materials* and whose
    ``n_pulses_used`` ≥ *min_pulses* are returned.  Metadata blocks (keys
    starting with ``_``) are silently skipped.

    Returns
    -------
    dict mapping material name → list of dicts with keys:
        name, tau_r, tau_f, t0   (all in ns, median values)
    """
    with json_path.open(encoding="utf-8") as fh:
        data = json.load(fh)

    by_material: dict[str, list[dict]] = {m: [] for m in materials}

    for key, entry in data.items():
        if key.startswith("_"):
            continue
        if not isinstance(entry, dict):
            continue

        material = entry.get("module_type", "Unknown")
        if material not in by_material:
            continue

        n_used = entry.get("n_pulses_used", 0)
        if n_used < min_pulses:
            continue

        try:
            tau_r = float(entry["tau_r_ns"]["median"])
            tau_f = float(entry["tau_f_ns"]["median"])
            t0    = float(entry["t0_ns"]["median"])
        except (KeyError, TypeError, ValueError):
            continue

        by_material[material].append({
            "name":  key,
            "tau_r": tau_r,
            "tau_f": tau_f,
            "t0":    t0,
            "n_pulses_used": int(n_used),
        })

    return by_material


# ---------------------------------------------------------------------------
# Per-material, per-parameter figure
# ---------------------------------------------------------------------------

_PARAM_META = {
    "tau_r": {
        "key":   "tau_r",
        "label": "τ_r (ns) — rise time",
        "cmap":  "viridis",
        "short": "τ_r",
    },
    "tau_f": {
        "key":   "tau_f",
        "label": "τ_f (ns) — fall time",
        "cmap":  "viridis",
        "short": "τ_f",
    },
    "t0": {
        "key":   "t0",
        "label": "t₀ (ns) — pulse onset",
        "cmap":  "viridis",
        "short": "t₀",
    },
}


def _resolve_clim(
    values: np.ndarray,
    vmin_override: float | None,
    vmax_override: float | None,
) -> tuple[float, float]:
    """Return (vmin, vmax) for the colour scale.

    Uses 2nd and 98th percentile by default; honour user overrides for either
    or both ends independently.
    """
    vmin = float(np.percentile(values, 2))  if vmin_override is None else vmin_override
    vmax = float(np.percentile(values, 98)) if vmax_override is None else vmax_override
    if vmin >= vmax:
        # Degenerate range — fall back to full range
        vmin = float(values.min())
        vmax = float(values.max())
    return vmin, vmax


def make_2d_map(
    material: str,
    param_name: str,
    channels: list[dict],
    geo_map: dict[str, tuple[float, float, float, float]],
    out_path: Path,
    vmin_override: float | None,
    vmax_override: float | None,
) -> None:
    """Render one 2D spatial map and save to *out_path*.

    Parameters
    ----------
    material : str
        Module material label (e.g. 'PbWO4').
    param_name : str
        One of 'tau_r', 'tau_f', 't0'.
    channels : list of dict
        Channel records for *material* (already filtered for min_pulses).
    geo_map : dict
        name → (x, y, sx, sy) from load_hycal_map().
    out_path : Path
        Destination PNG path.
    vmin_override, vmax_override : float or None
        CLI overrides for the colour scale.
    """
    meta = _PARAM_META[param_name]
    param_key = meta["key"]

    # --- Match channels to geometry ----------------------------------------
    xs:     list[float] = []
    ys:     list[float] = []
    sxs:    list[float] = []
    sys_:   list[float] = []
    vals:   list[float] = []
    missing = 0

    for ch in channels:
        name = ch["name"]
        geo = geo_map.get(name)
        if geo is None:
            print(f"  [WARN] {name}: not found in HyCal map — skipping",
                  file=sys.stderr)
            missing += 1
            continue
        x, y, sx, sy = geo
        xs.append(x)
        ys.append(y)
        sxs.append(sx)
        sys_.append(sy)
        vals.append(ch[param_key])

    n_plotted = len(vals)
    if n_plotted == 0:
        print(f"  [WARN] {material}/{param_name}: no modules to plot after "
              f"geometry matching — skipping figure.", file=sys.stderr)
        return

    xs_arr   = np.asarray(xs,   dtype=np.float64)
    ys_arr   = np.asarray(ys,   dtype=np.float64)
    sxs_arr  = np.asarray(sxs,  dtype=np.float64)
    sys_arr  = np.asarray(sys_, dtype=np.float64)
    vals_arr = np.asarray(vals, dtype=np.float64)

    vmin, vmax = _resolve_clim(vals_arr, vmin_override, vmax_override)

    # --- Build rectangle patches -------------------------------------------
    patches = [
        Rectangle((x - sx / 2.0, y - sy / 2.0), sx, sy)
        for x, y, sx, sy in zip(xs_arr, ys_arr, sxs_arr, sys_arr)
    ]
    pc = PatchCollection(patches, cmap=meta["cmap"], linewidths=0)
    pc.set_array(vals_arr)
    pc.set_clim(vmin, vmax)

    # --- Figure layout -------------------------------------------------------
    # Determine axis limits from the geometry extent, with ~7 % padding.
    x_lo = float((xs_arr - sxs_arr / 2.0).min())
    x_hi = float((xs_arr + sxs_arr / 2.0).max())
    y_lo = float((ys_arr - sys_arr / 2.0).min())
    y_hi = float((ys_arr + sys_arr / 2.0).max())
    pad_x = 0.07 * (x_hi - x_lo)
    pad_y = 0.07 * (y_hi - y_lo)

    fig, ax = plt.subplots(figsize=(9, 9))
    ax.add_collection(pc)

    cbar = fig.colorbar(pc, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label(meta["label"], fontsize=11)

    ax.set_xlim(x_lo - pad_x, x_hi + pad_x)
    ax.set_ylim(y_lo - pad_y, y_hi + pad_y)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)", fontsize=11)
    ax.set_ylabel("y (mm)", fontsize=11)

    med_val = float(np.median(vals_arr))
    title = (
        f"{material}  —  {meta['short']} 2D spatial map\n"
        f"{n_plotted} modules plotted"
        + (f"  ({missing} missing from map)" if missing else "")
        + f"   colour range [{vmin:.2f}, {vmax:.2f}] ns"
        + f"   median = {med_val:.2f} ns"
    )
    ax.set_title(title, fontsize=10)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Wrote {out_path}")


# ---------------------------------------------------------------------------
# n_pulses 2D map (log colour scale)
# ---------------------------------------------------------------------------

def make_2d_map_n_pulses(
    material: str,
    channels: list[dict],
    geo_map: dict[str, tuple[float, float, float, float]],
    out_path: Path,
    vmin_override: float | None,
    vmax_override: float | None,
) -> None:
    """Render a 2D spatial map colour-coded by n_pulses_used and save to *out_path*.

    The colour scale is logarithmic (``LogNorm``) because pulse counts span
    orders of magnitude across the detector face (thousands near the beam axis,
    hundreds or fewer at the periphery).  The lower bound is clamped to at
    least 1 so that LogNorm never receives a non-positive value.

    Parameters
    ----------
    material : str
        Module material label (e.g. 'PbWO4').
    channels : list of dict
        Channel records for *material* (already filtered for min_pulses).
        Each record must contain the key 'n_pulses_used'.
    geo_map : dict
        name → (x, y, sx, sy) from load_hycal_map().
    out_path : Path
        Destination PNG path.
    vmin_override, vmax_override : float or None
        CLI overrides for the colour-scale bounds (--n-vmin / --n-vmax).
        When None the 2nd and 98th percentiles of the plotted counts are used,
        with vmin clamped to ≥ 1.
    """
    cmap = "viridis"
    cbar_label = "n_pulses_used per channel"

    # --- Match channels to geometry ----------------------------------------
    xs:   list[float] = []
    ys:   list[float] = []
    sxs:  list[float] = []
    sys_: list[float] = []
    vals: list[float] = []
    missing = 0

    for ch in channels:
        name = ch["name"]
        geo = geo_map.get(name)
        if geo is None:
            print(f"  [WARN] {name}: not found in HyCal map — skipping",
                  file=sys.stderr)
            missing += 1
            continue
        x, y, sx, sy = geo
        xs.append(x)
        ys.append(y)
        sxs.append(sx)
        sys_.append(sy)
        vals.append(float(ch["n_pulses_used"]))

    n_plotted = len(vals)
    if n_plotted == 0:
        print(f"  [WARN] {material}/n_pulses: no modules to plot after "
              f"geometry matching — skipping figure.", file=sys.stderr)
        return

    xs_arr   = np.asarray(xs,   dtype=np.float64)
    ys_arr   = np.asarray(ys,   dtype=np.float64)
    sxs_arr  = np.asarray(sxs,  dtype=np.float64)
    sys_arr  = np.asarray(sys_, dtype=np.float64)
    vals_arr = np.asarray(vals, dtype=np.float64)

    # Colour-scale bounds — default to 2nd/98th percentile, but clamp vmin ≥ 1
    if vmin_override is None:
        vmin = max(1.0, float(np.percentile(vals_arr, 2)))
    else:
        vmin = max(1.0, float(vmin_override))

    if vmax_override is None:
        vmax = float(np.percentile(vals_arr, 98))
    else:
        vmax = float(vmax_override)

    if vmin >= vmax:
        # Degenerate range — fall back to full range (still clamped ≥ 1)
        vmin = max(1.0, float(vals_arr.min()))
        vmax = float(vals_arr.max())

    norm = LogNorm(vmin=vmin, vmax=vmax)

    # --- Build rectangle patches -------------------------------------------
    patches = [
        Rectangle((x - sx / 2.0, y - sy / 2.0), sx, sy)
        for x, y, sx, sy in zip(xs_arr, ys_arr, sxs_arr, sys_arr)
    ]
    pc = PatchCollection(patches, cmap=cmap, norm=norm, linewidths=0)
    pc.set_array(vals_arr)

    # --- Figure layout -------------------------------------------------------
    x_lo = float((xs_arr - sxs_arr / 2.0).min())
    x_hi = float((xs_arr + sxs_arr / 2.0).max())
    y_lo = float((ys_arr - sys_arr / 2.0).min())
    y_hi = float((ys_arr + sys_arr / 2.0).max())
    pad_x = 0.07 * (x_hi - x_lo)
    pad_y = 0.07 * (y_hi - y_lo)

    fig, ax = plt.subplots(figsize=(9, 9))
    ax.add_collection(pc)

    cbar = fig.colorbar(pc, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label(cbar_label, fontsize=11)

    ax.set_xlim(x_lo - pad_x, x_hi + pad_x)
    ax.set_ylim(y_lo - pad_y, y_hi + pad_y)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)", fontsize=11)
    ax.set_ylabel("y (mm)", fontsize=11)

    med_val = float(np.median(vals_arr))
    title = (
        f"{material}  —  n_pulses_used 2D spatial map\n"
        f"{n_plotted} modules plotted"
        + (f"  ({missing} missing from map)" if missing else "")
        + f"   colour range [{vmin:.0f}, {vmax:.0f}] (log scale)"
        + f"   median = {med_val:.0f}"
    )
    ax.set_title(title, fontsize=10)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Wrote {out_path}")


# ---------------------------------------------------------------------------
# Per-material summary
# ---------------------------------------------------------------------------

def print_material_summary(
    material: str,
    channels: list[dict],
    n_map_modules: int,
    min_pulses: int,
) -> None:
    """Print a concise per-material parameter summary to stdout."""
    tau_r_arr = np.array([c["tau_r"] for c in channels])
    tau_f_arr = np.array([c["tau_f"] for c in channels])
    t0_arr    = np.array([c["t0"]    for c in channels])
    n_pulses_arr = np.array([c["n_pulses_used"] for c in channels], dtype=np.int64)
    n = len(channels)
    print(
        f"\n=== {material}: {n}/{n_map_modules} map modules had templates"
        f" (min_pulses={min_pulses}) ==="
    )
    print(f"    τ_r range: {tau_r_arr.min():.1f} ns to {tau_r_arr.max():.1f} ns"
          f" (median {np.median(tau_r_arr):.1f} ns)")
    print(f"    τ_f range: {tau_f_arr.min():.1f} ns to {tau_f_arr.max():.1f} ns"
          f" (median {np.median(tau_f_arr):.1f} ns)")
    print(f"    t₀ range: {t0_arr.min():.1f} ns to {t0_arr.max():.1f} ns"
          f" (median {np.median(t0_arr):.1f} ns)")
    print(f"    n_pulses: min = {n_pulses_arr.min()}, max = {n_pulses_arr.max()},"
          f" median {int(np.median(n_pulses_arr))}, total {n_pulses_arr.sum()}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _default_hcmap_path(script_dir: Path) -> Path:
    """Resolve the default hycal_map.json path.

    Prefers $PRAD2_DATABASE_DIR/hycal_map.json if the env-var is set,
    otherwise falls back to <script_dir>/../../database/hycal_map.json.
    """
    env_dir = os.environ.get("PRAD2_DATABASE_DIR", "")
    if env_dir:
        return Path(env_dir).resolve() / "hycal_map.json"
    return (script_dir / ".." / ".." / "database" / "hycal_map.json").resolve()


def build_parser() -> argparse.ArgumentParser:
    script_dir = Path(__file__).resolve().parent
    default_hcmap = _default_hcmap_path(script_dir)

    p = argparse.ArgumentParser(
        description=(
            "Produce 2D spatial maps of HyCal pulse-template parameters "
            "(τ_r, τ_f, t₀) on the physical detector face to test whether "
            "observed multi-modal PbWO4 distributions organize spatially "
            "(radiation-damage / crystal-batch clustering) or are randomly "
            "distributed (per-PMT / electronics variation)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "template_json",
        type=Path,
        help="Path to the pulse_templates_*.json produced by fit_pulse_template.py.",
    )
    p.add_argument(
        "--hc-map-file",
        type=Path,
        default=default_hcmap,
        metavar="PATH",
        help="Path to hycal_map.json. Defaults to $PRAD2_DATABASE_DIR/hycal_map.json "
             "if the env-var is set, else <repo>/database/hycal_map.json.",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        metavar="DIR",
        help="Directory for output PNG files. Defaults to the directory of template_json.",
    )
    p.add_argument(
        "--min-pulses",
        type=int,
        default=50,
        metavar="N",
        help="Skip channels with n_pulses_used < N.",
    )
    p.add_argument(
        "--materials",
        type=str,
        default="PbWO4,PbGlass",
        metavar="MAT1,MAT2,...",
        help="Comma-separated list of module_type values to plot.",
    )
    p.add_argument(
        "--vmin",
        type=float,
        default=None,
        metavar="VAL",
        help="Override colour-scale minimum (ns) for the τ_r/τ_f/t₀ maps. "
             "Default: 2nd percentile of the plotted values.",
    )
    p.add_argument(
        "--vmax",
        type=float,
        default=None,
        metavar="VAL",
        help="Override colour-scale maximum (ns) for the τ_r/τ_f/t₀ maps. "
             "Default: 98th percentile of the plotted values.",
    )
    p.add_argument(
        "--n-vmin",
        type=float,
        default=None,
        metavar="VAL",
        help="Override colour-scale minimum for the n_pulses_used map "
             "(log scale; clamped to ≥ 1). "
             "Default: 2nd percentile of the plotted counts (clamped to ≥ 1).",
    )
    p.add_argument(
        "--n-vmax",
        type=float,
        default=None,
        metavar="VAL",
        help="Override colour-scale maximum for the n_pulses_used map (log scale). "
             "Default: 98th percentile of the plotted counts.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    # --- Validate inputs -----------------------------------------------------
    json_path: Path = args.template_json.resolve()
    if not json_path.exists():
        print(f"ERROR: template JSON not found: {json_path}", file=sys.stderr)
        return 1

    hcmap_path: Path = args.hc_map_file.resolve()
    if not hcmap_path.exists():
        print(f"ERROR: HyCal map not found: {hcmap_path}", file=sys.stderr)
        return 1

    out_dir: Path = (
        args.out_dir.resolve()
        if args.out_dir is not None
        else json_path.parent
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    min_pulses: int = args.min_pulses
    materials: list[str] = [
        m.strip() for m in args.materials.split(",") if m.strip()
    ]

    # --- Load inputs ---------------------------------------------------------
    print(f"Loading HyCal map:  {hcmap_path}")
    geo_map = load_hycal_map(hcmap_path)
    print(f"  {len(geo_map)} modules loaded from map.")

    print(f"\nLoading templates:  {json_path}")
    by_material = load_template_channels(json_path, min_pulses, materials)

    # --- Count map modules per material for the summary ----------------------
    map_modules_by_material: dict[str, int] = {}
    for name, (x, y, sx, sy) in geo_map.items():
        # Identify material by name prefix: 'W' → PbWO4, 'G' → PbGlass.
        # We don't store type in geo_map so we use a simple heuristic based
        # on the prefix that matches the naming convention in hycal_map.json.
        # A more robust approach uses the 't' field — we re-read it below.
        pass

    # Re-read the map JSON once more to get a per-module type count.
    with hcmap_path.open(encoding="utf-8") as fh:
        hcmap_records = json.load(fh)
    n_map_by_type: dict[str, int] = {}
    for rec in hcmap_records:
        t = rec.get("t", "Unknown")
        n_map_by_type[t] = n_map_by_type.get(t, 0) + 1

    # --- Process each material -----------------------------------------------
    for material in materials:
        channels = by_material.get(material, [])
        if not channels:
            print(
                f"\nWARNING: material '{material}' has no channels in the "
                f"template JSON that pass min_pulses={min_pulses}. Skipping."
            )
            continue

        # Quick check: how many of these channels appear in the map?
        n_in_map = sum(1 for ch in channels if ch["name"] in geo_map)
        if n_in_map == 0:
            print(
                f"\nWARNING: material '{material}' — none of the {len(channels)} "
                f"template channels were found in the HyCal geometry map. Skipping."
            )
            continue

        print(f"\nProcessing {material}: {len(channels)} channels "
              f"({n_in_map} found in map) …")

        for param_name in ("tau_r", "tau_f", "t0"):
            out_png = out_dir / f"template_2d_{material}_{param_name}.png"
            make_2d_map(
                material=material,
                param_name=param_name,
                channels=channels,
                geo_map=geo_map,
                out_path=out_png,
                vmin_override=args.vmin,
                vmax_override=args.vmax,
            )

        out_png_n = out_dir / f"template_2d_{material}_n_pulses.png"
        make_2d_map_n_pulses(
            material=material,
            channels=channels,
            geo_map=geo_map,
            out_path=out_png_n,
            vmin_override=args.n_vmin,
            vmax_override=args.n_vmax,
        )

        # Summary
        n_map_total = n_map_by_type.get(material, len(channels))
        print_material_summary(material, channels, n_map_total, min_pulses)

    return 0


if __name__ == "__main__":
    sys.exit(main())
