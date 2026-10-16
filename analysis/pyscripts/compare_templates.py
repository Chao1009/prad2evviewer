"""compare_templates.py — Compare per-material pulse-template parameters across
runs, height cuts, and t0 cuts.

Reads ``pulse_templates_<RUN>_<TRIGGER>_h<HMIN>_t0min<T0>.json`` files
produced by ``fit_pulse_template.py`` and prints side-by-side tables of τ_r,
τ_f, χ²/dof, and n_pulses for every (height_min, t0_min) combination, with
one column per run/energy.

Optionally generates four matplotlib figures:
  compare_tau_r_vs_energy.png   — τ_r vs. beam energy, one curve per cut combo
  compare_tau_f_vs_energy.png   — τ_f vs. beam energy, one curve per cut combo
  compare_tau_r_vs_hmin.png     — τ_r vs. height_min, one curve per energy
  compare_tau_r_vs_t0min.png    — τ_r vs. t0_min,    one curve per energy

No prad2py dependency — pure JSON + stdlib + numpy + (optionally) matplotlib.
"""

import argparse
import json
import pathlib
import sys

import numpy as np  # used for NaN sentinels; no heavy computation required

# ---------------------------------------------------------------------------
# Built-in run → energy mapping (MeV)
# ---------------------------------------------------------------------------
_RUN_ENERGY: dict[str, float] = {
    "024590": 728.918,
    "025320": 2239.51,
    "024512": 3488.43,
}


# ---------------------------------------------------------------------------
# File resolution
# ---------------------------------------------------------------------------

def _find_template(
    dir_path: pathlib.Path,
    run: str,
    trigger: str,
    hmin: str,
    t0: str,
) -> pathlib.Path | None:
    """Locate the JSON for one (run, trigger, hmin, t0) combination.

    The t0 value may be stored with 'p' as the decimal separator
    (e.g. ``t0min22p0``) or with a literal dot (e.g. ``t0min22.0``).
    Try the 'p' format first; fall back to the literal-dot format.
    """
    t0_p = str(t0).replace(".", "p")
    candidates = [
        dir_path / f"pulse_templates_{run}_{trigger}_h{hmin}_t0min{t0_p}.json",
        dir_path / f"pulse_templates_{run}_{trigger}_h{hmin}_t0min{t0}.json",
    ]
    for c in candidates:
        if c.is_file():
            return c
    return None


# ---------------------------------------------------------------------------
# JSON loading
# ---------------------------------------------------------------------------

def _load_material(path: pathlib.Path, material: str) -> dict | None:
    """Return the ``_by_type[material]`` block, or None on any error."""
    try:
        with path.open() as fh:
            data = json.load(fh)
        return data.get("_by_type", {}).get(material)
    except Exception as exc:  # noqa: BLE001
        print(f"  [warn] could not read {path.name}: {exc}", file=sys.stderr)
        return None


def _extract(block: dict | None, key: str, sub: str) -> float | None:
    """Extract block[key][sub] as float, or None if missing."""
    if block is None:
        return None
    try:
        return float(block[key][sub])
    except (KeyError, TypeError, ValueError):
        return None


def _extract_scalar(block: dict | None, key: str) -> float | None:
    """Extract block[key] as a plain scalar float."""
    if block is None:
        return None
    try:
        return float(block[key])
    except (KeyError, TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

def collect_results(
    dir_path: pathlib.Path,
    runs: list[str],
    trigger: str,
    hmins: list[str],
    t0s: list[str],
    material: str,
    energy_map: dict[str, float],
) -> dict:
    """Return nested dict: results[run][hmin][t0] = {tau_r, tau_f, chi2, n_pulses, path}."""
    results: dict = {}
    for run in runs:
        results[run] = {}
        for hmin in hmins:
            results[run][hmin] = {}
            for t0 in t0s:
                path = _find_template(dir_path, run, trigger, hmin, t0)
                if path is None:
                    results[run][hmin][t0] = None
                    continue
                block = _load_material(path, material)
                entry = {
                    "path": str(path),
                    "tau_r_med": _extract(block, "tau_r_ns", "median"),
                    "tau_r_mad": _extract(block, "tau_r_ns", "mad"),
                    "tau_f_med": _extract(block, "tau_f_ns", "median"),
                    "tau_f_mad": _extract(block, "tau_f_ns", "mad"),
                    "chi2_med":  _extract(block, "chi2_per_dof", "median"),
                    "chi2_mean": _extract(block, "chi2_per_dof", "mean"),
                    "n_pulses":  _extract_scalar(block, "n_pulses_total"),
                    "n_channels":      _extract_scalar(block, "n_channels"),
                    "n_channels_good": _extract_scalar(block, "n_channels_good"),
                }
                results[run][hmin][t0] = entry
    return results


# ---------------------------------------------------------------------------
# Table printing helpers
# ---------------------------------------------------------------------------

_MISSING = "      —      "


def _fmt_pm(med: float | None, mad: float | None) -> str:
    if med is None:
        return _MISSING
    if mad is None:
        return f"   {med:6.3f}       "
    return f" {med:6.3f} ± {mad:5.3f} "


def _fmt_chi2(med: float | None, mean: float | None) -> str:
    if med is None:
        return _MISSING
    if mean is None:
        return f"  med={med:5.3f}        "
    return f" med={med:5.3f}/mn={mean:5.3f} "


def _fmt_n(val: float | None) -> str:
    if val is None:
        return _MISSING
    return f"   {int(val):>10d}   "


def _col_header(run: str, energy_map: dict[str, float], width: int) -> str:
    energy = energy_map.get(run)
    if energy is not None:
        label = f"{run} ({energy/1000:.1f} GeV)"
    else:
        label = run
    return label.center(width)


# ---------------------------------------------------------------------------
# Print tables
# ---------------------------------------------------------------------------

COL_W = 20  # cell width in characters


def _print_table(
    title: str,
    runs: list[str],
    hmins: list[str],
    t0s: list[str],
    results: dict,
    energy_map: dict[str, float],
    cell_fn,          # callable(entry) -> str
) -> None:
    print(f"\n=== {title} ===")
    # header
    col_labels = [_col_header(r, energy_map, COL_W) for r in runs]
    header_run = " | ".join(col_labels)
    print(f"{'h_min':>6}  {'t0_min':>6} | {header_run}")
    sep_run = "-+-".join(["-" * COL_W] * len(runs))
    print(f"{'------':>6}  {'------':>6}-+-{sep_run}")

    for hmin in hmins:
        for t0 in t0s:
            cells = []
            for run in runs:
                entry = results.get(run, {}).get(hmin, {}).get(t0)
                cells.append(cell_fn(entry).center(COL_W))
            row = " | ".join(cells)
            print(f"{hmin:>6}  {t0:>6} | {row}")


def print_all_tables(
    runs: list[str],
    hmins: list[str],
    t0s: list[str],
    results: dict,
    energy_map: dict[str, float],
    material: str,
) -> None:
    _print_table(
        f"{material} tau_r (ns) [median ± MAD]",
        runs, hmins, t0s, results, energy_map,
        lambda e: _fmt_pm(
            e["tau_r_med"] if e else None,
            e["tau_r_mad"] if e else None,
        ),
    )
    _print_table(
        f"{material} tau_f (ns) [median ± MAD]",
        runs, hmins, t0s, results, energy_map,
        lambda e: _fmt_pm(
            e["tau_f_med"] if e else None,
            e["tau_f_mad"] if e else None,
        ),
    )
    _print_table(
        f"{material} chi2/dof [median / mean]",
        runs, hmins, t0s, results, energy_map,
        lambda e: _fmt_chi2(
            e["chi2_med"]  if e else None,
            e["chi2_mean"] if e else None,
        ),
    )
    _print_table(
        f"{material} n_pulses_total",
        runs, hmins, t0s, results, energy_map,
        lambda e: _fmt_n(e["n_pulses"] if e else None),
    )


# ---------------------------------------------------------------------------
# JSON output
# ---------------------------------------------------------------------------

def write_json(out_path: pathlib.Path, results: dict) -> None:
    """Serialise results to JSON (replace non-serialisable Nones gracefully)."""
    def _clean(obj):
        if isinstance(obj, dict):
            return {k: _clean(v) for k, v in obj.items()}
        if obj is None:
            return None
        if isinstance(obj, float) and (obj != obj):  # NaN
            return None
        return obj

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        json.dump(_clean(results), fh, indent=2)
    print(f"\nJSON written → {out_path}")


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def make_plots(
    results: dict,
    runs: list[str],
    hmins: list[str],
    t0s: list[str],
    energy_map: dict[str, float],
    material: str,
    out_dir: pathlib.Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    energies_mev = [energy_map.get(r, float("nan")) for r in runs]

    # ---- helper: collect (x, y, yerr) arrays over runs for one cut combo ----
    def _series_vs_energy(hmin, t0, key_med, key_mad):
        xs, ys, errs = [], [], []
        for run, e_mev in zip(runs, energies_mev):
            entry = results.get(run, {}).get(hmin, {}).get(t0)
            if entry is not None and entry.get(key_med) is not None:
                xs.append(e_mev)
                ys.append(entry[key_med])
                errs.append(entry.get(key_mad) or 0.0)
        return np.asarray(xs), np.asarray(ys), np.asarray(errs)

    # -----------------------------------------------------------------------
    # Plot 1 & 2: tau_r / tau_f vs. energy, one curve per (hmin, t0) combo
    # -----------------------------------------------------------------------
    for param, key_med, key_mad, ylabel, fname in [
        ("tau_r", "tau_r_med", "tau_r_mad", r"$\tau_r$ (ns)", "compare_tau_r_vs_energy.png"),
        ("tau_f", "tau_f_med", "tau_f_mad", r"$\tau_f$ (ns)", "compare_tau_f_vs_energy.png"),
    ]:
        fig, ax = plt.subplots(figsize=(8, 5))
        for hmin in hmins:
            for t0 in t0s:
                xs, ys, errs = _series_vs_energy(hmin, t0, key_med, key_mad)
                if len(xs) == 0:
                    continue
                label = f"h{hmin} t0={t0}"
                ax.errorbar(xs, ys, yerr=errs, marker="o", capsize=3, label=label)
        ax.set_xlabel("Beam energy (MeV)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{material} {param} vs. beam energy")
        ax.legend(fontsize=7, ncol=2)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        out_path = out_dir / fname
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"Plot → {out_path}")

    # -----------------------------------------------------------------------
    # Plot 3: tau_r vs. height_min, one curve per energy (fixed middle t0)
    # -----------------------------------------------------------------------
    t0_fixed = t0s[len(t0s) // 2]   # middle t0 value
    fig, ax = plt.subplots(figsize=(7, 5))
    for run, e_mev in zip(runs, energies_mev):
        xs, ys, errs = [], [], []
        for hmin in hmins:
            entry = results.get(run, {}).get(hmin, {}).get(t0_fixed)
            if entry is not None and entry.get("tau_r_med") is not None:
                xs.append(int(hmin))
                ys.append(entry["tau_r_med"])
                errs.append(entry.get("tau_r_mad") or 0.0)
        if xs:
            label = f"{run} ({e_mev:.0f} MeV)"
            ax.errorbar(xs, ys, yerr=errs, marker="o", capsize=3, label=label)
    ax.set_xlabel("Height min (ADC counts)")
    ax.set_ylabel(r"$\tau_r$ (ns)")
    ax.set_title(f"{material} tau_r vs. height cut  [t0_min = {t0_fixed}]")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_path = out_dir / "compare_tau_r_vs_hmin.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Plot → {out_path}")

    # -----------------------------------------------------------------------
    # Plot 4: tau_r vs. t0_min, one curve per energy (fixed height_min = 500
    # if available, else middle value)
    # -----------------------------------------------------------------------
    hmin_fixed = "500" if "500" in hmins else hmins[len(hmins) // 2]
    fig, ax = plt.subplots(figsize=(7, 5))
    for run, e_mev in zip(runs, energies_mev):
        xs, ys, errs = [], [], []
        for t0 in t0s:
            entry = results.get(run, {}).get(hmin_fixed, {}).get(t0)
            if entry is not None and entry.get("tau_r_med") is not None:
                xs.append(float(t0))
                ys.append(entry["tau_r_med"])
                errs.append(entry.get("tau_r_mad") or 0.0)
        if xs:
            label = f"{run} ({e_mev:.0f} MeV)"
            ax.errorbar(xs, ys, yerr=errs, marker="o", capsize=3, label=label)
    ax.set_xlabel("t0 min (ns)")
    ax.set_ylabel(r"$\tau_r$ (ns)")
    ax.set_title(f"{material} tau_r vs. t0 cut  [h_min = {hmin_fixed}]")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_path = out_dir / "compare_tau_r_vs_t0min.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Plot → {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Compare template tau_r/tau_f across runs, height cuts, t0 cuts.")
    ap.add_argument(
        "--dir", default="",
        help="Directory with pulse_templates_*.json (default: "
             "output/trigger_scan relative to this script).")
    ap.add_argument(
        "--runs", default="024590,025320,024512",
        help="Comma-separated run numbers to compare.")
    ap.add_argument(
        "--height-mins", default="100,300,500,700",
        help="Comma-separated height-min values.")
    ap.add_argument(
        "--t0-mins", default="20.0,22.0,24.0",
        help="Comma-separated t0-min values.")
    ap.add_argument(
        "--material", default="PbWO4",
        help="Module type to compare (default PbWO4).")
    ap.add_argument(
        "--trigger", default="SSP_RawSum",
        help="Trigger event type in the filename (default SSP_RawSum).")
    ap.add_argument(
        "--out", default="",
        help="Optional JSON output path.")
    ap.add_argument(
        "--plot", action="store_true",
        help="Produce tau_r/tau_f vs. energy plots (one curve per cut).")
    return ap


def main(argv=None) -> None:
    args = _build_parser().parse_args(argv)

    # Resolve directory
    script_dir = pathlib.Path(__file__).parent
    if args.dir:
        dir_path = pathlib.Path(args.dir).expanduser().resolve()
    else:
        dir_path = script_dir / "output" / "trigger_scan"

    if not dir_path.is_dir():
        sys.exit(f"ERROR: directory not found: {dir_path}")

    runs   = [r.strip() for r in args.runs.split(",") if r.strip()]
    hmins  = [h.strip() for h in args.height_mins.split(",") if h.strip()]
    t0s    = [t.strip() for t in args.t0_mins.split(",")    if t.strip()]
    material = args.material
    trigger  = args.trigger

    # Energy map — built-in with any run not in the table left as unknown
    energy_map = dict(_RUN_ENERGY)

    print(f"\nComparing {material} templates in: {dir_path}")
    print(f"  runs       : {runs}")
    print(f"  height_mins: {hmins}")
    print(f"  t0_mins    : {t0s}")
    print(f"  trigger    : {trigger}")

    results = collect_results(dir_path, runs, trigger, hmins, t0s, material, energy_map)

    print_all_tables(runs, hmins, t0s, results, energy_map, material)

    if args.out:
        write_json(pathlib.Path(args.out).expanduser().resolve(), results)

    if args.plot:
        make_plots(results, runs, hmins, t0s, energy_map, material, dir_path)


if __name__ == "__main__":
    main()
