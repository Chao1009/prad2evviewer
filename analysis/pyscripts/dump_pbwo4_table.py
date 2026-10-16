"""dump_pbwo4_table.py — Extract per-module PbWO4 tau_r / tau_f statistics
from two pulse-template JSON files and write a comparison table.

Default inputs
--------------
  output/pulse_templates_025308_h50_RawSum_max20modules.json   → run #25308
  output/pulse_templates_025320_SSP_RawSum_h50_max20modules.json → run #25320

Default output
--------------
  output/pulse_templates_h50.dat  (tab-separated, human-readable)

Columns
-------
  module   tau_r_mean_25308  tau_r_sigma_25308  tau_f_mean_25308  tau_f_sigma_25308
           tau_r_mean_25320  tau_r_sigma_25320  tau_f_mean_25320  tau_f_sigma_25320
           d_tau_r/(tau_r_25308)  d_tau_f/(tau_f_25308)

"mean" maps to the JSON key "median"; "sigma" maps to "mad".
Modules present in only one file get NaN for missing values.
The last two columns are relative differences: (tau_A - tau_B) / tau_A.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import re
import sys

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_SCRIPT_DIR = pathlib.Path(__file__).parent
_OUTPUT_DIR = _SCRIPT_DIR / "output"

_DEFAULT_FILE_25308 = _OUTPUT_DIR / "pulse_templates_025308_h50_RawSum_max20modules.json"
_DEFAULT_FILE_25320 = _OUTPUT_DIR / "pulse_templates_025320_SSP_RawSum_h50_max20modules.json"
_DEFAULT_OUTPUT    = _OUTPUT_DIR / "pulse_templates_h50.dat"

_NAN = float("nan")

# Keys extracted from _meta for the header comment block
_META_KEYS = [
    "model",
    "chi2_max",
    "height_min",
    "height_rms_mult",
    "min_pulses",
    "t0_min",
    "t0_max",
    "model_err_floor",
    "trigger_event_type_spec",
    "trigger_event_type_accept",
    "clk_ns",
    "pre_samples",
    "post_samples",
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _natural_key(name: str):
    """Sort key that orders 'W1' < 'W2' < ... < 'W10' < 'W11' ..."""
    m = re.match(r'^([A-Za-z]+)(\d+)$', name)
    return (m.group(1), int(m.group(2))) if m else (name, 0)


def _load_json(path: pathlib.Path) -> dict:
    """Load and return the full JSON document."""
    try:
        with path.open() as fh:
            return json.load(fh)
    except Exception as exc:
        sys.exit(f"ERROR reading {path}: {exc}")


def _load_pbwo4(data: dict) -> dict[str, dict]:
    """Return only the PbWO4 module entries (keyed by module name)."""
    return {
        k: v
        for k, v in data.items()
        if isinstance(v, dict) and v.get("module_type") == "PbWO4"
    }


def _get(mod: dict | None, param: str, stat: str) -> float:
    """Return mod[param][stat] as float, or NaN if absent."""
    if mod is None:
        return _NAN
    try:
        return float(mod[param][stat])
    except (KeyError, TypeError, ValueError):
        return _NAN


def _fmt(val: float, width: int = 10, prec: int = 6) -> str:
    """Format a float value; use 'NaN' for NaN."""
    if math.isnan(val):
        return "NaN".rjust(width)
    return f"{val:{width}.{prec}f}"


def _fmt_meta_value(val) -> str:
    """Format a metadata value as a compact string."""
    if val is None:
        return "None"
    if isinstance(val, list):
        return "[" + ", ".join(str(x) for x in val) + "]"
    return str(val)


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def build_table(
    mods_a: dict[str, dict],
    mods_b: dict[str, dict],
    label_a: str,
    label_b: str,
) -> list[dict]:
    """Return a list of row dicts (one per PbWO4 module, union of both files).

    The last two fields are relative differences:
        delta_tr = (tau_r_a - tau_r_b) / tau_r_a
        delta_tf = (tau_f_a - tau_f_b) / tau_f_a
    NaN is used when either value is missing or tau_a is zero.
    """
    all_modules = sorted(set(mods_a) | set(mods_b), key=_natural_key)

    rows: list[dict] = []
    for name in all_modules:
        mod_a = mods_a.get(name)
        mod_b = mods_b.get(name)

        tr_mean_a  = _get(mod_a, "tau_r_ns", "median")
        tr_sigma_a = _get(mod_a, "tau_r_ns", "mad")
        tf_mean_a  = _get(mod_a, "tau_f_ns", "median")
        tf_sigma_a = _get(mod_a, "tau_f_ns", "mad")

        tr_mean_b  = _get(mod_b, "tau_r_ns", "median")
        tr_sigma_b = _get(mod_b, "tau_r_ns", "mad")
        tf_mean_b  = _get(mod_b, "tau_f_ns", "median")
        tf_sigma_b = _get(mod_b, "tau_f_ns", "mad")

        # Relative differences: (a - b) / a  — NaN on missing or zero denominator
        if math.isnan(tr_mean_a) or math.isnan(tr_mean_b):
            delta_tr = _NAN
        elif tr_mean_a == 0.0:
            delta_tr = _NAN
        else:
            delta_tr = (tr_mean_a - tr_mean_b) / tr_mean_a

        if math.isnan(tf_mean_a) or math.isnan(tf_mean_b):
            delta_tf = _NAN
        elif tf_mean_a == 0.0:
            delta_tf = _NAN
        else:
            delta_tf = (tf_mean_a - tf_mean_b) / tf_mean_a

        rows.append({
            "module":        name,
            "tr_mean_a":     tr_mean_a,
            "tr_sigma_a":    tr_sigma_a,
            "tf_mean_a":     tf_mean_a,
            "tf_sigma_a":    tf_sigma_a,
            "tr_mean_b":     tr_mean_b,
            "tr_sigma_b":    tr_sigma_b,
            "tf_mean_b":     tf_mean_b,
            "tf_sigma_b":    tf_sigma_b,
            "delta_tr":      delta_tr,
            "delta_tf":      delta_tf,
        })
    return rows


def _meta_comment_block(meta_a: dict, meta_b: dict, label_a: str, label_b: str) -> str:
    """Build a multi-line comment string listing key metadata/cuts from both files."""
    run_a = label_a.lstrip("#")
    run_b = label_b.lstrip("#")

    lines = ["#", "# --- Fit parameters and cuts ---"]
    lines.append(f"#  {'Parameter':<28}  {('Run ' + run_a):<22}  {'Run ' + run_b}")
    lines.append(f"#  {'-'*28}  {'-'*22}  {'-'*22}")

    for key in _META_KEYS:
        val_a = _fmt_meta_value(meta_a.get(key, "N/A"))
        val_b = _fmt_meta_value(meta_b.get(key, "N/A"))
        lines.append(f"#  {key:<28}  {val_a:<22}  {val_b}")

    lines.append("#")
    return "\n".join(lines) + "\n"


def write_table(
    rows: list[dict],
    out_path: pathlib.Path,
    label_a: str,
    label_b: str,
    file_a: pathlib.Path,
    file_b: pathlib.Path,
    meta_a: dict,
    meta_b: dict,
) -> None:
    """Write the comparison table to *out_path* as a tab-separated file."""

    # Column widths
    W_MOD  = 8   # module name
    W_VAL  = 12  # numeric columns
    PREC   = 6   # decimal places

    sep = "\t"

    def _hcol(text: str, w: int) -> str:
        return text.ljust(w)

    run_a = label_a.lstrip("#")
    run_b = label_b.lstrip("#")

    header1_parts = [
        _hcol("# module", W_MOD),
        _hcol(f"tau_r_mean_{run_a}", W_VAL),
        _hcol(f"tau_r_sigma_{run_a}", W_VAL),
        _hcol(f"tau_f_mean_{run_a}", W_VAL),
        _hcol(f"tau_f_sigma_{run_a}", W_VAL),
        _hcol(f"tau_r_mean_{run_b}", W_VAL),
        _hcol(f"tau_r_sigma_{run_b}", W_VAL),
        _hcol(f"tau_f_mean_{run_b}", W_VAL),
        _hcol(f"tau_f_sigma_{run_b}", W_VAL),
        _hcol(f"d_tau_r/(tau_r_{run_a})", W_VAL + 4),
        _hcol(f"d_tau_f/(tau_f_{run_a})", W_VAL + 4),
    ]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        # Preamble comments
        fh.write(f"# PbWO4 per-module tau_r / tau_f comparison\n")
        fh.write(f"# Run {label_a}: {file_a}\n")
        fh.write(f"# Run {label_b}: {file_b}\n")
        fh.write(f"# 'mean' = median, 'sigma' = MAD  (from per-module fit statistics)\n")
        fh.write(f"# All times in ns.\n")
        fh.write(f"# Modules present in only one run have NaN for missing values.\n")
        fh.write(f"# Last two columns: relative difference = (tau_{run_a} - tau_{run_b}) / tau_{run_a}\n")
        fh.write(_meta_comment_block(meta_a, meta_b, label_a, label_b))
        fh.write(sep.join(header1_parts) + "\n")

        for row in rows:
            parts = [
                row["module"].ljust(W_MOD),
                _fmt(row["tr_mean_a"],  W_VAL, PREC),
                _fmt(row["tr_sigma_a"], W_VAL, PREC),
                _fmt(row["tf_mean_a"],  W_VAL, PREC),
                _fmt(row["tf_sigma_a"], W_VAL, PREC),
                _fmt(row["tr_mean_b"],  W_VAL, PREC),
                _fmt(row["tr_sigma_b"], W_VAL, PREC),
                _fmt(row["tf_mean_b"],  W_VAL, PREC),
                _fmt(row["tf_sigma_b"], W_VAL, PREC),
                _fmt(row["delta_tr"],   W_VAL + 4, PREC),
                _fmt(row["delta_tf"],   W_VAL + 4, PREC),
            ]
            fh.write(sep.join(parts) + "\n")

    print(f"Table written → {out_path}  ({len(rows)} PbWO4 modules)")


# ---------------------------------------------------------------------------
# Summary statistics to stdout
# ---------------------------------------------------------------------------

def print_summary(rows: list[dict], label_a: str, label_b: str) -> None:
    """Print a brief summary of relative-delta distributions to stdout."""
    import statistics

    deltas_tr = [r["delta_tr"] for r in rows if not math.isnan(r["delta_tr"])]
    deltas_tf = [r["delta_tf"] for r in rows if not math.isnan(r["delta_tf"])]

    n_only_a = sum(1 for r in rows if math.isnan(r["tr_mean_b"]))
    n_only_b = sum(1 for r in rows if math.isnan(r["tr_mean_a"]))
    n_both   = sum(1 for r in rows if not (math.isnan(r["tr_mean_a"]) or math.isnan(r["tr_mean_b"])))

    print(f"\n{'='*60}")
    print(f"  PbWO4 per-module comparison: {label_a} vs {label_b}")
    print(f"{'='*60}")
    print(f"  Total modules in table : {len(rows)}")
    print(f"  In both runs           : {n_both}")
    print(f"  Only in {label_a:<8}    : {n_only_a}")
    print(f"  Only in {label_b:<8}    : {n_only_b}")

    if deltas_tr:
        tr_med  = statistics.median(deltas_tr)
        tr_mean = statistics.mean(deltas_tr)
        tr_std  = statistics.stdev(deltas_tr) if len(deltas_tr) > 1 else _NAN
        print(f"\n  relative Δτ_r = (τ_r({label_a}) − τ_r({label_b})) / τ_r({label_a})")
        print(f"    median = {tr_med:+.6f}   mean = {tr_mean:+.6f}   σ = {tr_std:.6f}")

    if deltas_tf:
        tf_med  = statistics.median(deltas_tf)
        tf_mean = statistics.mean(deltas_tf)
        tf_std  = statistics.stdev(deltas_tf) if len(deltas_tf) > 1 else _NAN
        print(f"\n  relative Δτ_f = (τ_f({label_a}) − τ_f({label_b})) / τ_f({label_a})")
        print(f"    median = {tf_med:+.6f}   mean = {tf_mean:+.6f}   σ = {tf_std:.6f}")

    print(f"{'='*60}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Compare per-module PbWO4 tau_r/tau_f between two runs.")
    ap.add_argument(
        "--file-a",
        default=str(_DEFAULT_FILE_25308),
        help="JSON template file for run A (default: 025308 h50 RawSum max20modules).",
    )
    ap.add_argument(
        "--file-b",
        default=str(_DEFAULT_FILE_25320),
        help="JSON template file for run B (default: 025320 SSP_RawSum h50 max20modules).",
    )
    ap.add_argument(
        "--label-a",
        default="#25308",
        help="Label for run A used in column headers (default: '#25308').",
    )
    ap.add_argument(
        "--label-b",
        default="#25320",
        help="Label for run B used in column headers (default: '#25320').",
    )
    ap.add_argument(
        "--out",
        default=str(_DEFAULT_OUTPUT),
        help="Output .dat file path (default: output/pulse_templates_h50.dat).",
    )
    ap.add_argument(
        "--no-summary",
        action="store_true",
        help="Suppress the summary statistics printed to stdout.",
    )
    return ap


def main(argv=None) -> None:
    args = _build_parser().parse_args(argv)

    file_a  = pathlib.Path(args.file_a).expanduser().resolve()
    file_b  = pathlib.Path(args.file_b).expanduser().resolve()
    out_path = pathlib.Path(args.out).expanduser().resolve()
    label_a = args.label_a
    label_b = args.label_b

    for p in (file_a, file_b):
        if not p.is_file():
            sys.exit(f"ERROR: file not found: {p}")

    print(f"Loading {label_a}: {file_a.name}")
    data_a = _load_json(file_a)
    meta_a = data_a.get("_meta", {})
    mods_a = _load_pbwo4(data_a)
    print(f"  → {len(mods_a)} PbWO4 modules")

    print(f"Loading {label_b}: {file_b.name}")
    data_b = _load_json(file_b)
    meta_b = data_b.get("_meta", {})
    mods_b = _load_pbwo4(data_b)
    print(f"  → {len(mods_b)} PbWO4 modules")

    rows = build_table(mods_a, mods_b, label_a, label_b)

    if not args.no_summary:
        print_summary(rows, label_a, label_b)

    write_table(rows, out_path, label_a, label_b, file_a, file_b, meta_a, meta_b)


if __name__ == "__main__":
    main()
