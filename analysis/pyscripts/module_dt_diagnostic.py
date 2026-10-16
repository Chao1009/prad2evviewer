#!/usr/bin/env python3
"""module_dt_diagnostic.py — per-module ΔT diagnostic for HyCal timing.

Consumes the TSV/CSV produced by study_hycal_timing.py and, for each
neighbour module, fits a Gaussian + flat-background model to the ΔT
distribution, producing:

  module_dt_summary_<run>[_<tag>].csv            — one row per module
  module_dt_heatmap_<run>[_<tag>].{png,pdf}      — 2-D (module_id × ΔT) heatmap
  module_dt_fits_<run>[_<tag>].pdf               — per-module histogram + fit pages
  module_dt_flags_<run>[_<tag>].md               — Markdown flag report
  hycal_secondary_peak_map_<run>[_<tag>].{png,pdf} — HyCal spatial map (§14.c)

Dependencies
------------
  Required : numpy, matplotlib (incl. matplotlib.backends.backend_pdf)
  Soft     : scipy.optimize.curve_fit  (plots + sideband counts produced even
             without scipy; all fit columns become NaN / "scipy-missing")
  Stdlib   : argparse, csv, dataclasses, datetime, hashlib, json, math,
             os, pathlib, re, sys, tempfile, time, typing

No pandas. No C++ build required. Runs standalone on any machine that has
the TSV and the above packages.

Fit convention (§7.a of the plan)
-----------------------------------
  Histogram ΔT into counts n_i, then convert to density y_i = n_i / dt_bin
  (counts/ns) before fitting.  Model:
    f(x) = A · exp(-½((x-μ)/σ)²) + B    [A, B in counts/ns]
  Signal integral: S = A · σ · √(2π) · erf(3/√2)
  Background in ±3σ: B_w = B · 6σ
  Significance: S / √max(B_w, 1)

See docs/plans/plan-module-dt-diagnostic.md for full design rationale.

Follow-ups (deferred per plan §12)
------------------------------------
  - Singles-rate accidentals model (R_seed · R_module · 2Δt).
  - Automated writing of a new hycal_time_offsets JSON.
  - Per-(seed, neighbour) 2-D μ map.
  - Cross-run trending.
  - Per-module-type thresholds (PbWO4 vs. PbGlass).
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import datetime
import hashlib
import json
import math
import os
import pathlib
import re
import sys
import tempfile
import time
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm

# ---------------------------------------------------------------------------
# Scipy — soft dependency (§5 of the plan)
# ---------------------------------------------------------------------------
try:
    from scipy.optimize import curve_fit
    from scipy.signal import find_peaks
    HAVE_SCIPY = True
except ImportError:  # pragma: no cover
    curve_fit = None   # type: ignore[assignment]
    find_peaks = None  # type: ignore[assignment]
    HAVE_SCIPY = False

# ---------------------------------------------------------------------------
# Constants (§2, §7.b)
# ---------------------------------------------------------------------------

TSV_COLUMNS: List[str] = [
    "event_num", "trigger_bits",
    "seed_module", "seed_id", "seed_time", "seed_energy",
    "neighbor_module", "neighbor_id", "neighbor_time", "neighbor_energy",
    "dt", "dx_q", "dy_q", "dist_q",
]

REQUIRED_COLUMNS: Set[str] = {
    "seed_id", "neighbor_id", "seed_time", "neighbor_time",
    "dt", "seed_energy", "neighbor_energy", "dist_q",
}

_NaN = float("nan")
_ERF3_SQRT2 = 0.9973002039  # erf(3/√2) — precomputed

# ---------------------------------------------------------------------------
# Dataclasses (§7.b)
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class ModuleStats:
    module_id: int
    n_entries: int
    n_total_seen: int
    mu: float = _NaN
    mu_err: float = _NaN
    sigma: float = _NaN
    sigma_err: float = _NaN
    A_per_ns: float = _NaN
    B_per_ns: float = _NaN
    signal_counts: float = _NaN
    significance: float = _NaN
    fit_status: str = "low-stats"
    sideband_bg_per_ns: float = _NaN
    sideband_counts: float = _NaN
    # §14.d secondary-peak fields
    secondary_mu: float = _NaN
    secondary_mu_err: float = _NaN
    secondary_sigma: float = _NaN
    secondary_sigma_err: float = _NaN
    secondary_amp_per_ns: float = _NaN
    secondary_significance: float = _NaN
    flags: dataclasses.field(default_factory=list) = dataclasses.field(
        default_factory=list
    )


@dataclasses.dataclass
class RunContext:
    input_path: pathlib.Path
    output_dir: pathlib.Path
    run_id: str
    tag: str
    dt_lo: float
    dt_hi: float
    dt_bin: float
    fit_window: float
    sideband_lo: float
    sideband_hi: float
    sig_thr: float
    mu_thr: float
    sigma_thr: float
    min_entries: int
    module_id_filter: Optional[Set[int]]
    min_seed_energy: float
    max_dist_q: float
    time_calib_path: Optional[pathlib.Path]
    time_calib_sha1: str
    offsets: Dict[int, float]
    default_offset: float
    scipy_available: bool
    chunk_rows: int
    progress_every: int
    max_rows: int
    verbose: bool
    no_pdf: bool
    # §14.d secondary-peak / spatial-map fields
    secondary_lo: float = 5.0
    secondary_hi: float = 25.0
    secondary_sig_thr: float = 3.0
    hycal_geometry_path: Optional[pathlib.Path] = None
    spatial_map_enabled: bool = True
    dt_edges: np.ndarray = dataclasses.field(default=None)        # type: ignore[assignment]
    sideband_edges: np.ndarray = dataclasses.field(default=None)  # type: ignore[assignment]

    def __post_init__(self) -> None:
        self.dt_edges = np.arange(
            self.dt_lo, self.dt_hi + self.dt_bin, self.dt_bin
        )
        self.sideband_edges = np.arange(
            self.sideband_lo, self.sideband_hi + self.dt_bin, self.dt_bin
        )


@dataclasses.dataclass
class ModuleAccumulator:
    module_id: int
    dt_hist: np.ndarray        # int64, len = N_dt_bins
    sideband_hist: np.ndarray  # int64, len = N_sideband_bins
    n_entries: int = 0         # post-filter, post-dt-window, post-offset
    n_total_seen: int = 0      # post-filter (excl. dt-window), pre-offset

    @classmethod
    def empty(cls, module_id: int, ctx: RunContext) -> "ModuleAccumulator":
        return cls(
            module_id=module_id,
            dt_hist=np.zeros(len(ctx.dt_edges) - 1, dtype=np.int64),
            sideband_hist=np.zeros(len(ctx.sideband_edges) - 1, dtype=np.int64),
        )


@dataclasses.dataclass
class ModuleGeom:
    """HyCal module geometry entry (§14.b)."""
    module_id:   int    # PrimEx id (G: 1-576, W: 1001-2152)
    name:        str    # "G<n>" or "W<n>"
    x_mm:        float  # module centre x on HyCal face (mm)
    y_mm:        float  # module centre y on HyCal face (mm)
    size_x_mm:   float  # module face width  (mm) — "sx" in JSON
    size_y_mm:   float  # module face height (mm) — "sy" in JSON
    module_type: str    # "PbWO4" | "PbGlass" | "LMS" | "Veto" | "Unknown"
    sector:      str    # "Top" | "Bottom" | "Left" | "Right" | "Center"


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_time_calib(
    path: pathlib.Path,
) -> Tuple[Dict[int, float], float, str]:
    """Parse a HyCalTimeCalib JSON file.

    Returns (offsets, default_offset, sha1_hex).
    offsets is keyed by PrimEx id (int).
    """
    raw = path.read_bytes()
    sha1 = hashlib.sha1(raw).hexdigest()
    data = json.loads(raw)

    default_offset = float(data.get("default", 0.0))
    offsets: Dict[int, float] = {}

    for entry in data.get("modules", []):
        offset_ns = float(entry.get("offset_ns", 0.0))
        # Prefer module_id if present
        if "module_id" in entry:
            primex_id = int(entry["module_id"])
            offsets[primex_id] = offset_ns
        elif "name" in entry:
            name = str(entry["name"])
            # Strip leading G or W and map to PrimEx id
            if name.startswith("G") or name.startswith("g"):
                num = int(name[1:])
                primex_id = num  # G: 1–576, no added offset
            elif name.startswith("W") or name.startswith("w"):
                num = int(name[1:])
                primex_id = 1000 + num  # W: 1001–2152
            else:
                # Numeric string — treat directly
                primex_id = int(name)
            offsets[primex_id] = offset_ns

    return offsets, default_offset, sha1


def load_hycal_geometry(path: Optional[pathlib.Path]) -> Dict[int, ModuleGeom]:
    """Load hycal_map.json and return dict[PrimEx_id -> ModuleGeom].

    Returns {} (empty dict) on missing file, unreadable JSON, or non-array
    top-level.  Individual bad entries are skipped with a [WARN] line.
    """
    if path is None or not path.exists():
        if path is not None:
            print(
                f"[WARN] load_hycal_geometry: file not found: {path}",
                file=sys.stderr,
            )
        return {}

    try:
        raw = path.read_bytes()
        data = json.loads(raw)
    except json.JSONDecodeError as exc:
        print(
            f"[WARN] load_hycal_geometry: cannot parse {path}: {exc}",
            file=sys.stderr,
        )
        return {}

    if not isinstance(data, list):
        print(
            f"[WARN] load_hycal_geometry: {path} is not a JSON array",
            file=sys.stderr,
        )
        return {}

    result: Dict[int, ModuleGeom] = {}
    for idx, entry in enumerate(data):
        # Required keys at top level
        for key in ("n", "t", "geo"):
            if key not in entry:
                print(
                    f"[WARN] load_hycal_geometry: skipping entry {idx}: missing key '{key}'",
                    file=sys.stderr,
                )
                break
        else:
            geo = entry.get("geo", {})
            for gkey in ("sx", "sy", "x", "y"):
                if gkey not in geo:
                    print(
                        f"[WARN] load_hycal_geometry: skipping entry {idx}: missing key 'geo.{gkey}'",
                        file=sys.stderr,
                    )
                    break
            else:
                name = str(entry["n"])
                # Derive PrimEx id
                if name.upper().startswith("G"):
                    try:
                        primex_id = int(name[1:])
                    except ValueError:
                        print(
                            f"[WARN] load_hycal_geometry: skipping entry {idx}: cannot parse id from name '{name}'",
                            file=sys.stderr,
                        )
                        continue
                elif name.upper().startswith("W"):
                    try:
                        primex_id = int(name[1:]) + 1000
                    except ValueError:
                        print(
                            f"[WARN] load_hycal_geometry: skipping entry {idx}: cannot parse id from name '{name}'",
                            file=sys.stderr,
                        )
                        continue
                else:
                    print(
                        f"[WARN] load_hycal_geometry: skipping entry {idx}: name '{name}' does not match G/W prefix",
                        file=sys.stderr,
                    )
                    continue

                sector = str(geo.get("sec", "Unknown"))
                result[primex_id] = ModuleGeom(
                    module_id=primex_id,
                    name=name,
                    x_mm=float(geo["x"]),
                    y_mm=float(geo["y"]),
                    size_x_mm=float(geo["sx"]),
                    size_y_mm=float(geo["sy"]),
                    module_type=str(entry.get("t", "Unknown")),
                    sector=sector,
                )

    return result


def apply_offsets_chunk(
    seed_ids: np.ndarray,
    neighbor_ids: np.ndarray,
    dt: np.ndarray,
    offsets: Dict[int, float],
    default: float,
) -> np.ndarray:
    """Vectorised per-chunk offset correction.

    dt_corr = dt_raw - (offset_neighbour - offset_seed)

    If offsets is empty, returns dt unchanged (no copy).
    """
    if not offsets:
        return dt
    off_s = np.array(
        [offsets.get(int(s), default) for s in seed_ids], dtype=np.float32
    )
    off_n = np.array(
        [offsets.get(int(n), default) for n in neighbor_ids], dtype=np.float32
    )
    return dt - (off_n - off_s)


def _get_file_size(path: pathlib.Path) -> int:
    try:
        return path.stat().st_size
    except OSError:
        return 0


def stream_tsv(
    path: pathlib.Path,
    csv_mode: bool,
    no_header: bool,
    ctx: RunContext,
) -> Dict[int, ModuleAccumulator]:
    """Single-pass streaming reader.  Returns dict[neighbor_id -> ModuleAccumulator]."""

    delimiter = "," if csv_mode else "\t"
    accumulators: Dict[int, ModuleAccumulator] = {}

    file_size = _get_file_size(path)
    t0 = time.monotonic()
    rows_read = 0
    rows_kept = 0
    bytes_consumed = 0

    with path.open(newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter=delimiter)

        # --- Header resolution ---
        if not no_header:
            try:
                header_row = next(reader)
            except StopIteration:
                sys.exit("[ERROR] Input file is empty.")

            # Strip whitespace from header names
            header = [h.strip() for h in header_row]
            missing = REQUIRED_COLUMNS - set(header)
            if missing:
                sys.exit(
                    f"[ERROR] Missing required columns in header: {sorted(missing)}"
                )
            col_idx = {name: i for i, name in enumerate(header)}
        else:
            col_idx = {name: i for i, name in enumerate(TSV_COLUMNS)}

        idx_seed_id     = col_idx["seed_id"]
        idx_neighbor_id = col_idx["neighbor_id"]
        idx_dt          = col_idx["dt"]
        idx_seed_energy = col_idx["seed_energy"]
        idx_dist_q      = col_idx["dist_q"]

        # Chunk buffers
        # Python-list buffering; peak memory per chunk is ~a few hundred MB at 1M rows; see plan §7.c
        buf_seed_id:     List[int]   = []
        buf_neighbor_id: List[int]   = []
        buf_dt:          List[float] = []
        buf_seed_energy: List[float] = []
        buf_dist_q:      List[float] = []

        _MODULE_ID_MAX = 100_000  # sanity bound: well above HyCal count (~2200), below garbage

        def _process_chunk() -> None:
            nonlocal rows_kept
            if not buf_seed_id:
                return

            seed_ids     = np.array(buf_seed_id,     dtype=np.int32)
            neighbor_ids = np.array(buf_neighbor_id, dtype=np.int32)
            dt_raw       = np.array(buf_dt,          dtype=np.float32)
            seed_energy  = np.array(buf_seed_energy, dtype=np.float32)
            dist_q       = np.array(buf_dist_q,      dtype=np.float32)

            # Filter mask — committed order: seed_energy → dist_q → module_id
            mask = seed_energy >= ctx.min_seed_energy
            if ctx.max_dist_q < np.inf:
                mask &= dist_q <= ctx.max_dist_q
            if ctx.module_id_filter is not None:
                id_arr = np.zeros(len(neighbor_ids), dtype=bool)
                for mid in ctx.module_id_filter:
                    id_arr |= (neighbor_ids == mid)
                mask &= id_arr

            # n_total_seen updated AFTER the first three masks, BEFORE dt-window
            surviving_ids = neighbor_ids[mask]
            seed_ids_m    = seed_ids[mask]
            dt_raw_m      = dt_raw[mask]

            for uid in np.unique(surviving_ids):
                acc = accumulators.get(int(uid))
                if acc is None:
                    acc = ModuleAccumulator.empty(int(uid), ctx)
                    accumulators[int(uid)] = acc
                count_here = int(np.sum(surviving_ids == uid))
                acc.n_total_seen += count_here

            # dt-window on RAW dt (committed order, §7.b/§7.c)
            dt_mask = (dt_raw_m >= ctx.dt_lo) & (dt_raw_m <= ctx.dt_hi)
            surviving_ids_dt = surviving_ids[dt_mask]
            seed_ids_dt      = seed_ids_m[dt_mask]
            dt_after_window  = dt_raw_m[dt_mask]

            if len(surviving_ids_dt) == 0:
                return

            # Apply offsets (§2.a) — no-op when ctx.offsets is empty
            dt_corr = apply_offsets_chunk(
                seed_ids_dt, surviving_ids_dt,
                dt_after_window,
                ctx.offsets, ctx.default_offset,
            )

            rows_kept += len(dt_corr)

            if ctx.verbose:
                print(
                    f"  [chunk] {len(buf_seed_id)} rows → {len(dt_corr)} kept",
                    file=sys.stderr,
                )

            # Accumulate per module
            unique_ids = np.unique(surviving_ids_dt)
            for uid in unique_ids:
                mid_int = int(uid)
                sel = surviving_ids_dt == uid
                dt_slice = dt_corr[sel]

                acc = accumulators.get(mid_int)
                if acc is None:
                    acc = ModuleAccumulator.empty(mid_int, ctx)
                    accumulators[mid_int] = acc

                counts_dt, _ = np.histogram(dt_slice, bins=ctx.dt_edges)
                acc.dt_hist += counts_dt.astype(np.int64)
                acc.n_entries += int(counts_dt.sum())

                counts_sb, _ = np.histogram(dt_slice, bins=ctx.sideband_edges)
                acc.sideband_hist += counts_sb.astype(np.int64)

            buf_seed_id.clear()
            buf_neighbor_id.clear()
            buf_dt.clear()
            buf_seed_energy.clear()
            buf_dist_q.clear()

        next_progress = ctx.progress_every * 1_000_000 if ctx.progress_every > 0 else 0

        for lineno, row in enumerate(reader, start=1):
            rows_read += 1

            # Approximate bytes consumed (row length + newline)
            bytes_consumed += sum(len(f) for f in row) + len(row)

            if ctx.max_rows > 0 and rows_read > ctx.max_rows:
                break

            try:
                sid = int(row[idx_seed_id])
                nid = int(row[idx_neighbor_id])
                dt_val  = float(row[idx_dt])
                se_val  = float(row[idx_seed_energy])
                dq_val  = float(row[idx_dist_q])
            except (ValueError, IndexError) as exc:
                raw_preview = ("\t".join(row) if row else "")[:80]
                print(
                    f"[ERROR] Cannot parse data row {lineno} in {path}: {exc}\n"
                    f"  Row preview: {raw_preview!r}",
                    file=sys.stderr,
                )
                sys.exit(2)

            # Sanity-check module IDs for garbage values
            if not (0 <= sid < _MODULE_ID_MAX):
                raw_preview = ("\t".join(row) if row else "")[:80]
                print(
                    f"[ERROR] seed_id={sid} out of valid range [0, {_MODULE_ID_MAX}) "
                    f"at data row {lineno} in {path}\n"
                    f"  Row preview: {raw_preview!r}",
                    file=sys.stderr,
                )
                sys.exit(2)
            if not (0 <= nid < _MODULE_ID_MAX):
                raw_preview = ("\t".join(row) if row else "")[:80]
                print(
                    f"[ERROR] neighbor_id={nid} out of valid range [0, {_MODULE_ID_MAX}) "
                    f"at data row {lineno} in {path}\n"
                    f"  Row preview: {raw_preview!r}",
                    file=sys.stderr,
                )
                sys.exit(2)

            buf_seed_id.append(sid)
            buf_neighbor_id.append(nid)
            buf_dt.append(dt_val)
            buf_seed_energy.append(se_val)
            buf_dist_q.append(dq_val)

            if len(buf_seed_id) >= ctx.chunk_rows:
                _process_chunk()

            # Progress reporting
            if next_progress > 0 and rows_read >= next_progress:
                elapsed = time.monotonic() - t0
                frac = bytes_consumed / file_size if file_size > 0 else 0.0
                eta_str = ""
                if frac > 0:
                    eta_s = elapsed / frac * (1.0 - frac)
                    eta_str = f"  ETA {eta_s:.0f}s"
                print(
                    f"[progress] {rows_read/1e6:.1f}M rows read, "
                    f"{rows_kept/1e6:.2f}M kept, "
                    f"{elapsed:.1f}s elapsed{eta_str}",
                    file=sys.stderr,
                    flush=True,
                )
                next_progress += ctx.progress_every * 1_000_000

        # Flush remaining rows
        _process_chunk()

    elapsed_total = time.monotonic() - t0
    print(
        f"[stream_tsv] Done: {rows_read:,} rows read, "
        f"{rows_kept:,} kept, {len(accumulators)} modules, "
        f"{elapsed_total:.1f}s",
        file=sys.stderr,
        flush=True,
    )

    return accumulators


# ---------------------------------------------------------------------------
# Fit model
# ---------------------------------------------------------------------------

def gaussian_plus_const(
    x: np.ndarray,
    A: float,
    mu: float,
    sigma: float,
    B: float,
) -> np.ndarray:
    """Density model: A·exp(-½((x-μ)/σ)²) + B  [counts/ns]."""
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + B


# ---------------------------------------------------------------------------
# fit_module (§9)
# ---------------------------------------------------------------------------

def fit_module(acc: ModuleAccumulator, ctx: RunContext) -> ModuleStats:
    """Fit Gaussian + constant to acc.dt_hist. Returns ModuleStats."""

    base_kwargs = dict(
        module_id=acc.module_id,
        n_entries=acc.n_entries,
        n_total_seen=acc.n_total_seen,
    )

    # 1. Low-stats short-circuit
    if acc.n_entries < ctx.min_entries:
        return ModuleStats(fit_status="low-stats", **base_kwargs)

    # 2. Scipy-missing short-circuit
    if not ctx.scipy_available:
        return ModuleStats(fit_status="scipy-missing", **base_kwargs)

    # 3. Density view
    counts  = acc.dt_hist.astype(np.float64)
    centres = 0.5 * (ctx.dt_edges[:-1] + ctx.dt_edges[1:])
    y       = counts / ctx.dt_bin
    y_err   = np.maximum(np.sqrt(counts), 1.0) / ctx.dt_bin

    # 4. Seed parameters
    mu0 = float(centres[np.argmax(counts)])

    # B0 = median of outer-fifth bins
    n_bins     = len(y)
    outer_cut  = n_bins // 5  # 20% from each end
    tail_bins  = np.concatenate([y[:outer_cut], y[n_bins - outer_cut:]])
    B0 = float(np.median(tail_bins)) if len(tail_bins) > 0 else 0.0
    B0 = max(B0, 0.0)

    A0 = max(0.0, float(y.max()) - B0)
    sigma0 = 1.5  # ns

    # 5. Fit window mask
    mask = (centres >= mu0 - ctx.fit_window) & (centres <= mu0 + ctx.fit_window)
    n_nonzero = int(np.sum(counts[mask] > 0))
    if n_nonzero < 5:
        return ModuleStats(fit_status="low-stats", **base_kwargs)

    # 6. curve_fit call
    bounds = (
        [0.0,       ctx.dt_lo, 0.1,             0.0],
        [np.inf,    ctx.dt_hi, ctx.fit_window,  np.inf],
    )
    try:
        popt, pcov = curve_fit(
            gaussian_plus_const,
            centres[mask], y[mask],
            p0=[A0, mu0, sigma0, B0],
            sigma=y_err[mask],
            absolute_sigma=True,
            bounds=bounds,
            maxfev=5000,
        )
    except Exception:
        return ModuleStats(fit_status="fit-failed", **base_kwargs)

    A_fit, mu_fit, sigma_fit, B_fit = popt

    # Check covariance health: non-finite OR negative variance → singular/ill-conditioned
    # Per plan §9.b: singular covariance maps to fit_status = "fit-failed"
    diag = np.diag(pcov)
    if np.any(~np.isfinite(diag)) or np.any(diag < 0):
        return ModuleStats(
            fit_status="fit-failed",
            mu=float(mu_fit),
            sigma=float(sigma_fit),
            A_per_ns=float(A_fit),
            B_per_ns=float(B_fit),
            mu_err=_NaN,
            sigma_err=_NaN,
            significance=_NaN,
            **base_kwargs,
        )

    mu_err    = float(math.sqrt(diag[1]))
    sigma_err = float(math.sqrt(diag[2]))

    # 7. At-bound detection (§9.a step 7)
    tol_mu    = max(0.5 * ctx.dt_bin, 1e-6)
    tol_sigma = max(0.5 * ctx.dt_bin, 1e-6)

    at_bound = (
        abs(mu_fit - ctx.dt_lo) < tol_mu
        or abs(mu_fit - ctx.dt_hi) < tol_mu
        or abs(sigma_fit - 0.1) < tol_sigma
        or abs(ctx.fit_window - sigma_fit) < tol_sigma
    )

    if at_bound:
        fit_status = "at-bound"
    else:
        fit_status = "ok"

    # 8. Derived quantities (§9.a step 8)
    S    = A_fit * sigma_fit * math.sqrt(2.0 * math.pi) * _ERF3_SQRT2
    B_w  = B_fit * 6.0 * sigma_fit
    sig  = S / math.sqrt(max(B_w, 1.0))

    return ModuleStats(
        module_id=acc.module_id,
        n_entries=acc.n_entries,
        n_total_seen=acc.n_total_seen,
        mu=float(mu_fit),
        mu_err=mu_err,
        sigma=float(sigma_fit),
        sigma_err=sigma_err,
        A_per_ns=float(A_fit),
        B_per_ns=float(B_fit),
        signal_counts=float(S),
        significance=float(sig),
        fit_status=fit_status,
    )


# ---------------------------------------------------------------------------
# Secondary-peak detection (§14.a)
# ---------------------------------------------------------------------------

def detect_secondary_peak(
    acc: ModuleAccumulator,
    stats: ModuleStats,
    ctx: RunContext,
) -> None:
    """Detect a secondary residual peak in acc.dt_hist.  Mutates stats in-place.

    All six secondary_* fields are set (NaN when no candidate or preconditions
    fail).  Appends 'secondary-peak' to stats.flags when significance clears
    ctx.secondary_sig_thr.
    """
    # Preconditions (§14.a)
    if (
        not HAVE_SCIPY
        or stats.fit_status != "ok"
        or abs(stats.mu) > 2.0 * ctx.mu_thr
        or stats.sigma > 1.5 * ctx.sigma_thr
    ):
        return  # leave all secondary_* fields at NaN

    centres = 0.5 * (ctx.dt_edges[:-1] + ctx.dt_edges[1:])

    # Step 1: Build residual density (counts/ns)
    y      = acc.dt_hist.astype(np.float64) / ctx.dt_bin
    y_fit  = gaussian_plus_const(
        centres, stats.A_per_ns, stats.mu, stats.sigma, stats.B_per_ns
    )
    residual = y - y_fit

    # Step 2: Restrict to secondary search window
    mask = (centres >= ctx.secondary_lo) & (centres <= ctx.secondary_hi)
    if int(mask.sum()) < 5:
        return  # too few bins in the search window

    # Step 3: Candidate search with find_peaks
    local_B_counts = max(stats.B_per_ns * ctx.dt_bin, 1.0)
    prominence_counts = max(3.0 * math.sqrt(local_B_counts), 5.0)
    prominence_density = prominence_counts / ctx.dt_bin  # convert to counts/ns

    peaks_in_mask, _ = find_peaks(residual[mask], prominence=prominence_density)
    if len(peaks_in_mask) == 0:
        return

    # Pick the tallest candidate
    tallest_idx = peaks_in_mask[int(np.argmax(residual[mask][peaks_in_mask]))]
    mu2_seed = float(centres[mask][tallest_idx])

    # Step 4: Fit a single bare Gaussian to the residual around the candidate
    A2_seed    = float(residual[mask][tallest_idx])
    sigma2_seed = max(2.0 * ctx.dt_bin, 0.5)

    # Error on residual: approximate Poisson on pre-subtraction counts
    y_err_mask = np.sqrt(acc.dt_hist[mask].astype(np.float64) + local_B_counts) / ctx.dt_bin

    bounds_sec = (
        [0.0,           ctx.secondary_lo, ctx.dt_bin],
        [np.inf,        ctx.secondary_hi, 5.0],
    )
    try:
        popt2, pcov2 = curve_fit(
            lambda x, A2, mu2, sigma2: A2 * np.exp(-0.5 * ((x - mu2) / sigma2) ** 2),
            centres[mask], residual[mask],
            p0=[A2_seed, mu2_seed, sigma2_seed],
            sigma=y_err_mask,
            absolute_sigma=True,
            bounds=bounds_sec,
            maxfev=5000,
        )
    except Exception:
        return  # fit failed — leave secondary_* at NaN

    A2, mu2, sigma2 = popt2

    # Check covariance
    diag2 = np.diag(pcov2)
    if np.any(~np.isfinite(diag2)) or np.any(diag2 < 0):
        return

    mu2_err    = float(math.sqrt(diag2[1]))
    sigma2_err = float(math.sqrt(diag2[2]))

    # Step 5: Derived quantities
    S2   = A2 * sigma2 * math.sqrt(2.0 * math.pi) * _ERF3_SQRT2
    B_w2 = stats.B_per_ns * 6.0 * sigma2
    sec_sig = S2 / math.sqrt(max(B_w2, 1.0))

    # Populate stats
    stats.secondary_mu          = float(mu2)
    stats.secondary_mu_err      = mu2_err
    stats.secondary_sigma       = float(sigma2)
    stats.secondary_sigma_err   = sigma2_err
    stats.secondary_amp_per_ns  = float(A2)
    stats.secondary_significance = float(sec_sig)

    # Flag if above threshold
    if sec_sig >= ctx.secondary_sig_thr:
        stats.flags.append("secondary-peak")


# ---------------------------------------------------------------------------
# Sideband check (§10)
# ---------------------------------------------------------------------------

def sideband_check(
    acc: ModuleAccumulator,
    ctx: RunContext,
    mu0: float,
) -> Tuple[float, float, bool]:
    """Return (sideband_counts, bg_per_ns, overlap_bool).

    overlap_bool=True means the sideband window overlaps the fit window
    centred on mu0; in that case counts and bg_per_ns are NaN.
    """
    sb_lo = ctx.sideband_lo
    sb_hi = ctx.sideband_hi
    fw_lo = mu0 - ctx.fit_window
    fw_hi = mu0 + ctx.fit_window

    # Overlap check: two intervals [a,b] and [c,d] overlap iff a<d and c<b
    overlap = (sb_lo < fw_hi) and (fw_lo < sb_hi)
    if overlap:
        return _NaN, _NaN, True

    counts = float(int(acc.sideband_hist.sum()))
    width  = sb_hi - sb_lo
    bg_per_ns = counts / width if width > 0 else _NaN
    return counts, bg_per_ns, False


# ---------------------------------------------------------------------------
# Flagging (§7.b)
# ---------------------------------------------------------------------------

def flag_module(stats: ModuleStats, ctx: RunContext) -> List[str]:
    """Return list of flag strings in fixed order.

    Fixed order: low-stats → no-peak → offset-peak → wide-peak →
    sideband-overlap → secondary-peak.
    sideband-overlap is appended by sideband_check; secondary-peak is
    appended by detect_secondary_peak.  Both are preserved here.
    """
    flags: List[str] = []

    if stats.n_entries < ctx.min_entries:
        flags.append("low-stats")

    no_peak = stats.fit_status != "ok" or stats.significance < ctx.sig_thr
    if no_peak:
        flags.append("no-peak")

    if stats.fit_status == "ok":
        if abs(stats.mu) > ctx.mu_thr:
            flags.append("offset-peak")
        if stats.sigma > ctx.sigma_thr:
            flags.append("wide-peak")

    # sideband-overlap: set by sideband_check; preserve if already in stats.flags
    if "sideband-overlap" in stats.flags and "sideband-overlap" not in flags:
        flags.append("sideband-overlap")

    # secondary-peak: set by detect_secondary_peak; preserve if already in stats.flags
    if "secondary-peak" in stats.flags and "secondary-peak" not in flags:
        flags.append("secondary-peak")

    return flags


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_heatmap(
    accumulators: Dict[int, ModuleAccumulator],
    ctx: RunContext,
    path_stem: pathlib.Path,
) -> None:
    """Save module_id × ΔT heatmap as PNG and PDF."""
    ids = sorted(accumulators.keys())
    n_mods  = len(ids)
    n_bins  = len(ctx.dt_edges) - 1

    counts2d = np.zeros((n_mods, n_bins), dtype=np.float64)
    for iy, mid in enumerate(ids):
        counts2d[iy, :] = accumulators[mid].dt_hist.astype(np.float64)

    fig, ax = plt.subplots(figsize=(12, max(6, n_mods * 0.05 + 2)))
    im = ax.imshow(
        counts2d,
        aspect="auto",
        norm=LogNorm(vmin=1, vmax=max(counts2d.max(), 2)),
        origin="lower",
        extent=[ctx.dt_lo, ctx.dt_hi, 0, n_mods],
        interpolation="nearest",
    )
    plt.colorbar(im, ax=ax, label="counts")
    ax.axvline(0, color="white", linewidth=0.8, linestyle="--", alpha=0.7)

    # Y-axis ticks: label with actual PrimEx ids, thin if >50
    tick_step = max(1, n_mods // 50) if n_mods > 50 else 1
    tick_positions = list(range(0, n_mods, tick_step))
    tick_labels    = [str(ids[i]) for i in tick_positions]
    ax.set_yticks([p + 0.5 for p in tick_positions])
    ax.set_yticklabels(tick_labels, fontsize=6)

    ax.set_xlabel("ΔT = t_neighbour − t_seed (ns)")
    ax.set_ylabel("Neighbour PrimEx ID")
    tag_str = f"_{ctx.tag}" if ctx.tag else ""
    ax.set_title(f"Per-module ΔT  (run {ctx.run_id}{tag_str})")
    fig.tight_layout()

    for ext in ("png", "pdf"):
        fig.savefig(str(path_stem) + f".{ext}", dpi=150)
    plt.close(fig)


def plot_per_module_pdf(
    accumulators: Dict[int, ModuleAccumulator],
    stats_by_id: Dict[int, ModuleStats],
    ctx: RunContext,
    path: pathlib.Path,
) -> None:
    """Write per-module histogram + fit pages to a multi-page PDF."""
    centres = 0.5 * (ctx.dt_edges[:-1] + ctx.dt_edges[1:])
    x_fit   = np.linspace(ctx.dt_lo, ctx.dt_hi, 500)

    with PdfPages(str(path)) as pdf:
        for mid in sorted(accumulators.keys()):
            acc   = accumulators[mid]
            stats = stats_by_id.get(mid)

            fig, ax = plt.subplots(figsize=(9, 5))
            density = acc.dt_hist.astype(np.float64) / ctx.dt_bin
            ax.step(centres, density, where="mid", color="steelblue",
                    linewidth=0.8, label="data")

            title_extra = ""
            if stats is not None and stats.fit_status == "ok":
                y_fit = gaussian_plus_const(
                    x_fit, stats.A_per_ns, stats.mu, stats.sigma, stats.B_per_ns
                )
                ax.plot(x_fit, y_fit, color="red", linewidth=1.0, label="fit")
                ax.axvline(stats.mu, color="red", linestyle="--", linewidth=0.8)
                for pm in (stats.mu - 3 * stats.sigma,
                           stats.mu + 3 * stats.sigma):
                    ax.axvline(pm, color="orange", linestyle=":", linewidth=0.7)
                title_extra = (
                    f"  μ={stats.mu:.3f}±{stats.mu_err:.3f} ns"
                    f"  σ={stats.sigma:.3f}±{stats.sigma_err:.3f} ns"
                )
                # Secondary fit overlay (if detected)
                if not math.isnan(stats.secondary_mu):
                    # Full composite: primary Gaussian + secondary Gaussian + background
                    y_composite = (
                        gaussian_plus_const(
                            x_fit, stats.A_per_ns, stats.mu, stats.sigma, stats.B_per_ns
                        )
                        + stats.secondary_amp_per_ns * np.exp(
                            -0.5 * ((x_fit - stats.secondary_mu) / stats.secondary_sigma) ** 2
                        )
                    )
                    ax.plot(x_fit, y_composite, color="purple", linewidth=1.0,
                            label="primary + secondary")
                    # Secondary Gaussian alone (over flat background)
                    y_sec = (
                        stats.B_per_ns
                        + stats.secondary_amp_per_ns * np.exp(
                            -0.5 * ((x_fit - stats.secondary_mu) / stats.secondary_sigma) ** 2
                        )
                    )
                    ax.plot(x_fit, y_sec, color="green", linewidth=1.0, linestyle="--",
                            label="secondary fit")
                    # Vertical marker at secondary μ₂
                    ax.axvline(stats.secondary_mu, color="green", linestyle="--", linewidth=0.8)
                    # ±3σ₂ markers
                    for pm in (stats.secondary_mu - 3 * stats.secondary_sigma,
                               stats.secondary_mu + 3 * stats.secondary_sigma):
                        ax.axvline(pm, color="lightgreen", linestyle=":", linewidth=0.7)
                    title_extra += (
                        f"  μ₂={stats.secondary_mu:.2f}±{stats.secondary_mu_err:.2f} ns"
                    )
            elif stats is not None:
                title_extra = f"  NO FIT ({stats.fit_status})"

            ax.set_xlabel("ΔT (ns)")
            ax.set_ylabel("counts / ns")
            ax.set_title(f"Module {mid}{title_extra}")
            ax.axvline(0, color="gray", linewidth=0.5, linestyle="--", alpha=0.5)

            # Text box
            if stats is not None:
                flags_str = "; ".join(stats.flags) if stats.flags else "none"
                textlines = [
                    f"N={stats.n_entries} (seen {stats.n_total_seen})",
                    f"fit_status: {stats.fit_status}",
                ]
                if stats.fit_status == "ok":
                    textlines += [
                        f"B={stats.B_per_ns:.2f} cts/ns",
                        f"S={stats.signal_counts:.0f}  sig={stats.significance:.1f}σ",
                    ]
                    if not math.isnan(stats.secondary_mu):
                        textlines += [
                            f"μ₂={stats.secondary_mu:.2f}±{stats.secondary_mu_err:.2f} ns",
                            f"σ₂={stats.secondary_sigma:.2f} ns  sig₂={stats.secondary_significance:.1f}σ",
                        ]
                textlines.append(f"flags: {flags_str}")
                ax.text(
                    0.97, 0.97, "\n".join(textlines),
                    transform=ax.transAxes,
                    va="top", ha="right",
                    fontsize=7,
                    bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", alpha=0.8),
                )

            ax.legend(fontsize=7)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


# ---------------------------------------------------------------------------
# HyCal spatial map (§14.c)
# ---------------------------------------------------------------------------

def plot_hycal_spatial_map(
    stats_by_id: Dict[int, ModuleStats],
    geometry: Dict[int, "ModuleGeom"],
    ctx: RunContext,
    path_stem: pathlib.Path,
) -> None:
    """Four-panel HyCal spatial map for secondary-peak results.

    Panels:
      (a) top-left  – binary flag map (red = flagged, grey = ok)
      (b) top-right – secondary μ₂ colormap
      (c) bot-left  – secondary significance colormap
      (d) bot-right – x-y scatter of flagged modules with name labels

    Writes both .png and .pdf.  Does nothing if geometry is empty or
    ctx.spatial_map_enabled is False.
    """
    if not ctx.spatial_map_enabled or not geometry:
        return

    import matplotlib.patches as mpatches
    import matplotlib.cm as mplcm
    import matplotlib.colors as mplcolors

    fig, _axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = _axes.flatten()  # index as axes[0..3]

    # Compute HyCal footprint bounds
    all_x = [g.x_mm for g in geometry.values()]
    all_y = [g.y_mm for g in geometry.values()]
    all_sx = [g.size_x_mm for g in geometry.values()]
    all_sy = [g.size_y_mm for g in geometry.values()]
    if all_x:
        x_lo = min(x - sx / 2.0 for x, sx in zip(all_x, all_sx)) - 10
        x_hi = max(x + sx / 2.0 for x, sx in zip(all_x, all_sx)) + 10
        y_lo = min(y - sy / 2.0 for y, sy in zip(all_y, all_sy)) - 10
        y_hi = max(y + sy / 2.0 for y, sy in zip(all_y, all_sy)) + 10
    else:
        x_lo, x_hi, y_lo, y_hi = -700, 700, -700, 700

    # Collect secondary significance values for colormap scaling
    sec_sigs = [
        s.secondary_significance
        for s in stats_by_id.values()
        if not math.isnan(s.secondary_significance)
    ]
    sig_vmax = max(10.0, float(np.nanmax(sec_sigs))) if sec_sigs else None

    # Colormaps
    cmap_viridis = matplotlib.colormaps["viridis"]
    cmap_plasma  = matplotlib.colormaps["plasma"]
    grey_color   = (0.7, 0.7, 0.7, 1.0)
    outline_kw   = dict(edgecolor="black", facecolor="none", linewidth=0.5)
    sec_lo       = ctx.secondary_lo
    sec_hi       = ctx.secondary_hi

    # Fix 1: signature uses `linewidth` (not `lw`) and `facecolor` defaults to "none"
    # so that **outline_kw (which has linewidth + facecolor + edgecolor) binds correctly.
    def _add_rect(ax, geom, facecolor="none", edgecolor="none", linewidth=0.3):
        rect = mpatches.Rectangle(
            (geom.x_mm - geom.size_x_mm / 2.0,
             geom.y_mm - geom.size_y_mm / 2.0),
            geom.size_x_mm, geom.size_y_mm,
            facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth,
        )
        ax.add_patch(rect)

    # Identify flagged modules for panel (d)
    flagged_geoms = [
        geometry[gid]
        for gid, stats in stats_by_id.items()
        if gid in geometry and "secondary-peak" in stats.flags
    ]

    for gid, geom in geometry.items():
        stats = stats_by_id.get(gid)

        # --- Panel (a): Binary flag map ---
        ax_a = axes[0]
        if stats is None:
            _add_rect(ax_a, geom, **outline_kw)
        elif "secondary-peak" in stats.flags:
            _add_rect(ax_a, geom, facecolor="red")
        else:
            _add_rect(ax_a, geom, facecolor=grey_color)

        # --- Panel (b): μ₂ map ---
        ax_b = axes[1]
        if stats is None:
            _add_rect(ax_b, geom, **outline_kw)
        elif not math.isnan(stats.secondary_mu):
            t = (stats.secondary_mu - sec_lo) / (sec_hi - sec_lo)
            t = max(0.0, min(1.0, t))
            _add_rect(ax_b, geom, facecolor=cmap_viridis(t))
        else:
            _add_rect(ax_b, geom, facecolor=grey_color)

        # --- Panel (c): Significance map ---
        ax_c = axes[2]
        if stats is None:
            _add_rect(ax_c, geom, **outline_kw)
        elif sig_vmax is not None and not math.isnan(stats.secondary_significance):
            t = stats.secondary_significance / sig_vmax
            t = max(0.0, min(1.0, t))
            _add_rect(ax_c, geom, facecolor=cmap_plasma(t))
        else:
            _add_rect(ax_c, geom, facecolor=grey_color)

        # --- Panel (d): grey outline context for scatter ---
        _add_rect(axes[3], geom, **outline_kw)

    # --- Panel (d): x-y scatter of flagged modules ---
    ax_d = axes[3]
    n_flagged = len(flagged_geoms)
    _LABEL_LIMIT = 30

    if n_flagged == 0:
        ax_d.text(
            0.5, 0.5, "No secondary peaks detected",
            transform=ax_d.transAxes,
            ha="center", va="center", fontsize=11, color="black",
        )
        ax_d.set_title("Secondary-peak modules (N=0)")
    else:
        show_labels = n_flagged <= _LABEL_LIMIT
        label_note = "" if show_labels else f", labels suppressed"
        ax_d.set_title(f"Secondary-peak modules (N={n_flagged}{label_note})")

        for geom in flagged_geoms:
            stats = stats_by_id.get(
                next(gid for gid, g in geometry.items() if g is geom), None
            )
            sig_val = (
                stats.secondary_significance
                if stats is not None and not math.isnan(stats.secondary_significance)
                else 0.0
            )
            dot_size = min(300.0, 30.0 + 10.0 * sig_val)
            ax_d.scatter(geom.x_mm, geom.y_mm, s=dot_size, color="red",
                         zorder=3, linewidths=0)
            if show_labels:
                ax_d.annotate(
                    geom.name,
                    (geom.x_mm, geom.y_mm),
                    xytext=(3, 3),
                    textcoords="offset points",
                    fontsize=6,
                    color="darkred",
                )

    # Axes decoration (all four panels share the same spatial extent)
    tag_str = f"_{ctx.tag}" if ctx.tag else ""
    for ax in axes:
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")

    axes[0].set_title("secondary-peak flag  (red = flagged)")
    axes[1].set_title("secondary μ₂ (ns)")
    axes[2].set_title("secondary significance  S₂/√(B·6σ₂)")
    # axes[3] title already set above

    # Colorbars for panels (b) and (c)
    sm_b = mplcm.ScalarMappable(
        cmap=cmap_viridis,
        norm=mplcolors.Normalize(vmin=sec_lo, vmax=sec_hi),
    )
    sm_b.set_array([])
    fig.colorbar(sm_b, ax=axes[1], label="secondary μ₂ (ns)", shrink=0.8)

    if sig_vmax is not None:
        sm_c = mplcm.ScalarMappable(
            cmap=cmap_plasma,
            norm=mplcolors.Normalize(vmin=0.0, vmax=sig_vmax),
        )
        sm_c.set_array([])
        fig.colorbar(sm_c, ax=axes[2], label="S₂ / √(B·6σ₂)", shrink=0.8)
    else:
        # All-NaN case — annotate panel (c)
        axes[2].text(
            0.5, 0.5, "No secondary peaks detected",
            transform=axes[2].transAxes,
            ha="center", va="center", fontsize=11, color="black",
        )

    fig.suptitle(
        f"HyCal secondary-peak spatial map  (run {ctx.run_id}{tag_str})",
        fontsize=13,
    )
    fig.tight_layout()

    for ext in ("png", "pdf"):
        fig.savefig(str(path_stem) + f".{ext}", dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# CSV / Markdown output
# ---------------------------------------------------------------------------

def _fmt(v: float, fmt: str = ".4f") -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "nan"
    return format(v, fmt)


CSV_COLUMNS = [
    "module_id", "n_entries", "n_total_seen",
    "mu_ns", "mu_err_ns", "sigma_ns", "sigma_err_ns",
    "A_per_ns", "B_per_ns", "signal_counts", "significance",
    "fit_status",
    "sideband_bg_per_ns", "sideband_counts",
    "secondary_mu", "secondary_sigma", "secondary_amp_per_ns", "secondary_significance",
    "flags",
]


def write_summary_csv(
    stats_list: List[ModuleStats],
    ctx: RunContext,
    path: pathlib.Path,
) -> None:
    """Write summary CSV with # comment header (§8.a column list)."""
    now = datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    calib_mode = "raw (uncorrected)"
    if ctx.time_calib_path is not None:
        calib_mode = f"offset-corrected  json={ctx.time_calib_path}  sha1={ctx.time_calib_sha1}"

    with path.open("w", newline="", encoding="utf-8") as fh:
        # Comment header
        fh.write(f"# module_dt_diagnostic — generated {now}\n")
        fh.write(f"# input: {ctx.input_path}\n")
        fh.write(f"# dt_mode: {calib_mode}\n")
        fh.write(f"# dt_window: [{ctx.dt_lo}, {ctx.dt_hi}] ns  dt_bin: {ctx.dt_bin} ns\n")
        fh.write(f"# fit_window: ±{ctx.fit_window} ns\n")
        fh.write(f"# sideband: [{ctx.sideband_lo}, {ctx.sideband_hi}] ns\n")
        fh.write(
            f"# thresholds: sig>{ctx.sig_thr}  |mu|<{ctx.mu_thr} ns  sigma<{ctx.sigma_thr} ns  min_entries={ctx.min_entries}\n"
        )
        fh.write(
            f"# secondary: search_window=[{ctx.secondary_lo}, {ctx.secondary_hi}] ns"
            f"  sig_threshold={ctx.secondary_sig_thr}\n"
        )
        fh.write(
            f"# filter: min_seed_energy={ctx.min_seed_energy} MeV  max_dist_q={ctx.max_dist_q}\n"
        )
        fh.write(f"# scipy_available: {ctx.scipy_available}\n")
        if ctx.hycal_geometry_path is not None:
            fh.write(f"# hycal_geometry: {ctx.hycal_geometry_path}\n")

        writer = csv.writer(fh)
        writer.writerow(CSV_COLUMNS)

        for s in stats_list:
            flags_str = ";".join(s.flags)
            writer.writerow([
                s.module_id,
                s.n_entries,
                s.n_total_seen,
                _fmt(s.mu),
                _fmt(s.mu_err),
                _fmt(s.sigma),
                _fmt(s.sigma_err),
                _fmt(s.A_per_ns),
                _fmt(s.B_per_ns),
                _fmt(s.signal_counts, ".2f"),
                _fmt(s.significance, ".3f"),
                s.fit_status,
                _fmt(s.sideband_bg_per_ns),
                _fmt(s.sideband_counts, ".0f") if not math.isnan(s.sideband_counts) else "nan",
                _fmt(s.secondary_mu),
                _fmt(s.secondary_sigma),
                _fmt(s.secondary_amp_per_ns),
                _fmt(s.secondary_significance, ".3f"),
                flags_str,
            ])


def _module_name(module_id: int) -> str:
    """Convert a PrimEx id to its G<n>/W<n> name string."""
    if module_id >= 1001:
        return f"W{module_id - 1000}"
    return f"G{module_id}"


def _md_module_line(s: ModuleStats, geometry: Optional[Dict[int, "ModuleGeom"]] = None) -> str:
    """Format one flagged module line per §8.b.

    Returns one or two lines: the primary line, and (if secondary-peak is
    flagged) an indented secondary-fit summary line.
    """
    name = _module_name(s.module_id)
    # Append sector tag if geometry is loaded
    if geometry is not None and s.module_id in geometry:
        geom = geometry[s.module_id]
        name = f"{name} ({geom.sector})"

    mu_str  = _fmt(s.mu, ".3f")  if not math.isnan(s.mu)    else "nan"
    dmu_str = _fmt(s.mu_err, ".3f") if not math.isnan(s.mu_err) else "nan"
    sg_str  = _fmt(s.sigma, ".3f") if not math.isnan(s.sigma) else "nan"
    dsg_str = _fmt(s.sigma_err, ".3f") if not math.isnan(s.sigma_err) else "nan"
    sig_str = _fmt(s.significance, ".1f") if not math.isnan(s.significance) else "nan"
    flags_str = " ".join(f"[{f}]" for f in s.flags)
    primary_line = (
        f"- {name} (id {s.module_id})  n={s.n_entries}/{s.n_total_seen}"
        f"  μ={mu_str} ± {dmu_str} ns"
        f"  σ={sg_str} ± {dsg_str} ns"
        f"  sig={sig_str}  {flags_str}"
    )
    # Secondary-peak follow-up line
    if "secondary-peak" in s.flags:
        s2_mu  = _fmt(s.secondary_mu, ".2f")    if not math.isnan(s.secondary_mu)           else "nan"
        s2_dmu = _fmt(s.secondary_mu_err, ".2f") if not math.isnan(s.secondary_mu_err)       else "nan"
        s2_sg  = _fmt(s.secondary_sigma, ".2f") if not math.isnan(s.secondary_sigma)         else "nan"
        s2_dsg = _fmt(s.secondary_sigma_err, ".2f") if not math.isnan(s.secondary_sigma_err) else "nan"
        s2_sig = _fmt(s.secondary_significance, ".2f") if not math.isnan(s.secondary_significance) else "nan"
        secondary_line = (
            f"    secondary μ={s2_mu} ± {s2_dmu} ns"
            f"  σ={s2_sg} ± {s2_dsg} ns"
            f"  sig={s2_sig}"
        )
        return primary_line + "\n" + secondary_line
    return primary_line


def write_markdown_report(
    stats_list: List[ModuleStats],
    ctx: RunContext,
    path: pathlib.Path,
    geometry: Optional[Dict[int, "ModuleGeom"]] = None,
) -> None:
    """Write Markdown flag report (§8.b, §14.c)."""
    now = datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    calib_mode = "raw (uncorrected)"
    if ctx.time_calib_path is not None:
        calib_mode = (
            f"offset-corrected (json=`{ctx.time_calib_path}`,"
            f" sha1=`{ctx.time_calib_sha1}`)"
        )

    all_flags_vocab = [
        "low-stats", "no-peak", "offset-peak", "wide-peak",
        "sideband-overlap", "secondary-peak",
    ]

    n_ok        = sum(1 for s in stats_list if s.fit_status == "ok" and not s.flags)
    n_flagged   = len(stats_list) - n_ok

    lines: List[str] = []
    lines.append("# module_dt_diagnostic — flag report")
    lines.append("")
    lines.append("## Metadata")
    lines.append("")
    lines.append(f"- **Generated:** {now}")
    lines.append(f"- **Input:** `{ctx.input_path}`")
    lines.append(f"- **Run ID:** `{ctx.run_id}`")
    lines.append(f"- **dt mode:** {calib_mode}")
    lines.append(f"- **ΔT window:** [{ctx.dt_lo}, {ctx.dt_hi}] ns, bin {ctx.dt_bin} ns")
    lines.append(f"- **Fit window:** ±{ctx.fit_window} ns around argmax")
    lines.append(f"- **Sideband:** [{ctx.sideband_lo}, {ctx.sideband_hi}] ns")
    lines.append(
        f"- **Thresholds:** significance > {ctx.sig_thr},"
        f" |μ| < {ctx.mu_thr} ns, σ < {ctx.sigma_thr} ns,"
        f" min_entries = {ctx.min_entries}"
    )
    lines.append(
        f"- **Secondary-peak:** search_window=[{ctx.secondary_lo}, {ctx.secondary_hi}] ns,"
        f" sig_threshold={ctx.secondary_sig_thr}"
    )
    if not ctx.scipy_available:
        lines.append("- **WARNING:** scipy not available — running in plots-only mode."
                     " All fit columns are NaN; fit_status = 'scipy-missing'.")
    lines.append("")
    lines.append("## Summary counts")
    lines.append("")
    lines.append(f"- Total modules analysed: {len(stats_list)}")
    lines.append(f"- Clean (no flags): {n_ok}")
    lines.append(f"- Flagged: {n_flagged}")
    for fl in all_flags_vocab:
        cnt = sum(1 for s in stats_list if fl in s.flags)
        if cnt:
            lines.append(f"  - `{fl}`: {cnt}")
    lines.append("")

    # -----------------------------------------------------------------------
    # §14.c "Modules with secondary peak" roll-up (inserted before per-flag
    # sections)
    # -----------------------------------------------------------------------
    lines.append("## Modules with secondary peak")
    lines.append("")
    lines.append(
        f"Sorted by secondary_significance (descending).  Threshold: "
        f"`--secondary-sig-threshold = {ctx.secondary_sig_thr}`; search window: "
        f"`--secondary-search-window = {ctx.secondary_lo} {ctx.secondary_hi}` ns."
    )
    lines.append("")
    sec_flagged = [s for s in stats_list if "secondary-peak" in s.flags]
    if sec_flagged:
        sec_flagged.sort(
            key=lambda s: -(s.secondary_significance if not math.isnan(s.secondary_significance) else 0.0)
        )
        lines.append("| name  | id   | μ₂ ± err (ns)   | σ₂ ± err (ns)  | sig   |")
        lines.append("|-------|------|-----------------|----------------|-------|")
        for s in sec_flagged:
            name = _module_name(s.module_id)
            mu2  = _fmt(s.secondary_mu, ".2f")      if not math.isnan(s.secondary_mu)           else "nan"
            dmu2 = _fmt(s.secondary_mu_err, ".2f")   if not math.isnan(s.secondary_mu_err)       else "nan"
            sg2  = _fmt(s.secondary_sigma, ".2f")   if not math.isnan(s.secondary_sigma)         else "nan"
            dsg2 = _fmt(s.secondary_sigma_err, ".2f") if not math.isnan(s.secondary_sigma_err)   else "nan"
            ssig = _fmt(s.secondary_significance, ".1f") if not math.isnan(s.secondary_significance) else "nan"
            lines.append(
                f"| {name} | {s.module_id} | {mu2} ± {dmu2} | {sg2} ± {dsg2} | {ssig} |"
            )
    else:
        lines.append("No modules exceeded the secondary-peak significance threshold.")
    lines.append("")

    # Per-flag sections — sort worst-first within each group
    for fl in all_flags_vocab:
        flagged = [s for s in stats_list if fl in s.flags]
        if not flagged:
            continue
        lines.append(f"## Flag: `{fl}`")
        lines.append("")
        # Worst-first sorting per flag type
        if fl == "no-peak":
            flagged.sort(key=lambda s: s.significance if not math.isnan(s.significance) else float("inf"))
        elif fl == "wide-peak":
            flagged.sort(key=lambda s: -(s.sigma if not math.isnan(s.sigma) else 0.0))
        elif fl == "offset-peak":
            flagged.sort(key=lambda s: -(abs(s.mu) if not math.isnan(s.mu) else 0.0))
        elif fl == "secondary-peak":
            flagged.sort(
                key=lambda s: -(s.secondary_significance if not math.isnan(s.secondary_significance) else 0.0)
            )
        else:
            flagged.sort(key=lambda s: s.module_id)

        for s in flagged:
            lines.append(_md_module_line(s, geometry))
        lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("## Interpretation note")
    lines.append("")
    if ctx.time_calib_path is None:
        lines.append(
            "μ reported here is the **raw (uncorrected)** peak position.  "
            "For a module with id `m`, μ ≈ offset_m − mean(offset_seed), "
            "where the seed-module mix varies by run conditions.  "
            "This is not a post-calibration residual."
        )
    else:
        lines.append(
            "μ reported here is the **post-calibration residual** "
            f"after subtracting offsets from `{ctx.time_calib_path}`.  "
            "A well-calibrated module should show μ ≈ 0 ns."
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

_RUN_PAT = re.compile(r"prad[_\-]?(\d+)", re.IGNORECASE)


def sniff_run_id(input_path: pathlib.Path) -> str:
    """Infer run ID from filename.  Returns 'run_<N>' or a timestamp."""
    m = _RUN_PAT.search(input_path.name)
    if m:
        return f"run_{m.group(1)}"
    ts = datetime.datetime.now(datetime.UTC).strftime("%Y%m%dT%H%M%S")
    return ts


def _validate_args(ns: argparse.Namespace) -> None:
    """Validate CLI arg combinations.  Raises SystemExit(2) on problems."""

    def _die(msg: str) -> None:
        print(f"[ERROR] {msg}", file=sys.stderr)
        raise SystemExit(2)

    dt_lo, dt_hi = ns.dt_window
    if dt_lo >= dt_hi:
        _die("--dt-window: lo must be < hi")

    sb_lo, sb_hi = ns.sideband
    if sb_lo >= sb_hi:
        _die("--sideband: lo must be < hi")
    if sb_lo < dt_lo or sb_hi > dt_hi:
        _die(
            f"--sideband [{sb_lo}, {sb_hi}] must be inside "
            f"--dt-window [{dt_lo}, {dt_hi}]"
        )

    if ns.dt_bin >= (dt_hi - dt_lo) / 10.0:
        _die("--dt-bin is too large (must be < (dt_hi - dt_lo) / 10)")

    for name in ("sig_threshold", "mu_threshold", "sigma_threshold"):
        val = getattr(ns, name)
        if val < 0:
            _die(f"--{name.replace('_','-')}: threshold must be non-negative")

    if ns.time_calib is not None and not pathlib.Path(ns.time_calib).exists():
        _die(f"--time-calib: file not found: {ns.time_calib}")

    if ns.chunk_rows < 1:
        _die("--chunk-rows must be >= 1")
    if ns.progress_every < 0:
        _die("--progress-every must be >= 0")
    if ns.max_rows < 0:
        _die("--max-rows must be >= 0")

    # §14.e secondary-peak validation
    sec_lo, sec_hi = ns.secondary_search_window
    if sec_lo >= sec_hi:
        _die("--secondary-search-window: lo must be < hi")
    if sec_lo < dt_lo or sec_hi > dt_hi:
        _die(
            f"--secondary-search-window [{sec_lo}, {sec_hi}] must be inside "
            f"--dt-window [{dt_lo}, {dt_hi}]"
        )
    if ns.secondary_sig_threshold < 0:
        _die("--secondary-sig-threshold must be non-negative")


# ---------------------------------------------------------------------------
# Self-test (§11.a)
# ---------------------------------------------------------------------------

def _self_test() -> int:
    """Run both synthetic validation cases.  Returns 0 if both pass."""
    import tempfile

    rng = np.random.default_rng(42)
    passed = 0
    failed = 0

    # -----------------------------------------------------------------------
    # Case 1 — Synthetic Gaussian + flat background
    # -----------------------------------------------------------------------
    print("[self-test] Case 1: Synthetic Gaussian + flat background ...")

    # Build a minimal RunContext with §13.a defaults
    ctx = RunContext(
        input_path=pathlib.Path("synthetic"),
        output_dir=pathlib.Path("."),
        run_id="test",
        tag="",
        dt_lo=-40.0,
        dt_hi=40.0,
        dt_bin=0.5,
        fit_window=8.0,
        sideband_lo=-40.0,
        sideband_hi=-20.0,
        sig_thr=5.0,
        mu_thr=2.0,
        sigma_thr=3.0,
        min_entries=200,
        module_id_filter=None,
        min_seed_energy=0.0,
        max_dist_q=float("inf"),
        time_calib_path=None,
        time_calib_sha1="",
        offsets={},
        default_offset=0.0,
        scipy_available=HAVE_SCIPY,
        chunk_rows=1_000_000,
        progress_every=0,
        max_rows=0,
        verbose=False,
        no_pdf=True,
    )

    N_sig = 5000
    N_bg  = 2000
    true_mu    = 0.3
    true_sigma = 1.1
    true_B_per_ns = N_bg / (ctx.dt_hi - ctx.dt_lo)  # = 25 counts/ns

    sig_samples = rng.normal(true_mu, true_sigma, N_sig)
    bg_samples  = rng.uniform(ctx.dt_lo, ctx.dt_hi, N_bg)
    all_samples = np.concatenate([sig_samples, bg_samples])

    acc = ModuleAccumulator.empty(9999, ctx)
    dt_hist_full, _ = np.histogram(all_samples, bins=ctx.dt_edges)
    acc.dt_hist = dt_hist_full.astype(np.int64)

    sb_mask = (all_samples >= ctx.sideband_lo) & (all_samples <= ctx.sideband_hi)
    sb_hist, _ = np.histogram(all_samples[sb_mask], bins=ctx.sideband_edges)
    acc.sideband_hist = sb_hist.astype(np.int64)

    acc.n_entries    = int(acc.dt_hist.sum())
    acc.n_total_seen = N_sig + N_bg

    ok1 = True

    if not HAVE_SCIPY:
        print(
            "  [SKIP] scipy not available — cannot run fit.  "
            "Install scipy to enable this case.",
            file=sys.stderr,
        )
        ok1 = False
    else:
        stats = fit_module(acc, ctx)
        checks = {
            "fit_status == ok":        stats.fit_status == "ok",
            "|mu - 0.3| < 3*mu_err":  abs(stats.mu - true_mu) < 3 * stats.mu_err,
            "|mu - 0.3| < 0.1 ns":    abs(stats.mu - true_mu) < 0.1,
            "|sigma-1.1| < 3*sig_err": abs(stats.sigma - true_sigma) < 3 * stats.sigma_err,
            "|sigma-1.1| < 0.1 ns":   abs(stats.sigma - true_sigma) < 0.1,
            "|B - 25|/25 < 0.20":     abs(stats.B_per_ns - true_B_per_ns) / true_B_per_ns < 0.20,
            "significance > 10":       stats.significance > 10.0,
            "n_entries == hist.sum()": acc.n_entries == int(acc.dt_hist.sum()),
        }
        for desc, result in checks.items():
            status = "PASS" if result else "FAIL"
            print(f"  {status}  {desc}")
            if not result:
                ok1 = False
                # Print actual values on failure
                if "B" in desc:
                    print(f"       B_per_ns={stats.B_per_ns:.3f}, truth={true_B_per_ns:.3f}")
                elif "mu" in desc:
                    print(f"       mu={stats.mu:.4f}, mu_err={stats.mu_err:.4f}")
                elif "sigma" in desc:
                    print(f"       sigma={stats.sigma:.4f}, sigma_err={stats.sigma_err:.4f}")

    if ok1:
        print("[self-test] Case 1: PASS")
        passed += 1
    else:
        print("[self-test] Case 1: FAIL")
        failed += 1

    # -----------------------------------------------------------------------
    # Case 2 — Minimal time-calib JSON parse + apply_offsets_chunk sign check
    # -----------------------------------------------------------------------
    print("[self-test] Case 2: Time-calib JSON parse + offset sign ...")

    calib_json = {
        "default": 0.5,
        "modules": [
            {"name": "W4", "offset_ns": 1.0},
            {"module_id": 42, "offset_ns": 3.0},
        ],
    }
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", delete=False, encoding="utf-8"
    ) as tf:
        json.dump(calib_json, tf)
        tf_path = pathlib.Path(tf.name)

    ok2 = True
    try:
        offsets, default, sha1 = load_time_calib(tf_path)

        # W4 → module_id = 1004, offset = 1.0
        # module_id 42 → offset = 3.0
        checks2 = {
            "W4 → id 1004 in offsets":  1004 in offsets,
            "offsets[1004] == 1.0":     abs(offsets.get(1004, _NaN) - 1.0) < 1e-9,
            "id 42 in offsets":         42 in offsets,
            "offsets[42] == 3.0":       abs(offsets.get(42, _NaN) - 3.0) < 1e-9,
            "default == 0.5":           abs(default - 0.5) < 1e-9,
            "sha1 is non-empty str":    isinstance(sha1, str) and len(sha1) > 0,
        }
        for desc, result in checks2.items():
            status = "PASS" if result else "FAIL"
            print(f"  {status}  {desc}")
            if not result:
                ok2 = False

        # apply_offsets_chunk sign check:
        #   row 0: seed=1004 (off=1.0), nbr=42  (off=3.0): corr = 10 - (3-1) = 8
        #   row 1: seed=42   (off=3.0), nbr=1004 (off=1.0): corr = 10 - (1-3) = 12
        #   sum = 20.0
        seed_ids     = np.array([1004, 42],   dtype=np.int32)
        neighbor_ids = np.array([42,   1004], dtype=np.int32)
        dt_in        = np.array([10.0, 10.0], dtype=np.float32)
        dt_corr = apply_offsets_chunk(seed_ids, neighbor_ids, dt_in, offsets, default)

        sign_ok = abs(float(dt_corr[0]) - 8.0) < 1e-4 and abs(float(dt_corr[1]) - 12.0) < 1e-4
        sum_ok  = abs(float(dt_corr.sum()) - 20.0) < 1e-4

        print(f"  {'PASS' if sign_ok else 'FAIL'}  dt_corr values: {dt_corr[0]:.4f}, {dt_corr[1]:.4f} (expect 8.0, 12.0)")
        print(f"  {'PASS' if sum_ok  else 'FAIL'}  dt_corr sum == 20.0 ({float(dt_corr.sum()):.4f})")
        if not sign_ok or not sum_ok:
            ok2 = False
    finally:
        tf_path.unlink(missing_ok=True)

    if ok2:
        print("[self-test] Case 2: PASS")
        passed += 1
    else:
        print("[self-test] Case 2: FAIL")
        failed += 1

    # -----------------------------------------------------------------------
    # Case 3 — Secondary-peak detection (§14.h)
    # -----------------------------------------------------------------------
    print("[self-test] Case 3: Secondary-peak detection ...")

    ok3 = True
    if not HAVE_SCIPY:
        print(
            "  [SKIP] scipy not available — cannot run secondary-peak detection.  "
            "Install scipy to enable this case.",
            file=sys.stderr,
        )
        ok3 = False
    else:
        # Build a fresh RunContext with §13.a defaults (same as ctx above)
        ctx3 = RunContext(
            input_path=pathlib.Path("synthetic"),
            output_dir=pathlib.Path("."),
            run_id="test",
            tag="",
            dt_lo=-40.0,
            dt_hi=40.0,
            dt_bin=0.5,
            fit_window=8.0,
            sideband_lo=-40.0,
            sideband_hi=-20.0,
            sig_thr=5.0,
            mu_thr=2.0,
            sigma_thr=3.0,
            min_entries=200,
            module_id_filter=None,
            min_seed_energy=0.0,
            max_dist_q=float("inf"),
            time_calib_path=None,
            time_calib_sha1="",
            offsets={},
            default_offset=0.0,
            scipy_available=HAVE_SCIPY,
            chunk_rows=1_000_000,
            progress_every=0,
            max_rows=0,
            verbose=False,
            no_pdf=True,
            secondary_lo=5.0,
            secondary_hi=25.0,
            secondary_sig_thr=3.0,
        )

        rng3 = np.random.default_rng(7)

        # Primary: A1=1000 counts/ns, μ1=0, σ1=1 ns
        # N_sig1 ≈ A1 * σ1 * √(2π) * (dt_hi - dt_lo) truncated;
        # to get ~5000 signal counts in the histogram use N_sig1=5000 samples
        N_sig1 = 5000
        N_sig2 = 1000   # secondary peak at μ2=15, σ2=1
        N_bg   = 2000   # flat background

        sig1_samples = rng3.normal(0.0, 1.0, N_sig1)
        sig2_samples = rng3.normal(15.0, 1.0, N_sig2)
        bg_samples3  = rng3.uniform(ctx3.dt_lo, ctx3.dt_hi, N_bg)
        all3 = np.concatenate([sig1_samples, sig2_samples, bg_samples3])

        acc3 = ModuleAccumulator.empty(8888, ctx3)
        dt_hist3, _ = np.histogram(all3, bins=ctx3.dt_edges)
        acc3.dt_hist     = dt_hist3.astype(np.int64)
        acc3.n_entries   = int(acc3.dt_hist.sum())
        acc3.n_total_seen = N_sig1 + N_sig2 + N_bg

        # Sideband histogram (needed for flag_module)
        sb3_mask = (all3 >= ctx3.sideband_lo) & (all3 <= ctx3.sideband_hi)
        sb3_hist, _ = np.histogram(all3[sb3_mask], bins=ctx3.sideband_edges)
        acc3.sideband_hist = sb3_hist.astype(np.int64)

        stats3 = fit_module(acc3, ctx3)

        # Sideband check before flag
        mu0_3 = stats3.mu if stats3.fit_status == "ok" else float(
            (0.5 * (ctx3.dt_edges[:-1] + ctx3.dt_edges[1:]))[np.argmax(acc3.dt_hist)]
        )
        sb3_cnt, sb3_bg, sb3_ov = sideband_check(acc3, ctx3, mu0_3)
        stats3.sideband_counts    = sb3_cnt
        stats3.sideband_bg_per_ns = sb3_bg
        if sb3_ov:
            stats3.flags.append("sideband-overlap")
        stats3.flags = flag_module(stats3, ctx3)

        # Run secondary detection
        detect_secondary_peak(acc3, stats3, ctx3)

        # Re-run flag_module to pick up secondary-peak flag
        stats3.flags = flag_module(stats3, ctx3)

        checks3 = {
            "primary fit_status == ok":       stats3.fit_status == "ok",
            "|mu − 0| < 0.1 ns":             abs(stats3.mu) < 0.1 if not math.isnan(stats3.mu) else False,
            "|sigma − 1| < 0.1 ns":          abs(stats3.sigma - 1.0) < 0.1 if not math.isnan(stats3.sigma) else False,
            "|secondary_mu − 15| < 0.5 ns":  abs(stats3.secondary_mu - 15.0) < 0.5 if not math.isnan(stats3.secondary_mu) else False,
            "secondary_significance > 10":    stats3.secondary_significance > 10.0 if not math.isnan(stats3.secondary_significance) else False,
            "'secondary-peak' in flags":      "secondary-peak" in stats3.flags,
        }
        for desc, result in checks3.items():
            status = "PASS" if result else "FAIL"
            print(f"  {status}  {desc}")
            if not result:
                ok3 = False
                if "secondary_mu" in desc:
                    print(f"       secondary_mu={stats3.secondary_mu:.4f}")
                elif "secondary_sig" in desc:
                    print(f"       secondary_significance={stats3.secondary_significance:.4f}")
                elif "mu" in desc:
                    print(f"       mu={stats3.mu:.4f}")
                elif "sigma" in desc:
                    print(f"       sigma={stats3.sigma:.4f}")
                elif "flags" in desc:
                    print(f"       flags={stats3.flags}")

    if ok3:
        print("[self-test] Case 3: PASS")
        passed += 1
    else:
        print("[self-test] Case 3: FAIL")
        failed += 1

    print(f"\n[self-test] Results: {passed} passed, {failed} failed.")
    return 0 if failed == 0 else 1


# ---------------------------------------------------------------------------
# CLI / main
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="module_dt_diagnostic.py",
        description=(
            "Per-module ΔT diagnostic for HyCal timing.\n"
            "Consumes the TSV/CSV from study_hycal_timing.py and produces\n"
            "per-module Gaussian fits, a 2-D heatmap, and a flag report."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    ap.add_argument("input_tsv", metavar="INPUT_TSV", nargs="?",
                    help="TSV/CSV from study_hycal_timing.py")
    ap.add_argument("output_dir", metavar="OUTPUT_DIR", nargs="?",
                    help="Directory to write plots + summary into (created if missing)")

    # Input shape
    grp_input = ap.add_argument_group("input shape")
    grp_input.add_argument("--csv", action="store_true",
                           help="Parse INPUT_TSV as CSV (default: TSV)")
    grp_input.add_argument("--no-header", action="store_true",
                           help="TSV/CSV has no header row; use fixed column order from §2")

    # Grouping / filtering
    grp_filt = ap.add_argument_group("grouping / filtering")
    grp_filt.add_argument("--module-id", type=int, nargs="+", metavar="ID",
                          default=None,
                          help="Restrict to these neighbour PrimEx IDs (default: all)")
    grp_filt.add_argument("--min-seed-energy", type=float, default=0.0, metavar="MEV",
                          help="Drop rows with seed_energy < MEV (default 0)")
    grp_filt.add_argument("--max-dist-q", type=float, default=float("inf"), metavar="F",
                          help="Drop rows with dist_q > F (default inf)")
    grp_filt.add_argument("--min-entries", type=int, default=200, metavar="N",
                          help="Modules with fewer than N rows are flagged low-stats (default 200)")

    # ΔT windowing
    grp_dt = ap.add_argument_group("ΔT windowing")
    grp_dt.add_argument("--dt-window", type=float, nargs=2, default=[-40.0, 40.0],
                        metavar=("LO", "HI"),
                        help="Histogram range in ns (default -40 40)")
    grp_dt.add_argument("--dt-bin", type=float, default=0.5, metavar="F",
                        help="Bin width in ns (default 0.5)")
    grp_dt.add_argument("--fit-window", type=float, default=8.0, metavar="F",
                        help="Half-width (ns) around argmax for Gaussian seed window (default 8.0)")
    grp_dt.add_argument("--sideband", type=float, nargs=2, default=[-40.0, -20.0],
                        metavar=("LO", "HI"),
                        help="Off-time sideband range in ns (default -40 -20)")

    # Time calibration
    grp_calib = ap.add_argument_group("time calibration (optional)")
    grp_calib.add_argument("--time-calib", metavar="PATH", default=None,
                           help="Apply per-module offsets from this JSON (HyCalTimeCalib schema)")

    # Flagging thresholds
    grp_thr = ap.add_argument_group("flagging thresholds")
    grp_thr.add_argument("--sig-threshold", type=float, default=5.0, metavar="F",
                         help="Flag 'no-peak' if significance < F (default 5.0)")
    grp_thr.add_argument("--mu-threshold", type=float, default=2.0, metavar="F",
                         help="Flag 'offset-peak' if |μ| > F ns (default 2.0)")
    grp_thr.add_argument("--sigma-threshold", type=float, default=3.0, metavar="F",
                         help="Flag 'wide-peak' if σ > F ns (default 3.0)")

    # Streaming / performance
    grp_stream = ap.add_argument_group("streaming / performance")
    grp_stream.add_argument("--chunk-rows", type=int, default=1_000_000, metavar="N",
                            help="Rows per read chunk (default 1000000; performance knob only)")
    grp_stream.add_argument("--progress-every", type=int, default=10, metavar="N",
                            help="Print progress every N million rows (default 10; 0=silent)")
    grp_stream.add_argument("--max-rows", type=int, default=0, metavar="N",
                            help="Hard cap on rows read (default 0 = no cap)")

    # Secondary-peak detection (§14.e)
    _script_dir = pathlib.Path(__file__).resolve().parent
    _env_db = os.environ.get("PRAD2_DATABASE_DIR", "")
    if _env_db:
        _default_geom = str(pathlib.Path(_env_db).resolve() / "hycal_map.json")
    else:
        _default_geom = str((_script_dir / ".." / ".." / "database" / "hycal_map.json").resolve())

    grp_sec = ap.add_argument_group("secondary-peak detection (§14)")
    grp_sec.add_argument(
        "--secondary-search-window", type=float, nargs=2,
        default=[5.0, 25.0], metavar=("LO", "HI"),
        help="ΔT region (ns) to search for residual secondary peak (default 5.0 25.0)",
    )
    grp_sec.add_argument(
        "--secondary-sig-threshold", type=float, default=3.0, metavar="T",
        help="Minimum S₂/√(B·6σ₂) for the secondary-peak flag (default 3.0)",
    )
    grp_sec.add_argument(
        "--hycal-geometry", metavar="PATH", default=_default_geom,
        help=(
            "Path to hycal_map.json for the spatial map "
            "(default: $PRAD2_DATABASE_DIR/hycal_map.json or <repo>/database/hycal_map.json; "
            "missing file → warn + skip spatial map)"
        ),
    )
    grp_sec.add_argument(
        "--no-spatial-map", action="store_true",
        help="Skip the HyCal spatial-map output.  Secondary detection still runs.",
    )

    # Output control
    grp_out = ap.add_argument_group("output control")
    grp_out.add_argument("--run-id", metavar="STR", default=None,
                         help="Appended to output filenames (default: sniffed from input filename)")
    grp_out.add_argument("--tag", metavar="STR", default="",
                         help="Extra tag appended after run-id (default '')")
    grp_out.add_argument("--no-pdf", action="store_true",
                         help="Skip the per-module PDF (useful for CI runs)")
    grp_out.add_argument("-v", "--verbose", action="store_true",
                         help="Log per-module fit outcomes to stderr")

    # Self-test
    ap.add_argument("--self-test", action="store_true",
                    help="Run built-in validation cases (§11.a) and exit")

    return ap


def main(argv: Optional[List[str]] = None) -> int:
    ap = _build_parser()
    ns = ap.parse_args(argv)

    # Emit scipy-missing warning at start of main() — not at import time (§5)
    if not HAVE_SCIPY:
        print(
            "[WARN] scipy not available — running in plots-only mode",
            file=sys.stderr,
        )

    # --self-test short-circuits before any TSV I/O
    if ns.self_test:
        return _self_test()

    # --- Positional args check ---
    if not ns.input_tsv or not ns.output_dir:
        ap.error("INPUT_TSV and OUTPUT_DIR are required (or use --self-test)")

    _validate_args(ns)

    input_path = pathlib.Path(ns.input_tsv)
    output_dir = pathlib.Path(ns.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    run_id = ns.run_id if ns.run_id else sniff_run_id(input_path)
    tag    = ns.tag

    # Load time calibration if requested
    offsets: Dict[int, float] = {}
    default_offset = 0.0
    time_calib_sha1 = ""
    time_calib_path: Optional[pathlib.Path] = None
    if ns.time_calib is not None:
        time_calib_path = pathlib.Path(ns.time_calib)
        offsets, default_offset, time_calib_sha1 = load_time_calib(time_calib_path)
        if ns.verbose:
            print(
                f"[calib] Loaded {len(offsets)} module offsets from {time_calib_path}"
                f" (sha1={time_calib_sha1[:8]}…)",
                file=sys.stderr,
            )

    # Resolve geometry path (§14.b, §13.a)
    hycal_geom_path: Optional[pathlib.Path] = None
    if ns.hycal_geometry:
        hycal_geom_path = pathlib.Path(ns.hycal_geometry).resolve()

    # Build RunContext
    ctx = RunContext(
        input_path=input_path,
        output_dir=output_dir,
        run_id=run_id,
        tag=tag,
        dt_lo=ns.dt_window[0],
        dt_hi=ns.dt_window[1],
        dt_bin=ns.dt_bin,
        fit_window=ns.fit_window,
        sideband_lo=ns.sideband[0],
        sideband_hi=ns.sideband[1],
        sig_thr=ns.sig_threshold,
        mu_thr=ns.mu_threshold,
        sigma_thr=ns.sigma_threshold,
        min_entries=ns.min_entries,
        module_id_filter=set(ns.module_id) if ns.module_id else None,
        min_seed_energy=ns.min_seed_energy,
        max_dist_q=ns.max_dist_q,
        time_calib_path=time_calib_path,
        time_calib_sha1=time_calib_sha1,
        offsets=offsets,
        default_offset=default_offset,
        scipy_available=HAVE_SCIPY,
        chunk_rows=ns.chunk_rows,
        progress_every=ns.progress_every,
        max_rows=ns.max_rows,
        verbose=ns.verbose,
        no_pdf=ns.no_pdf,
        secondary_lo=ns.secondary_search_window[0],
        secondary_hi=ns.secondary_search_window[1],
        secondary_sig_thr=ns.secondary_sig_threshold,
        hycal_geometry_path=hycal_geom_path,
        spatial_map_enabled=not ns.no_spatial_map,
    )

    # Single streaming pass over the TSV
    accumulators = stream_tsv(input_path, ns.csv, ns.no_header, ctx)

    if not accumulators:
        print("[ERROR] Zero rows survived filtering.", file=sys.stderr)
        return 1

    # Per-accumulator: fit → sideband → flag → secondary detection (§14.g step 5a)
    stats_list: List[ModuleStats] = []
    stats_by_id: Dict[int, ModuleStats] = {}

    for mid in sorted(accumulators.keys()):
        acc   = accumulators[mid]
        stats = fit_module(acc, ctx)

        # Determine mu0 for sideband overlap check
        if stats.fit_status == "ok":
            mu0 = stats.mu
        else:
            centres = 0.5 * (ctx.dt_edges[:-1] + ctx.dt_edges[1:])
            mu0 = float(centres[np.argmax(acc.dt_hist)])

        sb_counts, sb_bg, sb_overlap = sideband_check(acc, ctx, mu0)
        stats.sideband_counts    = sb_counts
        stats.sideband_bg_per_ns = sb_bg
        if sb_overlap:
            stats.flags.append("sideband-overlap")

        # Primary flagging
        stats.flags = flag_module(stats, ctx)

        # Step 5a: secondary-peak detection (§14.g)
        detect_secondary_peak(acc, stats, ctx)
        # Re-run flag_module to pick up secondary-peak flag in fixed order
        stats.flags = flag_module(stats, ctx)

        if ns.verbose:
            print(
                f"  module {mid:5d}  n={stats.n_entries:8,}  "
                f"fit={stats.fit_status:15s}  flags={stats.flags}",
                file=sys.stderr,
            )

        stats_list.append(stats)
        stats_by_id[mid] = stats

    # Step 5b: load HyCal geometry (§14.g)
    geometry: Dict[int, ModuleGeom] = {}
    if ctx.spatial_map_enabled:
        geometry = load_hycal_geometry(ctx.hycal_geometry_path)
        if not geometry:
            # Empty dict means missing or malformed file; spatial map skipped
            ctx.spatial_map_enabled = False

    # --- Output filenames ---
    tag_suffix = f"_{tag}" if tag else ""
    stem = f"{run_id}{tag_suffix}"

    csv_path         = output_dir / f"module_dt_summary_{stem}.csv"
    heatmap_stem     = output_dir / f"module_dt_heatmap_{stem}"
    fits_pdf         = output_dir / f"module_dt_fits_{stem}.pdf"
    md_path          = output_dir / f"module_dt_flags_{stem}.md"
    spatial_map_stem = output_dir / f"hycal_secondary_peak_map_{stem}"

    write_summary_csv(stats_list, ctx, csv_path)
    print(f"[out] {csv_path}", file=sys.stderr)

    plot_heatmap(accumulators, ctx, heatmap_stem)
    print(f"[out] {heatmap_stem}.png / .pdf", file=sys.stderr)

    if not ns.no_pdf:
        plot_per_module_pdf(accumulators, stats_by_id, ctx, fits_pdf)
        print(f"[out] {fits_pdf}", file=sys.stderr)

    write_markdown_report(stats_list, ctx, md_path, geometry=geometry if geometry else None)
    print(f"[out] {md_path}", file=sys.stderr)

    # Spatial map (§14.g step 6, §14.c)
    if ctx.spatial_map_enabled and geometry:
        plot_hycal_spatial_map(stats_by_id, geometry, ctx, spatial_map_stem)
        print(f"[out] {spatial_map_stem}.png / .pdf", file=sys.stderr)

    # Summary to stdout
    n_flagged   = sum(1 for s in stats_list if s.flags)
    n_no_peak   = sum(1 for s in stats_list if "no-peak" in s.flags)
    n_offset    = sum(1 for s in stats_list if "offset-peak" in s.flags)
    n_wide      = sum(1 for s in stats_list if "wide-peak" in s.flags)
    n_ls        = sum(1 for s in stats_list if "low-stats" in s.flags)
    n_secondary = sum(1 for s in stats_list if "secondary-peak" in s.flags)
    print(
        f"{len(stats_list)} modules analysed, {n_flagged} flagged "
        f"(no-peak={n_no_peak}, offset={n_offset}, wide={n_wide}, "
        f"low-stats={n_ls}, secondary-peak={n_secondary})"
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
