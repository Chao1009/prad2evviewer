#!/usr/bin/env python3
"""
_common.py — shared helpers for analysis/pyscripts/.

The analysis scripts do the same setup (load DAQ config, runinfo,
HyCal+GEM systems, discover EVIO splits), run the same EVIO physics-event
loop and share the per-event boilerplate (waveform → cluster, GEM
ProcessEvent + Reconstruct, lab-frame transform, HyCal↔GEM matching).
This module factors that out so each script only differs in its
per-event accumulation + output.

The heavy wiring lives in `prad2::PipelineBuilder` (prad2det) and is
reached here through `det.PipelineBuilder` — so the same C++ code that
the analysis scripts and the live server use also drives the Python
analyses.

Requires:
  prad2py    (built from python/, exposes dec.* + det.*)
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterator, Optional, Sequence, Tuple

try:
    from prad2py import dec, det
except ImportError as e:
    raise SystemExit(
        f"[ERROR] cannot import prad2py: {e}\n"
        "        Build the python bindings (cmake -DBUILD_PYTHON=ON) and "
        "ensure the install directory is on PYTHONPATH."
    )


# ---- Path / run-number helpers ----

_RUN_PAT = re.compile(r"(?:prad|run)_0*(\d+)", re.IGNORECASE)


def extract_run_number(path: str) -> int:
    """Sniff the run number out of 'prad_NNNNNN.evio.*'-style names. -1 if none."""
    if not path:
        return -1
    m = _RUN_PAT.search(path)
    if not m:
        return -1
    try:
        return int(m.group(1))
    except ValueError:
        return -1


def resolve_db_path(p: str) -> str:
    """Resolve a possibly-relative database path against PRAD2_DATABASE_DIR."""
    if not p:
        return p
    if os.path.isabs(p):
        return p
    db = os.environ.get("PRAD2_DATABASE_DIR")
    if db is None:
        return p
    return os.path.join(db, p)


def hycal_pos_resolution(A: float, B: float, C: float, energy_mev: float) -> float:
    """sigma(E) at the HyCal face (mm), mirroring HyCalSystem::PositionResolution."""
    E_GeV = energy_mev / 1000.0 if energy_mev > 0 else 1e-6
    a = A / math.sqrt(E_GeV)
    b = B / E_GeV
    return math.sqrt(a * a + b * b + C * C)


def discover_split_files(any_path: str) -> list[str]:
    """Three modes by input shape:
      * '*' in path  → literal shell-style glob: expanded by Python's glob
        module so users can pick a subset (e.g. 'prad_024236.evio.0000*'
        gets splits .00000–.00009).  Quote it on the shell to keep the
        shell from expanding it first.
      * directory    → enumerate every prad_<run>.evio.<digits> in the dir,
        sniff run from dir name, warn (stderr) on gaps.
      * anything else → return [any_path] unchanged (single-file mode)."""
    if not any_path:
        return []
    p = Path(any_path)
    wants_glob = "*" in any_path
    is_dir = p.is_dir()

    if not wants_glob and not is_dir:
        return [any_path]

    # ---- glob mode: honor the literal pattern ------------------------------
    if wants_glob:
        import glob as _glob
        matches = sorted(_glob.glob(any_path))
        if not matches:
            sys.stderr.write(
                f"[WARN] discover_split_files: glob {any_path!r} matched "
                f"no files.\n"
            )
            return [any_path]
        return matches

    # ---- directory mode: enumerate by run number ---------------------------
    directory = p
    run = extract_run_number(p.name)

    if run < 0 or not directory.is_dir():
        sys.stderr.write(
            f"[WARN] discover_split_files: cannot resolve run/dir from "
            f"{any_path!r} — passing through as a single file.\n"
        )
        return [any_path]

    pat = re.compile(rf"^prad_0*{run}\.evio\.(\d+)$", re.IGNORECASE)
    matched: list[tuple[int, str]] = []
    for entry in directory.iterdir():
        m = pat.match(entry.name)
        if m:
            try:
                matched.append((int(m.group(1)), str(entry)))
            except ValueError:
                pass
    matched.sort()

    if matched:
        last = matched[-1][0]
        seen = {idx for idx, _ in matched}
        missing = [i for i in range(0, last + 1) if i not in seen]
        if missing:
            miss_str = " ".join(f".{i:05d}" for i in missing)
            sys.stderr.write(
                f"[WARN] split-file gaps in run {run} (found "
                f"{len(matched)} file(s), max suffix .{last:05d}): "
                f"missing {miss_str}\n"
            )

    if not matched:
        sys.stderr.write(
            f"[WARN] discover_split_files: no files matched 'prad_{run}.evio.*' "
            f"in {directory}\n"
        )
        return [any_path]

    return [path for _, path in matched]


def project_to_z(x: float, y: float, z: float, target_z: float
                 ) -> tuple[float, float, float]:
    """Straight-line target→hit projection to z = target_z.  Mirrors
    analysis::GetProjection (single-vertex assumption).  z must be > 0."""
    if z == 0.0:
        return x, y, target_z
    s = target_z / z
    return x * s, y * s, target_z


# ---- Pipeline state + setup ----

@dataclass
class Pipeline:
    """Initialized HyCal + GEM systems plus geometry & EVIO file list.  One
    call to setup_pipeline() produces a ready-to-loop bundle.

    The detector-side fields (`cfg`, `geo`, `hycal`, `gem_sys`,
    `hycal_xform`, `gem_xforms`, `crate_map`) come from the C++
    PipelineBuilder so the wiring is identical to what the live server
    does.  The clusterers (`hc_clusterer`, `gem_clusterer`,
    `wave_ana`) are constructed here because they hold per-event scratch
    state and stay outside the C++ Pipeline."""
    cfg:           "dec.DaqConfig"       = None
    crate_map:     dict[int, int]        = field(default_factory=dict)
    geo:           "det.RunConfig"       = None
    hycal:         "det.HyCalSystem"     = None
    hc_clusterer:  "det.HyCalCluster"    = None
    wave_ana:      "dec.WaveAnalyzer"    = None
    gem_sys:       "det.GemSystem"       = None
    gem_clusterer: "det.GemCluster"      = None
    evio_files:    list[str]             = field(default_factory=list)
    hycal_xform:   "det.DetectorTransform" = None
    gem_xforms:    list                  = field(default_factory=list)

    # Matching parameters from reconstruction_config.json:matching, copied
    # off the underlying C++ Pipeline so callers can read them without a
    # second JSON parse.
    hycal_pos_res:  list[float]          = field(default_factory=lambda: [2.6, 0.0, 0.0])
    gem_pos_res:    list[float]          = field(default_factory=lambda: [0.1, 0.1, 0.1, 0.1])
    target_pos_res: list[float]          = field(default_factory=lambda: [1.0, 1.0, 20.0])

    # The underlying det.Pipeline (kept alive so the borrowed `hycal`,
    # `gem_sys`, `hycal_xform`, `gem_xforms` references stay valid).
    _core: "det.Pipeline"                = None


def setup_pipeline(*,
                   evio_path: str,
                   run_num: int = -1,
                   gem_ped_file: str = "",
                   gem_cm_file: str = "",
                   hc_calib_file: str = "",
                   daq_config: str = "",
                   gem_map_file: str = "",
                   hc_map_file: str = "",
                   ) -> Pipeline:
    """Wire up HyCal + GEM detectors via the C++ PipelineBuilder, then add
    the per-event scratch (HyCalCluster + GemCluster + WaveAnalyzer) and
    file-discovery bits the analysis loop needs.

    `_file` args accept "" to fall back to runinfo / database defaults.
    Run number defaults to a sniff from the EVIO basename when -1."""
    b = det.PipelineBuilder()
    if daq_config:    b.set_daq_config(daq_config)
    if hc_calib_file: b.set_hycal_calib(hc_calib_file)
    if gem_ped_file:  b.set_gem_pedestal(gem_ped_file)
    if gem_cm_file:   b.set_gem_common_mode(gem_cm_file)
    if hc_map_file:   b.set_hycal_map(hc_map_file)
    if gem_map_file:  b.set_gem_map(gem_map_file)
    if run_num > 0:   b.set_run_number(run_num)
    elif evio_path:   b.set_run_number_from_evio(evio_path)
    core = b.build()

    p = Pipeline()
    p._core         = core
    p.cfg           = core.daq_cfg
    p.crate_map     = {int(r.tag): int(r.crate) for r in core.daq_cfg.roc_tags}
    p.geo           = core.run_cfg
    p.hycal         = core.hycal
    p.hc_clusterer  = det.HyCalCluster(core.hycal)
    p.hc_clusterer.set_config(core.hycal_cluster_cfg)
    p.wave_ana      = dec.WaveAnalyzer(core.daq_cfg.wave_cfg)
    p.gem_sys       = core.gem
    p.gem_clusterer = det.GemCluster()
    p.hycal_xform   = core.hycal_transform
    p.gem_xforms    = list(core.gem_transforms)
    p.hycal_pos_res  = list(core.hycal_pos_res)
    # Detectors missing from matching.gem_pos_res get 0.1 mm — the same
    # fallback as AppState::gemPosRes in the server — so callers can
    # index 0..3.
    gem_res = list(core.gem_pos_res)
    p.gem_pos_res    = gem_res + [0.1] * (4 - len(gem_res))
    p.target_pos_res = list(core.target_pos_res)

    p.evio_files = discover_split_files(evio_path or "")
    print(f"[setup] EVIO       : {len(p.evio_files)} split file(s) for "
          f"input {evio_path or '(null)'}", flush=True)
    for f in p.evio_files:
        print(f"           {f}", flush=True)
    return p


def setup_pipeline_from_args(args: argparse.Namespace) -> Pipeline:
    """setup_pipeline() driven by the options of add_common_args()."""
    return setup_pipeline(
        evio_path     = args.evio_path,
        run_num       = args.run_num,
        gem_ped_file  = args.gem_ped_file,
        gem_cm_file   = args.gem_cm_file,
        hc_calib_file = args.hc_calib_file,
        daq_config    = args.daq_config,
        gem_map_file  = args.gem_map_file,
        hc_map_file   = args.hc_map_file,
    )


def load_matching_config(p: Pipeline
        ) -> tuple[tuple[float, float, float],
                   list[float],
                   tuple[float, float, float]]:
    """Return ((A, B, C), gem_pos_res, (sx, sy, sz)) from the matching
    parameters resolved by setup_pipeline(); gem_pos_res has at least 4
    entries."""
    A, B, C = p.hycal_pos_res[0], p.hycal_pos_res[1], p.hycal_pos_res[2]
    gem = list(p.gem_pos_res)
    tgt = (p.target_pos_res[0], p.target_pos_res[1], p.target_pos_res[2])
    return (A, B, C), gem, tgt


# ---- Per-event reconstruction helpers ----
# Scripts call these inside their own event loop to keep the pipeline
# idempotent — the user passes the already-decoded EventData / SspEventData
# in (so trigger filtering happens before we pay the reco cost).

def iter_fadc_channels(fadc_evt) -> Iterator[tuple]:
    """Yield (roc_tag, slot, channel, ChannelData) for every present FADC
    channel that carries samples."""
    for ri in range(fadc_evt.nrocs):
        roc = fadc_evt.roc(ri)
        if not roc.present:
            continue
        for s in roc.present_slots():
            slot = roc.slot(s)
            for c in slot.present_channels():
                cd = slot.channel(c)
                if cd.nsamples > 0:
                    yield roc.tag, s, c, cd


def iter_hycal_peaks(p: Pipeline, fadc_evt,
                     window: Optional[Tuple[float, float]] = None
                     ) -> Iterator[tuple]:
    """Yield (module, peaks) for every HyCal channel of one decoded
    EventData.  `peaks` is the list of WaveAnalyzer peaks, restricted to
    lo < time < hi when window=(lo, hi) is given."""
    for tag, s, c, cd in iter_fadc_channels(fadc_evt):
        crate = p.crate_map.get(tag)
        if crate is None:
            continue
        mod = p.hycal.module_by_daq(crate, s, c)
        if mod is None or not mod.is_hycal():
            continue
        _, _, peaks = p.wave_ana.analyze(cd.samples)
        if window is not None:
            lo, hi = window
            peaks = [pk for pk in peaks if lo < pk.time < hi]
        yield mod, peaks


def add_hycal_peaks(cl, mod, peaks, multi_pulse: bool) -> None:
    """Feed one module's peaks to a HyCalCluster: every peak with its time
    stamp when multi_pulse, else only the largest-integral one."""
    if multi_pulse:
        for pk in peaks:
            cl.add_hit(mod.index, mod.energize(pk.integral), float(pk.time))
        return
    best = None
    best_i = -1.0
    for pk in peaks:
        if pk.integral > best_i:
            best = pk
            best_i = pk.integral
    if best is not None:
        cl.add_hit(mod.index, mod.energize(best.integral), float(best.time))


def reconstruct_hycal(p: Pipeline, fadc_evt) -> list:
    """Run HyCal waveform → energy → cluster on one decoded EventData.
    Returns a list of det.ClusterHit (detector frame).

    Two pulse-selection modes, switched by `seed_time_window` on the
    HyCal cluster config:

    * `seed_time_window > 0` (multi-pulse): every analyzer-detected peak
      is fed into the clusterer with its time stamp.  HyCalCluster groups
      them via seed-anchored timing coincidence — one or more clusters
      per event, each with time-coherent constituents.
    * `seed_time_window <= 0` (default): one ModuleHit per module — the
      largest-integral detected peak — matches `bestPeak()` in
      viewer_utils.h."""
    multi_pulse = p.hc_clusterer.get_config().seed_time_window > 0.0
    p.hc_clusterer.clear()
    for mod, peaks in iter_hycal_peaks(p, fadc_evt):
        add_hycal_peaks(p.hc_clusterer, mod, peaks, multi_pulse)
    p.hc_clusterer.form_clusters()
    return p.hc_clusterer.reconstruct_hits()


def reconstruct_gem(p: Pipeline, ssp_evt) -> None:
    """Run GEM ProcessEvent + Reconstruct.  After this the per-detector
    hit lists are accessible via p.gem_sys.get_hits(d)."""
    p.gem_sys.clear()
    p.gem_sys.process_event(ssp_evt)
    p.gem_sys.reconstruct(p.gem_clusterer)


def hycal_to_lab(p: Pipeline, h) -> tuple[float, float, float]:
    """Lab-frame (x, y, z) of a HyCal ClusterHit, z at the shower-max depth
    (det.shower_depth — same calc as Replay.cpp)."""
    return p.hycal_xform.to_lab(h.x, h.y, det.shower_depth(h.center_id, h.energy))


def best_gem_matches(hc_lab: Sequence[Sequence[float]],
                     gem_lab: Sequence[Sequence[Sequence[float]]],
                     abc: Tuple[float, float, float],
                     gem_pos_res: Sequence[float],
                     nsigma: float) -> list[tuple]:
    """Closest GEM hit per (HyCal cluster, GEM detector) pair, projected
    through the target at the origin — mirrors gem_hycal_matching.C.

    hc_lab items start with lab (x, y, z, energy); gem_lab[d] items start
    with lab (x, y, z), and the detector's plane z is taken from its first
    hit.  sigma_total = sqrt((sigma_face·z_gem/z_hc)² + gem_pos_res[d]²)
    with sigma_face = hycal_pos_resolution(*abc, energy); a hit matches
    within nsigma·sigma_total, the last of equally close hits winning.

    Returns [(k, d, gi, proj_x, proj_y, residual, sigma_total, sigma_face)]
    ordered by HyCal index k, then detector d."""
    A, B, C = abc
    out = []
    for k, h in enumerate(hc_lab):
        hx, hy, hz, he = h[0], h[1], h[2], h[3]
        if hz <= 0.0:
            continue
        sig_face = hycal_pos_resolution(A, B, C, he)
        for d, gl in enumerate(gem_lab):
            if not gl:
                continue
            z_gem = gl[0][2]
            if z_gem <= 0.0:
                continue
            scale         = z_gem / hz
            proj_x        = hx * scale
            proj_y        = hy * scale
            sig_hc_at_gem = sig_face * scale
            sig_gem       = gem_pos_res[d]
            sig_total     = math.sqrt(
                sig_hc_at_gem * sig_hc_at_gem + sig_gem * sig_gem)

            best_gi = -1
            best_dr = nsigma * sig_total
            for gi, g in enumerate(gl):
                dx = g[0] - proj_x
                dy = g[1] - proj_y
                dr = math.sqrt(dx * dx + dy * dy)
                if dr <= best_dr:
                    best_dr = dr
                    best_gi = gi
            if best_gi >= 0:
                out.append((k, d, best_gi, proj_x, proj_y,
                            best_dr, sig_total, sig_face))
    return out


# ---- Argparse helpers ----

def add_common_args(ap: argparse.ArgumentParser) -> None:
    """Register the args every analysis script accepts.  Mirrors the C++
    signatures of gem_hycal_matching / plot_hits_at_hycal."""
    ap.add_argument("evio_path",
                    help="EVIO input.  glob ('prad_NNNNNN.evio.*'), directory, "
                         "or single split (prad_NNNNNN.evio.00000).")
    ap.add_argument("out_path",
                    help="Output table path (.tsv / .csv).")
    ap.add_argument("--max-events", type=int, default=0,
                    help="Stop after N raw physics events (0 = all).")
    ap.add_argument("--run-num", type=int, default=-1,
                    help="Run number override (default -1 = sniff from filename).")
    ap.add_argument("--gem-ped-file",  default="",
                    help='GEM pedestal file (default "" = via runinfo).')
    ap.add_argument("--gem-cm-file",   default="",
                    help='GEM common-mode file (default "" = via runinfo).')
    ap.add_argument("--hc-calib-file", default="",
                    help='HyCal calibration file (default "" = via runinfo).')
    ap.add_argument("--daq-config",    default="",
                    help='DAQ config (default "" = installed default).')
    ap.add_argument("--gem-map-file",  default="",
                    help='GEM map (default "" = database/gem_map.json).')
    ap.add_argument("--hc-map-file",   default="",
                    help='HyCal map (default "" = database/hycal_map.json).')
    ap.add_argument("--csv", action="store_true",
                    help="Emit CSV instead of TSV.")
    ap.add_argument("--no-header", action="store_true",
                    help="Skip the column-name header row.")


# Physics trigger bits accepted for reconstruction.  Mirrors the server's
# `physics.accept_trigger_bits` set in monitor_config.json — currently
# {SSP0, SSP1, SSP2, SSP3} = bits {8,9,10,11} = mask 0xf00.  The check is
# bitwise-AND (an event passes if ANY accepted bit is set), matching
# `TriggerFilter::operator()` in src/app_state.h.
PHYSICS_TRIGGER_MASK = 0xf00


def passes_physics_trigger(trigger_bits: int) -> bool:
    """C++-equivalent physics trigger gate (`bits & PHYSICS_TRIGGER_MASK`)."""
    return bool(trigger_bits & PHYSICS_TRIGGER_MASK)


# SSP0-only trigger word.  An exact `bits == PHYSICS_TRIGGER_BITS` gate keeps
# only pure-SSP0 events, a subset of what passes_physics_trigger() (and the
# server) accepts: the cluster-counting SSP1-3 triggers are excluded.
PHYSICS_TRIGGER_BITS = 0x100


# ---- EVIO event loop ----

@dataclass
class LoopStats:
    """Counters maintained by iter_physics_events()."""
    n_files_open: int   = 0
    n_read:       int   = 0     # EVIO records read
    n_phys:       int   = 0     # decoded physics events (before the trigger gate)
    n_kept:       int   = 0     # physics events that passed the trigger gate
    cur_file:     str   = ""    # EVIO file being read
    t0:           float = field(default_factory=time.monotonic)

    @property
    def elapsed(self) -> float:
        return time.monotonic() - self.t0


def iter_physics_events(p: Pipeline, stats: LoopStats, *,
                        with_ssp: bool,
                        max_events: int = 0,
                        accept: Optional[Callable[[int], bool]] = passes_physics_trigger,
                        progress: Optional[Callable[[LoopStats], None]] = None,
                        progress_every: int = 5000) -> Iterator[tuple]:
    """Yield (fadc_evt, ssp_evt) for every decoded physics event in
    p.evio_files whose trigger bits pass `accept` (None = no gate).
    ssp_evt is None unless with_ssp.

    Stops right after the `max_events`-th physics event, counted before
    the trigger gate (<= 0 = all).  After each physics record,
    progress(stats) is called whenever n_phys has crossed the next
    multiple of `progress_every`.  Iterate it directly in a `for` so the
    EVIO channel is closed as soon as the loop is left."""
    ch = dec.EvChannel()
    ch.set_config(p.cfg)
    step = max(1, int(progress_every))
    next_progress = step
    for fpath in p.evio_files:
        if ch.open_auto(fpath) != dec.Status.success:
            print(f"[WARN] skip (cannot open): {fpath}", flush=True)
            continue
        stats.n_files_open += 1
        stats.cur_file = fpath
        print(f"[file {stats.n_files_open}/{len(p.evio_files)}] {fpath}",
              flush=True)
        try:
            while ch.read() == dec.Status.success:
                stats.n_read += 1
                if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
                    continue
                for i in range(ch.get_n_events()):
                    decoded = ch.decode_event(i, with_ssp=with_ssp)
                    if not decoded["ok"]:
                        continue
                    stats.n_phys += 1
                    fadc_evt = decoded["event"]
                    if accept is None or accept(int(fadc_evt.info.trigger_bits)):
                        stats.n_kept += 1
                        yield fadc_evt, decoded["ssp"]
                    if max_events > 0 and stats.n_phys >= max_events:
                        return
                if progress is not None and stats.n_phys >= next_progress:
                    progress(stats)
                    while next_progress <= stats.n_phys:
                        next_progress += step
        finally:
            ch.close()


# ---- Output helpers ----

def open_table_writer(out_path: str, csv_mode: bool):
    """Open `out_path` for writing.  Returns (file_handle, write_row callable)
    where write_row(seq) writes one row.  Caller closes the file."""
    import csv as _csv
    fh = open(out_path, "w", encoding="utf-8", newline="")
    if csv_mode:
        w = _csv.writer(fh, lineterminator="\n")
        return fh, w.writerow
    sep = "\t"
    def write_row(row):
        fh.write(sep.join("" if v is None else str(v) for v in row))
        fh.write("\n")
    return fh, write_row


def print_summary(p: Pipeline, stats: LoopStats, kept_label: str,
                  rows: Sequence[Tuple[str, object]], out_path: str) -> None:
    """Print the end-of-run '--- summary ---' block: loop counters (the
    trigger-passed count labelled `kept_label`), the script's own
    (label, value) rows, elapsed time and the output path."""
    print("--- summary ---", flush=True)
    for label, value in [
            ("EVIO files opened", f"{stats.n_files_open} / {len(p.evio_files)}"),
            ("EVIO records",      stats.n_read),
            ("physics events",    stats.n_phys),
            (kept_label,          stats.n_kept),
            *rows,
            ("elapsed (s)",       f"{stats.elapsed:.2f}"),
            ("wrote",             out_path)]:
        print(f"  {label:<22}: {value}")


def import_pyplot():
    """Headless (Agg) matplotlib.pyplot, or None if matplotlib is unavailable."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except ImportError:
        return None
