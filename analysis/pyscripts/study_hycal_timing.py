#!/usr/bin/env python3
"""
study_hycal_timing.py — sample (seed, neighbour, dt) pairs from real HyCal
waveform data so you can pick a value for HyCalClusterConfig.seed_time_window.

Per event:
  1. Decode FADC, run WaveAnalyzer per HyCal channel.
  2. Push EVERY detected peak (within an optional pre-window) into the
     clusterer with its time stamp — multiple pulses per module land as
     separate ModuleHits.
  3. Call HyCalCluster.collect_neighbor_timing(max_qdist) — for each seed
     candidate (largest pulse in a region passing min_center_energy that
     hasn't already seeded another cluster in this scan), emit one row
     per neighbouring pulse within `max_qdist` module units.  No timing
     cut is applied; no neighbour pulses are consumed across seeds.

Output is a flat TSV / CSV — one row per (seed, neighbour) pair — meant
to be plotted as histograms (dt; dt vs. dist_q; dt vs. neighbour energy)
in your tool of choice.  Once a sensible cut shows up, set
`seed_time_window` in your reconstruction config (or via the C++ API)
and the production clusterer will start gating on it.

Output columns
--------------
  event_num, trigger_bits,
  seed_module, seed_id, seed_time, seed_energy,
  neighbor_module, neighbor_id, neighbor_time, neighbor_energy,
  dt, dx_q, dy_q, dist_q

`*_id` are PrimEx IDs (G: 1-576, W: 1001-2152).  Times are ns;
energies are MeV; dx_q / dy_q / dist_q are in module units.

Usage
-----
  # full run, default 5-module radius:
  python analysis/pyscripts/study_hycal_timing.py \\
      /data/stage6/prad_023867/prad_023867.evio.* hcal_dt_023867.tsv

  # tighter spatial radius, narrower pre-window, cap at 50k events:
  python analysis/pyscripts/study_hycal_timing.py input.evio.* dt.tsv \\
      --max-qdist 3 --pre-window 80 250 --max-events 50000
"""

from __future__ import annotations

import argparse
import math
import sys

import _common as C


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    C.add_common_args(ap)
    ap.add_argument("--max-qdist", type=float, default=5.0,
                    help="Spatial radius around the seed module in quantized "
                         "module units (default 5.0 — captures ~9x9 around the "
                         "seed in PbWO4, more in PbGlass).")
    ap.add_argument("--pre-window", type=float, nargs=2, default=None,
                    metavar=("LO", "HI"),
                    help="Reject peaks outside [LO, HI] ns BEFORE feeding the "
                         "clusterer (useful to suppress out-of-trigger noise "
                         "pulses).  Default: no pre-window — all detected "
                         "peaks are kept so the dt landscape is unbiased.")
    ap.add_argument("--require-trigger", type=lambda s: int(s, 0),
                    default=C.PHYSICS_TRIGGER_BITS,
                    help=f"Trigger bits gate (default 0x{C.PHYSICS_TRIGGER_BITS:x}).  "
                         "Pass 0 to accept any trigger.")
    args = ap.parse_args(argv)

    p = C.setup_pipeline_from_args(args)
    print(f"[setup] max_qdist  : {args.max_qdist} module units", flush=True)
    if args.pre_window is not None:
        print(f"[setup] pre_window : [{args.pre_window[0]}, {args.pre_window[1]}] ns",
              flush=True)
    else:
        print("[setup] pre_window : (none — all detected peaks kept)", flush=True)

    cluster_cfg = p.hc_clusterer.get_config()
    print(f"[setup] HC cluster : min_mod_E={cluster_cfg.min_module_energy} "
          f"min_ctr_E={cluster_cfg.min_center_energy} "
          f"min_cl_E={cluster_cfg.min_cluster_energy}", flush=True)

    cols = [
        "event_num", "trigger_bits",
        "seed_module", "seed_id", "seed_time", "seed_energy",
        "neighbor_module", "neighbor_id", "neighbor_time", "neighbor_energy",
        "dt", "dx_q", "dy_q", "dist_q",
    ]
    fh, write_row = C.open_table_writer(args.out_path, args.csv)
    if not args.no_header:
        write_row(cols)

    n_seeds = n_rows = 0

    stats = C.LoopStats()
    try:
        for fadc_evt, _ in C.iter_physics_events(
                p, stats, with_ssp=False, max_events=args.max_events,
                accept=((lambda b: b == args.require_trigger)
                        if args.require_trigger else None),
                progress=lambda st: print(f"[progress] {st.n_phys} physics events  "
                                          f"seeds={n_seeds}  rows={n_rows}",
                                          flush=True)):
            trigger_bits = int(fadc_evt.info.trigger_bits)
            event_num = int(fadc_evt.info.event_number)

            # Push every detected peak into the clusterer.  We do NOT
            # call form_clusters() — collect_neighbor_timing() works
            # directly on the accumulated hits and applies its own
            # seed-finding logic without consuming pulses across seeds.
            p.hc_clusterer.clear()
            for mod, peaks in C.iter_hycal_peaks(p, fadc_evt, args.pre_window):
                C.add_hycal_peaks(p.hc_clusterer, mod, peaks, True)

            rows = p.hc_clusterer.collect_neighbor_timing(args.max_qdist)
            if not rows:
                continue

            # Track the seed-module set so we can count distinct seeds
            # — rows is one entry per (seed, neighbour), not per seed.
            seen_seeds = set()
            for r in rows:
                seed_mod = p.hycal.module(r.seed_module)
                nbr_mod  = p.hycal.module(r.neighbor_module)
                seen_seeds.add(r.seed_module)
                dist_q = math.sqrt(r.dx_q * r.dx_q + r.dy_q * r.dy_q)
                write_row([
                    event_num, trigger_bits,
                    r.seed_module, seed_mod.id,
                    f"{r.seed_time:.3f}", f"{r.seed_energy:.4f}",
                    r.neighbor_module, nbr_mod.id,
                    f"{r.neighbor_time:.3f}", f"{r.neighbor_energy:.4f}",
                    f"{r.dt:.3f}",
                    f"{r.dx_q:.4f}", f"{r.dy_q:.4f}", f"{dist_q:.4f}",
                ])
                n_rows += 1
            n_seeds += len(seen_seeds)
    finally:
        fh.close()

    C.print_summary(p, stats, "trigger-passed events", [
        ("total seeds",  n_seeds),
        ("rows written", n_rows),
    ], args.out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
