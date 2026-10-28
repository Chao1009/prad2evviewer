#!/usr/bin/env python3
"""
gem_hycal_matching.py — Python counterpart of analysis/scripts/gem_hycal_matching.C

Same pipeline (HyCal reco → GEM reco → straight-line target-vertex matching),
same trigger filter (trigger_bits == 0x100).  Difference vs. the ROOT script:

  * No ROOT.  Output is a flat per-match TSV / CSV table — one row per
    (event, HyCal cluster, GEM detector) tuple, with `det_id` distinguishing
    which GEM (0..3) won the match.

  * Best-match rule (HyCal cluster as baseline):
      For each (HC cluster, GEM detector) pair, keep at most ONE row —
      the GEM hit with the smallest 2D residual that's still inside the
      `--match-nsigma · σ_total` window.  A given GEM hit can win against
      multiple HC clusters (no GEM-side exclusivity).

Matching geometry (lab frame, target at origin, beam along +z):

    σ_hc_face = sqrt((A/sqrt(E_GeV))^2 + (B/E_GeV)^2 + C^2)  [mm at HyCal face]
    σ_hc@gem  = σ_hc_face · (z_gem / z_hc)
    σ_gem     = gem_pos_res[det_id] mm                       (per detector)
    σ_total   = sqrt(σ_hc@gem² + σ_gem²)
    cut       = nsigma · σ_total

    (A, B, C) and gem_pos_res come from reconstruction_config.json:matching.

Output columns (one row per matched pair):

  event_num, trigger_bits,
  hc_idx, hc_x, hc_y, hc_z, hc_energy, hc_center, hc_nblocks, hc_sigma,
  det_id,
  gem_x, gem_y, gem_z,                 # lab/target-centered (mm)
  gem_x_local, gem_y_local,            # detector-frame (mm)
  gem_x_charge, gem_y_charge,          # X/Y cluster total ADC
  gem_x_peak,   gem_y_peak,            # X/Y cluster max-strip ADC
  gem_x_max_tb, gem_y_max_tb,          # time sample of max-ADC strip (int)
  gem_x_size,   gem_y_size,            # X/Y cluster strip count
  proj_x, proj_y, residual, sigma_total

Coordinates labelled "lab" are target-centered (mm); hc_z includes
shower-depth.  The "_local" coords are the GEM detector frame (no
rotation/translation), useful for per-detector hit maps.  Convert
gem_*_max_tb to ns by multiplying by the cluster config's ts_period
(default 25 ns).

Usage
-----
  # full run (glob):
  python analysis/pyscripts/gem_hycal_matching.py \\
      /data/stage6/prad_023867/prad_023867.evio.* match_023867.tsv

  # single split (debugging):
  python analysis/pyscripts/gem_hycal_matching.py \\
      /data/stage6/prad_023867/prad_023867.evio.00000 match_023867_seg0.tsv

  # CSV output, tighter cut, cap at 50k events:
  python analysis/pyscripts/gem_hycal_matching.py input.evio.* out.csv \\
      --csv --match-nsigma 2.0 --max-events 50000
"""

from __future__ import annotations

import argparse
import sys

import _common as C


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    C.add_common_args(ap)
    ap.add_argument("--match-nsigma", type=float, default=3.0,
                    help="Matching window in σ_total (default 3.0).")
    args = ap.parse_args(argv)

    p = C.setup_pipeline_from_args(args)
    print(f"[setup] Match cut  : {args.match_nsigma:.2f} · sigma_total",
          flush=True)

    (pr_A, pr_B, pr_C), gem_pos_res, _ = C.load_matching_config(p)
    print(f"[setup] HC sigma(E)= sqrt(({pr_A:.3f}/sqrt(E_GeV))^2"
          f"+({pr_B:.3f}/E_GeV)^2+{pr_C:.3f}^2) mm", flush=True)
    print(f"[setup] GEM sigma  : {gem_pos_res} mm", flush=True)

    cols = [
        "event_num", "trigger_bits",
        "hc_idx", "hc_x", "hc_y", "hc_z", "hc_energy",
        "hc_center", "hc_nblocks", "hc_sigma",
        "det_id",
        "gem_x", "gem_y", "gem_z",
        "gem_x_local", "gem_y_local",
        "gem_x_charge", "gem_y_charge",
        "gem_x_peak",   "gem_y_peak",
        "gem_x_max_tb", "gem_y_max_tb",
        "gem_x_size",   "gem_y_size",
        "proj_x", "proj_y", "residual", "sigma_total",
    ]
    fh, write_row = C.open_table_writer(args.out_path, args.csv)
    if not args.no_header:
        write_row(cols)

    n_match = 0
    total_hc = 0
    total_gem = 0
    gem_per_det = [0, 0, 0, 0]

    stats = C.LoopStats()
    try:
        for fadc_evt, ssp_evt in C.iter_physics_events(
                p, stats, with_ssp=True, max_events=args.max_events,
                accept=lambda b: b == C.PHYSICS_TRIGGER_BITS,
                progress=lambda st: print(f"[progress] {st.n_phys} physics events",
                                          flush=True)):
            event_num    = int(fadc_evt.info.event_number)
            trigger_bits = int(fadc_evt.info.trigger_bits)

            # ---- HyCal: waveform → energy → cluster --------------
            hc_raw = C.reconstruct_hycal(p, fadc_evt)

            # Lab-frame HyCal hits with shower-depth applied to z.
            hc_lab = [(*C.hycal_to_lab(p, h), float(h.energy),
                       int(h.center_id), int(h.nblocks)) for h in hc_raw]
            total_hc += len(hc_lab)

            # ---- GEM: pedestal → CM → ZS → 1D + 2D --------------
            C.reconstruct_gem(p, ssp_evt)

            # Per-detector lab-frame hit lists, plus the raw GEMHit
            # for charge / size / peak / timing lookup at write time.
            # Tuple layout (positional, frozen):
            #   0: x_lab     1: y_lab     2: z_lab
            #   3: x_local   4: y_local
            #   5: x_charge  6: y_charge
            #   7: x_peak    8: y_peak
            #   9: x_max_tb 10: y_max_tb
            #  11: x_size   12: y_size
            gem_lab: list[list[tuple]] = [[], [], [], []]
            n_dets = min(p.gem_sys.get_n_detectors(), 4)
            for d in range(n_dets):
                raw = p.gem_sys.get_hits(d)
                gem_per_det[d] += len(raw)
                total_gem      += len(raw)
                xform = p.gem_xforms[d]
                for g in raw:
                    x, y, z = xform.to_lab(g.x, g.y)
                    gem_lab[d].append((
                        x, y, z,
                        float(g.x), float(g.y),
                        float(g.x_charge), float(g.y_charge),
                        float(g.x_peak),   float(g.y_peak),
                        int(g.x_max_timebin), int(g.y_max_timebin),
                        int(g.x_size),     int(g.y_size),
                    ))

            # ---- best match per HC × GEM detector ---------------
            for k, d, gi, proj_x, proj_y, residual, sig_total, sig_face in \
                    C.best_gem_matches(hc_lab, gem_lab, (pr_A, pr_B, pr_C),
                                       gem_pos_res, args.match_nsigma):
                hx, hy, hz, he, hc_center, hc_nblocks = hc_lab[k]
                g = gem_lab[d][gi]
                write_row([
                    event_num, trigger_bits,
                    k,
                    f"{hx:.4f}", f"{hy:.4f}", f"{hz:.4f}",
                    f"{he:.4f}",
                    hc_center, hc_nblocks,
                    f"{sig_face:.4f}",
                    d,
                    f"{g[0]:.4f}", f"{g[1]:.4f}", f"{g[2]:.4f}",
                    f"{g[3]:.4f}", f"{g[4]:.4f}",
                    f"{g[5]:.4f}", f"{g[6]:.4f}",
                    f"{g[7]:.4f}", f"{g[8]:.4f}",
                    g[9], g[10],
                    g[11], g[12],
                    f"{proj_x:.4f}", f"{proj_y:.4f}",
                    f"{residual:.4f}", f"{sig_total:.4f}",
                ])
                n_match += 1
    finally:
        fh.close()

    C.print_summary(p, stats, "passed trig cut 0x100", [
        ("total HyCal clusters", total_hc),
        ("total GEM 2D hits",
         f"{total_gem}  (det0={gem_per_det[0]} det1={gem_per_det[1]} "
         f"det2={gem_per_det[2]} det3={gem_per_det[3]})"),
        ("matched rows written", n_match),
    ], args.out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
