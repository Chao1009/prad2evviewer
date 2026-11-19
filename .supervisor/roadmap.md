# Roadmap

1. Energy resolution improvement program (PRad-II HyCal waveform analysis)
   - Start here next session: analysis/pyscripts/notes/week_plan_20261117.md (next-week plan; may still be being written).
   - Plan: analysis/pyscripts/notes/plan_energy_resolution_improvement_20261112.md (15 work packages; the plan stays in notes/, not copied here).
   - Research input: analysis/pyscripts/notes/energy_resolution_improvement_options_20261112.md
   - Background: P0 and P13 are done (2026-11-12). As of 2026-11-14 control v3 is the current chain; P16 and P10-M1 are adopted; the P7 U2 gap clause is waived for ring-scoped constants; the P11 outermost-layer cut is not applied; run 24937 is dropped. The amplitude-calibration track (analysis/pyscripts/notes/plan_amplitude_calibration_20261114.md) is parked after P1b, resume at P2 with cfd7/sum5. P14 snake scan: farm pass approved, plan analysis/pyscripts/notes/plan_p14_snake_scan_20261114.md. The uncommitted `pileup_rates.py` diff was reviewed as safe-additive.
   - Operating note: token economy; read only the week plan and the plan's Todo and Log tails; use Mode 2 for bookings and lookups; Review: light for mechanical work.
   - Decisions:
     - .supervisor/memory/decisions/2026-11-14-control-v3-is-the-energy-resolution-chain.md
     - .supervisor/memory/decisions/2026-11-14-shape-factor-refit-adopted-at-all-h2-types.md
     - .supervisor/memory/decisions/2026-11-14-i-phi-correction-m1-adopted-per-energy.md
     - .supervisor/memory/decisions/2026-11-14-amplitude-calibration-track-parked.md
     - .supervisor/memory/decisions/2026-11-14-snake-scan-uses-lms-ratio-no-tagger.md

2. Exotic waveforms: open follow-ups
   - Plan: analysis/pyscripts/notes/plan_pileup_fit_validation_20261116.md (open items T3a gate decision on per-module shape error, T4 to T9).
   - Background: the exotic study is closed with the gate stored below. Open: full ifarm scan with PILEUP_EXOTIC=1 (fetch recipe analysis/pyscripts/notes/rcdb_and_fetch_20261117.md; ifarm login broken); Ta block-H gains need a calibration that frees the inner ring; X17 number to quote (ep 1.98 % vs Møller 2.09 %); amplitude-carrying H2 full rescan; user compile and page check of the chapter.
   - Decisions:
     - .supervisor/memory/decisions/2026-11-18-exotic-waveforms-study-2026-11-17-18.md
