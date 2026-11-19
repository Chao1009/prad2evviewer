---
domain: science
summary: Exotic waveforms stay event markers; only clusters with a clipped seed holding >0.95 of the energy are vetoed, on 2.24 and 3.49 GeV H2; channel drop rejected.
affects: [analysis/pyscripts/waveform_analysis/pileup_analysis/exotic_events/, analysis/pyscripts/waveform_analysis/writeup/waveform_analysis_prad2_x17.tex]
verified: 2026-11-18
source: analysis/pyscripts/notes/plan_pileup_fit_validation_20261116.md#exotic-events-impact-packages-approved-2026-11-18
supersedes:
---
Four user-picked types (tail below baseline, wide smooth pulse, wide pulse fitted as close pile-up, clipped pulse with secondary pulses) plus saturation become exclusive flags, priority f4>f5>f1>f2>f3, one threshold set for all run types. f1 baseline shift >8σ and tail <−5σ; f2 w50>3 (template FWHM 23.5 ns), unclipped; f3 k≥2 and stretched-single BIC; f4 clipped and (secondary ≤8 ns from clip or E_sec>0.5 E_beam); f5 clipped, not f4. PbWO4 only, run 24937 excluded. Any-flag rate 0.069/0.090/0.073/1.114 % per kept channel (0.73/2.24/3.49 GeV H2/Ta); saturation ≤0.1 % of full-energy pulses, statistics too small.
No σ/E or peak shift beyond 2σ under veto, correction or drop. The gate vetoes only f4/f5-seeded clusters with seed fraction >0.95, on 2.24/3.49 GeV H2 (cost 0.0044/0.0067 %, 37 fakes removed); BASE_SLOPE and CLIP_SINGLE stay off; f5 stays for the full scan.
Rejected: wider veto (more cost); the user's channel drop (no gain, costs 0.03–0.10 % H2, 0.87 % Ta, shifts Ta tail and Møller peak).
