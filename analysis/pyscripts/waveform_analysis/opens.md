# Open questions and pending items (as of 2026-11-11, end of day)

State of the code: everything under `waveform_analysis/` and `crosstalk_study/` is standalone,
untracked and uncommitted (also `mf_pulse.py` and `notes/PbWO4_energy_resolution_20261108.md`).
Nothing has been committed all week.

## Running / just launched (results not yet read)

1. Energy benchmark comparison with the old chain (`benchmark_hycal_timing_test_20261021.py`):
   reference numbers from the notes, item-by-item difference table, old chain rerun on run 25320
   (gain refinement off, position correction off) for an apples-to-apples number.
2. Time-gate scan of the benchmark, 2–24 ns in 0.5 ns steps (report-only rerun on the existing
   scans); driver `run_energy_resolution.sh`, `GATES` line.

## Energy resolution benchmark (`waveform_analysis/energy_resolution_benchmark.py`)

3. Resolutions are worse than the previous studies (2.24 GeV: 2.004 % legacy; run 25320 alone gave
   1.959 % yesterday before the 10 ADC hit threshold, the 1 sigma baseline clip, the noise
   reference and the sigma(t, A) scaling). Suspects, in order: per-crystal gain refinement and the
   position correction are not implemented (parked as a later package because the reference
   output trees were deleted); no amplitude correction f(A) (energy scale +0.6 %); the hit
   threshold; the baseline clip. Decide what to add after item 1 reports.
   11-11 late: per-crystal gains (on by default, `GAIN_REFINE`) and the position correction added to
   `report` (gate 8/12/16 check, scratch output): legacy sigma/E CB none -> gains + position 2.004 ->
   1.878 % (2.24 GeV H2; 25320 alone 1.894 %), 2.320 -> 1.700 % (3.49), 3.744 -> 3.385 % (0.73), 2.408 ->
   2.132 % (Ta); mu/E' - 1 after the gains -0.01 to +0.14 %; still above the 11-05 reference (1.794 /
   1.591 % at 2.24 / 3.49 GeV, with Moller rejection). The X17 runs need a rescan (new physics-any gate).
4. Energy scale sits 2–4 % high in every run type (mu/E' - 1): calibration basis / amplitude
   correction, not a pulse-chain effect.
5. At 0.73 GeV the Gaussian +-2 sigma width (8–10 %) disagrees with the Crystal Ball (3.7 %): the
   Gaussian fit is unstable there; check the fit range.
6. The benchmark's hit list has no clipped-sample count and no chi-squared, so the saturated
   (>= 3 clipped samples) and bad-fit (chi-squared per degree of freedom > 1000) exclusions of the
   pile-up study are not applied there yet: a scan-level change.
7. Benchmark has no 10 ADC hit floor (hits down to the matched filter's own 5 ADC); decided to
   leave it, revisit if the Moller-energy clusters need it.

3a. Comparison result (2026-11-11 night): old chain on run 25320 with the same selection
    (no gains, no position correction, no Moller rejection) gives 2.015 +- 0.009 %, new chain
    2.004 %: the pulse chain is equivalent. The 1.794 % reference had gains + position correction
    + Moller rejection. Gap decomposition at 2.24 GeV: gains ~ +0.18 pp, position correction
    ~ +0.05 pp (+0.14 pp at 3.49 GeV), Moller rejection -0.01 pp. Both corrections being added
    (gains enabled, replacing f(A)); evidence /tmp/sup-bench/research-energy-resolution-gap.md.
3b. X17 benchmark scans are starved (55-90 k events per run vs ~1 M for H2): the benchmark gates
    on the exact physics bit 0x100 while X17 triggers on other physics bits (template extraction
    uses "physics-any" for X17). Make the gate run-type dependent, rescan the three X17 runs.
3c. 3.49 GeV: run 24937 (4 splits) yields no usable modules in the benchmark; runs 24698 and
    24911 differ by 2.3 % in peak position (gain epoch, absorbed by the gain refinement); the
    reference used the old runs 24386/24503/24512, so the numbers are not directly comparable.
3d. Non-linearity per-member correction (`--nl-members on` in the old chain) matters only at
    0.73 and 3.49 GeV (June7 calibration file).  The new script's clusterer already applied it
    (HyCalCluster non_linear_corr is on by default; no separate per-member step is missing); since
    2026-11-12 the `NONLINEAR` switch and the arms of `report` show its effect, and `MOLLER_REJECT`
    adds the old chain's kinematic Moller rejection.

3e. Full benchmark rerun (2026-11-12 02:38, X17 gate fixed: 2.52 M of 2.54 M X17 events now
    used; gains + position correction; 45 gates): Crystal Ball sigma/E, both corrections,
    legacy -> best gate: 0.73 GeV 3.385 -> 3.351 % (12.5 ns); 2.24 GeV 1.878 -> 1.857 % (15.5 ns);
    3.49 GeV 1.700 -> 1.679 % (20.5 ns, but yield 96.9 %); X17 2.364 -> 2.087 % (7.5 ns, yield
    107 %). Energy scale within +-0.2 % (legacy); the gated arms shift mu by +0.1 to -0.5 %.
    Caveats: yield is reported for the uncorrected arm only; "best gate" = minimum sigma/E, no
    significance test; the remaining gap to the old references (1.794 / 1.591 %) is the Moller
    rejection, not applied here.

## Pile-up study (`waveform_analysis/pileup_analysis/`)

8. The free-baseline offset hits the +-1 sigma clip for 13–24 % of channels (expected ~16 % from
   statistics, so partly real offsets, e.g. tails of pulses before the window); a flag instead of
   a clip may be the better treatment. Map of the clipped fraction not yet looked at.
9. H2 cut rates (0.96 / 0.55 / 0.28 %) are upper limits: the surviving extras sit 10–14 ns after
   the main pulse with r just above 0.1, the tail of the pulse-to-pulse shape fluctuation.
10. Synthetic efficiency at r = 0.05 is 52–82 % at |dt| = 12 ns and 0 at 8 ns (per-pulse 5 sigma
    significance); r = 0.1 at +8 ns is 77 %. Accepted for now; quote 12 ns as the first resolved
    gap for small secondaries.
11. W335 (0.73 GeV runs 24590/24603) is a ringing channel: dead-list candidate for those runs
    (database edit, user decision pending).
12. X17 beam-hole ring: 9 % truncated main pulses, 1.2 % bad waveforms (W527/W528/W563); the
    hot_modules.pdf ring example with chi-squared 839 (saturated pulse at 330 ns) — now excluded
    by the saturation flag; worth a look at what remains.
13. The generator shares one smear realisation between main and secondary pulse (the secondary's
    relative scatter is that of A1, not A2); the matched filter is conservative on the secondary.

## Templates (`waveform_analysis/template_extract.py`, `pulse_templates/`)

14. Pool selection order (elastic >= 20 -> general 200–1000 ADC band >= 20 -> small elastic ->
    single largest -> run-type average) "still needs some thought" (user). Measured effect: pools
    change the fall time by ~1 ns and the 20–45 ns flank by <= 1–2 % of peak; rise time unchanged;
    the user accepts a 1 % difference for now.
15. Single-pulse templates (50 / 5 / 29 / 118 modules per run type) are noisy by construction; the
    run-type average would be smoother for those.
16. Borrowed sigma(t) for ~500 outer modules at 0.73 GeV and X17 (fewer than 50 band pulses):
    data-limited.
17. The 4 ns shoulder in W493/W661/W663/W664 is real (confirmed by the sample-phase split); the
    wiggle check is skipped for pools >= 100.
18. Per-pulse rise time on raw 4 ns samples is 6.6 ns vs 5.35 ns for the template (sampling
    interpolation); the rise/fall jitter maps carry that effect.
19. A full critic review of the final `template_extract.py` has not been run (each step was
    reviewed lightly); a last concision pass before freezing/committing.

## Cross-talk study (`crosstalk_study/`)

20. Established: the late component (constant-fraction time lagging the matched-filter onset by
    > 8 ns) is an induced signal from the electronically adjacent digitizer channel (same crate
    and slot, channel +-1), ~1.1–1.3 % of the aggressor, peaking ~20 ns after its onset, ~15 ns
    wide; two boards (crate 2 slot 5, crate 4 slot 8) carry 33–40 % of late hits with fixed
    partner pairs on a crate-3 board (W667<-W664, W693<-W661, W761<-W729, W829<-W797,
    W660<-W662): coupling upstream of the digitizers (cables/patch panel), inferred, not tested.
21. 24–40 % of late hits have no identified partner: a full partner scan over all module pairs,
    then check the pairs against the cable/patch-panel map.
22. Possible veto for the analysis: partner (adjacent digitizer channel) > ~500 ADC and victim
    delay 10–19 ns; affects 10–20 ADC hits and isolated 10–100 MeV clusters, not elastic seeds.
23. The 10–50 MeV seed chi-squared second population at 2.24 / 3.49 GeV (2.7 / 6.6 % of seeds,
    4–5 ns late, adjacent channel firing 3–4x more often) matches the cross-talk signature; raw
    waveforms of these seeds next to their adjacent channel would confirm it.
24. The two one-off analyses (`seed_chi2_two_peaks_10_50MeV.*`, `late_component.*`) are not
    produced by `crosstalk_study.py`; fold them in if they should be reproducible.
25. Slices above the beam energy (few hundred entries, chi-squared 6–135) are artefacts; add a
    cut or a note in the figure.
26. The lone-channel veto is nearly inert at the 10 ADC threshold (2 085 / 10 490 / 20 438 / 3
    channels per run type); the burst modules are now caught by other flags — check that the
    hot-module list no longer needs it.

## Reference

- Beam hole = 2x2 PbWO4 positions (rows 17-18, cols 17-18 of the 34x34 grid).
  Ring 1 (12): W526 W527 W528 W529 W560 W563 W594 W597 W628 W629 W630 W631.
  Ring 2 (20): W491-W496, W525, W530, W559, W564, W593, W598, W627, W632, W661-W666.
  No dead module in either ring (source: database/hycal_map.json, general.json).

## Other

27. time_study/ rebuild parked (its inputs `output/module_timing_20261030` were deleted).
28. Cleanup block of 2026-11-11 evening (scratch dirs under /tmp/sup-*, `output/*_prev*`,
    chain scripts and watch logs) not yet run.
29. Commit everything (user to say when).
