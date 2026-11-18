# Fitted integral vs raw sum for single monotonic pulses (item 16)

**Date:** 2026-11-17. **Status:** closed for now. **Scope:** offline check on cached waveforms; no fitter change. Paths are relative to `analysis/pyscripts/`; `PA` = `waveform_analysis/pileup_analysis/`.

## Question

Asked by the user: for waveforms identified as a single monotonic peak, is the fitted integral A·I_m consistent with the raw pedestal-subtracted sum? If yes, a large chi2/dof on such a pulse is acceptable for the energy.

## Data and definitions

- 178,253 cached single pulses (k = 1, no overflow, A 50-3000 ADC) from the 2-file sample: 13 runs × 2 splits, PbWO4 first 1500 events per split, lead glass 60,000 (`output/resid_study/resid_cache_N60000_P1500.npz`). Default fit, MONO_EDGE and NOISE_FLAG off.
- Monotonic: the PILEUP_MONO_EDGE definition applied offline (rise over P−3..P−1, fall over P..P+9, steps within 2√2 σ_ref).
- I_m = `iint_wa` of the module's bank row: the phase-averaged walk-out sum of the unit template (0.1 × peak, two samples below end the walk).
- Raw sum, def. 1: sum(Y − ped) over the same walk-out window on the fitted pulse. Def. 2: fixed window [t0 − 20, t0 + 200] ns, which covers ~85 % of the template (ratio ≈ 0.853 baseline); not meaningful for lead glass.
- Ratio = A·I_m / raw sum. Complete window = the 32-sample fit window lies inside the record; the rest are the truncated class.

## Results

Monotonic, complete window, chi2/dof < 5, def. 1 (`PA/integral_vs_rawsum_20261117.tsv`):

| run type | n | median | rms | within ±2 % |
|---|---|---|---|---|
| 0.73 GeV H2 | 12,906 | 0.9954 | 0.011 | 89.8 % |
| 2.24 GeV H2 | 63,923 | 0.9941 | 0.013 | 86.6 % |
| 3.49 GeV H2 | 35,139 | 0.9944 | 0.012 | 87.6 % |
| 2.24 GeV Ta X17 | 51,993 | 0.9961 | 0.013 | 87.6 % |

All run types by chi2/dof (monotonic, complete window; def. 2 median for reference, baseline 0.853):

| chi2/dof | n | median def. 1 | rms def. 1 | within ±2 % | median def. 2 |
|---|---|---|---|---|---|
| < 1 | 117,293 | 0.9955 | 0.010 | 92.9 % | 0.853 |
| 1-2 | 40,448 | 0.9926 | 0.016 | 76.0 % | 0.853 |
| 2-5 | 6,220 | 0.9906 | 0.027 | 57.8 % | 0.842 |
| 5-20 | 583 | 0.9909 | 0.037 | 49.4 % | 0.838 |
| 20-100 | 46 | 0.9991 | 0.017 | 78.3 % | 0.849 |
| > 100 | 5 | 0.9910 | 0.030 | 40 % | 0.854 |

- chi2/dof > 20 with complete windows: only 51 pulses, all in 2.24 GeV H2 (n 13, median 1.0033, rms 0.024) and Ta X17 (n 38, median 0.9978, rms 0.017, 86.8 % within ±2 %). Their def. 2 ratio is 0.848-0.851, against the 0.853 baseline.
- Amplitude dependence: the def. 1 median drifts from 0.993 to 0.998 over 50-3000 ADC.
- Contrast: non-monotonic pulses at chi2/dof 5-20 sit at a median of 1.040 (n 317), so the agreement depends on the monotonic selection.
- Examples (`PA/integral_vs_rawsum_20261117.pdf`, page 6): 20 monotonic pulses at chi2/dof 20-107 with complete windows (only 2.24 GeV H2 and Ta exist) are clean single peaks, A 300-2400, with the fit following the data up to a small peak/tail misfit; ratio 0.97-1.03. Two lead-glass examples (G535, G455, A ≈ 70) show the fitted tail above the data, from the weak LG template statistics.

## Truncated pulses

6,131 pulses (5,954 monotonic) have their 32-sample fit window past the record end. Those with chi2/dof > 20 (0.73 GeV 9, 2.24 GeV H2 36, 3.49 GeV 175, Ta 622) have a ratio of 1.25-1.28 (rms 0.09-0.14), because A·I_m extrapolates past the record. Truncated pulses with chi2/dof < 5 have per-run-type medians of 1.00-1.02 (the 2-5 bin alone is at 1.09). The high chi2 there is a record-end effect, not a shape problem.

## Conclusions and decisions

1. For monotonic single peaks, the fitted integral equals the raw charge within ~0.5 % (median) and 1.3 % (rms) irrespective of chi2/dof (chi2/dof < 5 pooled; the 51 complete-window pulses above 20 agree as well; the sparse 5-20 bin is wider, rms 0.037). A high chi2/dof on a clean single peak is a shape detail, and the amplitude-based energy is safe. chi2 still matters for pile-up detection: the same mismatch produces fake secondaries, handled under T3a of `notes/plan_pileup_fit_validation_20261116.md`.
2. Late pulses truncated by the record end need no special treatment. They are not the trigger pulse, and the time gate (in-time window 140-156 ns) removes them from the energy measurement. A pulse starting ≳ 200 ns after the in-time pulse does not overlap its ~130 ns template window.
3. Pile-up rates: "pile-up" is defined as pulses whose onset lies within about [t_trig − 40 ns, t_trig + 130 ns] of the in-time pulse (overlapping column). All other pulses form an occupancy column, useful as rate × overlap window ≈ random pile-up probability. To be implemented at T6 of the plan.

Item 16 closed 2026-11-17.

## Files

- Script: `PA/integral_vs_rawsum.py` (reads the cache only; no new decode pass).
- Report: `PA/integral_vs_rawsum_20261117.pdf` (plots, table pages, examples including the truncated class).
- Tables: `PA/integral_vs_rawsum_20261117.tsv` (by class, run type, chi2/dof bin); `PA/integral_vs_rawsum_20261117_waveforms.tsv` (per waveform, 178,253 rows).
