# HyCal energy resolution: where we stand and how to improve it

**Date:** 2026-11-12
**Branch:** `waveform_jz` (not committed)
**Scope:** PbWO4 first; lead glass where the material allows. This note covers what `waveform_analysis/`, `crosstalk_study/`, `time_study/`, `mf_pulse.py` and `notes/` (10-21 to 11-08) have measured, and how other calorimeter experiments calibrate energy. It ranks the options for further improvement and is the evidence base for a later plan.
**Paths** are relative to `analysis/pyscripts/` unless stated.

**Tags.** **[F]** = read from a file or a cited paper (source given). **[I]** = my inference or arithmetic. **[L]** = an option that comes from the literature survey (section 5).

**Abbreviations.**
- CB: Crystal Ball core sigma/E. All widths are in %, differences in percentage points (pp).
- Pool: the benchmark's region. Seed radius < 200 mm, seed outside rows 15-20 x columns 15-20, modules with at least 30 elastic events.
- oos: out-of-sample, an even/odd split cross fit.
- R5nr-ri: the old-chain reference (11-05 note).
- d_pu: CB without pile-up-flagged clusters minus CB with all clusters.

## 0. Summary

**Where we stand.**
- **The pool core is close to the beam-test level of comparable crystals.** The template-fit benchmark gives a CB core sigma of 3.39 / 1.88 / 1.59 % at 0.73 / 2.24 / 3.49 GeV (legacy, gains + position correction). The PrimEx HyCal three-term fit predicts 3.36 / 1.80 / 1.55 %, and the GlueX CCAL fit predicts 3.45 / 1.90 / 1.54 % (section 1).
  - Those references are Gaussian or CB fits; on our Gaussian widths the excess is larger, 1.95 / 0.71 / 0.54 % (section 1).
- **The pool-core gains left are small.** In the pool core every lever left is worth 0.01-0.08 pp. Most of the remaining physics value lies outside the pool, where nothing is measured: outer PbWO4, the ring around the beam hole, the transition, and the lead glass. About 870 modules keep a prior gain, and every lead-glass calibration factor is 0.

**Ranked options.** Gains are on CB sigma/E. "Pool" means the current benchmark number; "coverage" means regions that have no number yet.

| # | Option | Budget term attacked | Expected gain | Cost | First test (data and scripts on disk) |
|---|---|---|---|---|---|
| 1 | Reconcile the new benchmark with the old chain at 2.24 GeV | Møller rejection, gain solver, peak fitter | up to -0.08 pp at 2.24 GeV (1.876 vs 1.794 on the same 42 splits) | low | `energy_resolution_benchmark.py report` on `output/energy_resolution/scan/2p24GeV_H2`, changing one item at a time: `MOLLER_REJECT`, the gain-window refit, the peak fitter |
| 2 | Region-resolved resolution, plus a calibration that reaches beyond the pool (ep + Møller, global solve) [L] | calibration constant term outside the pool | unknown; the resolution in those regions is unmeasured | medium | add per-ring output to `report`, then a ring-wise E/E' and E/E_ee spectrum from the existing scans |
| 3 | Secondary pulses in members: 24-34 % of pool clusters carry one (whether they are spurious is [I]) | amplitude taken by fitted secondaries | ≤ 0.04-0.05 pp in H2 legacy, ≤ 0.03-0.04 after the gate; Ta ≤ 0.11 after the gate | low-medium | dt and A2/A1 breakdown of flagged secondaries from the scan npz; then MIN_SEP 12 vs 8 ns on run 025320 (the open 11-08 to-do) |
| 4 | Better gain solve: seed-indexed refinement, or a global least-squares fit with split-sample validation [L] | residual per-module scale | -0.02 to -0.04 pp e-p; -0.07 to -0.09 pp Ta (old-basis measurement) | low-medium | port the refinement into `report` (no switch exists), judge it on oos |
| 5 | Edge-aware position correction: ring-dependent maps, a width term, a shower-shape regression, log-weight retune [L] | position and edge width | upper bound about 0.1 pp at 3.49 GeV [I]; realistic gain unknown | medium | edge-vs-centre sigma per run type from the existing pool; log_weight scan |
| 6 | Lead-glass calibration and non-linearity [L] | LG has no calibration at all | from "no number" to the about 6 %/sqrt(E) reference | high | LG E/E' spectra at 0.73 / 2.24 GeV from the existing scans (WaveAnalyzer hits) with a flat gain |
| 7 | Ta X17: per-run calibration, refinement, gate already adopted | pile-up and gain epochs | refinement -0.07 to -0.09 pp; the gate already gives -0.27 pp | low | report rerun on `2p24GeV_Ta_X17` |
| 8 | Port the onset-phase integral correction I(phi) to production | per-module phase bias | about -0.04 pp [I], untested | low | report with I(phi) off vs on |
| 9 | Per-module non-linearity from Møller vs ep | energy scale vs energy (Møller) | about 0 at the ep peak; matters for Møller and for the scale | medium | `NONLINEAR` arm; Møller-vs-ep peak table per ring |
| 10 | Common-mode noise subtraction, fixed-pedestal fit [L] | electronic noise | ≤ 0.01 pp at ≥ 2.24 GeV, ≤ 0.05 pp at 0.73 GeV [I] | medium | noise-only matched-filter fit on baseline windows, converted to MeV |
| 11 | Template refinements (amplitude-indexed shapes, tail-ripple fix) | shape scatter | ≤ 0.01-0.03 pp | medium | none needed first; low priority |

**Dropped as resolution levers** (section 4.3): crosstalk correction (< 0.03 pp; keep only as a veto for fake low-energy clusters), saturation (0.0003-0.014 % of clusters), within-run LMS drift correction (measured null), integral vs matched filter (already won), tighter time gates, chi2 cuts, leakage and Geant4-profile corrections.

**Next three steps** (section 6):
1. Reconcile the 2.24 GeV gap with report-only reruns.
2. Add region-resolved output and extend the calibration beyond the pool.
3. Test the secondary-pulse hypothesis: a breakdown of the flagged secondaries first, then MIN_SEP 12 vs 8 ns on run 025320.

## 1. Current state

**Template-fit benchmark** [F] (`waveform_analysis/energy_resolution/<RT>/summary.md`, default arm = gains + position correction, Møller rejection and non-linearity off: `energy_resolution_benchmark.py:153,155`):

| run type | legacy CB (in / oos) | best gate, oos CB | Gauss ±2σ | pile-up flag | d_pu (legacy) | arm "none" (no gains, no position correction) |
|---|---|---|---|---|---|---|
| 0.73 GeV H2 | 3.389 ± 0.014 / 3.420 | 3.404 (11.5 ns) | 3.884 | 9.8 % | -0.047 | 4.103 |
| 2.24 GeV H2 | 1.876 ± 0.005 / 1.889 | 1.871 (13 ns) | 1.935 | 24.3 % | -0.043 | 2.004 |
| 3.49 GeV H2 (new runs) | 1.585 ± 0.006 / 1.614 | 1.601 (13 ns) | 1.637 | 33.7 % | -0.054 | 2.184 |
| 3.49 GeV H2 old runs | 1.559 ± 0.006 / 1.579 | 1.568 | 1.616 | 25.4 % | -0.035 | 2.810 |
| 2.24 GeV Ta X17 | 2.308 ± 0.035 / 2.332 | 2.061 (7.5 ns) | 2.285 | 36.4 % | -0.270 | 2.597 |

- **Old-chain references** [F]. R5nr-ri (`notes/PbWO4_energy_resolution_20261105.md` §4) gives 3.417 / 1.794 / 1.591 % (0.73 / 2.24 / 3.49 GeV) and 1.800 % on C. It is the legacy path, with Møller rejection and Møller-rejected gains, gain refinement off, seed 30 MeV and log_weight 5.2 (job line `Commands.md`:149-156, with `reconstruction_config_r5.json` via `--hycal-cfg`). Its result tables are no longer on disk.
- **Refined R5** [F] (old basis, `notes/PbWO4_energy_resolution_20261029.md` §3.1): 1.472 / 1.773 % at 3.49 / 2.24 GeV.
- **Comparability.** The 3.49 GeV references used runs 24386/24503/24512, so they are not directly comparable with the new runs (`waveform_analysis/opens.md` 3c).
- **The 3-term fit is unconstrained** [F] (`energy_resolution/resolution_fit.tsv`). sqrt(a² + b²/E + c²/E²) gives a = 1.10, b = 1.99 % sqrt(GeV), c = 1.62 % GeV. It has three points and three parameters, so 0 dof. The 2-term form gives a = 0.36, b = 2.82 with chi2 212 for 1 dof, which rejects it.
  - [I] c = 1.62 % GeV means 16 MeV per cluster. The electronic noise estimate is 1.4-4 MeV per 3 x 3 cluster (section 2), so c is not electronic noise. It absorbs the 0.73 GeV specifics: position correction off there, Møller overlap, and the CB/Gaussian instability. A fourth energy point is needed; no 1.1 or 2.1 GeV runs are on disk.
- **References** [F, section 5]. The PRad-II proposal (arXiv:2009.10510, HyCal section) gives 2.6 %/sqrt(E), from a Gaussian fit of a 6 x 6 array in a tagged-photon beam. That is 3.05 / 1.74 / 1.39 % at our energies, stochastic only [I]. The PrimEx fit (X17 proposal PR12-21-003, arXiv:2108.13276 Fig. 14, which shows PrimEx data; term assignment inferred) is 1.13 ⊕ 1.73/sqrt(E) ⊕ 1.77/E, i.e. 3.36 / 1.80 / 1.55 %. The GlueX CCAL (PbWO4, 3 x 3) is 0.53 ⊕ 2.63/sqrt(E) ⊕ 1.07/E, i.e. 3.45 / 1.90 / 1.54 % (draft of NIM A 1013 (2021) 165683; coefficients unverified against the published version).
  - [I] The pool's CB core sits 0.35-0.53 % (in quadrature) above the PrimEx curve, and the Gaussian ±2σ widths 1.95 / 0.71 / 0.54 %. The old chain's 1.794 % at 2.24 GeV is on the curve. The references used Gaussian fits (2.6 %/sqrt(E): "Gaussian fit of the line-shape", arXiv:2009.10510) or CB fits (CCAL); our CB core is the narrower estimator.
  - [I] These references are beam tests or in-beam tests with tagged photons, while ours is in-beam elastic E/E'. Treat them as targets, not as floors.
- **Lead glass** [F]. All 576 G modules have calibration factor 0 in every calibration file (`../../database/calibration/*.json`, checked 11-12: G nonzero = 0 in each file). No lead-glass resolution exists. The published reference is 6.1-6.2 %/sqrt(E) (PRad, arXiv:2302.13818 §3).

## 2. Resolution budget at the ep peak (pool, legacy, gains + position)

No budget in the notes is current. The R4 2-term fit (a 2.48, b 0.70, `notes/PbWO4_energy_resolution_summary.md` §3.6) and the R2 3-term fit predate the runinfo basis. The R5nr-ri 2-term fit is unphysical (b² < 0, chi2 678). The table below is assembled from the measured pieces.

| term | 3.49 GeV (1.585 %) | 2.24 GeV (1.876 %) | status | source |
|---|---|---|---|---|
| photostatistics + shower (stochastic) | 1.39 % (if 2.6 %/sqrt(E)) | 1.74 % | literature, not measured here | arXiv:2009.10510 |
| excess over 2.6 %/sqrt(E), in quadrature | 0.76 % | 0.71 % | [I] arithmetic | — |
| pile-up flag / secondary pulses (d_pu as quadrature equivalent) | ≤ 0.41 % | ≤ 0.40 % | measured as a removal difference (an upper bound, includes selection) | `<RT>/summary.md` |
| gain refinement residual (refinement-off cost +0.021 / +0.027 pp) | about 0.25 % | about 0.32 % | measured on the old basis | `20261029.md` §3.1 [I conversion] |
| template shape scatter | about 0.15 % | about 0.2 % | [I] from the synthetic 0.30 % amplitude width at 3000 ADC | `pileup_analysis/efficiency.tsv` |
| crosstalk (members late by > 8 ns) | 0.08 % (rms 2.8 MeV) | 0.10 % (rms 2.2 MeV) | [I] unsaved diagnostic, 11-12 (see note below) | `output/crosstalk_study/scan/run_{24698,25228}.npz` |
| electronic noise, 3 x 3 | 0.04 % (uncorrelated) to 0.12 % (fully correlated) | 0.06-0.18 % | [I] from raw baseline std 1.0 ADC and template norm | `output/noise_diag_20261028/report.md` |
| position (after correction) | profile flat within ±0.1 %; edge-bin sigma 1.85 vs 1.55 % centre | 2.08 vs 1.81 % | measured (old basis) | `PbWO4_energy_resolution_summary.md` §3.5 |
| clipping | 0.0076 % of clusters | 0.014 % | measured | `pileup_rates_out/summary.txt` |
| within-run gain drift (LMS) | ≤ 0.07 % common | — | measured, correction null | `output/lms_gain_20261028/summary.md` |
| unattributed remainder | about 0.55 % | about 0.4 % | [I] excess ⊖ quadrature sum of the rows above (0.54-0.56 % at 3.49, 0.39-0.43 % at 2.24) | — |

**Crosstalk diagnostic.** Clusters were built with the `crosstalk_study.py` report clustering and gains. For each cluster in the elastic energy slice (1800-2600 MeV on run 25228, 2900-3800 MeV on run 24698, 550-900 MeV on run 24590; no e-p selection), the energy of non-seed members with t_cfd - t_mf > 8 ns was summed and its rms taken over all clusters. To regenerate: run the report clustering on the scan npz above and apply that cut.

[I] Reading of the table:
- The remainder is the candidate home for four effects: the edge-width excess, calibration of the prior-dominated neighbour modules (42-63 % of the modules listed in `gains.tsv` have prior_frac ≥ 0.5), longitudinal non-uniformity, and leakage.
- None of these is measured as a sigma term.
- Against the PrimEx curve (section 1), the pile-up-flag and refinement terms alone would cover the excess.

## 3. Inventory of what was done

| study | result (CB unless stated) | outcome |
|---|---|---|
| Integral (WaveAnalyzer rawsum) → matched filter (MF) amplitude, 10-22 | -0.061 / -0.040 / -0.073 (C) / -0.306 pp at 3.49 / 2.24 / C / 0.73 GeV; Ta -0.11 / -0.20 / 0 (`20261022.md` §4.1) | adopted; predates gains and position |
| Per-module spline templates + refitted gains + 3-band f(A), 10-25/26 | -0.114 / -0.119 / -0.074 / -0.108 pp; Ta -0.082 (`summary.md` S:108) | adopted; largest single step |
| Per-crystal gains (ridge solve) | -0.048 to -0.091 pp at 3.49 / 2.24 / C; free-module sigma_g 0.17-0.34 % (`gains.tsv`) | adopted |
| Gain refinement (seed-indexed) | on is better by 0.021-0.042 pp e-p, 0.07-0.09 pp Ta (`20261029.md` §3.1) | switched off by user decision; no rerun on the new basis |
| Møller-rejected gains | +0.008 ± 0.005 pp vs R4 (`20261029.md` §6). Unmeasured in the new chain (switch off) | kept in the old chain; off in the new benchmark |
| Position correction, new chain (legacy, in / oos, `poscorr_comparison.tsv`) | hybrid1000 vs none: -0.146 / -0.137 pp at 3.49, -0.048 / -0.046 at 2.24, +0.114 / +0.106 at 0.73, Ta -0.046 / -0.062; vs map12 at 3.49: -0.090 / -0.075 | adopted (hybrid1000 at 2.24 / 3.49 and Ta, none at 0.73) |
| even6 / even4 forms | old chain: per-module even6 +0.042 pp vs poly8 at 3.49 (`poscorr_form_choice_20261106.md`); new chain 11-12: global even4 / even6 +0.106 / +0.091 pp vs hybrid1000 at 3.49 (in); per-module even4 / even6 +0.054 / +0.044 pp vs hybrid1000 at 3.49 (in) | rejected |
| Best-50 / one-pulse / f(A) off / p1 shape cut | ≤ 0.02 pp, mostly ≤ 0.004 (`20261027.md` §2.2, `20261105.md` §3-4) | null |
| R7ri energy-dependent W0 | -0.033 / -0.012 / -0.016 pp (`20261105.md` §4) | not adopted |
| Time gate | -0.011 to -0.020 pp e-p; Ta -0.26 to -0.27 pp; < 6-8 ns loses yield (`<RT>/summary.tsv`) | adopt 12-16 ns e-p, 8-12 ns Ta |
| Pile-up rates (`pileup_rates.py`) | channel 0.28-0.96 %, cluster ≤ 1.46 % in H2; Ta 7.3 % channel (`pileup_rates_out/summary.txt`) | small in H2 |
| Synthetic two-pulse study; MIN_SEP 12 → 8 ns | ≥ 97.5 % efficiency at dt ≥ 12 ns, r ≥ 0.1; 0.6-3.4 % spurious at A = 3000 with smearing (`20261108.md`) | uncommitted; real-data check open |
| Noise diagnosis | raw std 1.0 ADC; 27 % tone, 33 % common, 40 % independent; subtraction gives no gain in sigma_t or chi2 (`output/noise_diag_20261028/report.md`) | not tested on energy |
| LMS drift | ≤ 0.07 % within a run; correction +0.000 pp; Ta crate-3 jumps between runs, up to tens of percent (`output/lms_gain_20261028/summary.md`) | null within a run |
| Non-linearity | only the June7 file has non-zero nl1/nl2 (1006 modules); `non_linear` is 0 in every file; round1_* and x17_* have nl1 = nl2 = 0; Møller vs ep gives \|NL\| ≲ 0.5 % over 0.5-2.5 GeV (`output/nonlin_20261104/step3`) | no sigma on/off number on disk |
| Crosstalk (`crosstalk_study/`, 10-30 and 11-07 notes) | electronic; adjacent fADC channel and fixed cross-crate pairs; k ≈ 1.0-1.4 %, delay 12-17 ns; rms 2-3 MeV per elastic cluster ([I], unsaved diagnostic, section 2) | negligible for sigma/E; veto untested |
| Clustering knobs | min_module_energy null; seed 30 ≤ 0.016 pp; log_weight 5.2 costs 0.013-0.041 pp; Geant4 profile ≤ 0.004; leakage inert; C++ energy_bias (R6) +0.026 to +0.051 (`summary.md`, `20261030.md` §8) | mostly null |
| Pulse shape / saturation | bumpy 0.09 %, no hot modules; clipped 0.0003-0.014 % of clusters (`report/pulse_shape_study.md`, `pileup_rates_out/summary.txt`) | negligible |
| Timing (`time_study/`) | seed-vs-veto sigma 0.35-0.42 ns; member-vs-seed 0.64-0.84 ns calibrated (`time_study/tables/pooled_vetos_allsplits.md`) | supports the gate; no energy lever |

**Unfinished** [F]:
- the 0.73 GeV CB-vs-Gaussian headline and the unattributed +0.273 pp (`20261030.md` §4);
- P3-P5 of the 11-05 plan (the 0.73 GeV Møller block, energy_bias and the Geant4 profile on the ri gains, NL steps 2-3);
- the MIN_SEP real-data check;
- the noise-model refit for in-time bad fits (chi2/dof > 3 for 33.9 % of 20-50 ADC pulses);
- the cross-crate crosstalk path;
- W335 at 0.73 GeV;
- the 2.6 % member_e mismatch;
- the benchmark has no clip count or chi2 exclusion (`opens.md` 6);
- `opens.md` items 20-26.

## 4. Options in detail

### 4.1 Pool core (current benchmark number)

**Option 1. Reconcile with the old chain at 2.24 GeV.**
- **Evidence** [F]: The new benchmark gives 1.876 %; R5nr-ri gives 1.794 % on the same 42 splits.
  - With no gains and no position correction, both chains agree on run 25320: 2.015 vs 2.004 % (`opens.md` 3a).
  - The two `opens.md` entries conflict on the Møller-rejection share of the gap. Entry 3a gives it as -0.01 pp; entry 3e attributes the remaining gap to Møller rejection. In the new chain it is unmeasured: the switch is off, so the "+ Møller rejection" arm equals the default (`2p24GeV_H2/summary.md`).
  - [I] Gains + position are worth 0.13 pp in the new chain against about 0.22 pp in the old one.
- **Candidate differences.** These come from the R5nr-ri job line (`Commands.md`:149-156, with `reconstruction_config_r5.json` via `--hycal-cfg`) compared with the benchmark constants (`energy_resolution_benchmark.py`:95-160). Seed 30 MeV, log_weight 5.2, the r < 200 mm pool with rows/columns 15-20 excluded, and the legacy path are already the same, and NL is inert at 2.24 GeV.
  - Møller rejection: `--moller-reject kin` plus `--moller-reject-gains on` in the old chain; `MOLLER_REJECT = False` in the new one.
  - Gain solver: `percrystal_calib.py fit --fixed-windows` in the old chain; the new solve refits the windows in its first passes (`GAIN_N_REFIT = 3`) and has its own prior strength (`GAIN_SEED_FRAC`).
  - Peak fitter: `--fit-select chi2 --cb-two-peak` in the old chain; the new two-peak CB uses a Møller tie, its own fit range and its own binning (`FIT_HALF`, `BIN`, `MOLLER_TIE_*`).
  - Pulse model and bank: they differ (`--mf-scale`, best-50 bank vs `pulse_templates/`), but the chains agree without corrections on 25320, so this item comes last.
- **Gain.** Up to 0.08 pp (4 %) at 2.24 GeV. At 3.49 GeV the comparison is blocked by the different runs.
- **Risk.** The reference tables are gone (`output/r5ri_20261104` is not on disk), so the 1.794 % cannot be re-read; it comes from the note.
- **First test.** The report takes minutes per run type (`energy_resolution/README.md` "Scanning the gate"). Rerun it on `output/energy_resolution/scan/2p24GeV_H2`, changing one item above at a time (the pulse model last); the switches are code constants, so `executor` makes the changes. Compare the paired in-sample CB with 1.794 % and report oos alongside.

**Option 3. Secondary pulses in cluster members (spurious or small real ones).**
- **Evidence** [F]:
  - 24.3 % / 33.7 % of pool clusters at 2.24 / 3.49 GeV have a member with two or more in-window pulses (`<RT>/summary.md`). This flag has no ratio or energy threshold.
  - The pile-up cut rate is 0.42-0.83 % per cluster (`pileup_rates_out/summary.txt`). It requires r > 0.1, |dt| ≥ 8 ns and E_main ≥ 0.1 E_beam, so it is a different quantity.
  - The raw channel k ≥ 2 rates (1.85 / 2.86 / 3.21 % at 0.73 / 2.24 / 3.49 GeV, same file), compounded over about 6-13 members, reproduce the flag fractions.
  - [I] So the flags are mostly small secondaries; whether they are spurious is not established.
  - The H2 extras sit 10-14 ns after the main pulse at r ≈ 0.12, which `opens.md`:9 reads as a shape-fluctuation tail.
  - An accepted secondary takes a median 6.3-7.7 % from the dominant pulse (`summary.md` S:348).
  - The synthetic study shows 0.6-3.4 % spurious seconds at A = 3000 with template smearing, at MIN_SEP 8 ns (`20261108.md`).
  - The benchmark scans already use 8 ns (worker check), and no 12 vs 8 ns real-data comparison exists.
- **Gain.** Bounded by d_pu: 0.043-0.054 pp in H2 legacy, 0.025-0.041 pp at g12. Ta: 0.27 pp legacy, 0.107 pp after the gate.
- **Risks.** A stricter secondary acceptance lowers the true pile-up efficiency, which matters for Ta. Summing all in-time pulses (the gated path) already recovers part of the loss.
- **First test.**
  1. Histogram the flagged secondaries' dt and A2/A1 from the existing scan npz (`n_pulses` and the extra flag are stored), and compute d_pu per bin.
  2. Rescan run 025320 at MIN_SEP 12 and at 8 ns. This takes several minutes, so it is a planned run, not a diagnostic. Compare the flag fraction, d_pu and CB.

**Option 4. Gain solve** [L: NPS global linear least-squares solve, CMS and GlueX split-sample validation].
- **Evidence** [F]:
  - Refinement off costs 0.021-0.042 pp e-p and 0.07-0.09 pp Ta (old basis).
  - The per-module absolute-constant spread is 2.3-2.5 %, 7-12 times the fit uncertainty, and the gains undo only about 30 % of the seed-file spread (`output/nonlin_20261104/step2/findings_step2.md`).
  - The in-sample/oos gap is 0.01-0.03 pp (splits.tsv), the same size as the expected gain.
- **Gain.** -0.02 to -0.04 pp e-p, -0.07 to -0.09 pp Ta.
- **Risks.** Overfitting (judge it only on oos), and the user's 10-29 decision to run with refinement off.
- **First test.** The benchmark has no refinement switch (no `refine` in `energy_resolution_benchmark.py`). Port the old chain's seed-indexed step into `report` (cost low-medium), then evaluate it oos.

**Option 5. Edge and position width.**
- **Evidence** [F]:
  - Edge-bin sigma/mu is 1.85 against 1.55 % at the centre (3.49 GeV) and 2.08 against 1.81 % (2.24 GeV); a scale map cannot remove this.
  - Ring and side dependence is 0.15-0.29 % RMS, with single bins up to 1.66 % (`PbWO4_energy_resolution_summary.md`:193, from the 10-25 bilinear-era maps).
  - log_weight 5.2 costs 0.013-0.041 pp (`20261029.md` §5), cause unresolved.
- **Literature** [L]:
  - CMS corrects with a BDT on shower shape and local crack coordinates. A centre-to-edge residual of about 1 % remains, worth 0.3-0.5 % RMS (arXiv:1306.2016 §4.4-4.5).
  - HPS uses an MC correction parametrised against the distance to the edge (arXiv:1610.04319 §3).
- **Gain** [I]. Upper bound about 0.1 pp at 3.49 GeV if the edge widths reached the centre width (two-population mixture with about 40 % edge clusters). The realistic gain is unknown, since part of the edge excess is physics (gaps, wrapping).
- **First test.** Edge-vs-centre CB per run type from the existing pool on the new basis, then a log_weight scan in report.

**Option 8. Phase correction I(phi).**
- **Evidence** [F]: 1.8 % p-p per module (`energy_resolution/README.md`). The benchmark corrects it; production `mf_pulse` uses the phase average (`ep_resolution_walkthrough_20261102.md`:63).
- **Gain** [I]: about 0.35 % per cluster, or about 0.04 pp. Untested.
- **First test.** Report with I(phi) off.

**Options 10-11. Noise and templates.**
- **Noise** [F]: 60 % of the baseline variance is common within a board. Subtracting it was tested only on time and chi2 (section 3).
- **Noise per cluster** [I]: 1.4 MeV (uncorrelated) to 4 MeV (correlated), so ≤ 0.01 pp at ≥ 2.24 GeV and ≤ 0.05 pp at 0.73 GeV.
- **Constraint** [F]: the readout is zero-suppressed (7-20 PbWO4 channels per event, `summary.md` S:150). That limits common-mode estimates to the channels read out.
- **Fixed pedestal** [L]: CMS multifit uses a fixed pedestal for better noise suppression, but a 1 ADC drift biases it by 0.3 ADC (arXiv:2006.14359 §5.3).
- **Templates**: shape scatter is about 0.15 % [I], and every variant tried moved CB by ≤ 0.03 pp.

### 4.2 Coverage outside the pool and other targets

**Option 2. Region-resolved resolution and calibration beyond the pool.**
- **Evidence** [F]:
  - The benchmark reports the pool only (`energy_resolution_benchmark.py:109-113`).
  - The gain solve frees only 267-284 modules, all at 53 < r < 200 mm. About 870 modules keep the prior (`nonlin_20261104/step2`). 257 modules in the 25320 table have no seed events.
  - Outside the pool, E/E' drifts by +0.5-0.7 % (step 3).
  - Rings at r < 85 mm show +0.4 / +1.4 % shifts.
- **Literature** [L]:
  - PRad-I calibrated in situ with ep and double-arm Møller, using about 400 M events (Weizhi 2016 slides, wiki.jlab.org pcrewiki).
  - NPS solves all 1080 gains as one linear system with a low-side radiative-tail cut (NPS 2025 note, Sec 3.5).
  - The 2023 HyCal calibration proposal plans cosmics + pi0 + Compton (jlab.org/primex HyCal_calibration_proposal.pdf).
- **Gain.** Unknown; it is the first number to measure.
- **Risks** [I]:
  - Møller-based gains need the E(θ) relation and a beam-offset correction: 1 mm gives 0.3-0.9 % in E/E_ee (step 3).
  - ep rates fall steeply with angle.
- **First test.** Per-ring E/E' and E/E_ee peak and width from the existing scans with the current gains. Report the prior fraction per ring.

**Option 6. Lead glass.**
- **Evidence** [F]: no calibration factor in any file. PRad-I measured an LG non-linearity of about -0.27/GeV, about 30 % over the range (APS 2017 slides), and the PMTs are more than 50 years old (PRad-II technical review rev14).
- **Note.** The PRad-II proposal planned to replace the LG with PbWO4 (arXiv:2009.10510), yet the data have 576 G modules. Is the LG in the physics acceptance? (Q1)
- **First test.** LG cluster E/E' at 0.73 / 2.24 GeV from the scans with a flat gain, to size the problem.

**Option 7. Ta X17.**
- The gate already gives -0.26 to -0.27 pp.
- Refinement would give -0.07 to -0.09 pp.
- LMS shows crate-3 slot shifts between runs (26230: slot 15 +110.9 %), so calibrate per run (runinfo already splits the files).
- The scans are statistics-starved.

**Option 9. Non-linearity.**
- **Evidence** [F]: |NL| ≲ 0.5 % over 0.5-2.5 GeV. The ep peak spans < 0.5 % in E', so it cannot constrain per-module NL. No sigma comparison with NL on vs off exists.
- [I] Little effect on the ep-peak sigma. It matters for the Møller energy scale and for the cross-section energy cut.
- **Literature** [L]: ALICE fits a piecewise NL by making the pi0 mass independent of pT (arXiv:1902.06145 §4).

### 4.3 Dropped, with reasons
- **Crosstalk correction**:
  - [I] Late members add rms 2.2 / 2.8 / 0.66 MeV per elastic cluster at 2.24 / 3.49 / 0.73 GeV (unsaved diagnostic, section 2).
  - [F] Seed chi2 does not change with a late member present.
  - The waveform is not a scaled copy of the aggressor (chi2 116-142 against 33 for real leakage).
  - The (t_cfd - t_mf) > 8 ns veto is useful only against fake 10-100 MeV clusters (`crosstalk_study/out/late_component.tsv`, `crosstalk_figures_20261030.md`).
- **Saturation**: 0.0003-0.014 % of clusters are clipped. The template fit already gives clipped samples zero weight.
- **Within-run LMS correction**: a measured +0.000 pp.
- **Integral vs MF**: settled in favour of MF.
- **Tighter gates and chi2 cuts**: they lose yield or energy (`20261030.md` §6).
- **Leakage, Geant4 profile, min_module_energy**: inert.

## 5. How other experiments calibrate calorimeter energy

Sources are labelled as in the URL list at the end. "Unverified" means the claim was not seen on a primary source page.

**PRad-I / PrimEx HyCal.**
- **Calibration** [P1-P3]:
  - a tagged-photon snake scan (PRad 0.25-1.05 GeV, PrimEx 1-5.5 GeV), iterating gain(n+1) = gain(n)/<E_meas/E_exp>;
  - then in situ ep + double-arm Møller, histogramming E_rec/E_exp for the central module of each cluster.
- **LMS and temperature** [P1, P3]: LED through an integrating sphere and fibres, with a PIN and an alpha-source reference PMT; it cannot correct radiation-damaged modules. Light yield changes by -2 %/°C, so the calorimeter was held at 14 ± 0.1 °C.
- **Non-linearity** [P2]: per module, <E_cl/E_γ> = α E + β. About 0.01/GeV for PbWO4; about -0.27/GeV for lead glass (preliminary).
- **Resolution**: 2.6 %/sqrt(E) PbWO4 [P5]; 2.4 %/sqrt(E) PbWO4 and 6.2 %/sqrt(E) LG [arXiv:2302.13818].
- **PRad-II plan** [P6]: cosmics (20-30 %), then pi0 (2.5-5 % per channel), then Compton.

**NPS (Hall C, PbWO4 + fADC250)** [N1-N4]:
- Amplitude: per-channel reference shapes from clean elastic pulses, then a full-waveform multi-pulse fit with a free baseline.
- Gains: one linear least-squares system for all 1080 blocks from H(e,e'p), five cycles, with cuts on the position residual, edge seeds and the low-side radiative tail.
- Threshold loss is absorbed by 1/(1+α(E-<E>)) with α = 0.05.
- Resolution: 1.2-1.3 % at 4.5-7.3 GeV (CB), degrading with irradiation.

**HPS ECal (PbWO4 + APD, fADC250)** [H1]:
- Calibration chain: cosmics (about 18 MeV per crystal), then elastic or full-energy-electron iteration with seed fraction > 60 %, two iterations; wide-angle bremsstrahlung for mid energies.
- LED stable to < 1 % over 3 days; the motherboards were redesigned against crosstalk.
- MC correction f = A/E + B/sqrt(E) + C as a function of the distance to the edge; log-weight w0 = 3.1.
- Resolution: 1.62/E ⊕ 2.87/sqrt(E) ⊕ 2.5 % (fiducial); the stochastic term worsens near edges (Fig. 12).

**CLAS12 FT-Cal (PbWO4 + APD)** [F1]:
- Calibration: cosmics, then a per-crystal ep seed-energy edge fit against Geant4 (seed fraction > 55 %), then pi0 iterated to ±0.5 MeV of the mass.
- Leakage is corrected per seed crystal with a 4th-order polynomial in E; leakage is 16 % at 0.5 GeV and 6.6 % at 4.5 GeV for a 10 MeV threshold.
- Resolution: 3.3 % at 2.2 GeV. It is about 1 % worse than design, blamed on the calibration spread and on thresholds.

**GlueX FCAL / CCAL / ECAL** [G1-G4, C2]:
- **FCAL amplitude**: fADC firmware integral over a fixed window chosen for signal-to-noise; run-averaged pedestals.
- **FCAL calibration**: iterative pi0-width minimisation, with a split-sample check that revealed fitting of background fluctuations, an LED-smoothness prior, and empirical NL from the pi0-mass energy dependence.
- **CCAL** (140 PbWO4): tagged-beam HV equalisation, then beam-energy and Compton gain refinement; LED with two reference PMTs (one with an Am-241 YAP pulser), stable to 1 % over 48 days.
  - An active-base amplifier caused percent-level non-linearity and about 10 % worse resolution. Resolution 2.63 %/sqrt(E) ⊕ 1.07 %/E ⊕ 0.53 % (draft).
- **ECAL** (1596 PbWO4): LMS HV slope, then cosmic MIP, then Compton; about 1.7 % at 6-7 GeV.

**CMS ECAL (PbWO4)** [M1-M3]:
- **Amplitude**:
  - Run 1 used a weights method (an optimal filter) with a per-event pedestal.
  - Multifit is an NNLS template fit with up to 10 amplitudes, a fixed pedestal and a measured pulse-shape covariance.
  - Multifit gains 10-15 % on the quadrature resolution at pT < 5 GeV.
- **Intercalibration**: phi-symmetry, pi0/eta, E/p and Z→ee, combined by precision. The precisions themselves come from the pairwise spread of three methods, with an odd/even statistical check. The result is 0.4 % in the central barrel.
- **Monitoring**: laser + PN diodes, correction (R0/R)^α, temperature held to 0.03 °C.
- **Cluster correction**: BDT regression with crack coordinates.
- **Resolution**: test beam 2.8 %/sqrt(E) ⊕ 12 %/E ⊕ 0.3 %. The constant term is driven by longitudinal uniformity, rear leakage and the intercalibration spread.

**ALICE PHOS (PbWO4, −25 °C)** [A1]:
- LED HV pre-equalisation, then an iterative pi0 update c = (m_pi0/m_i)^n. A toy study found n = 1.6 best (n = 2 oscillates), needing 2-3 iterations, with the correction applied to the central cell only.
- Run-by-run pi0 mass scale correction.
- Piecewise NL fixed by a pT-independent pi0 mass; E/p checks the scale to 0.2 %.
- Test beam: 1.3 %/E ⊕ 3.58 %/sqrt(E) ⊕ 1.12 %.

**PANDA EMC (PWO-II)** [D1]:
- Calibration: cosmics (about 10 %), then pi0/eta iteration (Crystal Barrel needed 1.2 M events for 1380 crystals), then a light pulser, with 0.1 °C stability.
- Waveform feature extraction by moving-window deconvolution and digital CFD (KVI slides).

**CsI, BGO, LKr and lead-glass families.** BaBar and L3 quote a/E^(1/4) ⊕ b, so their "a" is not comparable with ours.

- **NA48 LKr** [K1]:
  - Amplitude: a 3-sample digital filter binned in pulse time, accurate to 0.1 %.
  - Calibration: Ke3 E/p iteration using only 25-40 GeV electrons, to keep non-linearity out of the gains. It lowered the cell RMS from 0.41 % to 0.15 % and the constant term from about 0.6 % to 0.42 %.
  - Resolution: 3.2 %/sqrt(E) ⊕ 0.09/E ⊕ 0.42 %.
  - It publishes a constant-term budget: intercalibration 0.15 %, gaps 0.1-0.2 %, pulse reconstruction 0.1 %, calibration vs physics pulse shape 0.15 %, MC about 0.2 %. This is the model for section 2.
- **KTeV CsI** (PMT) [K2, K3]: laser drift correction before clustering; about 600 M Ke3 events give per-channel scales of 0.03 % and per-channel non-linearity corrections (0.5 % on average).
  - Cluster corrections: a 5 x 5 sub-cell transverse map, a partial-cluster correction for crystals below threshold, and out-of-cone and beam-hole leakage. Resolution 2 %/sqrt(E) ⊕ 0.4 %.
- **BaBar CsI(Tl)** [K4, K5]: a digital filter from the measured pulse and noise spectrum (equivalent noise 440 → 230 keV).
  - Bhabha events solved as one linear system for all crystals (about 200 e per crystal for 0.35 %). Seed 10 MeV and neighbour 1 MeV thresholds degrade the shower-edge resolution; 2-4 % crosstalk gave a residual non-linearity.
- **Belle / Belle II ECL** [K6, K7]: an FPGA chi2 template fit with a measured noise covariance (photon + hadron components); pile-up noise falls about 1.5-fold and a time cut rejects fake clusters 7-fold.
  - Calibration from the e+e- → γγ upper edge of the highest-energy crystal plus Bhabha; Belle II details come from slides and software documentation only.
- **COMPASS ECAL1/2 (lead glass)** [K8]: electron-beam calibration, then iterative pi0 with α = 1/(1+ΔM/M)², 8-10 iterations (pi0 sigma 7.6 → 4.6 MeV); a laser/LED run-by-run factor L_i(t)/L_i(0).
- **Other families** [K9-K13]:
  - PHENIX PbSc/PbGl: MIP calibration (about 3 % per tower) with an electron-beam scale and laser drift.
  - KOTO: cosmics, then a 3π0 width minimisation.
  - Mu2e: multi-pulse waveform fits beat window integration.
  - GRAAL BGO: per-crystal LED-attenuator linearity curves improved the eta width by about 20 %.
  - Crystal Ball/TAPS: iterative pi0; it left its PbWO4 out of analysis because its "energy resolution is not understood well".
- **Unverified (primary text not reached)**: the PRad Nature 2019 Methods (doi:10.1038/s41586-019-1721-2); the L3 BGO Xenon and RFQ texts (abstracts only); KOTO and Mu2e fit details; the PHENIX PbGl resolution.

**Optimal filtering** [O1]: Cleland and Stern weights minimise the noise variance of an unbiased amplitude, given the pulse shape and the noise autocorrelation (σ_N ≈ σ_n/‖h‖). A wrong pulse shape gives a systematic error and a raised quality factor. A generic noise-reduction factor against integration does not exist, because it depends on the shape and the noise correlation.

**Comparison table.**

| technique | used by | this analysis now | gap | applicability to PRad-II |
|---|---|---|---|---|
| template / multi-pulse waveform fit | NPS, CMS multifit, Belle II | MF with per-module spline templates and shape covariance | secondaries in 24-34 % of pool clusters (option 3) | direct; tune the secondary acceptance |
| fixed pedestal vs free baseline | CMS (fixed) | free offset clipped ±1σ, clip hit 19-24 % | untested on energy | small gain; needs pedestal stability per run |
| global linear gain solve over all blocks | NPS, PRad-I (ep + Møller) | ridge solve on the pool, ep only | no Møller, prior outside the pool | ep + Møller on disk; option 2 |
| split-sample validation of calibration | GlueX FCAL, CMS | even/odd oos for gains and position | refinement not judged oos | cheap; already built |
| several independent intercalibrations combined | CMS | ep only | Møller and LMS not used as independent constants | ep + Møller + LMS on disk |
| pi0 / Compton / cosmic calibration | FCAL, PHOS, CCAL, ECAL, FT-Cal | none | no pi0 or cosmic data known on disk (Q4) | only if such runs exist |
| LED / reference-PMT gain monitoring | PRad-I, CCAL, CMS, HPS | LMS checked, null within a run | between-run jumps on Ta crate 3 | per-run gains already; LMS for epochs |
| temperature control or correction (-2 %/°C) | all PbWO4 | not studied | unknown | needs slow-control logs (Q5) |
| per-module NL from an energy lever arm | PRad-I snake scan, ALICE pi0, NPS | only the June7 file; \|NL\| ≲ 0.5 % globally | per module unknown | Møller gives 0.5-2.5 GeV |
| MC leakage / edge correction | HPS, FT-Cal, CMS BDT | data-driven poly8 + map | width excess at the edge; no MC | Geant4 profile exists; regression needs MC |
| seed-fraction and radiative-tail cuts in the calibration sample | FT-Cal, HPS, NPS | elastic window ±2σ | seed-fraction cut untested | cheap report change |
| calibrate in a narrow energy band, keep NL separate | NA48 (25-40 GeV e), Belle II (single-crystal edge) | ep peak only for gains (good) | Møller mixes energies if added to the gains | use ep for gains, Møller for NL / cross-check |
| below-threshold (partial-cluster) correction | KTeV, BaBar (thresholds), NPS α term | none; zero-suppressed readout | threshold loss unquantified | needs the ZS threshold (Q6) and the Geant4 profile |
| published constant-term budget | NA48, CMS | none current (section 2 here) | budget mostly inferred | needs the runs named in section 4 |

## 6. Recommended sequence

1. **Close the 2.24 GeV gap** (option 1, then option 4).
   - Data: existing scans `output/energy_resolution/scan/2p24GeV_H2` and `3p49GeV_H2`; report only.
   - Change one item at a time (section 4.1, option 1): Møller rejection, the gain-window refit and prior, the peak fitter; then port the refinement (option 4).
   - Success: in-sample CB at 2.24 GeV at or below 1.80 % (oos reported alongside), with every step attributed. A paired gain of more than 2 SE counts as real.
2. **Region-resolved resolution and extended calibration** (option 2).
   - Add per-ring output (hole ring r < 85 mm, 85-120, 120-160, 160-200, 200-300, > 300 mm, transition) to `report`.
   - Measure E/E' (ep) at 2.24 / 3.49 GeV and E/E_ee (Møller) at 0.73 / 2.24 GeV per ring, with the prior fraction per ring.
   - Then add Møller clusters to the gain solve where ep is too sparse, with the NL fixed first (Møller spans several energies), validated on oos.
   - Success: sigma/E per ring for inner, mid and outer rings at ep and Møller energies, and the prior-dominated fraction reduced.
3. **Secondary-pulse test** (option 3).
   - First break the flagged secondaries down by dt and A2/A1 in the existing npz, with d_pu per bin. Then rescan run 025320 at MIN_SEP 12 and 8 ns.
   - Success: either the flags are shown to be real small pulses (and option 3 closes), or a stricter acceptance takes d_pu to 0 and improves CB at 2.24 / 3.49 GeV by up to 0.04 pp without losing Ta pile-up efficiency.

## 7. Open questions and gaps

**Questions for the user.**
- Q1. Is the lead glass in the PRad-II physics acceptance? The proposal planned an all-PbWO4 HyCal, but the data have 576 G modules, all with calibration factor 0.
- Q2. The 10-29 refinement-off decision: may refinement come back if it wins out-of-sample?
- Q3. The Møller rejection and NL switches are "pending confirmation" (`energy_resolution_benchmark.py:153,155`). Which default should the plan assume?
- Q4. Do PRad-II cosmic, snake-scan or pi0/Compton runs exist?
- Q5. Are HyCal temperature and HV logs available, given the -2 %/°C light yield?
- Q6. What is the fADC zero-suppression threshold per channel? Threshold loss feeds b/E and the NL (FT-Cal, NPS).
- Q7. Which metric matters for the cross section: the CB core, or the tail fraction (11-45 % outside 3σ)?

**Gaps.**
- The numbers below are not on disk:
  - the R5nr-ri tables;
  - sigma with NL on vs off;
  - a refine-on rerun on the new basis;
  - a 12 vs 8 ns real-data comparison;
  - the energy effect of the phase correction or of common-mode subtraction;
  - any per-ring resolution;
  - a lead-glass result;
  - a fourth beam energy.
- The PRad-I in-beam ep peak width at 1.1 / 2.2 GeV was not found in a primary source.
- The budget rows marked [I] are estimates and need the runs named in section 4.

## Sources (literature)
- P1 Kubantsev et al., PrimEx HyCal, https://arxiv.org/abs/physics/0609201
- P2 Li Ye, APS 2017, https://wiki.jlab.org/pcrewiki/images/6/61/APS_2017_LiYe.pdf (preliminary)
- P3 Weizhi Xiong, PRad calibration slides, https://wiki.jlab.org/pcrewiki/images/d/d7/Weizhi_10262016_.pdf
- P4 = P1
- P5 PRad-II proposal, https://arxiv.org/abs/2009.10510; technical review https://wiki.jlab.org/pcrewiki/images/c/c3/PRad_II_Technical_Review_rev14.pdf
- P6 HyCal calibration proposal, https://www.jlab.org/primex/weekly_meetings/tff/docs/HyCal_calibration_proposal.pdf
- Xiong and Peng review, https://arxiv.org/abs/2302.13818; PrimEx fit, https://arxiv.org/abs/2108.13276 (Fig. 14; term assignment inferred)
- N1 NPS note, https://mailman.jlab.org/pipermail/nps/attachments/20250821/6399aa15/attachment-0001.pdf; N2 https://indico.jlab.org/event/946/contributions/16506/attachments/12587/20068/20250505_pi0_Calibration_Hao_Huang.pdf; N3 https://hallcweb.jlab.org/wiki/images/archive/2/22/20230710133446%21NPS_ElasticCalibrations.pdf; N4 https://indico.jlab.org/event/758/contributions/13811/attachments/10628/16056/NPS_Performance_2024.pdf
- H1 HPS ECal, https://arxiv.org/abs/1610.04319
- F1 CLAS12 FT, NIM A 959 (2020) 163475, https://eprints.whiterose.ac.uk/157258/1/main.pdf
- G1 GlueX, https://arxiv.org/abs/2005.14272; G2 https://halldweb1.jlab.org/wiki/images/8/83/Fcal_cal_12_Feb_2015.pdf; G4 GlueX ECAL, https://arxiv.org/abs/2605.25301; C2 CCAL draft, https://halldweb.jlab.org/detectors/fcal2/conferences/ccal_nps/ccal_nps_ver2.pdf; B1 BCAL, https://arxiv.org/abs/1801.03088
- M1 CMS ECAL Run 1, https://arxiv.org/abs/1306.2016; M2 multifit, https://arxiv.org/abs/2006.14359; M3 Run 2, https://arxiv.org/abs/2403.15518
- A1 ALICE PHOS, https://arxiv.org/abs/1902.06145
- D1 PANDA EMC TDR, https://arxiv.org/abs/0810.1216
- O1 Cleland and Stern, NIM A 338 (1994) 467, doi:10.1016/0168-9002(94)91332-3 (via Delmastro, https://arxiv.org/abs/0812.3486)
- K1 NA48 LKr, https://arxiv.org/abs/hep-ex/0012011; K2 KTeV, https://arxiv.org/abs/hep-ex/0208007; K3 KTeV, https://arxiv.org/abs/1011.0127; K4 BaBar, https://arxiv.org/abs/hep-ex/0105044; K5 BaBar EMC, https://arxiv.org/abs/physics/0312128 and https://arxiv.org/abs/physics/0601138; K6 Belle II TDR, https://arxiv.org/abs/1011.0352; K7 Belle II ECL pulse-shape, https://arxiv.org/abs/2007.09642; K8 COMPASS, https://arxiv.org/abs/1410.1797; K9 PHENIX EMCal, https://arxiv.org/abs/nucl-ex/0202009; K10 KOTO, https://arxiv.org/abs/1509.03386; K11 Mu2e, https://arxiv.org/abs/1802.06341; K12 GRAAL, https://arxiv.org/abs/physics/9709027; K13 A2 at MAMI, https://arxiv.org/abs/1805.10504 and https://arxiv.org/abs/1506.08849
- PDG PbWO4: https://pdg.lbl.gov/2024/AtomicNuclearProperties/HTML/lead_tungstate.html (X0 0.89 cm, R_M 1.96 cm; HyCal crystal 18 cm = 20.2 X0)
