# Position-correction form per run type: decision and the resulting energy resolution

**Date:** 2026-11-06
**Context:** the r5ri chain (`output/r5ri_20261104/`, runinfo calibration basis) on all seven points; form study extended from 3.49 GeV H2 to every run type in `output/poscorr_form_20261105/`. Algorithm and definitions: `notes/position_correction_20261105.md`.
**Source of the numbers:** `output/poscorr_form_20261105/README.md` and `tables/` (held-out χ² and paired σ tables, cited per table below); `output/r5ri_20261104/tables/ri_table.md` and `ri_splits.tsv` for the resolution.

Paths are relative to `analysis/pyscripts/`; `$F` = `output/poscorr_form_20261105`, `$R` = `output/r5ri_20261104`. σ_E/E is the Crystal Ball (CB) core width of the inner-ring elastic peak, mean ± SE over splits, in %; "[G]" = iterative ±2σ Gaussian; Δ in pp, paired per split. Held-out χ²/ndf is the 5-fold per-module figure of merit of the form study (6 × 6 grid, bins with ≥ 30 events, event-mean normalisation, no parameters subtracted; `$F/tables/<slug>_heldout.md`, footnote). All numbers are on the runinfo basis; none is compared with a result from before 2026-11-04.

## 1. Decision

The forms already in production in the r5ri chain are kept. Per run type:

- **3.49 GeV H2, 2.24 GeV H2, 2.24 GeV C: hybrid poly8.** A per-module eight-coefficient polynomial for every seed module with ≥ 3000 clusters, the pooled 12 × 12 bilinear map for the rest. poly8 has the lowest held-out χ²/ndf of the four forms in all three run types and the lower σ_E/E in the benchmark (sections 2 and 3).
- **0.73 GeV H2: per-run-type bilinear.** Bilinear and poly8 tie in the held-out test (section 2), and no σ scan exists for this point; the simpler form stays.
- **Ta 25800, 26138, 26230 (each run separately, never pooled): the 2.24 GeV H2 bilinear map, by symlink.** No Ta module reaches 3000 clusters, so a per-module form cannot be fitted; a Ta-trained pooled poly8 is not better than the H2 map where it can be tested (section 2).
- **even6 and even4** (c6 = c7 = 0; c3 to c7 = 0) beat poly8 in no run type and are not used.

The R7ri variant (energy-dependent W0) is a separate decision, not taken here (section 4).

## 2. Held-out form comparison

Median held-out χ²/ndf over the modules with ≥ 3000 clusters (`$F/tables/all_runtypes_summary.md`; per-type tables `3p49_heldout.md`, `2p24_H2_heldout.md`, `2p24_C_heldout.md`, `0p73_heldout.md`):

| run type | modules (≥ 3000) | bilinear | poly8 | even6 | even4 | poly8 best in |
|---|---|---|---|---|---|---|
| 3.49 GeV H2 | 52 | 6.587 | 1.792 | 3.444 | 4.124 | 49 / 52 |
| 2.24 GeV H2 | 76 | 3.233 | 1.633 | 2.434 | 2.494 | 66 / 76 |
| 2.24 GeV C | 28 (27 usable) | 2.633 | 1.705 | 2.117 | 2.198 | 22 / 27 |
| 0.73 GeV H2 | 16 | 1.759 | 1.646 | 1.871 | 1.859 | 7 / 16 |

At 0.73 GeV the paired bilinear − poly8 difference is median −0.080, mean +0.051 ± 0.106, bilinear lower in 9 of 16 modules (`0p73_heldout.md`): a tie. At 2.24 GeV C one module-fold poly8 fit failed in the held-out test, so 27 of the 28 modules enter the comparison (`2p24_C_heldout.md`); the production map has 28 own fits and no failure (`stats.failed` of the production JSON, section 5).

**Ta** (`$F/tables/Ta_pooled_heldout.md`). The largest module holds 917 / 204 / 549 clusters in 25800 / 26138 / 26230, so the per-module test is impossible and the per-module forms cannot be fitted. The substitute pools all modules' 6 × 6 bins and compares the production H2 bilinear map with a pooled poly8 trained on the Ta train folds:

| sample | clusters (r window) | held-out bins | bilinear (2.24 H2 map) | pooled poly8 (Ta-trained) |
|---|---|---|---|---|
| Ta 25800 | 16 289 | 150 | 1.382 | 1.440 |
| Ta 26138 | 2 215 | none | n/a: no held-out 6 × 6 bin with ≥ 30 events | n/a |
| Ta 26230 | 3 389 | 63 | 1.648 | 1.785 |

The pooled poly8 is not better than the symlinked map on either testable run, and 26138 cannot be tested. χ² per bin shrinks with the cluster count (Ta has 50 to 400 times fewer clusters than the 2.24 GeV H2 reference in the same table), so these values are not comparable in absolute terms with the 2.24 GeV H2 ones.

Caveats from `$F/README.md`: for 0.73 GeV and Ta the bilinear column is the production map as fitted on the full sample, not refitted per fold, so its χ² is slightly optimistic; the 3.49 and 2.24 GeV bilinear columns carry the same bias (footnote of each `*_heldout.md`). The Ta per-module tables (`Ta25800_heldout.md` etc.) list only cluster counts.

## 3. Benchmark σ_E/E by form

Paired per-split comparison of the production hybrid (poly8) against the bilinear map run through the frozen benchmark on the same splits, legacy path (`$F/tables/3p49_sigma.md`, `2p24_H2_sigma.md`, `2p24_C_sigma.md`). even6 was run at 3.49 GeV H2 only (1.6328 ± 0.0067, Δ +0.0418 ± 0.0037); even4 was never run through the benchmark; 0.73 GeV and Ta have no σ scan.

| run type | splits | CB poly8 | CB bilinear | Δ(bilinear − poly8) | bilinear worse in | gated poly8 / bilinear | Δ gated |
|---|---|---|---|---|---|---|---|
| 3.49 GeV H2 | 27 | 1.5910 ± 0.0071 | 1.6985 ± 0.0073 | +0.1075 ± 0.0055 | — | 1.5907 / 1.6977 | +0.1070 ± 0.0055 |
| 2.24 GeV H2 | 42 | 1.7938 ± 0.0037 | 1.8353 ± 0.0038 | +0.0415 ± 0.0021 | 42 / 42 | 1.7917 / 1.8335 | +0.0418 ± 0.0021 |
| 2.24 GeV C | 8 | 1.7995 ± 0.0078 | 1.8193 ± 0.0087 | +0.0198 ± 0.0033 | 7 / 8 | 1.7985 / 1.8182 | +0.0198 ± 0.0035 |

The 3.49 GeV table carries no per-split count column. The Gaussian deltas agree in sign and size (+0.1013, +0.0380, +0.0194 pp legacy). Each map was fitted on all splits of its run type, so the absolute σ values are in-sample; the paired differences between forms are the comparison.

## 4. Energy resolution with the chosen forms

The chosen forms are the production configuration of the r5ri chain, so the R5nr-ri reference arm is the result; nothing was rerun. Basis: the calibration file that `database/runinfo/general.json` assigns to each run with nl1/nl2 applied, benchmark sha1 `3293ab96870cef4580e886eb2dda5ec37c8821d1`, `--nl-members on`, Ta per run, never pooled (`$R/README.md`). Numbers from `$R/tables/ri_table.md`; split counts checked against `ri_splits.tsv` (236 rows = 118 splits × legacy/gated; 25 Ta splits = 8 + 4 + 13).

| run type | splits | legacy CB | legacy [G] | gated CB | gated [G] |
|---|---|---|---|---|---|
| 3.49 GeV H2 | 27 | 1.591 ± 0.007 | 1.638 | 1.591 ± 0.007 | 1.638 |
| 2.24 GeV H2 | 42 | 1.794 ± 0.004 | 1.837 | 1.792 ± 0.004 | 1.836 |
| 0.73 GeV H2 | 16 | 3.417 ± 0.010 | 3.448 | 3.417 ± 0.010 | 3.448 |
| 2.24 GeV C | 8 | 1.800 ± 0.008 | 1.839 | 1.798 ± 0.008 | 1.839 |
| Ta 25800 | 8 | 2.051 ± 0.024 | 2.066 | 1.961 ± 0.033 | 2.014 |
| Ta 26138 | 4 | 2.111 ± 0.078 | 2.118 | 2.062 ± 0.076 | 2.088 |
| Ta 26230 | 13 | 2.142 ± 0.060 | 2.137 | 1.962 ± 0.060 | 1.998 |

`ri_table.md` gives no SE for the Gaussian means. The benchmark's preferred fitter is the Gaussian at 2.24 and 0.73 GeV H2; CB is the headline here as in the rest of the chain (`$F/README.md`, deviation 5).

**R7ri.** The energy-dependent W0 (per-type `log_weight_thres`) is a separate decision not taken in this note. From the same table, legacy CB: 1.558 ± 0.009 at 3.49 GeV H2, 1.782 ± 0.004 at 2.24 GeV H2, 1.784 ± 0.010 at 2.24 GeV C (Δ vs R5nr-ri −0.033 ± 0.004, −0.012 ± 0.002, −0.016 ± 0.004 pp); 0.73 GeV H2 and the three Ta runs are within their errors (+0.004 ± 0.005, +0.004 ± 0.010, −0.022 ± 0.041, +0.002 ± 0.018 pp).

## 5. Pointers

- Algorithm, definitions, constants and decision history: `notes/position_correction_20261105.md`.
- Form study: `$F/README.md` (method, commands, deviations), `$F/tables/` (held-out, per-module TSV and σ tables), `$F/tools/` (`heldout.py`, `heldout_rt.py`, `heldout_ta.py`, `sigma_table.py`, `sigma_table_rt.py`), bilinear benchmark arms under `$F/<slug>/bilinear/`.
- Production maps, `$R/poscorr/maps/` (`created` 2026-11-05):
  - `3p49GeV_H2/poscorr_r5ri_3p49GeV_H2_all.json` (hybrid, 52 own fits of 170 seeding modules), fallback `poscorr_r5ri_bilinear_3p49GeV_H2_all.json`, per-module `poscorr_r5ri_3p49GeV_H2_all_modules.tsv`
  - `2p24GeV_H2/poscorr_r5ri_2p24GeV_H2_all.json` (hybrid, 76 own fits of 164), fallback `poscorr_r5ri_bilinear_2p24GeV_H2_all.json`, `poscorr_r5ri_2p24GeV_H2_all_modules.tsv`
  - `2p24GeV_C/poscorr_r5ri_2p24GeV_C_all.json` (hybrid, 28 own fits of 191), fallback `poscorr_r5ri_bilinear_2p24GeV_C_all.json`, `poscorr_r5ri_2p24GeV_C_all_modules.tsv`
  - `0p73GeV_H2/poscorr_r5ri_0p73GeV_H2_all.json` → symlink to `poscorr_r5ri_bilinear_0p73GeV_H2_all.json` in the same directory
  - `2p24GeV_Ta_X17/poscorr_r5ri_2p24GeV_Ta_X17_all.json` → symlink to `../2p24GeV_H2/poscorr_r5ri_bilinear_2p24GeV_H2_all.json`; the Ta jobs of all three runs use it
- Own-fit counts and `stats.failed` read from the `stats` block of each hybrid JSON; the symlink targets from `ls -l`.
