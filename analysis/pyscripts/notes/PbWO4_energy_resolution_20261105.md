# PbWO4 energy resolution: runinfo calibration basis, R5nr-ri / R7ri / fAri, template shape cut and tail ripple, one-pulse templates, Ta timing rerun, position-correction form

**Date:** 2026-11-04 and 2026-11-05 (one note for both days; the position-correction note and the start of its extension, finished early 2026-11-06, are the last items)
**Branch:** `waveform_jz` (HEAD `2cc8d11`; the scripts of §2, §3 and §6 landed in `469dc71`, the notes and `Commands.md` in `2cc8d11`, both 2026-11-06 06:15)
**Previous notes:** `notes/PbWO4_energy_resolution_20261030.md` (R5nr, Møller, MF fit quality), `notes/rerun_inventory_20261104.md` (inventory under the calibration rule), `notes/position_correction_20261105.md` (algorithm and figures)

Paths are relative to `analysis/pyscripts/` in `/Users/L00338853/work/PRad/prad2evviewer` unless stated; `$R` = `output/r5ri_20261104`, `$T` = `output/templates_spline_best50_p1r3_20261104`. σ_E/E is the CB core width of the inner-ring elastic peak, mean ± SE over splits, in %; Δ in pp, paired per split; "[G]" = iterative ±2σ Gaussian. All numbers in this note are on the runinfo calibration basis; by the standing rule none is compared with a result from before 2026-11-04.

## 1. Status

- **Calibration basis switched (§2).** Standing rule: every job uses the calibration file that `database/runinfo/general.json` assigns to its run, with nl1/nl2 applied. `calib_resolve.py` resolves it; the benchmark gained `--nl-members on|off` (default off, byte-identical default path) and is frozen at sha1 `3293ab96870cef4580e886eb2dda5ec37c8821d1`.
- `notes/rerun_inventory_20261104.md` grades every earlier study against the rule; `plans/plan_ri_20261104.md` defines the rerun P0–P7. P1 and P2 are done; P6 is done except the RF e100 tier (§9).
- **Templates (§3).** The rising ≥ 3 shape cut has no material effect. The wavy tails are real per-module structure plus a junction feature at the data end, inferred to be a build artefact, in both arms. Wassim's one-pulse templates ripple more, not less.
- **Resolution on the new basis (§4).** R5nr-ri is the reference, all seven points. R7ri (energy-dependent W0) improves 3.49 GeV H2, 2.24 GeV H2 and 2.24 GeV C by −0.033, −0.012 and −0.016 pp and is neutral elsewhere. Dropping f(A) is neutral (e-p/C within ±0.004 pp, Ta within its errors). a/√E ⊕ b does not describe the H2 points. The one-pulse bank gives 1.584 ± 0.005 at 3.49 GeV with f = 1.
- **P6 Ta timing outputs (§5)** exist on the new basis: module timing, seed-time dt and pile-up rates for 26138/26230.
- **Position correction (§6).** Held-out form study at 3.49 GeV H2: per-module poly8 beats bilinear, even6 and even4. Explanatory note with figures written. The extension to the other run types started 2026-11-06 and is in progress.
- **Nothing adopted beyond the basis rule.** R7ri and the poly8 form are findings, not decisions.

## 2. Calibration basis

- **Resolver.** `calib_resolve.resolve_calib(run)` applies the C++ `RunInfoConfig` rule (largest `from_run` ≤ run, chained inheritance; `calib_resolve.py:1-4`). It exits on a missing file or a file with fewer than 100 modules with factor > 0; `check_pipeline` exits when under half of the PbWO4 modules carry `cal_factor > 0`.
- **Run → file** (`rerun_inventory_20261104.md` §1; sha1 of each file in `$T/calib_table.tsv`): 24386–24512 and, by inheritance from 24551, 24590/24603 → `calibration_factor_3p5_June7_nonLinear.json`; 25228 → `PRad2_2p2_round1_1`, 25257/25260 → `round1_2`, 25299/25308 → `round1_3`, 25320 → `round1_4`; 25800/26138/26230 → `x17_25825/26094/26214.json`.
- **nl1/nl2.** Only the June7 file carries them (nl1 in 1006 modules, nl2 in 986; median factor 1.06 at 0.73 GeV, 1.008 at 2.24 GeV). The `round1_*` and `x17_*` files have nl = 0, so the rule changes the 3.49 and 0.73 GeV numbers only (inventory §1).
- **Benchmark edit.** `--gain-file` rebuilt the cluster energy from member energies and lost nl (inventory §4.1). `benchmark_hycal_timing_test_20261021.py` was edited in place, +9/−1 (`git show 469dc71`): `--nl-members on` scales the member energies of every cluster record by `r.hit.energy / Σ h.energy`, the C++ non-linearity factor, so the gain fit and the `--gain-file` passes carry it. The `nl_members=on` header token is written only when on.
- **Gates** (`plans/plan_ri_20261104.md` Log, `$R/README.md`): G0, default path byte-identical to the pre-edit run (result_window.tsv, result_gains.tsv, dump arrays); G1, nl identity at 0.73 GeV, max deviation 3.9e-7, f ≠ 1 in 96.4 % of clusters, median f 1.062; G1b, f = 1 exactly at 2.24 GeV. `run1_ri.sh` refuses any other benchmark sha (rc 97; `Commands.md:141`).
- **Inventory.** No study ran with zero calibration. The energy-dependent studies of 10-21 to 11-04 used the three `*_fadc` files and the `x17_*_fadc` twins, none of which general.json names. Verdicts: (i) rerun, (ii) regenerate later, (iii) no rerun, (iv) invalid as recorded.
- The (iv) rows: the 2.24 GeV scans run with the 3.5 GeV file, the pre-window sweep on a private file, the Ta 26138/26230 jobs on the 25825 period, and the claim "non_linear_corr off = exact null" as a general statement (inventory §2).
- **Rerun list** (inventory §5): P1 R5nr-ri; P2 R7 and f(A) on it; P3 0.73 GeV Møller block; P4 energy_bias / Geant4 profile / R6; P5 non-linearity steps 2–3; P6 Ta timing; P7 low priority.
- **Ta file variant.** runinfo names `x17_*.json`; the `_fadc` twins are about 6 % larger per module (median ratio 1.062–1.064, per-module std 0.035–0.049; `plan_ri_20261104.md` Open Questions 2). Run as written; whether runinfo should name the `_fadc` files is open.

## 3. Template extraction

**Shape cut** (`$T/report.md`): arm B (`--p1-shape-cut --p1-min-rising 3`) against arm A (baseline), five run types, 238 select splits, 2 h 43 min, 2.2 GB.

- The cut rejects 80–88 % of candidates. Stored pools shrink from 339 914 to 249 298 (0.73 GeV), 1 899 227 to 1 522 107 (2.24 GeV H2), 261 225 to 200 606 (C), 595 279 to 408 519 (Ta), 789 258 to 679 279 (3.49 GeV).
- Run-type-average fallbacks rise by 65 modules at 0.73 GeV and by 5–25 elsewhere; the median χ²/dof moves by at most 5.5 %; G1-5 passes in both arms for every run type and G1-8 is identical in A and B.
- The B template is bit-identical to A in 850/1054 modules (2.24 GeV H2) and 617/866 (3.49 GeV), because the cut almost never touches the 50 highest-energy pulses that form a pool (`$T/tail_ripple_diag.md`). No material effect.

**Tail ripple** (`$T/tail_ripple_diag.md`; 2.24 GeV H2 and 3.49 GeV H2, fallback-0 modules):

- Ripple of 0.93 % (2.24 GeV) and 0.92 % (3.49 GeV) of peak (half p-p, 30–150 ns, median), the same in A and B. Dominant periods 30–90 ns, no module below 10 ns: neither the 4 ns sampling phase nor the knot spacing. Independent of the selected λ (0.95 / 0.99 / 0.95 % at 1e2 / 1e3 / 1e4).
- It is real: the tail ripple of two random halves of a pool correlates with median r = 0.91 (IQR 0.79–0.96; 58 modules).
- A second component is a junction feature where the data end (T_data ≈ 245–250 ns) and the exponential pseudo-points begin; 69 % of modules have a rising segment within T_data − 20 … + 40 ns.
- The log-y tail panel makes a 0.1 %-of-peak bump at y ≈ 0.01 look like a 10 % wiggle.
- Fix options, not implemented: stronger penalty or fewer knots beyond ~100 ns, or start the fitted exponential ~40 ns earlier.

**One-pulse (Wassim) templates** (`wassim_one_pulse_template.py`; `output/templates_wassim1_20261104/3p49GeV_H2/README.md`, `summary.tsv`):

- The `best_ref_shape_fit.C` candidate and ranking rules, one pulse per module over the 28 splits of 24386/24503/24512, interpolating cubic spline through the 100 samples, no pooling. About 50 min wall.
- 1148/1152 modules get a pulse (982 tier 2, 166 tier 1; 22 418 candidates). W69, W1055, W1123 and W1124 have none and fall back to the run-type average.
- The pipeline ripple counter flags 183/1148 modules against 0/1052 in the p1shape_r3 bank (built with `--osc-escalate`): the single pulse's noise is carried into the tail.

**Template PDF.**

- `spline_template_extract.py plot` re-renders the module pages of a bank without a rebuild and writes `onset_abs_ns.json` (absolute onset in the readout window = mean mf_t0 of the module's select pulses + bank onset shift).
- `output/templates_spline_best50_20261026/3p49GeV_H2/templates_3p49GeV_H2.pdf` was remade on 2026-11-05 with a linear-y 0–400 ns panel showing the template at its onset (dashed) and T_data (dotted); the previous version is kept as `templates_3p49GeV_H2_prev.pdf`.
- The module pages carry no h50 two-tau overlay (`spline_template_extract.py:2546-2576`); the one-pulse PDF omits it by request.

## 4. Energy resolution on the new basis

Tables: `$R/tables/ri_table.md` (per point, both paths, Δ per variant, a/√E ⊕ b), `ri_splits.tsv` (per split), `$R/README.md` (basis, commands, gates). Three chains of 119/118/118/1 jobs (dump, map, final, repro); Ta per run, never pooled.

- **R5nr-ri reference** (legacy path): 3.49 GeV H2 1.591 ± 0.007 (n 27), 2.24 GeV H2 1.794 ± 0.004 (42), 0.73 GeV H2 3.417 ± 0.010 (16), 2.24 GeV C 1.800 ± 0.008 (8), Ta 25800 2.051 ± 0.024 (8), Ta 26138 2.111 ± 0.078 (4), Ta 26230 2.142 ± 0.060 (13).
- Gated path within 0.002 pp of legacy for e-p/C; Ta gated 1.961–2.062.
- **R7ri** (energy-dependent W0 via per-type `log_weight_thres` configs, `Commands.md:170-174`): Δ −0.033 ± 0.004 (3.49 GeV H2), −0.012 ± 0.002 (2.24 GeV H2), −0.016 ± 0.004 (C); 0.73 GeV and Ta within their errors. The Gaussian deltas agree.
- **fAri** (banks without the `ampcorr_*` keys, f = 1): e-p/C Δ within ±0.004 pp on the legacy path; Ta −0.014 ± 0.013 to +0.011 ± 0.019, none beyond its SE. Neutral.
- **a/√E ⊕ b on the three H2 points:** R5nr-ri CB a = 2.858 ± 0.011 %·√GeV, b = −0.504 ± 0.030 %, χ² 677.71 for 1 dof; R7ri 408.27; the Gaussian series 759.81 and 585.28. b² < 0 and χ² of several hundred for one degree of freedom: the form does not describe the points.
- **One-pulse bank, 3.49 GeV H2** (`output/w1ri_20261104/tables/w1_table.md`, `README.md`): 1.584 ± 0.005 legacy, 1.582 ± 0.006 gated [G 1.622 ± 0.006], 27 splits, the r5ri chain with the bank swapped. Two things differ from R5nr-ri, the template shape and f = 1, because `ampcorr` refuses a bank without `ampcorr_meta`. Runtime 08:00–12:07 at P = 4.

**Chain mechanics** (`$R/tools/`):

- `ri_jobs.py` writes the job lists; `chain_ri.sh` runs one stage with gates G2–G6 (14 gain tables converged, 5 maps normalised and populated, 118 final rc = 0, repro twin byte-identical); `drive_ri.sh` runs the stages in order, re-reading P from `logs/P`; `run1_ri.sh` checks the benchmark sha.
- `launch_after_templates.sh` waited for the template chain's CHAIN_DONE (20:43 on 11-04), then started r7ri and fAri at P = 5 and raised r5ri to P = 6.
- r5ri dump started 18:56 on 11-04 and the repro finished 14:28:41 on 11-05; r7ri and fAri finished 15:03 (`logs/chain.out`). 5.2 GB per tree.
- **Overnight stall.** The gains gate compared `wc -l` output with `14`; macOS pads it with spaces, so the gate failed after the dump stage (22:54 on 11-04) and the chain idled until the fix (`| tr -d " "`, `chain_ri.sh:34`) and the resume at 08:00:42 on 11-05: nine hours lost.

## 5. Ta timing and pile-up on the new basis (P6)

- `output/module_timing_ri_20261104/`: `module_timing_dump.py` on the 17 splits of 26138 (4) and 26230 (13), started 18:31 on 11-04, 0 failed (`dump_run.log`).
- Veto-referenced summary, time-calib off: 301 modules fitted, σ median 2.476 ns, 16–84 % range [1.728, 2.995] ns (`module_timing_veto_summary.txt`); tcal, typealign, nogate and vs-integral variants alongside.
- `output/tseed_dt_ri_20261104/tables/tseed_dt.md`: T_module − T_seed for the 17 splits, raw and calibrated, gated and legacy; Ta 26138 gated robust σ 2.721 ns raw, 0.824 ns calibrated. Queued after the dump by `tools/chain.sh`.
- `output/pileup_rates_ri_20261104/report.md`: three splits per run, 60 000 physics events each. Events with ≥ 2 MF pulses 79.14 ± 0.10 % (26138) and 89.50 ± 0.07 % (26230); channels with a near extra pulse 4.080 / 6.233 %.
- `pileup_rates.py` now resolves the calibration through `calib_resolve` and checks that it loaded. The RF-comb footnote crashed on an empty H2/C pair list in a Ta-only run and now prints `n/a` (`git show 469dc71`).

## 6. Position correction

**Form study, 3.49 GeV H2** (`output/poscorr_form_20261105/README.md`, `tables/`):

- `position_corr.py fit --per-module --form {poly8,even6,even4}` added (`HYB_FORMS`, zero-filled higher coefficients). The default output is unchanged; the poly8 refit equals the stored r5ri hybrid exactly.
- Part A, 5-fold held-out χ²/ndf over the 52 modules with ≥ 3000 clusters (`tables/3p49_heldout.md`): medians bilinear 6.587, poly8 1.792, even6 3.444, even4 4.124; poly8 best in 49 of 52 modules. The bilinear map is the production one, not refitted per fold, so its χ² is optimistic.
- Part B, the benchmark on the 27 paired splits (`tables/3p49_sigma.md`, CB legacy): poly8 1.5910 ± 0.0071, even6 1.6328 ± 0.0067, bilinear 1.6985 ± 0.0073; Δ vs poly8 +0.0418 ± 0.0037 (even6) and +0.1075 ± 0.0055 (bilinear). even4 was not run through the benchmark.

**Explanatory note** `notes/position_correction_20261105.md`, figures in `notes/figures/poscorr_20261105/` (made by `poscorr_explain_plots.py` from the r5ri chain; numbers in `facts.md` / `facts_data.json`):

- Definitions, the bowl profile, the per-module poly8 and its fallback, the JSON schemas, the commands, a constants table and the decision history.
- Its §6 finds no recorded rationale for most constants, among them the 3000-cluster threshold, the accepted range [0.9, 1.1], the clamp 0.6, the 12 × 12 map binning and the 6 × 6 module grid.

**Extension to the other run types** (2.24 GeV H2, 2.24 GeV C, 0.73 GeV H2, Ta) started on 2026-11-06: held-out tables and the 2.24 GeV bilinear arms are being produced (`tables/*_heldout.md`, `2p24GeV_*/bilinear/`, `run_2p24.log`). In progress; no numbers here.

## 7. Documentation

- `Commands.md` was rewritten as the command reference for the chains: environment, template extraction, matched filter and timing dump, the energy-resolution chain with its R7ri / fAri / w1ri variants, timing and pile-up tools, non-linearity tools. Every line is taken from a script or driver; reconstructed lines are marked.
- The earlier running log (through 2026-11-03) is preserved unchanged as `Commands_log_pre20261105.md`.

## 8. Process

- Agent reports were lost in transit several times and had to be recovered from the task output files. The session's record of a step is the output tree and its logs, not the report.
- The gate bug of §4: a gate that compares a number should normalise it (`tr -d " "`), and a failed gate should be visible in the morning without reading `chain.out`.

## 9. Open points

- **P3** 0.73 GeV Møller block (kinematic rejection, two-peak CB, over-subtraction, CB-excess attribution) on R5nr-ri; nl at 0.73 GeV may change it.
- **P4** energy_bias / Geant4 profile / R6 on the R5nr-ri gains.
- **P5** non-linearity steps 2–3 (seed ratios, beam offsets) from R5nr-ri clusters.
- **P7** low-priority regenerations: gate scans, CFD bands, clustering seed 30, best-50 Ta cap, pile-up E columns, mfstd if pursued.
- The matched filter runs in each of the four benchmark stages of a chain (dump, map, final, repro) on the same waveforms. Proposal: cache the MF output per split once and let the later stages read it; not designed yet.
- Ta crosstalk: the crosstalk study has no Ta split (inventory §3); the Ta pile-up maps of §5 show the modules to start from.
- RF per-run calibration (e100 tier) for 26138/26230 on the runinfo file is still to be rerun (inventory P6 row); the module-timing dump of §5 is its input.
- The position-correction constants without a recorded rationale (§6); the form study on the other run types; whether runinfo should name the `x17_*_fadc` files (§2).
- `pileup_rates_ri_20261104/report.md` carries the 10-27 title and command lines; fix the template string before the next run.
