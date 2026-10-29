# PbWO4 energy resolution: amplitude dependence of the pulse shape, best-50 Ta degradation, R4 promotion

**Date:** 2026-10-28
**Branch:** `waveform_jz` (HEAD `3566071`; everything since is uncommitted, as in 20261026 §1.2)
**Plan:** `analysis/pyscripts/plans/plan_best50_pipeline_20261026.md` (S2–S4)
**Previous notes:** `analysis/pyscripts/notes/PbWO4_energy_resolution_20261027.md` (§2 best-50 rescan, §11 doublet audit), `analysis/pyscripts/notes/PbWO4_energy_resolution_summary.md` (§4)

Paths are relative to `analysis/pyscripts/` in `/Users/L00338853/work/PRad/prad2evviewer` unless stated. σ/E, preferred fit, legacy/gated, Δ in pp and the run-type names are as in the summary note. R3 = `output/resolution_fit_20261026/r3_splits`; B50 = `output/templates_spline_best50_20261026/`. Ta runs are reported per run (25800, 26138, 26230), never pooled.

## 1. Status

User decisions (10-28): (1) resolve the Ta degradation of the best-50 bank and promote best-50 as R4; (3) check the amplitude dependence of the pulse shape with amplitude-binned templates. No C++ changes. Outcome:

- Part A (amplitude dependence of the pulse shape, §2): weak. FWHM −0.17 ± 0.02 ns from 490 to 1560 ADC (0.7 %), rise and w99 flat, f(A) flat in 300–2500 ADC; it explains ≲ 10 % of the per-module Ta − H2 scatter and not the 0.73 GeV bank offset. Banks unchanged.
- Part B (best-50 on Ta, §3): not a module, region, gain or selection effect; ascribed (inference) to the 50-pulse pools of the 8 central modules; the Ta bank rebuilt from the same pulses with `--n-max 200` is within 1 SE of R3 on all three Ta runs (legacy −0.003 ± 0.009 / −0.001 ± 0.006 / −0.002 ± 0.031 pp; n200 − R4 legacy −0.032 ± 0.022 / −0.056 ± 0.039 / −0.027 ± 0.027 pp).
- R4 promoted (§4): best-50 banks and gains for H2 and C, the n200 bank for Ta; equal to R3 on e-p; summary and plan updated.

## 2. Amplitude dependence of the pulse shape

- Question: pool amplitude mixes differ between run types (mean ln A of the shipped pools 7.19 for 2.24 H2, 6.47 Ta, 5.89 at 0.73 GeV; fraction < 600 ADC 0.12 / 0.45 / 0.71). Does the template shape depend on amplitude, and does that explain the per-module Ta − H2 scatter and the narrower 0.73 GeV bank? Script `template_amp_binned.py` (numpy/scipy; 36 s, 8 workers); output `output/templates_spline_best50_20261026/amp_binned/` (`report.md`, `modules.tsv`, `ta_h2_modules.tsv`, `h073_h2_modules.tsv`, `amp_binned.pdf/png`; one-run cross-check in `run25320/`).
- Method: 2.24 GeV H2 select pools (all 5 runs; tier ≥ 1, band b0 300–2500 ADC) split into 300–600 / 600–1200 / 1200–2500 ADC (median A 487 / 721 / 1558). Per module and bin 40 random pulses (+ 12 same-bin hold-out) go through `process_module` with the shipped options; the shipped pool and hold-out are excluded so bins are independent of the shipped template. 193 modules have ≥ 52 pulses in every bin. Noise: two disjoint 40-pulse draws of one bin (σ of a lo−hi difference, MAD): rise / w99 / FWHM 0.17 / 0.47 / 0.22 ns. Reproduction: `process_module` on the whole module and on the rebuilt pool + shipped hold-out gives the shipped template, max|Δf| = 0 (4 modules).
- Results (1200–2500 minus 300–600, 193 modules; median [16, 84 %], SE of the median; n beyond ±2σ, noise alone gives 4.5 per side):

| metric | Δ (ns) | n > +2σ / < −2σ | excess scatter | slope per ln A (ns) |
|---|---|---|---|---|
| rise 10–90 | +0.02 [−0.18, +0.46], SE 0.02 | 39 / 25 | 0.28 | +0.016 ± 0.017 |
| w99 | −0.01 [−0.44, +0.59], SE 0.05 | 10 / 7 | 0.21 | +0.008 ± 0.037 |
| FWHM | −0.17 [−0.41, +0.13], SE 0.02 | 4 / 24 | 0.16 | −0.112 ± 0.021 (run 25320: −0.148 ± 0.031) |
| t_pk | +0.10 [−0.32, +0.51], SE 0.04 | 12 / 9 | 0.16 | +0.095 ± 0.027 |
| max abs Δf (%) | 4.5 [2.4, 7.3] vs null 3.1 | | | |

- Answer, shape: weakly. The median FWHM is narrower by 0.17 ± 0.02 ns (0.7 %) from 490 to 1560 ADC (low amplitude wider in 24 modules, narrower in 4); median rise and w99 do not move (|Δ| ≤ 0.02 ns). The module-dependent part is larger than the median: per-module rise changes beyond noise by 0.28 ns of either sign (39 / 25 modules beyond ±2σ), max|Δf| 4.5 % vs 3.1 % for noise. The bank's own band templates (`shape_lo/top`, 100–400 vs 1000–2500) show the same sign and a smaller size (FWHM −0.12, rise +0.02 ns).
- Selection: the bank's rule (top mf_energy) inside a bin picks wider pulses at low amplitude (energy-top minus random FWHM +0.20 / +0.27 / −0.06 ns in the three bins, SE 0.01–0.03), so the energy-top lo−hi FWHM is −0.43 ns and its slope −0.38 ns per ln A.
- Answer, f(A): a shape-implied integral factor (`_fa_ratios` recipe: noise-free bin-template waveforms fitted with the module template, analytic walk-out integral; reproduces the stored r_lo/r_mid/r_top to 0.0003 and the f(A) steps 0.9946 / 1.0004) falls by 0.35 ± 0.06 % from the low to the high bin (run 25320 alone: 0.65 ± 0.09 %), core area −0.5 %, tail −1.1 %. The f(A) table is flat there (+0.03 %): its −0.54 % step sits at 223 → 505 ADC, below the b0 bins. So f(A) does not capture the 490 → 1560 ADC change; it is below 0.7 % and untested under 300 ADC.
- Answer, Ta − H2 (1082 modules, shipped banks): observed median rise / w99 / FWHM −0.01 / −0.04 / +0.03 ns, predicted from the amplitude slope × Δ⟨ln A⟩ (−0.72) −0.01 / −0.01 / +0.09 (energy-top slopes −0.03 / +0.15 / +0.32): no median offset to explain, FWHM over-predicted. Spearman ρ of Ta − H2 with the Ta pool fraction < 600 ADC: rise −0.16 (p 1e-7), w99 −0.01 (p 0.7), FWHM −0.09 (p 3e-3). With each module's own slope × Δ⟨ln A⟩ (n 193, bins independent of the H2 template): FWHM ρ +0.37 (p 1e-7), rise +0.15 (p 0.04), w99 +0.07 (p 0.4); predicted σ 0.04 vs observed 0.19 ns for FWHM, so amplitude accounts for ≲ 10 % of the per-module scatter and none of the rise / w99 scatter.
- Answer, 0.73 GeV: bank-average offset vs 2.24 H2 rise / w99 / FWHM observed −0.07 / −0.02 / −0.28 ns; predicted (Δ⟨ln A⟩ = −1.30) −0.02 / −0.01 / +0.15 (energy-top +0.49 for FWHM). The FWHM sign is wrong (lower amplitude gives wider, the bank is narrower), so amplitude does not explain the 0.73 offset; the conclusion rests on the FWHM sign (rise: predicted −0.02, energy-top slope −0.05, of −0.07 observed, so a rise contribution is not excluded).
- Caveats: bins start at 300 ADC (b0), so < 300 ADC is untested; one 40-pulse draw per module and bin.
- How to run: `MPLBACKEND=Agg python3 template_amp_binned.py` (cross-check: `--run 25320 --n-pool 25 --n-ho 8 --out output/templates_spline_best50_20261026/amp_binned/run25320`, 22 s). Figure: `amp_binned.png`.

## 3. Why best-50 is worse on Ta

Script `ta_b50_diagnose.py` (`modules`, `templates`, `gains`, `select`, `selvar`, `runs`; `--help`), output `output/ta_b50_diag_20261028/`. Ta is always per run (25800 / 26138 / 26230). R4 = B50 bank + B50 gains. Per-run Δ(R4−R3) legacy CB: +0.029 ± 0.015 / +0.055 ± 0.035 / +0.026 ± 0.036 pp (1.9 / 1.6 / 0.7 SE).

**(i) Where.** R3 and R4 events paired per event (seed, cluster position, raw energy; 99.8 % aligned), post energies, CB fit per module (n ≥ 150: 31 / 8 / 9 modules), paired bootstrap SE.
- Not a few modules and not a region: per-run paired Δ on the pooled events +0.019 ± 0.024 / +0.086 ± 0.058 / +0.055 ± 0.049 pp; r < 80 mm holds 70–100 % of the events (+0.023 ± 0.022 in 25800; 80–100 mm +0.104 ± 0.062, 25800 only). No module has |Δ| > 2 SE (per-module SE 0.06–0.37 pp); weighted mean module Δ +0.001 ± 0.018 / +0.027 ± 0.061 / +0.024 ± 0.042; 55 / 62 / 78 % of the fitted modules are positive.
- Of the 7 modules fitted (n ≥ 150) in all three runs, three are positive in all three: W599 (+0.26 ± 0.14, +0.23 ± 0.29, +0.17 ± 0.18), W592 (+0.07 ± 0.21, +0.18 ± 0.20, +0.00 ± 0.19) and W459 (+0.03 ± 0.10, +0.00 ± 0.14, +0.07 ± 0.18); W565 is positive in 25800 and 26138 (+0.17 ± 0.17, +0.05 ± 0.23) and −0.05 ± 0.16 in 26230 (`modules_legacy.tsv`). The 8 module-median modules: mean Δ +0.064 ± 0.026 (SD 0.073); the per-module Δ SD in H2 is 0.04–0.07 (2.24 GeV H2 mean +0.005 ± 0.004, 122 modules), so 8 modules cannot average it out.
- Single events: R4 − R3 post energy median ≤ 0.03 %, IQR/1.35 0.35–0.40 %. Pre arm (pooled gains only, no seed-indexed step) Δ (R4 − R3, legacy) is +0.012 / +0.031 / +0.000 pp, post +0.029 / +0.055 / +0.026.

**(ii) Templates.** Inner modules (r < 200 mm, 296 valid in both banks), R3 production vs B50, normalised shapes.
- Ta: max|Δf| median 4.73 % [p16–p84 3.2–6.1], Δrise −0.19 ns, Δt_pk −0.06 ns (|Δ| median 0.22). 2.24 GeV H2 (no degradation): 4.73 %, −0.21 ns, −0.02 ns (0.20). The change is generic, not Ta-specific.
- Spearman, module Δ vs |max|Δf||: −0.29 / −0.10 / +0.50 (n 31 / 8 / 9, p ≥ 0.11); vs |Δrise|, |Δw99|, |Δt_pk|, |ΔI_int|: no significant value. (Single-event energy shift SD vs |Δt_pk| and |ΔI_int|: ρ 0.41, p 0.01 in 25800.)
- Central 8 modules: R3 pool 200 pulses, B50 pool 50 (top by energy). Sub-sample phase concentration phase_R mean 0.12 (R3) → 0.27 (B50); I_int phase peak-to-peak of W565 / W592 / W599 0.023 / 0.024 / 0.018 → 0.035 / 0.029 / 0.027.

**(iii) Gains.** `percrystal_best50` / `percrystal_spl`, free modules (138 / 54 / 105): median ratio 1.0056 / 1.0029 / 1.0051, RMS 0.43 / 0.38 / 0.42 % (H2: 0.54–0.66 %), max |dev| 1.4–1.8 %, one outlier > 4 MAD per run. Spearman of |seed gain deviation| vs module Δ: +0.02 / −0.33 / +0.58 (p 0.92 / 0.42 / 0.10); 3×3 RMS: +0.22 / −0.05 / −0.03. The gains are not the cause.

**(iv) Selection.** Ta select logs: `--mf-strict` rejects 21.3 / 17.9 / 21.8 % of c3b survivors, `--max-modules-per-event 20` 0.27 / 0.20 / 2.52 % of events, the 1.2 E_beam cap 13 / 3 / 29 694 candidates (the open 29 710: 26230 pile-up).
- Build pool (tier ≥ 1 band b0, top 50 by energy, all runs): < 600 ADC is 32.8 / 40.4 / 32.4 % (source run) against 7.2–7.5 % for 2.24 GeV H2. By region: r < 200 mm 0.0 % < 600 in both (600–1200 ADC: Ta 5.3 / 10.1 / 8.2 %, H2 0 %); r ≥ 200 mm 56.0 / 56.6 / 49.1 % vs H2 10.0–10.7 %. The low-amplitude pool sits outside the Ta cluster region, where Ta has ~12 pulses ≥ 300 ADC per module; it is rate, not a cut.
- Relaxing one cut (25800.00100 / 26230.00100, stored top-50 per module): `--max-modules-per-event 40` +170 / +1198 candidates (60 / 64 % < 600 ADC); no `--mf-strict` +4001 / +3746 (73 / 69 % < 600); no cap swaps in 1 / 687, all > 1200 ADC at r < 200 mm (pile-up above 1.2 E_beam, which the cap should keep out). The r < 200 mm stored < 600 fraction moves by ≤ 0.9 pp. No cut restores a high-amplitude pool; no selection change is made.

**Cause (ascribed; inference, not shown directly).** The Ta result rests on 8 central modules (r 73–79 mm) whose B50 templates come from 50 top-energy pulses of ~1250 (R3: 200). That pool is noisier and its sub-sample phases are concentrated (phase_R 0.27); H2 and C average over ≥ 78 modules and do not see it. Test below: same pulses, `--n-max 200`; n200 − R4 legacy −0.032 ± 0.022 / −0.056 ± 0.039 / −0.027 ± 0.027 pp (same events), which does not by itself separate a noisier bank from other changes.

**Fix tested (no code change).** Rebuild of the Ta bank from the B50 Ta pulses with the B50 options except `--n-max 200` (`n200`; control `phase50`: `--pool-order energy-phase`, n 50), fresh Ta gains, 25 Ta splits with the R3 `job.cmd` (stem, Ta `--mf-templates`, `--gain-file` changed). Paired mean ± SE over the splits of each run, CB post:

| Run (splits) | R3 | R4 Δ legacy | n200 Δ legacy | within 1 SE | R4 Δ gated | n200 Δ gated | phase50 Δ legacy | phase50 Δ gated |
|---|---:|---:|---:|---|---:|---:|---:|---:|
| 25800 (8) | 1.986 | +0.029 ± 0.015 | −0.003 ± 0.009 | yes | +0.026 ± 0.018 | −0.030 ± 0.021 | +0.018 ± 0.013 (no) | +0.033 ± 0.013 |
| 26138 (4) | 1.950 | +0.055 ± 0.035 | −0.001 ± 0.006 | yes | +0.065 ± 0.024 | −0.029 ± 0.028 | +0.002 ± 0.048 (yes) | +0.013 ± 0.036 |
| 26230 (13) | 1.983 | +0.026 ± 0.036 | −0.002 ± 0.031 | yes | +0.115 ± 0.037 | +0.061 ± 0.038 | −0.006 ± 0.047 (yes) | +0.089 ± 0.043 |

- n200 − R4 (same events): −0.032 ± 0.022 / −0.056 ± 0.039 / −0.027 ± 0.027 legacy, −0.055 ± 0.031 / −0.094 ± 0.027 / −0.054 ± 0.043 gated. Event-level fits per run: n200 − R3 −0.000 ± 0.019 / +0.041 ± 0.059 / +0.023 ± 0.048 (R4 +0.019 / +0.086 / +0.055).
- Central phase_R 0.09 and I_int phase p-p 0.015–0.019 (W565 / W592 / W599) with n200; yields: events 0.9957 and peak 0.992 of R3.
- Verdict: the n200 Ta bank is within 1 SE of R3 (legacy CB) in 3 of 3 runs; gated 26230 stays +0.061 ± 0.038. Phase balance alone (phase50) is not enough in 25800. Per-module checks are below the noise (W599 stays positive in 25800 and 26230).
- Limits: the evidence rests on SEs of 0.006–0.036 pp per run; 26138 has 4 splits. H2 and C keep n 50 (equal to R3 already). The G1-7 gate still fails as for B50.

### 3.1 How to run

```bash
cd /Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts && source setup_env.sh && export MPLBACKEND=Agg
python3 ta_b50_diagnose.py modules --boot 100; python3 ta_b50_diagnose.py modules --boot 100 --path gated   # ~1-2 min each
python3 ta_b50_diagnose.py templates; python3 ta_b50_diagnose.py templates --run-type 2p24GeV_H2           # ~70 s each
python3 ta_b50_diagnose.py gains; python3 ta_b50_diagnose.py select; python3 ta_b50_diagnose.py selvar
# select variants: select_variants/jobs.txt (name run split options), cat jobs.txt | xargs -P 6 -L 1 select_variants/run.sh  (5.5 min)
F=output/templates_spline_ta_fix_20261028; B=output/templates_spline_best50_20261026; RT=2p24GeV_Ta_X17
python3 spline_template_extract.py build --pulses "$B/$RT/pulses/*.npz" --run-type $RT --pool-order energy --n-max 200 --fit-points direct --shape-sigma 0.012 \
  --lambda-select heldout --lambda-chi2-sigma full --lambda-min 1e2 --lambda-fallback 1e3 --converge-tol 0.01 --converge-df 5e-3 --max-iterations 10 --osc-escalate --no-pdf --out $F/n200/templates_$RT   # 274 s
python3 spline_template_extract.py ampcorr --out $F/n200/ampcorr.json --templates $F/n200/templates_$RT.npz   # one file: no change
python3 spline_template_extract.py check --pulses "$B/$RT/pulses/*.npz" --templates $F/n200/templates_$RT.npz --out $F/n200/check/   # 83 s
# gains: B50 dump job.cmd (percrystal_best50_20261026/dumps/scan_{025800,026138,026230}) with the dumps stem and --mf-templates moved to
#   output/percrystal_ta_fix_20261028/n200/dumps and $F/n200; jobs: percrystal_ta_fix_20261028/logs/dump_jobs_n200.txt, runner logs/run1.sh, xargs -P 8 (25 jobs, ~8 min)
G=output/percrystal_ta_fix_20261028; for R in 025800 026138 026230; do python3 percrystal_calib.py fit --dumps "$G/n200/dumps/scan_$R/split_*/result_dump.npz" --splits all --fixed-windows --out $G/n200/$R/; done
# scan: R3 job.cmd (r3_splits) with stem, Ta --mf-templates and --gain-file substituted: resolution_fit_20261028_ta/logs/scan_jobs_n200.txt, xargs -P 8 (25 jobs, ~10 min)
S=output/resolution_fit_20261028_ta
python3 compare_scan_dirs.py --scans R3=output/resolution_fit_20261026/r3_splits R4=output/resolution_fit_20261027_best50/r4_splits Fix_n200=$S/n200/splits --gauss-window 2sigma --out $S/compare/n200_vs_r3
python3 ta_b50_diagnose.py runs; python3 ta_b50_diagnose.py modules --cand $S/n200/splits --tag n200
```

`phase50` is the same with `--pool-order energy-phase --n-max 50`. One job (phase50 25800.00500) hung in plotting after writing its tables; it was killed and rerun (identical `result.tsv`), the hung directory is kept as `split_00500.hung_killed`.

## 4. Promotion

**Decision rule and outcome.** User rule (10-28): if the fixed Ta bank is within 1 SE of R3 on all three Ta runs, promote best-50 H2/C plus the fixed Ta bank as R4. The n200 bank is within 1 SE of R3 (legacy CB, paired over the splits of each run) on 25800, 26138 and 26230 (§3), so R4 is promoted.

**R4.**
- 3.49 / 2.24 / 0.73 GeV H2 and 2.24 GeV C: banks `output/templates_spline_best50_20261026/<runtype>/`, gains `output/percrystal_best50_20261026/<run>/gains_<run>_all.tsv`, splits `output/resolution_fit_20261027_best50/r4_splits` (each R3 `job.cmd` with only the stem, `--mf-templates`, `--gain-file` changed).
- Ta X17 (25800, 26138, 26230): the best-50 Ta pulses rebuilt with `--n-max 200` (`output/templates_spline_ta_fix_20261028/n200/templates_2p24GeV_Ta_X17.npz`), gains `output/percrystal_ta_fix_20261028/n200/<run>/gains_<run>_all.tsv`, splits `output/resolution_fit_20261028_ta/n200/splits`.
- Everything else as R3: position maps per run type (hybrid at 3.49 / 2.24 H2 / C, bilinear at 0.73 and Ta), W = 16 ns, centred rectangle, ±2σ Gaussian, 118 splits.

**Headline** (σ_E/E %, preferred fit, mean over splits; Δ = R4 − R3 paired per split, mean ± SE in pp; Ta per run, never pooled). Sources: `output/resolution_fit_20261027_best50/compare/r4_vs_r3.{md,tsv}`, `output/resolution_fit_20261027_best50/compare_vs_r3.md` (H2, C); `output/resolution_fit_20261028_ta/compare/n200_vs_r3.{md,tsv}`, `n200_vs_r3_splits.tsv`, arm `Fix_n200` (Ta).

| Run type | fit | splits | R3 legacy / gated | R4 legacy / gated | Δ legacy | Δ gated |
|---|---|---:|---|---|---|---|
| 3.49 GeV H2 | CB | 27 | 1.492 / 1.491 | 1.494 / 1.493 | +0.002 ± 0.005 | +0.003 ± 0.004 |
| 2.24 GeV H2 | Gauss ±2σ | 42 | 1.854 / 1.853 | 1.855 / 1.854 | +0.001 ± 0.002 | +0.001 ± 0.002 |
| 0.73 GeV H2 | Gauss ±2σ | 16 | 3.718 / 3.714 | 3.714 / 3.711 | −0.004 ± 0.011 | −0.003 ± 0.010 |
| 2.24 GeV C | CB | 8 | 1.797 / 1.797 | 1.789 / 1.786 | −0.008 ± 0.007 | −0.011 ± 0.008 |
| Ta 25800 | CB | 8 | 1.986 / 1.938 | 1.983 / 1.909 | −0.003 ± 0.009 | −0.030 ± 0.021 |
| Ta 26138 | CB | 4 | 1.950 / 1.920 | 1.949 / 1.890 | −0.001 ± 0.006 | −0.029 ± 0.028 |
| Ta 26230 | CB | 13 | 1.983 / 1.828 | 1.981 / 1.889 | −0.002 ± 0.031 | +0.061 ± 0.038 |

R3 and R4 are equal on e-p (legacy Δ ≤ 0.008 pp, within 1.2 SE at every run type). The Ta legacy Δ is within 1 SE per run; gated 26230 stays +0.061 ± 0.038 (1.6 SE).

**Module medians** (per-module σ_M/μ_M post, %, legacy, `--group-runs`, N > 1000; modules in brackets), R3 → R4: `output/resolution_fit_20261027_best50/module_resolution/{r3_ref,r4}/summary/types_summary.tsv`; Ta `output/resolution_fit_20261028_ta/module_resolution/n200/summary/types_summary.tsv` (new today: `plot_module_resolution.py --scan-dir output/resolution_fit_20261028_ta/n200/splits --out output/resolution_fit_20261028_ta/module_resolution/n200 --group-runs`, 2 s).

| Run type | R3 → R4 | Δ (pp) |
|---|---|---:|
| 3.49 GeV H2 | 1.519 → 1.524 (111 → 111) | +0.005 |
| 2.24 GeV H2 | 1.842 → 1.845 (122 → 123) | +0.003 |
| 0.73 GeV H2 | 3.497 → 3.502 (98 → 98) | +0.006 |
| 2.24 GeV C | 1.810 → 1.797 (78 → 78) | −0.013 |
| 2.24 GeV Ta (pooled over runs, 8 modules; no per-run fit) | 1.908 → 1.929 (8 → 8) | +0.021 |

The Ta median with the 50-pulse bank was 2.001 (+0.093 pp); the n200 bank gives 1.929 (pooled over runs, 8 modules; no SE, no per-run fit).

**Decomposition** (`output/resolution_fit_20261027_best50/out_r4.txt`, H2 only; the Ta per-run points are excluded, so the Ta fix does not enter): CB for all, 2-term b 0.68 → 0.70 ± 0.03 %, a 2.49 → 2.48 %·√GeV, χ² 3.4 → 3.6; the 3-term fit stays unphysical (c² = −0.40). Preferred (±2σ) 2-term b −0.98 → −1.05 ± 0.02 %.

**Caveats.** As §3 Limits; in-sample gains (as R3); the hybrid and bilinear maps are R3's, fitted on R3 energies and not refitted; timing/pile-up were measured with the best-50 banks on 10-26 (summary §3.10, 20261026 §9), not with the Ta n200 bank. `output/resolution_fit_20261027_best50/compare_vs_r3.md` carries the 10-27 verdict "R3 stays" (written before the Ta fix).

**Edits.** Summary §1.4, §3.1, §3.6, §4, §5.1 (§3.9 unchanged); plan S2–S4 ticked with output paths (`gains_vs_spl.md`, `pregain_peak.tsv`, `modres.md` were not written under those names; the plan lines say where the content is).


## 5. U5 real-data two-pulse fit

- What: plan `plans/plan_ml_pileup_20261027.md` U5. The U0 production two-pulse fit (`twopulse_fit.py`, eps 0.02, search 80 ns, best-50 bank of the run type) on the raw waveforms of 120 k PbWO4 channels with an MF pulse 1 per pileup_rates split (25 splits, 3.00 M fits, 13-68 k events per run type; EVIO re-decoded with the process_split selection, event numbers asserted equal to the dump). Script `twopulse_realdata.py` (`decode`, `fit`, `report`; 63 s + 16.8 min with 10 workers + 45 s); outputs `output/ml_pileup_20261027/U5/` (`report.md`, `threshold_scan.tsv`, `mf_agreement.tsv`, `rf_phase.tsv`, `operating_points.tsv`, `fake_by_amplitude.tsv`, `dt_ratio_<group>.tsv`, `u5_realdata.pdf`).
- Fake rate: the U0 thresholds (107 / 98 / 175) do not transfer. Detections with |dt| < 16 ns are 3.4 / 4.3 / 1.1 / 4.5 % of fired channels (0.73 / 2.24 / 3.49 H2, 2.24 C) against 0.015-0.045 % genuine expected from the MF density at 24-60 / 16-60 ns (70-170x); Ta 4.1 / 3.7 / 4.1 % (25800 / 26138 / 26230) against 1.1 / 1.1 / 1.7 %. The excess sits at dt = +8 ... +50 ns after pulse 1 (tail ripple), 6-7 % of channels at 20-100 ADC. The U0 pool is cut-selected (pedestal rms, tail, no MF secondary), the real sample is not.
- Operating points (lowest thr with detections < 2 x exp_lo AND symmetric fake < exp_lo, i.e. fake < MF-expected real rate by both fake estimators; the 0.73 H2 OP sample is still ~60 % fake by the symmetric estimate (sym fake 0.013 % > sym genuine 0.008 %)): e-p 474-632 (0.73 H2 632, 2.24 H2 474, 3.49 H2 548; C 500), Ta 173 / 150 / 130. Efficiency at 8 ns, ratio 1 / 0.3: 98.5-98.9 / 22-46 % (e-p, C), 99.8-99.9 / 90-94 % (Ta); at 12 ns 99.95 / 87-95 % and 99.95 / 99 %. For A1 >= 100 ADC the e-p thresholds rise to 730-2000.
- MF agreement (24-80 ns, ratio > 0.3, +-4 ns): at the production thr 71-73 % of the fit detections are in the MF (Ta 79-80 %), 94-98 % of the MF extras are found (Ta 83-88 %); 0.73 H2 only 46 / 45 %.
- RF: fits cannot reach > 260 ns; comb == grid at |dt| <= 80. e-p shows no comb excess (23-25 %; 3.49 H2 ~3.5 sigma below random, unexplained); Ta on-comb 29 % (|dt| < 16) and 32-33 % (16-80) at the OP vs 36.0 % for the MF-complete reference, i.e. f_flat 0.60-0.66 (rate-based fake share 0.36-0.48; f_flat(< 16) is biased high, f_ref from 24-80 ns pairs, upper bound) and 0.3. Secondary phase mod 4 ns follows pulse 1. This does not separate RF correlation from grid locking.
- Close-pair fraction per fired channel, |dt| < 16 ns, fake-subtracted (pre-pulse symmetry): e-p 0.008 / 0.033 / 0.017 %, C 0.032 %; Ta 1.21 / 1.31 / 1.87 % (MF extrapolation 1.06-1.70 %); the production MF sees 1/5-1/8 (e-p) and 1/10-1/15 (Ta) of it. Fraction of real pile-up in |dt| < 16 = close (fake-subtracted, OP) / (close + MF extra pulses at 16-60 ns, same channels): 12 / 19 / 24 % (0.73 / 2.24 / 3.49 H2), C 21 %, Ta 23 / 25 / 23 % (detected, not efficiency-corrected; sign of net bias not determined).
- Caveats: expectation ignores ratio < 0.3 and MF inefficiency at 24-30 ns; efficiency is from the U0 clean pool at absolute dchi2 (Ta uses 2.24 H2 cells); e-p flat fractions have n = 103-309. Open: the real-data dchi2 needs an amplitude-dependent threshold or a better shape model before the fit can replace the MF.


## 7. LMS gain variation, X17 runs

- Question (user, 10-28): is the LMS gain variation within a run very different between the three X17 runs; plot it every 3 min. Script `lms_gain_variation.py` (numpy/scipy/matplotlib; 50 s with 6 jobs), output `output/lms_gain_20261028/`.
- Coverage caveat. Only every 100th split of each run is on disk (25800: 8 splits, 26138: 4, 26230: 13; 25320: 8 as e-p reference), each about 5-7 s of beam and about 10 min apart. The grid is 3 min, but every populated bin is one split (a split is never cut by the grid), so the plot has one point per ~10 min, 4-13 points per run. True 3-min resolution needs all splits (`replay_gainCorr`, ifarm; 20261021 §6).
- Quantity (as the production correction, `analysis/src/GainCorrCompute.cpp`): WaveAnalyzer `peak_integral`, npeaks == 1 and integral > 0, PbWO4 channels in LMS events (trigger bit 24); LMS1-3 (crate 3 slot 13) in LMS and alpha (bit 25) events. Bin value = median (production: histogram peak fits). `norm` = W · mean over LMS2, LMS3 of 1/R_j, R_j = median LMS_j(LMS events) / median LMS_j(alpha events) per bin (the form of `load_lms_gain_corr`; not W / mean R_j). Relative gain = bin / first populated bin − 1.
- Errors: 1.2533 · 1.4826 MAD / √n per module and bin; event bootstrap (200) of the module-median curve with the first bin included, and a χ² against flat with the bootstrap covariance; split-half (random halves of the events per bin) for the module-to-module noise of the drift, after removing the common mode.
- Time: TI timestamp × 4 ns (250 MHz), minutes since run start. Check: event_number vs time is linear and monotonic per run (R² 1.0000 / 0.980 / 0.9997 for 25800 / 26138 / 26230); the mean rate from run start to the first split is within 0.6-1.1 × of the rate inside that split (25800 started at lower current), so t = 0 is the run start; the LMS pulser runs at 8.5-9.1 Hz in this time.
- Amplitudes are unimodal: one mode in event / bin-median amplitude in all runs, 0 % of the modules multimodal, robust σ 2.0 % per event (one LMS intensity setting).

| Run | Splits → bins | LMS / alpha events | raw median / p90 \|drift\| (%) | raw n > 1 % / > 2 % (noise only) | stat. error of a drift (%) | raw module-median max \|dev\| (%) | norm module-median rms (%) |
|---|---|---|---|---|---|---|---|
| 25800 | 8 → 8, 20-109 min | 734 / 1353 | 0.31 / 0.79 | 48 / 3 (24 / 0.0) | 0.43 | 0.074 | 0.88 |
| 26138 | 4 → 4, 11-52 min | 220 / 405 | 0.34 / 0.86 | 65 / 6 (46 / 0.0) | 0.48 | 0.034 | 0.79 |
| 26230 | 13 → 13, 9-117 min | 564 / 1036 | 0.38 / 1.00 | 116 / 3 (91 / 0.5) | 0.55 | 0.039 | 0.84 |
| 25320 (e-p) | 8 → 8, 7-63 min | 563 / 990 | 0.30 / 0.75 | 45 / 0 (21 / 0.0) | 0.42 | 0.050 | 0.58 |

- Answer. No, the three X17 runs are not very different: the LMS gain is stable within each run and the runs are indistinguishable. The module-median raw gain stays within 0.07 % of the first bin at every bin (last bin −0.02 / +0.03 / −0.03 %); 26138 and 26230 are flat within noise (χ² p 0.57 / 0.94), 25800 has a 0.03 % rms wiggle with no trend (p ≈ 1e-5). Per-module drifts have median |drift| 0.3-0.4 %, equal to their statistical error (0.4-0.6 %): robust scatter 0.44 / 0.50 / 0.57 % against a split-half noise of 0.43 / 0.49 / 0.57 %, so the module-specific part is 0.09 / 0.12 / 0.06 % (consistent with zero). The ordering 25800 < 26138 < 26230 in median |drift| follows the LMS events per bin (median 67 / 55 / 44), not the gain. The n > 1 % counts are 1-2 × the noise-only expectation. Common (light source / HV) drift ≤ 0.07 %; no module-specific (PMT gain) drift above ~0.2 %.
- Reference. The normalised curves move by 0.6-0.9 % rms (bootstrap error 1-1.7 % per bin, χ² p 0.10-0.53, flat within noise): the R_j factor from 39-69 LMS and 75-130 alpha events per bin is noisier than any gain change, so the reference does not absorb a drift (there is none) and adds noise at this statistics. The first-bin normalisation turns the first-bin R_j noise into a common offset (25800: +1.5 % module-median at the last bin). The reference is usable only with production batches (4000 LMS events).
- Groups (crate 3 slots 15 / 8 / 18, central 26 modules; `lms_gain_groups.tsv`). Raw last-bin excess against all modules: central 26: −0.15 ± 0.11 (25800), +0.42 ± 0.12 (26138), −0.42 ± 0.14 % (26230); slots 15 / 8 / 18 at 26230: −0.33 ± 0.25, −0.46 ± 0.18, −0.28 ± 0.18 %. Up to 3.5 σ with the noise taken as independent per module (a contiguous block can be correlated), and the sign of the central block flips between 26138 and 26230: no repeatable gain shift tied to the 26230 timing anomaly (20261026 §3). The 12 modules at the beam hole have a larger within-run RMS (0.56-0.66 % vs 0.27-0.36 % outside |x|,|y| < 110 mm in the X17 runs).
- Cross-check, 25320 against `database/gain_factor/gain_correction/prad_025320_gain_corr.root` (9 batches × 4000 LMS events, read-only): R_j mean ours 0.8716 / 2.3546 / 0.4630 vs production 0.8697 / 2.3602 / 0.4641; module-median peak-to-peak raw 0.09 % (ours, 8 bins) vs 0.04 % (production fitted LMS mean), normalised 2.0 % (ours, reference noise) vs 0.23 % (production gain_W). Production has no file for the X17 runs.
- Benchmark. The LMS-corrected elastic benchmark ran with LMS off for 25800, 26138 and 26230 (no `gain_corr` file: 20261021 §4.1 note ², §6; 20261025 §5 (c), drift bound ≤ 0.06 % for the sampled splits of all runs; summary §3.2). This measurement agrees: the raw LMS gain is flat to ≤ 0.07 % in all three runs, so LMS off costs nothing at the sampled splits. A fast excursion between splits cannot be excluded.
- Files: `lms_gain_variation.pdf` (13 pages: one per run, the three-run overlay, per-run drift and RMS histograms with face maps, per-run group medians), `summary.md`, `lms_gain_runs.tsv` (per-run), `lms_gain_modules.tsv` (per module: drift, error, RMS, raw and norm), `lms_gain_bins.tsv` (per bin: references, R_j, module-median curves and errors), `lms_gain_modules_bins.tsv` (module × bin long table), `lms_gain_groups.tsv`, `run.log`, `cache/` (one npz per split).
- Run: `MPLBACKEND=Agg python3 lms_gain_variation.py --jobs 6` (default output `output/lms_gain_20261028/`; `--stage analyse` re-uses `cache/`; `--runs 25800 26230` for a subset).

**§7 update (round 2).**
- Deliverable is now one plot: `output/lms_gain_20261028/lms_gain_variation.pdf` (1 page) and `.png`; the 25800 > 26138 > 26230 runs are concatenated on one 3-min axis (x = 0 is the start of 25800), relative to the first bin of 25800, with all modules, 16-84 % band, module median (reference-normalised, bootstrap error), raw median (dashed, no reference) and the slots 15/8/18 + central 26 median. The per-run, histogram, face-map and group pages are removed; their numbers stay in the TSVs and `summary.md` (new `lms_gain_xrun_bins.tsv`).
- Run gaps from SYNC/EPICS unix_time (run start = unix − TI × 4 ns, splits agree within 1 s): 25800 starts 2026-07-12 09:33 UTC, 26138 601.2 h later, 26230 141.1 h after 26138 (≤ 599.4 h and ≤ 140.3 h after the last local split; run ends are not on disk). Module-median jump at the run boundaries: 25800 → 26138 raw +0.48 ± 0.03 %, norm +1.18 ± 1.13 %; 26138 → 26230 raw −0.48 ± 0.03 %, norm −0.05 ± 1.27 %.
- Noise-only counts are now per module, Σ_m erfc(thr / (√2 · err_m)): raw n > 1 % / > 2 % expected 46 / 2.4, 71 / 4.3, 113 / 8.1 (25800 / 26138 / 26230) against 48 / 3, 65 / 6, 116 / 3 observed; they match (pull rms 1.05 / 1.02 / 1.03). The earlier 24 / 46 / 91 and "1-2 ×" above used one global σ and are superseded.
- Group errors now use the group's own per-module errors: central 26 excess +0.42 ± 0.25 % (26138, 1.7σ) and −0.42 ± 0.27 % (26230, 1.6σ); the largest group deviation is 1.9σ, so the "up to 3.5σ" above is superseded.
- Typos in the bullets above: "median 67 / 55 / 44" is 68 / 54.5 / 43 LMS events per bin; "75-130 alpha events" is 70-130; "module-specific part ≤ 0.1 %" is ≤ 0.12 %.

**§7 update (round 3).**
- "Indistinguishable" in §7 holds for the within-run behaviour only. Between runs the groups shift (raw median per-module rel at the first bin of the run vs bin 1 of 25800, %; rest = outside slots 15/8/18 and central 26): 26138 all +0.5, slot 15 / 8 / 18 −14.5 / −12.5 / −13.6, central 26 −14.4, rest +0.5, 100 modules with |rel| > 5 %; 26230 all +0.1, slot 15 +110.9, slot 8 +0.6, slot 18 −5.1, central 26 +89.5, rest +0.1, 124 modules > 5 % (norm values are within 3 % of these for the groups). Consistent with an HV/gain change between runs; not established here. Table `lms_gain_xrun_groups.tsv` and `summary.md`; the plot marks the off-scale group median (±8 % axis).
- The quoted boundary jumps (+0.48 / −0.48 % raw, +1.18 / −0.05 % norm) are medians of the per-module jumps, not the jump of the module median; the plot, summary and this note now say "median per-module jump".
- Corrections to the round-2 block: central 26 at 26230 is 1.5σ (not 1.6σ). Superseded §7 body values: slot 15 / 8 / 18 errors ±0.25 / 0.18 / 0.18 → 0.39 / 0.29 / 0.31; central 25800 ±0.11 → 0.21; 25320 noise expectation (21 / 0.0) → (40 / 1.6).
- Fix to the first bullet: the norm group levels differ from the raw ones by 1.9-5.6 points, not "within 3 %" (26138 slot 15 −12.6 norm vs −14.5 raw; 26230 slot 15 +116.5 vs +110.9, central 26 +94.6 vs +89.5); see `lms_gain_xrun_groups.tsv`.

## 6. Baseline noise diagnostic
`noise_diag.py` -> `output/noise_diag_20261028/` (`noise_diag.pdf`, `noise_*.tsv`, `numbers.json`). The HyCal readout is zero-suppressed (16 e-p / 43 Ta PbWO4 channels per event, each with a pulse), so the baseline is samples 0..29 of every present channel whose earliest MF pulse starts >= 132 ns (99 % e-p, 75-79 % Ta), 59.4-59.9k events for each of 0.73 / 2.24 / 3.49 GeV H2 (24590/10, 25228/50, 24503/50) and Ta 26230/100, 25800/100, plus 48-sample late-only windows and 1087 LMS events (complete boards). e-p numbers are primary; Ta windows carry pile-up tails.
- **0.40 vs 0.92 is definitional.** e-p per-module raw std 0.96-0.98 (pooled ddof 1 1.01-1.04; raw/stored 2.20 / 2.54 / 3.08 at 5 / 50 / 95 %). The stored sig_ped is the sigma-clipped rms of the 3-tap triangular smoothed buffer (replica 0.407-0.410 vs stored 0.403; raw 1.4826 MAD is 1.483 in 70 % of windows, integer quantisation). The raw acf is lag1 -0.35 with an alternating tail to lag 29; the kernel then predicts 0.460 smoothed vs 0.452 measured (white 0.564).
- **Not white.** 28 % of the variance sits in the Nyquist bin (125 MHz = fs/2; alternating amplitude rms 0.52 ADC, random sign, kurtosis 2.4), plus a weak line at 62.5 MHz (fs/4, +57 % over neighbours); the other bins vary +-34 %; no excess within 8-125 MHz (4.2 % below 10 MHz, white 6.9 %; 88.6 % above 30 MHz, white 79 %), but see the per-event offset below. RF aliases (249.52 MHz -> 0.48 MHz) fall in the first bin; no RF-phase modulation (<= 0.009 ADC vs 0.52).
- **Coherent.** Same-sample corr of baseline windows (pairs with >= 300 common samples): within board 0.59, same crate other board 0.55, other crate 0.34 (alternating part removed: 0.46 / 0.41 / 0.16); the alternating amplitudes of one event correlate 0.88 across channels. LOO board-mean subtraction: sigma 1.00 -> 0.71, common-variance fraction 0.60 (0.47 after the tone). Budget: tone 27 %, other common (blue, 0.14 ADC^2/bin at 8 MHz to 0.65 at 60-117 MHz) ~33 %, independent ~40 % (0.63 ADC). Per-event window mean wanders 0.30 ADC (white 0.19); mean-removed windows cannot constrain it, so the excess is attributed to a per-event offset (~9 % of the variance).
- **No gain from a linear fix.** The board mean cannot be formed at pulse samples (zero-suppression); subtracting the alternating component (own pre-pulse fit, or the mean of the event's other channels) from real pulses leaves sigma_t vs RF (E >= 100 MeV) at 0.1445 / 0.1460 / 0.1459 ns (2.24 GeV; raw / own / common; 0.73 and 3.49 GeV within 1 %) and single-pulse dchi2 median / 97th 31.4 / 271 -> 31.0 / 271 (-1.2 % / -0.3 %; the tone-only test is insensitive by construction, the dchi2 evidence is U0's noise-only 0.14 / 16 vs real singles 8.0 / 96). Noise-only onset error at 301 ADC is 7-17 ps against 145-185 ps measured.
- **Verdict.** Structured and non-diagonal, yet timing and dchi2 are not noise limited, so ML denoising of the baseline has no target at E >= 100 MeV (only below ~60 ADC, where sigma_p^2 exceeds (0.02 y)^2).
- **sigma definition.** For chi2 on raw samples use the raw baseline std (~1.0 ADC per module) with the measured acf; there is no single diagonal equivalent (0.6-1.3 ADC: 30-sample windows do not fix the acf beyond lag 29 or the offset/pedestal handling). sig_ped (0.40) is right for WaveAnalyzer thresholds only; do not regenerate the table, add a `sig_raw` column instead.

## 8. Clustering parameter scan
`output/cluster_scan_20261028/` (`report.md`, `tables/*.tsv`, `cluster_scan_20261028.pdf`). The HyCal island-clusterer parameters come from `database/reconstruction_config.json` ("hycal"); they were overridden at run time from copied configs by `tools/cluster_cfg_wrap.py`, which wraps `_common.setup_pipeline` (no edit of database/, C++ or the benchmark). R4 values: seed (min_center_energy) 10 MeV, neighbour (min_module_energy) 1.0 MeV, min_cluster_size 1, corner_conn false, W 16 ns. There is no radius knob: the cluster energy is the island sum, and corner_conn sets the island extent. The baseline through the wrapper reproduces R4 bit-identically on the 20 reference splits.
- **One-at-a-time** (4 even-rank splits per run type; 60k events per e-p split vs a 60k baseline; Ta and 0.73 GeV seed points at full stat). Neighbour 0.5 / 2 MeV are identical to R4 (the MF pulse finder is the real floor, ~5 MeV), 5 MeV is null, and 10 MeV is worse (2.24 C +0.033 ± 0.013 pp). corner_conn 1 and min size 2 are null. W 8 / 32 change e-p by ≤ 0.012 pp (W 32: Ta gated +0.080 ± 0.025). The seed threshold gains monotonically up to a plateau at 30–40 MeV (3.49 H2 −0.046 ± 0.010, 2.24 H2 −0.036 ± 0.009, 2.24 C −0.037 ± 0.011 pp at 30 MeV), while 5 MeV gives no gain (mixed sign, < 2 SE). The seed × W grid shows W adds nothing on top of the seed change.
- **Validated** (seed 30 MeV chosen on the even subset; 60 odd-rank splits, full stat, vs R4, legacy / gated): 3.49 H2 −0.059 ± 0.010 / −0.059 ± 0.010; 2.24 H2 −0.046 ± 0.003 / −0.046 ± 0.003; 2.24 C −0.033 ± 0.005 / −0.033 ± 0.006; 0.73 H2 −0.040 ± 0.011 / −0.040 ± 0.010 pp. The single-cluster count changes by ≤ 0.03 % and the inner-ring count by −0.1 to +0.6 %. Events move from both tails into the peak (`tables/tail_migration.tsv`: core +0.3–0.7 %, high tail −3 to −19 %). The benchmark's peak_yield (−0.04 to −0.49 %) uses each configuration's own σ, so it is not an acceptance measure. For Ta, Δσ is inconclusive (25800 legacy −0.019 ± 0.007; 26138 / 26230 consistent with 0 at SE 0.1–0.19), and the high tail grows by +16 / +2 / +4 events (+6–8 %).
- **Proposal for R5 (not promoted):** min_center_energy 30 MeV, everything else as R4. Three conditions apply before promotion. (1) Re-derive the percrystal gains and the position correction at seed 30, since both were solved at seed 10. (2) Check secondary-shower merging for the Møller / two-cluster / X17 analyses; the Ta high tail is a weak warning sign. (3) Rerun the 118 splits for the headline. The summary reference is unchanged.

## 9. Merge of origin/main (6e638c6 → 2daedc9)
`output/merge_check_20261028/` (`report.md`, `tables/{summary,splits}.tsv`, `tools/merge_{driver,cfg_wrap,analyze}.py`, `logs/driver.log`: 123 jobs, 0 failed, 3257 s). HEAD 2daedc9 merges origin/main 6e638c6 into 1cd8e7d; the benchmark engine is unchanged by the merge. Jobs use the R4 `job.cmd` flags; only the database paths below are rewritten. Pre-merge files are `git show 1cd8e7d:` extracts in scratch; nothing was written to `database/`.
- **Build.** `cmake --build build --target prad2py -j 14` on the existing `build/` (Release, BUILD_PYTHON=ON, python 3.13): exit 0, `build/python/prad2py.cpython-313-darwin.so`, `dec.WaveConfig().clk_ns` = 4.0, no shim. A stale `prad2py.cpython-314-darwin.so` (Sep 14) stays in `build/python`; python 3.13 does not load it.
- **Variants.** ri: pre-merge "hycal" block and pre-merge `runinfo/general.json` (isolates the code). (a) R4-config: pre-merge "hycal" values patched in-process (leakage off, log_weight_thres 3.6, energy_bias off), merged runinfo. (b) main-config: merged `reconstruction_config.json` and runinfo, with the pre-merge LMS files and R4's explicit `--hc-calib-file`.
- **R4 reproduced bit-identically.** ri gives `result.tsv` byte-identical to R4 on all 41 splits (16 e-p/C reference splits, 4 per run type, and the 25 Ta splits). (a) moves 2 of 164 split × path × fit entries by > 0.005 pp (24590.00080 legacy Gauss +0.0067, 24590.00040 gated CB +0.0056), through main's runinfo, not the code: HyCal z 6267.0 → 6265.887 mm, W1019 and W1020 added to the dead list (578 → 580); runs ≥ 25308 also get a beam-target offset (0.805, 1.486) mm. Mean Δ(a) ≤ 0.0022 pp on e-p; Ta (a) equals R4.
- **Main's config is worse; only log_weight reaches the benchmark.** `_make_cluster` builds a fresh `HyCalClusterConfig()` and copies 8 fields; of main's "hycal" changes only log_weight_thres 5.2 gets through. leakage_correction and energy_bias_correction stay false (C++ default) in every variant. log_weight changes only (x, y), so it acts through pos-corr, the inner radius and the exclusions: (b) +0.014 to +0.023 pp at 3.49 / 2.24 GeV H2 and C (2–3.5 SE), +0.033 pp at 0.73 GeV (1.7 SE); Ta mixed within errors.

σ_E/E (%), mean over splits; Δ vs R4 paired per split (pp, mean ± SE). e-p/C: the 4 reference splits per run type, R4 = their stored `result.tsv` (not the 118-split headline). Ta per run, never pooled, all n200 splits:

| Run type (fit) | n | R4 legacy / gated | Δ(a) legacy / gated | (b) legacy / gated | Δ(b) legacy / gated |
|---|---:|---|---|---|---|
| 3.49 GeV H2 (CB) | 4 | 1.507 / 1.507 | +0.0002 / +0.0001 | 1.531 / 1.530 | +0.023 ± 0.007 / +0.023 ± 0.006 |
| 2.24 GeV H2 (Gauss) | 4 | 1.826 / 1.826 | −0.0004 / +0.0002 | 1.840 / 1.842 | +0.014 ± 0.006 / +0.015 ± 0.006 |
| 0.73 GeV H2 (Gauss) | 4 | 3.673 / 3.671 | +0.0022 / +0.0015 | 3.706 / 3.705 | +0.033 ± 0.020 / +0.034 ± 0.019 |
| 2.24 GeV C (CB) | 4 | 1.798 / 1.797 | −0.0007 / −0.0005 | 1.816 / 1.815 | +0.018 ± 0.009 / +0.018 ± 0.008 |
| Ta 25800 (CB) | 8 | 1.983 / 1.909 | 0.0000 / 0.0000 | 2.000 / 1.933 | +0.018 ± 0.006 / +0.025 ± 0.026 |
| Ta 26138 (CB) | 4 | 1.949 / 1.890 | 0.0000 / 0.0000 | 1.954 / 1.861 | +0.005 ± 0.012 / −0.030 ± 0.030 |
| Ta 26230 (CB) | 13 | 1.981 / 1.889 | 0.0000 / 0.0000 | 1.986 / 1.860 | +0.005 ± 0.043 / −0.030 ± 0.035 |

CB σ of the Gauss-preferred types, legacy / gated, R4 → (a) → (b): 2.24 GeV H2 1.759 / 1.759 → 1.759 / 1.760 → 1.765 / 1.764; 0.73 GeV 2.952 / 2.955 → 2.953 / 2.956 → 2.937 / 2.938.

- **Database paths.** (1) Main moved `database/calibration/{calibration_factor_fadc_3p5, calibration_factor_fadc_0p7, calibration_factor_2p2_fadc}.json` to `calibration/old/` (identical, git R100), so R4's `--hc-calib-file` paths no longer exist; the driver passes `calibration/old/<name>`. The Ta `x17_*_fadc.json` are unchanged. (2) LMS blobs: `--lms-corr auto` is used by 25308 and 25320; main deleted `gain_correction/prad_025308_gain_corr.root` and modified `prad_025320_gain_corr.root`. All variants use the pre-merge (1cd8e7d) blobs, so (b) does not include main's new 25320 table. The 71 Physics_calib deletions and the other gain_correction changes are not referenced by our jobs or `_common.py`. (3) Main's runinfo assigns `PRad2_2p2_round1_*` calibrations from run 25200 on; our jobs pass `--hc-calib-file`, so these defaults are not used.
- **max_events fix.** The merged `_common.setup_pipeline` no longer accepts `max_events`; the benchmark passed it (line 2140, unused before the merge) and raised a TypeError. Fixed the same day by removing `max_events=args.max_events` from that call (uncommitted at the time of the check).
- **Verdict.** Push of `waveform_jz` cleared: the merged code is functionally equivalent for this pipeline. Running on main without the wrapper needs the `max_events` fix and the `calibration/old/` paths. The merged `reconstruction_config.json` carries log_weight_thres 5.2 (R4: 3.6); adopting it is a physics decision, open.

## 10. Møller electrons in the per-module spectra

**(i) Diagnostic.** `output/twopeak_diag_20261028/` (`report.md`, `drift_cases.tsv`, `iteration_trace.tsv`, `scan_*.tsv`, `figures/fig_A_*`, `fig_B_*`). Question: two peaks per module in one split, gain drift or refinement over-iteration? Cases (2.24 GeV H2; 2-Gauss + flat fit to the seed-module raw cluster energy): C1–C3 25320.00100 / 00200 / 00300 (W460, W565, W592; `scan_20261025_ref`, LMS on), C4 25228.00050 (W459, R4, LMS off), C5 25320.00050 (W593, 10/20 bimodality scan, old engine).
- **Not gain drift.** LMS within the split −1.3 ± 1.3, +2.3 ± 1.9, −0.7 ± 1.8, −1.1 ± 1.6 % (C1–C4; C4 with 25320 as proxy), split to split ≤ 1.6 %; low-peak fraction flat in time (χ² p 0.26–0.96); both peaks in every time bin.
- **Not over-iteration.** Every refinement scales a seed module's whole spectrum by one gain: cal/raw is constant per seed module (relative spread < 1e-6) in all stored dumps. Of 553 (R4, 2.24 GeV H2) and 417 (10/20 scan) bimodal module × split spectra, 0 become unimodal after refinement. Elastic-window gain: pass 1 1.001 (C1) / 0.999 (C3), passes 2–5 exactly 1.000.
- **Low peak = single-arm Møller electrons.** High peak 0.997–1.000 × E'ep, low 0.997–1.010 × E_ee(cluster r); measured ratio 1.25 vs predicted 1.26–1.27 (C1–C4). The low population sits 4.0–5.7 mm toward the beam (high 1.6–1.9 mm; KS p < 1e-10). Near the beam hole E_ee is within 25 % of E'ep. The same 15-module ring is bimodal in every 2.24 GeV H2 split; 3.49 GeV H2 has no bimodal module; Ta was not searched. C1: 0 of 694 low-peak events inside the elastic window [2142, 2324] MeV (Møller gap 11.4 / 12.9 σ).
- **Old flip artefact (C5).** The retired cb-peak took the highest 20-MeV bin; with two near-equal peaks the gain came out 0.979 (high) or 1.132 (low) by noise (1.132 / 0.979 = 1.156 = the separation), and splits 00200 / 00250 alternated with each extra iteration (ghost at ~2590 MeV, not a merge).
- **The "merge" is in the pooled spectrum only.** The post pool drops uncalibrated modules (cal = NaN, including the inner `--exclude-rowcol` box; C1: 78 % of the inner events, 98k in 32 box modules). Pooled pre is bimodal (1892 / 2221 MeV), post unimodal; the calibrated low-shoulder fraction is 0.188 pre and post.
- Caveats: `event_num` is the running accepted-event index; 25228 has no LMS; partner test impossible (one cluster > 500 MeV per event); low-population slope 30–60 % of the Møller prediction.

**(ii) Rejection and two-peak fit.** `output/moller_reject_20261028/` (`report.md`, `tables/`, `figures/q_kin_0p73.png`, `cbfit_0p73_compare.png`). R4 settings, pre-merge configuration (§9 ri). New benchmark options, uncommitted, default off: `--moller-reject {off,pair,kin,any}`, `--moller-reject-gains {on,off}`, `--moller-kin-n` (3), `--moller-kin-m` (2), `--cb-two-peak`, `--cb-sigma-tie [FRAC]`; `plot_module_resolution.py` drops `moller_rej==1` rows. Flag off = R4: max|Δ| 0.00000 pp over 164 entries (41 reference splits).
- **Quantification, 0.73 GeV.** Modules with σ/μ > 2× median (3.50 %) are exactly a 15-module ring at r 79–98 mm outside the 6×6 box: W524 W457 W462 W695 W700 W490 W497 W660 W667 W425 W426 W566 W591 W600 W427. High peak at E'ep(θ) (+0.5 %), low at E_ee(θ) (618–653 MeV), within 1 % of prediction, 3.5–4.2σ apart; Møller is ~48 % of their events. The ring enters the gain refinement (217 of 240 module × split entries, ~10 % Møller in window), the module median (3.50 % with, 3.45 % without) and the headline pool (33 % of it; Møller ~9 % of the pool, 1.4 % inside ±2σ). 2.24 GeV H2 (16 modules), C (15) and 3.49 GeV (2) (pooled-spectrum count) show a Møller bump 10–11σ below the elastic peak, none above 2× median; Ta none.
- **Criterion, `--moller-reject kin`.** Reject a cluster if |E − E_ee(θ)| < 3σ_E and |E − E'(θ)| > 2σ_E; θ from the centroid, E after gain file and pos-corr, σ_E/E = √((2.5 %/√E'[GeV])² + 0.7 %²). Offline on the ring: ~90 % of the Møller and 2.5 % of the elastic events removed, < 1 % sculpting of the ±2σ window. The pair tag finds a partner for < 0.4 % at 0.73 GeV (partner at r ≈ 550–750 mm, mostly outside HyCal) and tags 12 % at 3.49 GeV with a 10 % tolerance; not used.

σ_E/E (%), legacy (gated on − off agrees within 0.02 pp); e-p/C on the 4 reference splits per run type, Ta all splits, per run (`tables/eval_summary.tsv`):

| Run type | fit | R4 (= off) | on | on, R4 gains | on − off (pp) | rejected (% of pool) |
|---|---|---|---|---|---|---|
| 3.49 GeV H2 | CB | 1.507 | 1.507 | 1.507 | +0.000 ± 0.000 | 1.7 |
| 2.24 GeV H2 | Gauss | 1.826 | 1.826 | 1.826 | −0.000 ± 0.001 | 13.1 |
| 2.24 GeV H2 | CB | 1.759 | 1.737 | 1.738 | −0.022 ± 0.005 | 13.1 |
| 0.73 GeV H2 | Gauss | 3.673 | 3.447 | 3.433 | −0.226 ± 0.013 | 28.0 |
| 0.73 GeV H2 | CB | 2.952 | 3.417 | 3.400 | +0.465 ± 0.029 | 28.0 |
| 2.24 GeV C | CB | 1.798 | 1.797 | 1.796 | −0.001 ± 0.001 | 3.6 |
| Ta 25800 | CB | 1.983 | 1.956 | 1.959 | −0.027 ± 0.014 | 2.3 |
| Ta 26138 | CB | 1.949 | 1.946 | 1.947 | −0.003 ± 0.003 | 1.9 |
| Ta 26230 | CB | 1.981 | 1.978 | 1.978 | −0.003 ± 0.006 | 2.2 |

- **Two-peak CB, 0.73 GeV** (`--cb-two-peak`, all 16 splits, rejection off; `tables/cbfit_0p73_splits.tsv`): legacy Gauss 3.714 %, old CB 2.974 %, two-peak elastic CB 3.265 % (gated the same). χ²/ndf 2.6 against ~50 for the old and the σ-tied CB; no fallback. The σ-tied option (±15 %: 3.13 %) tracks FRAC (3.31 at 0.10, 3.00 at 0.25). Fitted Møller fraction 0.34–0.48 per split (overstates the event fraction). Rejection on with the old single CB: 3.39–3.47 %.
- **Per module, 0.73 GeV, rejection on** (16 splits; `module_resolution_on/`): modules above 2× median 15 → 0; ring maximum σ/μ 11.6 → 3.7 %; fitted modules 98 → 85 (some fall below N > 1000; W524 now fails the Møller-gap check); median 3.447 % (`summary/types_summary.tsv`).
- **R4 at 0.73 GeV is Møller-biased.** Gauss 3.673 (4 reference splits) / 3.714 (16 splits, the headline) and CB 2.952 / 2.974 include the Møller peak; the single CB is a poor fit there (χ²/ndf ~50). Møller-free: ≈ 3.42–3.45 % (CB 3.417 / Gauss 3.447, 4 reference splits; 3.400 / 3.433 with the R4 gains). Cross-check: two-peak elastic CB 3.265 % (16 splits, Møller kept). Re-derived and R4 gains differ by ≤ 0.017 pp (`eval_summary.tsv`; the report text says ≤ 0.015).
- PDFs: two-peak fit `cbfit_0p73/benchmark_hycal_timing_scan_mf_pc_spl_all_calibrated_plots.pdf` (16 pages); rejection on `on/benchmark_hycal_timing_scan_mf_pc_spl_all_calibrated_plots.pdf` (16 pages); per-module spectra `module_resolution_on/0p73_H2/0p73_H2_spectra.pdf` (12 pages).
- Standing rule (2026-10-28): Møller rejection (`--moller-reject kin`, n=3, m=2) is proposed for R5; defaults unchanged pending user promotion. `--cb-two-peak` is the fit to use if Møller events are kept; not needed with rejection on.

## 11. main's energy_bias, leakage and Geant4 cluster profile
`output/energy_bias_20261028/` (`report.md`, `tables/{summary,splits,contrasts}.tsv`, `tools/eb_wrap.py`). Splits and control as §9 (ri = R4). `eb_wrap.py` enables main's corrections in-process: a recon-config copy in scratch with the flag on, plus a carried builder config, which also brings main's Geant4 cluster profile (hence the h / g controls); member energies are rescaled so the `--gain-file` factor is kept. Variants: f pos-corr off; h / g Geant4 profile, no correction, pos-corr off / on; c E-bias, pos-corr off; d E-bias + pos-corr; e leakage + pos-corr.
- **Where the corrections live** (HEAD 2daedc9). E-bias: `HyCalCluster::ReconstructHits` (`prad2det/src/HyCalCluster.cpp:452`), E × 1/(1 + bias[module][5×5 cell of (x, y)]); ep table if E > E_beam(1 − 0.12/√(E_beam/GeV)), else ee; tables `database/energy_bias/{0p7,2p2,3p5}GeV_{ee,ep}.json` (28800 cells each; 3.49 H2 → 3p5, 2.24 C / H2 / Ta → 2p2, 0.73 → 0p7), loaded only by `PipelineBuilder` (`PipelineBuilder.cpp:430`) when `hycal.energy_bias_correction` is true. Leakage: `apply_leakage_correction` (`HyCalCluster.cpp:113, 485`), only seeds flagged kLeakage with ≥ 4 hits. The benchmark's `_make_cluster` uses a fresh `HyCalClusterConfig` (SimpleProfile, no table, no leakage).

Contrasts, legacy, pp, paired per split (gated in `tables/contrasts.tsv`); Ta per run, never pooled:

| Run type (fit, n) | pos-corr removed (f − ri) | E-bias alone (c − h) | E-bias on pos-corr (d − g) | leakage (e − g) | Geant4 profile (g − ri) |
|---|---|---|---|---|---|
| 3.49 GeV H2 (CB, 4) | +0.182 ± 0.009 | −0.108 ± 0.006 | +0.119 ± 0.007 | −0.001 ± 0.001 | −0.047 ± 0.004 |
| 2.24 GeV H2 (Gauss, 4) | +0.066 ± 0.005 | +0.011 ± 0.006 | +0.118 ± 0.008 | +0.000 ± 0.000 | −0.037 ± 0.003 |
| 0.73 GeV H2 (Gauss, 4) | +0.012 ± 0.016 | +0.022 ± 0.021 | +0.033 ± 0.020 | −0.000 ± 0.002 | −0.011 ± 0.011 |
| 2.24 GeV C (CB, 4) | +0.062 ± 0.015 | +0.008 ± 0.005 | +0.126 ± 0.028 | +0.000 ± 0.000 | −0.029 ± 0.002 |
| Ta 25800 (CB, 8) | +0.003 ± 0.012 | −0.015 ± 0.033 | +0.050 ± 0.023 | +0.000 ± 0.001 | −0.024 ± 0.029 |
| Ta 26138 (CB, 4) | +0.036 ± 0.029 | +0.080 ± 0.116 | +0.123 ± 0.098 | −0.002 ± 0.002 | +0.092 ± 0.077 |
| Ta 26230 (CB, 13) | +0.022 ± 0.039 | −0.032 ± 0.067 | +0.058 ± 0.064 | +0.009 ± 0.014 | −0.023 ± 0.093 |

- **energy_bias: partial duplicate of pos-corr; conflicts when stacked.** Both correct the same in-cell response. E-bias is weaker: alone it recovers about 60 % of pos-corr at 3.49 GeV (−0.108 of +0.182 pp) and nothing at 2.24 GeV H2 (+0.011). Stacked on pos-corr it double-corrects (+0.12 pp at 3.49 and 2.24 GeV H2, +0.126 on C). Keep pos-corr; leave `energy_bias_correction` off.
- **Leakage inert.** |Δ| ≤ 0.002 pp on e-p and C; Ta within 1 SE.
- **Geant4 profile.** −0.011 … −0.047 pp on e-p; Ta noise (−0.024 ± 0.029 / +0.092 ± 0.077 / −0.023 ± 0.093). Mechanism unverified: the profile arrives with the carried builder config, not as a separate switch; a clean test needs a `_make_cluster` option that selects the profile.

## 12. Housekeeping
- **Presentation.** `output/presentation_20261028/`: `presentation.md` (collaboration-meeting status talk, 2026-10-28, 28 `##` sections including backup and open items), `figures/` (f01–f11 PNG, `make_presentation_figures.py`, `plotted_values.tsv`). Open items for the speaker (its last section):
  - TBD: main-branch `energy_bias` comparison (no source when the deck was written; now §11).
  - TBD: funding acknowledgment line; talk length and slot (deck sized for 10–15 min).
  - f08 R5-projected bars = R4 118-split mean + odd-split Δ (seed 30 MeV, §8), errors in quadrature, no Ta bar. The R5 table's "σ R4 on subset" is the odd-split mean (2.24 H2 1.866 %), not the headline (1.855 %).
  - The best-resolution table mixes R3 and R4 rows (summary §1.4: the choice is not a result); consider R4 throughout.
  - Pile-up rates are not comparable: summary §3.8 0.15–0.41 % (e-p, E_dom > 50 MeV, after burst veto) vs the deck's 10-27 per-channel 0.17–0.56 % (any extra pulse).
  - f05: one representative split per panel; 0.73 GeV shows a shoulder near 0.6 GeV (E_ee there is 618–653 MeV, §10).
  - Conflicts settled in the deck: cluster regression on hold, not dead (V1 HGB/MLP, credible gain ≤ 0.026 pp, 20261027 §9); Wassim − best-50 at 3.49 GeV is −0.020 ± 0.004 pp on 4 splits, a possible gain that needs the full scan.
- **a/b (corrected).** R4 2-term fit a = 2.48 ± 0.02 %·√GeV, b = 0.70 ± 0.03 % (R3: 2.49, 0.68); the deck and f06 use these (`compare_vs_r3.md`, summary §3.6).
- **LMS between-run shift (corrected).** At 26230 only slot 15 (+110.9 %) and central 26 (+89.5 %) jump; slot 8 +0.6 %, slot 18 −5.1 %. At 26138 all four groups are −12.5 … −14.5 %. §7 round 3 already has these values; the wrong version ("+90 to +111 % for slots 15/8/18 and central 26") was in the presentation brief only.
- **Regenerated per-split PDFs.** `output/resolution_fit_20261027_best50/r4_splits/benchmark_hycal_timing_scan_mf_pc_spl_all_calibrated_plots.pdf` (118 pages; its Ta splits use the 50-pulse Ta bank, not R4's n200); `output/resolution_fit_20261028_ta/n200/splits/benchmark_hycal_timing_scan_mf_pc_spl_all_calibrated_plots.pdf` (25 pages, R4 Ta).
- **Per-module spectra.** `output/resolution_fit_20261027_best50/module_resolution/r4/<runtype>/<runtype>_spectra.pdf`: 3p49_H2 14, 2p24_H2 14, 0p73_H2 12, 2p24_C 16 pages (files dated 10-27 09:16). Ta per run with the n200 bank: `output/resolution_fit_20261028_ta/module_resolution/n200/{25800,26138,26230}/<run>_spectra.pdf` (4 / 2 / 2 pages); the `Ta_X17` files in both directories pool the runs.
- **Stray working-tree change.** `notes/figures/methods_20261027/dataflow.mmd` is deleted and `notes/figures/methods_20261027/dataflow.md` is untracked (a Mermaid `flowchart TD`). Origin unknown; left as is for the user.
- **Edits.** Summary header (merge line), §1.4 (0.73 GeV footnote), §3.11–§3.13, §4 (R5 candidates). R4 reference numbers unchanged.
