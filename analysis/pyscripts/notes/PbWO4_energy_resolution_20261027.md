# PbWO4 energy resolution: best-50 template bank on the full R3 scan (not promoted), CFD-fraction scan

**Date:** 2026-10-27
**Branch:** `waveform_jz` (HEAD `3566071`; everything since is uncommitted, as in 20261026 §1.2)
**Plan:** `analysis/pyscripts/plans/plan_best50_pipeline_20261026.md` (S2, S3, S4 module medians; outputs moved to the paths below)
**Previous notes:** `analysis/pyscripts/notes/PbWO4_energy_resolution_20261026.md` (§8 bank, §9 reduced pipeline), `analysis/pyscripts/notes/PbWO4_energy_resolution_summary.md` (§4 R3)

Paths are relative to `analysis/pyscripts/` in `/Users/L00338853/work/PRad/prad2evviewer` unless they start with `prad2dec/` (repo root). σ/E, preferred fit, legacy/gated, Δ in pp and the run-type names are as in the summary note; ±2σ Gaussian. R3 = `output/resolution_fit_20261026/r3_splits` (summary §4). B50 = the best-50 bank `output/templates_spline_best50_20261026/` (20261026 §8).

## 1. Status

User decision (10-27): adopt B50 as the production reference, subject to the full rescan. Done today:
- Per-crystal gains with B50 for the remaining 9 runs (75 dump jobs incl. the 24386.00060 extra dump, 0 failures, 44 min at 12 parallel; 9 fits converged in 2 iterations). The 5 runs of 10-26 are reused.
- The 118-split rescan with the R3 configuration, changing only the bank and the gains (118 / 118, 0 failures, 85 min at 12 parallel, alongside another job).
- Comparison to R3 (§2). **Result: same at the four H2/C run types, worse at Ta. R3 stays the reference; the summary note is not changed by §2 (its 10-27 CFD rows come from §5).**

No code changes for §2. Timing was not rerun (§2.5). Also finished today, by a separate executor: the CFD-fraction scan (§5; Python only, `cfd_fraction_scan.py`, outputs `output/cfd_scan_20261027/`).

## 2. Best-50 adoption test

### 2.1 Gains

`output/percrystal_best50_20261026/<run>/gains_<run>_all.{tsv,json}`, all 14 runs. Procedure as for the 5 runs of 10-26: each 10-24 dump `job.cmd` (`output/scan_20261024/spl_percrystal/scan_<run>/split_*/job.cmd`) with the output stem moved to `output/percrystal_best50_20261026/dumps/scan_<run>/split_<n>/` and `--mf-templates` pointed at B50 (the substitution reproduces the 10-26 `job.cmd` of 25800.00100 exactly), then `percrystal_calib.py fit --splits all --fixed-windows`. 24386 also takes `dumps_extra/024386_00060_dump.npz`, as `percrystal_spl` does.

Consistency against `percrystal_spl` (legacy path, modules used in both; `logs/gains_vs_spl.txt`): n_free equal in 13 runs (25800: 135 → 134); median g_B50 / g_spl 1.0025–1.0096 (3.49 GeV 1.004, 2.24 GeV H2 and C 1.008–1.010, 0.73 GeV 1.003, Ta 1.003–1.006); per-module RMS of the ratio 0.37–0.70 %. The 9 new runs match the 5 reused ones.

### 2.2 σ/E per run type

Candidate scan `output/resolution_fit_20261027_best50/r4_splits/`: every R3 `job.cmd` with only the output stem, `--mf-templates` (production → B50) and `--gain-file` (`percrystal_spl` → `percrystal_best50_20261026`) changed. Paired Δ = B50 − R3 (mean ± SE over splits). Source `output/resolution_fit_20261027_best50/compare_vs_r3.md`.

| Run type (fit) | Splits | R3 | B50 | Δ legacy | Δ gated |
|---|---:|---:|---:|---:|---:|
| 3.49 GeV H2 (CB) | 27 | 1.492 | 1.494 | +0.002 ± 0.005 | +0.003 ± 0.004 |
| 2.24 GeV H2 (Gauss) | 42 | 1.854 | 1.855 | +0.001 ± 0.002 | +0.001 ± 0.002 |
| 0.73 GeV H2 (Gauss) | 16 | 3.718 | 3.714 | −0.004 ± 0.011 | −0.003 ± 0.010 |
| 2.24 GeV C (CB) | 8 | 1.797 | 1.789 | −0.008 ± 0.007 | −0.011 ± 0.008 |
| 2.24 GeV Ta X17 (CB) | 25 | 1.978 | 2.010 | **+0.031 ± 0.020** | **+0.078 ± 0.022** |

ε_yield 0.9993–1.0006 at H2 and C, 0.9944 at Ta; Ta ε_peak 1.0138 (R3 1.0135, so not a bank effect).

### 2.3 Module medians and decomposition

Per-module median σ_M/μ_M post (`--group-runs`, N > 1000; R3 recomputed on the R3 tree, which reproduces the summary's R3 medians):

| Run type | R3 (modules) | B50 (modules) | Δ |
|---|---|---|---:|
| 3.49 GeV H2 | 1.519 (111) | 1.524 (111) | +0.005 |
| 2.24 GeV H2 | 1.842 (122) | 1.845 (123) | +0.003 |
| 0.73 GeV H2 | 3.497 (98) | 3.502 (98) | +0.006 |
| 2.24 GeV C | 1.810 (78) | 1.797 (78) | −0.013 |
| Ta X17 | 1.908 (8) | 2.001 (8) | +0.093 |

2-term fit σ/E = a/√E ⊕ b, H2 only, legacy (`out_r4.txt`, copy of `phase0.py` with scan key r4; it reproduces the R3 rows of `out_r3.txt`):

| Scan / estimator | σ/E at 0.73 / 2.24 / 3.49 (%) | a (%·√GeV) | b (%) | χ² (1 dof) |
|---|---|---|---|---:|
| R3, CB for all | 2.974 / 1.806 / 1.492 | 2.49 ± 0.02 | 0.68 ± 0.03 | 3.4 |
| B50, CB for all | 2.974 / 1.807 / 1.494 | 2.48 ± 0.02 | 0.70 ± 0.03 | 3.6 |
| R3, preferred (±2σ) | 3.718 / 1.854 / 1.492 | 3.18 ± 0.01 | −0.98 ± 0.02 | 374 |
| B50, preferred (±2σ) | 3.714 / 1.855 / 1.494 | 3.23 ± 0.01 | −1.05 ± 0.02 | 346 |

### 2.4 Per split and verdict

Per-split uncertainty σ_Δ = the split-to-split SD of the paired Δ within the run type (R3 and B50 use the same events). 8 of 118 splits have |Δ| > 2 σ_Δ, against 5.4 expected by chance, 3 up and 5 down; the largest are 24386.00080 +0.078 (z +3.2), 25257.00200 +0.035 (+3.2) and 24590.00060 −0.125 (−2.9). All eight and the tables per run type are in `compare_vs_r3.md`.

**Verdict.** At 3.49 GeV, 2.24 GeV H2, 0.73 GeV and C, B50 is the same as R3 (|Δ| ≤ 1.2 SE legacy, ≤ 1.4 SE gated; module medians within 0.013 pp; b 0.68 → 0.70 ± 0.03 %). At Ta it is worse: +1.6 SE legacy, +3.5 SE gated, module median +0.09 pp. Under the adoption rule (no run type worse beyond its uncertainty), B50 is not promoted. R3 remains the reference in the summary note.

### 2.5 Timing (reference only, not rerun)

B50 σ_t per module and per cluster: 20261026 §9.4 (`output/best50_20261026/compare/`); the per-cluster arm there used the `percrystal_spl` gains.

## 3. Caveats and open items

| Item | Detail |
|---|---|
| Ta weight | Ta is "reported, not interpreted" in the summary (≈ 16 calibrated modules per split, 8 modules in the median). The non-promotion rests on Ta alone. If Ta is excluded from the grading, B50 qualifies as equal at the four H2/C run types; that is a user decision. |
| Ta gated vs legacy | The Ta Δ is larger on the gated path (+0.078) than on legacy (+0.031); B50 Ta gated is 1.956 vs legacy 2.010, R3 1.878 vs 1.978. Not investigated. |
| In-sample | B50 gains and the R3 hybrid maps are both fitted on the splits they are applied to; the maps were fitted on R3 (production-bank) energies and not refitted for B50. |
| Not changed (§2) | Summary note, `database/`, code, `plans/`. The plan checkboxes of `plan_best50_pipeline_20261026.md` S2–S4 are not ticked (outputs live under the paths of this note, not `output/best50_20261026/`). |

## 4. How to run

```bash
cd /Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts && source setup_env.sh && export MPLBACKEND=Agg
G=output/percrystal_best50_20261026; D=output/resolution_fit_20261027_best50

# §2.1 dumps: the job list as run (dir + command per line) is $G/logs/dump_jobs_9runs.txt; each job ran as
#   mkdir -p <dir>; echo <cmd> > <dir>/job.cmd; <cmd> > <dir>/benchmark.log 2>&1   (12 parallel; log $G/logs/dumps_9runs.out)
for R in 024386 024512 024603 025228 025257 025260 025299 026138 026230; do
  X=""; [ $R = 024386 ] && X="$G/dumps_extra/024386_00060_dump.npz"
  python3 percrystal_calib.py fit --dumps "$G/dumps/scan_$R/split_*/result_dump.npz" $X --splits all --fixed-windows --out $G/$R/ > $G/fit_$R.log 2>&1; done

# §2.2 scan: job list $D/logs/scan_jobs.txt (R3 job.cmd with stem, --mf-templates, --gain-file substituted), same runner, log $D/logs/scan.out
python3 compare_scan_dirs.py --scans R3=output/resolution_fit_20261026/r3_splits R4=$D/r4_splits --gauss-window 2sigma --out $D/compare/r4_vs_r3
# §2.3
python3 plot_module_resolution.py --scan-dir $D/r4_splits --out $D/module_resolution/r4 --group-runs
python3 plot_module_resolution.py --scan-dir output/resolution_fit_20261026/r3_splits --out $D/module_resolution/r3_ref --group-runs
(cd $D && python3 phase0.py --scans spl,r2,r3,r4 --also-fixed150 > out_r4.txt)
# §2.4 tables (verdict added by hand into compare_vs_r3.md)
python3 $D/compare_vs_r3.py
```

## 5. CFD-fraction scan

Question (10-26 §7): does a CFD fraction f ≠ 0.5 improve σ_t, and how does a CFD time compare with the MF time? The fraction is hard-coded at 0.5 (`prad2dec/src/WaveAnalyzer.cpp:523`) and enters only the pulse time. For pulses above 10 ped_rms (all pulses with E ≥ 100 MeV) the production time is a log-normal fit evaluated at f (`fit_log_normal_cfd`), with the linear CFD as fallback.

### 5.1 Method

- `cfd_fraction_scan.py` (new) reads one split per run type: 24590_00010 (0.73 GeV H2), 25320_00050 (2.24 GeV H2), 24512_00010 (3.49 GeV H2), 25800_00100 (X17 Ta).
- MF time: best-50 banks (`output/templates_spline_best50_20261026`), run as in `module_timing_dump.py`; m_t and m_E are bit-identical to the 10-26 dump (§5.2).
- For the dominant pulse with m_E ≥ 50 MeV, a numpy port of the WaveAnalyzer pickoff re-evaluates the time at f = 0.2, 0.3, …, 0.7:
  - `lin_f`: the linear CFD on the smoothed waveform only.
  - `cfd_f`: the full production chain. The log-normal fit is done once per pulse, t(f) = t0 + exp(μ − σ√(−2 ln f)), and the f-dependent gate and linear-CFD fallback are re-applied. cfd_0.5 is the production time.
- The smoothing is a numpy port of `WaveAnalyzer::smooth`: the `prad2py` binding `wave_ana.smooth()` segfaults (`request()` is called with the GIL released).
- σ_t is the iterative ±2σ Gaussian-core width of the RF-folded phase per module (80 bins, `fit_unwrapped`), as tier e100 of `rf_time_calib.py`. MF-fallback rows (`m_fb`) are dropped. All 13 methods (MF, six lin_f, six cfd_f) use the same rows, and a module enters only if all 13 fits are ok.
- Bands: e100 (E ≥ 100 MeV, N ≥ 300 rows per module) and e500 (E ≥ 500 MeV, N ≥ 150). Values are medians over modules.

### 5.2 Validation

Source: `output/cfd_scan_20261027/validation.txt`. The C++ comparison uses 25320_00050, E ≥ 100 MeV, no MF fallback, WA peak matched within 8 ns (N 588 339).

| Check | Result |
|---|---|
| MF reproduction, 4 splits | same rows (1 525 499 to 3 606 737 per split); max \|Δm_t\| 0 ns, max \|Δm_E\| 0 MeV; event, index and fallback flags equal |
| f = 0.5, linear-CFD rows (code 1, N 22 173) vs C++ `peak.time` | 99.995 % within 1e-3 ns; RMS 8.0e-3 ns, from one row where C++ accepted the fit and the port did not |
| f = 0.5, fit rows (code 2, N 566 130, 96.2 % of compared rows) | 74.7 % within 1e-3 ns; RMS 0.145 ns, max 2.06 ns |
| f = 0.5, fallback rows (code 0, N 36) | RMS 5.3 ns: pulses rising at sample 0, where C++ reads `buf[-1]` (out of bounds); the port restricts j ≥ 0 |
| σ_t of the C++ time vs the port cfd_0.5, same rows and modules | e100 0.2093 vs 0.2078 ns, e500 0.1675 vs 0.1654 ns (median per-module ratio 1.005 / 1.007) |
| MF σ_t vs `rf_calib_modules.tsv` (10-26, e100, same modules) | median ratio 0.989–0.994 per run |

The fit-row mismatch is in the C++ fit itself. A scratch harness (a C ABI shim around `WaveAnalyzer::fit_log_normal_cfd`, compiled at -O3 in the session scratchpad, not in the repo) was run on 3062 waveforms of 25320_00050. The shim reproduces `peak.time` to 1e-3 ns on all of them, and its fit acceptance equals the port's. Port vs shim on the rows both accept (N 2745): 79.7 % within 1e-3 ns, RMS 0.108 ns. The same shim compiled with `-ffp-contract=off` vs the default shim: 78.8 % within 1e-3 ns, RMS 0.108 ns. The mismatching rows have t0 on its lower bound and χ²/dof of 10–1000; the 60-iteration LM does not converge there and stops at different points of a flat valley. The cfd_f values include this noise, as would a C++ rerun.

### 5.3 Results

Median σ_t (ns) vs f. Source: `output/cfd_scan_20261027/sigma_t_medians.md`; n = modules.

E ≥ 100 MeV (e100):

| f | 0.73 H2 lin | 0.73 H2 cfd | 2.24 H2 lin | 2.24 H2 cfd | 3.49 H2 lin | 3.49 H2 cfd | Ta lin | Ta cfd |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.2 | 0.1799 | 0.2648 | 0.1771 | 0.2653 | 0.1767 | 0.2633 | 0.1946 | 0.2793 |
| 0.3 | 0.1962 | 0.2549 | 0.1928 | 0.2563 | 0.1907 | 0.2503 | 0.2057 | 0.2703 |
| 0.4 | 0.1765 | 0.2862 | 0.1677 | 0.3122 | 0.1656 | 0.2981 | 0.1858 | 0.2969 |
| 0.5 | 0.1595 | 0.2073 | 0.1501 | 0.2078 | 0.1497 | 0.2068 | 0.1731 | 0.2263 |
| 0.6 | 0.1536 | 0.1874 | 0.1395 | 0.1800 | 0.1388 | 0.1782 | 0.1668 | 0.2037 |
| 0.7 | 0.2236 | 0.2248 | 0.2084 | 0.2236 | 0.2107 | 0.2199 | 0.2312 | 0.2374 |
| MF | 0.1727 | | 0.1474 | | 0.1468 | | 0.1914 | |
| n | 108 | | 163 | | 274 | | 133 | |

E ≥ 500 MeV (e500):

| f | 0.73 H2 lin | 0.73 H2 cfd | 2.24 H2 lin | 2.24 H2 cfd | 3.49 H2 lin | 3.49 H2 cfd | Ta lin | Ta cfd |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.2 | 0.1445 | 0.2317 | 0.1385 | 0.2319 | 0.1406 | 0.2370 | 0.1498 | 0.2416 |
| 0.3 | 0.1681 | 0.1970 | 0.1604 | 0.2059 | 0.1602 | 0.2081 | 0.1669 | 0.2078 |
| 0.4 | 0.1425 | 0.1837 | 0.1329 | 0.1946 | 0.1328 | 0.2059 | 0.1364 | 0.1955 |
| 0.5 | 0.1208 | 0.1650 | 0.1111 | 0.1654 | 0.1134 | 0.1691 | 0.1160 | 0.1729 |
| 0.6 | 0.1122 | 0.1450 | 0.0981 | 0.1360 | 0.1002 | 0.1363 | 0.1023 | 0.1428 |
| 0.7 | 0.1849 | 0.1819 | 0.1744 | 0.1761 | 0.1744 | 0.1681 | 0.1820 | 0.1788 |
| MF | 0.1163 | | 0.0976 | | 0.1064 | | 0.1058 | |
| n | 55 | | 145 | | 232 | | 86 | |

The per-module ratios (cfd_f/cfd_0.5, cfd_f/MF, lin_f/MF) are in the same file; the per-module values in `sigma_t_modules.tsv`; the curves in `cfd_scan.pdf`.

### 5.4 Discussion

(a) **The log-normal refinement degrades the time.** At f = 0.2–0.6, cfd_f is worse than lin_f in every split and band. At the production f = 0.5, cfd is 0.207–0.226 ns vs lin 0.150–0.173 ns (e100). At f = 0.7 the gap is 0.5–7 % at e100, and at e500 cfd is 1.6–3.6 % better in three of four splits. A likely contributor is the non-converging fit of §5.2.

(b) **f = 0.6 is the best fraction without a phase correction.** It is the minimum of both lin_f and cfd_f in all eight split-band combinations. For the production chain, cfd_0.6 lowers σ_t from cfd_0.5 by 10–14 % at e100 and 12–19 % at e500 (ratio of medians). The linear CFD at 0.6 beats MF by 6–12 % at e100 (median per-module lin_0.6/MF 0.88–0.94) and ties it at e500 (0.95–0.98; 2.24 GeV H2 median 0.0981 vs 0.0976 ns). The production chain does not beat MF at any f: cfd_0.6/MF is 1.07–1.23 at e100 and 1.22–1.36 at e500.

(c) **σ_t(f) is not smooth.** lin_f is poor at 0.3 and 0.7, cfd_f at 0.4 and 0.7. Much of the lin_f structure comes from the sub-sample-phase non-linearity (DNL) of the linear interpolation between samples; cfd_f adds the fit instability. A one-split in-sample check on 25320 (sub-sample phase mapped through its empirical CDF) makes lin σ_t monotone in f, 0.117–0.144 ns for f = 0.2–0.7 (`report.md`; not saved as a table). f = 0.6 is therefore optimal only without a phase correction.

(d) **Next steps (study items, not scheduled):** a sub-sample phase correction on all four splits, held-out; a finer scan at f = 0.55–0.65; a check of the walk vs amplitude.

**Decision.** No C++ change was made. A production change (fraction 0.6, dropping or fixing the log-normal refinement) awaits the user's decision; `report.md` recommends linear CFD at f = 0.6 without the refinement.

**Caveats.** One split per run type; e500 at 0.73 GeV has 55 modules.

### 5.5 Files and how to run

- `cfd_fraction_scan.py` (new): subcommands `dump`, `analyze`, `validate`.
- `output/cfd_scan_20261027/`: `dump/*.npz`, `sigma_t_modules.tsv`, `sigma_t_medians.{md,tsv}`, `cfd_scan.pdf`, `validation.txt`, `report.md`, `analyze.log`, `dump_driver.log`.

```bash
source setup_env.sh; export MPLBACKEND=Agg
python3 cfd_fraction_scan.py dump --parallel 4 --out output/cfd_scan_20261027/dump --templates-root output/templates_spline_best50_20261026
python3 cfd_fraction_scan.py analyze --dump-dir output/cfd_scan_20261027/dump --rf-dir output/module_timing_20261025/rf --out-dir output/cfd_scan_20261027
python3 cfd_fraction_scan.py validate --out-dir output/cfd_scan_20261027
```

## 6. Wassim single-pulse criterion vs best-50

Outputs in `output/templates_wassim_20261027/` (W). The criterion comes from `Waveform_Wassim/Hycal_prad_wf_analysis.pdf` §4 and `best_ref_shape_fit.C`, with page references in `W/criterion.md`. Per module, the template is ONE measured pulse:
- **Selection.** Peak in 150–170 ns, threshold 2 ADC, at least 20 pre-pulse samples, at least 4 rising samples and at least 10 strictly falling samples, no second peak above 0.15 A outside ±6 samples.
- **Ranking.** Strict-background tier first, then the highest amplitude, then the smoother background.
- **Shape.** Reference = (W − pre-mean)/max, a cubic interpolating spline through the 100 samples. No averaging, alignment or smoothing.

Readings taken here, because the PDF does not fix them:
- One pulse per module per run type, scanned over the same 119 splits as best-50, behind the production select trigger gate.
- Clipped waveforms are vetoed.
- Bank-format adaptations: an exponential tail beyond the readout (bank extent 400 ns); no f(A), so f ≡ 1; modules without a valid pulse take the run-type average (fallback 2).

### 6.1 Code (`spline_template_extract.py`; defaults byte-identical)

- `select --wassim-ref`: clip veto plus the PASS 1 `valid_p1` of the existing port. It keeps the single best candidate per module and split, and runs no mf.
- `build --single-pulse`: takes the best candidate over all files and writes the same npz/ROOT/PDF/summary.tsv. The npz has no ampcorr and no holdout_keys.
- Regression, OLD vs NEW: default select, default build and the best-50 extras select give identical arrays and identical summary.tsv.
- Review: critic PASS. Two fixes were applied: the onset-shift search went from ±8 to ±16 ns, and the spline end condition changed from natural to not-a-knot, as TSpline3 uses.

### 6.2 Banks and G1 gates

Modules with a strict-tier / fallback-tier reference / run-type average:

| Run type | Strict / fallback / average |
|---|---|
| 0.73 GeV H2 | 689 / 435 / 28 |
| 2.24 GeV C | 888 / 243 / 21 |
| 2.24 GeV H2 | 1089 / 62 / 1 |
| Ta X17 | 765 / 359 / 28 |
| 3.49 GeV H2 | 982 / 166 / 4 |

`check` ran on the best-50 held-out pulses. A check-only copy of each bank carried the best-50 `holdout_keys` (copies in the session scratchpad, not in W). Results are in `W/<slug>/check/gates.tsv`:
- **G1-4 (held-out χ²/dof) passes in all 5 run types**, but the χ² is 2–4 times the best-50 value on the same pulses. At 2.24 GeV H2, 300–1000 ADC: 155.8 against 74.8. At 1000–2500 ADC: 960.9 against 257.7.
- **G1-5 (sub-sample rise) fails in all 5.** The [1, 2) ns rise residual is −3.2 % against −0.3 % for best-50.
- **G1-6 (I_int vs the C++ integral) fails in all 5.** At 2.24 GeV H2, 25 modules have a phase-mean deviation over 0.5 % and 253 have no peak; best-50 has 0 and 0. This is expected from a single pulse sampled every 4 ns.
- **G1-8 (coverage):** 3.49 GeV passes (1148 modules); best-50 fails here (1065).

### 6.3 Pile-up (`W/compare/pileup_3way.md`)

Same production pseudo-data and reconstruction as 20261026 §9.1 (`benchmark_pileup_mf.py --tolerance-ns 4 --window-rule 40:0.16`); only the bank changes. Values are for the headline block, r_model 0.3. Eff_match is in % over 60 cells.

| | R3 | best-50 | Wassim |
|---|---:|---:|---:|
| Eff_match, pooled over 5 run types | 35.05 | 35.43 | 35.49 |
| Fake rate per event | 0.032 | 0.026 | 0.027 |

- **Gains.** The Wassim bank gains most at Δt 20 ns, ratio 1: 28.8 / 40.8 / 63.2 % pooled. Single run types reach +38 (0.73 GeV) and +36 (3.49 GeV) points over best-50.
- **Losses.** At 2.24 GeV H2 it loses 6–14 points in single cells (Δt 24 ratio 1: 100 → 86; Δt 32 ratio 0.5: 100 → 90). The fake rate there doubles (0.044 against 0.021 per event).
- **Well-separated cells.** These lose 0.4–1.2 points pooled.
- **Per run type, Wassim − best-50:** +0.60 (0.73), +0.33 (C), −1.03 (2.24 H2), −0.13 (Ta), +0.57 (3.49) points.

### 6.4 σ_E/E (`W/compare/energy_3way.md`)

The 20 splits of `best50_vs_r3.md` with the R3 configuration; only `--mf-templates` and `--gain-file` change. Gains come from `W/gains/<run>/` (percrystal recipe, all splits of the five runs, all converged; n_free within 1 of best-50). Legacy path, preferred fit, post arm, mean over 4 splits ± SE, in %. The R3 and best-50 columns reproduce `best50_vs_r3.md`.

| Run type | R3 | best-50 | Wassim | Δ Wassim − best-50 (pp) |
|---|---|---|---|---|
| 3.49 GeV H2 | 1.505 ± 0.006 | 1.507 ± 0.004 | 1.488 ± 0.006 | −0.020 ± 0.004 |
| 2.24 GeV H2 | 1.824 ± 0.004 | 1.826 ± 0.006 | 1.828 ± 0.006 | +0.003 ± 0.003 |
| 0.73 GeV H2 | 3.706 ± 0.037 | 3.673 ± 0.010 | 3.652 ± 0.029 | −0.021 ± 0.023 |
| 2.24 GeV C | 1.801 ± 0.023 | 1.798 ± 0.027 | 1.790 ± 0.023 | −0.008 ± 0.010 |
| 2.24 GeV Ta (X17) | 2.005 ± 0.042 | 2.015 ± 0.038 | 2.035 ± 0.055 | +0.020 ± 0.029 |

The gated path gives the same picture (3.49 GeV: −0.019 ± 0.004 pp).

### 6.5 Verdict

The single-pulse bank is not worse than best-50 for σ_E/E. It is 0.02 pp (about 5 SE) better at 3.49 GeV and equal within errors elsewhere, on 4 splits per type.
- **Shape quality is worse.** Held-out χ² is 2–4 times higher, and G1-5 and G1-6 fail.
- **Pile-up is mixed.** It is better at Δt 20 ns for equal amplitudes, worse at 2.24 GeV H2 (fake rate doubled).

Best-50 stays the reference bank. The 3.49 GeV gain would need the full 118-split scan before any claim.

### 6.6 How to run

```bash
W=output/templates_wassim_20261027; BK=output/templates_spline_best50_20261026   # RT = run-type slug; trigger as in the best-50 select (Ta: physics-any)
python3 spline_template_extract.py select $E/prad_025320.evio.00050 --trigger physics --wassim-ref --out $W/2p24GeV_H2/pulses/025320_00050.npz   # x119 splits, 5 parallel
python3 spline_template_extract.py build --single-pulse --run-type $RT --pulses "$W/$RT/pulses/*.npz" --out $W/$RT/templates_$RT
python3 spline_template_extract.py check --pulses "$BK/$RT/pulses/*.npz" --templates <bank copy + best-50 holdout_keys> --out $W/$RT/check/
python3 benchmark_pileup_mf.py output/pileup_eff_20261024_fix/$RT/synth --templates $W/$RT/templates_$RT.npz --out $W/pileup/W_prodsynth_wassim/$RT --tolerance-ns 4 --window-rule 40:0.16
python3 best50_compare.py pileup --ref output/pileup_eff_20261024_fix --new output/best50_20261026/pileup/B_prodsynth_newbank --arms W=$W/pileup/W_prodsynth_wassim --out $W/compare/pileup_3way   # renamed to pileup_3way_full.*; pileup_3way.md is a summary of pileup_3way_cells.tsv
# gains: the 44 dump job.cmd of output/percrystal_best50_20261026/dumps/scan_<run>/ with templates -> $W, stem -> $W/gains/dumps/...; then
python3 percrystal_calib.py fit --dumps "$W/gains/dumps/scan_<run>/split_*/result_dump.npz" --splits all --fixed-windows --out $W/gains/<run>/
# scan: the job.cmd of output/scan_20261026_best50/scan_<run>/split_<n>/ with --mf-templates -> $W, --gain-file -> $W/gains/<run>/, stem -> $W/scan/
python3 compare_scan_dirs.py --scans R3=output/resolution_fit_20261026/r3_splits best50=output/scan_20261026_best50 wassim=$W/scan --gauss-window 2sigma --out $W/compare/energy_3way_raw
```

## 7. Best energy resolution per run type (2026-10-27)

Numbers: 118-split scan, mean ± SE over splits, R3 reference chain (equal to or better than best-50 at every run type; best-50 kept for timing/pile-up only). Sources: `output/resolution_fit_20261026/out_r3.txt` (pooled rows), per-run Ta from `output/resolution_fit_20261027_best50/compare/r4_vs_r3_splits.tsv` (R3 columns). For the e-p runs the gated and legacy paths agree to ≤ 0.004 pp; the 16 ns gate matters only for Ta (≈ 5 % better). Ta runs are reported per run, never pooled (rule of 2026-10-27).

| Run type | splits | best σ_E/E (%) | fit / path | CB core width (%) | Methods applied |
|---|---|---|---|---|---|
| 3.49 GeV H2 | 27 | **1.491 ± 0.007** | CB, gated (= legacy 1.492) | 1.491 | common chain + hybrid bowl position correction |
| 2.24 GeV H2 | 42 | **1.853 ± 0.004** | Gauss ±2σ, gated (= legacy 1.854) | 1.804 ± 0.005 | common chain + hybrid bowl |
| 2.24 GeV C | 8 | **1.797 ± 0.012** | CB, gated = legacy | 1.797 | common chain + hybrid bowl |
| 0.73 GeV H2 | 16 | **3.714 ± 0.011** | Gauss ±2σ, gated (= legacy 3.718) | 2.973 ± 0.024 | common chain + bilinear position correction |
| 2.24 GeV Ta, run 25800 | 8 | **1.938 ± 0.030** | CB, gated (legacy 1.986) | — | common chain + bilinear; 16 ns gate essential |
| 2.24 GeV Ta, run 26138 | 4 | **1.920 ± 0.111** | CB, gated (legacy 1.950) | — | same |
| 2.24 GeV Ta, run 26230 | 13 | **1.828 ± 0.056** | CB, gated (legacy 1.983) | — | same |

Common chain (R3): FADC waveform → pedestal → matched-filter pulse fit with the production spline template bank (per module, per run type) → MF-amplitude energy with per-crystal spline gains (`percrystal_spl`, elastic peak vs E_beam) → island clustering with seed-anchored coincidence (|t_hit − t_seed| ≤ 16 ns, pre-window 100–200 ns) → position correction (hybrid bowl at 3.49/2.24 H2 and C, bilinear at 0.73 and Ta) → 2-ring edge exclusion → elastic window → fit (CB for 3.49 H2, 2.24 C, Ta; iterative ±2σ Gaussian for 2.24 and 0.73 H2) → mean ± SE over splits. "CB core width" gives the Crystal Ball σ where the Gaussian is the preferred fit; the two estimators are not interchangeable (CB excludes the radiative tail). Methods tried and not adopted: see `notes/PbWO4_methods_summary_20261027.md` §7.


## 8. Two-pulse χ² ceiling (U0, 2026-10-27)

- Question: is the MF's 0 % at pileup_3way label ≤ 16 ns (onset separation ≤ 9.8 ns) an information limit or an algorithm limit? Script `twopulse_fit.py`; output `output/ml_pileup_20261027/U0/` (`report.md`, `cells.tsv`, `eff_vs_dt.pdf`).
- Method: best-50 spline template at continuous time; 0.5 ns onset grid (±80 ns) with A ≥ 0 linear solves for 1 and 2 pulses, then `least_squares`. Pulse 2 is the later one; the decision is symmetric in time order; separation floor 1 ns. σ² = σ_ped(module)² + (0.02 y)². N = 2 iff χ²₁ − χ²₂ > threshold, the (1 − fake) quantile on even-base singles, applied to odd-base events.
- Data: one split per type (024603.00100, 025228.00100, 024386.00040). Clean real pulses (100–4000 ADC) plus a noise-free best-50 template (`real`). Order random; Δt = onset separation ± 2 ns phase; 2000 events per cell; run types never pooled. Arms: `tpl_noise`, `nf_same`, `nf_xbank` (Wassim bank fitted with best-50), `real2` (two real pulses).
- **Convention:** the pileup_3way Δt label = onset separation + 6.19 ns (label 20 = 13.8 ns, label 24 = 17.8 ns). The MF on the production synth with best-50 reproduces the published numbers (label 20 r1: 34 / 52 / 40 %; label 24: 92 / 100 / 94 %).
- **Δt 8 ns, r1 (target > 80 %):** fit 100 / 100 / 100 % (0.73 / 2.24 / 3.49 H2) on `real`, `tpl_noise`, `nf_same`, `real2`; 100 / 97 / 100 % on `nf_xbank`; MF 0 %. Evaluation fake rate 0.021 / 0.028 / 0.024 at the ≤ 0.03 operating point. In label terms, label 8 (1.8 ns onset) r1 gives fit 12 / 36 / 4 %, MF 0 %.
- The fit is 100 % at Δt ≥ 12 for every ratio ≥ 0.3, and at Δt 8 for r ≥ 0.5 (r 0.3: 96–100 %). The MF reaches ≥ 99 % only at Δt ≥ 20 (r 1), ≥ 16 (r 2) and 32 (r 0.5; 91–95 % at 24), and gives 0 % at r 0.3 through 24 ns. The MF's limit is algorithmic (12 ns floor, greedy argmax-then-stop), not informational.
- Below Δt 8 ns: Δt 4 ns (2–6 ns) r1 `real` 86 / 83 / 75 %, `tpl_noise` and `nf_same` 100 %, `nf_xbank` 84 / 36 / 84 %. Noise and the minimiser are not the limit; template-shape error is (real pulses deviate from the template by 1–2 % rms, ≈ 5 % on the leading edge). The 36 % comes from a few bad Wassim modules at 2.24 GeV.
- Resolution (r1, detected events): A bias < 1 % at Δt ≥ 8. RMS(A/A_true) is 4–6 % at 8 ns, 2.3–2.9 % at 16 and 1.7–2.4 % at 24 (MF 1.7–2.8 % at 24–32). Δt RMS is 0.2 ns at 8 and 0.15 ns at 16.
- Kill/pass (plan §6 U0): PASS. On the production synth the fit gives 100 % at label 20 r1 and label 16 r1, all types.
- Caveats: single-module waveforms only. The threshold is calibrated on clean pulses (mf_npk 1, ratio < 0.05), so the real-data fake rate is unmeasured (U5). The injected pulse is exact (`real`) or a real shape (`real2`, which agrees to ≤ 5 points). Amplitudes are 100–4000 ADC.
- Consequence: a network can at most approach this ceiling. At Δt ≥ 8 ns there is no room left for ML; U1–U4 would matter only below 8 ns (template-error-limited) or for robustness on real data. σ_E/E stays bounded by §7 (|Δ| ≤ 0.02 pp on e-p and C).

## 9. Cluster-level regression (V0/V0b, 2026-10-27)

- Script `cluster_regression.py` (numpy/scipy); output `output/ml_cluster_20261027/V0/`. Gate: E_R3 is rebuilt from the R3 dumps + pooled gains + poscorr_20261025 + the seed-indexed g_M. It reproduces R3 σ_E/E to ≤ 0.0002 pp on all 93 e-p splits that have an R3 value (limit 0.005 pp on the 4 reference splits per type).
- Ridge, augment arm: y = E′/E_R3, Huber-IRLS, plan §3 features (no absolute energy), trained on even split rank and evaluated on odd. λ by grouped CV (1e-4). Δσ_E/E is paired against untouched R3, mean ± SE, preferred fit, with the seed-indexed step refitted on the new energies.
- Δ (pp): 3.49 H2 −0.053 ± 0.008 (13 splits), 2.24 H2 −0.085 ± 0.003 (20), 2.24 C −0.071 ± 0.013 (4). The 0.73 H2 result is uninterpretable: the shuffle control gives −0.054 ± 0.014 and the pool is 59 % background.
- Shuffled-truth control: y is permuted within seed module (E′ is ~constant per module, so permuting E′ would be a no-op). Results: 3.49 +0.003 ± 0.004, 2.24 H2 −0.005 ± 0.001, C +0.002 ± 0.001 (vs R3; refitting g_M on R3 itself shifts σ/E by +0.005 / −0.004 / +0.003 pp, so against R3_gm the control is ≈ −0.002 pp, 0 within ~2 SE). The in-sample oracle bin lookup (module, u, v, E_seed/E_cl) gives −0.038 / −0.022 / −0.019 pp.
- The gain is an energy proxy. ln n_members alone gives −0.047 / −0.080 / −0.054 pp. Without it the gain is −0.026 ± 0.009 / −0.035 ± 0.002 / −0.014 ± 0.010. With the remaining features residualised on E_R3 and E_R3², it is +0.008 ± 0.004 / −0.005 ± 0.001 / +0.003 ± 0.001 vs R3, i.e. zero against R3_gm. A narrower training window (1 ± 3σ) increases the apparent gain (−0.072 at 3.49), which also points to sculpting.
- Kill rule (|Δ| < 0.02 pp at 3.49 for both ridge and oracle): not triggered by the letter. The credible gain is ≤ 0.026 pp (shape-only) and 0 when E-orthogonal. Decision: hold V1 (HGB/MLP), because a nonlinear model on the same features would learn the same energy proxy.
- Caveats: the R3 gains, maps and g_M are in-sample, which favours R3 by 0.01–0.02 pp. E-orthogonalisation also removes real leakage-fluctuation information, so the 0 is a lower bound on the shape-only gain.

## 10. Real-data pile-up rates (matched filter, 2026-10-27)

**What.** Every MF pulse of every PbWO4 channel, from 25 EVIO splits (first 60 000 physics events each, 3-4 splits per type spread over early/middle/late splits; Ta: 3 splits per run, not pooled), best-50 bank, MF nsigma 5 / min_amp 5 / kmax 4 / r_model 0.3, trigger (bits & 0xf00, no LMS/alpha). Script `pileup_rates.py` (`dump` per split, `analyse` from the npz); results in `output/pileup_rates_20261027/` (`report.md`, `rates_by_type.tsv`, `cluster_rates.tsv`, `module_rates_<g>.tsv`, `module_map_<g>.png`, `dt_ratio_<g>.tsv`, `dt_ratio_distributions.pdf`, `rf_comb.tsv`).

**Definitions.** Pulse 1 = the in-time pulse (largest-E pulse with 100 < t < 200 ns, the clusterer's hit); every other MF pulse of the channel is an extra pulse, dt = t_k - t_1, ratio = A_k / A_1 (can exceed 1). "any" = >= 2 pulses anywhere in the readout window; "near" = an extra pulse within 60 ns of pulse 1. WaveAnalyzer-fallback channels are excluded (rows per group: 3.49 GeV H2 2, 2.24 GeV H2 33, 2.24 GeV C 28, Ta 26138 1, 0.73 GeV H2 and the other Ta runs 0). Clusters: legacy `HyCalCluster`, E >= 50 MeV.

| group | events | events with a >= 2-pulse channel, any / near (%) | channels >= 2 pulses, any (%) | same, E1 >= 100 MeV (%) | extra pulses / event | cluster any member / seed (%) |
|---|---|---|---|---|---|---|
| 0.73 GeV H2 | 239 623 | 0.72 / 0.21 | 0.169 | 0.126 | 0.013 | 0.45 / 0.12 |
| 2.24 GeV H2 | 239 005 | 3.48 / 1.00 | 0.559 | 0.454 | 0.091 | 2.06 / 0.42 |
| 3.49 GeV H2 | 208 346 | 1.74 / 0.36 | 0.244 | 0.281 | 0.048 | 1.18 / 0.40 |
| 2.24 GeV C | 239 284 | 2.94 / 0.84 | 0.519 | 0.457 | 0.082 | 2.01 / 0.44 |
| Ta 25800 | 178 456 | 76.4 / 34.0 | 11.98 | 13.95 | 4.03 | 35.0 / 10.3 |
| Ta 26138 | 178 733 | 79.1 / 35.5 | 12.20 | 14.47 | 4.15 | 34.7 / 10.8 |
| Ta 26230 | 178 934 | 89.5 / 46.5 | 17.87 | 26.50 | 7.49 | 47.1 / 18.0 |

**Findings.**
- The event-level "any" rate for the e-p types (0.7-3.5 %) is of the same order as the 3.6-4.0 % "multi-peak channels" figure quoted in `notes/PbWO4_energy_resolution_20261022.md` line 64 for 2.24 H2 (3.48 %; that denominator was not checked, so no equality is claimed). Per channel the rate is 0.17-0.56 %. It is highest at 2.24 GeV and lower at 3.49 and 0.73 GeV (it does not rise monotonically with beam energy), and rises with cluster energy (2.24 H2: 1.3 % at 50-200 MeV, 3.0 % above 2 GeV). Seed radius (2.24 H2, any-member): 0-100 mm 2.9 %, 150-400 mm 0.4-0.6 %, and the outer edge > 400 mm is the highest bin (4.3 ± 0.6 %, 1 295 clusters). The seed channel itself is piled up in 0.4 % of clusters (2.24 H2).
- Ta runs have 12-18 % of channels with >= 2 pulses, mostly out-of-time background (near rate 4-6 %); the rate grows with seed energy (26230: 25 % at 50-200 MeV, 83 % above 2 GeV for any-member).
- Extra pulses with |dt| 16-60 ns, ratio 0.3-1 (2.24 H2): 0.090 % per channel; ratio < 0.3: 0.028 %; ratio > 1: 0.023 %. Extra pulses come after pulse 1 about twice as often as before (2.24 H2: 13 769 / 7 197).
- RF comb: for |dt| < ~260 ns the 4.000 ns FADC sample grid lies within 0.5 ns of the 4.0077 ns comb, so the pulse-1 on-comb fraction (31-37 % over the groups, random 25 %) equals the 4.000 ns grid fraction (31-36 %) and cannot be separated from sampling-grid locking. Pulse times show no grid locking (t mod 1 ns is flat to 0.5 %). Pulse-1 |dt| never exceeds ~238 ns, so the 260-400 ns numbers belong to within-channel pulse pairs (separation from the channel's earliest pulse; `rf_comb.tsv` rows "pair sep"), not to pulse-1 dt. The grid-decorrelated test is the pairs at |sep| 260-400 ns: Ta on-comb 29-31 % vs 26-27 % on the grid (25800: 31.2 / 27.2, 26138: 30.8 / 27.2, 26230: 29.3 / 26.4 %), a real but smaller RF excess; H2 and C have only 0-32 such pairs per group, so no test there.

**Caveat (mark everything as a lower limit for |dt| < 16-20 ns and ratio < 0.3).** The MF secondary search applies the 0.3 A_parent gate only after the parent onset, has a 12 ns MIN_SEP and stops after the argmax candidate (`mf_pulse.py` ~1054-1091; `output/templates_wassim_20261027/why_hard_cells.md`); cells in those regions carry `*` in `dt_ratio_<g>.tsv` and are shaded in the PDF. Errors are binomial; channels within an event are correlated. Event caps: 60 000 physics events per split; the splits are listed in `rates_by_type.tsv`. Top modules (N fired >= 2000) are the ring around the beam hole, e.g. 2.24 H2: W563 1.9 %, W560 1.9 %, W594 1.8 %.

Commands: `python3 pileup_rates.py dump --sample --parallel 8 --out output/pileup_rates_20261027/dump --max-events 60000` (about 7-12 min wall at `--parallel 8`, 25 splits), `python3 pileup_rates.py analyse --dump-dir output/pileup_rates_20261027/dump --out output/pileup_rates_20261027`. Figures and report in that directory are the regenerated versions: module maps colour modules with N fired >= 200, the top-5 lists keep N fired >= 2000, and the dt marginal is drawn over all bins with denominators named in the captions.

## 11. Hidden doublets in template pools (2026-10-27)

- Question: the best-50 spline templates are built from pulses that may contain a second pulse. Does that widen the Ta (X17) templates relative to e-p? Script `template_doublet_audit.py` (numpy/scipy; 111 s, 18 workers); output `output/templates_spline_best50_20261026/doublet_audit/` (`report.md`, per-bank pulse and module TSVs, rebuilt templates, PDF).
- Pool: the selected indices are not stored. They are rebuilt with `spline_template_extract`'s own loader, argv and hold-out (top 50 by `mf_energy` of the p0 or p1 chain pool). `n_pool` matches for every module of all 5 banks. `process_module` re-run on 8 modules per bank reproduces the shipped shape with max|Δf| = 0, also for the pool-restricted call used for the rebuilds.
- Fit: `twopulse_fit.fit_wave` (eps 0.02) of every pool pulse with its own template. Flag = Δχ² > thr, Δt ≥ 2 ns, secondary/main amplitude ≥ 0.05. Thr = U0 per type (0.73 H2 107, 2.24 H2 98, 3.49 H2 175); 2.24 C and Ta have no U0 value and use 98.
- Flagged fraction (|Δt| < 16 ns in brackets; strict = Δt ≥ 8 ns and ratio ≥ 0.2): 0.73 H2 0.38 (0.38) %, 2.24 C 0.86 (0.85), 2.24 H2 0.96 (0.93), 2.24 Ta 0.75 (0.73), 3.49 H2 0.46 (0.45). Strict: 0.010, 0.044, 0.079, 0.029, 0.002 %. Extrapolated 0.1 % (e-p) and 3 % (Ta) are not seen: Ta is at the level of 2.24 e-p and C. Most flags are Δt 2–8 ns with ratio < 0.5, and 396–4992 pulses per bank above thr are degenerate (Δt < 2 ns or ratio < 0.05) and excluded. With the literal A2/A1 (later/earlier) the fractions are 0.57 / 3.67 / 4.83 / 3.79 / 2.51 % (Ta |Δt| < 16 ns 2.91 %, 2.24 H2 4.64 %); the extra flags are ~11–14 ADC leading-edge precursors (A1/A2 ≈ 0.007) seen equally in e-p and Ta. All counts are upper bounds. Pool pulses are high energy, so template mismatch also lifts Δχ².
- Modules with ≥ 3 flagged of 50: 4 / 20 / 32 / 14 / 4 (same bank order). Selection bias, flagged % top-10 / bottom-10 amplitude (ratio): 0.42 / 0.62 (0.7), 1.73 / 0.66 (2.6), 1.30 / 0.98 (1.3), 0.91 / 0.88 (1.0, Ta), 0.82 / 0.26 (3.1).
- Shape impact (top-20 flagged modules, pool minus flagged, same build routine): Ta median |Δ| rise / w99 / FWHM 0.03 / 0.10 / 0.07 ns, max|Δf| median 1.6 % vs 1.5 % for removing the same number of random pulses (σ of one 50-pulse template rise / w99 0.11 / 0.23 ns; split-half max|Δf| 3.5 %, xtype 2.6–2.9 %). 2.24 H2: 0.10 / 0.18 / 0.09 ns, 1.9 % vs 2.3 %. About a third of the rebuilds switch the held-out λ, which also happens under random removal.
- Ta vs 2.24 H2, shipped templates, 1082 modules (Ta − H2, median [16, 84 %]): rise −0.010 [−0.29, +0.16] ns, w99 −0.04 [−0.53, +0.44], FWHM +0.03 [−0.29, +0.33]; noise σ of the difference 0.14 / 0.29 / 0.14 ns. Ta is wider by > 2σ in 121 (rise) and narrower in 176 modules. Run-type averages rise / FWHM: Ta 5.49 / 23.41, 2.24 H2 5.52 / 23.38, 0.73 H2 5.45 / 23.10, 3.49 H2 5.54 / 23.46. W325: Δrise −0.12 ns, Δw99 +0.59 ns (2.1σ), ΔFWHM 0.00, flagged pulses 0/50 (Ta), 1/50 (H2).
- Answer: no. Hidden doublets in the pools are rare and no more frequent in Ta than in 2.24 e-p, removing them moves the templates by the random-removal amount, and the shipped Ta templates are not wider than 2.24 H2 (median differences within noise; per-module scatter 1.5–2× split-half noise, both signs, uncorrelated with n_flag). Open: the Ta pool is lower-amplitude (45 % below 600 ADC vs 12 % for 2.24 H2), so an amplitude-dependent shape is the next suspect.
- Caveats: Δχ² threshold is calibrated on typical clean pulses (U0), not on the pool amplitude mix; pool pulses are template build inputs (lower Δχ², biases the count low); pool = pre-outlier-rejection top 50 (the build keeps about 40).
