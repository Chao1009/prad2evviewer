# PbWO4 energy resolution: ±2σ Gaussian, amplitude correction, position correction, LMS, clustering checks, Møller veto

**Date:** 2026-10-25 (work of 2026-10-24 15:10 onward)
**Branch:** `waveform_jz` (HEAD `80e19b6`; everything after it is uncommitted, §1.3)
**Plan:** `analysis/pyscripts/plans/plan_20261025.md` (Phases 0–5; "Decisions (user, 2026-10-24)"; Progress log)
**Previous notes:** `analysis/pyscripts/notes/PbWO4_energy_resolution_20261023.md` (per-crystal gains, centred rectangle, sampled templates), `analysis/pyscripts/notes/PbWO4_energy_resolution_summary.md` (state as of 2026-10-24)

Paths are relative to `analysis/pyscripts/` in `/Users/L00338853/work/PRad/prad2evviewer` unless they start with `database/`, `docs/`, `prad2det/` or `python/` (repo root), or are absolute. σ/E, preferred fit (CB for 3.49 GeV H2, C, Ta; Gaussian for 2.24 and 0.73 GeV H2), legacy/gated, "paired" (mean ± SE of per-split differences), Δ in pp and the run-type names are as in the summary note; legacy path unless noted. **Every Gaussian number in this note uses the ±2σ definition of §2** unless marked "fixed150". "Reference" = `output/scan_20261024_spline/` (banks `output/templates_spline/`, gains `output/percrystal_spl/`, κ = 1.0, centred rectangle). R1 = the reference configuration rerun with the centroid columns (§3.6); R2 = R1 + position correction (§4.8).

## 1. Status

Complete (`plans/plan_20261025.md` progress log to 10-25 01:42): all scans of P0–P4 and both add-ons are finished; nothing running at hand-off.

### 1.1 Scope and decisions

| Date (2026-10-24) | Decision | Effect |
|---|---|---|
| plan approval | Items in the order 2 → 1 → 3 → 4; Q1, Q3–Q7 plan defaults; Q2: position correction after the per-crystal gains, before the seed-indexed step; Q8: this note | §3–§6 |
| evening (during P2) | Add-on: software Møller veto, tested at 0.73 GeV with a 3.49 GeV control | §7 |
| ≈ 19:12 | Gaussian fits use an iterative ±2σ window everywhere | §2 |
| 19:30 | No checkpoint after P2; P3 → P4 → P5 run without stopping; the lead takes the R2/R1 decision from the §5.4 rule of the plan | §4.8 |
| ≈ 19:35 | Deadline 2026-10-25 06:00, no new jobs after 05:00; P4 at full scale in the order `cc1`, `mme2`, `mme0p5`, `mme5`, switching to the S_dev subset (36 splits) if P4 + P5 would end after 05:00 | §6 |

Deviation: every earlier output directory is read-only, so study-pulse dumps and tables go to `output/mf_fit_study_20261025/` (dumps `dump/<arm>/<run>_<split>/pulses.npz`, read with `--arm-dir`) and the reference-bank candidate diagnosis to `output/templates_spline_20261025/diag_ref*/`.

### 1.2 Code changes

| Script | Change (defaults byte-identical unless stated) |
|---|---|
| Benchmark `benchmark_hycal_timing_test_20261021.py` | `--write-cluster-tsv` forces `rec_ext` and appends `cluster_x_mm, cluster_y_mm` (C++ log-weighted centroid) to both paths' rows; without `--gain-file` it now also writes the existing `# member_identity` header line. `--pos-corr FILE` (off): E → E/R(u, v) after `--gain-file`, before the window fit, both paths; header `# pos_corr=… sha1=… form=… norm=… n_applied=… n_clamped=…`. `--veto-second-cluster EMIN_MEV` (off), `--veto-dphi DEG` (refused without EMIN), `--moller-gap-nsigma X` (3; 0 re-admits the Møller-gap modules) (§7). `--gauss-window {2sigma,fixed150}` (**default 2sigma**: changes the reported Gaussian columns, §2). `--cluster-min-module-energy MEV`, `--cluster-corner-conn {0,1}`, `--cluster-non-linear-corr {0,1}` (default: config value; `# cluster_cfg` line only when one is given). Neighbourhood path passes the module index to `bank.f_amp` (no-op for 1-D banks). |
| `spline_template_extract.py` | `build --bands E0,…` (generalised npz keys `band_edges, shape_band, r_band, n_band, A_med_band`; old keys kept); `ampcorr` with n knots; `ampcorr --per-module-f OUTDIR` (copies with a per-module `ampcorr_f[module, ln A grid]`); `check --amp-bins`. Defaults reproduce the reference banks byte for byte (2.24 GeV C build, ampcorr on copies, check). |
| `mf_pulse.py` | `from_npz` reads both band layouts; 2-D `ampcorr_f`, `f_amp(amp, midx)`. Two-τ default untouched. |
| `mf_fit_study.py` | `intraw --bins/--bank/--arm/--arm-dir`, `diag --out`. |
| `compare_scan_dirs.py` | `--gauss-window {2sigma,fixed150,stored}` (default 2sigma: refit from the stored inner-ring histograms; `stored` reproduces existing tables byte for byte). |
| `plot_module_resolution.py` | `--gauss-window`. |
| New: `position_corr.py` | `map` (12 × 12 bins, harness even/odd split sets via `G.r_pairs`, bootstrap errors, ring/side tables, `--gauss-window`), `fit --form {bilinear,poly8,both}`, `fit --global` (§4). |
| New: `mf_ampshape_diag.py` | Option (ii) offline refit diagnostic (§3.4). |
| New: `regauss_tables.py` | Old-vs-new Gaussian table from a scan's stored histograms (§2). |
| Unchanged | Driver `benchmark_hycal_timing_scan.py`, harness `gate_scan_20261022.py`, `percrystal_calib.py`, C++ and bindings, `database/`. |

### 1.3 Identity checks and git

| Check | Result |
|---|---|
| G0 after the centroid + `--pos-corr` edit (025320.00050, 20k events) | rawsum and mf: all TSVs identical except `elapsed_s`. REF0 + cluster TSV + dump: first 10 columns identical (32 760 rows); centroid vs dump `x_cl/y_cl` 16 027/16 027 matched, max \|Δ\| 0.0005 mm (3-decimal rounding). Identity map: all TSVs identical except the `# pos_corr` line (n_applied 32 760, n_clamped 132). |
| Full-scan identity, R1 vs reference | `output/scan_20261025_ref/` (119/119): `result.tsv`, `result_window.tsv`, `result_lms.tsv` identical except `elapsed_s` in all 119 splits; `result_clusters.tsv` first 10 columns identical. |
| G0 after the veto flags | Default byte-identical (on a copy, before install). |
| G0 after the Gaussian add-on | `--gauss-window fixed150` identical to the previous live benchmark; default changes only the Gaussian columns, two headline rows and one `# gauss_fit` line; `_window`, `_lms`, `_clusters` byte-identical; in-job 2σ equals the refit from the stored histogram. |
| G0 after the clustering flags | Identical to base. `--cluster-non-linear-corr 0`: all TSVs identical except the `# cluster_cfg` line (§6.1). |
| Install of the add-on code | Veto flags first installed at 19:14 (the veto matrix started on them). The fixed veto (true sub-50 MeV threshold), the ±2σ Gaussian and the clustering flags installed in the live files at 19:50 by atomic replace, after the review fixes (pre-install copies `$S/review/preinstall/`, diff `$S/review/install_addons.diff`). Jobs started before 19:50 ran on the previous file. |
| Configs and run scripts | Copied to `output/scan_20261025/configs/` (`cfgs_*.json`, `run_*.sh`, `form_cmp.py`). |
| Git | Modified since `80e19b6`: the benchmark, `compare_scan_dirs.py`, `mf_fit_study.py`, `mf_pulse.py`, `plot_module_resolution.py`, `spline_template_extract.py`, `notes/PbWO4_energy_resolution_summary.md`. Untracked: `mf_ampshape_diag.py`, `position_corr.py`, `regauss_tables.py`, `plans/plan_20261025.md`, this note. `output/` is gitignored. Nothing committed. |

## 2. Gaussian definition change (2026-10-24 ≈ 19:12)

**Definition.** Start from the old fit (Gaussian around the histogram argmax in a fixed ±150 MeV × kS window), refit in [μ − 2σ, μ + 2σ], iterate until μ and σ change by < 1 %, at most 5 refits; the final window and the refit count are reported (`# gauss_fit`, `/5*` when capped). `--gauss-window fixed150` restores the old fit.

**All Gaussian numbers before this change** (every note and plan up to 2026-10-24 19:12, and the stored Gaussian columns of every scan on disk before then) use the fixed ±150 MeV × kS window. CB values, gains, windows and cluster TSVs are unchanged by definition; the elastic-window μ_global fit, the CB seeds and the peak-yield s_ref stay on the fixed fit. Existing scans are recomputed offline from the inner-ring pre/post histograms stored in each `result.tsv`, the exact input of the in-job fit. Check: the `fixed150` refit reproduces the stored σ within 0.005 MeV on all 952 values of `scan_20261024_spline`.

Source: `output/regauss_20261025/{r1,spline,pc}.{md,tsv}` (`regauss_tables.py`). R1 = `output/scan_20261025_ref` (values identical to `spline`). Post σ/E (%), mean over splits; half-width = final (hi − lo)/2 in MeV; refits = mean count; cap = splits at the 5-refit cap.

| Run type (preferred fit), splits | Path | fixed150 | ±2σ | Δ (pp) | Half-width fixed150 / ±2σ | Refits | Cap | CB (unchanged) |
|---|---|---|---|---|---|---|---|---|
| 3.49 GeV H2 (CB), 28 | legacy | 1.783 | 1.745 | −0.038 | 187 / 122 | 2.0 | 0 | 1.709 |
| | gated | 1.783 | 1.745 | −0.039 | 187 / 122 | 2.0 | 0 | 1.710 |
| **2.24 GeV H2 (Gauss), 42** | legacy | **1.952** | **1.908** | −0.044 | 150 / 85 | 2.0 | 0 | 1.875 |
| | gated | **1.950** | **1.907** | −0.044 | 150 / 85 | 2.0 | 0 | 1.873 |
| **0.73 GeV H2 (Gauss), 16** | legacy | **4.366** | **3.713** | −0.653 | 86 / 54 | 4.0 | 5 | 2.954 |
| | gated | **4.364** | **3.712** | −0.652 | 86 / 54 | 4.1 | 5 | 2.953 |
| 2.24 GeV C (CB), 8 | legacy | 1.953 | 1.905 | −0.048 | 150 / 85 | 2.0 | 0 | 1.858 |
| | gated | 1.952 | 1.904 | −0.048 | 150 / 85 | 2.0 | 0 | 1.856 |
| Ta X17 (CB), 25 | legacy | 2.109 | 2.055 | −0.054 | 150 / 92 | 2.0 | 0 | 2.000 |
| | gated | 2.037 | 1.981 | −0.056 | 150 / 89 | 1.8 | 0 | 1.909 |

| Preferred-fit values that change | fixed150 → ±2σ, legacy (gated) |
|---|---|
| 2.24 GeV H2, reference = R1 (`scan_20261024_spline`) | 1.952 → **1.908** (1.950 → 1.907) |
| 0.73 GeV H2, reference = R1 | 4.366 → **3.713** (4.364 → 3.712) |
| 2.24 GeV H2, two-τ + pooled gains (`scan_20261024_pc`) | 2.071 → 2.042 (2.069 → 2.041) |
| 0.73 GeV H2, `scan_20261024_pc` | 4.440 → 3.771 (4.439 → 3.772) |

- Peak position (legacy post, `r1.tsv`) moves +0.54 MeV (3.49 GeV), +0.68 MeV (2.24 GeV H2), +0.67 MeV (C), +0.42 MeV (Ta) and +2.83 MeV at 0.73 GeV (723.10 → 725.93 MeV).
- **Non-convergence at 0.73 GeV.** 5 of 16 post spectra per path (8 of 16 pre) end at the 5-refit cap with σ still shrinking; the refit alternates between two bin sets on the tailed peak. The review added the converged flag; the user's definition with the cap is kept (no mean-of-last-two estimator). Capped spectra are flagged.
- CB results at 3.49 GeV, C and Ta, and every Δ quoted on a CB run type, are unaffected.

## 3. Amplitude correction (item 2)

### 3.1 Diagnosis

Study pulses (§6.1 population of the 20261023 note; 025320 × 3, 024503 × 3, 024590 × 1 splits, 20k events each). Median Q_template/Q_raw per acceptance bin (ADC), SE 0.0001–0.0007; spread = max/min − 1 over the bins. Source: `output/mf_fit_study_20261025/diag_{fixed,sampled}_{acc,fine}/intraw_bins_*.md`.

| Run | Arm | 100–250 | 250–500 | 500–1000 | 1000–2000 | 2000–4000 | Spread (%) |
|---|---|---|---|---|---|---|---|
| 25320 | two-τ fixed | 0.9778 | 0.9782 | 0.9813 | 0.9929 | 0.9942 | 1.68 |
| | sampled (reference banks) | 0.9975 | 0.9921 | 0.9909 | 0.9975 | 1.0014 | 1.06 |
| 24503 | two-τ fixed | 0.9843 | 0.9828 | 0.9833 | 0.9926 | 0.9965 | 1.40 |
| | sampled | 1.0002 | 0.9910 | 0.9898 | 0.9990 | 1.0032 | 1.35 |
| 24590 | two-τ fixed | 0.9804 | 0.9913 | 0.9962 | – | – | 1.61 |
| | sampled | 0.9901 | 0.9944 | 0.9972 | – | – | 0.72 |

- Sampled arm: a dip at 250–1000 ADC and a rise above 1000 ADC at 2.24 and 3.49 GeV; at 0.73 GeV a low first bin rising to 1000 ADC. Worst pairwise pulls 19–90 SE: not flat. Fine bins (24503): 1000–1500 0.9962, 1500–2000 1.0010, 2000–3000 1.0030, 3000–4000 1.0091.
- Fallback class: class 0 holds > 99 % of the pulses and equals "all"; class 1 deviates by up to 0.9 % (≈ 4 SE) in single bins (e.g. 25320 100–250: 1.0064 ± 0.0021 vs 0.9975), with no effect on "all".
- Template candidates (G1-7, held-out, reference banks, five acceptance bins): spread 0.34 / 0.63 / 0.94 / 0.45 / 0.97 % (0.73 / 2.24 H2 / C / Ta / 3.49).

### 3.2 Option (i): finer bands

Banks `output/templates_spline_20261025/<slug>/`, `build --bands 100,250,500,1000,2000,4000` from the stored pools; ampcorr scope per run type (run types differ by up to 9.7 SE).

| Run type | f at 100–250 / 250–500 / 500–1000 / 1000–2000 / 2000–4000 | G1-7 spread, new (reference) |
|---|---|---|
| 2.24 GeV H2 | 1.0097 / 0.9996 / 0.9943 / 0.9890 / 0.9842 | 0.73 (0.63) % |
| 3.49 GeV H2 | 1.0074 / 0.9993 / 0.9922 / 0.9964 / 0.9877 | 0.54 (0.97) % |
| 2.24 GeV C | 1.0066 / 1.0000 / 0.9938 / 0.9869 / 0.9819 | 0.92 (0.94) % |
| Ta X17 | 1.0044 / 1.0029 / 0.9983 / 0.9942 / 0.9913 | 0.18 (0.45) % |
| 0.73 GeV H2 | 1.0054 / 0.9979 / 0.9934 / – / – | 0.17 (0.34) % |

- G1-7 not flat in any run type. Per-module RMS of the G1-7 ratio at A > 300 (quadrature-corrected) unchanged: 0.86 / 0.45 / 0.26 / 0.22 / 0.31 %, so the option (iii) trigger fires.
- Study pulses (arm `fa`, `output/mf_fit_study_20261025/fa_acc/`): 25320 0.9989 / 0.9904 / 0.9886 / 0.9930 / 0.9920, spread 1.05 % (reference 1.06); 24503 1.0017 / 0.9924 / 0.9879 / 0.9969 / 0.9947, 1.39 % (1.35); 24590 0.9910 / 0.9933 / 0.9946, 0.37 % (0.72). The 250–1000 dip stays and deepens by 0.2 %: the noise-free band-template f does not track the in-window pulse population.

### 3.3 Option (iii): per-module f(A)

Built as `ampcorr --per-module-f` on the option (i) banks (post-processing; `r_band, n_band, A_med_band` are stored) → `output/templates_spline_20261025_pm/`. A module keeps its own chained value where its statistical SE < 0.3 %; other cells take the run-type value.

| Item | Detail |
|---|---|
| Review fix (round 2, blocking C1) | The first version set non-own cells to the run-type knot's absolute f, so a module's row dropped from its own level to the run-type level above its last own knot (median 0.9 %, max 10 %, in the seed amplitude range). The `fa_percrystal` dump pass (14 of 93 jobs) was stopped and its outputs deleted. |
| Rebuilt banks (v4) | Non-own cells chain from the nearest own knot; SE model on the capped count (c = 0.96 %, s0 = 0.57 %). The 0.3 % rule then keeps every cell with a ratio (all but 1 of 2068 at 0.73 GeV): non-binding on these banks. Round 3: PASS (one latent KeyError fixed; banks unchanged). |
| Outlier guard | Tried and dropped: W1014 at f = 0.886 is real on held-out pulses; the guard raised the 2.24 GeV H2 module RMS 0.08 → 0.39 %. |
| f range | 0.878–1.231; 0.878–1.087 without W335 (§9). |
| G1-7 (v4) | Spread 0.18 / 0.23 / 0.37 / 0.11 / 0.27 % (0.73 / 2.24 H2 / C / Ta / 3.49); Ta flat by 2 SE, the others not. Module RMS at A > 300 (corrected) 0.40 / 0.20 / 0.23 / 0.16 / 0.37 %. |

Study pulses, arm `fa_pm` on the v4 banks (`output/mf_fit_study_20261025/fa_pm_{acc,fine}/`):

| Run | 100–250 | 250–500 | 500–1000 | 1000–2000 | 2000–4000 | Spread (reference) | Module RMS, A > 300 (reference) |
|---|---|---|---|---|---|---|---|
| 25320 | 0.9973 | 0.9974 | 0.9970 | 0.9989 | 0.9990 | **0.19 %** (1.06) | 0.31 % (0.41) |
| 24503 | 0.9987 | 0.9976 | 0.9974 | 0.9978 | 0.9991 | **0.17 %** (1.35) | 0.33 % (0.47) |
| 24590 | 0.9951 | 0.9979 | 0.9986 | – | – | **0.36 %** (0.72) | 0.24 % (0.24) |

Not flat by 2 SE (pulls 9–13). 24503 fine bins: spread 0.66 %, 3000–4000 ADC bin at 1.0030.

### 3.4 Option (ii): amplitude-indexed templates, offline diagnostic

`mf_ampshape_diag.py` refits the stored study waveforms (≤ 200 per module) with (a) the module template + run-type f(A) of option (i) and (b) the shape interpolated linearly in ln A between the neighbouring band templates, no f. (a) reproduces the dump energies exactly. Source: `output/templates_spline_20261025/ampshape_diag/ampshape_diag.md`.

| Run | Spread (a) → (b), all pulses | Worst pull (b) | χ²/dof (a) → (b), 2000–4000 ADC | \|Median t_b − t_a\| per bin |
|---|---|---|---|---|
| 24503 | 0.58 → 0.16 % | 6.6 | 1.26 → 0.68 | ≤ 0.045 ns |
| 25320 | 0.58 → 0.23 % | 8.8 | 1.14 → 0.75 | ≤ 0.041 ns |
| 24590 | 0.19 → 0.23 % | 4.8 | – (500–1000: 1.02 → 0.99) | ≤ 0.026 ns |

(b) is not flat by 2 SE anywhere; χ²/dof is lower with (b) in every bin. This stored-wave subset shows a smaller (a) spread than the full in-window population of `intraw` (0.58 vs 1.05 %). Implementation in the matched filter needs the user's separate go (Q1); not done.

### 3.5 Choice and production

**Choice (Q1 default): option (iii)**, the smallest study-pulse spread. Option (i) does not qualify (spread not smaller than the reference at 2.24 and 3.49 GeV) and was not scanned.

| Step | Result |
|---|---|
| Dump pass `fa_percrystal` (pm banks) | Stage R 93/93, RX 25/25, no failures; + `output/percrystal_fa/dumps_extra/024386_00060_dump.npz`. |
| Gains `output/percrystal_fa/<run>/` | 14/14 converged in 2 iterations. vs `percrystal_spl` (free modules in both): median ratio 1.0004–1.0009 (3.49), 0.9989–0.9992 (2.24 H2), 0.9978 (C), 0.9964–0.9965 (0.73), 0.9995–1.0002 (Ta); per-module RMS 0.56–0.63 / 0.40–0.42 / 0.32 / 0.23–0.25 / 0.23–0.27 %. |
| Driver scan `output/scan_20261025_fa/` | 119/119, no failures. |

Production vs reference, σ/E (%), paired over 119 splits. Computed before §2: the Gaussian rows (2.24 and 0.73 GeV H2) are **fixed150**. Source: `output/scan_20261025_fa/compare/fa_vs_ref.{md,tsv}`.

| Run type (fit) | Reference | fa (pm banks) | Δ legacy ± SE [SD] | Δ gated | ε_yield | ε_peak |
|---|---|---|---|---|---|---|
| 3.49 GeV H2 (CB) | 1.709 | 1.753 | **+0.043 ± 0.004** [0.021] | +0.041 ± 0.004 | 0.9979 | 1.0001 |
| 2.24 GeV H2 (Gauss, fixed150) | 1.952 | 1.983 | **+0.031 ± 0.001** [0.008] | +0.031 ± 0.001 | 0.9983 | 1.0004 |
| 0.73 GeV H2 (Gauss, fixed150) | 4.366 | 4.360 | −0.006 ± 0.025 [0.099] | +0.008 ± 0.029 | 0.9992 | 0.9999 |
| 2.24 GeV C (CB) | 1.858 | 1.903 | **+0.044 ± 0.003** [0.008] | +0.046 ± 0.003 | 0.9993 | 1.0003 |
| Ta X17 (CB), report only | 2.000 | 1.996 | −0.003 ± 0.009 [0.047] | +0.008 ± 0.010 | 0.9942 | 1.0133 |

- Redone with the ±2σ Gaussian (`output/scan_20261025_fa/compare/fa_vs_ref_2sigma.md`): 2.24 GeV H2 +0.029 ± 0.002, 0.73 GeV −0.003 ± 0.007 pp; CB rows unchanged.
- **Rejected.** Every H2/C run type at 2.24 and 3.49 GeV is worse by > 2 SE (plan §4.5 item 3); Q1's fallback keeps the current f(A). Per-module medians (`output/module_resolution_20261025_fa/`): 3.49 GeV 1.61 → 1.63 %, 2.24 GeV H2 1.92 → 1.93 %, C 1.83 → 1.90 %, 0.73 GeV 3.60 → 3.60 %.
- Reading (hypothesis, not tested): flattening Q_template/Q_raw per module (study-pulse spread 1.1–1.4 → 0.17–0.19 %) does not improve the energy; the per-module band ratios (≤ 200 pulses per band) likely imprint template noise on every pulse, and Q_raw is not an unbiased truth.

### 3.6 R1

**R1 := `output/scan_20261025_ref`**: driver rerun of the reference configuration (banks `output/templates_spline/`, gains `output/percrystal_spl/`) with `--write-cluster-tsv`, so the cluster TSVs carry the centroid; identity-checked against `scan_20261024_spline` (§1.3). It was run before the acceptance test, during machine idle time.

## 4. Position correction (item 1)

### 4.1 Coordinates and pool

u = (`cluster_x_mm` − `seed_x_mm`)/20.77, v = (`cluster_y_mm` − `seed_y_mm`)/20.75 (reconstructed log-weighted coordinates). r = `cluster_E_cal_MeV`/E′(seed module), E′ = `eprime_MeV` of `result_window.tsv`. Pool = `plot_module_resolution.py`'s pool: legacy path, seed r < 200 mm, seed E ≥ 200·kE, outside the 1-based `15 20 15 20` rectangle, finite `cluster_E_cal` (Møller-gap and low-statistics seed modules already out). Ring k = max(|row − 17.5|, |col − 17.5|) − 2.5 (1-based); sides as in 20261023 §8a.4.

### 4.2 Maps on R1 and gate

`position_corr.py map`, all/even/odd, 12 × 12 bins over u, v ∈ [−0.6, 0.6], preferred fit of r per bin, `--min-n 300`, bootstrap errors (30 refits). The maps used for the fit were made before §2 (Gaussian bins fixed150); the ±2σ re-view is `output/poscorr_20261025/before_2sigma/`. R RMS with the statistical part removed. Source: `output/poscorr_20261025/map_summary.tsv`, `before_2sigma/map_summary.tsv`.

| Run type (fit) | Events | Bins | R peak-to-peak (%) | R RMS, fixed150 maps (%) | R RMS, ±2σ maps (%) | σ/μ edge / centre (%) | Gate < 0.3 % |
|---|---|---|---|---|---|---|---|
| 3.49 GeV H2 (CB) | 580 708 | 132/144 | 1.91 | **0.44** | 0.44 | 1.87 / 1.54 | fail |
| 2.24 GeV H2 (Gauss) | 1 159 603 | 140 | 1.31 (±2σ 1.34) | 0.17 | 0.17 | 2.12 / 1.85 (±2σ 2.08 / 1.83) | pass |
| 2.24 GeV C (CB) | 332 725 | 110 | 1.15 | 0.21 | 0.21 | 1.97 / 1.74 | pass |
| 0.73 GeV H2 (Gauss) | 465 118 | 132 | 3.81 (±2σ 4.89) | 0.24 (stat 0.24) | 0.37 | 4.72 / 3.78 (±2σ 4.35 / 3.43) | pass (±2σ: fail) |
| Ta X17 (CB) | 23 870 | 11 | 0.70 | 0.10 | – | – | report only |

- Edge = max(|u|, |v|) > 0.4; centre < 0.2.
- **Gate: FAIL** (3.49 GeV above 0.3 %), so the correction was fitted and scanned. At 0.73 GeV the ±2σ map RMS rises from 0.24 to 0.37 % because the narrower window follows the peak more closely.
- Even vs odd pull RMS 1.03 / 0.94 / 1.12 / 0.97 (0.73 / C / 2.24 H2 / 3.49), 94–95 % of bins within 2 (`eo_check.tsv`): the bootstrap errors are calibrated.

### 4.3 Fit, global vs per run type

`fit --form both` per run type and set → `output/poscorr_20261025/<slug>/poscorr_<slug>_<set>{,_poly8}.json` (bilinear default, normalised to a population-weighted mean of 1). poly8 residual vs the binned map 0.05–0.11 % at 2.24, 3.49 GeV and C, 0.25 % at 0.73 GeV. Global test (leave-one-out, `global_all/consistency_global.tsv`), fraction of bins within 2 SE: 3.49 GeV 0.31, 2.24 GeV H2 0.41, 0.73 GeV 0.67, C 0.91 (rule: ≥ 0.95 in every H2/C type) → **per run type** (Q3). Ta takes the 2.24 GeV H2 files (symlinks), reported only.

### 4.4 Ring and side dependence

Ring maps vs the run-type map (6 × 6, normalised over common bins, SE = se_g·√(1 − n_g/n_all)), n-weighted RMS / max of the difference:

| Run type | Ring 1 | Ring 2 | Rings 3–4 | Ring-1 sides |
|---|---|---|---|---|
| 3.49 GeV H2 | 0.15 / 0.39 % | 0.16 / 0.46 % | 0.29 / 0.99 % | 0.21–0.55 % RMS (bottom worst, max 1.66 %) |
| 2.24 GeV H2 | 0.11–0.27 % RMS over rings and sides, max 0.81 % | | | |
| 2.24 GeV C | 0.11–0.32 % RMS, max 1.06 % | | | |
| 0.73 GeV H2 | – | 0.13–0.17 % (rings 2, 3–4) | | 1.5–3 % RMS on 3–12 sparse bins |

Max |ΔR|/SE is 6–17 at 2.24 and 3.49 GeV (`map.log`): the ring dependence is significant and 3–10× smaller than σ/E. A per-type map captures the average only.

### 4.5 Bilinear vs poly8 (offline)

First-order test: E_cal/R(u, v) applied after the seed-indexed step, pooled σ/μ of E/E′ on the R1 pool. Δ in pp. Source: `output/poscorr_20261025/form_compare_offline.{py,tsv}`.

| Run type | Bilinear in-sample / e→o / o→e | poly8 in-sample / e→o / o→e |
|---|---|---|
| 3.49 GeV H2 | −0.096 / −0.088 / −0.100 | −0.096 / −0.091 / −0.101 |
| 2.24 GeV H2 | −0.010 / −0.010 / −0.009 | −0.010 / −0.011 / −0.009 |
| 2.24 GeV C | −0.021 / −0.021 / −0.022 | −0.021 / −0.024 / −0.019 |
| 0.73 GeV H2 | +0.016 / +0.004 / +0.011 | +0.013 / +0.001 / +0.001 |

Bilinear and poly8 agree within 0.003 pp at 2.24, 3.49 GeV and C. At 0.73 GeV the bilinear map (per-bin SE ≈ 0.25 %) adds a little noise. Only bilinear is scanned (Q4).

### 4.6 Held-out validation

Harness tags `pos_even` (even-split maps) on the 45 odd splits and `pos_odd` on the 48 even splits, stage R, R1 configuration + `--pos-corr`; no failures. The job lists equal the complementary map split sets exactly. Paired vs R1, ±2σ Gaussian, both halves pooled. Source: `output/scan_20261025/pos_{even,odd}/compare/heldout.{md,tsv}`, `output/scan_20261025/heldout_pooled.md` (`pooled_delta.py`; ε_peak = peak-yield ratio cmp/ref).

| Run type (fit) | n | Path | R1 | + pos | Δ ± SE [SD] | ε_yield | ε_peak |
|---|---|---|---|---|---|---|---|
| 3.49 GeV H2 (CB) | 27 | legacy | 1.709 | 1.611 | **−0.098 ± 0.005** [0.026] | 1.0000 | 0.9939 |
| | | gated | 1.710 | 1.610 | −0.100 ± 0.005 [0.028] | 1.0000 | 0.9939 |
| 2.24 GeV H2 (Gauss) | 42 | legacy | 1.908 | 1.899 | −0.009 ± 0.001 [0.007] | 1.0000 | 0.9992 |
| | | gated | 1.907 | 1.898 | −0.009 ± 0.001 [0.007] | 1.0000 | 0.9992 |
| 2.24 GeV C (CB) | 8 | legacy | 1.858 | 1.837 | −0.021 ± 0.005 [0.014] | 1.0000 | 0.9995 |
| | | gated | 1.856 | 1.834 | −0.022 ± 0.005 [0.013] | 1.0000 | 0.9995 |
| 0.73 GeV H2 (Gauss) | 16 | legacy | 3.713 | 3.707 | −0.006 ± 0.008 [0.033] | 1.0000 | 0.9852 |
| | | gated | 3.712 | 3.705 | −0.007 ± 0.008 [0.032] | 1.0000 | 0.9851 |

Per half (legacy, even-map / odd-map): 3.49 GeV −0.094 ± 0.007 / −0.102 ± 0.008; 2.24 GeV H2 −0.007 ± 0.001 / −0.011 ± 0.002; 0.73 GeV −0.004 ± 0.013 / −0.007 ± 0.011. Harness ε_yield and ε_peak ≥ 0.9997 in both halves. **Success rule met at 3.49 GeV** (Δ ≤ −0.02 pp, > 2 SE); no H2/C run type worse.

### 4.7 Full scan, module maps, after-map

`pos_all`: stage R 93/93, RX 25/25, no failures (118 splits, cluster TSVs). Paired vs R1, ±2σ Gaussian. Source: `output/scan_20261025/pos_all/compare/pos_vs_r1.{md,tsv}`, `pos_vs_r1_pooled.md` (ε_peak = peak-yield ratio cmp/ref).

| Run type (fit) | n | Path | R1 | R2 | Δ ± SE [SD] | ε_yield | ε_peak |
|---|---|---|---|---|---|---|---|
| 3.49 GeV H2 (CB) | 27 | legacy | 1.709 | 1.610 | **−0.099 ± 0.005** [0.025] | 1.0000 | 0.9937 |
| | | gated | 1.710 | 1.609 | −0.101 ± 0.005 [0.025] | 1.0000 | 0.9937 |
| 2.24 GeV H2 (Gauss) | 42 | legacy | 1.908 | 1.898 | −0.010 ± 0.001 [0.007] | 1.0000 | 0.9992 |
| | | gated | 1.907 | 1.897 | −0.010 ± 0.001 [0.008] | 1.0000 | 0.9992 |
| 2.24 GeV C (CB) | 8 | legacy | 1.858 | 1.831 | −0.027 ± 0.003 [0.009] | 1.0000 | 0.9993 |
| | | gated | 1.856 | 1.830 | −0.026 ± 0.003 [0.007] | 1.0000 | 0.9993 |
| 0.73 GeV H2 (Gauss) | 16 | legacy | 3.713 | 3.718 | +0.005 ± 0.008 [0.030] | 1.0000 | 0.9863 |
| | | gated | 3.712 | 3.714 | +0.002 ± 0.007 [0.028] | 1.0000 | 0.9863 |
| Ta X17 (CB), report only | 25 | legacy | 2.000 | 1.978 | −0.021 ± 0.021 [0.107] | 1.0000 | 1.0011 |
| | | gated | 1.909 | 1.878 | −0.031 ± 0.024 [0.120] | 1.0000 | 1.0013 |

Verdicts (harness rule): 3.49 GeV success; 2.24 GeV H2 better (|Δ| < 0.02 pp); C better (not graded); 0.73 GeV neutral; Ta report only. Harness ε_yield 0.9993–0.9999, ε_peak 1.0000–1.0003. Full-scan and held-out Δ (§4.6) agree within 0.006 pp at 3.49 GeV, 2.24 GeV H2 and C; at 0.73 GeV +0.005 vs −0.006 pp, both within 1 SE of zero.

Per-module σ_M/μ_M, median [16–84 %] (%), R1 → R2, ±2σ Gaussian / CB as preferred, `--group-runs`, N > 1000 (`output/module_resolution_20261025_{ref,pos}/summary/types_summary.tsv`):

| Run type | Modules fitted | Median [16–84 %] | Pooled σ/E |
|---|---|---|---|
| 3.49 GeV H2 | 111 → 111 | 1.609 [1.47, 1.85] → 1.582 [1.44, 1.75] | 1.710 → 1.616 |
| 2.24 GeV H2 | 124 → 123 | 1.873 [1.81, 1.96] → 1.879 [1.80, 1.95] | 1.908 → 1.898 |
| 2.24 GeV C | 78 → 78 | 1.833 [1.76, 1.98] → 1.827 [1.74, 1.95] | 1.862 → 1.839 |
| 0.73 GeV H2 | 99 → 98 | 3.515 [3.16, 8.43] → 3.497 [3.18, 5.70] | 3.724 → 3.717 |
| Ta X17 | 8 → 8 | 1.991 → 1.908 | 1.998 → 1.994 |

The pooled gain at 3.49 GeV (−0.10 pp) is larger than the per-module median gain (−0.03 pp): the correction acts mostly on the spread of module means and on the impact-position mix inside the pool, and little on single-module widths.

After-map (`output/poscorr_20261025/after/map_summary.tsv`, ±2σ; before = `before_2sigma/`), R RMS with the statistical part removed, and what remains (σ/μ edge vs centre, `after/sigma_*.png`):

| Run type | R RMS before → after (%) | R peak-to-peak before → after (%) | σ/μ edge / centre after (%) |
|---|---|---|---|
| 3.49 GeV H2 | 0.44 → **0.05** | 1.91 → 0.45 | 1.85 / 1.55 |
| 2.24 GeV H2 | 0.17 → 0.00 | 1.34 → 0.54 | 2.08 / 1.81 |
| 2.24 GeV C | 0.21 → 0.00 | 1.15 → 0.38 | 1.95 / 1.76 |
| 0.73 GeV H2 | 0.37 → 0.20 | 4.89 → 4.55 | 4.37 / 3.41 |

The gate (< 0.3 %) passes everywhere on the after-map. The edge resolution stays 0.2–0.3 pp worse than the centre (≈ 1 pp at 0.73 GeV); this is a width effect, which a scale map cannot remove.

### 4.8 Decision

**Position correction adopted; R2 := `output/scan_20261025/pos_all`** (lead's decision under the plan §5.4 rule; checkpoint skipped by the user's 19:30 instruction). Configuration: R1 + `--pos-corr output/poscorr_20261025/{runtype}/poscorr_{runtype}_all.json` (bilinear, per run type; Ta with the 2.24 GeV H2 map). R2 holds 118 splits (no 24386.00060); comparisons pair on common splits. P3 and P4 run on R2.

## 5. LMS for all runs (item 3)

Scope (fixed 2026-10-24, before this package):

| Item | Detail |
|---|---|
| Local files | `output/gain_corr_20261024/` holds usable files for 24503 (2 batches) and 24512 (1 batch) only, at `-b 1600`; the local splits are a sparse sample with 220–3600 LMS events per run, and the tool drops the trailing partial batch. |
| (a) Two-run check | Harness tag `lms_2run` = R2 + `--lms-corr output/gain_corr_20261024/prad_{run}_gain_corr.root`, runs 24503 and 24512 (17 splits). 24503 carries one real gain step; 24512 (one batch, `first-batch` gives f ≡ 1) is a pipeline identity check. |
| (b) Drift bound | The `-b 200` diagnostic over the sampled splits of all 12 runs (`output/gain_corr_20261024/b200_diagnostic/`) puts the raw per-module LMS means within 0.9996–1.0006: in-run drift ±0.05 % (≤ 0.06 %), two orders of magnitude below the constant term. The larger b200 jumps are failed alpha fits on the reference PMTs: **W1124** in every run, **W781** in the Ta runs. On 25320, where a full database file exists, LMS on/off leaves σ/E unchanged (20261021 §3). |
| (c) Other 10 runs | LMS off, as in every scan so far. |
| Deferred (ifarm) | Full-run gain_corr files need the complete runs: `replay_gainCorr -b 1600` (or 4000) over every split of each run on ifarm, the 12 files copied into `database/gain_factor/gain_correction/` there. Never from local partial data. Not part of this package. |

Results: `lms_2run`, 17 jobs, no failures, `output/scan_20261025_p3/lms_2run/` (separate harness out dir: the P2 harness was writing `output/scan_20261025/manifest.json`). Paired vs R2 (`lms_2run/compare/lms_vs_r2.{md,tsv}`, ±2σ Gaussian):

| Run type (fit) | Splits | R2 | + LMS | Δ legacy ± SE [SD] | Δ gated | ε_yield |
|---|---|---|---|---|---|---|
| 3.49 GeV H2 (CB), 24503 + 24512 | 17 | 1.619 | 1.619 | **+0.000 ± 0.000** [0.002] | +0.000 ± 0.001 | 0.9999 |

| Run | Detail |
|---|---|
| 24503 | Two batches; batch 1 starts at event 8 231 291, 25.7 min after batch 0. Per-module factors of batch 1: median 1.0005, p05–p95 0.9990–1.0021, min 0.9966. Splits 00010–00040 lie entirely in batch 0 (f ≡ 1), 00050 straddles, 00060–00100 lie in batch 1. |
| 24512 | One batch, f ≡ 1: `result.tsv` and `result_window.tsv` identical to R2 except `elapsed_s` and the `lms_corr` header in all 7 splits (`result_lms.tsv` differs by construction: R2 has no LMS file for this run). |

LMS tracking leaves σ/E unchanged where it can be applied; the drift bound (b) and the deferred ifarm item stand.

## 6. Clustering checks (item 4)

### 6.1 Constants and the non-linearity identity job

| Constant (pipeline JSON) | Value | Override flag | Variant tags |
|---|---|---|---|
| min_module_energy | 1.0 MeV | `--cluster-min-module-energy` | `mme0p5`, `mme2`, `mme5` (0.5 / 2 / 5 MeV) |
| corner_conn | false (4-connected) | `--cluster-corner-conn` | `cc1` |
| non_linear_corr | true (struct default, not copied from the JSON) | `--cluster-non-linear-corr` | identity job only |
| min_center_energy / min_cluster_energy | 10 / 50 MeV | none | – |
| log_weight_thres | 3.6 | none | – |

| Calibration files in use | nl keys | `linear_corr` |
|---|---|---|
| H2/C `_fadc` files (0.7, 2.2, 3.5 GeV) | only the legacy key `non_linear` = 0, which `LoadCalibration` never reads | 1 |
| The three X17 `_fadc` files | nl1 = nl2 = 0 | 1 |

`linear_corr` = 1 everywhere, and with `--gain-file` the cluster energy bypasses it. Identity job (025320.00050, 20k events, REF0 ± `--cluster-non-linear-corr 0`): all TSVs identical except the `# cluster_cfg` line. **`non_linear_corr` off is an exact null.**

### 6.2 Variants

Harness, R2 configuration + one flag, full scale (stages R + RX), paired vs R2. Winner rule: Δ ≤ −0.02 pp with |Δ| > 2 SE at 2.24 and 3.49 GeV H2, no H2/C type worse by > 2 SE, ε_yield and ε_peak ≥ 0.99. Pooled gains are not refit per variant; a winner gets a confirmation pass with refit gains (≈ 2.75 h) before adoption.

`cc1` (corner_conn on), full scale (R 93/93, RX 25/25, no failures). Paired vs R2 (`cc1/compare/cc1_vs_r2.{md,tsv}`, pooled `cc1_vs_r2_pooled.md`, ±2σ Gaussian); σ/E and ε on the legacy path:

| Run type (fit) | Splits | R2 | `cc1` | Δ legacy ± SE [SD] | Δ gated ± SE [SD] | ε_yield | ε_peak |
|---|---|---|---|---|---|---|---|
| 3.49 GeV H2 (CB) | 27 | 1.610 | 1.611 | +0.001 ± 0.003 [0.018] | +0.002 ± 0.003 [0.016] | 0.9993 | 1.0002 |
| 2.24 GeV H2 (Gauss) | 42 | 1.898 | 1.893 | −0.005 ± 0.001 [0.008] | −0.005 ± 0.001 [0.008] | 0.9996 | 0.9988 |
| 2.24 GeV C (CB) | 8 | 1.831 | 1.817 | −0.014 ± 0.005 [0.015] | −0.013 ± 0.005 [0.015] | 1.0001 | 0.9983 |
| 0.73 GeV H2 (Gauss) | 16 | 3.718 | 3.640 | **−0.077 ± 0.009** [0.035] | −0.071 ± 0.009 [0.035] | 0.9997 | 0.9954 |
| 2.24 GeV Ta (X17) (CB) | 25 | 1.978 | 2.010 | +0.031 ± 0.027 [0.133] | +0.033 ± 0.032 [0.158] | 0.9986 | 0.9918 |

CB on the Gaussian run types, legacy, paired (`cc1_vs_r2_splits.tsv`): 0.73 GeV 2.974 → 2.943 %, Δ −0.031 ± 0.017; 2.24 GeV H2 1.858 → 1.849 %, Δ −0.009 ± 0.003.

**`cc1` does not meet the winner rule**: 3.49 GeV H2 is neutral (+0.001 ± 0.003) and the 2.24 GeV H2 gain (−0.005 ± 0.001) stays below 0.02 pp. The gain is at 0.73 GeV (Gaussian −0.077 ± 0.009) and on C (−0.014 ± 0.005).

`mme2` (2 MeV), full scale (118/118): **exact null**. Δ = 0.000 ± 0.000 [0.000] and ε = 1.0000 in every run type and path (`mme2/compare/mme2_vs_r2_pooled.md`). `result.tsv`, `result_window.tsv` and `result_lms.tsv` are identical to R2 in all 118 splits except the `# cluster_cfg` line and, in 11 splits, the floating-point `member_identity maxrel` diagnostic (10⁻⁷). The reconstructed pulse energies never fall below ≈ 2.2 MeV at 2.24 and 3.49 GeV and ≈ 3.7 MeV at 0.73 GeV (sampled study dumps: mf minimum 2.17–3.67 MeV, fraction below 2 MeV 0.0000), so any `min_module_energy` ≤ 2 MeV removes no hit.

`mme0p5` (0.5 MeV, below the current 1 MeV), full scale (118/118, the first split resumed, no failures): **exact null**, as the energy floor predicts. Δ = 0.000 ± 0.000 [0.000] and ε = 1.0000 in every run type and path (`mme0p5/compare/mme0p5_vs_r2_pooled.md`). `result.tsv` and `result_window.tsv` are identical to R2 in all 118 splits apart from the `# cluster_cfg` line and the floating-point `member_identity` diagnostic.

`mme5` (5 MeV), full scale (118/118, no failures): **neutral** in every H2/C run type (`mme5/compare/mme5_vs_r2_pooled.md`). Δ legacy: 3.49 GeV H2 −0.001 ± 0.001, 2.24 GeV H2 −0.000 ± 0.000, 0.73 GeV H2 +0.001 ± 0.001 (Gaussian) and +0.001 ± 0.001 (CB), C −0.000 ± 0.001; ε_yield 1.0000 in every H2/C run type. Ta −0.027 ± 0.014 (report only).

**No P4 variant meets the winner rule**, so no confirmation pass is triggered. `cc1` gains only at 0.73 GeV and on C; `mme2` and `mme0p5` are exact nulls (energy floor); `mme5` is neutral.

## 7. Software Møller veto at 0.73 GeV (add-on)

| Item | Detail |
|---|---|
| Rule | `--veto-second-cluster EMIN_MEV`: reject a single-cluster event if any other cluster (PbWO4 or lead glass) has energy > EMIN, per path (legacy and gated use their own clusters). `--veto-dphi DEG`: the vetoing cluster must also lie within ±DEG of back-to-back in azimuth (\|Δφ − 180°\| ≤ DEG, centroid x, y). Veto count in the result header. Default off, byte-identical. |
| Re-admission | `--moller-gap-nsigma X` (default 3): seed modules with an elastic–Møller gap < X σ_M are excluded from the window fit; 0 re-admits them. |
| 50 MeV issue and fix | Review finding: EMIN < 50 MeV was silently clamped to the clusterer's `min_cluster_energy` (50 MeV), so the first v30 and v30p20 rows were effectively 50 MeV. Fix (installed 19:50): sub-50 MeV clusters taken from `get_clusters()`; EMIN ≥ 50 unchanged; `--veto-dphi` without EMIN refused. The clamped v30/v30p20 outputs were moved to `$S/v30_clamped/` and rerun with the fixed veto (20:38–20:50). |
| Matrix as run | `output/scan_20261025_veto/`, harness stage R, R1 configuration (no position correction) + `--max-events 40000 --write-cluster-tsv`, `--parallel 6`. 0.73 GeV: 24590 + 24603, 16 splits × 10 tags: `vb` (base), `v30`, `v50`, `v100`, `v30p20`, `v50p20`, `v100p20` (±20°), `vb_rm`, `v50_rm`, `v50p20_rm` (`--moller-gap-nsigma 0`). 3.49 GeV accidental-veto control: 24503, 10 splits × `vb`, `v30`, `v50p20`. |
| Full-statistics confirmation | `vf30` (EMIN 30, no φ) and `vf50p20` (EMIN 50, ±20°) on all 16 full 0.73 GeV splits (32 jobs, no failures), R1 configuration, paired with R1 (`cfgs_vetofull.json`, `run_vetofull.sh`). |

Screening matrix, 40k events per split, legacy path, Δ vs `vb` paired. CB σ/E stored; r < 0.90 = fraction of the pooled E_cal/E′ below 0.90. Source: `output/scan_20261025_veto/analysis/veto_table.{py,md,tsv}`.

| Tag | Vetoed (%) | CB σ/E | ΔCB ± SE (pp) | r < 0.90 (%) | r < 0.90 kept (%) | CB peak μ − 1 (%) |
|---|---|---|---|---|---|---|
| `vb` (0.73 GeV) | 0 | 3.061 | – | 57.0 | 100 | +0.69 |
| `v30` | 24.4 | 3.120 | +0.060 ± 0.031 | 42.5 | 54.3 | +0.61 |
| `v50` | 16.7 | 3.080 | +0.019 ± 0.025 | 46.7 | 65.3 | +0.65 |
| `v100` | 6.8 | 3.079 | +0.018 ± 0.013 | 52.8 | 84.5 | +0.67 |
| `v30p20` | 2.6 | 3.072 | +0.012 ± 0.011 | 55.8 | 94.9 | +0.66 |
| `v50p20` | 1.6 | 3.071 | +0.010 ± 0.009 | 56.2 | 96.9 | +0.67 |
| `v100p20` | 0.6 | 3.069 | +0.008 ± 0.006 | 56.7 | 98.9 | +0.67 |
| `vb_rm` | 0 | 4.147 | +1.086 ± 0.106 | 51.6 | 202.7 | +0.62 |
| `v50_rm` | 16.7 | 4.173 | +1.113 ± 0.097 | 41.3 | 131.9 | +0.58 |
| `v50p20_rm` | 1.6 | 4.147 | +1.086 ± 0.105 | 50.9 | 196.7 | +0.60 |
| `vb` (3.49 GeV, 24503) | 0 | 1.671 | – | 11.3 | 100 | +0.02 |
| `v30` | 41.6 | 1.628 | −0.043 ± 0.009 | 5.7 | 38.2 | +0.00 |
| `v50p20` | 22.7 | 1.670 | −0.002 ± 0.002 | 8.5 | 72.6 | +0.02 |

- No setting improves the 0.73 GeV CB σ/E (veto-only rows: all ΔCB positive and ≤ 2 SE; the `_rm` rows are ≈ 10 SE worse). Gated within 0.006 pp of legacy.
- Møller-gap modules: 12.8 excluded per split (union 16) at 0.73 GeV, 0 at 3.49 GeV. Re-admission (`_rm`) doubles the pool and raises the CB σ/E by 1.09 pp; restricted to seeds outside `vb`'s Møller-gap set, the `_rm` pool reproduces the `vb` pool exactly, so the excess comes from the re-admitted modules.
- At 40k events per split only ≈ 27 modules are calibrated at 0.73 GeV and the ±2σ Gaussian is unstable (`vb` 9.31 % vs 3.71 % at full statistics), so the screening Gaussian Δs (−1.38 ± 0.41 for `v30`) are not used.
- At 3.49 GeV the vetoed fraction includes real second clusters (radiative photon, Møller): an upper bound on the accidental rate.

Full statistics, 0.73 GeV, 16 splits, legacy, paired with R1 (`analysis/cmpfull_{vf30,vf50p20}.{md,tsv}`):

| Tag | Vetoed (%) | ε_yield | Gaussian σ/E R1 → tag | ΔGauss ± SE | CB σ/E R1 → tag | ΔCB ± SE |
|---|---|---|---|---|---|---|
| `vf30` | 24.4 | 0.756 | 3.713 → 3.629 | **−0.083 ± 0.009** | 2.954 → 2.993 | +0.039 ± 0.022 |
| `vf50p20` | 1.55 | 0.984 | 3.713 → 3.713 | +0.001 ± 0.004 | 2.954 → 2.993 | +0.038 ± 0.008 |

The benchmark without the veto reproduces R1 on 024590.00010 exactly (CB, windows), so the CB shifts are real: removing tail events changes the CB tail parameter (α 0.331 → 0.375 on that split) and raises the core σ.

**Recommendation: no veto in the 0.73 GeV reference configuration.**
- The back-to-back veto (the Møller signature) fires on 0.6–2.6 % of events and keeps 95–99 % of the E/E′ < 0.90 population: the Møller partner is rarely seen as a separate cluster.
- The no-φ veto lowers the Gaussian σ/E by 0.08 pp at full statistics and discards 24 % of the events (42 % at 3.49 GeV, mostly non-Møller second clusters).
- No setting improves the CB σ/E; re-admitting the Møller-gap modules is worse (CB +1.09 ± 0.11 pp).
- For Møller-tail studies, EMIN 30 MeV without φ is the setting that changes the tail.

## 8. Resolution decomposition refit

Source: `output/resolution_fit_20261025/phase0.py` and `out.txt` (`python3 phase0.py --also-fixed150`; scans spl 119, r1 119, r2 118 splits; ±2σ Gaussian refit from the stored histograms, CB stored). `--gauss fixed150` reproduces `output/resolution_fit_20261024/out.txt` (χ² 574.7 vs 575.1, stored-value rounding). σ² = a²/E + b² (2-term, 1 dof) and + c²/E² (3-term, exactly determined); a in %·√GeV, b in %, c in %·GeV. Legacy; gated within 0.01.

| Scan / estimator | σ/E % at 0.73 / 2.24 / 3.49 | 2-term a, b (χ²) | 3-term a, b, c |
|---|---|---|---|
| R1, CB for all | 2.954 / 1.875 / 1.709 | 2.20, 1.20 (58) | 1.55, 1.44, 1.33 |
| **R2, CB for all** | 2.974 / 1.858 / 1.610 | **2.37, 0.98 (3.6)** | **2.25, 1.05, 0.64** |
| R1, preferred (±2σ Gaussian at 0.73, 2.24) | 3.713 / 1.908 / 1.709 | 3.07, −0.58 (1150) | 0.40, 1.55, 2.43 |
| R2, preferred (±2σ) | 3.718 / 1.898 / 1.610 | 3.12, −0.78 (616) | 1.77, 1.16, 2.08 |
| R2, preferred (fixed150, continuity) | 4.338 / 1.942 / 1.610 | 2.91, 0.11 ± 0.29 (414) | 1.17, 1.24, 2.85 |

- **The preferred (mixed) estimator is unusable for the decomposition with the ±2σ Gaussian**: the 0.73 GeV value is the core of a tailed peak (54 MeV half-width), so b² < 0 in the 2-term fit. CB for all is the consistent estimator.
- R2, CB for all: the 2-term fit is acceptable (χ² 3.6 for 1 dof) and the constant term drops from 1.20 to 0.98 % (2-term) and from 1.44 to 1.05 % (3-term). 3-term shares at 0.73 / 2.24 / 3.49 GeV: a/√E 2.64 / 1.51 / 1.21, b 1.05, c/E 0.88 / 0.29 / 0.18 %.

## 9. Caveats and open items

| Item | Detail |
|---|---|
| W335, 0.73 GeV | Broken module template in the option (i)/(iii) 0.73 GeV bank: χ²/ndf 3375 (≤ 119 for every other 0.73 GeV module), I_int 4.52 vs 7.0–7.3 for the same module in the other run types and the two-τ value, lohi_iint_ratio 0.65. Its per-module f (1.231 / 1.133 / 1.100) partly offsets the bad integral. One module of 1103; not rebuilt. Also the 0.73 GeV pile-up hot channel. The reference bank `output/templates_spline/0p73GeV_H2/` has the same defect (χ²/ndf 3375, I_int 4.52, lohi_iint_ratio 0.65; the option (i) build reused the module template), so R1 and R2 carry it. |
| In-sample gains and maps | The pooled gains are fitted on the splits they are applied to. `pos_all` uses the all-split map in-sample; the held-out halves (§4.6) give the same size at 3.49 GeV. |
| pm f(A) rejection | The reading in §3.5 (template-noise imprint, Q_raw not an unbiased truth) is a hypothesis; not tested. |
| Gains not refit | After the position correction and per clustering variant (first-order argument: the map has mean 1; member-sum changes ≲ 0.3 % absorbed by the seed-indexed step). |
| Ring dependence | Not corrected: one map per run type; ring and side differences up to 0.8–1.7 % in single bins (§4.4). |
| ε_peak with the position correction | Peak-yield ratio 0.986 at 0.73 GeV and 0.994 at 3.49 GeV in the full scan (§4.7; 0.985 / 0.994 held-out); not investigated. |
| Veto screening | 40k events per split, ≈ 27 calibrated modules at 0.73 GeV; its Gaussian values are unstable. Only the full-statistics rows support the recommendation (§7). |
| Decomposition | Two estimators disagree on b (0.98 % CB for all vs a negative b² with the preferred set); errors are SEM-only (§8). |
| ±2σ Gaussian at 0.73 GeV | Non-convergence within 5 refits in 5 of 16 post spectra per path (§2). |
| Mixed stored Gaussian columns | Scans started across the 19:50 install store fixed150 (jobs before) and ±2σ (jobs after) Gaussian columns: the rest of `pos_odd`, all of `pos_all`, the rest of the veto matrix. Never used in tables: every Gaussian number is refit from the stored histograms with `compare_scan_dirs.py --gauss-window 2sigma`. |
| Harness resume | `job_done` compares the command string, pre-window, W, pulse and `bench_args` (`gate_scan_20261022.py:828-858`); it does not check the content (sha1) of files the command names. A JSON or gain table rewritten under the same path is not detected. |
| Ta | Maps 11/144 bins populated; never fitted; report only. |
| Deferred | Full-run LMS files on ifarm (§5); option (ii) in the matched filter (needs a separate go); pile-up and efficiency tables were not redone (f(A) unchanged after the rejection). P4 variants: done (§6.2). |

## 10. How to run

Configs and run scripts are in `output/scan_20261025/configs/`: `cfgs_20261025.json` (`fa_percrystal`), `cfgs_20261025_p2.json` (`pos_*`), `cfgs_p34_withpos.json` (`lms_2run`, `cc1`, `mme*`), `cfgs_veto.json`, `cfgs_vetofull.json`, `run_*.sh`. The copied run scripts set `S=` to this directory. The tools do not log argv; commands marked (r) are reconstructed from plan §13 and `--help`.

```bash
cd /Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts && source setup_env.sh
E=/Users/L00338853/work/PRad/data/evio; DB=$PRAD2_DATABASE_DIR/calibration
WP="--calibrate --exclude-rowcol 15 20 15 20 --pre-window 100 200 --inner-radius 200 --min-per-module 30 --window 16"
S=output/scan_20261025/configs     # configs and run_*.sh
T0=output/templates_spline; T1=output/templates_spline_20261025; PM=output/templates_spline_20261025_pm
G="python3 gate_scan_20261022.py"; OUT=output/scan_20261025; P="--out $OUT --parallel 12 --timeout-seconds 1800"

# (a) P0 G0: pre-edit copies (+ output symlink) in a scratch dir W, rawsum/mf diff, centroid vs dump; repeated after every benchmark edit
W=<scratch>; mkdir -p $W/before && cp benchmark_hycal_timing_test_20261021.py mf_pulse.py _common.py percrystal_calib.py spline_template_extract.py $W/before/ && ln -s $PWD/output $W/before/output
for M in rawsum mf; do for V in old new; do B=benchmark_hycal_timing_test_20261021.py; [ $V = old ] && B=$W/before/$B
  python3 $B $E/prad_025320.evio.00050 $W/g0_${V}_$M $WP --lms-corr auto --hc-calib-file $DB/calibration_factor_2p2_fadc.json --max-events 20000 --pulse-method $M; done
  for f in $W/g0_old_$M*.tsv; do g=${f/g0_old/g0_new}; diff <(grep -v elapsed_s $f) <(grep -v elapsed_s $g) >/dev/null || echo "DIFF $g"; done; done
REF0="--pulse-method mf --mf-templates $T0/2p24GeV_H2/templates_2p24GeV_H2.npz --mf-scale 1.0 --gain-file output/percrystal_spl/025320/gains_025320_all.tsv"
python3 benchmark_hycal_timing_test_20261021.py $E/prad_025320.evio.00050 $W/g0_xy $WP --lms-corr auto \
  --hc-calib-file $DB/calibration_factor_2p2_fadc.json --max-events 20000 $REF0 --write-cluster-tsv --cluster-dump $W/g0_xy_dump.npz

# (b) banks: option (i) (build, ampcorr, check), option (iii) per-module copies (r)
for R in 0p73GeV_H2 2p24GeV_H2 2p24GeV_C 2p24GeV_Ta_X17 3p49GeV_H2; do
  python3 spline_template_extract.py build --pulses "$T0/$R/pulses/*.npz" --run-type $R --bands 100,250,500,1000,2000,4000 --out $T1/$R/templates_$R; done
python3 spline_template_extract.py ampcorr --out $T1/ampcorr.json --templates $T1/*/templates_*.npz
python3 spline_template_extract.py ampcorr --out $PM/ampcorr.json --templates $T1/*/templates_*.npz --per-module-f $PM
for B in $T1 $PM; do for R in 0p73GeV_H2 2p24GeV_H2 2p24GeV_C 2p24GeV_Ta_X17 3p49GeV_H2; do
  python3 spline_template_extract.py check --pulses "$T0/$R/pulses/*.npz" --templates $B/$R/templates_$R.npz --amp-bins 100,250,500,1000,2000,4000 --out $B/$R/check/; done; done

# (c) study dumps (7 x 20k events, one arm per bank set) and intraw / option (ii) diagnostic
zsh $S/run_fa_dumps.sh fa $T1; zsh $S/run_fa_dumps.sh fa_pm $PM     # -> output/mf_fit_study_20261025/dump/<arm>/
FS=output/mf_fit_study_20261023/uniform; MS=output/mf_fit_study_20261025                                       # (r)
python3 mf_fit_study.py intraw --root $FS --arm sampled --bins 100,250,500,1000,2000,4000 --out $MS/diag_sampled_acc   # fine: 100,250,500,1000,1500,2000,3000,4000
python3 mf_fit_study.py intraw --root $FS --arm fa_pm --arm-dir $MS/dump/fa_pm --bank "$PM/{runtype}/templates_{runtype}.npz" \
  --bins 100,250,500,1000,2000,4000 --out $MS/fa_pm_acc
python3 mf_ampshape_diag.py --dumps "$MS/dump/fa/*/pulses.npz" --templates "$T1/{runtype}/templates_{runtype}.npz" --out $T1/ampshape_diag/

# (d) production of the pm f(A): dump pass, 24386.00060 extra dump, 14 gain fits, driver scan, compare, module maps
zsh $S/run_fa_dumppass.sh 8       # harness tag fa_percrystal, stages R and RX, + the extra dump
zsh $S/run_fa_prod.sh             # percrystal_calib.py fit --splits all --fixed-windows; driver --mf-templates-map "$PM/{runtype}/..." --gain-dir output/percrystal_fa
# R1 (r): driver rerun of the reference configuration with the centroid columns
python3 benchmark_hycal_timing_scan.py --pulse-method mf --window 16 --pre-window 100 200 --mf-templates-map "$T0/{runtype}/templates_{runtype}.npz" \
  --mf-scale 1.0 --gain-file-map auto --gain-dir output/percrystal_spl --write-cluster-tsv --parallel 8 --timeout-seconds 1800 --output-dir output/scan_20261025_ref

# (e) position maps, fit, held-out, full scan (R1 = output/scan_20261025_ref)
R1=output/scan_20261025_ref; PC=output/poscorr_20261025
for SET in all even odd; do MPLBACKEND=Agg python3 position_corr.py map --scan-dir $R1 --splits $SET --gauss-window fixed150 --out $PC/; done  # as made (pre-§2)
MPLBACKEND=Agg python3 position_corr.py map --scan-dir $R1 --splits all --out $PC/before_2sigma/                                          # ±2σ view
for SL in 0p73GeV_H2 2p24GeV_H2 2p24GeV_C 3p49GeV_H2; do for SET in all even odd; do
  python3 position_corr.py fit --map $PC/map_${SL}_$SET.tsv --form both --out $PC/$SL/poscorr_${SL}_$SET.json; done; done
for SET in all even odd; do python3 position_corr.py fit --global --maps "$PC/map_*_$SET.tsv" --out $PC/global_$SET/; done
mkdir -p $PC/2p24GeV_Ta_X17 && for SET in all even odd; do ln -sf ../2p24GeV_H2/poscorr_2p24GeV_H2_$SET.json $PC/2p24GeV_Ta_X17/poscorr_2p24GeV_Ta_X17_$SET.json; done
python3 $PC/form_compare_offline.py > $PC/form_compare_offline.tsv   # bilinear vs poly8, offline (r)
zsh $S/run_p2.sh                                        # pos_even (odd splits), pos_odd (even splits), compares, pos_all R+RX, compare, module maps, after-map
python3 $OUT/pooled_delta.py pos $OUT/pos_even/compare/heldout_splits.tsv $OUT/pos_odd/compare/heldout_splits.tsv > $OUT/heldout_pooled.md   # (r)

# (f) Gaussian definition: old-vs-new tables, and any comparison with an explicit window
for X in spline:scan_20261024_spline pc:scan_20261024_pc r1:scan_20261025_ref; do
  python3 regauss_tables.py --scan-dir output/${X#*:} --out output/regauss_20261025/${X%%:*}; done                               # (r)
python3 compare_scan_dirs.py --scans ref=$R1 pos=$OUT/pos_all --gauss-window 2sigma --out $OUT/pos_all/compare/pos_vs_r1      # fixed150 | stored for old tables

# (g) veto matrix (+ v30 rerun), LMS two-run check, veto full statistics, clustering variants
zsh $S/run_veto.sh; zsh $S/run_v30_p3.sh                # veto -> output/scan_20261025_veto/; lms_2run -> output/scan_20261025_p3/lms_2run/
python3 compare_scan_dirs.py --scans ref=$OUT/pos_all lms=output/scan_20261025_p3/lms_2run --out output/scan_20261025_p3/lms_2run/compare/lms_vs_r2   # (r)
python3 output/scan_20261025_veto/analysis/veto_table.py                                                                                         # (r)
zsh $S/run_vetofull.sh                                  # vf30, vf50p20 on the 16 full 0.73 GeV splits, cmpfull_* vs R1
(cd output/resolution_fit_20261025 && python3 phase0.py --also-fixed150 > out.txt)                                                              # §8
zsh $S/run_p4.sh cc1 mme2 mme0p5 mme5                   # all four run at full scale; R + RX per tag, compare vs R2 -> $OUT/<tag>/compare/<tag>_vs_r2
#   S_dev fallback (r): add --runs 025320,024503,024590,025308 (36 splits) to both harness calls in run_p4.sh
```
