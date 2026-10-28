# PbWO4 energy resolution: ±2σ Gaussian, amplitude correction, position correction, LMS, clustering checks, Møller veto, module timing

**Date:** 2026-10-25 (work of 2026-10-24 15:10 onward)
**Branch:** `waveform_jz` (HEAD `3566071`, 10-25 14:14; later work uncommitted, §1.3)
**Plan:** `analysis/pyscripts/plans/plan_20261025.md` (Phases 0–5; "Decisions (user, 2026-10-24)"; Progress log)
**Previous notes:** `analysis/pyscripts/notes/PbWO4_energy_resolution_20261023.md` (per-crystal gains, centred rectangle, sampled templates), `analysis/pyscripts/notes/PbWO4_energy_resolution_summary.md` (state as of 2026-10-24)

Paths are relative to `analysis/pyscripts/` in `/Users/L00338853/work/PRad/prad2evviewer` unless they start with `database/`, `docs/`, `prad2det/` or `python/` (repo root), or are absolute. σ/E, preferred fit (CB for 3.49 GeV H2, C, Ta; Gaussian for 2.24 and 0.73 GeV H2), legacy/gated, "paired" (mean ± SE of per-split differences), Δ in pp and the run-type names are as in the summary note; legacy path unless noted. **Every Gaussian number in this note uses the ±2σ definition of §2** unless marked "fixed150". "Reference" = `output/scan_20261024_spline/` (banks `output/templates_spline/`, gains `output/percrystal_spl/`, κ = 1.0, centred rectangle). R1 = the reference configuration rerun with the centroid columns (§3.6); R2 = R1 + position correction (§4.8).

## 1. Status

Complete (`plans/plan_20261025.md` progress log to 10-25 01:42): all scans of P0–P4 and both add-ons are finished; nothing running at hand-off. Addendum 10-25 23:33, decision 10-26: the hybrid per-module position map of Study B was confirmed by a full scan and **adopted per run type as R3**: hybrid at 3.49 GeV, 2.24 GeV H2 and C, the bilinear R2 map at 0.73 GeV and Ta (0.73 GeV degrades with the hybrid; Ta has no own fits); no rescan (§4.9).

Addendum 10-25 14:30–23:00, module timing (§9): per-module t_module − t_ref of the PbWO4 matched-filter pulses against the veto counters, the event median, the seed module and the CEBAF RF, with six new scripts (§1.2). Outputs in `output/module_timing_20261025/` and `output/module_energy_20261025/`. The per-module timing resolution above 100 MeV is 0.18–0.30 ns (RF reference; diagnostic 10/25). Accidentals arrive in beam-bunch steps of 4.00755 ± 0.00013 ns, and no gate or veto-timing cut helps e-p or C. Open items in §10.

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
| Benchmark `benchmark_hycal_timing_test_20261021.py` | `--write-cluster-tsv` forces `rec_ext` and appends `cluster_x_mm, cluster_y_mm` (C++ log-weighted centroid) to both paths' rows; without `--gain-file` it now also writes the existing `# member_identity` header line. `--pos-corr FILE` (off): E → E/R(u, v) after `--gain-file`, before the window fit, both paths; header `# pos_corr=… sha1=… form=… norm=… n_applied=… n_clamped=…`. `--veto-second-cluster EMIN_MEV` (off), `--veto-dphi DEG` (refused without EMIN), `--moller-gap-nsigma X` (3; 0 re-admits the Møller-gap modules) (§7). `--gauss-window {2sigma,fixed150}` (**default 2sigma**: changes the reported Gaussian columns, §2). `--cluster-min-module-energy MEV`, `--cluster-corner-conn {0,1}`, `--cluster-non-linear-corr {0,1}` (default: config value; `# cluster_cfg` line only when one is given). Neighbourhood path passes the module index to `bank.f_amp` (no-op for 1-D banks). 10-25: `--pos-corr` also reads `poscorr_hybrid_v1` by schema; the header line gains `n_own=… own_modules=…` for hybrid maps only; poscorr_v1 output unchanged (§4.9). |
| `spline_template_extract.py` | `build --bands E0,…` (generalised npz keys `band_edges, shape_band, r_band, n_band, A_med_band`; old keys kept); `ampcorr` with n knots; `ampcorr --per-module-f OUTDIR` (copies with a per-module `ampcorr_f[module, ln A grid]`); `check --amp-bins`. Defaults reproduce the reference banks byte for byte (2.24 GeV C build, ampcorr on copies, check). |
| `mf_pulse.py` | `from_npz` reads both band layouts; 2-D `ampcorr_f`, `f_amp(amp, midx)`. Two-τ default untouched. |
| `mf_fit_study.py` | `intraw --bins/--bank/--arm/--arm-dir`, `diag --out`. |
| `compare_scan_dirs.py` | `--gauss-window {2sigma,fixed150,stored}` (default 2sigma: refit from the stored inner-ring histograms; `stored` reproduces existing tables byte for byte). |
| `plot_module_resolution.py` | `--gauss-window`. |
| New: `position_corr.py` | `map` (12 × 12 bins, harness even/odd split sets via `G.r_pairs`, bootstrap errors, ring/side tables, `--gauss-window`), `fit --form {bilinear,poly8,both}`, `fit --global` (§4). 10-25: `fit --per-module` writes a `poscorr_hybrid_v1` map (per-module poly8 for seeds with ≥ `--min-clusters` training clusters, the run-type poscorr_v1 embedded as fallback; §4.9). |
| New: `mf_ampshape_diag.py` | Option (ii) offline refit diagnostic (§3.4). |
| New: `regauss_tables.py` | Old-vs-new Gaussian table from a scan's stored histograms (§2). |
| New (10-25 timing): `module_timing_dump.py`, `plot_module_timing.py`, `plot_module_timing_vs_integral.py`, `plot_module_energy.py`, `rf_time_dump.py`, `plot_module_timing_rf.py` | Per-module timing dump, timing PDFs against the veto, event-median, seed and RF references, timing against the integral, per-module R2 energy spectra (§9.1). Read the EVIO, the dumps and the scans only. |
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
| Git, 10-25 23:40 | HEAD `3566071` (10-25 14:14, "energy_cali_benchmarking updates") committed the morning package: the benchmark, `compare_scan_dirs.py`, `mf_ampshape_diag.py`, `mf_fit_study.py`, `mf_pulse.py`, `plot_module_resolution.py`, `position_corr.py`, `regauss_tables.py`, `spline_template_extract.py`, this note, the summary and `plans/plan_20261025.md`. Modified since: the benchmark and `position_corr.py` (hybrid, §4.9), this note, the summary, the plan. Untracked: `module_timing_dump.py`, `plot_module_energy.py`, `plot_module_timing.py`, `plot_module_timing_rf.py`, `plot_module_timing_vs_integral.py`, `plot_template_overlay.py`, `rf_time_dump.py`. |

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
| f range | 0.878–1.231; 0.878–1.087 without W335 (§10). |
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

### 4.9 Hybrid per-module map (Study B form (d)): confirmation scan (10-25 22:30–23:33), adopted per run type as R3 (10-26)

**Configuration.** R2's command line with only the map changed: `--pos-corr output/poscorr_20261025/hybrid/{runtype}/poscorr_hybrid_{runtype}_all.json` (harness tag `pos_hybrid`, `output/scan_20261025/configs/{cfgs_hybrid.json,run_hybrid.sh}`). Stage R 93/93 and RX 25/25, no failures, 61 min at 12 jobs.

**Map.** `position_corr.py fit --per-module --min-clusters 3000` (schema `poscorr_hybrid_v1`):
- Seed modules with ≥ 3000 training clusters get their own poly8, c0[1 + c1u² + c2v² + c3u²v² + c4u⁴ + c5v⁴ + c6u + c7v], with u, v clamped to ±0.6.
- The estimator is Study B's: 6 × 6 iterative ±2σ truncated means (≥ 30 per bin, ≥ 10 bins), weighted LSQ, normalised to the module's event mean, range check 0.9–1.1.
- Training: r = E_cal/E′ in [0.85, 1.12] from the R1 cluster TSVs (`output/scan_20261025_ref`, all splits), the input of the bilinear maps and of Study B.
- All other seeds use the run-type bilinear `_all.json`, embedded verbatim. Ta is fitted on its own clusters with the 2.24 GeV H2 fallback that R2 uses.
- Check: even/odd fits reproduce Study B's own-module sets and per-module n exactly. Coefficients differ by ≤ 0.035 (R by 0.01 % RMS) because the TSV centroids give exact u = 0 for single-module clusters.

| Run type | Seeding modules | ≥ 3000 | Own fits | Share of modules | Share of training events |
|---|---|---|---|---|---|
| 3.49 GeV H2 | 170 | 52 | 52 | 30.6 % | 78.3 % |
| 2.24 GeV H2 | 165 | 76 | 76 | 46.1 % | 87.9 % |
| 2.24 GeV C | 191 | 28 | 28 | 14.7 % | 55.6 % |
| 0.73 GeV H2 | 141 | 26 | 25 | 17.7 % | 55.6 % (W463 failed the range check, 0.969–1.183) |
| Ta X17 | 40 | 0 | 0 | 0 | 0 (22 356 training events in total) |

The all-split maps give own fits to about twice as many modules as Study B's half-sample maps (24–28 / 41–45 / 13). The header `n_own` counts all calibration events of both paths and is a different fraction: for example, 30 893 of 366 609 on 24590.00010.

**In-sample Δ vs R2** (paired, legacy; gated within 0.003 pp). Sources: `output/scan_20261025/pos_hybrid/compare/hybrid_vs_r2.{md,tsv}`, `_splits.tsv`, `hybrid_vs_r2_per_run.md` and `hybrid_vs_r2_pooled.md` (ε_peak is the peak-yield ratio cmp/ref).

| Run type (fit) | n | R2 | Hybrid | Δ preferred ± SE [SD] | Δ CB | Δ ±2σ Gauss | ε_yield | ε_peak |
|---|---|---|---|---|---|---|---|---|
| 3.49 GeV H2 (CB) | 27 | 1.610 | 1.492 | **−0.118 ± 0.007** [0.035] | −0.118 ± 0.007 | −0.103 ± 0.007 | 1.0000 | 0.9936 |
| 2.24 GeV H2 (Gauss) | 42 | 1.898 | 1.854 | **−0.044 ± 0.002** [0.011] | −0.053 ± 0.003 | −0.044 ± 0.002 | 1.0000 | 0.9974 |
| 2.24 GeV C (CB) | 8 | 1.831 | 1.797 | **−0.034 ± 0.005** [0.014] | −0.034 ± 0.005 | −0.029 ± 0.002 | 1.0000 | 0.9983 |
| 0.73 GeV H2 (Gauss) | 16 | 3.718 | 3.909 | **+0.191 ± 0.040** [0.158] | +0.127 ± 0.014 | +0.191 ± 0.040 | 1.0000 | 1.0192 |
| Ta X17 (CB), report only | 25 | 1.978 | 1.978 | 0.000 (identical) | 0.000 | 0.000 | 1.0000 | 1.0000 |

- Harness ε_yield ≥ 0.9993 and ε_peak ≥ 0.9999 in H2 and C.
- Per run (preferred): 3.49 GeV −0.134 ± 0.016 / −0.108 ± 0.007 / −0.111 ± 0.005 (24386 / 24503 / 24512); 2.24 GeV H2 −0.038 to −0.053 in all five runs; 0.73 GeV +0.135 ± 0.043 (24590) and +0.285 ± 0.065 (24603), with CB +0.134 ± 0.015 / +0.115 ± 0.031.

**Held-out Δ.** Study B protocol: test half evaluated offline with the map fitted on the other half, through the benchmark's own reader. Source: `output/offline_20261025/hybrid_heldout/heldout.{md,py}`, `heldout_splits.tsv`.

| Run type | Δ vs R2, preferred | Δ CB | Δ ±2σ Gauss | vs Study B (d) |
|---|---|---|---|---|
| 3.49 GeV H2 | −0.097 ± 0.007 | −0.097 ± 0.007 | −0.081 ± 0.007 | +0.001 ± 0.001 |
| 2.24 GeV H2 | −0.035 ± 0.002 | −0.042 ± 0.003 | −0.035 ± 0.002 | −0.000 ± 0.000 |
| 2.24 GeV C | −0.016 ± 0.004 | −0.016 ± 0.004 | −0.014 ± 0.002 | +0.001 ± 0.001 |

0.73 GeV and Ta were not part of Study B and have no held-out number.

**Per-module σ_M/μ_M**, median (N > 1000, `--group-runs`; `output/module_resolution_20261025_hybrid/summary/types_summary.tsv`), R2 → hybrid:
- 3.49 GeV H2: 1.582 → 1.519
- 2.24 GeV H2: 1.879 → 1.842
- C: 1.827 → 1.810
- 0.73 GeV H2: 3.497 → 3.505, with the 84th percentile 5.70 → 10.04
- Ta: unchanged (1.908)

**Decision (user, 2026-10-26): adopted as R3, per run type.** R3 uses the hybrid map at 3.49 GeV H2, 2.24 GeV H2 and C and keeps the bilinear map (`poscorr_v1`, identical to R2) at 0.73 GeV H2 and Ta: at 0.73 GeV the hybrid worsens σ/E (+0.191 ± 0.040 pp Gaussian, +0.127 ± 0.014 pp CB), and Ta has 0 of 40 own fits. R3 σ/E, legacy, preferred fit: 1.492 / 1.854 / 1.797 / 3.718 / 1.978 % at 3.49 GeV H2 / 2.24 GeV H2 / C / 0.73 GeV H2 / Ta (R2 1.610 / 1.898 / 1.831 / 3.718 / 1.978). The map failed the rule, which required improvements at 3.49 GeV H2, 2.24 GeV H2 and C consistent with Study B within about 2 SE, and no degradation by more than 1 SE at 0.73 GeV H2 or Ta:
- **0.73 GeV H2 fails:** +0.191 ± 0.040 pp (4.8 SE; CB +0.127 ± 0.014), in both runs.
- **In-sample gains exceed Study B (d) vs R2:** −0.118 vs −0.098 ± 0.008 (1.9 SE), −0.044 vs −0.035 ± 0.002 (3.2 SE), −0.034 vs −0.017 ± 0.005 (2.4 SE). The all-split maps cover twice as many modules and are applied to the splits they were fitted on. The held-out numbers above (−0.097 / −0.035 / −0.016) are the unbiased estimate and agree with Study B.
- **Ta carries no information:** no module reaches 3000 clusters, so the hybrid equals R2 by construction.

**Reading (hypothesis, not tested).** At 0.73 GeV the per-module linear terms scatter about 3× more across modules: SD of c6, c7 0.028 / 0.025, against 0.006–0.010 in the other run types. The 3000-cluster threshold chosen at 2.24–3.49 GeV gives noisy or tail-biased per-module fits at 0.73 GeV, where the peak is ≈ 3.7 % wide with a Møller tail.

**R3 configuration.** Every job uses only its own run type's map, so R3 exists without a new scan: `pos_hybrid` for the 3.49, 2.24 H2 and C splits, and `pos_all` for the 0.73 and Ta splits. Its in-sample and held-out Δ vs R2 are the rows above.

### 4.10 Best-50 spline templates (10-26): re-extraction, oscillation fix, validation. Not adopted; awaiting the user's decision

**Request (user).** Per module, keep only the 50 best pulses of a run type: highest energy below 1.2 E_beam, exactly one matched-filter (MF) pulse, at most 20 channels per event (LMS veto), then normalise, build the template, and iterate until it is stable and smooth. The production banks (`output/templates_spline/`) and the benchmark pipeline are unchanged. Outputs are in `output/templates_spline_best50_20261026/`, `output/percrystal_best50_20261026/` and `output/scan_20261026_best50/`.

**Configuration.**
- `select` (per split): production trigger (physics; physics-any for Ta), with `--mf-bank output/templates_spline/<slug>/templates_<slug>.npz --per-module-best 50 --emax-frac 1.2 --mf-strict --max-modules-per-event 20 --exclude-trigger-bits 0x1000000`.
  - E = MF A · I_int · calibration factor. The calibration file comes from `module_timing_dump.RUN_CFG`; C uses `--calib-file database/calibration/calibration_factor_2p2_fadc.json`, because 25308 is not in `RUN_CFG`.
- `build` (final, "(a) + escalation"): `--pool-order energy --n-max 50 --fit-points direct --shape-sigma 0.012 --lambda-select heldout --lambda-chi2-sigma full --lambda-min 1e2 --lambda-fallback 1e3 --converge-tol 0.01 --converge-df 5e-3 --max-iterations 10 --osc-escalate`.
  - Fit pool: the top 50 by energy of a module over all splits.
  - Held-out set: the next 13 by energy.
  - λ: chosen per module from {1e2, 1e3, 1e4} by the mean held-out χ²/dof. That χ² uses σ² = rms² + (A · 0.15 ns · f′)² + (A · 0.012 · f)². With fewer than 5 held-out pulses, λ = 1e3.
  - Escalation: a module template with any spurious extremum is rebuilt with λ × 10 until it has none. A spurious extremum is a maximum or minimum with prominence ≥ 0.2 % of the main peak, in 2–150 ns, more than 3 ns from the peak, where the peak is the largest value in 3–30 ns. A template still flagged at 1e4 gets the run-type average (fallback 2, `osc_fallback = 1` in `summary.tsv`).
  - Then `ampcorr` over the five banks and `check`.
- New `build` options, all off by default: `--lambda-select`, `--converge-df`, `--shape-sigma`, `--lambda-chi2-sigma`, `--lambda-min`, `--lambda-fallback`, `--osc-escalate`, and `--pool-order energy-phase` (comparison only). The default build is still byte-identical to the unmodified script on the Ta production candidates (npz arrays; `summary.tsv` identical to production).

**S1, selection (pass 1, MF reference = production bank).** 119/119 splits, 62 min at 12 parallel jobs. Multiplicity counts channels with A ≥ 50 ADC per event, before the quality cuts. The MF-strict fraction is mf_npk ≠ 1 among the survivors of the other cuts.

| Bank | Events (gated) | Multiplicity veto (> 20) | Mean / p99 multiplicity | Bit-24 veto | MF-strict rejected | Above 1.2 E_beam | Stored per module, median | Modules with ≥ 50 |
|---|---|---|---|---|---|---|---|---|
| 0.73 GeV H2 | 3 133 150 | 25 (0.001 %) | 2.30 / 5 | 0 | 15.6 % | 17 | 146 | 1091 |
| 2.24 GeV H2 | 6 702 649 | 70 (0.001 %) | 4.75 / 9 | 0 | 12.9 % | 323 | 1986 | 1151 |
| 2.24 GeV C | 1 309 752 | 17 (0.001 %) | 4.46 / 9 | 0 | 11.9 % | 101 | 205 | 1009 |
| Ta | 2 520 014 | 34 555 (1.37 %) | 10.18 / 22 | 2 | 21.0 % | 29 710 | 274 | 1097 |
| 3.49 GeV H2 | 3 117 324 | 136 (0.004 %) | 6.09 / 10 | 0 | 9.9 % | 0 | 582 | 1089 |

The multiplicity veto and bit 24 act almost only on Ta (physics-any trigger).

**Pass 1 oscillated: the λ grid {0 … 1e4} with noise-only σ.**
- The first build used the held-out λ grid {0, 1e1, 1e2, 1e3, 1e4}, the noise-only point σ (rms/A plus the timing term), the noise-only held-out χ² and `--lambda-fallback 1e2`. It is kept in `pass1_osc/`.
- 33–46 % of the fallback-0 templates showed spurious maxima at 0.5 % prominence, against 1.6–2.2 % in production (2.24 H2 and Ta).
- The held-out χ² was nearly flat across λ (median χ²/min 1.12 → 1.02 from λ = 0 to 1e4; χ²/dof ≈ 200–300 at high amplitude), and the split-half (even vs odd splits) max |Δf| was 8.6 % of peak where both halves chose λ ≤ 1e2, against 2.6 % at λ ≥ 1e3.
- The sub-sample-phase spread of I_int (`iint_phase_pp`) was 4.2–5.4 % at λ ≤ 1e2, against 1.6–2.0 % at λ ≥ 1e3 (fallback-0 medians, 2.24 H2 and Ta).
- The cause: the noise-only σ lets the largest pulses dominate the direct point fit, and a low λ leaves the spline free between the populated sub-sample phases.

**Grid (2.24 GeV H2 / Ta).** Common settings: shape σ 0.012, full-σ selection, fallback 1e3. The held-out χ²/dof (full σ) is computed on the same next-best-by-energy pulses for every bank (11 067 / 7415 pulses). "osc" = share of fallback-0 modules with a spurious maximum at 0.5 % prominence.

| Bank | osc | w99, ns | n(w99 > 2.5 ns) | Rise 10–90 %, ns | Held-out χ²/dof | Converged / hit cap |
|---|---|---|---|---|---|---|
| production | 1.6 / 2.2 % | 1.66 / 1.67 | 16 / 23 | 5.70 / 5.66 | 5.71 / 5.81 | – |
| pass 1 | 33.2 / 46.2 % | 1.18 / 1.05 | 38 / 29 | 5.75 / 5.92 | 4.36 / 5.42 | 275/876, 232/877 |
| (a) direct, λ ≥ 1e2 | 3.3 / 4.4 % | 1.76 / 1.68 | 22 / 22 | 5.54 / 5.52 | 4.22 / 5.33 | 778/373, 776/333 |
| (b) direct, λ ≥ 1e3 | 2.0 / 2.1 % | 1.80 / 1.76 | 23 / 27 | 5.50 / 5.49 | 4.25 / 5.34 | 796/355, 775/334 |
| (c) bins, λ ≥ 1e2 | 15.1 / 23.5 % | 1.50 / 1.38 | 39 / 41 | 5.57 / 5.54 | 4.39 / 5.47 | 471/680, 306/803 |
| (d) = (a) with the phase-stratified pool | 3.5 / 4.1 % | 1.76 / 1.69 | 24 / 21 | 5.53 / 5.52 | 4.15 / 5.32 | 778/373, 785/324 |

- No configuration met the first acceptance test: osc ≤ 5 %, plus w99 and n(w99 > 2.5) no larger than production.
- The user relaxed the w99 test to "median within +0.15 ns, n within +10" (the remaining width is the N = 50 λ effect) and chose (a).
- (d) is within noise of (a) and is not used.
- The user then asked for the per-module escalation, because some (a) modules still showed wiggles.

**Phase check (2.24 GeV H2, 806 modules with ≥ 60 eligible candidates; phase = mf_t0 mod 4 ns).**
- The candidates are near-uniform in phase: median Rayleigh R = 0.10.
- The top-50 pools have median R = 0.169, against 0.153 for random 50-pulse draws from the same candidates.
- In 11.2 % of modules the pool R exceeds the module's own 95th percentile of random draws (5 % expected). The circular-linear correlation between energy rank and phase has median 0.05.
- So the energy ranking induces a weak phase clustering. Most of the pool R is the statistics of 50 draws.
- **The MF energy depends on phase.** With the production bank as MF reference, E/⟨E⟩ in eight 0.5 ns phase bins (median over modules) is 1.0004, 1.0005, 0.9975, 0.9971, 0.9963, 1.0007, 1.0006, 1.0045. That is 0.8 % peak to peak, against 0.2 % for phase-shuffled data. Reading (not tested): a ±0.4 % bias of the MF energy itself; the sample is energy-selected (top 50 per split, MF-strict) (§10).

**S2, final banks: (a) + escalation, against production.**
- "flagged" = modules with a spurious extremum at 0.2 % before escalation.
- "osc 0.2 %" = maxima and minima at 0.2 % prominence over modules with an own template (after escalation / production).
- w99 and rise are medians over the fallback-0 modules.
- χ²_ho = held-out χ²/dof with full σ, on the common next-best set.
- Pool size: production median 58 / 200 / 84 / 60 / 188.

| Bank | Fallback 0/1/2 (production) | Flagged / fixed at +1 / +2 steps / osc_fallback | λ_final 1e2 / 1e3 / 1e4 | Converged / hit cap | osc 0.2 % | w99, ns | n(w99 > 2.5) | Rise, ns | iint_phase_pp | χ²_ho |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.73 GeV H2 | 607/481/64 (577/543/32) | 175 / 147 / 1 / 27 | 86 / 704 / 298 | 811/277 | 0 / 11.2 % | 1.67 (1.66) | 15 (12) | 5.46 (5.59) | 1.74 (1.71) % | 5.99 (6.14) |
| 2.24 GeV H2 | 1059/68/25 (1062/89/1) | 96 / 69 / 3 / 24 | 148 / 507 / 472 | 773/354 | 0 / 4.7 % | 1.77 (1.66) | 24 (16) | 5.53 (5.70) | 1.69 (1.95) % | 4.26 (5.71) |
| 2.24 GeV C | 753/281/118 (742/310/100) | 132 / 112 / 1 / 19 | 104 / 565 / 365 | 729/305 | 0 / 8.0 % | 1.73 (1.68) | 14 (14) | 5.51 (5.68) | 1.69 (1.77) % | 4.71 (5.84) |
| Ta | 840/249/63 (812/297/43) | 156 / 136 / 0 / 20 | 110 / 621 / 358 | 767/322 | 0 / 9.5 % | 1.71 (1.67) | 25 (23) | 5.52 (5.66) | 1.69 (1.71) % | 5.36 (5.81) |
| 3.49 GeV H2 | 869/196/87 (858/268/26) | 184 / 140 / 2 / 42 | 113 / 552 / 400 | 718/347 | 0 / 8.8 % | 1.75 (1.67) | 20 (9) | 5.58 (5.70) | 1.67 (1.85) % | 5.48 (6.17) |

- Fewer than 5 held-out pulses (flag 1, λ = 1e3) holds for 283–683 modules per bank.
- The `osc_fallback` modules are listed in each `summary.tsv` (column `osc_fallback`).

**Gates, final against production (`<slug>/check/`).**

| Gate | Result |
|---|---|
| G1-4 | PASS in all five, as production; sampled/two-τ 0.05–0.08 at 300–1000 ADC and 0.02 at 1000–2500 ADC. The 100–300 ADC bin is empty, because the held-out set is the next-best by energy. |
| G1-5 | PASS in all five; largest quartile 3.49 GeV −0.47 ± 0.03 % (limit 0.5 %; production +0.27 %). |
| G1-6 | **PASS in all five; production FAILs.** 0 modules with a phase-mean > 0.5 % (production 3–5); no-peak 0 (production 33–73). |
| G1-7 | FAIL as production (0.73 GeV N/A: only the 300–1000 ADC bin is filled). |
| G1-8 | PASS at 2.24 GeV H2 (1127) and C (1034). **FAIL at 3.49 GeV: 1065 < 1100**, from the 42 osc fallbacks (production 1126). |
| G1-9 | FAIL as production. lo/hi max-diff median 1.92–2.39 % (production 1.46–1.90 %), noise floor 1.41–1.53 % (production 1.09–1.52 %); I_hi/I_lo median 0.991–1.000 (production 0.980–0.990). |

`check` now applies the build's `fit_points` and `shape_sigma` in the G1-9 noise-floor rebuilds (review fix).

**S3, outer loop (pass 2 = `select` with the best-50 bank as MF reference, 2.24 GeV H2 and Ta, rebuilt as S2).**
- Pass 2 selected with the (a) bank before escalation. Escalation changes the flagged 8–17 % of own templates, by one or two smoothing steps.
- Split-half noise is the final configuration built separately from the even and odd split files.

| Quantity | 2.24 GeV H2 | Ta |
|---|---|---|
| Stored candidates in common, pass 1 / pass 2 | 95.3 % / 95.0 % | 94.4 % / 95.1 % |
| MF-strict rejection, pass 1 → pass 2 | 12.9 → 20.1 % | 21.0 → 22.9 % |
| max \|Δf\| per module, median / p90 / max (fallback 0 in both) | 0.48 % / 6.5 % / 44 % (1053 modules) | 0.47 % / 6.6 % / 13 % (834) |
| Δ rise, median (16–84 %) | 0.000 (−0.029, +0.035) ns | +0.000 (−0.032, +0.033) ns |
| Δ w99, median (16–84 %) | −0.001 (−0.11, +0.11) ns | +0.001 (−0.09, +0.11) ns |
| Δ held-out χ²/dof at the chosen λ, median (16–84 %) | +0.002 (−0.52, +0.54) | −0.002 (−0.52, +0.42) |
| Same λ chosen | 66 % | 67 % |
| Split-half: max \|Δf\| median / p90; \|Δ rise\|; \|Δ w99\| | 2.6 % / 8.8 %; 0.08 ns; 0.28 ns | 2.9 % / 8.1 %; 0.10 ns; 0.30 ns |

- **Verdict: converged.** The pass-1 → pass-2 changes have zero median and are 5–7 times smaller than the split-half noise.
- The 6.5 % p90 comes from modules whose final λ changes between passes: 27 % at 2.24 GeV H2 and 22 % at Ta. There, max |Δf| has median 2.0 / 6.1 % and p90 9.5 %. With the same λ the median is 0.27–0.29 % and p90 1.2 %. That is statistical wobble of the 13-pulse λ choice.
- Pass-2 gates: as pass 1, plus G1-7 PASS at 2.24 GeV H2.
- The MF-strict rejection rises with the best-50 bank as MF reference (2.24 GeV H2 +7.2 pp): the new templates find more secondary pulses.

**S4, energy resolution (4 splits per run type against R3).**
- Gains were re-derived with the final banks, following the production recipe (dump pass with `--cal-method percrystal --cluster-dump auto`, then `percrystal_calib.py fit --splits all --fixed-windows`). They cover all splits of the five runs that hold the chosen splits, into `output/percrystal_best50_20261026/<run>/`.
- All five fits converge in 2 iterations, with the same number of free modules as production (138–284), except Ta (138 vs 139).
- Gain ratio new/production (legacy, prior fraction < 0.5), median (MAD-RMS): 24503 1.0041 (0.63 %), 25320 1.0098 (0.49 %), 25308 1.0085 (0.39 %), 24590 1.0029 (0.33 %), 25800 1.0056 (0.39 %).
- Benchmark: the R3 `job.cmd` of each split, with only `--mf-templates` and `--gain-file` changed, into `output/scan_20261026_best50/scan_<run>/split_<n>/`.
- Before the comparison, the current benchmark code was rerun on one R3 split per type with the R3 inputs (`output/scan_20261026_best50/control/`). `result.tsv` is identical to R3 apart from `elapsed_s` in all five.
- Δ = new − R3, paired over the 4 splits, ±2σ Gaussian refit (`compare_scan_dirs.py --gauss-window 2sigma`).

| Run type (fit) | Splits | Δ per split, legacy (pp) | Legacy Δ, mean ± SE | Gated Δ, mean ± SE | σ/E R3 → new, legacy, 4-split mean (%) | Peak shift Δμ/μ, legacy (%) |
|---|---|---|---|---|---|---|
| 3.49 GeV H2 (CB) | 24503.20/40/60/80 | +0.003, +0.007, −0.002, +0.003 | +0.003 ± 0.002 | +0.007 ± 0.003 | 1.505 → 1.507 | −0.004 ± 0.005 |
| 2.24 GeV H2 (Gauss) | 25320.100/200/300/400 | +0.011, +0.008, −0.002, −0.011 | +0.001 ± 0.005 | +0.001 ± 0.005 | 1.824 → 1.826 | −0.005 ± 0.004 |
| 2.24 GeV C (CB) | 25308.100/200/300/400 | +0.019, +0.011, −0.017, −0.024 | −0.003 ± 0.011 | −0.004 ± 0.012 | 1.801 → 1.798 | −0.002 ± 0.005 |
| 0.73 GeV H2 (Gauss) | 24590.20/40/60/80 | +0.005, −0.028, −0.125, +0.015 | −0.033 ± 0.032 | −0.035 ± 0.031 | 3.706 → 3.673 | +0.002 ± 0.014 |
| Ta (CB) | 25800.200/400/600/800 | +0.026, −0.012, +0.028, +0.001 | +0.011 ± 0.010 | +0.024 ± 0.034 | 2.005 → 2.015 | +0.008 ± 0.009 |

- The peak shift is for the preferred fit, legacy path.
- The refitted gains absorb the template scale, so the elastic-peak positions are unchanged within 0.01 %.
- Files: `output/scan_20261026_best50/compare/best50_vs_r3{.md,.tsv,_splits.tsv}`.

**Outcome.**
- The best-50 banks with escalation have no spurious extrema at 0.2 % in any own template; production has them in 5–11 % of modules.
- They fit the high-energy held-out pulses better than production (χ²_ho −2 to −25 %; optimistic for the best-50 banks, see §10).
- They pass G1-6, which production fails.
- They are 0.01–0.11 ns wider at 99 % of the peak.
- They fail G1-8 at 3.49 GeV.
- The outer loop has converged.

- The energy resolution is unchanged within the 4-split precision: Δ σ/E = −0.033 … +0.011 pp legacy, each within 1.4 SE of zero. The gated path agrees, except 3.49 GeV, where it is +0.007 ± 0.003 pp (2.4 SE, worse).

**Not adopted; awaiting the user's decision.** The benchmark pipeline, `output/templates_spline/` and `output/percrystal_spl/` are unchanged.

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

## 9. Module timing (10-25 14:30–23:00)

Per-module time distributions of the PbWO4 matched-filter pulses against four references (veto counters, event median, seed module, CEBAF RF), on all 111 splits of the 13 H2 and Ta X17 runs. Outputs are in `output/module_timing_20261025/` (PDFs referenced, none embedded); `dump/` (4.9 GB) and `rf/` (582 MB) are kept out of git.

### 9.1 Scripts, data and references

| Script | Role |
|---|---|
| `module_timing_dump.py` | One npz per split (`dump/`, 111 splits): per event `ev_num` and trigger bits; per PbWO4 channel the dominant matched-filter pulse (largest Eu; time `p_time` = fitted onset + template t50, reference spline banks, R2 pulse settings, no per-module offsets; energy without LMS or per-crystal gains, used only for selection); every V1–V4 pulse with WaveAnalyzer CFD time, integral and height. No cuts. |
| `plot_module_timing.py` | Per-module t_module − t_ref histograms and fits. `--ref veto\|evt`, `--veto-align counter\|centroid\|type` (default counter, §9.2), `--veto-gate 10\|none`, `--time-calib`, `--min-n 10000`. Selection (trig & 0xf00) ≠ 0 and (trig & `TRIG_EXCL`) = 0, `TRIG_EXCL` = 0x3000000: LMS flasher and alpha pulser events leaked through the 0xf00 mask and made a 0.2 ns spike (diagnostic 10/25). |
| `plot_module_timing_vs_integral.py` | t_module − t_veto (X17: t_seed) against the fitted integral; `--veto-counter any\|V1…V4` (§9.4). |
| `plot_module_energy.py` | `output/module_energy_20261025/module_energy.pdf`: 123 modules, four run-type panels of calibrated R2 spectra (`output/scan_20261025/pos_all`) with the preferred fit; 98 / 123 / 111 / 8 panels with N > 1000 (0.73 / 2.24 / 3.49 GeV / Ta). |
| `rf_time_dump.py` | RF TDC hits per event → `rf/` (111 splits, 38 s wall). |
| `plot_module_timing_rf.py` | RF phase per module, joined to `dump/` by event number (§9.5). |

| Reference | Definition | Available |
|---|---|---|
| t_veto | V1–V4 pulse with integral > 2000 and \|Δ − D0\| < 10 ns, largest integral wins; D0 = peak of t_veto − t_event per run type | H2 only (§9.2) |
| t_event | Leave-one-out median of the other PbWO4 pulses with E > 10 MeV (≥ 3 pulses) | All |
| t_seed | Highest-energy PbWO4 pulse of the event (E > 10 MeV); on the seed's own page t_seed − t_2nd | All; used for X17 |
| t_RF | Nearest channel-A RF hit; t_module − t_RF folded into one bunch period | All |

Median σ over modules with N > 10 000, ns (modules in parentheses):

| Reference | 0.73 GeV H2 | 2.24 GeV H2 | 3.49 GeV H2 | Ta X17 | Source |
|---|---|---|---|---|---|
| t_veto, counter-aligned | 0.630 (93) | 0.622 (39) | – (0) | – | `module_timing_veto_summary.txt` |
| t_event, raw | 1.301 (203) | 0.975 (858) | 0.994 (472) | 1.961 (412) | `module_timing_evtref_summary.txt` |
| t_event, time-calibrated | 0.579 (203) | 0.535 (857) | 0.533 (472) | 0.785 (412) | `module_timing_evtref_tcal_summary.txt` |
| t_seed, raw | – | – | – | 2.640 (374) | `module_timing_veto_summary.txt` |
| RF phase, runs aligned | 0.500 (205) | 0.461 (871) | 0.583 (473) | 0.588 (465) | `module_timing_rf.tsv` |

- Event-median outliers (σ > 2 × median, 2.24 and 3.49 GeV): W290–W292, W318–W326, W351, W353–W357, W386–W387, W532, W557, W560, W591, W625. W464 has mean +9.0 ns at 3.49 GeV with σ 0.94 ns.
- The veto reference reaches N > 10 000 in 93 modules at 0.73 GeV and 39 at 2.24 GeV, none at 3.49 GeV (best W630, N = 2960) or X17.

### 9.2 Veto reference and counter offsets

V1–V4 are FADC crate 3 slot 3 ch 0–3, 50 × 50 mm at x = −630 to −730 mm (`database/hycal_map.json:1733-1736`). They are read out in the H2 runs; the X17 Ta data have no veto pulse. Fractions per physics event from `veto_summary.tsv` (D0 from `run_veto.log`):

| Run type | D0 (σ), ns | Any veto pulse | Integral > 2000 | In time | Runs used | Modules N > 10 000 |
|---|---|---|---|---|---|---|
| 0.73 GeV H2 | −33.42 (1.58) | 26.8 % | 25.9 % | 25.3 % | 24590, 24603 | 93 |
| 2.24 GeV H2 | −33.32 (0.85) | 21.6 % | 8.5 % | 1.4 % (3.0 % in 25299, 25320) | 25299, 25320; no peak in 25228, 25257, 25260 | 39 |
| 3.49 GeV H2 | −33.70 (1.28) | 8.8 % | 3.9 % | 0.37 % | 24503, 24512; 24386 peaks at +18.49 ns | 0 |
| Ta X17 | – | 0 | 0 | 0 | none | 0 |

At 3.49 GeV the veto pulse is in time in only 11 % of the events with a pulse above 2000 (24503 + 24512; 98 % at 0.73 GeV): 11.5 k in-time events over the three runs. This is why the veto reference has no statistics there.

**Counter offsets.** The four counters have no entries in `database/hycal_time_offsets/time_calib.json` and differ by cable delays. Picking the largest-integral counter without alignment made a ≈ 2.4 ns double peak (W492: V1 at 36.5 ns, 42 %; V4 at 34.1 ns, 53 %; W457 and W526 similar; diagnostic 10/25). Run drift is excluded (24590 vs 24603 and 25299 vs 25320 agree to 0.1 ns), as are 4 ns quantisation and pile-up (diagnostic 10/25).

`--veto-align counter` fits the per-(module, counter) peak positions as a_module + b_counter and shifts each counter by b_counter (`run_veto.log`, ns):

| Run type | V1 | V2 | V3 | V4 | V1 − V4 | Residual RMS (modules) |
|---|---|---|---|---|---|---|
| 0.73 GeV H2 | −0.91 | −0.74 | +0.52 | +1.45 | −2.36 | 0.069 (395) |
| 2.24 GeV H2 | −1.48 | −1.36 | −0.05 | +0.85 | −2.33 | 0.074 (532) |
| 3.49 GeV H2 | −0.82 | −0.67 | +0.49 | +1.31 | −2.13 | 0.115 (149) |

- `--veto-align centroid` moves the centroid of v_t − t_event per counter to D0. Without time calibration t_event carries the raw offsets of the HyCal region each counter covers, so this is wrong: W628 went from 0.75 to 1.65 ns (diagnostic 10/25). Kept for the record.
- Effect, σ(t_module − t_veto) unaligned → counter-aligned (`module_timing_veto_valign_type{.tsv,_summary.txt}` → `module_timing_veto{.tsv,_summary.txt}`): median 0.730 → 0.630 ns at 0.73 GeV and 1.016 → 0.622 ns at 2.24 GeV; W457 1.55 → 0.67 ns (0.73 GeV), W526 1.50 → 0.73 ns (2.24 GeV), W492 0.74 → 0.52 ns (0.73 GeV).
- V1 has two timing populations at 0.73 GeV (≈ −37.0 and −33.8 ns against the event), the early one dominating when W492 fires (diagnostic 10/25). Origin not established (§10).
- Regenerated with the alignment: `module_timing.pdf` (381 pages, 374 module pages including the X17 t_seed pages), `module_timing_vs_integral.pdf` and `_V4.pdf`, `module_timing_nogate.pdf` (375 pages each); same module order in all. The unaligned copy is `module_timing_valign_type.pdf`; `module_timing_tcal.pdf` was regenerated on 10-26 with the aligned code (`--ref veto --time-calib`; 381 pages; companions module_timing_veto_tcal.tsv, module_timing_veto_tcal_summary.txt, veto_summary_tcal.tsv, log run_veto_tcal.log). For the veto-referenced runs the database offsets change nothing: median σ 0.630 → 0.623 ns at 0.73 GeV (93 modules with N > 10⁴) and 0.622 → 0.618 ns at 2.24 GeV (39 modules), as expected since the module-referenced alignment already absorbs them. For X17 (t_module − t_seed, 374 modules) the median σ drops from 2.640 to 1.036 ns and the number of modules with |mean − median| > 3 ns from 74 to 0, i.e. the seed-referenced offset blocks of §9.5 are the database offsets themselves; 19 modules remain mean outliers (|mean − median| > 5 × MAD-σ, 0.60–1.26 ns).

### 9.3 Beam-bunch accidentals at 4.008 ns

- At 3.49 GeV the t_module − t_veto main peak has satellites at ±4 and ±8 ns: 10–20 % of the main peak each, pooled over 163 modules; W630 counts in ±0.3 ns are 458 at 0, 69 / 81 at −4 / +4 and 42 / 94 at −8 / +8 (diagnostic 10/25).
- The same events show no teeth in t_module − t_event, and t_veto − t_event alone has a flat accidental floor at 3.49 GeV.
- The no-gate all-pairs Δt (`--veto-gate none`) shows the accidental floor as a comb with 60–90 % modulation at all three energies (diagnostic 10/25).

**Period.** Tooth centroids for k = ±2 … ±24 over ±100 ns, fitted as kP + c (`accidental_comb_period.{png,txt}`); P_RF = 128.247 / 32 = 4.0077 ns from the RF TDC (§9.5):

| Set | Teeth | P (ns) | (P − P_RF)/σ | (P − 4.000)/σ | Tooth σ (ns) |
|---|---|---|---|---|---|
| 2.24 GeV H2 (25299 + 25320) | 46 | **4.00755 ± 0.00013** | −1.3 | +59.5 | 0.553 |
| 0.73 GeV H2 (24590 + 24603) | 45 | 4.00751 ± 0.00031 | −0.7 | +23.9 | 0.569 |
| 3.49 GeV H2 (24503 + 24512) | 46 | 4.00761 ± 0.00022 | −0.5 | +35.0 | 0.558 |

Per run, per counter and for the three 2.24 GeV runs without a veto peak, P lies within 2.1σ of P_RF. The tooth σ is module ⊕ veto ⊕ bunch length for pairs in different bunches.

**Reading.**
- Accidental veto and module hits come from neighbouring beam bunches. At 3.49 GeV the veto pulse above 2000 is accidental in 89 % of the events that have one, and the ±10 ns gate admits ±2 bunches: five teeth.
- The teeth do not appear against the event median. At 3.49 GeV the median is polluted by low-energy pulses from other bunches (−5 to +8 ns shoulders in t_module − t_event; diagnostic 10/25), and low-energy accidental module pulses have poor time resolution.
- The matched-filter time is continuous (no one-sample ambiguity), and `MIN_SEP_NS = 12` (`mf_pulse.py:62`) rules out a second MF pulse 4 ns away. The FADC sampling clock is not RF-locked: the period reads 4.0077 in FADC ns.

No-gate pairs, pooled (`module_timing_nogate.tsv`, `run_nogate.log`):

| Run type | Events | Pairs / event | Accidentals / ns / event | In-time fraction | S/B in ±10 ns |
|---|---|---|---|---|---|
| 0.73 GeV H2 | 3 133 153 | 2.19 | 2.4 × 10⁻⁴ | 0.96 | 448 |
| 2.24 GeV H2 (25299, 25320) | 2 568 914 | 1.84 | 3.3 × 10⁻³ | 0.28 | 8.1 |
| 3.49 GeV H2 (24503, 24512) | 2 658 931 | 0.77 | 1.8 × 10⁻³ | 0.066 | 1.4 |

**Consequence for the energy resolution.** In-time pile-up and accidentals arrive in 4.008 ns steps. A window narrower than one bunch spacing cuts genuine shower members: the two-module time difference has σ ≈ 0.8 ns even after time calibration (diagnostic 10/25). This agrees with the gate table in §9.6.

### 9.4 Integral dependence

Source `module_timing_vs_integral.{pdf,tsv}` (+ `_V4`). Horizontal axis: fitted integral Eu = A · I_int = m_E / cal_factor (ADC·samples).

- Flat above ≈ 500 ADC·samples (≈ 70 MeV with the median PbWO4 factors 0.139–0.142, `run_rf.log`). Slopes for Eu > 500 are −0.019 / −0.002 / −0.045 ± 0.38 ns per 1000 (0.73 / 2.24 / 3.49 GeV; diagnostic 10/25).
- Below 500 a walk-like rise of 0.3–0.9 ns (diagnostic 10/25); W492 rises by 0.86 ns at 0.73 GeV (`dt_high_minus_low` −0.865).
- Full-range linear slopes, median: −0.053 (93 modules) and −0.034 ns per 1000 (36) at 0.73 and 2.24 GeV with N ≥ 10⁴, −0.454 (28) at 3.49 GeV with N ≥ 2000. The low-energy step dominates them; a straight line is a poor model.
- Outliers at 0.73 GeV (`run_vs_integral.log`): W430, W500, W521, W568, W569, W693, slopes −0.79 to −1.17 ns per 1000.
- At 2.24 and 3.49 GeV few high-integral pulses coincide with a veto, so the upper axis is sparse.

### 9.5 X17 without veto: seed and RF references

**t_module − t_seed.** 374 modules, median σ 2.64 ns raw. 74 modules have means more than 3 ns from the common median (+0.57 ns), in blocks that follow the raw crate and cable offsets (`module_timing_veto_summary.txt`). Early: W219–W327, W353–W367, W560, W593, W598 and others (−2.4 to −7.8 ns). Late: W464, W559, W589–W590, W773–W775 and others (+3.6 to +6.8 ns).

**RF in the data.** The CEBAF RF divided by 32 is recorded in TDC ROC 0x40 slot 16 ch 0 and 8 (`database/daq_config.json:52`, `prad2dec/include/TdcData.h:40-43`, Python `RfTimeData`). Every physics event in X17 and H2 has RF hits (`rf_dump_run.log`), 5.61 channel-A hits per event, spaced 128.247 ns = 32 × 4.00771 ns (`rf_summary.tsv`). `docs/REPLAYED_DATA.md:201`, `test/evio_dump.cpp:685` and the comment at `TdcData.h:38` quote ≈ 131 ns; the measured spacing is 128.25 ns. Existing recipe: `docs/analysis_notes/rf_time_reconstruction_plan.md:52-61`; cluster-level `cl_dt_rf` in `prad2det/include/EventData.h:289-304`, offsets `database/hycal_rf_offsets/24840.json`.

Other timing sources:
- LMS reference PMTs (crate 3 slot 13 ch 10–12, `database/hycal_map.json:1730-1732`) are present in the X17 data (all LMS-trigger events), so LMS gain tracking is possible for X17.
- GEM: 6 × 25 ns APV25 samples and only a charge-weighted time of the seed strip (`prad2det/src/GemCluster.cpp:273-290`); no timing reference.
- TI, FADC and VTP timestamps: 4 ns ticks only.

**RF-referenced timing.** φ = t_module − t_RF folded into (−P/2, P/2], with t_RF the nearest channel-A hit (as `RfTimeData::nearest_a` / `ClusterDeltaRf`). P per split = mean channel-A spacing / 32 = 4.00771 ns in all 13 runs (`run_rf.log`). The join by event number matches 100 % (`rf_summary.tsv`). The fold is sharp: the trigger is synchronous enough that no 4 ns box appears. Runs are aligned onto the type's largest run before summing. Outputs: `module_timing_rf.pdf` and `module_timing_rf_tcal.pdf` (874 pages each), `module_timing_rf_vs_integral.{pdf,tsv}`, `module_timing_rf.tsv`, `rf_summary.tsv`, `run_rf.log`.

- Medians in the §9.1 table; run 24386 degrades 3.49 GeV.
- On the 93 veto modules at 0.73 GeV: median σ_RF 0.462 against σ_veto 0.630 ns, and σ_RF is smaller in each of the 93.
- **Above 100 MeV the per-module width is 0.18–0.30 ns** (pooled 0.21 ns at 2.24 GeV, 0.25 ns at X17; diagnostic 10/25). This is the intrinsic timing resolution.
- Spot values, 0.73 / 2.24 / 3.49 GeV / X17 (ns): W492 0.28 / 0.39 / 0.60 / 0.35; W457 0.50 / 0.66 / 0.89 / 0.53; W630 0.57 / 0.59 / 0.80 / 0.59; W527 0.52 / 0.46 / 0.65 / 0.50.
- Integral slopes of φ, median over N > 10 000 (ns per 1000 ADC·samples, `module_timing_rf_vs_integral.tsv`): −0.001 / −0.021 / −0.024 / **+0.105**; at X17 the 16–84 % range is 0.01–0.24. The X17 slope is largest on outer modules, W493 +0.3 ns over its range (diagnostic 10/25). The low-energy tail widens the ±2σ fit everywhere.

**Findings.**
1. The FADC–RF phase changes run to run, partly per FADC board. Between 25800 and 26230, 1139 of 1150 modules shift by more than 0.5 ns (`rf_summary.tsv`), some boards by ≈ P/2 (diagnostic 10/25). The split-to-split span within a run is 1–4 ps in 9 of 13 runs (24512 0.06, 26230 0.04, 24503 0.23 ns; 24386 below). RF time calibration must be per run.
2. 26 central modules have X17 time offsets 1.0–2.0 ns from the X17 median (`offset_tcal`, `module_timing_rf.tsv`) with normal widths (0.25–0.62 ns): W492–W495, W525–W529, W559, W560, W563, W564, W593, W594, W597, W598, W627–W631, W663–W665, W697. X17 needs per-run-period offsets.
3. **Run 24386 is anomalous.** Its RF phase drifts between splits (0.88 ns span; 175 of 728 modules differ by > 0.5 ns from reference run 24503, `rf_summary.tsv`). Within most splits the events sit 50/50 in two phases 1.33 ns apart, with no link to the TI timestamp (diagnostic 10/25). Its veto coincidence peak is at +18.49 ns (§9.2). Excluded from timing studies and flagged for the DAQ/run log.

Pulses from neighbouring bunches fold onto the same phase. The RF fold measures resolution and cannot separate pile-up; t_module − t_seed does that.

### 9.6 Gate provenance, timing cuts and the Ta-vs-W table

The working point (raw times, W = 16 ns, pre-window 100–200 ns) comes from the 10-21 scan on the rawsum path (`plans/plan_20261021.md` §8; 20261022 §7) and was confirmed on h50 mf. It was not re-selected on the spline banks or on R2; R1 and R2 ran at it (`WP` in §11). The Ta gate gain was measured with the 26138-only bank.

Δgate = gated − ungated σ/E, pp, preferred fit, half A unless noted. Rawsum rows: `output/gate_scan_20261021/gate_scan_runtype.tsv`; mf rows: `output/gate_scan_20261022/gate_scan_runtype.tsv`.

| Times, W (ns) | 3.49 GeV H2 | 2.24 GeV H2 | 2.24 GeV C | 0.73 GeV H2 | Ta X17 |
|---|---|---|---|---|---|
| raw, 4 | +0.36 | +0.60 | +0.39 | +0.64 | +0.27 |
| raw, 8 | +0.04 | +0.03 | +0.02 | +0.02 | −0.03 (B −0.16 ± 0.04) |
| raw, 16 (working point) | 0.00 | 0.00 | 0.00 | 0.00 | −0.05 (B −0.09 ± 0.05) |
| calibrated, 2 | +0.13 | +0.17 | +0.15 | +0.10 | +0.12 |
| calibrated, 3 | +0.06 | +0.09 | +0.08 | +0.08 | −0.02 |
| calibrated, 4 | +0.04 | +0.05 | +0.04 | +0.06 | −0.06 (B, pre-window 130–200: −0.00 ± 0.06) |
| calibrated, 6–12 | 0.00 to +0.03 | 0.00 to +0.01 | 0.00 to +0.01 | −0.02 to +0.02 | −0.05 to −0.07 |
| mf raw, 8 | 0.00 | +0.01 | 0.00 | +0.02 | −0.12 |
| mf raw, 16 | 0.00 | 0.00 | 0.00 | 0.00 | −0.16 (B −0.10 ± 0.02) |

- No gate helps e-p or C. Ta gains 0.05–0.16 pp, already in the working point. Calibrated windows of 2–3 ns cost 0.06–0.17 pp in H2 and C.
- A cut on t_module − t_veto or t_cluster − t_veto cannot help: the veto exists only in H2, is in time in 0.4–25 % of events (§9.2), and tags non-elastic particles.
- Correction of an earlier statement: no benchmark applies a veto-based timing cut or a "non-linearity" cut. R2 ran with the second-cluster veto off (§7) and the clusterer non-linearity correction at its default, on (§6.1).

### 9.7 RF time calibration per run (10-26)

Plan `plans/plan_rf_calib_20261026.md` (executed 10-26, S0–S7). Script `rf_time_calib.py`, output `output/rf_calib_20261026/`, full run 80 s with 12 jobs (`run.log`). Inputs: the 111 dump and RF splits of §9.1. All numbers below come from `rf_calib_validation.txt`, `rf_calib_runs.tsv`, `rf_calib_modules.tsv` and `run.log` of that directory.

**Method** (plan §3).
- Rows: the §9.5 rows, φ = fold(t_module − t_RF) per (run, module), with `RF:prep` and `RF:fit_unwrapped` unchanged (S0 only exposes the event mask `sel`). `m_fb` rows are dropped: 660 of 272 356 802 rows (the plan expected none).
- Offset o_{m,r}, three tiers: `e100` (E ≥ 100 MeV, N ≥ 300); else `allE` (all pulses, N ≥ 300, flag `allE`); else `pred` = fold(τ_m + Φ_r) from the `time_calib.json` prior τ and the run phase Φ_r; else unresolved (W1123, W1124 have no prior). Estimator: Gaussian core of the 80-bin circular histogram (0.050 ns) unwrapped at the smoothed argmax; 34 argmax fallbacks in 49 368 fits (0.07 %).
- ρ = fold(o − τ); Φ_r = circular median of ρ over the measured modules; ρ^r = fold(ρ − Φ_r).
- Branch of the continuous correction D = fold(o − τ − G_r) + kP: G_r = circular median of fold(o − τ − e) over modules, e = leave-one-out event-median reference (in-window mates, E > 10 MeV, ≥ 3 mates), k = argmin |D − e|, `firm` if the margin ≥ 1 ns. Pass B recomputes e with the firm D¹ applied to the mates (D²). It moves G_r by a run constant, 0.047–0.148 ns in the usable runs (`G1_ns`, `G2_ns`).
- Held-out: dhalf = fold(c_even − c_odd) over the even/odd split halves (≥ 150 rows each). Bimodal module: residual of the `allE` histogram after the Gaussian core ≥ 0.35 × amplitude beyond 0.75 ns; the run is unusable above 30 % flagged modules. Drift: span of the per-split pooled centres > 0.10 ns.
- JSON per run in the `24840.json` schema (`default` 0.0, `modules[name, offset_ns, …]`, extra keys ignored by `HyCalRfOffsets.h`), no NaN (`grep -c NaN` is 0 in all 13 files).

**Per run** (`rf_calib_runs.tsv`; σ = fit width of the `e100` tier, dhalf over the modules where it is defined):

| Run | Type | e100 | allE | pred | median σ (ns) | median \|dhalf\| (ns) | \|ρ^r\| < 0.3 | Quality, span (ns) |
|---|---|---|---|---|---|---|---|---|
| 24386 | 3.49 | 438 | 205 | 507 | 0.177 | 0.020 | 0.827 | **bimodal**, 0.993; 78.2 % bimodal modules, unusable |
| 24503 | 3.49 | 627 | 429 | 94 | 0.184 | 0.014 | 0.928 | **drift**, 0.133 |
| 24512 | 3.49 | 583 | 416 | 151 | 0.183 | 0.013 | 0.930 | ok, 0.044 |
| 24590 | 0.73 | 326 | 807 | 17 | 0.183 | 0.026 | 0.924 | ok, 0.002 |
| 24603 | 0.73 | 243 | 779 | 128 | 0.190 | 0.030 | 0.945 | ok, 0.003 |
| 25228 | 2.24 | 977 | 173 | 1 | 0.190 | 0.011 | 0.932 | ok, 0.003 |
| 25257 | 2.24 | 935 | 191 | 24 | 0.188 | 0.013 | 0.933 | ok, 0.003 |
| 25260 | 2.24 | 952 | 180 | 18 | 0.188 | 0.013 | 0.936 | ok, 0.003 |
| 25299 | 2.24 | 955 | 179 | 16 | 0.189 | 0.014 | 0.937 | ok, 0.001 |
| 25320 | 2.24 | 954 | 180 | 16 | 0.189 | 0.013 | 0.937 | ok, 0.006 |
| 25800 | X17 | 327 | 815 | 8 | 0.202 | 0.024 | 0.905 | ok, 0.006 |
| 26138 | X17 | 235 | 841 | 74 | 0.203 | 0.031 | 0.899 | ok, 0.003 |
| 26230 | X17 | 386 | 764 | 1 | 0.203 | 0.020 | 0.876 | ok, 0.004 |

- (i) Median σ per type (e100 tier, usable runs): 0.187 (0.73 GeV, 568 fits), 0.189 (2.24, 4773), 0.183 (3.49, 1210), 0.203 ns (X17, 948). Fraction with σ < 0.3 ns: 0.974 (26230) to 1.000 per run.
- (ii) Gate PASS: median \|dhalf\| 0.0167 ns, 0.9985 below 0.30 ns, over 13 011 of 13 254 measured (run, module) pairs with dhalf defined (usable runs). Counting undefined dhalf as failing the fraction is 0.9802 (not gated; the plan's §4.5 wording 'of measured' gives this reading and it is below 0.99; which denominator applies is open for the user). 24503 alone: 0.0138 / 0.9981, max 0.318.
- Fallback tier: `pred` counts 94 / 151 / 128 at 24503 / 24512 / 24603 and 74 at 26138 (the plan expected "a handful" outside 26138). Run 24512 has 999 measured modules (e100 + allE), one below the S2 acceptance of 1000.
- (iii) Prior consistency, \|ρ^r\| < 0.3 ns: 0.876–0.945 per usable run, below the expected ≳ 0.95. The `e100` tier alone is 0.969–0.984 (26230: 0.894), the `allE` tier 0.721–0.933 (the shortfall is mainly in the `allE` tier). Φ_r: +1.93 / +1.92 (24503 / 24512), −1.74 / −1.78 (0.73 GeV), −1.69 (2.24), −1.88 (25800), +0.079 (26138), +0.138 ns (26230).
- Walk bias (median c_allE − c_e100): −0.067 … +0.049 ns in the usable runs (24386: +0.467).

**Global 25800 → 26230 jump** ((iv), `rf_calib_validation.txt`). 1142 modules are measured in both; the circular median of ρ_25800 − ρ_26230 is +1.990 ns. Of 83 boards with ≥ 4 modules measured in both, 80 are within ±0.1 ns of it. The three outside (Δ_b = circular median of ρ per board):

| Board (crate, slot) | n | Δ_b(25800) | Δ_b(26230) | jump − global | mad_b 25800 / 26230 |
|---|---|---|---|---|---|
| (3,15) | 8 | +1.983 | −1.585 | +1.578 | 0.036 / 0.307 |
| (3,8) | 16 | +1.903 | −0.410 | +0.322 | 0.037 / 0.347 |
| (3,18) | 15 | +1.943 | −0.191 | +0.143 | 0.033 / 0.224 |

(3,15) and (3,8) were expected; (3,18) was not (+0.143). The block spans six boards; the three deviating boards are those where the shifted modules are about half or more of the measured modules (8/8, 8/16, 7/15); (3,6), (3,14), (3,5) are within ±0.1. The run-to-run phase change is therefore global, not per board (correction to "partly per FADC board" in §9.5 finding 1); a board deviates only where the 26230 central block is about half of its modules or more.

**The 33 modules** ((v), gate PASS): 26 central X17 modules (W492–W495, W525–W529, W559, W560, W563, W564, W593, W594, W597, W598, W627–W631, W663–W665, W697) and the 7 halo modules (W491, W496, W530, W632, W666, W661, W662). In 25320, 25800, 26138 and 26230 all 33 are measured and firm (margins 3.83–4.01 ns), max \|D² − e²\| 0.089 / 0.074 / 0.065 / 0.074 ns; no `branch_changed` among the pass-A-firm modules (all usable runs).

| Run | D² of the 33 (ns) |
|---|---|
| 2.24 GeV (25228–25320) | −0.51 … −0.18; \|D²\| ≤ 0.4 for 24 of 33 (expected 33); none beyond 0.7. The nine beyond 0.4: W493, W494, W495, W526, W527, W559, W629, W664, W697 (−0.41 … −0.51 at 25320) |
| 25800 / 26138 | 29 modules −0.31 … +0.04 (both runs); W525 −0.97 / −0.96, W528 −0.96 / −1.01, W597 −0.91 / −0.89, W697 −1.33 / −1.35 |
| 26230 | 26 central −2.10 … −0.85, 7 halo −0.76 … −0.34; W525 −1.39, W528 −2.00, W597 −1.58, W697 −1.10 |

- **The central block is shifted in 26230 only.** All 33 have \|ρ^r\| ≥ 0.68 ns there (min); in 25800 and 26138 the other 29 agree with the prior (D² −0.31 … +0.04), and the four exceptions follow. Six modules (W493, W494, W527, W528, W629, W664) have ρ^r = +1.57 … +2.00 ns at 26230, which the event median resolves to D² = −2.10 … −1.67 ns (k = −1); the prior alone would take the wrong branch.
- **W525, W528, W597 and W697 are early in all three X17 runs** (D² −0.89 … −1.35 ns at 25800 and 26138, −1.10 … −2.00 ns at 26230); at 2.24 GeV the same four are −0.26 … −0.51 ns (25228–25320).
- Against the plan's expectations: 26230 D² ∈ [−2.10, −0.34] (expected [−2.6, −0.4], upper end 0.06 ns outside); 25800 / 26138 as expected; the 2.24 GeV runs miss the expectation \|D²\| ≤ 0.4 for all 33 (24 of 33).

**24386 and 24503.** 24386 is flagged bimodal (78.2 % of measured modules, bimodal ratio 0.35, every other run ≤ 5.5 %); `24386.json` has `usable: false`, `quality: "bimodal"`, `modules: []`; its per-module rows stay in the TSV with `bimodal` / `run_unusable` and are excluded from every statistic. 24503 gets one run offset with `quality: "drift"`, `split_span_ns` 0.1333; it is the only `drift` run (the other 11 usable runs: span ≤ 0.044).

**(vi) Consistency and cross-check** (gate PASS). circ_mad of fold(o − τ − e): pass A 0.026–0.049 ns, pass B 0.020–0.040 ns per usable run (limit 0.15). The pooled tcal histogram fit reproduces `rf_summary.tsv:run_offset_tcal` in all 13 runs (max difference 0.00000 ns, limit 0.002).

**Caveats.**
- Time definition: the offsets centre the matched-filter `p_time` of the Python pipeline (`time_definition` in each JSON, with the template path, here `output/templates_spline/2p24GeV_Ta_X17/templates_2p24GeV_Ta_X17.npz` for 26230). `Replay.cpp:1093` applies `hycal_rf_offsets` to the WaveAnalyzer CFD time of `ClusterHit`: format-compatible, values not transferable without a constant per module type.
- The 26230-only shift of the central block may belong to the run (rate, DAQ) or to the X17 bank; a bank change invalidates the X17 files.
- `allE`-tier offsets carry the walk bias (−0.067 … +0.049 ns, both signs; above).
- Adoption into `database/hycal_rf_offsets/` is not done (plan Open Question 1).

**Output** (`output/rf_calib_20261026/`): `<run>.json` ×13, `rf_calib_modules.tsv`, `rf_calib_runs.tsv`, `rf_calib_validation.txt`, `rf_calib_<run>.pdf` ×13, `run.log`.

## 10. Caveats and open items

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
| RF time calibration | Done 10-26 (§9.7, `output/rf_calib_20261026/`, `plans/plan_rf_calib_20261026.md`): 12 usable runs, 24386 unusable. Caveats in §9.7: time definition (Python matched-filter `p_time`, not the `Replay.cpp` CFD time), template-bank dependence, `allE` walk bias; not adopted into `database/`. Correction: the run-to-run phase change is global, not per board (the three boards holding most of the 26230 central block, 23 of 33 modules, excepted). |
| Run 24386 | RF phase drift between splits, two phases 1.33 ns apart within splits, veto peak at +18.49 ns; excluded from timing studies; DAQ/run-log question (§9.5). |
| Veto hardware | V1 two timing populations at 0.73 GeV; 24386 peak at +18.49 ns; no coincidence peak in 25228, 25257, 25260 (§9.2). Run-log question. |
| RF period in the docs | `docs/REPLAYED_DATA.md:201`, `test/evio_dump.cpp:685` and `prad2dec/include/TdcData.h:38` quote ≈ 131 ns; the measured spacing is 128.25 ns. Fix open. |
| `module_timing_tcal.pdf` | Regenerated 10-26 on the aligned code (§9.2); the X17 seed-referenced blocks disappear with the database offsets. |
| Gate | Ta-only re-scan on the final pipeline deferred until more X17 data are available (§9.6). |
| R3 | Adopted 10-26 (§4.9). |
| Best-50 templates (§4.10) | Not adopted; awaiting the user's decision. Banks in `output/templates_spline_best50_20261026/` (pass 1, grid and the (a) banks kept in `pass1_osc/`, `grid/`, `a_noesc/`); gains for five runs only; 4 splits per run type benchmarked. |
| MF energy vs sub-sample phase | With the production bank as MF reference, E/⟨E⟩ varies by ±0.4 % over the sub-sample phase of mf_t0 (0.8 % peak to peak against 0.2 % for shuffled phases; 2.24 GeV H2 best-50 candidates, §4.10). Reading (not tested; measured on energy-selected candidates): a bias of the MF energy itself, which would also affect the production pipeline. Open: not measured with other banks, at other energies or in the cluster energy. |
| Best-50 held-out numbers | λ is chosen on the same next-best-by-energy pulses that `check` and the χ²_ho column hold out, so those held-out numbers are optimistic for the best-50 banks. The held-out set covers only the top of the energy range: the G1-4 and G1-7 bins below 300 ADC are empty, and the low-amplitude behaviour rests on f(A). |
| Best-50 3.49 GeV H2 | G1-8 FAIL (1065 < 1100 modules with an own template): the 42 modules still wiggly at λ = 1e4 use the run-type average. n(w99 > 2.5 ns) = 20 against 9. Recorded as caveats (user, 10-26). |
| Best-50 pass 2 | Selected with the (a) bank before escalation; escalation changes 8–17 % of the own templates (flagged modules) by one or two smoothing steps. |
| Best-50 S4 scope | Gains in-sample, from all splits of the five runs that hold the 20 benchmark splits; position maps were fitted with the production bank. Ta: 29 710 candidates above 1.2 E_beam (calibration `x17_25825_fadc.json`), against ≤ 323 elsewhere; not investigated. |

## 11. How to run

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
# hybrid per-module map (§4.9): fits (all; even/odd for the held-out check), scan + compare + module maps, held-out, per-run table
for SL in 3p49GeV_H2 2p24GeV_H2 2p24GeV_C 0p73GeV_H2; do python3 position_corr.py fit --per-module --scan-dir $R1 --splits all --slug $SL --fallback $PC/$SL/poscorr_${SL}_all.json --min-clusters 3000 --out $PC/hybrid/$SL/poscorr_hybrid_${SL}_all.json; done
python3 position_corr.py fit --per-module --scan-dir $R1 --splits all --slug 2p24GeV_Ta_X17 --fallback $PC/2p24GeV_Ta_X17/poscorr_2p24GeV_Ta_X17_all.json --out $PC/hybrid/2p24GeV_Ta_X17/poscorr_hybrid_2p24GeV_Ta_X17_all.json
#   even/odd: the same with --splits even|odd for 3p49GeV_H2, 2p24GeV_H2, 2p24GeV_C (fallback poscorr_<slug>_<set>.json)
zsh $S/run_hybrid.sh 12                                 # pos_hybrid R + RX, compare vs R2 (hybrid_vs_r2), module maps -> output/module_resolution_20261025_hybrid
(cd output/offline_20261025/hybrid_heldout && python3 heldout_hybrid.py 3)   # held-out, Study B protocol
python3 $OUT/pooled_delta.py hyb $OUT/pos_hybrid/compare/hybrid_vs_r2_splits.tsv > $OUT/pos_hybrid/compare/hybrid_vs_r2_pooled.md

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

# (h) module timing (§9): argv as recorded at the end of run_*.log; (r) = reconstructed from --help (dump_run.log, rf_dump_run.log, run_rf.log, run_veto_tcal.log and run_veto_valign_type.log record none)
M=output/module_timing_20261025
python3 module_timing_dump.py --all --parallel 12 --out $M/dump                                    # (r) 111 splits
MPLBACKEND=Agg python3 plot_module_timing.py --ref veto                                         # module_timing.pdf, _veto.tsv, _veto_summary.txt, veto_summary.tsv
python3 plot_module_timing.py --ref evt; python3 plot_module_timing.py --ref evt --time-calib   # module_timing_evtref{,_tcal}.*
MPLBACKEND=Agg python3 plot_module_timing.py --ref veto --veto-gate none                        # module_timing_nogate.{pdf,tsv}
python3 plot_module_timing.py --ref veto --veto-align type --time-calib                         # (r) superseded: module_timing_tcal.pdf is the 10-26 regeneration on the next line; with the current code this command writes *_tcal_typealign.*
MPLBACKEND=Agg python3 plot_module_timing.py --ref veto --time-calib                            # module_timing_tcal.pdf (+ _tcal companions), 76 s
python3 plot_module_timing.py --ref veto --veto-align type                                     # (r) now writes *_typealign.*; the on-disk *_valign_type.* files are the earlier unaligned run, renamed
MPLBACKEND=Agg python3 plot_module_timing_vs_integral.py                                        # module_timing_vs_integral.*
MPLBACKEND=Agg python3 plot_module_timing_vs_integral.py --veto-counter V4                      # (r) _V4.* (the log records only "V")
python3 plot_module_energy.py                                                                   # (r) defaults: --scan-dir output/scan_20261025/pos_all --out output/module_energy_20261025
python3 rf_time_dump.py --all --parallel 12                                                     # (r) -> $M/rf
python3 plot_module_timing_rf.py                                                                # (r) defaults; run_rf.log: min-e 0.0, min-n 10000, min-run-n 200
#   RF time calibration (§9.7): S0 acceptance (RF: outputs unchanged by the prep/pass1 edit; SC = scratch dir), then the full run
SC=<scratch>
MPLBACKEND=Agg python3 plot_module_timing_rf.py --runs 26138 --jobs 4 --out-dir $SC/rf_before   # before the edit
MPLBACKEND=Agg python3 plot_module_timing_rf.py --runs 26138 --jobs 4 --out-dir $SC/rf_after    # after
for f in module_timing_rf_test.tsv rf_summary_test.tsv module_timing_rf_vs_integral_test.tsv; do diff -q $SC/rf_before/$f $SC/rf_after/$f; done
MPLBACKEND=Agg python3 rf_time_calib.py --jobs 12 2>&1 | tee output/rf_calib_20261026/run.log        # 80 s, gates (ii), (v), (vi)
#   accidental_comb_period.{png,txt}: session scratch script comb_period.py, not kept in the repo; method in the header of the .txt
# (i) best-K template pool (10-26; options default off): per split, then per run type
python3 spline_template_extract.py select $E/prad_025320.evio.00050 --out <dir>/025320_00050.npz --trigger physics --per-module-best 50 --emax-frac 1.2 --mf-strict --mf-bank $T0/2p24GeV_H2/templates_2p24GeV_H2.npz --max-modules-per-event 20 --exclude-trigger-bits 0x1000000
python3 spline_template_extract.py build --pulses "<dir>/*.npz" --run-type 2p24GeV_H2 --pool-order energy --n-max 50 --out <dir>/templates_2p24GeV_H2   # + --fit-points direct, --converge-tol 0.01 optional
# (j) best-50 spline templates (§4.10; not adopted)
B=output/templates_spline_best50_20261026
SEL="--per-module-best 50 --emax-frac 1.2 --mf-strict --max-modules-per-event 20 --exclude-trigger-bits 0x1000000"
#   select on all 119 splits, 12 parallel; Ta: --trigger physics-any; C: + --calib-file database/calibration/calibration_factor_2p2_fadc.json (25308 is not in module_timing_dump.RUN_CFG)
python3 spline_template_extract.py select $E/prad_025320.evio.00050 --trigger physics --mf-bank $T0/2p24GeV_H2/templates_2p24GeV_H2.npz $SEL --out $B/2p24GeV_H2/pulses/025320_00050.npz
BX="--pool-order energy --n-max 50 --fit-points direct --shape-sigma 0.012 --lambda-select heldout --lambda-chi2-sigma full --lambda-min 1e2 --lambda-fallback 1e3 --converge-tol 0.01 --converge-df 5e-3 --max-iterations 10 --osc-escalate"
for R in 0p73GeV_H2 2p24GeV_H2 2p24GeV_C 2p24GeV_Ta_X17 3p49GeV_H2; do python3 spline_template_extract.py build --pulses "$B/$R/pulses/*.npz" --run-type $R $BX --out $B/$R/templates_$R; done    # 7-8 min each, 5 parallel
python3 spline_template_extract.py ampcorr --out $B/ampcorr.json --templates $B/{0p73GeV_H2,2p24GeV_H2,2p24GeV_C,2p24GeV_Ta_X17,3p49GeV_H2}/templates_*.npz
for R in 0p73GeV_H2 2p24GeV_H2 2p24GeV_C 2p24GeV_Ta_X17 3p49GeV_H2; do python3 spline_template_extract.py check --pulses "$B/$R/pulses/*.npz" --templates $B/$R/templates_$R.npz --out $B/$R/check/; done
#   pass 2 (2.24 H2, Ta): select as above with --mf-bank $B/a_noesc/<slug>/templates_<slug>.npz (the (a) bank, before escalation) --out $B/pass2/<slug>/pulses/...; build + check as above in $B/pass2/<slug>/ (no ampcorr)
#   pass 1 ($B/pass1_osc/): BX without --shape-sigma, --lambda-chi2-sigma, --lambda-min, --lambda-fallback, --osc-escalate; grid $B/grid/{a,b,c,d}/<slug>/: --no-pdf, (a) = BX without --osc-escalate, (b) --lambda-min 1e3, (c) --fit-points bins, (d) --pool-order energy-phase
# gains: the 10-24 dump job.cmd of every split of 24503, 25320, 25308, 24590, 25800 (output/scan_20261024/spl_percrystal/scan_<run>/split_*/job.cmd), output stem -> output/percrystal_best50_20261026/dumps/scan_<run>/split_<n>/result, --mf-templates -> $B/<slug>/templates_<slug>.npz; 44 jobs, 10 parallel, 26 min; then per run
python3 percrystal_calib.py fit --dumps "output/percrystal_best50_20261026/dumps/scan_025320/split_*/result_dump.npz" --splits all --fixed-windows --out output/percrystal_best50_20261026/025320/
# benchmark: the R3 job.cmd (output/scan_20261025/pos_hybrid for 3.49, 2.24 H2, C; pos_all for 0.73, Ta) with --mf-templates and --gain-file replaced -> output/scan_20261026_best50/scan_<run>/split_<n>/result; control/ = R3 inputs unchanged
python3 compare_scan_dirs.py --scans R3=output/resolution_fit_20261026/r3_splits best50=output/scan_20261026_best50 --gauss-window 2sigma --out output/scan_20261026_best50/compare/best50_vs_r3
```
