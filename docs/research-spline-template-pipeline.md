# Research: per-run-type spline pulse-template pipeline (adapting Wassim's HyCal macros)

**Date:** 2026-10-23 · **Branch:** waveform_jz · **Status:** research only. This document contains no plan or code.

## Scope

This document collects the facts needed to plan a new pulse-template pipeline. The pipeline picks one "best" waveform per HyCal module for each run type and spline-interpolates it into a template. Those templates then feed the deconvolution and pile-up benchmarks (pile-up rate heat maps, 10 pile-up examples per run, a Table-4-style efficiency) and the energy-calibration harness.

The run types are:

| Beam energy | Target | Runs |
|---|---|---|
| 0.73 GeV | H2 | 24590, 24603 |
| 2.24 GeV | H2 | 25228, 25257, 25260, 25299, 25320 |
| 2.24 GeV | C | 25308 |
| 2.24 GeV | Ta | 25800, 26138, 26230 |
| 3.49 GeV | H2 | 24386, 24503, 24512 |

Out of scope: writing code, running scans, designing the plan.

Path abbreviations used below:

| Short | Full path |
|---|---|
| `WW/` | /Users/jingyi.zhou/research/PRad-II/Waveform_Wassim/ |
| `R/` | /Users/L00338853/work/PRad/prad2evviewer/ |
| `PY/` | `R/analysis/pyscripts/` |

---

## 1. Access (verified)

`WW/` exists and can be read from account L00338853. The files are owned by jingyi.zhou with mode `-rw-------`, but an ACL grants read access: `0: user:L00338853 inherited allow read,execute`. Every file in the directory was readable.

| File | Size | Date | What it is |
|---|---|---|---|
| best_ref_shape_fit.C | 35 444 | Aug 16 | PASS 1: picks the best reference waveform per module (§2) |
| hycal_mf_fit_premodule_visual.C | 38 519 | Aug 16 | PASS 2: per-module matched filter plus multi-pulse spline fit (§3) |
| fit(2).C | 54 933 | Aug 13 | Original NPS (Hall C) macro by Hamdi and Mazouz that PASS 2 derives from (§4) |
| Hycal_prad_wf_analysis.pdf | 300 533 | Aug 16 | Hamdi's 20-page write-up of PASS 1 and PASS 2 for PRad HyCal (§4) |
| Waveform_Analysis.pdf | 4 823 016 | Aug 13 | Hamdi's NPS collaboration-meeting slides, 18 Jul 2024 |
| 1-s2.0-S0168900226000495-main.pdf | 2 248 567 | Aug 13 | NPS resolution paper, NIM A 1086 (2026) 171323 |
| 1712.06041v2.pdf, PhysRevD.60.022002.pdf, An_introduction_to_matched_filters.pdf | — | Sep | Background papers on matched filtering (gravitational-wave and Turin 1960) |
| ref.txt | 1 362 | Aug 13 | One normalised reference pulse: 110 samples, peak = 1 at sample 44.5. Neither macro reads or writes it. |

**Input format (verified).** Both macros read the repo's own `replay_rawdata` output, `prad_<run>.<split>_raw.root` (`R/analysis/README.md:30`). They use tree `events` and the leaves `hycal.nch`, `hycal.module_id`, `hycal.module_type`, `hycal.nsamples`, `hycal.samples` (a flat array), `hycal.ped_mean` and `hycal.ped_rms` (best_ref_shape_fit.C:550-564).

`module_id` is filled from `Replay::moduleID` (`R/analysis/src/Replay.cpp:98-114, 511-517`). The encoding is G<n> → n and W<n> → n + 1000 (`R/analysis/include/GainCorrCompute.h:23`). `p.hycal.module_by_id` decodes it.

**Local raw ROOT files (verified).** `/Users/L00338853/work/PRad/data/raw_root/` holds only run24386: 11 splits, 9.5 GB. Every other target run exists locally only as EVIO in `/Users/L00338853/work/PRad/data/evio/`.

---

## 2. best_ref_shape_fit.C: PASS 1, best reference per module (verified, 1102 lines)

**Entry point.** `best_ref_shape_fit(input="prad_025308.00000_raw.root", output_root="hycal_references_best_timecut.root", output_ps, output_txt, max_events=10000, adc_to_mV=1.0, threshold_mV=2.0, min_rising_samples=4, min_descending_samples=10, background_samples=20, background_max_abs_mV=2.0, background_max_rms_mV=1.0, background_max_roughness_mV=1.0, second_peak_fraction=0.15, second_peak_exclusion=6, max_amplitude_mV=-1)` (L521-537).

Because `adc_to_mV=1` by default, the values called "mV" are actually ADC counts. Module keying uses `module_id` only, and the sample period is fixed at `SAMPLE_DT_NS=4.0` (L95).

**Candidate cuts** (`evaluateCandidate`, L372-489):
- **Pedestal:** the branch value `ped_mean` is subtracted (L194).
- **Height:** the peak must exceed 2.0 (L410).
- **Peak time:** must lie in 150–170 ns, i.e. samples about 38–42 (L96-97, L417-421). The window is hard-coded from the ~157 ns peak population of run 25308.
- **Pre-pulse room:** at least 20 samples must precede the first sample above threshold (L267).
- **Pulse shape:** at least 4 strictly ascending samples above threshold (L225-251), and at least 10 strictly descending samples after the peak (L314-331).
- **Pile-up veto:** reject the waveform if any local maximum outside ±6 samples of the peak exceeds both 2 and 0.15·A_main (L339-370).
- **Not applied:** there is no saturation or clipping cut (the max-amplitude cut is off by default), no χ² or shape cut, and no trigger or cluster requirement.

**Best-waveform choice.** Exactly one waveform per module is kept, with no averaging and no alignment (`betterCandidate`, L491-516, applied at L673). The ranking is:
1. A "strict background" candidate beats a fallback candidate. Strict means max|pre| ≤ 2, RMS ≤ 1 and roughness ≤ 1 over the 20 pre-pulse samples (L299-309, L448-455).
2. Higher amplitude wins.
3. Lower roughness wins.
4. Lower RMS wins.

**Template construction** (L870-946):
- The baseline is re-subtracted using the mean of the 20 pre-pulse samples.
- The waveform is divided by its maximum *sample* (not the spline maximum).
- A `TSpline3` (default end conditions, not "natural") is threaded through every sample, at x = bin centres 0.5, 1.5, … in sample units (bin centres at L924, TSpline3 built at L939-940).
- No resampling, smoothing, CFD, onset alignment or sub-sample alignment is applied.

**Outputs:**
- **ROOT file**, holding:
  - `refs/ref_<module_id>` (TH1D, nsamples bins over [0, nsamples));
  - `splines/spline_<module_id>` (TSpline3);
  - `reference_canvases/c_ref_<id>` (TCanvas);
  - the `ref_summary` TTree (module_id, module_type, nused, polarity, valid, …; L751-788);
  - TNamed and TParameter objects recording the cuts.
- **Plot file:** a `.ps` with one page per module (data, spline, the 11-point MF window).
- **Text file:** a `.txt` table with one row per module (L792).

**Pitfalls when adapting (verified in source):**
- **Clipped pulses win.** "Highest amplitude wins" with no overflow check selects clipped pulses. The Python decoder exposes `pk.overflow`.
- **ADC plateaus rejected.** The strict `<` in the 10-falling-samples test rejects real pulses that have ADC plateaus.
- **Silent drops.** If `maxval<=0` (L880-887), the module is written with valid = 0 and gets no `.txt` row.
- **Wrong default.** `max_events` defaults to 10000, which contradicts the TNamed selection note ("scan all events"). The example call at L51 passes −1.
- **Hard-coded timing window.** The 150–170 ns window was tuned on run 25308 only.

## 3. hycal_mf_fit_premodule_visual.C: PASS 2, per-module MF and multi-pulse fit (verified, 1077 lines)

**Entry point and inputs:**
- The function is named `hycal_mf_fit_permodule(...)` (:411) and does not match the file name.
- Defaults: `reference_file="hycal_references_best.root"`, `max_events=10000`, `mf_nsigma=5`, `tspectrum_fraction=0.02`, `min_physics_amplitude=1000` ADC, `amplitude_bound_fraction=0.12` (:412-423).
- The reference-file default differs from PASS 1's output name (`..._timecut.root`).
- `loadReferences` (:288) reads `ref_summary` and the `refs/ref_<id>` TH1D. It ignores the TSpline3 and rebuilds a cubic interpolator from the histogram with `ROOT::Math::Interpolator(kCSPLINE)` on the integer sample grid x = k (:329, :343-345). That is half a sample (2 ns) off PASS 1's x = k + 0.5 knots, and it uses natural end conditions.
- After the last bin the template continues as an exponential tail, fitted log-linearly to the last 8 bins with τ between 1 and 100 samples (:347-392).

**Matched filter and peak finding:**
- The matched filter correlates with the module's 11-sample reference slice centred on the peak, normalised by ΣR² (:165-189).
- `TSpectrum::Search(σ=2, fraction 0.02)` finds candidates on the MF output (:685-687).
- A candidate is kept only if its MF value is at least 5·ped_rms (:702) and it sits at least 5 samples from either edge (:701).
- It must also pass the raw-waveform gate `hasFourAscendingSamples` (:207-252): the last three of the four rising samples must exceed the local baseline + 2.
- Peaks closer than 3 samples are merged; at most 6 pulses are kept (:722-733).

**Fit:**
- **Model:** `B + Σ_p A_p·R_m(x − t_p)`, fitted with TF1/Minuit chi2 over the *whole* waveform (`hwf->Fit(&ffit,"Q0R")`, :864).
- **Time limits:** t_p is held to the MF seed ±4 samples.
- **Amplitude limits:** A_p is held to the raw local height ±12%.
- **Baseline limits:** B is held to B15 ± 3·max(σ15, ped_rms, 0.5), where B15 is the mean of the first 15 samples (:635-639, :765-771).
- **Errors:** σ_i = √|S_i·4.096/2| / 4.096, floored at 1 (:255-267). This is inherited from NPS.
- **Quality:** χ²/ndf = TF1 chi2 / NDF (:867).
- **Skipped channels:** channels with peak − B15 < 1000 ADC are skipped (fit_status = −2, :661).

**Pile-up handling.** No explicit pile-up flag exists. Every accepted MF peak simply becomes another pulse in the fit. Downstream code can use `nfound > 1` or χ²/ndf as a flag.

**Outputs:**
- **ROOT file:** `hycal_mf_permodule.root`, tree `fit` (:522-543). Branches are event, channel, module_id, nfound, fit_status, and the arrays fit_amp[6], fit_time[6] (ns) and fit_integral[6], plus chi2ndf and the seeds.
- **Per-waveform visuals:** `visuals/module_<ID>/...` objects.
- **PDF:** `Plots_per_module/module_<ID>_fits_nps_errors.pdf` with 2 pages per fitted waveform: the waveform with fit and label box, then the MF trace with threshold and peak markers (:978-1029).

**Python port notes (inference):**
- The interpolator maps to `scipy.interpolate.CubicSpline`.
- TSpectrum maps to `scipy.signal.find_peaks`. TSpectrum's σ = 2 smoothing and relative "fraction" do not map one-to-one.
- TF1 with SetParLimits maps to `least_squares` with bounds.
- Tree reading maps to uproot.

## 4. Other useful material

- **Hycal_prad_wf_analysis.pdf** (Hamdi, 16 Aug 2026, 20 pp). It is the authoritative text for §2–3 and confirms the cuts above:
  - pp.4-7: PASS 1;
  - p.7: the spline is "only an interpolator of the measured reference";
  - pp.7-10: PASS 2;
  - p.11: the reference matters more than the error scale.

  It gives no efficiency, resolution or pile-up ROC numbers. **Reuse:** the cut definitions and rationale. The constants (150–170 ns, the 1000 ADC cut, thresholds of 2) were tuned on a single run (25308).
- **Waveform_Analysis.pdf** (NPS slides). Reference criteria (p.3), spline (p.4), pulse finding (p.5). NPS results: σ_t 0.4–0.6 ns, 33% better energy resolution. **Reuse:** background only.
- **NIM A 1086 (2026) 171323** (NPS paper, sec. 3.1). This is the citable reference for the method.
- **fit(2).C.** The NPS original: up to 12 pulses, amplitude limits 0.2–5× the seed, reference files `ref_wf_<i>.txt` in the ref.txt format, hard-coded NPS paths and geometry. **Reuse:** none directly.
- **Sibling directories:**
  - `/Users/jingyi.zhou/research/PRad-II/project_jinst_paper_plan.md` plans the template extractor, pile-up generator, NNLS deconvolution and pile-up ROC, so it is directly relevant context.
  - `time_study/` holds only `prad_025308.00000-00099_raw.root` (100 raw ROOT splits of the carbon run, no scripts). This is the only `_raw.root` source for run 25308 found on the machine.

## 5. Table 4 and the pile-up generator (verified)

**Table 4 itself.** Table 4 is at `R/docs/technical_notes/waveform_analysis/fit_pulse_template.md:508-522`, titled "WaveAnalyzer detection efficiency":

| ΔT (ns) | r=0.1 | r=0.3 | r=0.5 | r=1.0 | r=2.0 |
|---|---|---|---|---|---|
| 8, 16 | 0 | 0 | 0 | 0 | 0 |
| 20 | 0 | 6 | 2 | 16 | 14 |
| 24 | 0 | 4 | 8 | 100 | 100 |
| 32 | 0 | 8 | 8 | 100 | 100 |
| 48 | 8 | 8 | 100 | 100 | 100 |
| 64 | 8 | 6 | 100 | 100 | 100 |
| 96 | 6 | 8 | 100 | 100 | 100 |
| 128 | 16 | 18 | 100 | 100 | 100 |

All values are in %.

**What "detection" means.** The footnote (:523) defines detection as the percentage of events where `findPeaks` returned 2 or more peaks. The 6–18% at low ratio "are noise-driven false positives".

**Generation settings** (:495-506):
- Run 025308 (C), SSP_RawSum, PbWO4 two-τ template (τ_r 2.46, τ_f 23.97 ns).
- Grid: ΔT ∈ {8, 12, 16, 20, 24, 32, 48, 64, 96, 128} ns × ratio ∈ {0.1, 0.2, 0.3, 0.5, 1, 2}, 50 events per cell. The table omits ΔT = 12 and r = 0.2. The efficiency resolution is 2%.
- Generator options: `--height-min 500 --t0-min 22 --residual-veto-sigma 999`.
- The doc also lists `--model-err-floor 0.03`, which is not an option of the current generator CLI.
- Noise is the real base-waveform noise only: `--noise-model` defaults to `none` (`PY/pileup_generator.py:1518`).

**How efficiency is computed:**
- `PY/benchmark_deconv.py:268`: `detection_eff = n_wa_2peak/n_events`, with the check at :420.
- The count does not require a match to the injected time.
- **No fake rate is defined anywhere, and there is no injection-free control sample.**

**Generator** (`PY/pileup_generator.py`, 1826 lines, from cd120a5):
- **CLI** (:1453-1545): positional EVIO paths, `--template <pulse_templates JSON>` (required), `--out-dir`, `--material`, `--n-per-config`, `--dt-grid`, `--ratio-grid`, `--random-dt`, `--height-min/max`, `--t0-min`, `--chi2-max`, `--noise-model {none,gaussian-iid}`, `--seed`, among others.
- **Shape:** analytic only. `load_template_shape` (:338-395) takes the median τ_r, τ_f and t0 of `_by_type[material]`. `eval_template` (:402-420) evaluates the two-τ form.
- **Injection:** `inject_pulse` (:424-480) adds `height·shape` onto a real single-pulse EVIO waveform.
- **Arbitrary sampled templates: not supported.** Support would need:
  - a sampled `TemplateShape` alternative (`TemplateShape` is defined at :164);
  - interpolation at `t − onset`;
  - `T_max` and `t_peak_offset_ns` taken from the sampled shape (:170-171, used at :454-455);
  - per-module keying, which is currently out of scope according to its docstring.
- **Outputs:** one `.npz` per cell with truth times, heights and onsets, plus `manifest.json` (`write_grid_point`, :965).

**Deconvolution and matched-filter benchmarks:**
- `benchmark_deconv.py` passes the template through `dec.PulseTemplateStore` (`_setup_template`, :118) into the C++ `WaveAnalyzer.deconvolve`. The Python binding exposes only τ_r, τ_f and is_global (`R/python/bind_dec.cpp:318-330`; `R/prad2dec/include/Fadc250Data.h:292`). **A sampled template cannot enter the C++ deconvolver without changing both the C++ code and the binding.**
- `benchmark_matched_filter.py:205` `build_kernel(tau_r, tau_f, clk_ns, kernel_samples=30)` is the only point where the shape enters that script, so a sampled kernel is a local swap there.

## 6. Current Python side (verified)

**`PY/mf_pulse.py` (992 lines, tracked; last changed in 5359d63):**
- **Shape model:** analytic two-τ throughout.
  - `template` (:64) gives a unit-peak shape that is 0 for t ≤ 0; `template_and_slope` (:70) adds the derivative for Gauss–Newton.
  - `load_h50` (:174) reads per-module (τ_r, τ_f) from `output/pulse_templates_h50.dat` (column 25320). Values outside the gates (1 ≤ τ_r ≤ 6, 15 ≤ τ_f ≤ 35 ns, or a fit marked not good) fall back to the PbWO4 median.
- **Grid:** `build_bank` (:329) builds kernels `K[m,φ,l] = p(4l+φ)`, with 4 phases per 4 ns sample (a 1 ns grid) and L = 24 (:277-285).
- **Time reference:** onset t0 is where the template becomes nonzero. The reported time is t0 + t50[m] (:784).
- **I_int:** the walk-out integral of the triangle-smoothed unit template (`walk_integral` :97, `unit_integral` :135), tabulated over phase (:296-311). Energy is E = A·I_int·κ (:793).
- **Interface a sampled bank must provide (inference from the code):**
  - per-module key, with an aggregate fallback;
  - unit peak, with t = 0 at onset;
  - an interpolant evaluable at arbitrary t, with a derivative;
  - T_max, t_pk and t50 measured from the sampled shape;
  - I_int obtained by walking the smoothed sampled shape. `unit_integral` already supports this once `template()` is swapped.
- **Code to change:** `template`, `template_and_slope`, `tmax_of`, `TemplateBank.__init__`, `build_iint`, `build_bank` and `load_h50`. Callers that use `bank.tau_r/tau_f/tmx` must also change: `_basis` (:403), `amp_at` (:434), `_gn` (:460), `_refine_multi` (:535), `_group` (:680) and `_synth` (:808).

**`PY/plot_pileup_examples.py` (520 lines, untracked):**
- **Selection:** trigger 0x100, unclipped PbWO4 channels, MF (nsigma 5, kmax 4).
- **Denominator:** E_dom > `--dominant-frac`·E_beam, default 0.5 (cut at :189, counted at :204, flag at :357). This is **not** "above 50 MeV".
- **Pile-up rule:** at least one other pulse with E/E_dom > 0.10, which matches the requested 10% rule. Pulses that are too early or at the edge are dropped.
- **Output:** one PDF per run, one page per candidate (data, components, residual, table).
- **Reusable for (ii)?** Yes, with three changes: loop over runs into one `PdfPages`, change the denominator line, and swap the template inside `fit_terms`/`plot_candidate`.

**Module positions for the heat map:**
- `p.hycal.module(i)` returns `.name/.id/.x/.y/.size_x/.size_y/.type` in mm, module centres (`R/python/bind_det.cpp:535-545`).
- The same information is in `R/database/hycal_map.json` and `hycal_module.txt`.
- An existing per-module patch heat map is `PY/plot_template_2d_map.py:105-137, 389-392, 533-536`. Modules are about 20.75 mm (W) and 38.15 mm (G).

**Existing extractor (`PY/fit_pulse_template.py`).** This script produced the h50 JSONs; `dump_pbwo4_table.py` then wrote them to the `.dat`.
- **Recorded settings:** height_min 50, height_rms_mult 10, chi2_max 3, min_pulses 50, SSP_RawSum trigger, pre/post 8/40 samples, max_modules_per_event 20, and t0_min 22 for 25320 and X17 (None for 25308). These are the values in the RawSum files' `_meta`; the no-cut files have max_modules_per_event 0 and no trigger.
- **Good-pulse cuts** (:1090-1117): exactly one peak from `wave_ana.analyze`, `quality==0`, `not overflow`, window inside the buffer, and height ≥ max(50, 10·rms).
- **Reusable pieces:**
  - `C.setup_pipeline` (`PY/_common.py:204`);
  - the `dec.EvChannel` read/decode loop;
  - `_resolve_trigger_event_types` (:136);
  - `hycal.module_by_daq`;
  - `wave_ana.analyze` (ped, rms, peak pos, height, quality, overflow);
  - the multiplicity pre-pass (:1010-1055).
- **Raw-pulse dumps:** with `--plot-dir` it writes per-channel `.npz` files of raw pulses, which are candidate inputs for a best-waveform picker.

## 7. ROOT I/O from Python (verified by running imports)

| Interpreter | PyROOT | uproot | prad2py |
|---|---|---|---|
| /Users/L00338853/Software/opt/miniforge3/bin/python (3.13, default) | no | 5.7.7 | yes |
| /opt/homebrew/bin/python3 (3.14) | 6.40.02 | no | yes |

**What uproot can write.** uproot 5.7.7 writers (`uproot.writing.identify.to_*`) cover TH1/2/3, TProfile, TTree, TObjString, TList and TArray. **They cover neither TGraph nor TSpline3.** gain_corr.root was read with uproot (`PY/benchmark_hycal_timing_test_20261021.py:663-672`).

**Storage recommendation (inference).** Write from the default env with uproot, using a layout that `hycal_mf_fit_permodule` already reads:
- one file per run type, e.g. `templates_<E>_<target>.root`;
- `refs/ref_<module_id>` as TH1D, the unit-peak sampled template, the same as PASS 1;
- a `ref_summary` TTree with module_id (Replay encoding), module_type, name, valid, polarity and nused. Add provenance (the source event) and the quantities mf_pulse.py needs (onset, t50, T_max, I_int);
- an explicit statement of the template grid's x origin (k or k + 0.5) and spline end condition (see Risk 6);
- a TObjString holding JSON metadata (run list, cuts, git hash).

How the two sides read it:
- **ROOT macros** rebuild a spline from the TH1D, exactly as PASS 2 already does at :343.
- **Python** rebuilds `CubicSpline` from the same arrays.

Also write a parallel `.npz` (array [n_mod, n_samp], module names, metadata) for fast Python loading. If a TSpline3 or TGraph is really needed, write it in a second step under homebrew PyROOT.

## 8. Run types and expected coverage (facts, then inference)

**Facts:**
- **Sparse readout.** The FADC readout is sparsified (`PY/plans/plan_20261023.md:339`, F11). About 7.1 (0.73 GeV), 15–16 (2.24 GeV H2 and C) and 19.9 (3.49 GeV) PbWO4 channels are read out per event.
- **Unread neighbours.** Of the 3×3 neighbours of a seed, 48% / 13% / 12% / 4% are unread (24590 / 25320 / 25308 / 24503).
- **Local data volume.** All 14 runs are on disk as EVIO splits (plan :92), each of about 160k events (plan :99). Split counts: 24590 10, 24603 6, 25228 11, 25257 7, 25260 8, 25299 8, 25320 8, 25308 8, 25800 8, 26138 4, 26230 13, 24386 11, 24503 10, 24512 7.
- **Module count.** HyCal has 1152 W and 576 G modules. W1141 is dead.
- **Existing h50 statistics** (`PY/output/*.json`, min_pulses 50, height ≥ 50, max 20 modules per event):

| Template file | Good W fits | Pulses per module | Lead glass |
|---|---|---|---|
| 25320, 12 splits | 1150 of 1151 | median 1306, min 55 | 4 of 43 have ≥ 50 pulses |
| 25308 (C) | 1078 of 1151 | median 384; 72 below 50 | 0 of 83 good |
| X17 / Ta combined | 377 of 1131 | median 16 | 15 listed |
| 25320 without the trigger and multiplicity cuts | — | — | 127 of 145 have ≥ 50 pulses |

**How strongly each cut filtered (verified from `_meta`):**
- The max-20-module cut removed almost nothing: `n_events_multiplicity_filtered` is 19 (25320), 17 (25308) and 25 (X17).
- The SSP_RawSum trigger cut barely thins 25320 (`n_events_trigger_filtered` 7 861 of 1.93M) or 25308 (3 995 of 1.31M). It does thin X17/Ta strongly: `n_events_after_trigger` = 173 299 of about 2.54M.
- The glass gap is unexplained by `_meta`. 25308 RawSum (t0_min None) has 8 of 83 G with ≥ 50 pulses, against 135 of 150 in the no-cut `025308_h50.json`, yet neither the trigger nor the multiplicity cut removed many events. The no-cut files also cap pulses at 500 per channel (`PY/fit_pulse_template.py:746`), and that cap is not recorded in `_meta`.

**Inference: modules with at least one waveform passing the strict PASS-1 cuts, per run type.** PASS 1 applies no trigger cut, so it sees every read-out channel.

| Run type | Expected W coverage | Notes |
|---|---|---|
| 2.24 GeV H2 | about 1100–1150 W | |
| 3.49 GeV H2 | about 1100–1150 W | |
| 2.24 GeV C | about 1050–1080 W | |
| 0.73 GeV H2 | perhaps 900–1000 W with any waveform, fewer passing the strict cuts | Few channels read out per event |
| 2.24 GeV Ta | likely better than the X17 JSON's one third, because that file kept only about 7% of events after the trigger cut | It is unknown how many of the non-RawSum events carry physics-like pulses in the 150–170 ns window |

The no-cut 25320 and 25308 files show that glass waveforms exist in the data; why the RawSum files lose them is not established. The likely gaps are the outer and corner PbWO4 modules, the W/G transition, and W1141 (fallback: Open question 2).

Because the selector keeps the single highest-amplitude waveform, a module with a handful of candidates gets an unrepresentative template. This matters most at 0.73 GeV and for whichever trigger selection is used for Ta.

---

## Risks

1. **Input format gap.**
   - Wassim's macros need `_raw.root`, which is mostly absent locally (see §1 and §4), while the Python stack reads EVIO.
   - Options: run `replay_rawdata` for each run, or port PASS 1 to Python on top of the `fit_pulse_template.py` EVIO loop.
2. **Saturation.** PASS 1's highest-amplitude ranking has no overflow cut, so it can pick clipped pulses.
3. **Run-specific timing window.**
   - The 150–170 ns window was tuned on 25308.
   - The Python side uses `t0_min 22` ns.
   - Peak timing per run type needs checking before the window is reused.
4. **The C++ deconvolver cannot take sampled templates.** This blocks the "deconvolution test with new templates" unless the C++ code and binding change, or the test runs on the Python MF/fit path instead.
5. **The Table-4 definition is weak.** "≥ 2 peaks found" has no time matching and no fake-rate control. A faithful "like Table 4" would reproduce that definition (see Open question 5).
6. **Normalisation and time-reference mismatch.**
   - PASS 1 normalises to the maximum sample and uses no onset.
   - `mf_pulse.py` needs onset t = 0, unit peak, t50 and I_int.
   - Converting needs an onset definition on the spline, such as the start of the rise or a fixed CFD fraction.
   - The x-grid conventions also differ. PASS 1 puts knots at k + 0.5 with TSpline3 default ends, while PASS 2 uses k with natural ends. `scipy` `CubicSpline` defaults to not-a-knot. The stored file must state its grid origin and end condition.

## Open questions for the user

1. **Data source.** Run `replay_rawdata` to make `_raw.root` for all 14 runs so the ROOT macros can be used as-is, or reimplement PASS 1 in Python on EVIO?
2. **Fallback.** For modules without a passing waveform (glass, Ta, 0.73 GeV corners), which fallback is acceptable: a neighbour, the type average, or the same module from another run type?
3. **Lead glass.** Should glass get templates at all, given that the current analysis path is PbWO4-only?
4. **Single waveform or average.** Should "one best sample" stay strictly single, or is an aligned average of the top N clean waveforms acceptable to reduce noise?
5. **Efficiency definition for (iii).** Keep the Table 4 definition (≥ 2 peaks found), or add time-matching and a fake-rate column?
6. **Deconvolution path.** For the deconvolution test, extend the C++ `PulseTemplate` to sampled grids, or benchmark only the Python MF/fit path?
7. **Pile-up rate denominator.** Is "events above 50 MeV in that module" measured on the dominant pulse energy (A·I_int·κ·cal) with the new template?

## Coverage gaps of this research

- No EVIO scan was run, so the coverage figures are inferences from the existing h50 JSONs.
- The figures in the PDFs were not viewed; the text was extracted with pypdf.
- The exact command line behind the h50 JSONs was not found. Its parameters, including `max_modules_per_event: 20`, are recorded in each JSON's `_meta`.
