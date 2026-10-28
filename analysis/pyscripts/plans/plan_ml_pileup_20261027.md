# Plan: ML pile-up rejection and recovery on single-module PbWO4 waveforms

**User decisions (2026-10-27)**

- CPU torch + scikit-learn in a separate environment: approved (not needed for U0).
- Target: > 80 % detection efficiency at Δt = 8 ns, ratio 1 (two FADC samples; CEBAF bunch spacing 4.0077 ns). The Δt = 20 ns cell is no longer the headline cell.
- Grid for all units in RF-bucket multiples: Δt ∈ {4, 8, 12, 16, 20, 24, 32} ns.
- U0 is also the replacement for the asymmetric MF secondary gate (see output/templates_wassim_20261027/why_hard_cells.md): symmetric two-pulse χ² fit, no 12 ns minimum-separation floor.
- e-p run types first; X17/Ta per run later.

**Date:** 2026-10-27. **Status:** draft for user approval; nothing implemented.
**Reads first:** `notes/PbWO4_energy_resolution_summary.md` §3.8, §4; `notes/PbWO4_energy_resolution_20261027.md` §6.3; `output/templates_wassim_20261027/compare/pileup_3way{,_full}.md`; `mf_pulse.py` (`MFReconstructor`, lines 833–1243); `pileup_generator.py` (docstring); `benchmark_pileup_mf.py`; `best50_compare.py pileup`.
Paths are relative to `analysis/pyscripts/`. Standing constraints: Python only, no C++ change without approval; no commits; `database/` never written; Ta runs 25800/26138/26230 reported per run, never pooled; all outputs under `output/ml_pileup_20261027/`; no new evio data (dumped waveforms only); no package installs without approval (§8).

## 1. Problem statement

- Waveform: FADC 250 MHz, 4 ns/sample, 100 samples (400 ns), pulse peak at 150–170 ns; PbWO4 onset-to-peak 5.5–6.4 ns (≈ 1.5 samples; `pileup_generator.py` docstring), fall ≈ 25 ns (two-τ τ_f validity range 15–35 ns, `mf_pulse.py`). Beam bunches every 4.0077 ns; the FADC clock is not RF-locked (`notes/daq_queries_20261026.md`).
- Baseline: per-module sampled-template matched filter (MF) with subtract-and-repeat deconvolution, kmax 4, 12 ns minimum separation, r_model 0.3, joint Gauss–Newton refit (methods note §2.2).
- Current numbers (production pseudo-data, 60 cells Δt {8,12,16,20,24,32,48,64,96,128} ns × ratio {0.1,0.2,0.3,0.5,1,2}, 50 events per cell per run type, eff_match = both onsets within ±4 ns; `pileup_3way.md`):
  - Δt ≤ 16 ns: 0 % in every cell (the 12 ns separation floor and one-sample rise make the second pulse invisible to the filter); fake rate 0.7/event at Δt 16 ratio 1 (the filter puts one pulse between the two).
  - Δt 20 ns: ratio 1 → 29 / 41 / 63 % (R3 / best-50 / Wassim bank), ratio 2 → 99 %, ratio ≤ 0.5 → 0 %.
  - Δt 24 ns: ratio 1 → 82 / 91 / 88 %, ratio 0.5 → 2–3 %.
  - Δt ≥ 32 ns: ≥ 97 % for ratio ≥ 0.5; ratio ≤ 0.3 stays at 0 % until Δt 96–128 ns.
  - Fake rate 0.024–0.055 per event pooled per run type.
- Real pile-up rate (E_dom > 50 MeV channels, after burst veto): 0.15–0.41 % on e-p and C, 10.3–11.0 % on Ta (summary §3.8). An accepted secondary takes 6–8 % of the dominant energy. σ_E/E is template-invariant to within 0.02 pp at 3.49 GeV (within errors elsewhere) although the Δt 20 ns efficiency moves by 34 points (20261027 §6.4–6.5), and the time gate changes σ_E/E by −0.003 to +0.004 pp on e-p and C and by −0.048 / −0.030 / −0.155 pp on Ta runs 25800 / 26138 / 26230 (legacy → gated, 20261027 §7). This bounds what any pile-up method can do to σ_E/E (§7).

## 2. Goal and approach

Goal: decide, with a frozen test set and a sourced ceiling, whether a learned waveform model can (a) resolve two pulses at Δt 12–24 ns where the MF cannot and (b) return the first-pulse amplitude with no bias when a secondary overlaps; and quantify the σ_E/E consequence (expected ≈ 0, §7).

Approach, cheapest first: (1) establish the information ceiling with a non-ML two-pulse χ² fit on the same templates; (2) only if the ceiling is far above the MF, train small networks on synthetic overlays of real clean pulses; (3) validate on real data through RF-phase and Ta statistics. Rejected alternatives: training on simulated pulses (no pulse simulation exists; real clean pulses with real noise are available in millions); a C++ implementation (constraint); a GPU cloud run (no data may leave the machine without a decision, and the problem is small enough for CPU).

## 3. Candidate models

| # | Model | Role | Output | Why |
|---|---|---|---|---|
| M0 | Denser-grid MF + two-pulse template fit (scipy `least_squares`), likelihood-ratio test single vs double, onset grid 0.5 ns, no 12 ns floor | non-ML ceiling and baseline | N ∈ {1, 2}, A₁, A₂, t₁, t₂ | Any template-based method is bounded by this fit; if it fails at Δt 16–20 ns, so will a network trained on the same templates |
| M1 | 1D CNN classifier + regressor: 4–6 conv layers (kernel 5–7, 16–64 channels), 100-sample input (pedestal-subtracted, scaled by amplitude), heads: N ∈ {1,2}, (A_k, t_k) for k ≤ 2 | main candidate | as M0 | Standard for pulse-train parameter regression; ≈ 30–100 k parameters, CPU-trainable |
| M2 | Unrolled deconvolution / learned matched filter: 5–10 ISTA-like iterations with a learned per-module kernel bank initialised from the spline bank | interpretable alternative | sparse amplitude train on a 1 ns grid | Keeps the template physics; learns only the thresholding and the phase handling |
| M3 | Small transformer (2 layers, 4 heads, d 64) over 100 tokens | listed per brief; after M2 only | as M1 | ≈ 3–5× M1 on CPU (estimate, unsourced) |

Decision rule between M1–M3: the one with the highest eff_match at Δt 16–24 ns, ratio 0.5–1, at a fake rate ≤ the MF's 0.03/event and a first-pulse amplitude bias < 0.2 % in every cell.

## 4. Training and test data

- Clean pool: the best-50 `select` outputs, `output/templates_spline_best50_20261026/<slug>/pulses/<run>_<split>.npz` (`samples (N,100) uint16`, `ped`, `rms`, `mf_t0`, `mf_amp`, `mf_npk`, `mf_energy`, `name`); keep `mf_npk == 1`, `mf_sec_ratio < 0.05`, clip veto, 100–4000 ADC. Sizes: 16 / 42 / 8 / 25 / 28 split files (0.73 / 2.24 H2 / C / Ta / 3.49), ≈ 10⁴–4·10⁴ pulses per file.
- Overlays: `pileup_generator.py` with `--template <best-50 bank>` (the module's own spline), `--dt-grid 8..40 step 2` and `--ratio-grid 0.1,0.2,0.3,0.5,0.7,1,1.5,2,3`, `--n-per-config 2000`, plus `--random-dt --dt-min 4 --dt-max 200` for a continuous set. Phase: the fixed grid puts the injected pulse at the base pulse's sub-sample phase (Δt is a multiple of 4 ns); the ML sets use a random sub-sample phase for the injected pulse (generator flag to add: `--phase-jitter`, ±2 ns uniform; one Python change, ≈ 20 lines). Noise: the base pulse carries real noise; the injected template is noise-free, so the composite noise is the real one (no `--noise-model`). Single-pulse events (N = 1, 30 % of the set) come straight from the clean pool so the classifier learns the null; the generator injects one pulse per event and that stays.
- Split by split index (even/odd = parity of the split's 0-based rank in the run's sorted split list, so every run keeps both halves): train on even split files, validate on odd; the five production-synth base splits (24590.00010, 25308.00050, 25320.00050, 24503.00010) and all 25800 splits are excluded from the clean pool and from every overlay set; the production synth (`output/pileup_eff_20261024_fix/<slug>/synth`, 50 per cell) stays the frozen test set, never used in training. Five run types trained jointly with a run-type input; Ta per run in the reports.
- Real-data validation (what exists): (a) RF-folded phase of every recovered secondary (`rf_time_calib.py`, `output/module_timing_20261025/rf/`): a genuine beam-correlated secondary has phase consistent with the bunch structure (σ_t 0.18–0.30 ns per module above 100 MeV, summary §3.9), a fake has a flat phase; the flat fraction measures the real fake rate. (b) Ta runs (10–11 % channel pile-up) per run: recovered-secondary Δt spectrum must show the 4.0077 ns comb; the comb contrast is the real-data efficiency proxy. (c) The 14 burst events and the isolated +23 ns bumps (`output/pileup_20261024/artefact/`) as labelled negatives/positives. LMS pulses are not a usable truth (single pulses only). No labelled real two-pulse set exists; state this in the paper.

## 5. Metrics

- The existing 60-cell table and the dense grid (17 Δt × 9 ratios): eff_match (±4 ns both onsets), fake rate per event, n2, each with Wilson intervals (`benchmark_pileup_mf.py` conventions, `--window-rule 40:0.16` arm kept for comparability).
- Recovered first-pulse amplitude: bias and RMS of A₁/A₁,true per cell (this is what feeds σ_E/E); the same for single-pulse events (must stay at the MF's ±0.4 % phase-dependent level or better).
- Onset resolution of both pulses; ROC of the N = 2 decision at Δt 16, 20, 24 ns.
- Domain-shift controls: train 2.24 H2, test 3.49 H2 and Ta; train on 100–1000 ADC, test on 1000–4000 ADC; modules with fallback templates excluded from training, reported separately; a leave-one-run-out check on the real RF-phase fake fraction.
- σ_E/E: the R3 118-split configuration with the ML pulse list read through the `--pulse-list` hook of U5b on 4 splits per run type first (`best50_compare.py` style paired Δ), full 118 only if |Δ| > 2 SE on any H2/C type.

## 6. Steps (bounded units, each with an acceptance criterion)

- U0, cheapest unit that proves or kills the idea (minutes of compute, numpy/scipy only): M0 on the Δt 12–24 ns × ratio 0.3–2 cells of the production synth, 5 run types. Accept: if the χ² ceiling reaches ≥ 80 % at Δt 20 ratio 1 and ≥ 50 % at Δt 16 ratio 1 at fake ≤ 0.05/event, continue; if the ceiling is < 45 % at Δt 20 ratio 1, stop (report that the rise time, not the method, limits separation) and keep only U5. User target (2026-10-27): > 80 % at Δt = 8 ns, ratio 1; the 20 ns thresholds remain as the first screen.
- U1: dense training set with `pileup_generator.py` (+ `--phase-jitter`), even/odd split; manifest with sha1 of pools and bank. Accept: `validate_output` passes; cell counts ≥ 2000; phase histogram of injected onsets flat within 3σ.
- U2: M1 training (needs torch, §8), 3 seeds; model card with parameters, epochs, wall time. Accept: validation eff_match within 2 points of training; fake rate ≤ 0.03/event on single-pulse validation events.
- U3: M1 vs M0 vs MF on the frozen 60-cell test set and the dense grid; amplitude bias table. Accept: the paper numbers exist with intervals; decision rule of §3 applied.
- U4 (only if U3 shows a gain ≥ 10 points at Δt 20 ratio 1): M2, then M3; same tables.
- U5: real-data validation (§4): RF-phase fake fraction and Ta comb contrast for MF and for the best model, Ta per run. Accept: a flat-phase fraction with its uncertainty per run type.
- U5b: `--pulse-list FILE` input in the benchmark `benchmark_hycal_timing_test_20261021.py` (Python; the MF now runs inside it via `--mf-templates`). Accept: the MF's own list through the hook reproduces R3 on one split, TSVs identical except the header.
- U6: σ_E/E on 4 splits per type with the best model (legacy and gated), paired vs R3; Ta per run. Accept: Δ with SE; decision per §5 last bullet.
- U7: paper material: one figure (eff_match vs Δt at ratio 0.5, 1, 2 for MF, M0, best ML, with the ceiling), one table (Δt 16/20/24 ns cells), two sentences on σ_E/E invariance and the real-data fake fraction.

Compute on the Mac (18 cores, 48 GB, no GPU; estimate, unsourced): M0 ≈ 10 ms per waveform → 60 cells × 50 × 5 types in minutes, the dense grid in ≈ 1 h; M1 with ≤ 10⁵ parameters on 1.5·10⁶ waveforms: ≈ 1–3 h per training run of 20 epochs, inference ≈ 10⁴ waveforms/s, so a 118-split rescan of the pulse list is ≈ 1–2 h; M3 ≈ 3–5× M1.

## 7. Improvement estimate from other experiments

Survey (researcher-worker, 2026-10-27, web). No published result gives a minimum resolvable two-pulse separation at a stated efficiency for a network against a matched filter or template deconvolution; the comparisons that exist report resolution under pile-up or amplitude error against "no correction". What was found:

| Source | Detector, sampling, rise | Model vs baseline | Reported gain | Train data |
|---|---|---|---|---|
| ATLAS LAr, Chiedde et al., arXiv:2111.08590 (JINST proc.; the journal paper Comput. Softw. Big Sci. 5, 19 (2021) was not retrieved) | LAr EM barrel cell, 25 ns/sample, peak between samples 2–3 (rise ≈ 1–2 samples), pulses overlap over ≤ 25 bunch crossings | 3/4-conv CNN (tagging + energy heads), RNN, LSTM vs optimal filter (OF, 5 samples + peak finder) | tagging efficiency "well above 90 %" vs ≈ 80 % for OF for deposits > 3σ (240 MeV); all networks beat OF in resolution, most at small gaps (Fig. 4, figure only) | AREUS simulation, ⟨μ⟩ 140 |
| ATLAS LAr, Aad et al., arXiv:2510.11469 (EPJC 2026) | same cell; networks of 89–431 parameters | RNN, Dense, CNN vs OF | energy resolution ≈ 80 MeV vs ≈ 90 MeV for OF at ⟨μ⟩ 200 (≈ 10 % relative); OF underestimates the scale | AREUS simulation, 10⁶ events |
| PMT waveforms, Xu et al., arXiv:2112.06913 | 1029-sample windows, SER τ 8 ns | 5-layer 1D CNN vs fast stochastic matching pursuit (FSMP), LucyDDM, Fourier deconvolution | CNN on par with the model-based FSMP, better than LucyDDM and Fourier; FSMP best overall; CNN charge distribution distorted | simulation |
| Micromegas, Kim et al., arXiv:2304.14496 (NIM A 1055, 168492) | TexAT, 40 ns/bin, 100-bin windows | 1D-CNN autoencoder vs no correction | peak-amplitude error 24 → 5 % and 259 → 15 %; timing error 2.7 → 0.4 bins | synthetic sums of labelled real pulses |
| Belle II ECL, Wermuth, BELLE2-MTHESIS-2021-072 (docs.belle2.org/record/2557) | CsI(Tl) | NN vs offline multi-template fit | better crystal energy resolution for photons, robust to beam-background pile-up; snippet-level only | not stated |
| Stilbene, Fu et al., Ann. Nucl. Energy 120 (2018) 410 (via MTV 2020 poster) | organic scintillator | ANN vs charge integration | recovers both pulses only at ≥ 60 ns separation; sampling rate not given; snippet-level | not stated |

Translation to our setting (all numbers in this block are estimates, unsourced, derived from the rows above and §1):
- Units: our Δt 16 / 20 / 24 ns is 4 / 5 / 6 samples and 2.7 / 3.3 / 4 rise times (rise 6 ns). The ATLAS gain sits at gaps of 1–2 samples ≈ 1 rise time, where the OF's single-pulse assumption breaks; our hardest cells are easier than that in rise-time units, and the second pulse (≥ 0.3 × 500–2000 ADC ≫ σ_ped ≈ 1–2 ADC, estimate, unsourced) is far above the 3σ threshold at which ATLAS tags > 90 %. The MF's 0 % at Δt ≤ 16 ns is therefore an algorithm limit (12 ns floor, greedy subtract-and-repeat, candidate finder), not an information limit; the PMT result (CNN ≈ model-based fit when the pulse model is exact) says the ceiling will be set by the two-pulse χ² fit of M0, and a network will approach, not beat, it.
- Expected range if M0 confirms the ceiling: Δt 20 ns ratio 1 from 29–63 % to 80–95 %; Δt 16 ns ratio 1 from 0 % to 50–85 %; Δt 24 ns ratio 0.5 from 2–3 % to 70–95 %; Δt 12 ns ratio 1 from 0 % to 10–50 % (template error 1.2 % and the one-sample rise limit it); ratio ≤ 0.3 cells are set by the r_model 0.3 policy, not by the method, and lowering it trades against the real +23 ns, ratio ≈ 0.18 bump population (summary §3.8). Driving assumptions: the overlay model (template + real noise) matches data; fake rate held ≤ 0.03/event; eff_match tolerance ±4 ns unchanged.
- σ_E/E: on e-p and C, pile-up touches 0.15–0.41 % of E > 50 MeV channels; the time gate (removes out-of-time pulses) moves σ_E/E by ≤ 0.004 pp and a 34-point change in the Δt 20 ns efficiency moves it by ≤ 0.02 pp (§1), so the expected effect of any pile-up method is |Δ| ≤ 0.02 pp. The ATLAS 8 % resolution gain is for pulses overlapping in every event (⟨μ⟩ 200) and does not transfer to a 0.3 % rate. On Ta (10–11 % rate; gate gain −0.048 / −0.030 / −0.155 pp on runs 25800 / 26138 / 26230, 20261027 §7) a recovery method could reach a comparable gain per run, −0.03 to −0.15 pp (estimate, unsourced), with ≈ 16 calibrated modules per split.

## 8. Risks and compatibility

- Package installs: `torch`, `scikit-learn`, `xgboost` are not installed; numpy 2.4.6 and scipy 1.18.0 are. U0 and U1 need nothing new. U2 needs a CPU torch wheel in a separate environment (`python-env` skill procedure); not done without approval.
- The overlay model is "template + real noise", so the real-data checks of U5 are mandatory before any claim.
- No C++: the production WaveAnalyzer stays; the ML pulse list enters only through the Python benchmark, as the MF does.

## 9. Acceptance criteria and validation commands

- `output/ml_pileup_20261027/u0/ceiling_<slug>.md` with the Δt 12–24 ns cells and the go/no-go line.
- `python3 pileup_generator.py --self-test` and `benchmark_pileup_mf.py` on one existing cell reproduce `pileup_3way_full.md` numbers before and after the `--phase-jitter` edit (the default stays byte-identical).
- `best50_compare.py pileup --ref output/pileup_eff_20261024_fix --new <M0 or ML dir> --out output/ml_pileup_20261027/compare/<arm>` for every arm.

## Open Questions

1. May a CPU `torch` be installed in a separate conda/uv environment for U2 (nothing in the base environment changes)?
2. Is U0 alone acceptable as the paper's statement if the χ² ceiling kills the idea?
3. Should Ta per-run tables use all three runs even where a run has < 2000 pile-up channels (26138 has ≈ 5 k channels in the (a) denominator)?

## Todo

- [ ] U0 χ² ceiling (M0) on Δt 12–24 ns cells, 5 run types, go/no-go recorded
- [ ] U1 dense overlay sets (even/odd), `--phase-jitter` added with regression check
- [ ] U2 M1 trained, model card written (after Open Question 1)
- [ ] U3 frozen-test comparison tables and amplitude-bias table
- [ ] U4 M2/M3 (conditional)
- [x] U5 real-data RF-phase fake fraction and Ta comb, per run (two-pulse fit, 2026-10-28: `output/ml_pileup_20261027/U5/`)
- [ ] U5b `--pulse-list` hook in the benchmark, R3 reproduced on one split
- [ ] U6 σ_E/E on 4 splits per type, paired vs R3
- [ ] U7 paper figure and table
