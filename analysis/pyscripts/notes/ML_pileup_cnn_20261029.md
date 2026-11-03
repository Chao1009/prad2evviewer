# Pile-up CNN on single-module PbWO4 waveforms: structure, training, evaluation, M1 → M2 → M2b → M3

**Date:** 2026-10-29
**Sources:** `output/ml_pileup_20261027/` (M1, M2, M2b, M3, U0, U5). Paths in parentheses are relative to that directory unless they start with `notes/`, `plans/` or a script name in `analysis/pyscripts/`.
**Plan:** `plans/plan_ml_pileup_20261027.md` (user decisions 2026-10-27; U2 line). Ta runs 25800, 26138, 26230 are reported per run, never pooled.

All training is local: Mac CPU, CPU torch in `~/.venvs/prad-ml-cpu`, 4 torch threads per job (M1/tools/m1_train.py). Research models only; nothing is used in production (M3/model_card.md).

## 1. Purpose and the classical baseline

**Task.** Label a 100-sample FADC waveform (4 ns/sample) as one pulse or two, and regress the onset separation dt and the amplitude ratio A_later/A_earlier (M1/tools/m1_model.py docstring). User target: > 80 % detection efficiency at Δt = 8 ns, ratio 1 (two FADC samples); grid Δt ∈ {4, 8, 12, 16, 20, 24, 32} ns; e-p first, Ta per run later (plans/plan_ml_pileup_20261027.md, user decisions).

**Production MF.** 0 % efficiency at pileup_3way label Δt ≤ 16 ns (onset separation ≤ 9.8 ns) in every cell of the production pseudo-data, set by the 12 ns minimum separation and the greedy argmax search (plans/plan_ml_pileup_20261027.md §1; notes/PbWO4_energy_resolution_20261027.md §8).

**Two-pulse χ² fit (U0, `twopulse_fit.py`).** Best-50 spline template at continuous time, 1- and 2-pulse fits; N = 2 iff χ²₁ − χ²₂ > thr, thr = (1 − fake) quantile on clean singles (notes/PbWO4_energy_resolution_20261027.md §8).

- **Ceiling on synthetic doublets:** 100 / 100 / 100 % at Δt 8 ns, r 1 (0.73 / 2.24 / 3.49 GeV H2), evaluation fake 0.021 / 0.028 / 0.024. 100 % at Δt ≥ 12 ns for r ≥ 0.3. Δt 4 ns r 1: 86 / 83 / 75 % on `real`, limited by template-shape error (notes/PbWO4_energy_resolution_20261027.md §8).
- **Consequence stated there:** at Δt ≥ 8 ns there is no room left for ML on synthetic data; ML matters below 8 ns or for behaviour on real data (same §8).

**Real-data fake problem (U5).** The U0 thresholds do not transfer (U5/report.md; notes/PbWO4_energy_resolution_20261028.md §5):

- Detections with |dt| < 16 ns: 3.41 / 4.32 / 1.10 / 4.46 % of fired channels (0.73 / 2.24 / 3.49 H2, 2.24 C) against 0.015–0.045 % genuine expected from the MF density, i.e. 70–170× (U5/report.md).
- The excess sits at dt = +8 … +50 ns after pulse 1 (tail ripple from template mismatch); 6–7 % of 2.24 H2 channels fire at 20–100 ADC (U5/report.md).
- The cure costs efficiency: at the e-p operating points (thr 474–632) the 8 ns efficiency is 98.5–98.9 % at r 1 and 22.2–45.6 % at r 0.3 (U5/report.md, operating points).
- Same effect on the CNN test set: χ² at 2 % synthetic fake fires 16.7 % on real in-domain singles; at 2 % real fire its 8 ns r 0.3 efficiency drops to 66.8 % on injected-later doublets (M1/tables/headline.md, M1/tables/headline_basefirst.md).

The CNN is meant to keep the χ² fit's small-ratio efficiency at a fake rate measured on real singles.

## 2. Data

Run type for M1–M2b: 2.24 GeV H2 only (M1/tools/m1_data.py, `SLUG`). M3 uses the three Ta runs (M3/model_card.md).

**Singles and doublets (M1/tools/m1_data.py).**

- Bases: clean real pulses from the best-50 pools: MF single, mf_sec_ratio < 0.05, no clip (max < 4095), in-bank non-fallback module, doublet-audit rows removed. M1 cuts: 100–4000 ADC, fitted onset ≥ 150 ns.
- Every base enters twice: as a single (label 0) and as the base of a doublet (label 1).
- Doublet draw law: dt uniform 2–40 ns with p 0.75, else 40–100 ns; ratio log-uniform 0.03–3.3; injected pulse first with p 0.5 unless it would start before 120 ns; A_base + A_inj ≤ 4000 ADC (up to 20 redraws).
- Injected pulse: the module's best-50 template, noise-free; the base carries the real noise.
- Rows: train 578 172, val 49 630 (M1/model/config.json).

**Hold-out (M2/model_card.md; M1/tools/m1_data.py).**

- Modules: 676 train / 113 val / 338 test (M1 `modules.json`).
- Splits: 36 train, 18 val; test 025228_00100, 025228_00200; real 025228_00050, 025260_00400, 025299_00700, 025320_00300. Test modules and the 6 test/real splits never enter training.
- Pulses: test singles are split by parity into calibration (half 0) and evaluation (half 1); real channels into half 0 / half 1.

**Test sets (M1/tools/m1_eval.py, M1/tools/m1_data.py).**

- **U0 grid** (`test_synth.npz`): dt {4, 8, 12, 16, 20, 24, 32} ns × ratio {0.05, 0.1, 0.2, 0.3, 0.5, 1, 2}, 2 ns jitter, 500 doublets per split, cell and arm (1000 per cell and arm). Arm `real`: template injection; arm `real2`: real partner pulse from the same module.
- **Injected-later subset** ("base-first"): leading pulse = base, onset ≥ 150 ns. Injected-first doublets carry an onset shortcut for the CNN (M1/tables/headline.md header).
- **Real 2.24 H2** (`real_2p24.npz`, U5 channels with an MF pulse 1; 480 033 channels, M1/tables/summary.json):
  - S_dom: MF single, fitted 100–4000 ADC, fitted onset ≥ 150 ns, no clip (39 197);
  - S100: MF single, MF A1 100–4000 ADC (94 357);
  - U5_all: every channel (480 033).

**Changes per step.**

| step | base onset | shift | amplitude floor | real negatives | rows train / val | source |
|---|---|---|---|---|---|---|
| M1 | t_b ≥ 150 ns | none | 100 ADC | none | 578 k / 50 k | M2/model_card.md |
| M2 | t_b ≥ 135 ns | k uniform in [−5, 5] samples (±20 ns), same k for a single and its doublet | 20 ADC; 25 % of bases (A_b ≤ 400) scaled to 20–50 ADC with correlated noise added back (pool floor 50 ADC) | 79 k unfitted pool rows, MF single, 100 < mf_t0 < 200 ns | 554 k / 59 k | M2/model_card.md, M2/tools/m2_data.py |
| M2b | as M2 | as M2 | 20 ADC, real pulses: scaled-base rows dropped (train 118 730, val 12 596) | M2's 79.2 k + 80 k real MF singles (mf_A1 ≥ 20, no clip, not audit-flagged) from 5 new 2.24 H2 splits, each also used as a doublet base (78 867 doublets) | 594 247 / 64 644 | M2b/model_card.md, M2b/model/config_s20261028.json |
| M3 (per Ta run) | as M2 | as M2 | 20 ADC | the run's own 80 k / 9 k real singles + doublets on them; injected pulse = Ta n200 template | ≈ 152 k / 17 k per run | M3/model_card.md |

The M2b real singles: 45 % at 20–50 ADC, 19 % at 50–100, 17 % at 100–300, 11 % at 300–1000, 8 % ≥ 1000; onset 1 / 50 / 99 % = 142.3 / 149.5 / 157.1 ns (M2b/model_card.md). The M2b test sets are M1's, sha1 asserted (M2b/model_card.md).

**Ta test data (M3/model_card.md).** Non-U5 Ta splits, first 30 000 physics events: 25800 three splits, 26138 one split (00199, the only non-U5 one on disk), 26230 three splits. Per run: grid bases 11.2 k / 3.0 k / 11.3 k; held-out channels 410 k / 133 k / 481 k (M1 test modules). Ta singles onset 1 / 50 / 99 % = 104–106 / 146–147 / 186–188 ns, against e-p 142–157 ns.

## 3. Model

Code: M1/tools/m1_model.py (`M1CNN`, `ARCH`, `preprocess`). M2, M2b and M3 import it unchanged (M2/tools/m2_train.py, M3/tools/m3_train.py).

**Input normalisation (preprocess v1).** y = raw − ped (MF pedestal); s = max(max_t y, 1); two channels of 100 samples:

- ch0 = y / s (peak-normalised shape);
- ch1 = log10(s) / 3.5, constant over the 100 samples (amplitude scale).

**Layers.** Each conv is `Conv1d` (padding d(k − 1)/2, length preserved) → `BatchNorm1d` → `ReLU`. Parameter counts are computed from `ARCH`.

| # | layer | in → out ch | kernel | dilation | length out | params |
|---|---|---|---|---|---|---|
| 1 | Conv1d + BN + ReLU | 2 → 32 | 7 | 1 | 100 | 480 + 64 |
| 2 | Conv1d + BN + ReLU | 32 → 32 | 5 | 1 | 100 | 5 152 + 64 |
| 3 | Conv1d + BN + ReLU, MaxPool1d(2) | 32 → 48 | 5 | 2 | 50 | 7 728 + 96 |
| 4 | Conv1d + BN + ReLU, MaxPool1d(2) | 48 → 48 | 5 | 2 | 25 | 11 568 + 96 |
| 5 | Conv1d + BN + ReLU | 48 → 32 | 3 | 1 | 25 | 4 640 + 64 |
| 6 | Flatten, Linear + ReLU + Dropout(0.1) | 800 → 64 | | | | 51 264 |
| 7a | cls head, Linear | 64 → 1 | | | | 65 |
| 7b | reg head, Linear | 64 → 2 | | | | 130 |

**Parameter count:** 81 411 (M1/model/config.json). The plan had foreseen 4–6 conv layers, kernel 5–7, 16–64 channels, ≈ 30–100 k parameters (plans/plan_ml_pileup_20261027.md §3).

**Heads.** `cls`: one logit, score = sigmoid(logit); the working points cut on the logit. `reg`: two outputs, decoded as dt_pred = 10 · reg₀ ns and ratio_pred = exp(reg₁) (m1_model.py `predict`).

## 4. Loss and training

**Yes, there is a loss function.** `loss_fn` in M1/tools/m1_model.py, used unchanged by M1, M2, M2b and the M3 fine-tunes:

L = L_cls + w_reg · L_reg

- **L_cls:** binary cross-entropy on the logit, `torch.nn.functional.binary_cross_entropy_with_logits`, mean over all rows; target label ∈ {0, 1}.
- **L_reg:** Smooth-L1 (Huber), `torch.nn.functional.smooth_l1_loss` with beta = 0.1, mean over the doublets in the batch and the two targets [dt/10 ns, ln(A_later/A_earlier)]. With dt/10 the beta corresponds to 1 ns.
- **Masking:** L_reg uses only rows with label > 0.5 (doublets; dt and ratio are NaN for singles). A batch with no doublets gives L_reg = 0. The mask is the truth label, independent of the classifier output.
- **w_reg:** CLI `--w-reg`, default 1.0 (m1_train.py, m2_train.py). Value used: 1.0 in every M1, M2, M2b and M3 run (config*.json).
- Weight decay 1e-4 is applied by AdamW (decoupled), outside L.

**Optimiser and schedule (M1/tools/m1_train.py).**

| item | M1 | M2, M2b | M3 fine-tune (per run) | M3 pooled-training variant |
|---|---|---|---|---|
| init | random, `torch.manual_seed(seed)` | random | M2b s20261028 checkpoint | M2b s20261028 |
| optimiser | AdamW, lr 2e-3, wd 1e-4 | same | AdamW, lr 5e-4, wd 1e-4 | same as fine-tune |
| schedule | OneCycleLR, max_lr = lr, pct_start 0.1, stepped per batch | same | same | same |
| batch | 512 | 512 | 512 | 512 |
| epochs | 10 (script default 30) | 10 | 8 | 4 |
| seeds | 20261028 | 20261028 (primary), 20261038, 20261048 | 20261130 / 20261131 / 20261132 (25800 / 26138 / 26230) | 20261133 |
| best epoch | 10 | M2 9 / 10 / 10; M2b 10 / 10 / 10 | 8 / 8 / 8 | 3 |
| wall | 32.0 min | M2 30 / 48 / 42 min; M2b 32.0 / 32.1 / 32.3 min | 7.3 / 6.5 / 6.5 min | 9.8 min |
| source | M1/model/config.json | M2/model_card.md, M2b/model_card.md, configs | M3/tools/m3_train.py, M3/model_card.md, configs | same |

- **Early stopping:** the checkpoint kept is the epoch with the lowest validation total loss L. Training stops after `--patience` 5 epochs without improvement, or before an epoch when elapsed time + the longest epoch so far would exceed `--max-minutes` (55 used for M1, M2 and M2b, 50 being the m1_train.py default; ≤ 25 for M3) (m1_train.py, m2_train.py, m3_train.py). Every recorded run stopped at max_epochs (`stop_reason` in all configs).
- **No per-sample augmentation** during training; the M2 shift is applied once at data build (M2/model_card.md).
- **Validation AUC (best epoch):** about 0.99 on e-p; Ta fine-tunes 0.941 / 0.943 / 0.932 for 25800 / 26138 / 26230 (M3/model_card.md).
- **Plan U2** asked for 3 seeds and a model card with parameters, epochs and wall time (plans/plan_ml_pileup_20261027.md U2). M1 has 1 seed and no model card; M2 and M2b have 3 seeds and cards; M3 has 1 seed per variant.

## 5. Evaluation protocol

Fire = logit > thr for the CNN, dchi > thr for χ² (U0 rule, no veto) (M1/tools/m1_eval.py).

**Working points.**

- **Synthetic 2 % ("syn"):** thr = 2 % fire on synthetic calibration singles (test half 0); synthetic fake quoted on evaluation singles (half 1).
- **Real 2 % ("real"):** thr = 2 % fire on real S_dom half 0; real rates quoted on half 1. This is the headline WP from M2 on.
- **Ta WP-2 %:** 2 % fire on the run's own S_dom half 0, rates on half 1 (M3/tools/m3_eval.py).
- **Ta WP-c:** fire on S_dom half 0 = 2 % + the run's empirical contamination: 10.1 / 10.2 / 17.8 % for 25800 / 26138 / 26230 (M3/tables/contamination.tsv, `wpc_target_fire_S_dom`). For M2b it reproduces the e-p threshold (M3/model_card.md).

**Metrics.**

- **Efficiency per cell:** fire fraction of grid doublets, arm real [real2], all doublets and injected-later only (M1/tables/headline_basefirst.md).
- **ROC at 8 ns:** AUC and TPR at FPR 0.1 / 0.5 / 1 / 2 / 5 % against the 5171 synthetic evaluation singles, per ratio (M2b/tables/roc8.tsv).
- **Fake vs onset:** every 5 ns onset bin in 130–170 ns with n ≥ 1000 within ±1 pp of the pooled 130–170 ns rate (M3/tools/m3_eval.py; on e-p these bins span 140–160 ns, M2b/model_card.md).
- **Fake vs amplitude:** bins 20–100, 100–300, 300–1000, ≥ 1000 ADC each within ±2 pp of the pooled singles rate (M3/tools/m3_eval.py, M2/model_card.md).
- **Efficiency criterion:** dt 8 ns, r 0.3, arm real, injected-later ≥ 80 % (M2/model_card.md).
- **Shift test:** synthetic evaluation singles shifted by ±4 / 8 / 12 ns, fire at the syn WP (M1/tables/shift_test.tsv, M2b/tables/shift_test.tsv).
- **Residuals:** median and 16–84 % half-width (hw68) of dt_pred − dt_true on fired doublets (M1/tools/m1_eval.py; M2b/tables/residuals.tsv).
- **MF agreement:** fire fraction on real MF pile-up candidates per (|dt|, ratio) cell with n ≥ 30, and the median dt offset from MF (M2b/tables/mf_agreement.tsv).
- **Close-pair excess per channel:** fraction of channels with CNN fire and dt_pred < 16 ns (χ²: U5 detection with |dt| < 16) minus the MF-expected genuine rate (U5 rule: density of ratio > 0.3 MF extras in 24–60 / 16–60 ns × 32 ns) (M3/tools/m3_eval.py; U5/report.md).

## 6. Results M1 → M2 → M2b (2.24 GeV H2, real WP)

Efficiencies on injected-later doublets, arm real. M1 is one seed; M2 and M2b are 3-seed mean ± std where the source gives one. Fake rates on half 1.

| model | 8 ns r 1 | 8 ns r 0.3 | 8 ns r 0.1 | 4 ns r 1 | syn fake | real fake S_dom / S100 / U5_all | onset flat (±1 pp) | amplitude flat (±2 pp) | shift −4 / −8 / −12 ns |
|---|---|---|---|---|---|---|---|---|---|
| CNN M1 | 100.0 | 97.8 | 45.6 | 96.1 | 0.7 % | 1.8 / 30.5 / 68.0 % | fails: 95 % at 140–145 ns, 4 % at 150–155 ns | fails: 51 / 25 / 40 / 30 % | 70 / 94 / 99 % |
| CNN M2 | 100 | 97.5 ± 0.2 | 39.8 ± 0.6 | 94.7 ± 0.8 | 0.7 % | 1.9 / 2.0 ± 0.1 / 4.1 ± 0.2 % | FAIL, max dev 1.0 / 1.5 / 1.3 pp (155–160 ns) | FAIL, max dev 2.7–3.0 pp (20–100 ADC 4.9–5.1 %) | 2.0 / 2.2 / 2.2 % |
| CNN M2b | 100 | 97.5 ± 0.3 | 42.3 ± 0.8 | 94.2 ± 0.7 | 0.67 ± 0.10 % | 1.94 / 2.12 ± 0.06 / 1.55 ± 0.12 % | PASS 3/3, max dev 0.39 / 0.28 / 0.22 pp | PASS 3/3, max dev 1.18 / 1.01 / 1.00 pp | 1.85 / 1.88 / 1.83 % |
| χ² fit | 100.0 | 66.8 | 2.1 | 63.4 | 0.2 % | 1.8 / 1.7 / 1.2 % | 1.3 / 1.3 / 2.2 % at 145–150 / 150–155 / 155–160 ns | 1.2 / 2.0 / 1.8 / 1.1 % | n/a |

Sources by row. M1: M1/tables/headline_basefirst.md (8 ns), M1/tables/cells_basefirst.tsv (4 ns r 1, real WP), M1/tables/headline.md (fakes), M2/model_card.md (M1 onset and amplitude bins), M1/tables/shift_test.tsv (shift, syn WP; unshifted 1.7 %). M2: M2/model_card.md. M2b: M2b/model_card.md, M2b/tables/seed_spread.tsv. χ²: M1/tables/headline_basefirst.md, M1/tables/headline.md, M2/model_card.md (amplitude), M2b/model_card.md (onset). Amplitude bins are 20–100 / 100–300 / 300–1000 / ≥ 1000 ADC. The M1 flatness entries are values; the criterion was defined at M2.

**At the syn WP** (same doublets): 8 ns r 0.3 / r 0.1 is M1 98.5 / 56.2, M2b 98.7 ± 0.2 / 52.7 ± 0.9, χ² 98.1 / 19.8 % (M1/tables/headline_basefirst.md, M2b/model_card.md). At this WP χ² fires 16.7 % on real S_dom, M2b 4.4 % (M2b/tables/headline.md).

**ROC at 8 ns, M2b primary seed, arm real** (M2b/tables/roc8.tsv):

| ratio | AUC CNN / χ² | TPR at FPR 0.1 % CNN / χ² | TPR at FPR 2 % CNN / χ² |
|---|---|---|---|
| 0.05 | 0.812 / 0.651 | 0.104 / 0.005 | 0.259 / 0.044 |
| 0.1 | 0.938 / 0.863 | 0.441 / 0.017 | 0.674 / 0.174 |
| 0.3 | 0.9994 / 0.9973 | 0.963 / 0.232 | 0.991 / 0.989 |

**Residuals** (dt 8 ns, r 0.3, arm real): M2b median +0.38 ns, hw68 0.51 ns; χ² +0.09 / 0.48 ns. At dt 4 ns, r 1: M2b +0.67 / 0.50 ns, χ² +0.07 / 0.28 ns (M2b/model_card.md, M2/model_card.md).

**MF agreement** (real MF candidates, cells n ≥ 30): ratio ≥ 0.3 M2b 83–95 %, M2 95–100 %, χ² 96–100 %; ratio < 0.3 M2b 77–86 %. The M2b loss sits in the 20–100 ADC candidates (M2b/model_card.md).

**Lesson per step.**

- **M1: onset shortcut.** Injected-first training doublets start as early as 120 ns while bases start at 150 ns, so an early onset meant "doublet" (M1/tools/m1_data.py `GEN`, `BASE_T_MIN`; M2/tools/m2_data.py). M1 fired on 95 % of real singles at 140–145 ns and on 68.0 % of U5_all; shifting singles 4 ns earlier raised the fire from 1.7 % to 70 % (M2/model_card.md, M1/tables/headline.md, M1/tables/shift_test.tsv). Its 8 ns efficiency on injected-later doublets was already 97.8 % at r 0.3.
- **M2: amplitude coverage.** The ±20 ns shift and the 135 ns onset cut removed the shortcut (U5_all 4.1 %, shift test flat). The 20–50 ADC bases were scaled pool pulses, and real singles at 20–100 ADC fired at 4.9–5.1 % against 1.0–1.1 % at ≥ 1000 ADC (M2/model_card.md, M2b/model_card.md).
- **M2b: real negatives.** Replacing the scaled bases with real singles down to 20 ADC brought 20–100 ADC to 1.5–1.7 % and U5_all to 1.55 ± 0.12 %, with 8 ns r 0.3 unchanged at 97.5 %. The cost is lower firing on low-amplitude MF pile-up candidates (ratio ≥ 0.3: 83–95 % against 95–100 %) (M2b/model_card.md).
- **The 155–160 ns bin** that set M2's onset failure is real pile-up: all-channel MF extras within 60 ns are 4–6× more frequent there than at 145–155 ns, and χ² fires 1.7× more often on the same singles. Under M2b the bin is +0.3 pp above the pool (M2b/model_card.md, M2b/tables/onset_edge_155.tsv).

## 7. Ta per run (M3)

**Contamination finding.** Ta MF singles carry genuine pile-up the MF misses, so a 2 % fire on Ta S_dom lands inside the pile-up population (M3/model_card.md). Per run (M3/tables/contamination.tsv, M3/model_card.md):

| run | MF multi-pulse fraction of kept channels | MF-density estimate (r > 0.3, dt < 24 ns), %/ch (lower limit) | empirical, M2b e-p WP − 2 %: S_dom CNN / χ² | same, singles A ≥ 20: CNN / χ² | WP-c target fire on S_dom |
|---|---|---|---|---|---|
| 25800 | 12.3 % | 0.98–1.01 | 8.1 / 10.4 % | 6.3 / 5.1 % | 10.1 % |
| 26138 | 12.3 % | 1.01–1.04 | 8.2 / 9.7 % | 6.2 / 4.5 % | 10.2 % |
| 26230 | 20.3 % | 1.65–1.71 | 15.8 / 17.4 % | 9.4 / 7.4 % | 17.8 % |

- The excess follows the run ordering of the channel pile-up rates, and CNN and Ta-bank χ² show the same U-shaped onset profile (M3/model_card.md).
- S_dom (onset ≥ 150 ns) is the late tail on Ta: 31 % (25800) and 42 % (26230) of S_dom sits at ≥ 155 ns (M3/model_card.md).
- The fine-tune learns part of the contamination as "single", since its negatives are about 6.3 / 6.2 / 9.4 % pile-up for 25800 / 26138 / 26230 (CNN, singles A ≥ 20; M3/tables/contamination.tsv, M3/model_card.md).

**Per-run results.** Held-out M1 test modules, half 1; eff on injected-later doublets, arm real, n200 injection; close = per-channel fraction with fire and dt < 16 ns, MF-expected genuine rate in brackets; excess = close − MF-expected high edge (M3/tables/headline_per_run.md).

| run | WP | method | thr | fire singles A ≥ 20 | eff 8 ns r 0.3 | close %/ch (MF-exp) | excess | MF agree 16–60 ns r ≥ 0.3 | onset / amp / 8 ns |
|---|---|---|---|---|---|---|---|---|---|
| 25800 | 2 % | M2b as-is | 9.78 | 1.07 % | 2.2 % | 0.12 (1.03–1.05) | −0.93 | 77.2 % | P / P / F |
| 25800 | 2 % | fine-tuned | 2.71 | 0.72 % | 84.2 % | 0.32 | −0.73 | 76.4 % | F / P / P |
| 25800 | 2 % | χ² | 2937 | 0.54 % | 1.8 % | 0.08 | −0.97 | 73.8 % | P / P / F |
| 25800 | c | M2b as-is | 0.63 | 8.37 % | 95.4 % | 2.64 | +1.59 | 97.9 % | F / F / P |
| 25800 | c | fine-tuned | 0.55 | 4.72 % | 94.9 % | 1.49 | +0.44 | 87.9 % | F / F / P |
| 25800 | c | χ² | 444 | 5.40 % | 62.7 % | 0.98 | −0.07 | 97.9 % | F / P / F |
| 26138 | 2 % | M2b as-is | 9.09 | 1.23 % | 5.7 % | 0.16 (1.01–1.05) | −0.90 | 78.3 % | P / P / F |
| 26138 | 2 % | fine-tuned | 2.83 | 0.58 % | 79.1 % | 0.25 | −0.81 | 71.9 % | F / P / F |
| 26138 | 2 % | χ² | 2491 | 0.62 % | 2.7 % | 0.10 | −0.95 | 78.3 % | F / P / F |
| 26138 | c | M2b as-is | 0.62 | 7.98 % | 93.7 % | 2.52 | +1.46 | 97.8 % | F / F / P |
| 26138 | c | fine-tuned | 0.64 | 4.26 % | 92.2 % | 1.28 | +0.22 | 85.8 % | F / F / P |
| 26138 | c | χ² | 395 | 5.37 % | 61.3 % | 0.93 | −0.13 | 98.4 % | F / P / F |
| 26230 | 2 % | M2b as-is | 14.43 | 0.65 % | 0.6 % | 0.01 (1.70–1.76) | −1.74 | 54.8 % | F / P / F |
| 26230 | 2 % | fine-tuned | 2.79 | 0.63 % | 66.2 % | 0.30 | −1.44 | 67.3 % | P / F / F |
| 26230 | 2 % | χ² | 5749 | 0.33 % | 1.1 % | 0.03 | −1.71 | 52.3 % | P / P / F |
| 26230 | c | M2b as-is | 0.65 | 11.22 % | 94.1 % | 3.58 | +1.84 | 98.0 % | F / F / P |
| 26230 | c | fine-tuned | −0.06 | 8.76 % | 92.4 % | 2.43 | +0.69 | 86.3 % | F / F / P |
| 26230 | c | χ² | 395 | 8.15 % | 58.2 % | 1.45 | −0.29 | 96.7 % | F / F / F |

**Reference, χ² at U5's Ta OP** (dchi 173 / 150 / 130), same channels: close 2.61 / 2.63 / 3.79 %, excess +1.56 / +1.57 / +2.05 % for 25800 / 26138 / 26230. U5's OP fake residual (symmetric, U5 splits): 0.69 / 0.77 / 1.14 % (M3/tables/headline_per_run.md, U5/report.md).

**Close-pair excess against target.** No numeric target is recorded (M3/model_card.md). Fine-tune excess at WP-c, +0.44 / +0.22 / +0.69 %/ch, is below U5's OP fake residual (0.69 / 0.77 / 1.14 %) and below χ² at U5's OP (+1.56 / +1.57 / +2.05 %) on every run. The excess assumes full efficiency on genuine close pairs; it is a rate-based fake estimate (M3/model_card.md).

**Verdicts (M3/model_card.md, M3/tables/headline_per_run.md).**

- **WP-2 %, 8 ns r 0.3 ≥ 80 %:** PASS only for the 25800 fine-tune (84.2 %). 26138 fine-tune 79.1 % FAIL; 26230 fine-tune 66.2 % FAIL. M2b as-is and χ² collapse to 0.6–5.7 %, because the 2 % quantile on contaminated S_dom lies in the pile-up tail (logit 9–14, dchi 2500–5700).
- **WP-c:** 8 ns passes for M2b as-is and the fine-tune on all runs (92.2–95.4 %); χ² gets 58.2–62.7 %.
- **Flatness on Ta:** fails for every CNN at WP-c, and for χ² in onset on every run; the U-shape and the amplitude rise follow the contamination's structure. On Ta singles the flatness criteria are not a fake-rate test.
- **MF agreement cost:** the fine-tune cuts the M2b close-pair excess by 2.7–6.6× at the same S_dom fire, at 10–12 pp lower MF agreement (85.8–87.9 % against 97.8–98.0 %).
- **Pooled-training variant** (three runs' train sets concatenated, evaluated per run): close pairs at WP-c 1.49 / 1.44 / 2.39 %/ch for 25800 / 26138 / 26230, against fine-tune 1.49 / 1.28 / 2.43; at WP-2 % 8 ns r 0.3 83.4 / 79.7 / 63.3 % (M3/tables/headline_per_run.md).

## 8. Current conclusion

**Established (e-p, 2.24 GeV H2).**

- At 2 % fire on real in-domain singles, M2b keeps 97.5 ± 0.3 % at 8 ns r 0.3 and 42.3 ± 0.8 % at r 0.1, against χ² 66.8 % and 2.1 % on the same doublets (M2b/model_card.md, M1/tables/headline_basefirst.md).
- Its fake rate is flat in onset (≤ 0.39 pp) and amplitude (≤ 1.18 pp) on real singles over 3 seeds; the test-module-only verdicts also pass (M2b/model_card.md).
- The two failure modes found on the way were data-construction effects: M1's onset prior and M2's scaled low-amplitude bases (§6).
- Loss (BCE + Smooth-L1, w_reg 1.0), architecture and trainer settings are unchanged from M1 to M2b (§4).

**Not established.**

- A close-pair operating point on e-p. At the real WP, M2b's fired channels with predicted |dt| < 16 ns are 0.42–0.43 % per channel (U5_all half 1), 10× the MF-expected genuine 0.042–0.045 % (M2b/model_card.md, U5/report.md). The acceptance proposed for this (< 2× genuine, < 0.09 %/ch) is not met at the 2 % WP and has not been scanned (M2b/model_card.md).
- Seed spread on Ta: one seed per variant (M3/model_card.md).
- Label noise: about 18 % of doublets on bases below 50 ADC carry an injected pulse below 3 ADC (M2b/model_card.md).

**What the e-p model is good for.** The e-p channel pile-up rate is 0.17–0.56 % per channel (MF, ≥ 2 pulses anywhere in the window; notes/PbWO4_energy_resolution_20261027.md §10), and genuine close pairs (|dt| < 16 ns) on 2.24 H2 are 0.042–0.045 % per channel (U5/report.md). At the 2 % WP at least ≈ 90 % of the CNN's e-p close-pair tags are fake (derived here: 0.43 against 0.045 %/ch). Any pile-up method's effect on σ_E/E is bounded on e-p and C: |Δ| ≤ 0.02 pp (notes/PbWO4_energy_resolution_20261027.md §8, referring to §7). On e-p, M2b serves as the validated recipe and the initial weights for Ta (M3 starts from M2b s20261028).

**What Ta needs next** (from M3/model_card.md known limitations and verdicts):

- Cleaner negatives: the Ta singles carry genuine pile-up (per run in the §7 contamination table), and the fine-tune learns some of it as "single".
- A per-run OP; WP-c is anchored on M2b's e-p threshold and assumes the e-p fake rate transfers.
- More than one seed per run, and more 26138 data (one split; 3.0 k grid bases, 1.5 k S_dom half-0 rows set its WP).
- A Ta fake-rate test other than flatness, which is not usable on contaminated singles.

**LCRC hand-off.** `plans/plan_ml_lcrc_20261028.md` exists and is on hold since 2026-10-29 (user decision, plan status line).
