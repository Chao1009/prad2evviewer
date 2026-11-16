# PbWO4 energy resolution: MIN_SEP_NS 12 → 8 ns

**Date:** 2026-11-08
**Branch:** `waveform_jz`
**Previous note:** `notes/PbWO4_energy_resolution_20261105.md`

Paths are relative to `analysis/pyscripts/` unless stated.

## 1. MIN_SEP_NS 12 → 8 ns (mf_pulse.py)

- **Change.** `mf_pulse.py:69`: `MIN_SEP_NS` 12.0 → 8.0. It is used at `:1013` (no new pulse within the separation of an accepted one), `:1035` (candidate acceptance) and `:1117` (merge of close peaks after the GN steps). No mf_pulse.py docstring states the value.
- **Why.** Synthetic two-pulse study (`waveform_analysis/pileup_analysis/pileup_generator_stage1.py`, `pileup_deconvolution_study.py`): 5 modules near the beam hole, 2p24GeV_H2 best-50 spline templates with shape_err smearing, 10k events per cell, no separation guard in the standalone MF port.
- Secondary pulses resolve with ~100 % efficiency from dt = 8 ns for A2/A1 ≥ 0.1 and from 12 ns for 0.05. At 4 ns (one sample) the pair is absorbed into one wider pulse for r ≤ 0.5; this is treated as a plateau, not pile-up, so the grid now starts at 8 ns.
- r = 0 controls: no false positives unsmeared; 0.6–3.4 % spurious second pulses at A1 = 3000 with K = 1 template smearing.
- The 12 ns guard was introduced 2026-10-23 against a −6 ns pre-pulse split on real data (+0.30 pp); 8 ns still blocks that split.
- **To do.** Verify on real data before relying on it: energy-resolution benchmark and merged-pair count at 8–12 ns on run 025320.
- **Results.** `output/pileup_stage1/pileup_stage1.npz`; `waveform_analysis/pileup_analysis/efficiency*.{tsv,pdf}`.
