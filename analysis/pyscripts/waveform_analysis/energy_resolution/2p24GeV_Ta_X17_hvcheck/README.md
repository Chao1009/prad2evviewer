# X17 Ta runs against the inner-layer HV change (08-10 to 08-14)

Question: is run 26230 (08-12) affected by the HV raise of the two innermost crystal layers [4540517], partly reversed on 08-14 because the resolution worsened [4543254], and is there a better Ta run? Source: stored v3 outputs (no rerun of the benchmark); `hvcheck.py` (one process, about 2 min for the tables plus about 15 min decoding the example waveforms) produces `summary.tsv`, `report.pdf` and `waveform_examples.tsv`. Standing rule PbWO4 only.

## Runs

| run | date | HV window | calibration file | events seen / physics | splits | ep events (gate 4 ns) |
|---|---|---|---|---|---|---|
| 25800 | about 07-11 to 07-13 (no logbook entry; between 25783 on 07-10 and 25833 on 07-14/15) | before | x17_25825 | 870618 / 862710 | 8 | 78270 |
| 26138 | 08-06 | before | x17_26094 | 415288 / 412340 | 4 | 27266 |
| 26230 | 08-12 | inside | x17_26214 | 1252541 / 1244962 | 13 | 72935 |

The EVIO file timestamps are the copy time (2026-11-14) and the physics events carry no unix time, so dates come from the elog summary (`notes/elog_pradlog_summary_20261117.md`). All three runs are after the 07-03 latency change [4519495]. No local Ta run is after the window (the run ended 08-14).

## Per-run numbers (gate 4 ns, the 1.987 % arm; out-of-sample Crystal Ball, mean ± standard error over the run's splits)

| run | σ/E oos, gate 4 ns (%) | σ/E oos, legacy (%) | fitted-g median ratio, layers 1-2 / layers ≥ 5 | calibration-file factor ratio | factor × g ratio | raw peak ADC per MeV, layer 1 (10 examples, median) | median hit time (ns) | hits in 140-156 ns |
|---|---|---|---|---|---|---|---|---|
| 25800 | 1.984 ± 0.032 | 2.145 ± 0.034 | 0.993 | 1.683 | 1.638 | 0.59 | 151.5 | 79.8 % |
| 26138 | 1.942 ± 0.078 | 2.331 ± 0.090 | 0.998 | 1.877 | 1.854 | 0.64 | 151.2 | 79.6 % |
| 26230 | 2.003 ± 0.040 | 2.389 ± 0.039 | 0.995 | 0.857 | 0.827 | 1.37 | 150.9 | 78.8 % |

Layers are PbWO4 rings around the beam hole (layer 1 = 12 modules, layers 1-2 = 32). The reference is layers ≥ 5 (211 to 238 modules with a fitted gain; 137 for 26138). Time statistics use matched-filter hits with E ≥ 200 MeV, 100 < t < 200 ns, t = fitted onset + template t50, raw (no per-module time offsets).

## Reading

- The HV raise is visible, but not in the fitted gains. The fitted g is a correction on top of each run's calibration file, so it stays at 1 (ratio 0.99 to 1.00) in all three runs. The step is in the calibration file: the inner/outer factor ratio drops from 1.68 and 1.88 (x17_25825, x17_26094) to 0.86 (x17_26214), a factor 2.2, matching "2-3x gain" [4540647]. The raw waveforms show the same: the layer-1 peak per reconstructed MeV rises from 0.6 to 1.4 ADC/MeV (about 2.2x).
- Pulse timing is the same in the three runs (median 151 ns, 79 to 80 % inside 140-156 ns, inner-layer median 150.5 to 151.5 ns). No latency-related shift.
- Resolution: 26230 is 2.003 ± 0.040, 25800 is 1.984 ± 0.032, 26138 is 1.942 ± 0.078. The differences (0.02 and 0.06 points) are inside one to two standard errors; no resolution degradation from the HV window is resolved. 26230 is the worst of the three at gate 4 ns and by a wider margin at legacy (2.39 against 2.15 for 25800), but 26138 has large errors and the legacy gap 25800 to 26230 is 0.24 ± 0.05.
- The two runs outside the window give 25800 + 26138 pooled (inverse-variance weights): 1.978 ± 0.029. The three-run mean in the summary is 1.987 ± 0.025.

## Recommendation

Two Ta runs lie before the window, so the preferred X17 number candidate is their pooled out-of-sample Crystal Ball σ/E = 1.98 ± 0.03 % (gate 4 ns; 25800 alone 1.984 ± 0.032, 26138 alone 1.942 ± 0.078). It agrees with the all-run 1.987 ± 0.025 % within a fraction of a standard error, so the 26230 inclusion changes the number by 0.01 points. Keeping 26230 costs nothing but carries the unresolved risk that 08-12 sat in the period reported as worse [4543254]. The split-to-split spread is large for 26138 (4 splits), so the pooled value is dominated by 25800. The elog note names no other Ta run beyond 26138 and 26230; X17 Ta production before the window runs 07-05 to 08-06 (runs about 25608 to 26138, X17_v1 to v4), and any of those would be fetchable, but the elog summary does not list their numbers, so the run list would need RCDB (the note's open question 1).

## Files

- `summary.tsv`: per-run table (all numbers above, plus per-layer time medians and the effective ratio).
- `report.pdf`: page 1 summary plots (σ/E per run, gain ratios, hit-time distributions with the 140-156 ns window shaded); page 2 fitted-gain maps per run; pages 3 to 5 ten in-time layer-1 waveforms per run (the 30 are in `waveform_examples.tsv`).
- `hvcheck.py`: reproduces everything; reads `output/energy_resolution/control_phase2_v3_20261114/v3_defaults/2p24GeV_Ta_X17/` and the scan files.
