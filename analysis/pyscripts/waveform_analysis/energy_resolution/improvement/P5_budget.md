# P5: resolution budget, separation table and phase gate (2026-11-13)

Inputs: `P1_reconcile`, `P2_regions`/`P2_results`, `P3_secondaries`, `P4_switches`, `P12_threshold`, `P15_tgate_calibrated` (this
directory). Tables: `P5_latest.tsv`, `P5_budget.tsv`, `P5_separation.tsv`, written by `tools/p5_tables.py` from the variant trees.
New runs: `output/energy_resolution/variants/P5/` (scripts `scan_lane.sh`, `run_combo.sh`, `run_ext.sh`; `logs/walltimes.txt`).
Rescans: `output/energy_resolution/scan_variants/{rmodel0p3,rebaseline}/`. Units: `res_cb` = Crystal Ball core σ/E of the ep pool
(%), mean ± SE over splits; Δ = variant − control, paired per (run, split), in pp; "oos" = even/odd cross-fit (D1).

**Tail columns (checked, BM:1326-1356).** `tail` (in `summary.tsv`, 11-49 %) and `tail_lo3` + `tail_hi3` = it are **shares of all
clusters in the fit window E' ± 36 % (± 800 MeV · E/2.24 GeV) that lie beyond μ ∓ 3σ_CB**, not a Crystal Ball tail parameter. The
window holds the Møller single clusters (0.73, 2.24 GeV), the radiative tail and, on Ta, quasi-elastic events, so the low side
(`tail_lo3`, 11-45 %) is mostly physics. Below they are labelled "window share below μ−3σ" / "above μ+3σ"; `tail_lo_f`, `tail_hi_f`
are the same shares beyond the fixed offset Δ_π = 145 MeV·E'/E. Only `tail_hi3` (0.05-0.16 % in the pool) is a detector tail.

## 1. Latest energy resolutions (pool, `res_cb` %)

Best gate = the P15 calibrated-time gate with the lowest oos `res_cb` (P15 final scan, 1-24 ns H2, 1-16 ns Ta; P15's one-SE rule
picks 12 / 16 / 12 / 12 ns H2 and 4 ns Ta, within one SE of these; H2 calibrated best = raw best within one SE). Offsets are applied
before the gate (user, 11-13).

| run type | legacy (in / oos) | best calibrated gate | in / oos | Δ in (oos) vs legacy | + log weight R7ri, same gate (measured) | + R_MODEL 0.3, same gate (measured) |
|---|---|---|---|---|---|---|
| 0.73 GeV H2 | 3.389 ± 0.014 / 3.420 | 20 ns | 3.365 ± 0.013 / 3.397 | −0.024 ± 0.007 (−0.023 ± 0.007) | null at 0.73 (P4b), not applied | not rescanned |
| 2.24 GeV H2 | 1.876 ± 0.005 / 1.889 | 16 ns | 1.856 ± 0.005 / 1.867 | −0.020 ± 0.002 (−0.021 ± 0.002) | 1.847 ± 0.005 / 1.858 | 1.833 ± 0.004 / 1.842 (bound only, §4 iii) |
| 3.49 GeV H2 | 1.585 ± 0.006 / 1.614 | 16 ns | 1.565 ± 0.006 / 1.595 | −0.020 ± 0.003 (−0.019 ± 0.002) | 1.538 ± 0.008 / 1.561 | 1.524 ± 0.006 / 1.552 (bound only) |
| 3.49 GeV H2 old runs | 1.559 ± 0.006 / 1.579 | 16 ns | 1.546 ± 0.006 / 1.564 | −0.013 ± 0.002 (−0.015 ± 0.002) | 1.509 ± 0.006 / 1.532 | not rescanned |
| 2.24 GeV Ta X17 | 2.308 ± 0.035 / 2.332 | 8 ns | 2.004 ± 0.024 / 2.019 | −0.304 ± 0.030 (−0.313 ± 0.028) | not run | Ta guard below |

- The H2 curve is flat from 12 to 24 ns within about 0.01 pp (P15); with log weight, 20 ns is again no better than 16 ns (+0.002,
  +0.002, +0.001 pp in; `variants/P5/toff_logw`). Log weight adds on top of the calibrated gate in full: Δ −0.010 ± 0.002 /
  −0.027 ± 0.003 / −0.037 ± 0.003 pp (2.24 / 3.49 / 3.49 old, g16 paired against P15 g16), the same as on g12 (P4b).
- **Projection (every real lever).** Measured and adoptable now (calibrated gate + log weight): 3.365 / 1.847 / 1.538 / 1.509 / 2.004 % (0.73 / 2.24 / 3.49 /
  3.49 old / Ta). Adding the two phase-2 levers that have a measured bound, assuming each removes its quadrature term in full and
  independently of the rest: P8 at the R_MODEL 0.3 gain on top of gate + log weight (−0.013 pp, 2.24 and 3.49 GeV only) and P9
  bringing every edge cluster to the centre resolution (q_edge 0.94 / 0.35 / 0.38 / 0.39 / 0.44 %, legacy, transferred unchanged
  to the best gate) gives **3.23 ± 0.03 / 1.80 ± 0.01 / 1.48 ± 0.01 / 1.46 ± 0.01 / 1.95 ± 0.06 %**. Not included (no size on the current basis):
  P6 (old basis −0.02 to −0.04 pp H2, −0.07 to −0.09 Ta), P16 (threshold), P17 (asymmetric window); P7 and P11 change coverage,
  not the pool number. Both bounds are optimistic (P8's setting must also pass the Ta guard; P9 cannot remove the crack loss).

## 2. Resolution budget (pool, legacy in-sample; quadrature sizes q in %, q = √(R² − R'²))

Reference curve 2.6 %/√E (PRad-II proposal, stochastic only): 3.05 / 1.74 / 1.39 % → excess over it 1.48 / 0.71 / 0.76 % at
0.73 / 2.24 / 3.49 GeV (0.70 % old runs). "Applied" terms are already in the benchmark number (size = cost of switching them off);
"lever" terms are available improvements (Δ, then q). Errors are paired SE; q errors propagate them. Full rows, incl. 3.49 old and Ta: `P5_budget.tsv`.

| term | 0.73 GeV | 2.24 GeV | 3.49 GeV | status | source |
|---|---|---|---|---|---|
| gains, applied (none → gains) | q 2.31 | q 0.56 (Δ −0.080) | q 1.33 | applied | `poscorr_comparison.tsv`, `res_cb_none` |
| position correction, applied | — (off at 0.73) | q 0.43 (Δ −0.048) | q 0.70 | applied | same |
| onset-phase correction, applied | — | q 0.47 ± 0.03 (Δ +0.058 ± 0.007 off; run 25320) | — | applied; production lacks it (P10) | P4a |
| time gate, calibrated, best | Δ −0.024 ± 0.007, q 0.40 ± 0.06 | Δ −0.020 ± 0.002, q 0.27 ± 0.02 | Δ −0.020 ± 0.003, q 0.25 ± 0.02 | lever, real | P15 |
| log weight R7ri (on g12) | Δ +0.006 ± 0.003, null | Δ −0.011 ± 0.002, q 0.21 ± 0.02 | Δ −0.026 ± 0.003, q 0.28 ± 0.01 (old −0.033) | lever, real at 2.24/3.49 | P4b |
| secondary pulses: R_MODEL 0.3 | not rescanned | Δ −0.025 ± 0.003 (q 0.30); g16 + log weight −0.013 ± 0.002; `tail_hi3` +0.010 ± 0.002 | Δ −0.030 ± 0.002 (q 0.30); −0.013 ± 0.003; `tail_hi3` +0.019 ± 0.004 | real, replicated; not adoptable (`tail_hi3` up, Ta g8 +0.083 ± 0.025) | P3B + P5 extension |
| secondary pulses: MIN_SEP 12 ns | — | Δ −0.000 ± 0.002, null | — | null | P3B |
| edge excess (pool → centre, bound) | q 0.94 (edge − centre +0.21 ± 0.04) | q 0.35 (+0.09 ± 0.01) | q 0.38 (+0.12 ± 0.01) | upper bound for P9 | P2 regions, paired |
| below-threshold loss fluctuation | q 1.31 (0.91-2.15) | q 0.54 (0.37-0.90) | q 0.42 (0.29-0.69) | model-dependent extrapolation (variance ∝ t², range t³-t) | P12 §3 |
| Møller rejection | +0.167 ± 0.013 (worse) | +0.009 ± 0.001 (worse) | 0 | no lever (U3: keep off) | P1a |
| non-linearity (clusterer factor) | −0.033 ± 0.011 (fails g12, replication) | identity | +0.095 ± 0.003 (worse; old +0.103) | no lever (U3: keep off) | P1g |
| gain solve vs old chain | — | 0.088 pp residual (`ref` 1.882 vs note 1.794), q 0.57 if real | — | unverified (old outputs gone); RA §2 refinement residual, old basis q 0.25 / 0.32 % → P6 | P1, RA §2 |
| noise / template shape / crosstalk | ≤ 0.05 pp / — / — | 0.06-0.18 / ≈ 0.2 / 0.10 (q) | 0.04-0.12 / ≈ 0.15 / 0.08 (q) | unmeasured, [I] of RA §2 | RA §2 |
| LMS drift, clipping, saturation | — | ≤ 0.07 % common; 0.014 % of clusters | | measured null (RA §2, P2 drift ≤ 0.05 pp) | RA §2, P2 |

Reading (quadrature, assuming independent terms; overlaps noted):
- 2.24 GeV: excess² 0.50 %². Measured levers combined (gate + log weight + R_MODEL, 1.876 → 1.833): 0.16 %²; edge bound
  0.12 %²; threshold term 0.29 %² (t² model). Together 0.57 %², so the remainder is 0 to 0.08 %² over the threshold-model range
  (0.14-0.81 %²).
- 3.49 GeV: excess² 0.57 %². Measured levers combined (1.585 → 1.524): 0.19 %²; edge bound 0.14 %²; threshold 0.17 %²
  (0.08-0.47). Remainder 0.07 %² (q ≈ 0.27 %) for the t² model, between 0.16 %² and none over the model range.
- 0.73 GeV: excess² 2.19 %²; gate 0.16, edge bound 0.88 and threshold 1.72 (t²) overlap (P12: at 0.73 GeV the loss acts at the
  centre and closes the edge gap), so they are not added; the threshold term alone can carry most of the excess.
- Overlaps: the time gate and R_MODEL both act on secondary pulses; the edge bound contains part of the position term. The
  measured combinations in §1 replace the sums where they exist.

## 3. Separation table (U7; pool and regions, legacy unless stated)

σ_E in MeV; `sep_pi` = Δ_π/σ_E; window shares in % (see the tail note); `eff_cut` = % above μ − c, `leak_cut` = % above μ + Δ_π − c
(c = 49 / 84 / 110 MeV). Mean ± SE over splits; regions from P2 (raw times), pool best gate from P15 / P5. Full table (also g12,
I1, I2, O, Ta, 3.49 old): `P5_separation.tsv`.

| run type | rows | n (k) | `res_cb` | σ_E | `sep_pi` | share < μ−3σ | share > μ+3σ | share < μ−Δ_π | `eff_cut` | `leak_cut` |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.73 | pool legacy | 329 | 3.389 ± 0.014 | 24.6 | 5.86 ± 0.02 | 44.7 ± 0.6 | 0.047 ± 0.006 | 21.1 ± 0.2 | 49.3 ± 0.6 | 0.01 ± 0.00 |
| 0.73 | pool, cal. 20 ns | 333 | 3.365 ± 0.013 | 24.5 | 5.90 ± 0.02 | 44.7 ± 0.6 | 0.049 ± 0.005 | 21.1 ± 0.3 | 49.4 ± 0.6 | 0.01 ± 0.00 |
| 0.73 | M | 24 | 4.149 ± 0.045 | 30.1 | 4.79 ± 0.05 | 7.1 ± 0.2 | 1.137 ± 0.076 | 4.4 ± 0.1 | 86.0 ± 0.2 | 1.04 ± 0.05 |
| 0.73 | T | 3 | 4.821 ± 0.181 | 35.0 | 4.21 ± 0.16 | 13.0 ± 0.9 | 0.736 ± 0.156 | 8.3 ± 0.4 | 75.2 ± 1.0 | 1.08 ± 0.14 |
| 0.73 | H | 1384 | 7.131 ± 0.042 | 51.8 | 2.79 ± 0.02 | 15.3 ± 0.2 | 0.250 ± 0.020 | 16.5 ± 0.1 | 60.9 ± 0.1 | 2.70 ± 0.02 |
| 2.24 | pool legacy | 1075 | 1.876 ± 0.005 | 41.9 | 3.45 ± 0.01 | 28.0 ± 0.1 | 0.160 ± 0.004 | 27.2 ± 0.1 | 68.6 ± 0.1 | 4.66 ± 0.03 |
| 2.24 | pool, cal. 16 ns | 1074 | 1.856 ± 0.005 | 41.5 | 3.48 ± 0.01 | 27.7 ± 0.1 | 0.126 ± 0.004 | 26.9 ± 0.0 | 69.1 ± 0.1 | 4.57 ± 0.03 |
| 2.24 | + log weight | 1074 | 1.847 ± 0.005 | 41.3 | 3.50 ± 0.01 | 27.7 ± 0.1 | 0.125 ± 0.004 | 26.9 ± 0.0 | 69.2 ± 0.1 | 4.48 ± 0.03 |
| 2.24 | + R_MODEL 0.3 | 1074 | 1.833 ± 0.004 | 41.0 | 3.53 ± 0.01 | 27.7 ± 0.1 | 0.131 ± 0.004 | 26.9 ± 0.0 | 69.2 ± 0.1 | 4.34 ± 0.03 |
| 2.24 | centre | 690 | 1.844 ± 0.006 | 41.2 | 3.51 ± 0.01 | 26.8 ± 0.1 | 0.197 ± 0.006 | 26.1 ± 0.1 | 70.2 ± 0.1 | 4.54 ± 0.04 |
| 2.24 | edge | 386 | 1.936 ± 0.009 | 43.3 | 3.34 ± 0.02 | 30.1 ± 0.1 | 0.098 ± 0.007 | 29.4 ± 0.1 | 65.6 ± 0.1 | 4.85 ± 0.05 |
| 2.24 | M | 67 | 2.043 ± 0.016 | 45.6 | 3.17 ± 0.03 | 11.2 ± 0.1 | 1.232 ± 0.073 | 10.8 ± 0.1 | 83.2 ± 0.2 | 8.58 ± 0.10 |
| 2.24 | O | 15 | 2.043 ± 0.033 | 45.6 | 3.20 ± 0.05 | 15.0 ± 0.5 | 1.458 ± 0.094 | 14.3 ± 0.4 | 79.5 ± 0.4 | 8.60 ± 0.27 |
| 2.24 | T | 10 | 2.326 ± 0.083 | 51.8 | 2.95 ± 0.12 | 33.4 ± 1.1 | 3.070 ± 0.239 | 33.3 ± 0.6 | 57.3 ± 0.7 | 10.27 ± 0.52 |
| 2.24 | H | 2938 | 2.365 ± 0.005 | 52.9 | 2.73 ± 0.01 | 58.4 ± 0.0 | 0.081 ± 0.002 | 59.2 ± 0.0 | 35.6 ± 0.0 | 4.44 ± 0.02 |
| 3.49 | pool legacy | 352 | 1.585 ± 0.006 | 55.2 | 2.62 ± 0.01 | 11.2 ± 0.1 | 0.137 ± 0.011 | 12.4 ± 0.1 | 84.3 ± 0.1 | 22.28 ± 0.09 |
| 3.49 | pool, cal. 16 ns | 352 | 1.565 ± 0.006 | 54.5 | 2.65 ± 0.01 | 11.0 ± 0.1 | 0.116 ± 0.009 | 12.1 ± 0.1 | 84.7 ± 0.1 | 22.15 ± 0.10 |
| 3.49 | + log weight | 352 | 1.538 ± 0.008 | 53.5 | 2.70 ± 0.01 | 11.1 ± 0.1 | 0.108 ± 0.011 | 12.1 ± 0.1 | 84.8 ± 0.1 | 21.77 ± 0.12 |
| 3.49 | + R_MODEL 0.3 | 352 | 1.524 ± 0.006 | 53.0 | 2.72 ± 0.01 | 11.1 ± 0.1 | 0.109 ± 0.011 | 12.0 ± 0.1 | 84.9 ± 0.1 | 21.59 ± 0.10 |
| 3.49 | centre | 227 | 1.540 ± 0.007 | 53.6 | 2.70 ± 0.01 | 9.9 ± 0.1 | 0.164 ± 0.014 | 10.7 ± 0.1 | 86.3 ± 0.1 | 22.21 ± 0.12 |
| 3.49 | edge | 125 | 1.662 ± 0.012 | 57.8 | 2.50 ± 0.02 | 13.5 ± 0.1 | 0.094 ± 0.015 | 15.6 ± 0.1 | 80.4 ± 0.2 | 22.16 ± 0.18 |
| 3.49 | M | 27 | 1.983 ± 0.035 | 68.9 | 2.10 ± 0.04 | 14.8 ± 0.3 | 0.597 ± 0.073 | 19.4 ± 0.2 | 76.4 ± 0.3 | 25.01 ± 0.31 |
| 3.49 | O | 6 | 1.951 ± 0.063 | 67.6 | 2.16 ± 0.07 | 20.7 ± 0.7 | 0.800 ± 0.123 | 25.2 ± 0.5 | 70.8 ± 0.6 | 23.22 ± 0.83 |
| 3.49 | T | 4 | 2.378 ± 0.111 | 82.3 | 1.80 ± 0.07 | 32.1 ± 1.8 | 1.553 ± 0.361 | 40.4 ± 1.3 | 55.8 ± 1.3 | 21.60 ± 1.00 |
| 3.49 | H | 937 | 2.411 ± 0.130 | 84.0 | 2.09 ± 0.41 | 49.0 ± 1.1 | 2.955 ± 2.843 | 55.4 ± 2.1 | 41.1 ± 2.0 | 16.75 ± 2.15 |
| Ta 2.24 | pool legacy | 56 | 2.308 ± 0.035 | 51.6 | 2.81 ± 0.04 | 44.6 ± 0.9 | 4.484 ± 0.224 | 44.9 ± 0.9 | 51.6 ± 0.8 | 11.29 ± 0.35 |
| Ta 2.24 | pool, cal. 8 ns | 55 | 2.004 ± 0.024 | 44.8 | 3.24 ± 0.04 | 45.8 ± 0.9 | 0.609 ± 0.038 | 45.4 ± 0.9 | 51.0 ± 0.9 | 5.01 ± 0.15 |

- The pool gains 0.08 in `sep_pi` at 3.49 GeV from gate + log weight (2.62 → 2.70) and loses 0.5 pp of `leak_cut`; the high-side
  window share falls 0.137 → 0.108 %. Region spreads are much larger: M, O, T and the hole ring H have `sep_pi` 1.8-2.2 at 3.49 GeV
  (pool 2.6), `tail_hi3` 4-11 × the pool (H 22 ×), and all of them run on prior gains (prior fraction 0.98-1.0; H 0.47-0.81).
- H (excluded from the pool) holds more ep clusters than the pool (2.9 M vs 1.1 M at 2.24 GeV), at `res_cb` 2.37 / 2.41 %.
- Edges: edge − centre `eff_cut` −4.4 to −5.9 pp and `tail_lo_f` +3.2 to +4.8 pp (≥ 25 SE at 2.24 and 3.49 GeV; 0.73: −4.0 / +2.0 pp, 9-20 SE); `tail_hi3` is lower at
  the edge. Edge clusters lose energy, they are not smeared up: an efficiency term at the fixed cut, not a leak term.

## 4. Rulings

(i) **P4b log weight R7ri: adopt at 2.24 and 3.49 GeV (keep 5.2 at 0.73).** The failed guard is an artefact of the σ-scaled
boundary: σ_E narrows by 0.9-1.3 MeV at 3.49 GeV, which moves μ − 3σ inward by 2.7-3.8 MeV and moves radiative-tail clusters
across it. At the fixed offset the low side improves: Δ`tail_lo_f` −0.084 ± 0.014 / −0.059 ± 0.012 pp (3.49 / old, g12), Δ`eff_cut`
+0.128 ± 0.022 / +0.120 ± 0.014; the figures of merit improve (`res_cb` −0.026 / −0.033, `leak_cut` −0.37 / −0.50, `tail_hi3`
−0.012 / −0.012, all ≥ 3 SE, replicated), oos unchanged, yield loss ≤ 0.2 %; the +0.007-0.012 % scale shift is absorbed by the
gains. On g16 calibrated it adds in full (§1). U7 rule: separation up, efficiency up → adopt. D3's `tail_lo3` guard should read
`tail_lo_f` (fixed offset) for every later variant that narrows σ.
(ii) **P12 KTeV-style correction: in, as phase-2 package P16** (§5); the term is visible on the plan's criterion (0.27-0.39 of
`res_cb` in amplitude, t² model) but its size rests on an extrapolation, so P16 starts with a measurement that can kill it.
(iii) **R_MODEL 0.3 (P8 gate): P8 go; R_MODEL 0.3 itself is not adopted.** The D5 extension replicates the core gain: legacy −0.025 ± 0.003 (2.24 GeV, 42 splits) and −0.030 ± 0.002 pp (3.49 GeV, 16), oos −0.028 / −0.034,
`leak_cut` −0.15 / −0.28 pp, U2 gap unchanged. But it fails two conditions: `tail_hi3` rises +0.010 ± 0.002 / +0.019 ± 0.004 pp
(legacy; at the adopted control g16 + log weight +0.005 ± 0.002 / +0.001 ± 0.003, the latter null), and the Ta guard (run 26138, 4 splits) worsens the gated `res_cb` by +0.083 ± 0.025 pp at g8
(`leak_cut` +0.28 ± 0.10): [I] a 0.3 amplitude floor merges real pile-up into the main pulse. The D7 P8 kill clause has therefore fired for R_MODEL 0.3
(3.3 SE, on 4 splits of one run). Most of the gain is already taken by the
calibrated gate (−0.013 pp left at g16 + log weight). P8's search is therefore for an acceptance that removes the tail-locus extras
(dt 8-48 ns after the main pulse, A2/A1 ≈ 0.05; P3A, exotic-waveform survey) and keeps pre-pulses and A2/A1 ≥ 0.1, judged on
`tail_hi3` and the Ta guard as well as `res_cb`; P8's go rests on that untested acceptance.
(iv) **P4a phase correction: P10 go.** +0.058 ± 0.007 pp (legacy) and +0.063 ± 0.007 (g12) when switched off, `leak_cut` +0.44 pp:
the correction the benchmark applies is worth 0.06 pp, and production (`mf_pulse.py`) does not apply it. Run 25320 only; the D5
multi-run extension is not needed for the budget (the term is already in every number above) and is left to P10's first slice.
(v) **P11 go, low priority.** Gate met: in T, 7.6 / 14.8 / 27.5 % (old 29.9 %) of clusters have a lead-glass coincidence and their
median E/E' is 0.804 / 0.877 / 0.923 against 1.001 / 0.979 / 0.981 without (P2). T is 0.7-1.1 % of the ep clusters; the cheap first
step is a report-level fiducial cut on T_lg against the lead-glass-gain option (U1).
(vi) **P1 residual 0.088 pp: close P1 with the residual named.** The steps that differ are the gain solve (−0.080 here vs −0.18 old)
and the Møller step (+0.006 vs −0.010). The D6 old-chain rebuild (about 3 CPU-h + a best-50 bank) is not run: it can only confirm a
note-only number, while P6 ports the old solver into BM and measures the same thing out of sample. Keep the D6 rebuild as P6's fallback.
(vii) **Asymmetric time window (tight early, loose late side): go, as P17, cheap and Ta-led.** The calibrated seed − member Δt
(`figures/P15_seed_member_dt.png`) has a core of 0.7-1.45 ns, a late-member tail beyond 8 ns that grows with energy (H2 0.06 / 1.75 / 3.06 %)
and sets the H2 plateau, and a flat accidental floor on both sides that calibration does not remove. A symmetric gate must open both
sides to keep the late tail, so it admits the early-side floor. Expected effect: small on H2 (the plateau moves ≤ 0.01 pp from 12 to
24 ns, so the early-side floor is worth ≤ 0.01 pp there), largest on Ta, where members later than 8 ns are 4.2 % (H2 0.06 / 1.75 / 3.06 %) and the early side beyond
8 ns 2.9 %, and the best symmetric gate (4-8 ns) cuts the late tail that H2 needs. The late tail may be time walk of small members (it grows with energy);
P17 also tests a walk-corrected member time against the loose late side.

## 5. Phase gate (D7) and phase-2 order (D13)

| package | gate metric (D7) | measured | ruling | order, deps |
|---|---|---|---|---|
| — adopt now | P15 gate, P4b | §1 | best calibrated gate (§1: 20 / 16 / 16 / 16 ns H2, 8 ns Ta) and log weight R7ri at 2.24 / 3.49 GeV become the phase-2 control (report constants; one P0-style default change, then `pre/` rerun) | 0, before any child plan |
| P8 secondary acceptance | P3 cell excess ≥ 2 SE and a rescan confirms | cells ≥ 2 SE (P3A); R_MODEL 0.3 replicated on 5 + 2 runs: legacy −0.025 / −0.030 pp, −0.013 on g16 + log weight; `tail_hi3` +0.010 / +0.019 (worse), Ta g8 +0.083 ± 0.025 | go on the P3 clause; the D7 kill fired for R_MODEL 0.3 (Ta g8 +0.083 ± 0.025, 3.3 SE, 4 splits), so P8 rests on an untested acceptance (ruling iii) | 1; deps P3, P5 |
| P9 edge correction | edge − centre ≥ 0.15 pp at 2.24 and 3.49, or an edge tail excess ≥ 2 SE | core +0.083 ± 0.010 / +0.117 ± 0.013 (g12; not met); `tail_lo_f` +3.2 / +4.7 pp, `eff_cut` −4.4 / −5.7 pp (≥ 25 SE; met) | go, on the tail clause; target `eff_cut` and `tail_lo_f` at the edge; core gain ≤ 0.03 / 0.05 pp | 2; deps P2, P4, P5; refit after P16 if both adopted |
| P6 gain solve (comparison arm, U2) | P1 closed | closed with the 0.088 pp residual named (ruling vi) | go; carries the D6 rebuild as fallback | 3; deps P1, P5 |
| P7 beyond the pool | a region with ≥ 2000 ep, σ − pool σ ≥ 0.15 pp, or `tail_hi3` above pool by 2 SE | M: 67 k / 27 k clusters, +0.17 / +0.40 pp, `tail_hi3` 1.23 vs 0.16 % (2.24 GeV); H: 2.9 M clusters, +0.49 pp | go; start with H and M (most clusters) | 4; deps P2, P5, P6 |
| P16 below-threshold correction (new) | P12: loss term ≥ 10 % of σ_E in quadrature | 0.27-0.39 (t²), 0.18-0.64 over models | go for slice 1 (measurement), slice 2 only if slice 1 confirms | 5; deps P12, P5 |
| P17 asymmetric time window (new) | P15: asymmetric Δt (late tail 1.8-4.2 % beyond 8 ns, flat early floor) | H2 plateau flat ≤ 0.01 pp over 12-24 ns; Ta best symmetric gate 4-8 ns, −0.30 pp vs legacy | go; slice 1 emulates [−a, +b] report-side at the loose gate (C++ `seed_t_win` is symmetric; `Amend parent:`); acceptance D3 on Ta and ≥ 2 H2 energies, kill if null on both | 2b; deps P15, P5; cheap (1 switch, about 10 reports) |
| P10 I(φ) in `mf_pulse.py` | P4a \|Δ\| ≥ 2 SE and ≥ 0.02 pp | +0.058 ± 0.007 (run 25320) | go; first slice = D5 extension (shared control: existing scans) | 6; deps P4, P5; production only, does not change the benchmark |
| P11 boundary and lead glass | ≥ 1 % of T with a coincidence and ≥ 0.5 % deficit | 7.6-29.9 %, deficit 6-20 % (relative) | go, low priority (ruling v) | 7; deps P2, P5 |
| P14 snake-scan / cosmic | P13 finds readable runs | runs found, data at JLab only | blocked on the user (data copy, scan log, HV history, tagger) | after the data arrive |

D13 order (effect on `leak_cut`, `tail_hi3`, `sep_pi` at 3.49 GeV): the adoptions first (`sep_pi` +0.08, `leak_cut` −0.5 pp
measured), then P8 (bound −0.013 pp, `leak_cut` −0.14 / −0.18 pp at the adopted gate), P17 (2b: no measured D13 effect yet, placed early because it is cheap and Ta-led), P9
(efficiency at a third of the pool), P6 → P7 (coverage: `sep_pi` 1.8-2.2 → pool's 2.6 over regions holding 2-3 × the pool's
clusters), P16, P10, P11. P9 and P6 touch disjoint code (position correction vs gain solve) and can
run in parallel; P7 waits for P6.

**P16 (new, `Amend parent:`): below-threshold energy correction (KTeV; U6).**
- Goal: recover the energy and the event-to-event width lost to the 9.5-9.7 MeV readout threshold (5-13 MeV across modules).
- Scope: slice 1 (measurement, report-level, no BM default change): split the pool by the members' effective threshold from
  `variants/P12/raw/edge_*.npz` (modules with edge < 8 vs > 11 MeV, matched in radius) and compare `res_cb`, `mu_none`, `tail_lo_f`;
  this tests the P12 extrapolation without a lower-threshold run. Slice 2 (only if slice 1 shows a threshold dependence ≥ 2 SE): a
  report switch adding the KTeV expected below-threshold energy per cluster, gains and position correction refitted after it.
- Acceptance: D3 adoption out of sample at ≥ 2 H2 energies (expected largest at 0.73 GeV), `eff_cut` and `tail_lo_f` not worse by
  > 2 SE, yield loss ≤ 1 %; U2 gap rule for any fitted constant. Kill: slice 1 null at every energy.
- Owns: a new BM switch (P0 `Amend parent:` route), `ER/improvement/P16_*`, `variants/P16/`. Deps: P12, P5. Cost: 3 reports, then low-medium.

## 6. Caveats

- **No re-baseline needed (D5).** The 2.24 GeV bank rebuild (15dfd867 → 7f125b75) leaves the hits unchanged: `h_e`, `h_t`, `h_a`,
  `h_extra`, `h_k` are bit-identical between the old scan and a current-code rescan for run 25320 (C0) and run 25228
  (`scan_variants/rebaseline/2p24GeV_H2/run_25228.npz`; 2 of 5 runs checked, on which D5's same-bank rule is relaxed); the 3.49 GeV and Ta banks are unchanged (294184b0, 8faf5d6f). `pre/` is the
  phase-2 control; the existing scans served as the extension's control arm (§7).
- **C0 vs pre/ on run 25320** (single-run report 1.882 / oos 1.906 vs the run's rows in the five-run report 1.867 / 1.884,
  legacy): a single-run report on the old scan reproduces C0 to every digit, so the difference is the training sample of gains and
  position correction (one run vs five), not the bank. Single-run Δ (P3B, P4a) are on that basis.
- **Run 24937** (3.49 GeV): four 2.1 GB files hold 24 k events each (157 k in 24698's files of the same size), only 389 clusters
  reach the gain solve, and every channel is anomalous (exotic-waveform survey). Its splits are NaN in every report; the 3.49 GeV
  numbers rest on 24698 and 24911 (16 splits). Recommend dropping it from D4.
- Quadrature sums assume independent terms; the measured combinations (§1) are the numbers to quote. Single-run and per-run-type
  training means split SEs understate run-to-run spread for fitted quantities (plan Risks); replication on 3.49 old is the guard.
- P15's offsets come from WaveAnalyzer times (cross-basis); lead glass keeps offset 0, so T and lead-glass rows are not compared
  across raw and calibrated times.

## 7. Runs made by P5 (banks: 2.24 rescans 7f125b75, existing 15dfd867; 3.49 294184b0; Ta 8faf5d6f)

- Reports `toff_logw` (calibrated times + log weight R7ri, existing scans, gates 5-20 ns); R_MODEL 0.3 rescans of 25228, 25257, 25260,
  25299, 24698, 24911, 26138 (+ 25320 from P3) with reports `ext_ctrl` / `ext_rmodel` (control = existing scans), `toff_logw_rmodel`,
  `toff_logw_2run`; identity rescan of 25228. Scripts and wall times: `variants/P5/{scan_lane,run_combo,run_ext}.sh`, `logs/`.
