# AMP P1b: where the raw-peak excess comes from (run 25320, 2026-11-14)

Follow-up of P1 (`AMP_P1_estimators.md`): the raw peak sample (`ENERGY_DEF=rawpk`, CFD time, no tables) is +0.456 +- 0.011 pp worse than the integral on run 25320 (2.284 against 1.827 % out of sample; quadrature excess 1.37 %). Question: noise, shape, pile-up, low-energy members, or calibration? Benchmark and every default untouched (sha1 86cf8bb2... as at the P1 scans; no new `OVERRIDABLE` entry, v3 byte-identical by construction, `check_identical.sh` not needed). Tools `tools/p1b_*.py` (benchmark wrappers, not edits), tables `AMP_P1b_tables/`, arms `output/energy_resolution/variants/AMP/p1b_25320/`, driver `p1b_reports.sh`, numbers `AMP_P1b_excess.tsv`, `AMP_P1b_arms.tsv`, `AMP_P1b_arms_contrasts.tsv`.

## Verdict

1. **The excess is not single-sample noise.** Pedestal noise explains about 1 % of it. The measured excess is a 1/sqrt(A) fluctuation of the peak sample (about 0.44 of 2.0 percent^2), an amplitude-proportional shape term (0.18), the sampling-phase curve (0.83, recoverable), and pile-up members (0.2); about 0.3 is unexplained and overlaps the cost of the gain solve with a noisier estimator (arms 2.00 against 1.56 in the gain-solve-free form).
2. **About 41 % of it is recoverable by calibration** (the module phase curve against the integral, 0.83 +- 0.09 of 2.00 +- 0.11 percent^2; an amplitude table adds 0.0 +- 0.04 in the report arms; the gain-solve-free form gives 46 % for both tables together, 0.65 + 0.07 of 1.56). With that table in sample the raw peak is still +0.28 pp behind the integral (2.146 against 1.866 % Crystal Ball, in sample, one split).
3. **The remaining 1.1-1.2 percent^2 is intrinsic to a single sample**: the same table on a multi-sample fixed sum (template free) recovers most of it: 5-sample sum from the peak sample (`sum5`, [p-1, p+3]) 1.982 %, 7-sample sum ([p-2, p+4]) 1.960 %, and the best, a 7-sample sum anchored at the CFD sample (`cfd7`, [jc-1, jc+5], jc the last sample below half the peak) 1.935 %, which is +0.070 pp behind the table-free integral (1.866). None reaches the integral, and v3 with its tables is another 0.07 pp better (P1: 1.827 -> 1.761 % out of sample), so the best variant is about +0.14 pp behind control v3 (whose tables are in sample for 25320, as are these).
4. **Recommendation: stop the raw-peak track; do not start the P2 rescan as planned.** If the user still wants a template-free amplitude chain for independence, use a 7-sample sum anchored at the CFD sample (or the peak-anchored 5-sample sum, same within errors) with a per-module phase table and expect to lose 0.14-0.18 pp against control v3 (to be confirmed out of sample; see Caveats). Without the table every fixed-window sum is far worse than the raw sample (`cfd7` 3.76 %), because the sums have a much larger phase dependence.

## Data and method

- A rescan of **one file** (run 25320, split 50: 159 596 physics events, 2.52 M hits) with `tools/p1b_scan.py`: the benchmark's own scan with `--set PEAK_COLUMNS=True` plus extra per-hit columns from the same subtracted waveform and peak sample p (fixed sums `h_s3/4/5/7` over [p-1,p+1], [p-1,p+2], [p-1,p+3], [p-2,p+4], CFD-anchored sums `h_c5/c7`, pedestal std and lag-1 correlation of the waveform, template sum of squares, peak-sample offset, 16-sample windows of 218 k hits with the fitted model). Every benchmark column except `h_apk_par0` (dropped for the hit id) is bit-identical to the P1 scan of that split (checked). Wall 555 s; 4-file scans were 3x slower under machine load and were stopped (see Deviations).
- Hit level: single-pulse hits (h_k = 1, not extra), h_a >= 30 ADC; reference h_a (an efficient estimator, so var(col / h_a - 1) is the variance in excess of the matched filter); per-module phase curve (16 bins, from >= 2000 hits, else pooled) removed; fit v(A) = (n/A)^2 + a/A + c^2 to the binned robust variances (`fit.tsv`; the 4-split P1 lane gives n 0.69, a 0.100, c 0.0062: same).
- Cluster level, two ways. **Arms**: the benchmark report (gain solve, position correction, gate 16 ns, no M1/P16 tables; the non-integral arms use the CFD time and the integral control the matched-filter time, a null effect of +0.002 pp in P1; per-module `AMP_SCALE` = the 8-split raw-peak solve of P1 for every column, `p1b_scale.py`) for the integral and each variant, compared on the 19 023 elastic clusters present in every arm (|e/E' - median| < 7 %): robust sigma (1.4826 MAD) of e/E', difference of variances, standard error from 200 resamplings of 40 event blocks (paired). **Direct** (`p1b_direct.py`, gain-solve free): for the members of the integral arm's elastic clusters look up each hit and form D = sum_m e_m (raw/integral ratio of the hit / module median - 1) / E; var(D) with the mean removed (a global offset is absorbed by the gain solve) is the excess variance, split into an independent part (sum of squared member terms) and cross terms. Total: arms 2.00 +- 0.11 percent^2 (robust sigma^2 difference; Crystal Ball: 2.321^2 - 1.866^2 = 1.91, P1 eight splits out of sample 1.88), direct 1.555 +- 0.014; the 0.44 difference is what the gain solve and position map add for the noisier estimator.
- The mix arms use the hit's integral energy as the selector (the seed is known only after clustering); "hits >= 1000 MeV" is the seed in nearly all elastic clusters (median seed amplitude 1604 ADC).

## Hypotheses

**H1, single-sample noise.** Pedestal std of one sample 0.95 ADC (module 10-90 % range 0.86-1.07; the benchmark's MAD reference 1.13), lag-1 correlation -0.37; the matched filter's white-noise error is 0.45 ADC (template sum of squares 4.47), so a raw sample carries 0.84 ADC of noise in excess. Median seed amplitude 1604 ADC, so sigma/A of the seed is 0.059 %. Propagated over the members of the elastic clusters (member energies, module scale): 0.15 % of cluster energy for all members, 0.04 % for the seed alone, **0.018 +- 0.004 percent^2 of the 2.0** (1 %). The hit-level fit finds n = 0.71 ADC against 0.82 predicted from the pedestals. Hit scatter against the white prediction (`scatter_rawpk.tsv`): A = 44 ADC 5.0 % measured against 1.8 % white; 192 ADC 2.3 against 0.42; 1213 ADC 1.0 against 0.07. Rejected.

**H1b/H2, scaling and the residual after the harmonic.** The scatter of raw/h_a at fixed A follows A^-1/2 (sigma x sqrt(A) = 0.33-0.38 over 45-1600 ADC; 0.49 in the sparse 2000-4000 ADC bin), not 1/A: a = 0.100 ADC (a/A), c = 0.60 % (A-independent), n small. The phase harmonic is not the residual: after removing each module's curve the scatter keeps both terms. Peak-sample position jitter is not it: 79-86 % of hits have the template-peak sample as the largest, and hits whose largest sample is one later scatter the same (3.6 against 3.4 % at 50-200 ADC; 0.93 against 0.99 % above 1000 ADC; `peakpos.tsv`). The a/A term is a fluctuation of the sample against the fitted template that grows as the pulse (photoelectron-count scale), is averaged out by a multi-sample window (a: raw 0.100, sum3 0.029, sum4 0.0087, sum5 0.0041 ADC; c: 0.60, 0.55, 0.39, 0.24 %) and so is intrinsic to one sample.

**H3, pile-up / U2.** Pile-up is not the cause. Clusters without a pile-up flag (no member with two in-window pulses, 73 % of clusters) have a Crystal Ball excess of 1.44 % (2.337 against 1.843), the same as all clusters (1.38 %); the numbers, with the threshold arms, are in `AMP_P1b_thr_pu.tsv`. The scatter of h_k >= 2 primaries at fixed amplitude is not larger than for h_k = 1 (`classes.tsv`). Members with a second pulse or extra pulse carry 0.20 +- 0.04 percent^2 (arms; raw only in those hits) or 0.12 (direct), 6-10 % of the total.

**H4, seed against members** (raw only for the hits named, integral for the rest, paired, percent^2 of 2.00 +- 0.11):

| raw only for | excess |
|---|---|
| hits >= 1000 MeV (the seed) | 1.43 +- 0.09 |
| hits < 1000 MeV | 0.61 +- 0.08 |
| hits >= 400 MeV | 1.66 +- 0.11 |
| hits < 400 MeV | 0.41 +- 0.05 |
| single-pulse hits | 1.88 +- 0.10 |
| pile-up hits | 0.21 +- 0.04 |

The pairs add up to the total within errors (2.04, 2.06, 2.09). Seventy percent of the excess enters at the seed, which is the opposite of a low-S/N story. In the direct form: members >= 1000 MeV 1.04, < 1000 MeV 0.43 (sum 1.47 against 1.56).

**H5, member threshold.** Removing hits below 14, 19 and 25 MeV in both arms leaves the gap unchanged: Dsigma 0.442, 0.439, 0.454 pp (0.453 at the default), excess variance 2.03, 2.05, 2.21 percent^2 (2.00); the resolution of both arms degrades (integral 1.98 -> 2.07, 2.12, 2.21 % robust). Low-energy members (< 100 MeV) carry 0.19 of 1.56 in the direct form (12 %), and a higher threshold costs as much energy as it saves noise. Rejected. (A 38 MeV threshold fails in the position-correction map and was dropped; the nominal threshold of the scan is not 9.5 MeV exactly, so the three values are absolute.)

## Budget (percent^2 of cluster energy, robust-sigma basis; total 2.00 +- 0.11, direct 1.555 +- 0.014)

| cause | percent^2 | recoverable by |
|---|---|---|
| 1 sampling-phase curve of the raw peak against the integral (per-module 16-bin table; arms) | 0.83 +- 0.09 (direct 0.65) | calibration (P5 M1-type table) |
| 2 amplitude-dependent bias raw/integral (0.9 % swing), after 1 | 0.0 +- 0.04 (direct 0.07) | calibration, small |
| 3 pile-up members (h_k >= 2 or extra pulse; arms, no tables) | 0.20 +- 0.04 (direct 0.12) | estimator or fitter |
| 4 single-sample fluctuation ~ 1/sqrt(A), predicted from the hit fit | 0.44 +- 0.07 | intrinsic to one sample |
| 5 amplitude-proportional shape term c = 0.6 % | 0.18 +- 0.03 | intrinsic to one sample |
| 6 white pedestal noise | 0.013-0.018 | intrinsic |
| sum of 1-6 | 1.67 +- 0.14 | |
| unexplained (total - sum) | 0.33 +- 0.18 | |
| cross terms of the member shifts after both tables (direct, intra-cluster correlation) | 0.007 +- 0.007 | none: errors are independent between members |
| arms total minus direct total (gain solve and position map with a noisier estimator; not isolated) | 0.44 +- 0.11 | more calibration statistics |

The uncertainties on 4 and 5 are 15 % modelling (fit weighting, binning; the 4-split P1 lane agrees to 1 % in a and 3 % in c); on 1 and 3 the paired bootstrap. The hit-level fit was made after the per-module phase curve was removed, so causes 4-6 do not contain cause 1. After both tables the direct remainder is 0.835, all of it independent member scatter: 0.63 is causes 4-6, 0.08 pile-up members, 0.12 other (tails, members below 30 ADC). Causes 4-6 were propagated over all members, pile-up members included, and cause 3 is measured without tables, so the sum of 1-6 double counts the pile-up members' share (about 0.1) and the phase part of their 0.2; the unexplained 0.33 is therefore a lower bound on the part that is not accounted for, and the 0.44 arms-minus-direct row (calibration cost) is its natural home, not an intra-cluster correlation.

**Calibration or intrinsic.** Recoverable by P3-P5 type calibration: 0.83 +- 0.09 percent^2 (41 %; direct 46 %), all of it the phase table; with it the raw peak sits at 1.17 +- 0.06 (arms), +0.28 pp Crystal Ball above the integral. Intrinsic to the single sample: 0.6-0.7 (32 %, causes 4-6), plus pile-up 0.1-0.2, plus the 0.3-0.4 that sits in the calibration solve (unexplained).

## Estimator variants (Crystal Ball sigma/E %, g16, in sample, one split; Dvar = variance excess over the integral, percent^2)

| estimator | no tables | Dvar | with phase + amplitude tables (reference: integral energy) | Dvar |
|---|---|---|---|---|
| integral (control, no tables) | 1.866 | 0 | | |
| raw peak sample | 2.321 | 2.00 +- 0.11 | 2.139 | 1.18 +- 0.07 |
| 3-sample sum [p-1, p+1] | 3.151 | 10.4 | 2.112 | 1.53 +- 0.07 |
| 4-sample sum [p-1, p+2] | 2.516 | 4.12 | 2.034 | 0.81 +- 0.06 |
| **5-sample sum [p-1, p+3]** | 2.146 | 1.25 +- 0.10 | **1.982** | **0.40 +- 0.06** |
| 7-sample sum [p-2, p+4] | 2.879 | 4.24 | 1.960 | 0.49 +- 0.06 |
| CFD-anchored sum [jc, jc+4] | 5.00 | 19 | 1.985 | 0.63 +- 0.05 |
| **CFD-anchored sum [jc-1, jc+5]** | 3.757 | 9.8 | **1.935** | **0.37 +- 0.05** |

The sums have a larger phase dependence (4.9 % peak to peak for `sum5` against the integral, 2.6 % for the raw sample), so they need the table more; their no-table values are not a verdict on the sum. `cfd7` and `sum5` with tables agree within errors (difference 0.03 +- 0.06). A table from the 8-split P1 data for the raw peak and from the one split itself for the sums (in sample; pooled phase curve for modules with < 2000 hits) was used. The raw-peak report arm with the h_a-referenced table (the P1 style reference, `corrphase`) recovers only 0.48 +- 0.06: the integral's own phase correction ph(x) is inside h_e, so the table must be referenced to it. The v3 integral with its tables is 1.761 % (P1, eight splits, tables in sample for 25320, gain and position out of sample) against 1.827 % without them, a 0.066 pp gain; adding it to the table-free gaps puts `cfd7`, `sum7` and `sum5` with tables at +0.14, +0.16 and +0.18 pp behind control v3 (one split against eight: only the differences are usable).

Without tables the CFD-anchored sums are poor (the window position steps by one sample with the sub-sample phase, so the pulse content of the window changes: a phase curve of several percent); with the E-referenced table `cfd7` is the best variant. That mechanism is a reading of the numbers, not a separate test.

## Recommendation

**Stop the raw-peak track (do not run P2 as planned).** After every calibration the table can give, the raw peak stays +0.28 pp behind the integral and about +0.35 pp behind v3 (adding the 0.07 pp that v3's tables gain on the integral); no member threshold, no pile-up cut and no timing choice changes that. The best template-free estimator found, the CFD-anchored 7-sample sum with a phase table, still loses 0.07 pp to the table-free integral and about 0.14 pp to v3. The deficit is the shape information that the matched filter uses and a fixed window does not (a and c in the hit scatter); the per-module phase curve and its amplitude dependence are the only calibration-recoverable parts.

If an amplitude chain is still wanted for independence (a template-free cross-check, not a resolution gain): carry `h_c7`, `h_s5` and `h_s7` (fixed sums) in the P2 rescan instead of the parabola, run P2 with the h_e-referenced phase table, and treat about +0.14-0.18 pp against control v3 as the expected cost; the decision is the user's.

## Caveats

- One split (50), in sample, no out-of-sample column; the 8-split P1 numbers (out of sample) agree with the totals (Crystal Ball 1.91 here against 1.88). Standard errors are paired bootstraps over 40 event blocks of that split, which do not include split-to-split variation (the P1 eight-split standard error on the raw-integral gap is 0.011 pp on a 0.456 pp gap).
- Tables of the sums are fitted on the same split they are tested on; the raw-peak table comes from 8 splits that include it. An out-of-sample test needs a second file (another 0.15 CPU-h) or the P2 data.
- The hit-level fit is on single-pulse hits with h_a >= 30 ADC; members below 30 ADC (5.4 % of the member energy, 0.98 MeV per ADC) are extrapolated.
- "Unexplained" 0.33 +- 0.18 is not distinguished from zero at 2 sigma; the direct form finds no intra-cluster correlation (cross terms 0.007 +- 0.007), so the arms-minus-direct 0.44 (gain solve and position map) is the candidate; it was not isolated.
- The excess is measured against the integral, which is itself 4.2 permille off in phase against E' (P1): a table referenced to the integral transfers that residual to the raw peak.

## Deviations from the brief

- CPU: the scan used 555 s (0.15 CPU-h) for the kept file, plus about 0.11 CPU-h for a four-file attempt stopped after six minutes (machine load 13, three times slower than the P1 scans) and 0.04 CPU-h of profiling: about 0.3 CPU-h in total for scans, at the cap, and about 43 single-split reports of about 20 s each (0.24 h wall, mostly one worker); so the P1b is on one split of eight, not on the P1 sample.
- Report arms ran on one split, so `res_cb_oos` does not exist; in-sample Crystal Ball and the robust sigma are quoted.
- Hit-level scatter was also run on the 8-split P1 scans (`AMP_P1b_tables/p1_lane_a`, no pedestal columns).
- Mix arms select by hit energy (the seed is not known before clustering); the member-threshold arms use absolute 14, 19, 25 MeV, not multiples of 9.5 MeV, and 38 MeV failed.
