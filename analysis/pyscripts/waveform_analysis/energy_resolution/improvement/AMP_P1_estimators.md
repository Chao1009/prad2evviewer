# AMP P1: peak estimators on run 25320 (2026-11-14)

Plan `notes/plan_amplitude_calibration_20261114.md`, package P1 (U1, U2, D2). Scan: run 25320, all 8 splits (50-400), bank 7f125b75, `--set PEAK_COLUMNS=True`, two lane files `scan_variants/amp_p1/2p24GeV_H2/run_25320_{a,b}.npz` (h_t/h_e/h_a/h_k/h_extra byte-identical to `scan_variants/c0`). Tables: `AMP_P1_tables/*.tsv` (tool `tools/amp_estimators.py`), report arms `output/energy_resolution/variants/AMP/p1_25320/<arm>/2p24GeV_H2` (driver `p1_reports.sh`), summary `AMP_P1_estimators.tsv`. Tests are in sample at this stage except the `oos` columns (even/odd cross fit of the report).

## Verdict

**Raw peak sample.** The 3-point parabola is the worse estimator on both criteria and no user action is needed (Q1 rule: the two do not disagree).

- Phase dependence against the true energy (clusters with 0.95 < E/E' < 1.05, first harmonic of E/E' against the seed onset phase, per unit member fraction, as P10): raw 11.9 permille, parabola 14.2, v3 integral with the phase table but no M1 4.2 (P10: 9.3 with no correction at all). Residual cluster rms of the fitted k = 1, 2 curve: raw 6.0 permille, parabola 7.4, integral 2.3.
- Phase dependence against the template-fit amplitude `h_a` (single-pulse hits, `h_a` > 200 ADC, each module divided by its median; median over 128 modules with >= 2000 hits): raw first harmonic 0.67 %, peak to peak 1.8 %; parabola 1.6 % and 3.9 %; the integral's own correction curve `ph(phi)` (`h_e / h_a`) 0.79 % and 2.1 %. `h_a` carries a phase bias of its own (7.9 permille for the integral here against 4.2 against E'), so these curves rank the estimators and the E/E' fit sizes them. Both references give the same ranking, and the parabola is about twice as bad as raw on the `h_a` curves (its binned curve is a saw tooth, -0.8 % to +2.5 %, a bias of the 3-point interpolation on a 4 ns sampled, sharp peak).
- Resolution, paired out of sample over 8 splits at the 16 ns gate (no M1, no P16, per-module `AMP_SCALE` from pass 1): `res_cb_oos` integral 1.827 +- 0.006 %, raw 2.284 +- 0.012, parabola 2.473 +- 0.011; parabola minus raw +0.189 +- 0.021 pp (paired, `d_res_cb_oos` in `delta.tsv`; in sample +0.183 +- 0.019); the guards also favour raw (`tail_lo_f` +0.106 +- 0.023 and `leak_cut` +1.33 +- 0.09 for the parabola); the reverse-sign metrics are the 3-sigma tails (`tail_hi3_oos` -0.022 +- 0.007, `tail_lo3` -0.34 +- 0.04), which are in units of each arm's own wider sigma.

Neither estimator is competitive before calibration: raw is +0.456 +- 0.011 pp worse than `integral_mf` (same selection, no tables) with `tail_lo_f` +1.67 +- 0.03, `eff_cut` -3.0, yield -1.4 %. About 6 permille of its 21.1 permille cluster width is the phase term that an M1 table (P5) can remove (integral: 4.2 permille harmonic; the v3 tables, M1 and P16 together, take the cluster width from 18.4 to 18.0 and `res_cb_oos` from 1.827 to 1.761); the rest of the excess is not phase. Whether an amplitude track can still win is a P5/P6 question, not settled here.

## Tables

Estimator table (`estimator_table.tsv`; ratio = column / `h_e` module medians, 128 modules; scale s = median of `h_e` / column over non-extra hits with `h_e` >= 200 MeV):

| estimator | s (MeV/ADC) | ratio median | IQR over modules | phase harmonic (pooled / module median) | peak to peak (pooled / module median) |
|---|---|---|---|---|---|
| integral (`h_e / h_a`) | 1 | 1 | 1-1 | 0.63 / 0.79 % | 1.5 / 2.1 % |
| raw (`h_apk_raw`) | 0.998 | 1.043 | 0.997-1.086 | 0.61 / 0.67 % | 1.4 / 1.8 % |
| parabola (`h_apk_par`) | 0.983 | 1.060 | 1.014-1.104 | 1.45 / 1.63 % | 3.2 / 3.9 % |

The 4-9 % module spread of column / `h_e` is the per-module `gain * iint` the global factor does not carry; pass 2 (`AMP_SCALE`) absorbs it but changes `res_cb_oos` by only -0.02 (raw 2.304 -> 2.284 at 16 ns), so it is not the cause of the gap. `energy_bins.tsv` gives the ratio in 8 energy bins.

Other arms (`AMP_P1_estimators.tsv`, `res_cb_oos` %, 16 ns gate): `integral_cfd` 1.830 (+0.002 +- 0.002 against `integral_mf`, null; legacy path identical by construction), pass 1 of raw 2.304 and parabola 2.513, `v3chain_integral_mf` (v3 default tables, run 25320 in the table's training set) 1.761. Offsets are `time_calib.json` (MF basis) for the CFD arms.

## U2 subtraction (`overlap.tsv`)

Primary hits with `h_k` >= 2 against the signed time to the module's largest extra pulse (`h_t` secondary minus primary), 254 k tabulated of 266 k such hits with `h_a` > 200 ADC in the 8 splits (4.5 % have no stored extra pulse; the table spans |dt| < 60 ns). Size of the subtraction `h_apk_par / h_apk_par0 - 1` (median) and column / `h_a` relative to the module's single-pulse median, with and without the subtraction:

| dt (ns) | n | subtraction | raw rel | par rel | par0 rel (no subtraction) |
|---|---|---|---|---|---|
| -60 .. -30 | 693 | -1 .. -2 % | 1.000 | 1.002 | 1.02 |
| -30 .. -14 | 868 | -3.5 .. -10 % | 1.006 | 1.006 | 1.05 .. 1.14 |
| -14 .. 0 | 1143 | -25 .. -28 % | 1.057 | 1.050 | 1.38 .. 1.44 |
| 0 .. 8 | 6175 | -3.8 % | 1.004 | 1.008 | 1.053 |
| 8 .. 14 | 78447 | -0.07 % (13 % of hits > 1 %) | 1.003 | 1.004 | 1.008 |
| >= 14 | 167 k | 0 | 0.990-0.997 | 0.990-0.998 | same |

Zero for `h_k` = 1 by construction (par equals par0 exactly). The subtraction removes the secondary's contribution to within about 5 % of the module's single-pulse value at close separation (an unsubtracted peak is 40 % high at -14 .. 0 ns, a preceding pulse). Primaries with a secondary 14-60 ns later sit 1 % below single-pulse hits for all three columns (a selection or `h_a` effect of joint-fit rows, not the subtraction). Cluster-level effect of the subtraction not run (no `parpk0` energy switch).

## CFD against MF time (`cfd_vs_mf.tsv`)

`h_tcfd - h_t` after removing each module's median (module medians: rms 0.13 ns across modules, range -0.47 .. +3.0 ns; one module with 3 ns). Robust sigma 0.23 ns (`h_a` 50-200), 0.15 (200-1000), 0.13 (> 1000); |dt| > 1 ns for 0.2-0.5 % of single-pulse hits, > 4 ns for 3e-7. `h_k` >= 2 primaries: sigma 0.14 ns, 0.4 % > 1 ns, 0.2 % without a CFD (no crossing). Extra pulses: sigma 0.44 ns, 12 % > 1 ns, 2.5 % > 4 ns, 9.6 % without CFD. CFD fallback (NaN `h_tcfd`) over PbWO4 hits: 5e-5. The distribution is narrow because both times are on the same pulse; it says the two definitions agree per hit, not which is better as a gate time (that is P3/P4, on the pair network).

## Implementation facts and deviations

- CFD of a hit is computed on the waveform with the other fitted pulses subtracted, with the peak sample of the same waveform; the plan (D2, XS) wrote the unsubtracted waveform. Reason (critic): on the unsubtracted one an extra pulse inside about 24 ns of its primary gets the primary's time (14.6 % of extras beyond 15 ns); for single pulses the two are identical.
- `PEAK_COLUMNS` switch (default False) controls the scan columns ("when enabled"); `TIME_OFFSETS_FILE` and `AMP_SCALE` default to `""` instead of None (the `--set` type check).
- Non-integral reports apply the v3 `BELOW_THR_BY_RT` and `PHASE_RESID_BY_RT` tables unless they are emptied: the P1 arms empty them (`--set "BELOW_THR_BY_RT={}"`), and P2 onward must do the same until P5 builds the amplitude tables.
- The phase test of the plan (`E_est / h_e` against onset phase) is run as column / `h_a` against onset phase plus the P10-style E/E' fit; `h_a` is the reference of the first.
- Baseline of the peak estimators is the matched filter's start baseline (pedestal, or `b_ref` for 3.0 sigma outliers), not the fitted free offset.
