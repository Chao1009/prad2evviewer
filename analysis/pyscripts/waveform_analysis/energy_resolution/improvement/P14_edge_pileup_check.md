# P14.2b: all-event cross-check of the endpoint edge against the pile-up rejection (2026-11-14)

Script `snake_scan/edge_pileup_check.py` (imports `fit_edge`, `make_hist`, `load_calib`, `neighbours` from `edge_pass.py`, which is unchanged). Table `P14_edge_pileup_check.tsv` (one row per module and setting, plus the transfer spreads as `#` lines). Run 23701, the 9 local files, same dwell ranges, map and calibration as wave 1. One decode per file (2800 CPU s in total) stores, for every dwell event with a central peak and for each of the 9 members, the number of peaks and the time, height and integral of up to 3 peaks; all settings are then applied offline and fitted with the wave-1 two-pass procedure (pass 2 divides the neighbours' W by their s_1). The `w1` setting reproduces `edge_23701.tsv` exactly (s, n_single, n_pileup of all nine modules).

## Settings

Rejection of edge_pass: `collect_events` flags an event as pile-up when any fired member of the 3x3 has `len(peaks) >= 2` from `WaveAnalyzer.analyze` (peak threshold `max(peak_nsigma x ped rms, min_peak_height)`, no height or time condition). At the scan rate, secondary peaks sit everywhere in time (histogram of signed time to the in-time peak, 20k W599 events: broad humps at -150..-75 ns and +25..+225 ns, none within +-25 ns) and 25 % of the central modules have a secondary above 200 ADC, so most flagged secondaries are separate, out-of-time bunches that do not touch the in-time integral. The two looser settings reject only secondaries that can overlap the in-time pulse. The in-time peak is the stored peak closest to T0 = 158 ns (median central peak time).

| setting | rejection | value of a member | rejected, 9-module range |
|---|---|---|---|
| `none_sum` | none | sum of all stored peak integrals (edge_pass convention for >= 2 peaks) | 0 % |
| `none_main` | none | in-time peak integral only | 0 % |
| `w1` | wave 1: any member with >= 2 peaks | its integral | 29-54 % |
| `loose1` | secondary peak >= 50 ADC within 100 ns of the in-time peak, any member | in-time peak integral | 7-17 % |
| `loose2` | secondary peak >= 200 ADC within 60 ns, any member | in-time peak integral | 4-9 % |

`loose2` is the one that matches the plan's expected rate of about 10 % (W599 sample: 10 %); `loose1` is the intermediate setting (150 and 50 ns windows and 20 ADC thresholds were tried in the W599 sample: 9-38 %). Only 3 peaks per member are stored; events with more are cut at the 3 earliest.

## Edge per module (3x3 convention, MeV of the production W; s = E_edge / 2108.85)

E_edge in MeV; d = shift of s relative to `w1` in %. Full numbers (s_err, window spread, sigma, n kept, fraction rejected) in the tsv.

| module | none_sum E (d) | none_main E (d) | w1 E, s, s_err | loose1 E (d) | loose2 E (d) | n kept w1 / none / loose2 (k) | sys spread % | wave-1 s_err % |
|---|---|---|---|---|---|---|---|---|
| W547 | 2107.0 (+1.19) | 2087.7 (+0.26) | 2082.3, 0.9874, 0.0031 | 2086.0 (+0.18) | 2086.8 (+0.22) | 283 / 552 / 502 | 1.19 | 0.32 |
| W552 | 2354.8 (+1.79) | 2322.4 (+0.39) | 2313.5, 1.0970, 0.0037 | 2320.1 (+0.29) | 2321.2 (+0.33) | 280 / 612 / 555 | 1.79 | 0.34 |
| W557 | 2185.2 (+1.22) | 2165.7 (+0.32) | 2158.8, 1.0237, 0.0042 | 2163.6 (+0.22) | 2164.3 (+0.26) | 217 / 388 / 355 | 1.22 | 0.41 |
| W564 | 2460.2 (+5.28) | 2344.5 (+0.33) | 2336.8, 1.1081, 0.0097 | 2333.8 (-0.13) | 2339.0 (+0.10) | 122 / 235 / 226 | 5.41 | 0.88 |
| W580 | 2045.4 (+1.09) | 2028.3 (+0.25) | 2023.2, 0.9594, 0.0031 | 2026.5 (+0.16) | 2027.4 (+0.21) | 342 / 622 / 567 | 1.09 | 0.33 |
| W585 | 2137.5 (+1.23) | 2115.7 (+0.20) | 2111.5, 1.0013, 0.0028 | 2114.6 (+0.14) | 2114.5 (+0.14) | 188 / 367 / 334 | 1.23 | 0.28 |
| W590 | 2180.2 (+1.00) | 2163.4 (+0.22) | 2158.7, 1.0236, 0.0032 | 2161.7 (+0.14) | 2162.5 (+0.18) | 103 / 189 / 172 | 1.00 | 0.32 |
| W594 | 2121.0 (+0.74) | 2106.6 (+0.05) | 2105.5, 0.9984, 0.0033 | 2106.3 (+0.04) | 2106.4 (+0.04) | 247 / 350 / 338 | 0.74 | 0.33 |
| W599 | 2197.3 (+1.27) | 2176.5 (+0.31) | 2169.8, 1.0289, 0.0040 | 2174.4 (+0.21) | 2175.3 (+0.25) | 352 / 662 / 600 | 1.27 | 0.39 |

(`sys spread` = (max - min of s over the five settings) / s_w1; it equals the none_sum shift in eight modules. For W564 it is 5.41 % against a none_sum shift of 5.28 %, because `loose1` lies 0.13 % below `w1`.)

## Systematic

- `none_sum` is not a usable estimate: summing out-of-time peaks into the member value smears the edge (sigma 84.5-277.7 MeV against 67.4-94.2 MeV for `w1`; eight of the nine modules are at 113 MeV or more), the fit error rises by 1.0-7.4 times (W564 0.071, W552 0.023; W590 unchanged) and s moves by +0.7 to +5.3 %. It is the "all events, no cut" literal reading and bounds the effect from above.
- Among the settings that keep the in-time integral (`none_main`, `w1`, `loose1`, `loose2`) s moves by +0.04 to +0.39 % relative to wave 1, always upward and in the same direction for all nine modules (W564 -0.13 % for `loose1`, inside its fit error). The spread of s across these four settings is 0.05-0.46 % per module (median 0.26 %; the 0.46 % is W564, where `loose1` is -0.13 %), about equal to the wave-1 `s_err` of 0.28-0.88 %, and the shift is much smaller than the window-spread term of `s_err`. So the all-event edge differs from the single-pulse edge by +0.1 to +0.4 %; pile-up does not bias the edge at the 1 % level once secondaries are not summed in.
- The shift is nearly common to the modules: about +0.2 to +0.3 % (mean of the four pool modules: +0.31 `none_main`, +0.21 `loose1`, +0.26 `loose2`; `loose2` per module +0.33, +0.26, +0.18, +0.25 %), and a common factor cancels in N_rt. The module-to-module part is the sd of the four pool shifts: 0.06-0.07 % for the three in-time settings and 0.33 % for `none_sum`; all nine modules 0.09-0.12 % (in-time) and 1.4 % (`none_sum`).
- Effect on the gate (`transfer.transfer` re-run with each setting's s, primary L_scan, N_rt refit, same gains): pool spread of r at 2.24 GeV 2.28 % (w1), 2.34 % (none_main), 2.33 % (loose1), 2.34 % (loose2), 2.69 % (none_sum); at 3.49 GeV 3.03, 3.15, 3.14, 3.15, 3.60 %; pool+H 5.60, 5.66, 5.38, 5.53, 6.18 % (2.24 GeV). The 2.29 % spread (excess over the expected 0.87 % is 2.1 %) moves by +0.05 pp for the in-time settings (+0.4 pp for `none_sum`, upward). The pile-up systematic therefore explains less than 3 % of the excess variance (0.07^2 / 2.1^2 for the module-to-module part; 0.33 % in the extreme `none_sum` case: 2.5 %) and none of the gate result: the failure of the transfer (pool spread 2.29 % against 2 %) is not caused by the pile-up rejection.

## Recommendation for the farm pass

Keep the wave-1 rejection (`edge_pass.collect_events` as is: reject when any fired member of the 3x3 has `len(peaks) >= 2`; member value = the single peak's integral; no parameter changes). Criterion: smallest fit error and window spread per module, and smallest spread of r over the 4 pool modules. `w1` has the smallest s_err of the four in-time settings in all nine modules and the smallest 4-module pool spread of r at both energies (2.28 % and 3.03 %); the looser settings differ by +0.05-0.12 pp, inside the noise, and keep 1.3-2.0 times as many events. The pool+H spread does not favour `w1`: at 2.24 GeV `loose1` is lowest (5.38 % against 5.60 % for `w1`), at 3.49 GeV `w1` (2.998 %) is the largest of the four in-time settings (2.985-2.998 %, a difference of 0.013 pp). The recommendation rests on the 4-pool-module spread and the smallest s_err, not on the pool+H spread. If statistics are needed in the farm (modules with few events), `loose2` (secondary >= 200 ADC within 60 ns, any member, in-time integral) is an equivalent alternative: its s is +0.04-0.33 % above wave 1; a farm pass that mixes the two between modules would add a 0.1 % module-dependent shift, so use one setting for all. Add 0.3 % (relative, common to all modules) to the s systematic if desired; it cancels in N_rt. Do not use `none_sum`.

## Caveats

- One run (23701), nine modules, one dwell per module; the fit is the same linear-plateau x erfc model, so a model systematic common to all settings (the window spread) is not tested here.
- Only 3 peaks per member are stored. The in-time peak is defined by time (T0 = 158 ns); secondaries are identified by the WaveAnalyzer peak finder, which can merge close pulses, so a pile-up pulse closer than the finder's `peak_pileup_gap` to the main pulse is invisible to every setting here (an irreducible high edge tail, part of the s systematic that this check cannot bound).
- `loose1` and `loose2` thresholds (50 and 200 ADC; 100 and 60 ns) are chosen to span 4-17 % rejection, not tuned; the central amplitude at the edge is about 1500-1600 ADC, so 200 ADC is 13 % of an edge pulse.
- W594's 3x3 has 7 members (containment differs); its shifts are the smallest (0.04-0.74 %) because 30 % of its events were rejected in wave 1.
- Event caches (400 MB, npz) are in the scratchpad `sup-edgecache/full`, not in the repository; rerun `edge_pileup_check.py collect` (about 300 CPU s per file) to regenerate.
