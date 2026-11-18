# energy_resolution

Output of `waveform_analysis/energy_resolution_benchmark.py report`: PbWO4 energy resolution of legacy and time-gated
island clustering, one directory per run type: `0p73GeV_H2`, `2p24GeV_H2`, `3p49GeV_H2`, `3p49GeV_H2_old` (the previous 3.49 GeV
runs 24386, 24503, 24512; beam energy 3488.43 MeV from `database/runinfo/general.json`, the new runs 3485.41 MeV) and
`2p24GeV_Ta_X17`.  Run `2p24GeV_H2` before `2p24GeV_Ta_X17`: the Ta position correction reads the 2.24 GeV H2 polynomials and map.

## Chain

1. `scan` (one run): physics events (trigger bits == 0x100; X17: trigger bits & 0xf00 without the light-monitoring bit, as
   `template_extract.py`); every PbWO4 channel fitted by the matched filter of
   `pileup_analysis/pileup_deconvolution_study.py` with the run type's template bank; hit time = onset + t50, energy =
   calibration x amplitude x walk-out integral x I(phi) / mean(I), with I(phi) the walk-out integral of the module's template at the
   onset phase phi = t0 mod 4 ns (256 phases, linear interpolation; 1.8 % peak to peak per module).  Lead glass and clipped
   PbWO4 channels without a usable pulse use WaveAnalyzer.  Only hits with 100 < t < 200 ns are stored, each with module,
   energy, amplitude, the number of pulses in its channel and an extra flag (not the module's largest in-window pulse).
   Calibration file, beam energy and dead modules come from `database/runinfo/general.json`.  `3p49GeV_H2_old` uses the
   `3p49GeV_H2` bank: the bank of its own runs (`pulse_templates/3p49GeV_H2_old`) has no `a_band`, which the matched filter needs.
2. `report` (per run type and split): HyCalCluster with the builder configuration, min_center_energy 30 MeV, the Geant4 profile
   and the builder's non-linearity switch on.  Legacy path: seed_time_window -1, one hit per module (the largest in-window
   pulse).  Gated path, one per gate: seed_time_window = gate, every in-window pulse a hit.
3. Selection as `benchmark_hycal_timing_test_20261021.py` (elastic-window method, seed-indexed refinement off): events with
   exactly one cluster above 500 kE; pool = seed radius < 200 mm, seed energy >= 200 kE, seed outside rows 15-20 x columns
   15-20; used modules = at least 30 events within 2 sigma of their Crystal Ball prefit and an elastic-Moller gap >= 3 sigma.
4. Fits of the used-module spectrum on E' +- 800 kE (20 kS bins): two-peak Crystal Ball (headline sigma/E), Gaussian within
   2 sigma (iterated), interquartile range / 1.349; peak position / E' - 1; peak yield within +- 2.5 s_ref of E' on the modules used
   by both paths.  The same fits without the clusters that have a member module with at least 2 in-window pulses give the
   pile-up flag systematic.
5. Arms.  The corrections are derived from the legacy pool (step 3 cuts, all splits) and applied identically to every path.  An
   arm is a set of switches; the default arm (all on) is fitted for legacy and every gate, the other arms for legacy and the
   gate `ARM_GATE` = 12 ns:

   | arm | gains + position | Møller rejection | non-linearity |
   |---|---|---|---|
   | none | no | no | no |
   | gains + position | yes | no | no |
   | + Møller rejection | yes | yes | no |
   | + non-linearity | yes | no | yes |
   | all on (default) | yes | yes | yes |

   `MOLLER_REJECT` and `NONLINEAR` (both False by default, pending confirmation, 2026-11-12) switch these two effects on in every arm that has them.

   - Per-crystal gains (`percrystal_calib.py`): one gain per member module and run, solved on all splits of the run: the
     sum of g_i E_i of the events inside their seed module's elastic window is matched to E', with a prior pulling every gain to
     g0 = mean E' / mu_global.  A cluster's energy is sum g_i E_i; modules never in a selected event keep g0.
   - Position correction (`position_corr.py`, notes/position_correction_20261105.md): u, v = (cluster centroid - seed module
     centre) / module size; corrected E = E / R(u, v), R fitted to r = E / E' of the pool (r in 0.85-1.12, mean 1) with the gains.
     Setting per run type (`PC_SETTING`, the best out-of-sample variant of `poscorr_comparison.tsv`): hybrid (own eight-coefficient
     polynomial per seed module with >= 1000 clusters, the 12 x 12 bilinear map elsewhere) at 3.49 and 2.24 GeV H2, no position
     correction (gains only) at 0.73 GeV, and the full 2.24 GeV H2 hybrid (own polynomials and map) for Ta.  Map bins are fitted with the Crystal Ball at 3.49 GeV and the Gaussian within 2 sigma elsewhere.
   - Møller rejection = the old chain's `--moller-reject kin` (`benchmark_hycal_timing_test_20261021.py:3297-3306`, N = 3, M = 2):
     per cluster, at its polar angle theta (cluster centroid; in the gain stage the seed module centre, as the old dump), reject if |E - E_ee(theta)| < 3 sigma and |E - E'(theta)| > 2 sigma, with
     sigma = E' sqrt((0.025 / sqrt(E' / GeV))^2 + 0.007^2) and E the cluster energy after gains and position correction.
     It removes the clusters from the event selection and, as the old chain's gains stage (`percrystal_calib.py:560-636`), from
     the gain fit: gains of the full pool, tag with them (no position correction), gains refitted on the rest, position correction
     fitted on that pool with the refitted gains.
   - Non-linearity = the old chain's `--nl-members on` (`benchmark_hycal_timing_test_20261021.py:2791-2792`): the clusterer's
     factor 1 / (1 + nl1 d + nl2 d^2) of the seed module (calibration file, `prad2det/src/HyCalCluster.cpp:457-473`; d = (E - base
     energy) / GeV; 1 above 3800 MeV or outside 0.7-1.3) multiplies the cluster energy and its member energies, so the gains carry
     it.  The clusterer applies it by default; the arms without it use the sum of the raw member energies.  Files without
     nl1, nl2 (all except `calibration_factor_3p5_June7_nonLinear.json`, that is 3.49 and 0.73 GeV) are unaffected.

kE = E_beam / 2239.51, kS = sqrt(kE).  All settings are named constants at the top of the script, grouped by step.  The gains and
the position correction are fitted on the splits they are applied to (in-sample, as in the old chain); with `OUT_OF_SAMPLE` (default True) the
default arm for legacy and every gate is also evaluated with corrections fitted on the other half of the splits (columns `*_oos`).

## Outputs (`<run type>/`)

- `splits.tsv`: one row per run, split, path and variant (`all`, `nopu` = pile-up-flagged clusters removed) of the default arm:
  counts, fit parameters, `res_cb` (Crystal Ball), `res_g2` (Gaussian within 2 sigma), `res_iqr` (interquartile range) in %,
  `tail` (% outside 3 sigma of the Crystal Ball), `mu_ep` (peak / E' - 1, %), `pu` and `clip` (% of the used pool), peak yield
  (`yield_n`, `yield_legacy_n` on `n_common` modules).  `n_single` counts events with one cluster, `n_pool` the pool after the
  Møller rejection.  The other arms add `res_cb_<arm>`, `res_g2_<arm>`, `mu_ep_<arm>` (arms `none`, `gp` gains + position, `gpm`
  + Møller rejection, `gpn` + non-linearity) for legacy and gate 12 ns, empty elsewhere.
  `res_cb_oos`, `res_g2_oos`, `mu_ep_oos`: the default arm with gains and position correction fitted on the other half of the splits
  (even/odd cross fit; half = ordinal 0, 2, 4, ... or 1, 3, 5, ... of the split within its run, as split numbers are almost all even),
  for legacy and every gate, empty for a run's lone split; X17 uses the same borrowed map in both halves.
- `modules.tsv`: the used modules per run, split, path and variant.
- `summary.tsv`: per path, mean and standard error over splits (`_se`) of the default arm for legacy and every gate; `d_cb`, `d_g2`
  = gated - legacy per split (percentage points, paired); `d_pu` = Crystal Ball without the flagged clusters minus all; the arm
  columns as in `splits.tsv` (including `*_oos`) plus `d_cb_<arm>`.
- `summary.md`: the default arm for legacy and gate 12 ns (with an out-of-sample Crystal Ball column), and one table with the five arms.
- `report.pdf`: page 1 compares legacy and every gate (resolution of the default arm, the other arms as markers at legacy and
  gate 12 ns, peak position, peak yield, paired difference per split, pile-up flag systematic, table); page 2 the spectra and
  the resolution per split.
- `spectra.npz`: E / E' histograms of the used-module clusters (all splits, bins `edges`), one array per path.
- `runs.tsv` (written by `summary`): per run, for legacy and the chosen gate, mean and standard error over the run's splits of the
  Crystal Ball sigma/E and peak / E' - 1 (in-sample and out-of-sample), used clusters, and a flag for runs more than 3 standard errors from the run-type mean.
- `gains.tsv`: per switch combination (`moller`, `nonlinear`), run and module: gain, its uncertainty, seeds and members in the
  selected events, prior fraction.
- `poscorr_map.tsv` (bilinear map per switch combination, `filled` = taken from neighbours) and `poscorr_modules.tsv` (own
  polynomials c0..c7 of R = c0 (1 + c1 u^2 + c2 v^2 + c3 u^2 v^2 + c4 u^4 + c5 v^4 + c6 u + c7 v)).

## Scanning the gate

The scans hold no clustering, so a new gate list needs only the report (minutes; `--out <dir>` writes elsewhere):

```bash
python3 waveform_analysis/energy_resolution_benchmark.py report output/energy_resolution/scan/2p24GeV_H2 4 8 12 16 24 32
```

`run_energy_resolution.sh` sets the default list (`GATES`: 45 values, 2-24 ns in 0.5 ns steps).  A scan may be split over several
files of one run (`run_<run>_a.npz`, ...): the report pools them per run.

## Summary

`python3 waveform_analysis/energy_resolution_benchmark.py summary waveform_analysis/energy_resolution` writes `summary.pdf` (and each
run type's `runs.tsv`) from the five reports and `poscorr_comparison.tsv`: overview table of the best configuration per run type
(best gate = the smallest gate whose out-of-sample Crystal Ball sigma/E is within one standard error of the minimum over the gates
with peak yield >= 99 % of legacy), resolution progression (none, gains, gains + position correction, best gate), gate scans,
resolution per run, position-correction variants, and one E / E' distribution page per run type.  `run_energy_resolution.sh` runs it last.
