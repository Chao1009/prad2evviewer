# 2p24GeV_Ta_X17_moller: Moller two-cluster calibration of the Ta X17 runs

Runs 25800, 26138, 26230 (2.24 GeV, Ta target, two-cluster trigger), PbWO4 only.  Produced by

    python3 waveform_analysis/energy_resolution_benchmark.py report --calib moller output/energy_resolution/scan/2p24GeV_Ta_X17

(`waveform_analysis/moller_calib.py`; `moller_calib.py test` runs the kinematics unit test).  Without `--calib` (or with `--calib ep`)
the benchmark is the control-v3 ep chain, unchanged.  Files: `summary.tsv` (all results, columns in its comment line), `cutflow.tsv`,
`geometry.tsv` (beam centre and vertex distance per run), `gains.tsv` (Moller gains), `compare.tsv` (third-cluster rules), `runs.tsv`
(per run), `report.pdf` (selection distributions, gain convergence, Crystal Ball fits, resolution versus energy, 20 example events,
per-split values, the isolation radius, the rule comparison and the per-run table).  The previous default, exactly two clusters above
300 MeV, is `--cl-min 300 --r-iso 0`.

## Existing X17 measurement (control v3, for comparison)

Single-cluster events (exactly one cluster above 500 MeV), seed in the ep pool (radius < 200 mm, outside the central 6 x 6 block),
target E' = elastic e-Ta at the seed module's polar angle.  Per-module gains per run fitted to E' on the legacy path, the 2.24 GeV H2
position correction borrowed, evaluated on the 4 ns gate path with gains from the other half of the splits: 1.987 +- 0.025 % at
E' = 2.237 GeV.  The H2 curve 0.90 (+) 2.05 / sqrt(E) (+) 1.35 / E gives 1.75 % there.

## Method

- Clustering as the benchmark (same hit preprocessing, legacy and 4 ns gate paths).  Candidates: the two largest clusters of an
  event, both PbWO4-seeded.
- Cuts, in order: isolation, no other cluster above 50 MeV (the lowest reconstructed cluster energy) within 60 mm of either pair
  cluster (next section); both seeds outside the central-hole block H and the outer two PbWO4 layers T (the benchmark's region
  definitions), seed energy >= 200 MeV; energy sum, calibration-file energies, within 20 % of Ebeam; coplanarity
  | |phi1 - phi2| - 180 | < 10 deg about the beam centre (core within ~5 deg, flat tail beyond 10 deg); Moller angle relation
  |tan th1 tan th2 (gamma + 1) / 2 - 1| < 0.15 (core sigma 0.056, about 2.7 sigma; flat background outside).  No global
  third-cluster threshold (`--cl-min` adds one).
- Geometry per run from the selected pairs: beam centre = point nearest to the pair lines (x0 = 0.78-0.87 mm, y0 = 1.42-1.84 mm), vertex
  distance z = median of sqrt((gamma + 1) / 2 R1 R2) (7583-7595 mm against hycal_z 7531 mm: shower depth).
- E_pred(theta) = Moller E'(theta) at 2239.5 MeV minus the target and window energy loss (`e_moller`, as `PhysicsTools`).  Unit test:
  E = (Ebeam + m_e) / 2 at the symmetric angle 1.2236 deg, E1 + E2 = Ebeam + m_e and momentum balance at 7 angles.
- Gains: one per member module and run (as the ep chain), both clusters of each pair, sum g_i E_i fitted to E_pred with the
  benchmark's prior, windows, refit passes and convergence test (`solve_gains` with per-cluster targets).  One departure, forced by
  the statistics (about 13 clusters per module and run): a module without its own window uses the global window instead of being
  dropped.  In sample: all splits of the run; out of sample: even / odd split ordinal within the run, gains from the other half.  All
  nine solves converge in 3-5 passes (in sample 149 / 93 / 178 free modules, 114 / 36 / 172 with their own window for runs 25800 /
  26138 / 26230).  One solve over all three runs (`--combine-runs`, previous selection) is worse: sum 2.16 +- 0.05 % per split against
  2.11, every cluster 4.24 against 4.06 %; the three runs use different calibration files (x17_25825, x17_26094, x17_26214.json).
- Resolution: Crystal Ball as the benchmark's `cb_single`, per split (mean +- standard error over the 25 splits, the control-v3
  protocol) and pooled (+- bootstrap over the splits).

## Cut flow (4 ns gate path; legacy in cutflow.tsv)

| cut | events |
|---|---|
| all events | 2 520 012 |
| hit sum > 0.8 Ebeam (prefilter) | 1 340 288 |
| two largest clusters PbWO4, loose (dphi < 30 deg, sum within 30 %) | 176 723 |
| isolation, 60 mm | 83 815 |
| fiducial (outside H and T, seed >= 200 MeV) | 38 744 |
| energy sum within 20 % | 36 223 |
| coplanar, 10 deg | 30 378 |
| Moller angle relation, 0.15 | 20 491 |

## Results (isolation 60 mm, 4 ns gate, Moller gains, out of sample)

| quantity | E (GeV) | per split, mean +- s.e. (%) | pooled (%) | H2 curve (%) |
|---|---|---|---|---|
| (ii) E1 + E2 (angle independent) | 2.236 | 2.09 +- 0.04 | 2.16 +- 0.04 | 1.75 (two clusters in quadrature) |
| symmetric pairs, sqrt(2) sigma_sum / E_sum = one cluster | 1.120 | 2.85 +- 0.11 | 3.03 +- 0.12 | 2.45 |
| (i) E / E_pred, every cluster | 1.119 | 4.10 +- 0.05 | 4.11 +- 0.05 | 2.45 |
| (iii) E / E_pred, E_pred 0.32-0.77 GeV | 0.612 | 5.99 +- 0.12 | 6.42 +- 0.09 | 3.54 |
| (iii) E_pred 0.77-1.12 GeV | 0.942 | 4.86 +- 0.07 | 4.98 +- 0.08 | 2.71 |
| (iii) E_pred 1.12-1.47 GeV | 1.295 | 3.90 +- 0.08 | 4.01 +- 0.08 | 2.27 |
| (iii) E_pred 1.47-1.90 GeV | 1.625 | 3.03 +- 0.07 | 3.06 +- 0.06 | 2.02 |

In sample the sum is 2.03 +- 0.02 %.

- The sum (ii) is the resolution measurement: it does not depend on the angles.  It is 19 % wider than the H2 curve predicts for the
  same two clusters.  The ep X17 value is 14 % above the curve at 2.24 GeV; the symmetric pairs, 2.85 +- 0.11 % at 1.12 GeV, are 16 %
  above it (pooled 3.03 %, 23 %).  The Moller pairs give the same picture as ep: the X17 data sit 15-25 % above the H2 curve.
- (i) and (iii) are not resolution measurements.  E_pred(theta) inherits the angle measurement: the Moller angle relation has a core
  width of 5.6 %, comparable to the measured widths, and the peaks move from -3.0 % in the low bin to 0 % in the high bin.

## Isolation radius

Report page 10: the distance from the pair to the nearest other cluster above 50 MeV (every other cut applied, no third-cluster cut;
25 812 pairs have one, 1 044 have none) and the raw pair sum versus that distance.  Below 50 mm the Crystal Ball peak of E1 + E2
sits 1.5-4.5 % below the plateau (split-off shower or radiated photon); from 50 mm on it is within 0.3 % of the plateau of the isolated pairs.  The
median keeps recovering until about 120 mm (the low tail, which the Crystal Ball tail absorbs).  60 mm is the default; 100 mm gives the
same width (table below) with 37 % fewer pairs.

## Third-cluster rules compared (E1 + E2, 4 ns gate, out of sample)

| rule | pairs | own windows (3 runs) | sum per split (splits fitted) | sum pooled | one cluster at 1.12 GeV: per split; pooled |
|---|---|---|---|---|---|
| isolation 60 mm (default) | 20 491 | 322 | 2.09 +- 0.04 (25) | 2.16 +- 0.04 | 2.85 +- 0.11; 3.03 +- 0.12 |
| isolation 100 mm | 12 924 | 175 | 2.09 +- 0.04 (25) | 2.13 +- 0.04 | 2.75 +- 0.14; 2.94 +- 0.13 |
| exactly two above 300 MeV (previous default) | 19 798 | 249 | 2.11 +- 0.05 (25) | 2.15 +- 0.06 | 2.75 +- 0.13; 3.00 +- 0.15 |
| exactly two above 100 MeV | 5 131 | 8 | 2.74 +- 0.12 (25) | 2.92 +- 0.12 | (2 splits fit); 4.04 +- 0.31 |
| exactly two above 50 MeV | 1 042 | 3 | 2.07 +- 0.22 (6) | 2.50 +- 0.15 | (none fit); 3.13 +- 0.87 |

The isolation cut and the 300 MeV threshold agree within the errors; the isolation cut removes the split-offs that a global threshold
keeps (100 MeV) or keeps the statistics a strict threshold loses (50 MeV).

## Per run (isolation 60 mm, 4 ns gate, out of sample)

| run | start | calibration file | pairs | splits | sum per split | sum pooled | inner layers 1-2: gain / calibration-file part | layer 3 |
|---|---|---|---|---|---|---|---|---|
| 25800 | 2026-07-12 | x17_25825 | 7 525 | 8 | 2.10 +- 0.04 | 2.11 +- 0.04 | 1 / 1 | 1 / 1 |
| 26138 | 2026-08-06 | x17_26094 | 3 742 | 4 | 1.96 +- 0.09 | 1.98 +- 0.10 | 0.880 / 0.884 | 0.885 / 0.901 |
| 26230 | 2026-08-12 | x17_26214 | 9 224 | 13 | 2.13 +- 0.06 | 2.17 +- 0.06 | 1.975 / 1.989 | 0.861 / 0.894 |

Gain = mean over the layer's modules of 1 / (calibration factor x Moller gain of the run's in-sample solve), relative to run 25800;
layers are rings around the beam hole (1-2 = the central block H, 32 modules; 3 = the next ring, 28 modules).  The layer 1-2 modules
are never Moller seeds (fiducial cut), so their value is essentially the calibration file.  Run dates: run start from the EVIO
unix time (notes/PbWO4_energy_resolution_20261028.md) and the elog summary (notes/elog_pradlog_summary_20261117.md).  Elog
context: gain files followed the radiation damage (x17_v4 on 08-08, inner modules 50-70 % damage); the innermost two layers' HV was
raised on 08-10 (2-3x gain) and lowered again on 08-14.  The calibration files used for these runs give: in 26138 (08-06) the inner gain 12 % below 25800 in layers 1-2 and in layer 3; this
file predates the x17_v4 damage file of 08-08, so the 12 % does not measure the 50-70 % loss.  In 26230 (08-12, after the HV raise)
layers 1-2 are 2.0x 25800 and layer 3 is 14 % below.  The pair-sum resolution does not follow these steps: 2.10, 1.96 and 2.13 % agree within 1.6 standard errors,
and the Moller fiducial region excludes layers 1-2.
