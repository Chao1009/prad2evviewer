# 3p49GeV_H2_old: PbWO4 energy resolution, legacy versus time-gated clustering

Runs 24386, 24503, 24512; mean ± standard error over 28 splits; default arm = gains and position correction.  The scan over the gates 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5, 21, 21.5, 22, 22.5, 23, 23.5, 24 ns is in summary.tsv; columns are defined in README.md.

## Default arm, legacy and gate 12 ns (differences gated - legacy in percentage points)

| path | Crystal Ball sigma/E (%) | out-of-sample (even/odd cross fit), Crystal Ball sigma/E (%) | Gaussian within 2 sigma (%) | interquartile range / 1.349 (%) | Crystal Ball peak / E' - 1 (%) | difference, Crystal Ball | difference, Gaussian | yield (%) | pile-up difference | pile-up flag (%) | clipped seed (%) | two-peak fits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| legacy | 1.559 ± 0.006 | 1.579 ± 0.006 | 1.616 ± 0.005 | 1.990 ± 0.008 | 0.164 ± 0.013 | - | - | 100.00 | -0.035 ± 0.004 | 25.367 ± 0.355 | 0.001 ± 0.001 | 0/28 |
| gate 12 ns | 1.548 ± 0.007 | 1.570 ± 0.006 | 1.600 ± 0.005 | 1.972 ± 0.008 | 0.203 ± 0.013 | -0.011 ± 0.002 | -0.016 ± 0.003 | 100.24 | -0.025 ± 0.004 | 25.322 ± 0.363 | 0.001 ± 0.001 | 0/28 |

## Correction arms

Every arm after "none" has the gains and the position correction; Moller rejection and non-linearity are added as named.  Crystal Ball sigma/E (%) unless stated.

| path | none | gains + position | + Moller rejection (switch off) | + non-linearity (switch off) | all on (switch off) |
|---|---|---|---|---|---|
| legacy | 2.810 ± 0.062 | 1.559 ± 0.006 | 1.559 ± 0.006 | 1.559 ± 0.006 | 1.559 ± 0.006 |
| gate 12 ns | 2.811 ± 0.066 | 1.548 ± 0.007 | 1.548 ± 0.007 | 1.548 ± 0.007 | 1.548 ± 0.007 |
| difference, gated - legacy (percentage points) | 0.002 ± 0.006 | -0.011 ± 0.002 | -0.011 ± 0.002 | -0.011 ± 0.002 | -0.011 ± 0.002 |
| legacy, Crystal Ball peak / E' - 1 (%) | 3.197 ± 0.162 | 0.164 ± 0.013 | 0.164 ± 0.013 | 0.164 ± 0.013 | 0.164 ± 0.013 |
