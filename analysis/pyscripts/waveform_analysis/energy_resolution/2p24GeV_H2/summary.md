# 2p24GeV_H2: PbWO4 energy resolution, legacy versus time-gated clustering

Runs 25228, 25257, 25260, 25299, 25320; mean ± standard error over 42 splits; default arm = gains and position correction.  The scan over the gates 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5, 21, 21.5, 22, 22.5, 23, 23.5, 24 ns is in summary.tsv; columns are defined in README.md.

## Default arm, legacy and gate 12 ns (differences gated - legacy in percentage points)

| path | Crystal Ball sigma/E (%) | out-of-sample (even/odd cross fit), Crystal Ball sigma/E (%) | Gaussian within 2 sigma (%) | interquartile range / 1.349 (%) | Crystal Ball peak / E' - 1 (%) | difference, Crystal Ball | difference, Gaussian | yield (%) | pile-up difference | pile-up flag (%) | clipped seed (%) | two-peak fits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| legacy | 1.876 ± 0.005 | 1.889 ± 0.005 | 1.935 ± 0.004 | 8.670 ± 0.105 | 0.065 ± 0.005 | - | - | 100.00 | -0.043 ± 0.002 | 24.290 ± 0.052 | 0.007 ± 0.001 | 42/42 |
| gate 12 ns | 1.864 ± 0.005 | 1.874 ± 0.004 | 1.922 ± 0.003 | 8.455 ± 0.103 | 0.120 ± 0.005 | -0.012 ± 0.002 | -0.013 ± 0.001 | 100.35 | -0.030 ± 0.002 | 24.221 ± 0.051 | 0.007 ± 0.001 | 42/42 |

## Correction arms

Every arm after "none" has the gains and the position correction; Moller rejection and non-linearity are added as named.  Crystal Ball sigma/E (%) unless stated.

| path | none | gains + position | + Moller rejection (switch off) | + non-linearity (switch off) | all on (switch off) |
|---|---|---|---|---|---|
| legacy | 2.004 ± 0.006 | 1.876 ± 0.005 | 1.876 ± 0.005 | 1.876 ± 0.005 | 1.876 ± 0.005 |
| gate 12 ns | 1.989 ± 0.005 | 1.864 ± 0.005 | 1.864 ± 0.005 | 1.864 ± 0.005 | 1.864 ± 0.005 |
| difference, gated - legacy (percentage points) | -0.015 ± 0.002 | -0.012 ± 0.002 | -0.012 ± 0.002 | -0.012 ± 0.002 | -0.012 ± 0.002 |
| legacy, Crystal Ball peak / E' - 1 (%) | 2.126 ± 0.007 | 0.065 ± 0.005 | 0.065 ± 0.005 | 0.065 ± 0.005 | 0.065 ± 0.005 |
