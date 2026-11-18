# 3p49GeV_H2: PbWO4 energy resolution, legacy versus time-gated clustering

Runs 24698, 24911, 24937; mean ± standard error over 20 splits; default arm = gains and position correction.  The scan over the gates 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5, 21, 21.5, 22, 22.5, 23, 23.5, 24 ns is in summary.tsv; columns are defined in README.md.

## Default arm, legacy and gate 12 ns (differences gated - legacy in percentage points)

| path | Crystal Ball sigma/E (%) | out-of-sample (even/odd cross fit), Crystal Ball sigma/E (%) | Gaussian within 2 sigma (%) | interquartile range / 1.349 (%) | Crystal Ball peak / E' - 1 (%) | difference, Crystal Ball | difference, Gaussian | yield (%) | pile-up difference | pile-up flag (%) | clipped seed (%) | two-peak fits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| legacy | 1.585 ± 0.006 | 1.614 ± 0.006 | 1.637 ± 0.006 | 2.029 ± 0.006 | 0.147 ± 0.008 | - | - | 100.00 | -0.054 ± 0.005 | 33.742 ± 0.223 | 0.003 ± 0.001 | 0/20 |
| gate 12 ns | 1.571 ± 0.006 | 1.602 ± 0.007 | 1.614 ± 0.006 | 2.007 ± 0.006 | 0.205 ± 0.008 | -0.014 ± 0.002 | -0.023 ± 0.004 | 100.33 | -0.041 ± 0.005 | 33.705 ± 0.230 | 0.003 ± 0.001 | 0/20 |

## Correction arms

Every arm after "none" has the gains and the position correction; Moller rejection and non-linearity are added as named.  Crystal Ball sigma/E (%) unless stated.

| path | none | gains + position | + Moller rejection (switch off) | + non-linearity (switch off) | all on (switch off) |
|---|---|---|---|---|---|
| legacy | 2.184 ± 0.032 | 1.585 ± 0.006 | 1.585 ± 0.006 | 1.585 ± 0.006 | 1.585 ± 0.006 |
| gate 12 ns | 2.177 ± 0.033 | 1.571 ± 0.006 | 1.571 ± 0.006 | 1.571 ± 0.006 | 1.571 ± 0.006 |
| difference, gated - legacy (percentage points) | -0.007 ± 0.003 | -0.014 ± 0.002 | -0.014 ± 0.002 | -0.014 ± 0.002 | -0.014 ± 0.002 |
| legacy, Crystal Ball peak / E' - 1 (%) | 4.184 ± 0.277 | 0.147 ± 0.008 | 0.147 ± 0.008 | 0.147 ± 0.008 | 0.147 ± 0.008 |
