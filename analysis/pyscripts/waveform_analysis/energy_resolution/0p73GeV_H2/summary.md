# 0p73GeV_H2: PbWO4 energy resolution, legacy versus time-gated clustering

Runs 24590, 24603; mean ± standard error over 16 splits; default arm = gains and position correction.  The scan over the gates 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5, 21, 21.5, 22, 22.5, 23, 23.5, 24 ns is in summary.tsv; columns are defined in README.md.

## Default arm, legacy and gate 12 ns (differences gated - legacy in percentage points)

| path | Crystal Ball sigma/E (%) | out-of-sample (even/odd cross fit), Crystal Ball sigma/E (%) | Gaussian within 2 sigma (%) | interquartile range / 1.349 (%) | Crystal Ball peak / E' - 1 (%) | difference, Crystal Ball | difference, Gaussian | yield (%) | pile-up difference | pile-up flag (%) | clipped seed (%) | two-peak fits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| legacy | 3.389 ± 0.014 | 3.420 ± 0.014 | 3.884 ± 0.035 | 14.389 ± 0.114 | 0.226 ± 0.015 | - | - | 100.00 | -0.047 ± 0.008 | 9.784 ± 0.146 | 0.000 ± 0.000 | 16/16 |
| gate 12 ns | 3.370 ± 0.013 | 3.401 ± 0.013 | 3.827 ± 0.029 | 14.379 ± 0.115 | 0.279 ± 0.013 | -0.019 ± 0.004 | -0.057 ± 0.022 | 100.23 | -0.027 ± 0.007 | 9.754 ± 0.147 | 0.000 ± 0.000 | 16/16 |

## Correction arms

Every arm after "none" has the gains and the position correction; Moller rejection and non-linearity are added as named.  Crystal Ball sigma/E (%) unless stated.

| path | none | gains + position | + Moller rejection (switch off) | + non-linearity (switch off) | all on (switch off) |
|---|---|---|---|---|---|
| legacy | 4.103 ± 0.077 | 3.389 ± 0.014 | 3.389 ± 0.014 | 3.389 ± 0.014 | 3.389 ± 0.014 |
| gate 12 ns | 4.112 ± 0.078 | 3.370 ± 0.013 | 3.370 ± 0.013 | 3.370 ± 0.013 | 3.370 ± 0.013 |
| difference, gated - legacy (percentage points) | 0.009 ± 0.074 | -0.019 ± 0.004 | -0.019 ± 0.004 | -0.019 ± 0.004 | -0.019 ± 0.004 |
| legacy, Crystal Ball peak / E' - 1 (%) | -1.734 ± 0.086 | 0.226 ± 0.015 | 0.226 ± 0.015 | 0.226 ± 0.015 | 0.226 ± 0.015 |
