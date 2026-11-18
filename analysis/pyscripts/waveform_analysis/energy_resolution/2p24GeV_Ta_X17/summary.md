# 2p24GeV_Ta_X17: PbWO4 energy resolution, legacy versus time-gated clustering

Runs 25800, 26138, 26230; mean ± standard error over 25 splits; default arm = gains and position correction.  The scan over the gates 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5, 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5, 20, 20.5, 21, 21.5, 22, 22.5, 23, 23.5, 24 ns is in summary.tsv; columns are defined in README.md.

## Default arm, legacy and gate 12 ns (differences gated - legacy in percentage points)

| path | Crystal Ball sigma/E (%) | out-of-sample (even/odd cross fit), Crystal Ball sigma/E (%) | Gaussian within 2 sigma (%) | interquartile range / 1.349 (%) | Crystal Ball peak / E' - 1 (%) | difference, Crystal Ball | difference, Gaussian | yield (%) | pile-up difference | pile-up flag (%) | clipped seed (%) | two-peak fits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| legacy | 2.308 ± 0.035 | 2.332 ± 0.036 | 2.285 ± 0.029 | 20.485 ± 0.464 | 0.009 ± 0.028 | - | - | 100.00 | -0.270 ± 0.058 | 36.449 ± 1.151 | 0.000 ± 0.000 | 13/25 |
| gate 12 ns | 2.036 ± 0.016 | 2.050 ± 0.019 | 2.090 ± 0.013 | 20.431 ± 0.463 | -0.184 ± 0.035 | -0.272 ± 0.027 | -0.195 ± 0.023 | 103.36 | -0.107 ± 0.053 | 34.959 ± 1.149 | 0.000 ± 0.000 | 19/25 |

## Correction arms

Every arm after "none" has the gains and the position correction; Moller rejection and non-linearity are added as named.  Crystal Ball sigma/E (%) unless stated.

| path | none | gains + position | + Moller rejection (switch off) | + non-linearity (switch off) | all on (switch off) |
|---|---|---|---|---|---|
| legacy | 2.597 ± 0.040 | 2.308 ± 0.035 | 2.308 ± 0.035 | 2.308 ± 0.035 | 2.308 ± 0.035 |
| gate 12 ns | 2.304 ± 0.029 | 2.036 ± 0.016 | 2.036 ± 0.016 | 2.036 ± 0.016 | 2.036 ± 0.016 |
| difference, gated - legacy (percentage points) | -0.293 ± 0.024 | -0.272 ± 0.027 | -0.272 ± 0.027 | -0.272 ± 0.027 | -0.272 ± 0.027 |
| legacy, Crystal Ball peak / E' - 1 (%) | 1.678 ± 0.327 | 0.009 ± 0.028 | 0.009 ± 0.028 | 0.009 ± 0.028 | 0.009 ± 0.028 |
