# P2: region-resolved resolution and separation

Run: `report --set EVAL_REGIONS=True --set MOLLER_PAIRS=True --set MOLLER_BEAM_OFFSET=<mean of the run type's runs in `output/nonlin_20261104/step3/beam_offsets.tsv`, (0,0) where absent>`, gates 8 12 (Ta 8). Values: mean ± SE over per-split fits, out-of-sample columns marked oos. sep_pi = 145 MeV / sigma_E.


## 0p73GeV_H2


### ep, path legacy

| region | status | n_range | modules | prior frac | res_cb (%) | res_cb oos | mu (%) | sigma_E (MeV) | sep_pi | tail_hi3 | leak_cut | eff_cut |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pool | ok | 329391 | 139 | 0.049 | 3.389 ± 0.014 | 3.420 | 0.226 ± 0.015 | 24.6 | 5.86 | 0.047 | 0.0 | 49.3 |
| I1 | ok | 269446 | 60 | 0.029 | 3.418 ± 0.019 | 3.449 | 0.238 ± 0.019 | 24.8 | 5.81 | 0.046 | 0.0 | 47.8 |
| I2 | ok | 59615 | 71 | 0.073 | 3.337 ± 0.040 | 3.397 | -0.032 ± 0.024 | 24.3 | 5.96 | 0.049 | 0.0 | 55.7 |
| I3 | ep-starved | 330 | 8 | 0.094 | 3.070 ± 0.184 | 3.069 | 0.348 ± 0.178 | 22.3 | 6.52 | 0.000 | 0.0 | 88.2 |
| pool_edge | ok | 169239 | 139 | 0.049 | 3.461 ± 0.030 | 3.484 | -0.067 ± 0.029 | 25.2 | 5.74 | 0.051 | 0.0 | 47.3 |
| pool_centre | ok | 160152 | 139 | 0.049 | 3.256 ± 0.023 | 3.326 | 0.524 ± 0.030 | 23.7 | 6.10 | 0.059 | 0.0 | 51.4 |
| H | ok | 1383894 | 32 | 0.807 | 7.131 ± 0.042 | 6.919 | -2.627 ± 0.078 | 51.8 | 2.79 | 0.250 | 2.7 | 60.9 |
| M | ok | 23955 | 348 | 0.984 | 4.149 ± 0.045 | 4.160 | 0.572 ± 0.043 | 30.1 | 4.79 | 1.137 | 1.0 | 86.0 |
| O | ok | 4868 | 252 | 1.000 | 4.081 ± 0.106 | 4.079 | 0.474 ± 0.095 | 29.6 | 4.92 | 0.567 | 0.4 | 87.8 |
| T | ok | 3279 | 253 | 1.000 | 4.821 ± 0.181 | 4.795 | 1.004 ± 0.215 | 35.0 | 4.21 | 0.736 | 1.1 | 75.2 |
| T_lg | not-fitted | 201 | 105 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |
| T_nolg | not-fitted | 3078 | 253 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |

### ep, path g12

| region | status | n_range | modules | prior frac | res_cb (%) | res_cb oos | mu (%) | sigma_E (MeV) | sep_pi | tail_hi3 | leak_cut | eff_cut |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pool | ok | 329779 | 140 | 0.049 | 3.370 ± 0.013 | 3.401 | 0.279 ± 0.013 | 24.5 | 5.89 | 0.041 | 0.0 | 49.4 |
| I1 | ok | 269524 | 60 | 0.030 | 3.405 ± 0.021 | 3.438 | 0.283 ± 0.018 | 24.7 | 5.83 | 0.038 | 0.0 | 47.9 |
| I2 | ok | 59925 | 72 | 0.074 | 3.315 ± 0.042 | 3.377 | 0.026 ± 0.026 | 24.1 | 6.00 | 0.051 | 0.0 | 55.8 |
| I3 | ep-starved | 330 | 8 | 0.094 | 3.084 ± 0.189 | 3.066 | 0.454 ± 0.188 | 22.4 | 6.49 | 0.000 | 0.0 | 88.2 |
| pool_edge | ok | 169906 | 140 | 0.049 | 3.440 ± 0.029 | 3.447 | 0.021 ± 0.025 | 25.0 | 5.78 | 0.041 | 0.0 | 47.5 |
| pool_centre | ok | 159873 | 140 | 0.049 | 3.255 ± 0.020 | 3.323 | 0.535 ± 0.029 | 23.7 | 6.10 | 0.051 | 0.0 | 51.4 |
| H | ok | 1383709 | 32 | 0.807 | 7.114 ± 0.040 | 6.899 | -2.576 ± 0.079 | 51.7 | 2.79 | 0.233 | 2.7 | 60.9 |
| M | ok | 23888 | 348 | 0.984 | 4.126 ± 0.045 | 4.121 | 0.625 ± 0.040 | 30.0 | 4.82 | 1.121 | 1.0 | 86.2 |
| O | ok | 4867 | 252 | 1.000 | 4.039 ± 0.107 | 4.050 | 0.536 ± 0.097 | 29.3 | 4.97 | 0.567 | 0.4 | 88.1 |
| T | ok | 3271 | 253 | 1.000 | 4.751 ± 0.160 | 4.763 | 1.091 ± 0.196 | 34.5 | 4.26 | 0.737 | 1.0 | 75.2 |
| T_lg | not-fitted | 201 | 104 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |
| T_nolg | not-fitted | 3070 | 253 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |

### moller2 ((E1+E2)/E_beam) and moller (E/E_ee), path legacy

| kind | region | status | n | n_range | sigma/E (%) | mu (%) | mu(moller) - mu(ep) (pp) |
|---|---|---|---|---|---|---|---|
| moller2 | all | ok | 61840 | 58936 | g2 4.100 | -3.182 ± 0.016 | nan |
| moller | pool | ok | 256514 | 256514 | g2 4.502 | -0.614 ± 0.011 | -0.840 |
| moller | I1 | ok | 227139 | 227139 | g2 4.432 | -0.592 ± 0.014 | -0.830 |
| moller | I2 | ok | 27157 | 27157 | g2 4.921 | -0.686 ± 0.043 | -0.654 |
| moller | I3 | ok | 2218 | 2218 | g2 9.013 | -0.908 ± 0.211 | -1.256 |
| moller | M | ok | 1337 | 1300 | g2 6.462 | -94.346 ± 88.311 (unreliable fit) | nan |
| moller | O | no-fit | 277 | 227 | g2 nan | nan ± nan | nan |
| moller | T | no-fit | 261 | 213 | g2 nan | nan ± nan | nan |

T clusters (legacy): with lead-glass coincidence 288 (7.6 %), median E/E' 0.8036 (deficit 19.64 %); without 3510, median 1.0006 (deficit -0.06 %).

### drift: mu of the pool, legacy, against split order per run

| run | splits | mu mean (%) | slope (pp / split) ± SE | mu first-last (pp) |
|---|---|---|---|---|
| 24590 | 10 | 0.240 | -0.0013 ± 0.0005 | -0.168 |
| 24603 | 6 | 0.202 | -0.0001 ± 0.0003 | 0.015 |

Continuity: pool legacy res_cb 3.3888 vs baseline pre 3.3888.

## 2p24GeV_H2


### ep, path legacy

| region | status | n_range | modules | prior frac | res_cb (%) | res_cb oos | mu (%) | sigma_E (MeV) | sep_pi | tail_hi3 | leak_cut | eff_cut |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pool | ok | 1075229 | 164 | 0.040 | 1.876 ± 0.005 | 1.889 | 0.065 ± 0.005 | 41.9 | 3.45 | 0.160 | 4.7 | 68.6 |
| I1 | ok | 972928 | 76 | 0.021 | 1.875 ± 0.005 | 1.886 | 0.033 ± 0.005 | 41.9 | 3.45 | 0.160 | 4.5 | 66.7 |
| I2 | ok | 101356 | 75 | 0.068 | 1.877 ± 0.012 | 1.906 | 0.002 ± 0.012 | 42.0 | 3.45 | 0.161 | 5.9 | 86.6 |
| I3 | ep-starved | 945 | 13 | 0.123 | 1.863 ± 0.079 | 1.717 | 0.101 ± 0.082 | 41.6 | 3.52 | 0.536 | 5.5 | 92.1 |
| pool_edge | ok | 385708 | 164 | 0.040 | 1.936 ± 0.009 | 1.956 | 0.080 ± 0.008 | 43.3 | 3.34 | 0.098 | 4.9 | 65.6 |
| pool_centre | ok | 689521 | 164 | 0.040 | 1.844 ± 0.006 | 1.850 | 0.060 ± 0.006 | 41.2 | 3.51 | 0.197 | 4.5 | 70.2 |
| H | ok | 2937778 | 32 | 0.468 | 2.365 ± 0.005 | 2.230 | -0.527 ± 0.011 | 52.9 | 2.73 | 0.081 | 4.4 | 35.6 |
| M | ok | 67147 | 348 | 0.984 | 2.043 ± 0.016 | 2.045 | 0.021 ± 0.016 | 45.6 | 3.17 | 1.232 | 8.6 | 83.2 |
| O | ok | 14664 | 252 | 1.000 | 2.043 ± 0.033 | 2.044 | -0.110 ± 0.037 | 45.6 | 3.20 | 1.458 | 8.6 | 79.5 |
| T | ok | 10283 | 253 | 1.000 | 2.326 ± 0.083 | 2.320 | -0.730 ± 0.090 | 51.8 | 2.95 | 3.070 | 10.3 | 57.3 |
| T_lg | not-fitted | 1774 | 154 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |
| T_nolg | not-fitted | 8509 | 253 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |

### ep, path g12

| region | status | n_range | modules | prior frac | res_cb (%) | res_cb oos | mu (%) | sigma_E (MeV) | sep_pi | tail_hi3 | leak_cut | eff_cut |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pool | ok | 1074544 | 164 | 0.040 | 1.864 ± 0.005 | 1.874 | 0.120 ± 0.005 | 41.7 | 3.47 | 0.118 | 4.5 | 68.9 |
| I1 | ok | 971890 | 76 | 0.021 | 1.861 ± 0.005 | 1.872 | 0.086 ± 0.005 | 41.6 | 3.47 | 0.117 | 4.4 | 66.9 |
| I2 | ok | 101712 | 75 | 0.068 | 1.871 ± 0.011 | 1.891 | 0.078 ± 0.011 | 41.8 | 3.46 | 0.130 | 5.8 | 87.0 |
| I3 | ep-starved | 942 | 13 | 0.124 | 1.837 ± 0.086 | 1.824 | 0.108 ± 0.074 | 41.1 | 3.57 | 0.414 | 5.9 | 92.6 |
| pool_edge | ok | 391092 | 164 | 0.040 | 1.915 ± 0.008 | 1.932 | 0.165 ± 0.008 | 42.8 | 3.38 | 0.071 | 4.7 | 66.0 |
| pool_centre | ok | 683452 | 164 | 0.040 | 1.832 ± 0.006 | 1.838 | 0.101 ± 0.005 | 41.0 | 3.53 | 0.147 | 4.4 | 70.4 |
| H | ok | 2937465 | 32 | 0.468 | 2.344 ± 0.005 | 2.213 | -0.465 ± 0.011 | 52.4 | 2.76 | 0.049 | 4.3 | 35.6 |
| M | ok | 66499 | 348 | 0.984 | 2.019 ± 0.015 | 2.022 | 0.086 ± 0.015 | 45.1 | 3.21 | 0.907 | 8.2 | 83.9 |
| O | ok | 14266 | 252 | 1.000 | 2.004 ± 0.031 | 2.008 | -0.031 ± 0.034 | 44.7 | 3.26 | 0.968 | 7.9 | 81.3 |
| T | ok | 9870 | 253 | 1.000 | 2.281 ± 0.071 | 2.275 | -0.660 ± 0.078 | 50.8 | 2.98 | 2.014 | 9.2 | 58.0 |
| T_lg | not-fitted | 1771 | 154 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |
| T_nolg | not-fitted | 8099 | 253 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |

### moller2 ((E1+E2)/E_beam) and moller (E/E_ee), path legacy

| kind | region | status | n | n_range | sigma/E (%) | mu (%) | mu(moller) - mu(ep) (pp) |
|---|---|---|---|---|---|---|---|
| moller2 | all | ok | 901109 | 885373 | g2 2.155 | -2.543 ± 0.004 | nan |
| moller | pool | ok | 224165 | 224165 | g2 2.879 | 0.478 ± 0.007 | 0.412 |
| moller | I1 | ok | 220994 | 220994 | g2 2.865 | 0.482 ± 0.007 | 0.448 |
| moller | I2 | ok | 2152 | 2152 | g2 7.343 | -0.962 ± 0.300 | -0.964 |
| moller | I3 | no-fit | 1019 | 1019 | g2 nan | nan ± nan | nan |
| moller | M | no-fit | 397 | 361 | g2 nan | nan ± nan | nan |

T clusters (legacy): with lead-glass coincidence 2191 (14.8 %), median E/E' 0.8773 (deficit 12.27 %); without 12650, median 0.9790 (deficit 2.10 %).

### drift: mu of the pool, legacy, against split order per run

| run | splits | mu mean (%) | slope (pp / split) ± SE | mu first-last (pp) |
|---|---|---|---|---|
| 25228 | 11 | 0.068 | 0.0001 ± 0.0000 | 0.067 |
| 25257 | 7 | 0.064 | 0.0000 ± 0.0001 | 0.047 |
| 25260 | 8 | 0.062 | -0.0000 ± 0.0000 | -0.055 |
| 25299 | 8 | 0.064 | -0.0000 ± 0.0001 | 0.006 |
| 25320 | 8 | 0.069 | 0.0002 ± 0.0001 | 0.041 |

Continuity: pool legacy res_cb 1.8759 vs baseline pre 1.8759.

## 3p49GeV_H2


### ep, path legacy

| region | status | n_range | modules | prior frac | res_cb (%) | res_cb oos | mu (%) | sigma_E (MeV) | sep_pi | tail_hi3 | leak_cut | eff_cut |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pool | ok | 352310 | 155 | 0.040 | 1.585 ± 0.006 | 1.614 | 0.147 ± 0.008 | 55.2 | 2.62 | 0.137 | 22.3 | 84.3 |
| I1 | ok | 308862 | 76 | 0.021 | 1.575 ± 0.007 | 1.603 | 0.099 ± 0.010 | 54.8 | 2.64 | 0.144 | 22.1 | 84.4 |
| I2 | ok | 43272 | 74 | 0.068 | 1.658 ± 0.015 | 1.682 | 0.031 ± 0.018 | 57.7 | 2.51 | 0.092 | 22.9 | 83.3 |
| I3 | ep-starved | 176 | 5 | 0.113 | 1.689 ± nan | nan | -0.168 ± nan | 58.7 | 2.46 | 0.000 | 25.7 | 88.6 |
| pool_edge | ok | 125400 | 155 | 0.040 | 1.662 ± 0.012 | 1.692 | 0.211 ± 0.015 | 57.8 | 2.50 | 0.094 | 22.2 | 80.4 |
| pool_centre | ok | 226910 | 155 | 0.040 | 1.540 ± 0.007 | 1.565 | 0.123 ± 0.008 | 53.6 | 2.70 | 0.164 | 22.2 | 86.3 |
| H | ok | 937161 | 32 | 0.572 | 2.411 ± 0.130 | 2.313 | -1.579 ± 0.796 | 84.0 | 2.09 | 2.955 | 16.8 | 41.1 |
| M | ok | 27085 | 348 | 0.988 | 1.983 ± 0.035 | 1.986 | -0.323 ± 0.050 | 68.9 | 2.10 | 0.597 | 25.0 | 76.4 |
| O | ok | 5729 | 252 | 1.000 | 1.951 ± 0.063 | 2.050 | -0.192 ± 0.085 | 67.6 | 2.16 | 0.800 | 23.2 | 70.8 |
| T | ok | 3892 | 253 | 1.000 | 2.378 ± 0.111 | 2.328 | -0.658 ± 0.114 | 82.3 | 1.80 | 1.553 | 21.6 | 55.8 |
| T_lg | not-fitted | 1137 | 178 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |
| T_nolg | not-fitted | 2755 | 240 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |

### ep, path g12

| region | status | n_range | modules | prior frac | res_cb (%) | res_cb oos | mu (%) | sigma_E (MeV) | sep_pi | tail_hi3 | leak_cut | eff_cut |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pool | ok | 352219 | 155 | 0.040 | 1.571 ± 0.006 | 1.602 | 0.205 ± 0.008 | 54.7 | 2.64 | 0.108 | 22.2 | 84.5 |
| I1 | ok | 308791 | 76 | 0.021 | 1.562 ± 0.007 | 1.593 | 0.154 ± 0.010 | 54.4 | 2.66 | 0.110 | 22.0 | 84.6 |
| I2 | ok | 43252 | 74 | 0.068 | 1.637 ± 0.016 | 1.665 | 0.116 ± 0.018 | 57.0 | 2.54 | 0.087 | 22.6 | 83.6 |
| I3 | ep-starved | 176 | 5 | 0.113 | 1.673 ± nan | 1.391 | -0.130 ± nan | 58.2 | 2.48 | 0.000 | 24.3 | 88.6 |
| pool_edge | ok | 126389 | 155 | 0.040 | 1.645 ± 0.010 | 1.678 | 0.271 ± 0.013 | 57.2 | 2.53 | 0.066 | 22.1 | 80.8 |
| pool_centre | ok | 225830 | 155 | 0.040 | 1.527 ± 0.009 | 1.552 | 0.178 ± 0.008 | 53.1 | 2.72 | 0.128 | 22.0 | 86.5 |
| H | ok | 936782 | 32 | 0.572 | 2.514 ± 0.068 | 2.450 | -2.480 ± 1.123 | 87.6 | 1.68 | 5.221 | 19.5 | 43.5 |
| M | ok | 27001 | 348 | 0.988 | 1.959 ± 0.036 | 1.960 | -0.256 ± 0.050 | 68.0 | 2.13 | 0.539 | 24.9 | 76.6 |
| O | ok | 5678 | 252 | 1.000 | 1.996 ± 0.055 | 2.040 | -0.207 ± 0.066 | 69.2 | 2.10 | 0.704 | 23.8 | 72.1 |
| T | ok | 3799 | 253 | 1.000 | 2.344 ± 0.138 | 2.394 | -0.565 ± 0.149 | 81.1 | 1.88 | 1.414 | 21.7 | 56.3 |
| T_lg | not-fitted | 1132 | 178 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |
| T_nolg | not-fitted | 2667 | 237 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |

### moller2 ((E1+E2)/E_beam) and moller (E/E_ee), path legacy

| kind | region | status | n | n_range | sigma/E (%) | mu (%) | mu(moller) - mu(ep) (pp) |
|---|---|---|---|---|---|---|---|
| moller2 | all | ok | 638237 | 626127 | g2 2.210 | -1.243 ± 0.360 | nan |
| moller | pool | ok | 7302 | 7301 | g2 4.683 | 0.308 ± 0.114 | 0.161 |
| moller | I1 | ok | 6521 | 6521 | g2 4.533 | 0.370 ± 0.128 | 0.271 |
| moller | I2 | no-fit | 568 | 568 | g2 nan | nan ± nan | nan |
| moller | I3 | no-fit | 213 | 212 | g2 nan | nan ± nan | nan |
| moller | M | no-fit | 75 | 74 | g2 nan | nan ± nan | nan |

T clusters (legacy): with lead-glass coincidence 1326 (27.5 %), median E/E' 0.9230 (deficit 7.70 %); without 3497, median 0.9808 (deficit 1.92 %).

### drift: mu of the pool, legacy, against split order per run

| run | splits | mu mean (%) | slope (pp / split) ± SE | mu first-last (pp) |
|---|---|---|---|---|
| 24698 | 8 | 0.148 | -0.0001 ± 0.0001 | 0.007 |
| 24911 | 8 | 0.146 | 0.0000 ± 0.0000 | -0.004 |
| 24937 | 0 | - | - | - |

Continuity: pool legacy res_cb 1.5854 vs baseline pre 1.5854.

## 3p49GeV_H2_old


### ep, path legacy

| region | status | n_range | modules | prior frac | res_cb (%) | res_cb oos | mu (%) | sigma_E (MeV) | sep_pi | tail_hi3 | leak_cut | eff_cut |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pool | ok | 548305 | 169 | 0.035 | 1.559 ± 0.006 | 1.579 | 0.164 ± 0.013 | 54.3 | 2.66 | 0.120 | 22.0 | 84.8 |
| I1 | ok | 484093 | 76 | 0.021 | 1.554 ± 0.006 | 1.571 | 0.112 ± 0.013 | 54.2 | 2.67 | 0.124 | 21.9 | 84.9 |
| I2 | ok | 62948 | 75 | 0.067 | 1.626 ± 0.014 | 1.683 | 0.040 ± 0.019 | 56.6 | 2.55 | 0.061 | 22.2 | 83.2 |
| I3 | ep-starved | 1264 | 18 | 0.104 | 1.598 ± 0.130 | 1.627 | 0.091 ± 0.084 | 55.6 | 2.78 | 0.824 | 24.2 | 87.7 |
| pool_edge | ok | 194028 | 169 | 0.035 | 1.630 ± 0.013 | 1.661 | 0.261 ± 0.018 | 56.8 | 2.55 | 0.083 | 21.8 | 81.4 |
| pool_centre | ok | 354277 | 169 | 0.035 | 1.510 ± 0.008 | 1.528 | 0.129 ± 0.014 | 52.6 | 2.75 | 0.136 | 21.7 | 86.6 |
| H | ok | 1316620 | 32 | 0.476 | 2.870 ± 0.084 | 3.051 | -1.750 ± 0.075 | 100.0 | 1.48 | 0.029 | 19.4 | 47.4 |
| M | ok | 34008 | 348 | 0.987 | 2.062 ± 0.044 | 2.100 | 0.743 ± 0.159 | 71.7 | 2.04 | 0.369 | 25.8 | 76.4 |
| O | ok | 6434 | 252 | 1.000 | 2.050 ± 0.072 | 2.046 | 1.596 ± 0.184 | 71.1 | 2.07 | 0.593 | 23.6 | 70.8 |
| T | ok | 4079 | 253 | 1.000 | 2.201 ± 0.117 | 2.192 | 1.329 ± 0.181 | 76.2 | 2.00 | 0.789 | 19.2 | 52.5 |
| T_lg | not-fitted | 1144 | 167 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |
| T_nolg | not-fitted | 2935 | 243 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |

### ep, path g12

| region | status | n_range | modules | prior frac | res_cb (%) | res_cb oos | mu (%) | sigma_E (MeV) | sep_pi | tail_hi3 | leak_cut | eff_cut |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pool | ok | 548224 | 170 | 0.035 | 1.548 ± 0.007 | 1.570 | 0.203 ± 0.013 | 53.9 | 2.68 | 0.101 | 21.8 | 84.9 |
| I1 | ok | 484089 | 76 | 0.021 | 1.544 ± 0.006 | 1.562 | 0.148 ± 0.013 | 53.8 | 2.69 | 0.105 | 21.8 | 85.1 |
| I2 | ok | 62806 | 75 | 0.066 | 1.610 ± 0.015 | 1.672 | 0.110 ± 0.018 | 56.1 | 2.58 | 0.051 | 22.0 | 83.4 |
| I3 | ep-starved | 1329 | 19 | 0.106 | 1.529 ± 0.122 | 1.608 | 0.183 ± 0.090 | 53.2 | 2.91 | 0.883 | 22.4 | 87.8 |
| pool_edge | ok | 195103 | 170 | 0.035 | 1.615 ± 0.014 | 1.650 | 0.300 ± 0.018 | 56.2 | 2.57 | 0.077 | 21.7 | 81.6 |
| pool_centre | ok | 353121 | 170 | 0.035 | 1.500 ± 0.008 | 1.514 | 0.168 ± 0.014 | 52.2 | 2.77 | 0.110 | 21.6 | 86.7 |
| H | ok | 1316278 | 32 | 0.476 | 2.852 ± 0.087 | 3.039 | -1.673 ± 0.077 | 99.4 | 1.49 | 0.018 | 19.3 | 47.3 |
| M | ok | 34002 | 348 | 0.987 | 2.057 ± 0.041 | 2.088 | 0.797 ± 0.157 | 71.5 | 2.04 | 0.349 | 25.7 | 76.5 |
| O | ok | 6422 | 252 | 1.000 | 2.046 ± 0.080 | 2.038 | 1.632 ± 0.185 | 71.0 | 2.09 | 0.656 | 23.4 | 71.0 |
| T | ok | 4067 | 253 | 1.000 | 2.250 ± 0.113 | 2.257 | 1.310 ± 0.180 | 77.9 | 1.95 | 0.913 | 20.2 | 53.6 |
| T_lg | not-fitted | 1140 | 167 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |
| T_nolg | not-fitted | 2927 | 243 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |

### moller2 ((E1+E2)/E_beam) and moller (E/E_ee), path legacy

| kind | region | status | n | n_range | sigma/E (%) | mu (%) | mu(moller) - mu(ep) (pp) |
|---|---|---|---|---|---|---|---|
| moller2 | all | ok | 656893 | 645505 | g2 2.270 | -2.153 ± 0.009 | nan |
| moller | pool | ok | 10534 | 10534 | g2 4.159 | 0.716 ± 0.113 | 0.552 |
| moller | I1 | ok | 9763 | 9763 | g2 4.105 | 0.741 ± 0.114 | 0.629 |
| moller | I2 | no-fit | 592 | 592 | g2 nan | nan ± nan | nan |
| moller | I3 | no-fit | 179 | 179 | g2 nan | nan ± nan | nan |
| moller | M | no-fit | 35 | 34 | g2 nan | nan ± nan | nan |

T clusters (legacy): with lead-glass coincidence 1342 (29.9 %), median E/E' 0.9214 (deficit 7.86 %); without 3143, median 0.9980 (deficit 0.20 %).

### drift: mu of the pool, legacy, against split order per run

| run | splits | mu mean (%) | slope (pp / split) ± SE | mu first-last (pp) |
|---|---|---|---|---|
| 24386 | 11 | 0.186 | 0.0010 ± 0.0009 | -0.046 |
| 24503 | 10 | 0.155 | -0.0007 ± 0.0007 | -0.051 |
| 24512 | 7 | 0.143 | -0.0006 ± 0.0006 | -0.011 |

Continuity: pool legacy res_cb 1.5591 vs baseline pre 1.5591.

## 2p24GeV_Ta_X17


### ep, path legacy

| region | status | n_range | modules | prior frac | res_cb (%) | res_cb oos | mu (%) | sigma_E (MeV) | sep_pi | tail_hi3 | leak_cut | eff_cut |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pool | ok | 56164 | 46 | 0.063 | 2.308 ± 0.035 | 2.332 | 0.009 ± 0.028 | 51.6 | 2.81 | 4.484 | 11.3 | 51.6 |
| I1 | ok | 56164 | 46 | 0.063 | 2.308 ± 0.035 | 2.332 | 0.009 ± 0.028 | 51.6 | 2.81 | 4.484 | 11.3 | 51.6 |
| I2 | ep-starved | 0 | 0 | nan | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |
| I3 | ep-starved | 0 | 0 | nan | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |
| pool_edge | ok | 21513 | 46 | 0.063 | 2.333 ± 0.060 | 2.396 | -0.050 ± 0.045 | 52.2 | 2.82 | 4.891 | 11.4 | 47.4 |
| pool_centre | ok | 34651 | 46 | 0.063 | 2.265 ± 0.036 | 2.283 | 0.052 ± 0.033 | 50.7 | 2.87 | 4.274 | 11.3 | 54.2 |
| H | ok | 387636 | 32 | 0.740 | 19.238 ± 1.813 | 19.118 | -35.618 ± 1.661 | 430.4 | 0.48 | 17.846 | 85.0 | 97.0 |
| M | ok | 2289 | 348 | 1.000 | 1.985 ± 0.125 | 2.052 | -0.481 ± 0.516 | 44.4 | 3.69 | 0.364 | 5.8 | 62.2 |
| O | ep-starved | 222 | 241 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |
| T | ep-starved | 95 | 199 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |
| T_lg | not-fitted | 13 | 55 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |
| T_nolg | not-fitted | 82 | 178 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |

### ep, path g8

| region | status | n_range | modules | prior frac | res_cb (%) | res_cb oos | mu (%) | sigma_E (MeV) | sep_pi | tail_hi3 | leak_cut | eff_cut |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pool | ok | 54835 | 47 | 0.063 | 2.035 ± 0.022 | 2.054 | -0.376 ± 0.044 | 45.5 | 3.18 | 0.583 | 5.2 | 51.4 |
| I1 | ok | 54835 | 47 | 0.063 | 2.035 ± 0.022 | 2.054 | -0.376 ± 0.044 | 45.5 | 3.18 | 0.583 | 5.2 | 51.4 |
| I2 | ep-starved | 0 | 0 | nan | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |
| I3 | ep-starved | 0 | 0 | nan | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |
| pool_edge | ok | 20857 | 47 | 0.063 | 2.087 ± 0.036 | 2.106 | -0.384 ± 0.051 | 46.7 | 3.12 | 0.520 | 5.2 | 47.4 |
| pool_centre | ok | 33978 | 47 | 0.063 | 2.008 ± 0.025 | 2.027 | -0.377 ± 0.043 | 44.9 | 3.23 | 0.617 | 5.2 | 54.0 |
| H | ok | 359024 | 32 | 0.740 | 13.908 ± 2.070 | 12.897 | -26.995 ± 3.620 | 311.2 | 1.15 | 15.460 | 64.6 | 78.8 |
| M | ok | 2307 | 348 | 1.000 | 2.116 ± 0.099 | 2.096 | -0.690 ± 0.549 | 47.3 | 3.17 | 0.124 | 6.3 | 62.1 |
| O | ep-starved | 223 | 241 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |
| T | ep-starved | 96 | 197 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |
| T_lg | not-fitted | 14 | 56 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |
| T_nolg | not-fitted | 82 | 178 | 1.000 | nan ± nan | nan | nan ± nan | nan | nan | nan | nan | nan |

### moller2 ((E1+E2)/E_beam) and moller (E/E_ee), path legacy

| kind | region | status | n | n_range | sigma/E (%) | mu (%) | mu(moller) - mu(ep) (pp) |
|---|---|---|---|---|---|---|---|
| moller2 | all | ok | 81127 | 54867 | g2 3.526 | -3.280 ± 0.093 | nan |
| moller | pool | ok | 16099 | 16099 | g2 4.520 | -1.630 ± 0.103 | -1.639 |
| moller | I1 | ok | 14252 | 14252 | g2 4.183 | -1.768 ± 0.125 | -1.777 |
| moller | I2 | ok | 1243 | 1243 | g2 6.463 | -0.867 ± 0.513 | nan |
| moller | I3 | no-fit | 604 | 604 | g2 nan | nan ± nan | nan |
| moller | M | no-fit | 820 | 805 | g2 nan | nan ± nan | nan |
| moller | O | no-fit | 77 | 47 | g2 nan | nan ± nan | nan |
| moller | T | no-fit | 18 | 5 | g2 nan | nan ± nan | nan |

T clusters (legacy): with lead-glass coincidence 85 (17.3 %), median E/E' 0.8137 (deficit 18.63 %); without 405, median 0.8817 (deficit 11.83 %).

### drift: mu of the pool, legacy, against split order per run

| run | splits | mu mean (%) | slope (pp / split) ± SE | mu first-last (pp) |
|---|---|---|---|---|
| 25800 | 8 | 0.022 | -0.0003 ± 0.0001 | -0.276 |
| 26138 | 4 | 0.045 | 0.0004 ± 0.0013 | 0.023 |
| 26230 | 13 | -0.010 | -0.0002 ± 0.0001 | -0.032 |

Continuity: pool legacy res_cb 2.3081 vs baseline pre 2.3081.
