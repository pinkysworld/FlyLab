**T13b_ablation_information.** How much of the full model's ordering of the library each ablated level recovers, per concentration. `reproduces_full_ordering` is applied to the **signed** rank correlation: an ordering that is a perfect inversion of the full model's reproduces nothing, and the previous absolute-value test would have credited it.

| conc_M | level | unit | spearman_rho_vs_full | pearson_r_vs_full | residual_rms_standardised | residual_rms_hz | reproduces_full_ordering | information_added_vs_previous |
|---|---|---|---|---|---|---|---|---|
| 1.00e-08 | A_receptor_only | dimensionless | -0.046 | -0.289 | 1.606 | -- | False | -- |
| 1.00e-08 | B_composition_only | excitation index | 0.906 | 0.998 | 0.068 | -- | True | 0.952 |
| 1.00e-08 | C_topology_only | Hz | 0.928 | 0.470 | 1.030 | 8.515 | True | 0.022 |
| 1.00e-08 | D_full_flylab | Hz | 1.000 | 1.000 | 0.000 | 0.000 | True | 0.072 |
| 1.00e-08 | C_topology_only_floor | Hz | -0.427 | -0.369 | 1.655 | 10.1 | False | -- |
| 1.00e-07 | A_receptor_only | dimensionless | -0.348 | -0.064 | 1.459 | -- | False | -- |
| 1.00e-07 | B_composition_only | excitation index | 0.989 | 0.717 | 0.752 | -- | True | 1.337 |
| 1.00e-07 | C_topology_only | Hz | 0.963 | 0.638 | 0.851 | 4.593 | True | -0.026 |
| 1.00e-07 | D_full_flylab | Hz | 1.000 | 1.000 | 0.000 | 0.000 | True | 0.037 |
| 1.00e-07 | C_topology_only_floor | Hz | 0.279 | 0.054 | 1.375 | 6.696 | False | -- |
| 1.00e-06 | A_receptor_only | dimensionless | -0.431 | 0.031 | 1.392 | -- | False | -- |
| 1.00e-06 | B_composition_only | excitation index | 0.988 | 0.756 | 0.699 | -- | True | 1.419 |
| 1.00e-06 | C_topology_only | Hz | 0.961 | 0.654 | 0.832 | 8.149 | True | -0.027 |
| 1.00e-06 | D_full_flylab | Hz | 1.000 | 1.000 | 0.000 | 0.000 | True | 0.039 |
| 1.00e-06 | C_topology_only_floor | Hz | 0.242 | -0.070 | 1.463 | 11.7 | False | -- |
| 1.00e-05 | A_receptor_only | dimensionless | -0.506 | 0.044 | 1.383 | -- | False | -- |
| 1.00e-05 | B_composition_only | excitation index | 0.987 | 0.754 | 0.702 | -- | True | 1.493 |
| 1.00e-05 | C_topology_only | Hz | 0.939 | 0.651 | 0.836 | 9.729 | True | -0.048 |
| 1.00e-05 | D_full_flylab | Hz | 1.000 | 1.000 | 0.000 | 0.000 | True | 0.061 |
| 1.00e-05 | C_topology_only_floor | Hz | 0.236 | -0.118 | 1.495 | 13.9 | False | -- |
