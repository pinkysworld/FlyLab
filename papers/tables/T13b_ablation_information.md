**T13b_ablation_information.** How much of the full model's ordering of the library each ablated level recovers, per concentration. `reproduces_full_ordering` is applied to the **signed** rank correlation: an ordering that is a perfect inversion of the full model's reproduces nothing, and the previous absolute-value test would have credited it.

| conc_M | level | unit | spearman_rho_vs_full | pearson_r_vs_full | residual_rms_standardised | residual_rms_hz | reproduces_full_ordering | information_added_vs_previous |
|---|---|---|---|---|---|---|---|---|
| 1.00e-08 | A_receptor_only | dimensionless | 0.095 | -0.033 | 1.437 | -- | False | -- |
| 1.00e-08 | B_composition_only | excitation index | 0.876 | 0.760 | 0.692 | -- | False | 0.781 |
| 1.00e-08 | C_topology_only | Hz | 0.924 | 0.268 | 1.210 | 1.619 | True | 0.048 |
| 1.00e-08 | D_full_flylab | Hz | 1.000 | 1.000 | 0.000 | 0.000 | True | 0.076 |
| 1.00e-08 | C_topology_only_floor | Hz | -0.371 | -0.043 | 1.444 | 1.958 | False | -- |
| 1.00e-07 | A_receptor_only | dimensionless | -0.519 | -0.136 | 1.507 | -- | False | -- |
| 1.00e-07 | B_composition_only | excitation index | 0.918 | 0.706 | 0.767 | -- | True | 1.438 |
| 1.00e-07 | C_topology_only | Hz | 0.962 | 0.620 | 0.872 | 4.584 | True | 0.044 |
| 1.00e-07 | D_full_flylab | Hz | 1.000 | 1.000 | 0.000 | 0.000 | True | 0.038 |
| 1.00e-07 | C_topology_only_floor | Hz | 0.175 | -0.014 | 1.424 | 6.690 | False | -- |
| 1.00e-06 | A_receptor_only | dimensionless | -0.564 | -0.008 | 1.420 | -- | False | -- |
| 1.00e-06 | B_composition_only | excitation index | 0.921 | 0.750 | 0.707 | -- | True | 1.486 |
| 1.00e-06 | C_topology_only | Hz | 0.968 | 0.645 | 0.842 | 8.151 | True | 0.046 |
| 1.00e-06 | D_full_flylab | Hz | 1.000 | 1.000 | 0.000 | 0.000 | True | 0.032 |
| 1.00e-06 | C_topology_only_floor | Hz | 0.161 | -0.116 | 1.494 | 11.7 | False | -- |
| 1.00e-05 | A_receptor_only | dimensionless | -0.579 | 0.042 | 1.385 | -- | False | -- |
| 1.00e-05 | B_composition_only | excitation index | 0.954 | 0.750 | 0.707 | -- | True | 1.533 |
| 1.00e-05 | C_topology_only | Hz | 0.950 | 0.636 | 0.853 | 9.731 | True | -0.004 |
| 1.00e-05 | D_full_flylab | Hz | 1.000 | 1.000 | 0.000 | 0.000 | True | 0.050 |
| 1.00e-05 | C_topology_only_floor | Hz | 0.151 | -0.165 | 1.527 | 13.9 | False | -- |
