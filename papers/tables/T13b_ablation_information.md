**T13b_ablation_information.** How much of the full model's ordering of the library each ablated level recovers, per concentration.

| conc_M | level | unit | spearman_rho_vs_full | pearson_r_vs_full | residual_rms_standardised | residual_rms_hz | reproduces_full_ordering | information_added_vs_previous |
|---|---|---|---|---|---|---|---|---|
| 1.00e-08 | A_receptor_only | dimensionless | -0.046 | -0.403 | 1.675 | -- | False | -- |
| 1.00e-08 | B_composition_only | excitation index | 0.906 | 0.998 | 0.068 | -- | True | 0.860 |
| 1.00e-08 | C_topology_only | Hz | -0.427 | -0.369 | 1.655 | 10.1 | False | -0.479 |
| 1.00e-08 | D_full_flylab | Hz | 1.000 | 1.000 | 0.000 | 0.000 | True | 0.573 |
| 1.00e-07 | A_receptor_only | dimensionless | -0.348 | 0.001 | 1.413 | -- | False | -- |
| 1.00e-07 | B_composition_only | excitation index | 0.989 | 0.717 | 0.752 | -- | True | 0.640 |
| 1.00e-07 | C_topology_only | Hz | 0.279 | 0.054 | 1.375 | 6.696 | False | -0.710 |
| 1.00e-07 | D_full_flylab | Hz | 1.000 | 1.000 | 0.000 | 0.000 | True | 0.721 |
| 1.00e-06 | A_receptor_only | dimensionless | -0.431 | 0.082 | 1.355 | -- | False | -- |
| 1.00e-06 | B_composition_only | excitation index | 0.988 | 0.756 | 0.699 | -- | True | 0.556 |
| 1.00e-06 | C_topology_only | Hz | 0.242 | -0.070 | 1.463 | 11.7 | False | -0.746 |
| 1.00e-06 | D_full_flylab | Hz | 1.000 | 1.000 | 0.000 | 0.000 | True | 0.758 |
| 1.00e-05 | A_receptor_only | dimensionless | -0.506 | 0.086 | 1.352 | -- | False | -- |
| 1.00e-05 | B_composition_only | excitation index | 0.987 | 0.754 | 0.702 | -- | True | 0.481 |
| 1.00e-05 | C_topology_only | Hz | 0.236 | -0.118 | 1.495 | 13.9 | False | -0.751 |
| 1.00e-05 | D_full_flylab | Hz | 1.000 | 1.000 | 0.000 | 0.000 | True | 0.764 |
