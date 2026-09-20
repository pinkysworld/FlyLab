**T16_uncertainty_budget.** Global uncertainty budget for the neighbourhood mean rate (imidacloprid at 1 µM, 11264 model evaluations, Jansen estimators on a Saltelli cross-sample). S1 is the variance share resolving that factor alone would remove; ST includes its interactions. The shares do not sum to one; the remainder is the interaction row. A null factor with no effect on the deterministic rate engine is included on purpose and reports the estimator's noise floor at this sample size.

| source | kind | share_of_variance_S1 | share_with_interactions_ST | variance_removed_hz2 | ci_low | ci_high |
|---|---|---|---|---|---|---|
| gain_transform | categorical | 0.410 | 0.626 | 1.290 | 0.313 | 0.500 |
| weight_threshold | continuous | 0.294 | 0.542 | 0.925 | 0.199 | 0.391 |
| drive | continuous | 0.019 | 0.141 | 0.058 | -0.072 | 0.112 |
| transmitter | continuous | 0.006 | 0.024 | 0.019 | -0.052 | 0.070 |
| expression | continuous | -0.006 | 1.28e-04 | -0.018 | -0.059 | 0.044 |
| potency | continuous | -0.006 | 6.78e-04 | -0.019 | -0.061 | 0.046 |
| lif_seed | seed | -0.007 | 0.000 | -0.021 | -0.062 | 0.043 |
| hill_n | continuous | -0.008 | 7.53e-04 | -0.025 | -0.062 | 0.040 |
| gain_coef | continuous | -0.046 | 0.024 | -0.144 | -0.107 | 0.028 |
| interactions (higher order) | residual | 0.344 | -- | 1.084 | -- | -- |
