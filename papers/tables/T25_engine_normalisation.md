**T25_engine_normalisation.** The composition-versus-full rank correlation under each of the rate engine's row normalisations, and under each of the two glutamate sign conventions the two levels of the ablation ladder use. The shipped engine divides every row of the signed weight matrix by its own total absolute input, which makes each cell's recurrent input a composition-weighted average of its presynaptic gains; that is a modelling choice, it was documented nowhere, and it is what the composition-dominance verdict is measured on. Levels are not comparable across normalisations -- the unnormalised operator is supercritical and the r_max clip shapes its rates -- only the orderings are.

| conc_M | normalisation | spearman_rho_b_vs_d | reproduces_full_ordering |
|---|---|---|---|
| 1.00e-06 | row_abs | 0.921 | True |
| 1.00e-06 | none | 0.935 | True |
| 1.00e-06 | degree | 0.551 | False |
| 1.00e-08 | row_abs | 0.876 | False |
| 1.00e-08 | none | 0.909 | True |
| 1.00e-08 | degree | -0.583 | False |
| 1.00e-08 | glutamate sign: wholens (B excitatory) | 0.876 | -- |
| 1.00e-08 | glutamate sign: rate_engine (B inhibitory) | 0.963 | -- |
| 1.00e-07 | glutamate sign: wholens (B excitatory) | 0.918 | -- |
| 1.00e-07 | glutamate sign: rate_engine (B inhibitory) | 0.992 | -- |
| 1.00e-06 | glutamate sign: wholens (B excitatory) | 0.921 | -- |
| 1.00e-06 | glutamate sign: rate_engine (B inhibitory) | 0.991 | -- |
| 1.00e-05 | glutamate sign: wholens (B excitatory) | 0.954 | -- |
| 1.00e-05 | glutamate sign: rate_engine (B inhibitory) | 0.992 | -- |
