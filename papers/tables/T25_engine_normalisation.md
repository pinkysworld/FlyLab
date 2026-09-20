**T25_engine_normalisation.** The composition-versus-full rank correlation under each of the rate engine's row normalisations, and under each of the two glutamate sign conventions the two levels of the ablation ladder use. The shipped engine divides every row of the signed weight matrix by its own total absolute input, which makes each cell's recurrent input a composition-weighted average of its presynaptic gains; that is a modelling choice, it was documented nowhere, and it is what the composition-dominance verdict is measured on. Levels are not comparable across normalisations -- the unnormalised operator is supercritical and the r_max clip shapes its rates -- only the orderings are.

| conc_M | normalisation | spearman_rho_b_vs_d | reproduces_full_ordering |
|---|---|---|---|
| 1.00e-06 | row_abs | 0.988 | True |
| 1.00e-06 | none | 0.906 | True |
| 1.00e-06 | degree | 0.666 | False |
| 1.00e-08 | row_abs | 0.906 | True |
| 1.00e-08 | none | 0.921 | True |
| 1.00e-08 | degree | -0.372 | False |
| 1.00e-08 | glutamate sign: wholens (B excitatory) | 0.906 | -- |
| 1.00e-08 | glutamate sign: rate_engine (B inhibitory) | 0.964 | -- |
| 1.00e-07 | glutamate sign: wholens (B excitatory) | 0.989 | -- |
| 1.00e-07 | glutamate sign: rate_engine (B inhibitory) | 0.987 | -- |
| 1.00e-06 | glutamate sign: wholens (B excitatory) | 0.988 | -- |
| 1.00e-06 | glutamate sign: rate_engine (B inhibitory) | 0.984 | -- |
| 1.00e-05 | glutamate sign: wholens (B excitatory) | 0.987 | -- |
| 1.00e-05 | glutamate sign: rate_engine (B inhibitory) | 0.979 | -- |
