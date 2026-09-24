**T14b_specification_diagnostics.** Per-specification diagnostics behind T14: the gains each specification produces for the two reference compounds at 1 uM, and the mean circuit-minus-receptor selectivity gap over the nicotinic and Nav/AChE sets. An empty gap means no compound of that set crossed the circuit threshold under that specification, so the index is undefined rather than zero.

| spec | g_ach_imidacloprid | g_gaba_fipronil | mean_si_gap_nicotinic | mean_si_gap_nav_ache | runtime_s |
|---|---|---|---|---|---|
| linear@0.50x | 1.198 | 0.507 | -- | -0.189 | 15.2 |
| linear@0.75x | 1.297 | 0.261 | -- | 0.040 | 15.7 |
| linear@1.00x | 1.396 | 0.050 | -1.734 | 0.198 | 16.0 |
| linear@1.25x | 1.495 | 0.050 | -1.042 | 0.325 | 17.5 |
| linear@1.50x | 1.595 | 0.050 | -0.790 | 0.424 | 16.6 |
| saturating@0.50x | 1.200 | 0.501 | -- | 0.503 | 18.4 |
| saturating@0.75x | 1.299 | 0.252 | -- | 0.738 | 16.9 |
| saturating@1.00x | 1.399 | 0.050 | -1.150 | 0.916 | 18.0 |
| saturating@1.25x | 1.499 | 0.050 | -0.443 | 1.018 | 17.5 |
| saturating@1.50x | 1.599 | 0.050 | -0.199 | 1.098 | 18.3 |
| weak_biphasic@0.50x | 0.805 | 0.265 | -- | -0.346 | 17.8 |
| weak_biphasic@0.75x | 0.708 | 0.050 | -- | -0.028 | 18.0 |
| weak_biphasic@1.00x | 0.611 | 0.050 | -- | 0.148 | 18.1 |
| weak_biphasic@1.25x | 0.514 | 0.050 | -2.407 | 0.279 | 17.6 |
| weak_biphasic@1.50x | 0.416 | 0.050 | -1.617 | 0.390 | 18.0 |
| flylab_biphasic@0.50x | 0.413 | 0.507 | -1.510 | -0.189 | 17.6 |
| flylab_biphasic@0.75x | 0.119 | 0.261 | -1.060 | 0.040 | 18.2 |
| flylab_biphasic@1.00x | 0.050 | 0.050 | -0.884 | 0.198 | 18.6 |
| flylab_biphasic@1.25x | 0.050 | 0.050 | -0.779 | 0.325 | 13.5 |
| flylab_biphasic@1.50x | 0.050 | 0.050 | -0.706 | 0.424 | 13.0 |
| strong_biphasic@0.50x | 0.050 | 0.050 | -0.751 | -1.135 | 11.7 |
| strong_biphasic@0.75x | 0.050 | 0.050 | -0.580 | -0.930 | 10.9 |
| strong_biphasic@1.00x | 0.050 | 0.050 | -0.507 | -0.792 | 11.2 |
| strong_biphasic@1.25x | 0.050 | 0.050 | -0.465 | 0.032 | 10.7 |
| strong_biphasic@1.50x | 0.050 | 0.050 | -0.424 | 0.208 | 10.9 |
