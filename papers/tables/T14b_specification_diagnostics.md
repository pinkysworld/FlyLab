**T14b_specification_diagnostics.** Per-specification diagnostics behind T14: the gains each specification produces for the two reference compounds at 1 uM, and the mean circuit-minus-receptor selectivity gap over the nicotinic and Nav/AChE sets. An empty gap means no compound of that set crossed the circuit threshold under that specification, so the index is undefined rather than zero.

| spec | g_ach_imidacloprid | g_gaba_fipronil | mean_si_gap_nicotinic | mean_si_gap_nav_ache | runtime_s |
|---|---|---|---|---|---|
| linear@0.50x | 1.198 | 0.507 | -- | -0.088 | 19.9 |
| linear@0.75x | 1.297 | 0.261 | -- | 0.117 | 21.0 |
| linear@1.00x | 1.396 | 0.050 | -1.734 | 0.284 | 21.7 |
| linear@1.25x | 1.495 | 0.050 | -1.042 | 0.413 | 22.6 |
| linear@1.50x | 1.595 | 0.050 | -0.790 | 0.504 | 23.5 |
| saturating@0.50x | 1.200 | 0.501 | -- | 0.588 | 24.1 |
| saturating@0.75x | 1.299 | 0.252 | -- | 0.837 | 27.1 |
| saturating@1.00x | 1.399 | 0.050 | -1.150 | 1.000 | 31.0 |
| saturating@1.25x | 1.499 | 0.050 | -0.443 | 1.091 | 28.7 |
| saturating@1.50x | 1.599 | 0.050 | -0.199 | 1.170 | 26.9 |
| weak_biphasic@0.50x | 0.805 | 0.265 | -- | -0.228 | 25.1 |
| weak_biphasic@0.75x | 0.708 | 0.050 | -- | 0.059 | 24.7 |
| weak_biphasic@1.00x | 0.611 | 0.050 | -- | 0.232 | 26.5 |
| weak_biphasic@1.25x | 0.514 | 0.050 | -2.407 | 0.372 | 26.1 |
| weak_biphasic@1.50x | 0.416 | 0.050 | -1.617 | 0.475 | 26.9 |
| flylab_biphasic@0.50x | 0.413 | 0.507 | -1.510 | -0.101 | 24.2 |
| flylab_biphasic@0.75x | 0.119 | 0.261 | -1.060 | 0.112 | 24.6 |
| flylab_biphasic@1.00x | 0.050 | 0.050 | -0.884 | 0.277 | 23.5 |
| flylab_biphasic@1.25x | 0.050 | 0.050 | -0.779 | 0.408 | 23.9 |
| flylab_biphasic@1.50x | 0.050 | 0.050 | -0.706 | 0.501 | 20.7 |
| strong_biphasic@0.50x | 0.050 | 0.050 | -0.751 | -1.008 | 20.8 |
| strong_biphasic@0.75x | 0.050 | 0.050 | -0.580 | -0.830 | 13.4 |
| strong_biphasic@1.00x | 0.050 | 0.050 | -0.507 | -0.489 | 11.9 |
| strong_biphasic@1.25x | 0.050 | 0.050 | -0.465 | 0.181 | 11.8 |
| strong_biphasic@1.50x | 0.050 | 0.050 | -0.424 | 0.312 | 11.3 |
