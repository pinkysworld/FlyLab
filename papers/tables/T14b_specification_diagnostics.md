**T14b_specification_diagnostics.** Per-specification diagnostics behind T14: the gains each specification produces for the two reference compounds at 1 uM, and the mean circuit-minus-receptor selectivity gap over the nicotinic and Nav/AChE sets. An empty gap means no compound of that set crossed the circuit threshold under that specification, so the index is undefined rather than zero.

| spec | g_ach_imidacloprid | g_gaba_fipronil | mean_si_gap_nicotinic | mean_si_gap_nav_ache | runtime_s |
|---|---|---|---|---|---|
| linear@0.50x | 1.198 | 0.507 | -- | -0.088 | 13.8 |
| linear@0.75x | 1.297 | 0.261 | -- | 0.117 | 13.3 |
| linear@1.00x | 1.396 | 0.050 | -1.734 | 0.284 | 10.3 |
| linear@1.25x | 1.495 | 0.050 | -1.042 | 0.413 | 9.690 |
| linear@1.50x | 1.595 | 0.050 | -0.790 | 0.504 | 9.780 |
| saturating@0.50x | 1.200 | 0.501 | -- | 0.588 | 9.780 |
| saturating@0.75x | 1.299 | 0.252 | -- | 0.837 | 9.710 |
| saturating@1.00x | 1.399 | 0.050 | -1.150 | 1.000 | 10.2 |
| saturating@1.25x | 1.499 | 0.050 | -0.443 | 1.091 | 10.3 |
| saturating@1.50x | 1.599 | 0.050 | -0.199 | 1.170 | 10.5 |
| weak_biphasic@0.50x | 0.805 | 0.265 | -- | -0.228 | 9.830 |
| weak_biphasic@0.75x | 0.708 | 0.050 | -- | 0.059 | 9.230 |
| weak_biphasic@1.00x | 0.611 | 0.050 | -- | 0.232 | 9.290 |
| weak_biphasic@1.25x | 0.514 | 0.050 | -2.407 | 0.372 | 9.420 |
| weak_biphasic@1.50x | 0.416 | 0.050 | -1.617 | 0.475 | 9.770 |
| flylab_biphasic@0.50x | 0.413 | 0.507 | -1.510 | -0.101 | 9.570 |
| flylab_biphasic@0.75x | 0.119 | 0.261 | -1.060 | 0.112 | 10.1 |
| flylab_biphasic@1.00x | 0.050 | 0.050 | -0.884 | 0.277 | 10.3 |
| flylab_biphasic@1.25x | 0.050 | 0.050 | -0.779 | 0.408 | 10.6 |
| flylab_biphasic@1.50x | 0.050 | 0.050 | -0.706 | 0.501 | 10.4 |
| strong_biphasic@0.50x | 0.050 | 0.050 | -0.751 | -1.008 | 10.7 |
| strong_biphasic@0.75x | 0.050 | 0.050 | -0.580 | -0.830 | 10.0 |
| strong_biphasic@1.00x | 0.050 | 0.050 | -0.507 | -0.489 | 9.620 |
| strong_biphasic@1.25x | 0.050 | 0.050 | -0.465 | 0.181 | 10.1 |
| strong_biphasic@1.50x | 0.050 | 0.050 | -0.424 | 0.312 | 9.820 |
