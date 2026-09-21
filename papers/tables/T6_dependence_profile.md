**T6_dependence_profile.** Connectome-dependence profile. `real_effect` is treated minus vehicle on the real MaleCNS cut; the null is the same contrast on `n` degraded copies of that cut, with the seed block and node order held fixed. `p_two_sided` is the empirical permutation probability (k+1)/(n+1) and is the statistic to read first; `p_resolution` is its floor, and `z` is a standardised distance from a usually non-normal null. `verdict` records `distinguishable` when the test rejects. The legacy API label `equivalent_within_tolerance` means only that the point estimate of the absolute real-effect to null-median gap is below the configured margin `delta`; it is a descriptive point-gap diagnostic, not a formal equivalence test. `indeterminate` means neither condition was met. No interval for the gap, TOST, distributional closeness criterion, or error-controlled equivalence decision is computed. `sign_permute` is a **joint** target-set-and-sign null and is off the ladder (`on_ladder = False`); the rank-3 rung is the weight-matched transmitter null, which holds each transmitter's share of total outgoing weight fixed (T23). The legacy within-tolerance label is a point-gap diagnostic, not a formal equivalence test. Rate-engine quantities labelled Hz in legacy fields or axes are rate-model units, not calibrated physiological firing rates.

| compound | assay | mode | on_ladder | real_effect | p_two_sided | verdict | abs_gap_from_null_median | delta |
|---|---|---|---|---|---|---|---|---|
| imidacloprid | subgraph | sign_permute | False | -6.173 | 0.275 | indeterminate | 2.652 | 0.332 |
| imidacloprid | subgraph | weight_permute | True | -6.173 | 0.586 | equivalent_within_tolerance | 0.126 | 0.332 |
| imidacloprid | subgraph | rewire_degree_preserving | True | -6.173 | 0.472 | equivalent_within_tolerance | 0.123 | 0.332 |
| imidacloprid | subgraph | erdos_renyi | True | -6.173 | 9.99e-04 | distinguishable | 6.015 | 0.332 |
| imidacloprid | subgraph | sign_permute_weight_matched | True | -6.173 | 0.389 | equivalent_within_tolerance | 0.236 | 0.332 |
| fipronil | subgraph | sign_permute | False | 0.900 | 0.113 | indeterminate | 0.419 | 0.332 |
| fipronil | subgraph | weight_permute | True | 0.900 | 0.006 | distinguishable | 0.310 | 0.332 |
| fipronil | subgraph | rewire_degree_preserving | True | 0.900 | 0.007 | distinguishable | 0.243 | 0.332 |
| fipronil | subgraph | erdos_renyi | True | 0.900 | 9.99e-04 | distinguishable | 0.900 | 0.332 |
| fipronil | subgraph | sign_permute_weight_matched | True | 0.900 | 0.210 | equivalent_within_tolerance | 0.060 | 0.332 |
| fipronil | taste_map | sign_permute | False | 0.413 | 0.417 | indeterminate | 0.197 | -- |
| fipronil | taste_map | weight_permute | True | 0.413 | 0.879 | indeterminate | 0.060 | -- |
| fipronil | taste_map | rewire_degree_preserving | True | 0.413 | 0.203 | indeterminate | 0.326 | -- |
| fipronil | taste_map | erdos_renyi | True | 0.413 | 0.245 | indeterminate | 0.345 | -- |
| fipronil | taste_map | sign_permute_weight_matched | True | 0.413 | 0.498 | indeterminate | 0.031 | -- |
