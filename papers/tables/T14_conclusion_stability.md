**T14_conclusion_stability.** Conclusion stability matrix: every principal conclusion re-derived under each of 25 prespecified engagement-to-gain specifications (biphasic, monotone-linear and saturating shapes at several coefficient scales). `n_undecidable` counts specifications under which the readout does not exist (for example a ratio whose denominator is silenced); the headline fraction counts those against the conclusion. The two topology rows are decided on 100 permutations per specification (resolution 0.0099) and on Benjamini-Hochberg adjusted probabilities across the 100 structural tests of the run, with the uncorrected count beside them. Until v0.6.1 the specification context rebound the gain function on the assay modules but not on the engine the permutation path resolves, so every specification fed the *default* gains to its nulls and the matrix's topology rows were empty; they are now computed per specification. For a conclusion shaped as a failure to reject, `n_equivalent_within_tolerance` is the part of the retention that is evidence of equivalence and `n_indeterminate` the part that is only a non-rejection.

| conclusion | n_specs | n_retained | n_lost | n_undecidable | n_retained_uncorrected | n_equivalent_within_tolerance | n_indeterminate |
|---|---|---|---|---|---|---|---|
| C1_nicotinic_suppression | 25 | 15 | 10 | 0 | -- | -- | -- |
| C5_nicotinic_buffering | 25 | 18 | 0 | 7 | -- | -- | -- |
| C6_nav_ache_amplification | 25 | 19 | 6 | 0 | -- | -- | -- |
| C2_rdl_disinhibition | 25 | 25 | 0 | 0 | -- | -- | -- |
| C3_imidacloprid_topology_not_distinguishable | 25 | 25 | 0 | 0 | 25 | 17 | 8 |
| C4_fipronil_topology_exceeds | 25 | 25 | 0 | 0 | 25 | 0 | 0 |
| C7_map_bitter_veto | 25 | 25 | 0 | 0 | -- | -- | -- |
