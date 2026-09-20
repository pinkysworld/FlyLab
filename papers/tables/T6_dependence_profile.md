**T6_dependence_profile.** Connectome-dependence profile. `real_effect` is treated minus vehicle on the real MaleCNS cut; the null is the same contrast on `n` degraded copies of that cut, with the seed block and node order held fixed. `p_two_sided` is the empirical permutation probability (k+1)/(n+1) and is the statistic to read first; `p_resolution` is its floor, and `z` is a standardised distance from a usually non-normal null. `n_stabilised` is the smallest permutation count from which the verdict no longer moved.

| compound | assay | mode | information_kept | n | real_effect | p_two_sided | z |
|---|---|---|---|---|---|---|---|
| imidacloprid | subgraph | sign_permute | the real edges, the real weights and the NT histogram | 1000 | -6.173 | 0.275 | -1.229 |
| imidacloprid | subgraph | weight_permute | the real edge list and the transmitters | 1000 | -6.173 | 0.586 | -0.502 |
| imidacloprid | subgraph | rewire_degree_preserving | every node's in/out degree and its transmitter | 1000 | -6.173 | 0.472 | -0.720 |
| imidacloprid | subgraph | erdos_renyi | N, E, the weight histogram and the transmitter census | 1000 | -6.173 | 9.99e-04 | -47.1 |
| fipronil | subgraph | sign_permute | the real edges, the real weights and the NT histogram | 1000 | 0.900 | 0.113 | 1.680 |
| fipronil | subgraph | weight_permute | the real edge list and the transmitters | 1000 | 0.900 | 0.006 | 3.108 |
| fipronil | subgraph | rewire_degree_preserving | every node's in/out degree and its transmitter | 1000 | 0.900 | 0.007 | 2.952 |
| fipronil | subgraph | erdos_renyi | N, E, the weight histogram and the transmitter census | 1000 | 0.900 | 9.99e-04 | 657 |
| fipronil | taste_map | sign_permute | the real edges, the real weights and the NT histogram | 300 | 0.413 | 0.417 | 0.105 |
| fipronil | taste_map | weight_permute | the real edge list and the transmitters | 300 | 0.413 | 0.879 | 0.083 |
| fipronil | taste_map | rewire_degree_preserving | every node's in/out degree and its transmitter | 300 | 0.413 | 0.203 | 0.220 |
| fipronil | taste_map | erdos_renyi | N, E, the weight histogram and the transmitter census | 300 | 0.413 | 0.245 | 0.198 |
