**T6_nullmodels.** Connectome null models. `real_effect` is treated minus vehicle on the real MaleCNS cut; the null distribution is the same contrast on `n` degraded copies of that cut. `p` is the empirical two-sided permutation p, (k+1)/(n+1), so its floor is 1/(n+1).

| compound | assay | readout | mode | n | real_effect | z | p_two_sided |
|---|---|---|---|---|---|---|---|
| imidacloprid | subgraph | mean_hz | sign_permute | 50 | -6.173 | -1.228 | 0.255 |
| imidacloprid | subgraph | mean_hz | weight_permute | 50 | -6.173 | -0.535 | 0.529 |
| imidacloprid | subgraph | mean_hz | rewire_degree_preserving | 50 | -6.173 | -0.794 | 0.373 |
| imidacloprid | subgraph | mean_hz | erdos_renyi | 50 | -6.173 | -45.7 | 0.020 |
| imidacloprid | taste_map | bitter_veto_ratio | sign_permute | 10 | -- | -- | -- |
| imidacloprid | taste_map | bitter_veto_ratio | weight_permute | 10 | -- | -- | -- |
| imidacloprid | taste_map | bitter_veto_ratio | rewire_degree_preserving | 10 | -- | -- | -- |
| imidacloprid | taste_map | bitter_veto_ratio | erdos_renyi | 10 | -- | -- | -- |
| fipronil | subgraph | mean_hz | sign_permute | 50 | 0.900 | 1.682 | 0.176 |
| fipronil | subgraph | mean_hz | weight_permute | 50 | 0.900 | 3.372 | 0.020 |
| fipronil | subgraph | mean_hz | rewire_degree_preserving | 50 | 0.900 | 3.016 | 0.039 |
| fipronil | subgraph | mean_hz | erdos_renyi | 50 | 0.900 | 1355 | 0.020 |
| fipronil | taste_map | bitter_veto_ratio | sign_permute | 10 | 0.413 | -0.201 | 0.636 |
| fipronil | taste_map | bitter_veto_ratio | weight_permute | 10 | 0.413 | 0.249 | 0.818 |
| fipronil | taste_map | bitter_veto_ratio | rewire_degree_preserving | 10 | 0.413 | 0.461 | 0.545 |
| fipronil | taste_map | bitter_veto_ratio | erdos_renyi | 10 | 0.413 | 0.711 | 0.455 |
