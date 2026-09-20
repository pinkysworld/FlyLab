**T27_scale_study.** The same dependence profile across the cuts the scaling study could afford, smallest first: the cut's structure, the permutation count it was run at (which falls with the edge count, so the resolution coarsens up the ladder), the per-mode empirical probabilities, the necessary information level and the class. `named` and `taste_motor` are built by a different recipe from the `scale_*` rungs, so the clean scale axis is the `scale_*` rows alone. Stability inside the range tested is not stability. The largest rung here is far short of the whole CNS, and the permutation count falls with scale, so 'settled' means 'did not move over the cuts this study could afford'.

| cut | n_nodes | n_edges | mean_degree | compound | n | p_degree | p_weight | p_sign_wm | necessary_information_level | class |
|---|---|---|---|---|---|---|---|---|---|---|
| scale_1k | 1000 | 22857 | 22.9 | fipronil | 894 | 0.001 | 0.002 | 0.666 | wiring_with_weighted_transmitter_balance | topology-dependent |
| scale_1k | 1000 | 22857 | 22.9 | imidacloprid | 894 | 0.493 | 0.007 | 0.180 | size_and_composition | topology-dependent |
| named | 1126 | 1360 | 1.208 | fipronil | 1000 | 0.007 | 0.006 | 0.210 | wiring_with_weighted_transmitter_balance | topology-dependent |
| named | 1126 | 1360 | 1.208 | imidacloprid | 1000 | 0.472 | 0.586 | 0.389 | degree_sequence | composition-dominated |
| taste_motor | 1841 | 19066 | 10.4 | fipronil | 1000 | 9.99e-04 | 9.99e-04 | 0.841 | wiring_with_weighted_transmitter_balance | topology-dependent |
| taste_motor | 1841 | 19066 | 10.4 | imidacloprid | 1000 | 9.99e-04 | 0.063 | 0.892 | topology_without_weight_pairing | topology-dependent |
| scale_5k | 5000 | 279845 | 56.0 | fipronil | 77 | 0.013 | 0.013 | 0.115 | wiring_with_weighted_transmitter_balance | topology-dependent |
| scale_5k | 5000 | 279845 | 56.0 | imidacloprid | 77 | 0.013 | 0.013 | 0.026 | real_connectome | topology-dependent |
| scale_10k | 10000 | 621600 | 62.2 | fipronil | 35 | 0.028 | 0.028 | 0.111 | wiring_with_weighted_transmitter_balance | topology-dependent |
| scale_10k | 10000 | 621600 | 62.2 | imidacloprid | 35 | 0.028 | 0.028 | 0.028 | real_connectome | topology-dependent |
| scale_25k | 25000 | 1364375 | 54.6 | fipronil | 20 | 0.048 | 0.048 | 0.143 | wiring_with_weighted_transmitter_balance | topology-dependent |
| scale_25k | 25000 | 1364375 | 54.6 | imidacloprid | 20 | 0.048 | 0.048 | 0.143 | wiring_with_weighted_transmitter_balance | topology-dependent |
| scale_50k | 50000 | 2216881 | 44.3 | fipronil | 20 | 0.048 | 0.048 | 0.095 | wiring_with_weighted_transmitter_balance | topology-dependent |
| scale_50k | 50000 | 2216881 | 44.3 | imidacloprid | 20 | 0.048 | 0.048 | 0.571 | wiring_with_weighted_transmitter_balance | topology-dependent |
