**T26_scale_frontier.** What a dependence analysis costs at each scale, from the measured per-shuffle and per-edge constants of `flylab.analysis.scale` on one core: a profile at n = 1000 permutations, an 84-cell landscape at n = 100, and the dense matrix the LIF engine materialises. `feasible` means one result in under a day on one core, which is a generous bar for a single result and a hopeless one for a study repeated over the specification family. Node and edge counts are measured on MaleCNS v1.0.

| scale | n_nodes | n_edges | dependence_profile_h | dependence_profile_feasible | landscape_days | lif_dense_gb | lif_feasible |
|---|---|---|---|---|---|---|---|
| named | 1126 | 1360 | 0.030 | True | 0.010 | 0.010 | True |
| taste_motor | 1841 | 19066 | 0.200 | True | 0.070 | 0.010 | True |
| scale_1k | 1000 | 22857 | 0.230 | True | 0.080 | 0.000 | True |
| scale_5k | 5000 | 279845 | 2.680 | True | 0.940 | 0.100 | True |
| scale_10k | 10000 | 621600 | 5.940 | True | 2.080 | 0.400 | True |
| scale_25k | 25000 | 1364375 | 13.0 | True | 4.560 | 2.500 | True |
| scale_50k | 50000 | 2216881 | 21.1 | True | 7.400 | 10.0 | True |
| whole_cns_w5 | 165122 | 6235682 | 59.4 | False | 20.8 | 109 | False |
| whole_cns_w1 | 165122 | 25563197 | 244 | False | 85.2 | 109 | False |
