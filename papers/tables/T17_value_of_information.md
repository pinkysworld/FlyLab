**T17_value_of_information.** Value of information. VOI_j = S_j x Var(Y) is the model variance that would disappear if assumption j were resolved exactly while everything else stayed as uncertain as it is; the upper bound uses the total-order index and is what resolving j *last* would buy. This is variance of a model output under assumed input ranges, not an expected gain in accuracy about a living fly.

| rank | factor | voi_fraction_of_var | voi_var_hz2 | experiment | cost |
|---|---|---|---|---|---|
| 1 | gain_transform | 0.410 | 1.290 | calibration of receptor engagement against synaptic gain | high |
| 2 | weight_threshold | 0.294 | 0.925 | no experiment: a reconstruction-confidence analysis | none (compute only) |
| 3 | drive | 0.019 | 0.058 | in vivo baseline firing rates of the driven cell classes | medium |
| 4 | transmitter | 0.006 | 0.019 | immunostaining or ground-truth transmitter labels | medium |
| 5 | expression | -0.006 | -0.018 | FISH / scRNA-seq for receptor expression in the named cell class | medium |
| 6 | potency | -0.006 | -0.019 | electrophysiological dose-response on the stated receptor | medium |
| 7 | lif_seed | -0.007 | -0.021 | none: increase the number of simulated replicates | none (compute only) |
| 8 | hill_n | -0.008 | -0.025 | the same dose-response, read for slope rather than midpoint | low (rides along with the potency ladder) |
| 9 | gain_coef | -0.046 | -0.144 | the same calibration, read for amplitude rather than shape | medium |
