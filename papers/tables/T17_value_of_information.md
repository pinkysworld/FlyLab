**T17_value_of_information.** Value of information. VOI_j = S_j x Var(Y) is the model variance that would disappear if assumption j were resolved exactly while everything else stayed as uncertain as it is; the upper bound uses the total-order index and is what resolving j *last* would buy. A value of information cannot be negative, so a negative first-order estimate is clipped to zero for the decision value and kept unclipped in `voi_fraction_raw`. `state` says whether the sample could resolve the factor at all: an unresolved factor carries no ranking claim, which is a statement about the sample size rather than about the factor. This is variance of a model output under assumed input ranges, not an expected gain in accuracy about a living fly. Rate-engine quantities labelled Hz in legacy fields or axes are rate-model units, not calibrated physiological firing rates.

| rank | factor | state | voi_fraction_of_var | voi_var_hz2 | experiment | cost |
|---|---|---|---|---|---|---|
| 1 | gain_transform | resolved | 0.410 | 1.290 | calibration of receptor engagement against synaptic gain | high |
| 2 | weight_threshold | resolved | 0.294 | 0.925 | no experiment: a reconstruction-confidence analysis | none (compute only) |
| 3 | drive | null control | 0.019 | 0.058 | in vivo baseline firing rates of the driven cell classes | medium |
| 4 | transmitter | null control | 0.006 | 0.019 | immunostaining or ground-truth transmitter labels | medium |
| 5 | expression | null control | 0.000 | 0.000 | FISH / scRNA-seq for receptor expression in the named cell class | medium |
| 6 | potency | null control | 0.000 | 0.000 | electrophysiological dose-response on the stated receptor | medium |
| 7 | lif_seed | null control | 0.000 | 0.000 | none: increase the number of simulated replicates | none (compute only) |
| 8 | hill_n | null control | 0.000 | 0.000 | the same dose-response, read for slope rather than midpoint | low (rides along with the potency ladder) |
| 9 | gain_coef | null control | 0.000 | 0.000 | the same calibration, read for amplitude rather than shape | medium |
