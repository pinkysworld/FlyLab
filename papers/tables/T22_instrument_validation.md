**T22_instrument_validation.** Ground-truth recovery and power for the dependence ladder. A recurrent cholinergic cycle of known strength is planted in a synthetic cut of the same size, density and transmitter composition as the `named` cut, and a gain patch that collapses `g_ach` removes an amplification that exists only while the cycle is intact -- so the drug *contrast*, not merely the rate, depends on the wiring. `loop_strength = 0` plants nothing and is the negative control, whose detection rate is the empirical false-positive rate. Nothing here touches the connectome or the compound library: the experiment is about the instrument.

| experiment | loop_strength | n | replicates | detection_rate | class |
|---|---|---|---|---|---|
| recovery | 0.000 | 200 | 1 | 0.000 | composition-dominated |
| recovery | 0.500 | 200 | 1 | 0.000 | composition-dominated |
| recovery | 1.000 | 200 | 1 | 1.000 | topology-dependent |
| recovery | 2.000 | 200 | 1 | 1.000 | topology-dependent |
| power | 0.000 | 50 | 6 | 0.000 | -- |
| power | 0.000 | 200 | 6 | 0.000 | -- |
| power | 0.000 | 1000 | 6 | 0.000 | -- |
| power | 0.250 | 50 | 6 | 0.167 | -- |
| power | 0.250 | 200 | 6 | 0.000 | -- |
| power | 0.250 | 1000 | 6 | 0.167 | -- |
| power | 0.500 | 50 | 6 | 0.667 | -- |
| power | 0.500 | 200 | 6 | 0.833 | -- |
| power | 0.500 | 1000 | 6 | 0.833 | -- |
| power | 1.000 | 50 | 6 | 1.000 | -- |
| power | 1.000 | 200 | 6 | 1.000 | -- |
| power | 1.000 | 1000 | 6 | 1.000 | -- |
