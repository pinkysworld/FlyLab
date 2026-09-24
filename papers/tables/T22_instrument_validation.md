**T22_instrument_validation.** Limited planted-cycle recovery check for the dependence ladder. A recurrent cholinergic cycle of known strength is planted in a synthetic 250-node cut with about 1500 background edges. Its transmitter mix was chosen close to the `named` cut, but size and connectivity were not matched. A gain patch that collapses `g_ach` removes an amplification that exists only while the cycle is intact, so the drug *contrast* depends on the wiring by construction. `loop_strength = 0` plants nothing and is the negative control. The small sample gives an observed control count, not a precise false-positive rate or power estimate for the paper's endpoint. Nothing here touches the connectome or compound library.

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
