**T23_transmitter_null_balance.** Share of total outgoing synaptic weight carried by each transmitter on the `named` cut, on the real graph and across the draws of each transmitter null. A gain patch acts *through* transmitter identity, so this share is what the drug sees. Plain label permutation moves it -- the real graph sits outside its own null on the cholinergic share -- which makes that mode a joint target-set-and-sign null rather than a null about transmitter identity. The weight-matched null is a constrained shuffle that holds every tracked share inside the tolerance band and is the ladder's rank-3 rung.

| mode | transmitter | real | null_mean | null_sd | percentile_of_real | within_tol |
|---|---|---|---|---|---|---|
| sign_permute | acetylcholine | 0.619 | 0.544 | 0.027 | 100 | False |
| sign_permute | gaba | 0.301 | 0.262 | 0.023 | 94.9 | False |
| sign_permute | glutamate | 0.044 | 0.050 | 0.011 | 32.5 | False |
| sign_permute_weight_matched | acetylcholine | 0.619 | 0.617 | 3.87e-04 | 99.4 | True |
| sign_permute_weight_matched | gaba | 0.301 | 0.299 | 4.72e-04 | 98.3 | True |
| sign_permute_weight_matched | glutamate | 0.044 | 0.043 | 6.54e-04 | 94.1 | True |
