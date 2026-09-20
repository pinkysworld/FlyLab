**T19_cut_census.** Structure of the two committed cuts. `share_onto_seeds` is the fraction of edges terminating on a seed cell and `share_edges_from_nodes_with_input` the fraction whose source itself receives input -- the recurrence budget, and the only edges that can carry a path longer than one hop. A cut where most edges point at a few hubs and most nodes have no input is an in-star, and a degree-preserving rewire of it destroys very little; the transmitter census is by cell, not by synapse, and a transmitter absent from a cut cannot be perturbed on it.

| cut | n_nodes | n_edges | mean_degree | share_onto_seeds | nodes_with_in_degree | share_edges_from_nodes_with_input | in_star |
|---|---|---|---|---|---|---|---|
| named | 1126 | 1360 | 1.208 | 0.849 | 184 | 0.249 | True |
| taste_motor | 1841 | 19066 | 10.4 | 0.096 | 1643 | 0.973 | False |
