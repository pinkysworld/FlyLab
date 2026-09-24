**T24_composition_reference.** Reference distribution for the composition-versus-full rank correlation. The composition level and the full model are not independent models: both are functions of the same gain vector, so the correlation has a large structural floor and the observed value must be compared with a matched reference. The latter draws pseudo-compounds with the shape the shipped mechanism rules produce (one receptor, one transmitter, one gain moved) and no pharmacology at all; shuffling the compound labels leaves the observed value unchanged, which is a consequence of moving the paired B and D values together.

| condition | spearman_rho | n_compounds | draws |
|---|---|---|---|
| observed (library gain vectors) | 0.921 | 21 | -- |
| library, floor-saturated compounds dropped | 0.826 | 16 | -- |
| matched reference: one gain moved per pseudo-compound | 0.847 | 21 | 50 |
| dense reference: every gain moved per pseudo-compound | 0.808 | 21 | 50 |
| compound labels shuffled | 0.921 | 21 | 1 |
