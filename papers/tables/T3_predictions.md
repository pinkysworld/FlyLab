**T3_predictions.** Pre-registered predictions H1-H7. `predicted_effect` and its CI are model-internal (teaching-EC50 jitter, drive jitter, RNG seed) and carry no biological variance; `suggested_n_per_group` therefore caps the standardised effect at d = 1.0 before the power calculation. `live_result` is null for every row and stays null until a real table is imported by hand.

| id | assay | compound | readout | predicted_direction | predicted_effect | suggested_n_per_group | status |
|---|---|---|---|---|---|---|---|
| H1 | occupancy | imidacloprid | occupancy_gap_insect_minus_vertebrate | increase | 0.854 | 16 | software_prediction |
| H2 | subgraph | imidacloprid | mean_hz | decrease | -6.145 | 23 | software_prediction |
| H3 | subgraph | fipronil | mean_hz | increase | 0.896 | 16 | software_prediction |
| H4 | subgraph | diazepam | mean_hz | none | 0.000 | -- | software_prediction |
| H5 | taste_map | -- | mn9_sugar_bitter_hz - mn9_sugar_hz | decrease | -0.296 | 23 | software_prediction |
| H6 | taste_map | fipronil | bitter_veto_ratio | increase | 0.413 | 16 | software_prediction |
| H7 | subgraph | picrotoxin | mean_hz | increase | 0.561 | 16 | software_prediction |
