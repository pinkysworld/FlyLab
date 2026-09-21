**T20_transformation_rule.** The validated evidence dispatch rule, materialised. The permitted transformation is a function of two facts: the parameter the source measured and the ordered distance E0-E4 between that source and this compound at this receptor in this species. A binding constant yields physical fractional occupancy only at E0; transferred, it becomes a labelled proxy carrying a warning that names the gap. Anything at E4 is `not_modelled` and returns N/A. The finite tests in `tests/test_pharm_evidence.py` enumerate every row in the shipped library and each listed engagement-to-circuit path: no number reaching a circuit readout originates from a row whose (param_type, distance) pair does not admit one. This is an audit of that finite implementation surface, not a proof over arbitrary future rows or paths.

| param_type | evidence_distance_label | engagement_model | model_strength |
|---|---|---|---|
| Kd | E0 (on-target) | binding_occupancy | 4 |
| Kd | E1 (cross-species) | binding_engagement_proxy | 3 |
| Kd | E2 (related receptor) | binding_engagement_proxy | 3 |
| Kd | E3 (class extrapolation) | functional_engagement_proxy | 1 |
| Kd | E4 (unsupported) | not_modelled | 0 |
| Ki | E0 (on-target) | binding_occupancy | 4 |
| Ki | E1 (cross-species) | binding_engagement_proxy | 3 |
| Ki | E2 (related receptor) | binding_engagement_proxy | 3 |
| Ki | E3 (class extrapolation) | functional_engagement_proxy | 1 |
| Ki | E4 (unsupported) | not_modelled | 0 |
| EC50 | E0 (on-target) | functional_engagement | 2 |
| EC50 | E1 (cross-species) | functional_engagement_proxy | 1 |
| EC50 | E2 (related receptor) | functional_engagement_proxy | 1 |
| EC50 | E3 (class extrapolation) | functional_engagement_proxy | 1 |
| EC50 | E4 (unsupported) | not_modelled | 0 |
| IC50 | E0 (on-target) | functional_engagement | 2 |
| IC50 | E1 (cross-species) | functional_engagement_proxy | 1 |
| IC50 | E2 (related receptor) | functional_engagement_proxy | 1 |
| IC50 | E3 (class extrapolation) | functional_engagement_proxy | 1 |
| IC50 | E4 (unsupported) | not_modelled | 0 |
| Kb | E0 (on-target) | functional_engagement | 2 |
| Kb | E1 (cross-species) | functional_engagement_proxy | 1 |
| Kb | E2 (related receptor) | functional_engagement_proxy | 1 |
| Kb | E3 (class extrapolation) | functional_engagement_proxy | 1 |
| Kb | E4 (unsupported) | not_modelled | 0 |
| relative_potency | E0 (on-target) | not_modelled | 0 |
| relative_potency | E1 (cross-species) | not_modelled | 0 |
| relative_potency | E2 (related receptor) | not_modelled | 0 |
| relative_potency | E3 (class extrapolation) | not_modelled | 0 |
| relative_potency | E4 (unsupported) | not_modelled | 0 |
| class_order | E0 (on-target) | not_modelled | 0 |
| class_order | E1 (cross-species) | not_modelled | 0 |
| class_order | E2 (related receptor) | not_modelled | 0 |
| class_order | E3 (class extrapolation) | not_modelled | 0 |
| class_order | E4 (unsupported) | not_modelled | 0 |
| unknown | E0 (on-target) | not_modelled | 0 |
| unknown | E1 (cross-species) | not_modelled | 0 |
| unknown | E2 (related receptor) | not_modelled | 0 |
| unknown | E3 (class extrapolation) | not_modelled | 0 |
| unknown | E4 (unsupported) | not_modelled | 0 |
