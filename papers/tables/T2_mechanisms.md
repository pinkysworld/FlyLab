**T2_mechanisms.** Mechanism to gain-patch rules. `th` is fractional occupancy at the named receptor; every gain is a dimensionless multiplier with 1.0 = vehicle and a floor of 0.05. This table is generated from `flylab.pharm.mechanisms.MECHANISM_TABLE`, the single source of truth every assay calls.

| receptor | direction | gain | formula |
|---|---|---|---|
| insect_nAChR | agonist | g_ach | max(0.05, 1 + 0.4*th - 1.6*th^2) |
| insect_nAChR | partial_agonist | g_ach | max(0.05, 1 + 0.4*th - 1.6*th^2) |
| insect_nAChR | positive_modulator | g_ach | max(0.05, 1 + 0.4*th - 1.6*th^2) |
| insect_nAChR | antagonist | g_ach | max(0.05, 1 - th) |
| insect_RDL | antagonist | g_gaba | max(0.05, 1 - th) |
| insect_RDL | agonist | g_gaba | max(0.05, 1 + 0.4*th) |
| insect_RDL | positive_modulator | g_gaba | max(0.05, 1 + 0.4*th) |
| insect_GluCl | agonist | g_glu | max(0.05, 1 + 0.8*th) |
| insect_GluCl | antagonist | g_glu | max(0.05, 1 - th) |
| insect_AChE | inhibitor | ach_tone | 1 + 2.0*th |
| insect_AChE | inhibitor | g_ach | g_ach * ach_tone * (1 + 0.4*th - 1.6*th^2) |
| insect_Nav | positive_modulator | g_nav | 1 + 1.5*th |
| insect_Nav | antagonist | g_nav | max(0.05, 1 - th) |
| insect_OctR | agonist | g_oct | 1 + 0.5*th |
| insect_OctR | antagonist | g_oct | max(0.05, 1 - th) |
