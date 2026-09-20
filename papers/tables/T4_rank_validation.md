**T4_rank_validation.** Retrospective rank validation against `data/literature/published_rank_orders.yaml`. `skipped` entries are those the teaching library cannot cover (missing compounds or placeholder-only rows); they are reported, not dropped.

| assay | id | species | n_compounds | skipped | spearman_rho | exact_match |
|---|---|---|---|---|---|---|
| occupancy | neonic_potency_dmel_larval | Drosophila melanogaster | 2 | False | 1.000 | True |
| occupancy | neonic_potency_dmel_adult | Drosophila melanogaster | 2 | False | 1.000 | True |
| occupancy | neonic_chronic_lifespan_dmel | Drosophila melanogaster | 4 | True | -- | -- |
| occupancy | neonic_resistance_ratio_myzus | Myzus persicae | 2 | False | 1.000 | True |
| occupancy | neonic_insect_vs_vertebrate_nachr | insect (Drosophila, Musca, Myzus) vs vertebrate (rat/chicken alpha4beta2, alpha7) | 5 | False | 0.707 | False |
| occupancy | neonic_alpha7_partial_agonism | Rattus norvegicus (recombinant alpha7) | 4 | False | 0.400 | False |
| occupancy | gaba_blocker_potency_insect_vs_vertebrate | Musca domestica (insect) vs Homo sapiens recombinant GABA-A | 3 | True | -- | -- |
| occupancy | fipronil_human_beta3_homomer | Homo sapiens (recombinant) | 3 | True | -- | -- |
| occupancy | fipronil_phenylpyrazole_selectivity_direction | Musca domestica vs Mus musculus | 1 | True | -- | -- |
| occupancy | rdl_antagonist_order_drosophila | Drosophila melanogaster | 3 | True | -- | -- |
| occupancy | benzodiazepine_inactive_at_insect_rdl | Drosophila melanogaster | 2 | True | -- | -- |
| occupancy | ivermectin_gluCl_vs_vertebrate | Haemonchus contortus GluCl vs mammalian recombinant GlyR alpha1 / GABA-A | 3 | False | 1.000 | True |
| subgraph | neonic_potency_dmel_larval | Drosophila melanogaster | 2 | True | -- | -- |
| subgraph | neonic_potency_dmel_adult | Drosophila melanogaster | 2 | True | -- | -- |
| subgraph | neonic_chronic_lifespan_dmel | Drosophila melanogaster | 4 | True | -- | -- |
| subgraph | neonic_resistance_ratio_myzus | Myzus persicae | 2 | False | 1.000 | True |
| subgraph | neonic_insect_vs_vertebrate_nachr | insect (Drosophila, Musca, Myzus) vs vertebrate (rat/chicken alpha4beta2, alpha7) | 5 | False | 1.000 | False |
| subgraph | neonic_alpha7_partial_agonism | Rattus norvegicus (recombinant alpha7) | 4 | False | 0.316 | False |
| subgraph | gaba_blocker_potency_insect_vs_vertebrate | Musca domestica (insect) vs Homo sapiens recombinant GABA-A | 3 | True | -- | -- |
| subgraph | fipronil_human_beta3_homomer | Homo sapiens (recombinant) | 3 | True | -- | -- |
| subgraph | fipronil_phenylpyrazole_selectivity_direction | Musca domestica vs Mus musculus | 1 | True | -- | -- |
| subgraph | rdl_antagonist_order_drosophila | Drosophila melanogaster | 3 | True | -- | -- |
| subgraph | benzodiazepine_inactive_at_insect_rdl | Drosophila melanogaster | 2 | True | -- | -- |
| subgraph | ivermectin_gluCl_vs_vertebrate | Haemonchus contortus GluCl vs mammalian recombinant GlyR alpha1 / GABA-A | 3 | True | -- | -- |
