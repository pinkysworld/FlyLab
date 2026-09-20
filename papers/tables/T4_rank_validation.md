**T4_rank_validation.** Literature concordance of the model's ordering with published orderings in `data/literature/published_rank_orders.yaml`. `independence` is computed by comparing DOIs/PMIDs: `shared-source` means the published ordering comes from a publication the library already cites for one of the compounds it orders, so agreement there is internal consistency and not out-of-sample validation. `skipped` entries are those the library cannot cover; they are reported, not dropped.

| assay | id | species | independence | n_compounds | skipped | spearman_rho | exact_match |
|---|---|---|---|---|---|---|---|
| occupancy | neonic_potency_dmel_larval | Drosophila melanogaster | source-disjoint | 2 | False | 1.000 | True |
| occupancy | neonic_potency_dmel_adult | Drosophila melanogaster | source-disjoint | 2 | False | 1.000 | True |
| occupancy | neonic_chronic_lifespan_dmel | Drosophila melanogaster | source-disjoint | 4 | True | -- | -- |
| occupancy | neonic_resistance_ratio_myzus | Myzus persicae | shared-source (not out-of-sample) | 2 | False | 1.000 | True |
| occupancy | neonic_insect_vs_vertebrate_nachr | insect (Drosophila, Musca, Myzus) vs vertebrate (rat/chicken alpha4beta2, alpha7) | source-disjoint | 5 | False | 0.707 | False |
| occupancy | neonic_alpha7_partial_agonism | Rattus norvegicus (recombinant alpha7) | source-disjoint | 4 | True | -- | -- |
| occupancy | gaba_blocker_potency_insect_vs_vertebrate | Musca domestica (insect) vs Homo sapiens recombinant GABA-A | source-disjoint | 3 | True | -- | -- |
| occupancy | fipronil_human_beta3_homomer | Homo sapiens (recombinant) | source-disjoint | 3 | True | -- | -- |
| occupancy | fipronil_phenylpyrazole_selectivity_direction | Musca domestica vs Mus musculus | source-disjoint | 1 | True | -- | -- |
| occupancy | rdl_antagonist_order_drosophila | Drosophila melanogaster | source-disjoint | 3 | True | -- | -- |
| occupancy | benzodiazepine_inactive_at_insect_rdl | Drosophila melanogaster | source-disjoint | 2 | True | -- | -- |
| occupancy | ivermectin_gluCl_vs_vertebrate | Haemonchus contortus GluCl vs mammalian recombinant GlyR alpha1 / GABA-A | source-disjoint | 3 | False | 1.000 | True |
| subgraph | neonic_potency_dmel_larval | Drosophila melanogaster | source-disjoint | 2 | True | -- | -- |
| subgraph | neonic_potency_dmel_adult | Drosophila melanogaster | source-disjoint | 2 | True | -- | -- |
| subgraph | neonic_chronic_lifespan_dmel | Drosophila melanogaster | source-disjoint | 4 | True | -- | -- |
| subgraph | neonic_resistance_ratio_myzus | Myzus persicae | shared-source (not out-of-sample) | 2 | False | 1.000 | True |
| subgraph | neonic_insect_vs_vertebrate_nachr | insect (Drosophila, Musca, Myzus) vs vertebrate (rat/chicken alpha4beta2, alpha7) | source-disjoint | 5 | False | 1.000 | False |
| subgraph | neonic_alpha7_partial_agonism | Rattus norvegicus (recombinant alpha7) | source-disjoint | 4 | False | 0.316 | False |
| subgraph | gaba_blocker_potency_insect_vs_vertebrate | Musca domestica (insect) vs Homo sapiens recombinant GABA-A | source-disjoint | 3 | True | -- | -- |
| subgraph | fipronil_human_beta3_homomer | Homo sapiens (recombinant) | source-disjoint | 3 | True | -- | -- |
| subgraph | fipronil_phenylpyrazole_selectivity_direction | Musca domestica vs Mus musculus | source-disjoint | 1 | True | -- | -- |
| subgraph | rdl_antagonist_order_drosophila | Drosophila melanogaster | source-disjoint | 3 | True | -- | -- |
| subgraph | benzodiazepine_inactive_at_insect_rdl | Drosophila melanogaster | source-disjoint | 2 | True | -- | -- |
| subgraph | ivermectin_gluCl_vs_vertebrate | Haemonchus contortus GluCl vs mammalian recombinant GlyR alpha1 / GABA-A | source-disjoint | 3 | True | -- | -- |
