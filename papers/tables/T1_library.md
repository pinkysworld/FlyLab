**T1_library.** The sourced compound library (schema v3). Every row carries an evidence tier, a named source, the KIND of parameter the source reported (`param_type`: Kd/Ki are binding constants, EC50/IC50/Kb are functional potencies) and how far that source sits from this compound/receptor/species (`relation`); `class_placeholder` rows carry no number at all and report N/A rather than a small response.

| compound | class | receptor | param_type | value_M | direction | evidence_tier |
|---|---|---|---|---|---|---|
| acetamiprid | neonicotinoid | insect_nAChR | EC50 | 5.00e-08 | agonist | literature_order |
| acetamiprid | neonicotinoid | vertebrate_nAChR_a7 | EC50 | 7.30e-04 | partial_agonist | literature_order |
| acetamiprid | neonicotinoid | vertebrate_nAChR_a4b2 | EC50 | 1.00e-04 | agonist | literature_order |
| acetamiprid | neonicotinoid | insect_RDL | unknown | -- | none | class_placeholder |
| acetamiprid | neonicotinoid | vertebrate_GABA_A | unknown | -- | none | class_placeholder |
| acetylcholine | endogenous transmitter | insect_nAChR | EC50 | 1.00e-05 | agonist | literature_order |
| acetylcholine | endogenous transmitter | vertebrate_nAChR_a4b2 | EC50 | 1.00e-05 | agonist | literature_order |
| acetylcholine | endogenous transmitter | vertebrate_nAChR_a7 | EC50 | 1.60e-04 | agonist | literature_order |
| acetylcholine | endogenous transmitter | insect_AChE | unknown | -- | none | class_placeholder |
| acetylcholine | endogenous transmitter | insect_RDL | unknown | -- | none | class_placeholder |
| acetylcholine | endogenous transmitter | vertebrate_GABA_A | unknown | -- | none | class_placeholder |
| caffeine | methylxanthine | vertebrate_GABA_A | IC50 | 0.001 | antagonist | literature_order |
| caffeine | methylxanthine | insect_nAChR | unknown | -- | none | class_placeholder |
| caffeine | methylxanthine | insect_RDL | unknown | -- | none | class_placeholder |
| caffeine | methylxanthine | vertebrate_nAChR_a4b2 | unknown | -- | none | class_placeholder |
| chlordimeform | formamidine | insect_OctR | EC50 | 3.00e-05 | partial_agonist | literature_order |
| chlordimeform | formamidine | insect_nAChR | unknown | -- | none | class_placeholder |
| chlordimeform | formamidine | insect_RDL | unknown | -- | none | class_placeholder |
| chlordimeform | formamidine | vertebrate_nAChR_a4b2 | unknown | -- | none | class_placeholder |
| chlordimeform | formamidine | vertebrate_GABA_A | unknown | -- | none | class_placeholder |
| chlorpyrifos_oxon | organophosphate (oxon metabolite) | insect_AChE | unknown | -- | none | class_placeholder |
| chlorpyrifos_oxon | organophosphate (oxon metabolite) | vertebrate_AChE | IC50 | 3.00e-07 | inhibitor | literature_order |
| chlorpyrifos_oxon | organophosphate (oxon metabolite) | insect_nAChR | unknown | -- | none | class_placeholder |
| chlorpyrifos_oxon | organophosphate (oxon metabolite) | vertebrate_nAChR_a4b2 | unknown | -- | none | class_placeholder |
| chlorpyrifos_oxon | organophosphate (oxon metabolite) | insect_RDL | unknown | -- | none | class_placeholder |
| chlorpyrifos_oxon | organophosphate (oxon metabolite) | vertebrate_GABA_A | unknown | -- | none | class_placeholder |
| clothianidin | neonicotinoid | insect_nAChR | EC50 | 3.00e-08 | agonist | literature_order |
| clothianidin | neonicotinoid | vertebrate_nAChR_a7 | EC50 | 7.40e-04 | partial_agonist | literature_order |
| clothianidin | neonicotinoid | vertebrate_nAChR_a4b2 | EC50 | 1.00e-04 | agonist | literature_order |
| clothianidin | neonicotinoid | insect_RDL | unknown | -- | none | class_placeholder |
| clothianidin | neonicotinoid | vertebrate_GABA_A | unknown | -- | none | class_placeholder |
| ddt | organochlorine | insect_Nav | EC50 | 6.50e-05 | positive_modulator | literature_order |
| ddt | organochlorine | vertebrate_Nav1_x | EC50 | 0.001 | positive_modulator | literature_order |
| ddt | organochlorine | insect_nAChR | unknown | -- | none | class_placeholder |
| ddt | organochlorine | insect_RDL | unknown | -- | none | class_placeholder |
| ddt | organochlorine | vertebrate_GABA_A | unknown | -- | none | class_placeholder |
| deltamethrin | type II pyrethroid | insect_Nav | EC50 | 4.30e-08 | positive_modulator | literature_order |
| deltamethrin | type II pyrethroid | vertebrate_Nav1_x | EC50 | 5.00e-06 | positive_modulator | literature_order |
| deltamethrin | type II pyrethroid | insect_nAChR | unknown | -- | none | class_placeholder |
| deltamethrin | type II pyrethroid | insect_RDL | unknown | -- | none | class_placeholder |
| deltamethrin | type II pyrethroid | vertebrate_GABA_A | unknown | -- | none | class_placeholder |
| diazepam | benzodiazepine | insect_nAChR | unknown | -- | none | class_placeholder |
| diazepam | benzodiazepine | vertebrate_nAChR_a4b2 | unknown | -- | none | class_placeholder |
| diazepam | benzodiazepine | vertebrate_GABA_A | EC50 | 1.00e-08 | positive_modulator | literature_order |
| diazepam | benzodiazepine | insect_RDL | unknown | -- | none | class_placeholder |
| dieldrin | cyclodiene organochlorine | insect_RDL | IC50 | 1.00e-07 | antagonist | literature_order |
| dieldrin | cyclodiene organochlorine | vertebrate_GABA_A | IC50 | 1.00e-06 | antagonist | literature_order |
| dieldrin | cyclodiene organochlorine | insect_nAChR | unknown | -- | none | class_placeholder |
| dieldrin | cyclodiene organochlorine | vertebrate_nAChR_a4b2 | unknown | -- | none | class_placeholder |
| fipronil | phenylpyrazole | insect_nAChR | unknown | -- | none | class_placeholder |
| fipronil | phenylpyrazole | vertebrate_nAChR_a4b2 | unknown | -- | none | class_placeholder |
| fipronil | phenylpyrazole | vertebrate_GABA_A | IC50 | 1.10e-06 | antagonist | literature_order |
| fipronil | phenylpyrazole | insect_RDL | IC50 | 3.00e-08 | antagonist | literature_order |
| gaba | endogenous transmitter | insect_RDL | EC50 | 5.00e-05 | agonist | literature_order |
| gaba | endogenous transmitter | vertebrate_GABA_A | EC50 | 1.00e-05 | agonist | literature_order |
| gaba | endogenous transmitter | insect_nAChR | unknown | -- | none | class_placeholder |
| gaba | endogenous transmitter | vertebrate_nAChR_a4b2 | unknown | -- | none | class_placeholder |
| imidacloprid | neonicotinoid | insect_nAChR | EC50 | 2.00e-08 | agonist | literature_order |
| imidacloprid | neonicotinoid | vertebrate_nAChR_a4b2 | EC50 | 1.00e-05 | agonist | literature_order |
| imidacloprid | neonicotinoid | vertebrate_GABA_A | unknown | -- | none | class_placeholder |
| imidacloprid | neonicotinoid | insect_RDL | unknown | -- | none | class_placeholder |
| imidacloprid | neonicotinoid | insect_nAChR_native_dmel | Kd | 2.00e-09 | agonist | literature_order |
| imidacloprid | neonicotinoid | insect_nAChR_beta1 | Kd | 8.30e-11 | agonist | literature_order |
| ivermectin | avermectin | insect_GluCl | EC50 | 2.00e-08 | agonist | literature_order |
| ivermectin | avermectin | insect_RDL | unknown | -- | none | class_placeholder |
| ivermectin | avermectin | vertebrate_GlyR | EC50 | 1.00e-06 | positive_modulator | literature_order |
| ivermectin | avermectin | vertebrate_GABA_A | EC50 | 1.00e-05 | positive_modulator | literature_order |
| ivermectin | avermectin | insect_nAChR | unknown | -- | none | class_placeholder |
| nicotine | alkaloid | insect_nAChR | EC50 | 2.00e-07 | agonist | literature_order |
| nicotine | alkaloid | vertebrate_nAChR_a4b2 | EC50 | 1.00e-06 | agonist | literature_order |
| nicotine | alkaloid | vertebrate_GABA_A | unknown | -- | none | class_placeholder |
| nicotine | alkaloid | insect_RDL | unknown | -- | none | class_placeholder |
| nitenpyram | neonicotinoid | insect_nAChR | EC50 | 5.00e-08 | agonist | literature_order |
| nitenpyram | neonicotinoid | vertebrate_nAChR_a4b2 | EC50 | 2.00e-05 | agonist | literature_order |
| nitenpyram | neonicotinoid | vertebrate_GABA_A | unknown | -- | none | class_placeholder |
| nitenpyram | neonicotinoid | insect_RDL | unknown | -- | none | class_placeholder |
| permethrin | type I pyrethroid | insect_Nav | EC50 | 4.00e-07 | positive_modulator | literature_order |
| permethrin | type I pyrethroid | vertebrate_Nav1_x | unknown | -- | none | class_placeholder |
| permethrin | type I pyrethroid | insect_nAChR | unknown | -- | none | class_placeholder |
| permethrin | type I pyrethroid | insect_RDL | unknown | -- | none | class_placeholder |
| permethrin | type I pyrethroid | vertebrate_GABA_A | unknown | -- | none | class_placeholder |
| picrotoxin | plant convulsant | insect_RDL | IC50 | 1.00e-06 | antagonist | literature_order |
| picrotoxin | plant convulsant | vertebrate_GABA_A | IC50 | 1.00e-06 | antagonist | literature_order |
| picrotoxin | plant convulsant | vertebrate_GlyR | unknown | -- | none | class_placeholder |
| picrotoxin | plant convulsant | insect_nAChR | unknown | -- | none | class_placeholder |
| picrotoxin | plant convulsant | vertebrate_nAChR_a4b2 | unknown | -- | none | class_placeholder |
| spinosad | spinosyn | insect_nAChR | unknown | -- | none | class_placeholder |
| spinosad | spinosyn | vertebrate_nAChR_a4b2 | unknown | -- | none | class_placeholder |
| spinosad | spinosyn | vertebrate_nAChR_a7 | unknown | -- | none | class_placeholder |
| spinosad | spinosyn | insect_RDL | unknown | -- | none | class_placeholder |
| spinosad | spinosyn | vertebrate_GABA_A | unknown | -- | none | class_placeholder |
| spinosad | spinosyn | insect_nAChR_alpha6 | unknown | -- | none | class_placeholder |
| sulfoxaflor | sulfoximine | insect_nAChR | unknown | -- | none | class_placeholder |
| sulfoxaflor | sulfoximine | vertebrate_nAChR_a4b2 | EC50 | 0.001 | agonist | literature_order |
| sulfoxaflor | sulfoximine | vertebrate_nAChR_a7 | unknown | -- | none | class_placeholder |
| sulfoxaflor | sulfoximine | insect_RDL | unknown | -- | none | class_placeholder |
| sulfoxaflor | sulfoximine | vertebrate_GABA_A | unknown | -- | none | class_placeholder |
| thiamethoxam | neonicotinoid | insect_nAChR | EC50 | 1.00e-04 | partial_agonist | literature_order |
| thiamethoxam | neonicotinoid | vertebrate_nAChR_a7 | unknown | -- | none | class_placeholder |
| thiamethoxam | neonicotinoid | vertebrate_nAChR_a4b2 | unknown | -- | none | class_placeholder |
| thiamethoxam | neonicotinoid | insect_RDL | unknown | -- | none | class_placeholder |
| thiamethoxam | neonicotinoid | vertebrate_GABA_A | unknown | -- | none | class_placeholder |
