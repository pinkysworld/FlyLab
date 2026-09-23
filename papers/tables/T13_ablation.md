**T13_ablation.** Ablation ladder. Each level is a prediction for the same compound-concentration cell from a model that has been denied one layer of information: receptor engagement only (dimensionless), mechanism gains on the graph's transmitter composition only (excitation index), the real cut with one generic multiplier instead of mechanism-specific gains (Hz), and the full model (Hz). The generic multiplier is direction-aware: generic in magnitude, with one bit of sign taken from the mechanism table, which is the minimum a connectome-without-pharmacology model needs to order a library containing disinhibitors. `C_topology_only_floor` is the historical depression-only rule, kept as a floor: every one of its entries is at most zero, so it cannot express disinhibition and is not a competitive baseline. Levels have different units, so they are compared by ordering, never by value. Rate-engine quantities labelled Hz in legacy fields or axes are rate-model units, not calibrated physiological firing rates.

| compound | class | conc_M | A_receptor_only | B_composition_only | C_topology_only | D_full_flylab |
|---|---|---|---|---|---|---|
| acetamiprid | neonicotinoid | 1.00e-08 | 0.127 | 0.543 | 0.867 | 0.191 |
| acetylcholine | endogenous transmitter | 1.00e-08 | 2.51e-04 | 0.002 | 0.002 | 7.61e-04 |
| caffeine | methylxanthine | 1.00e-08 | 0.000 | 0.000 | 0.000 | 0.000 |
| chlordimeform | formamidine | 1.00e-08 | 3.33e-04 | 0.000 | 0.002 | 0.000 |
| chlorpyrifos_oxon | organophosphate (oxon metabolite) | 1.00e-08 | -0.500 | 47.8 | 3.582 | 41.9 |
| clothianidin | neonicotinoid | 1.00e-08 | 0.211 | 0.286 | 1.458 | 0.100 |
| ddt | organochlorine | 1.00e-08 | 1.54e-04 | 0.002 | 0.001 | 0.003 |
| deltamethrin | type II pyrethroid | 1.00e-08 | 0.189 | 3.059 | 1.300 | 4.405 |
| diazepam | benzodiazepine | 1.00e-08 | 0.000 | 0.000 | 0.000 | 0.000 |
| dieldrin | cyclodiene organochlorine | 1.00e-08 | -0.091 | 1.048 | 0.620 | 0.071 |
| fipronil | phenylpyrazole | 1.00e-08 | -0.211 | 2.433 | 1.458 | 0.168 |
| gaba | endogenous transmitter | 1.00e-08 | 2.83e-06 | -1.30e-05 | -1.91e-05 | -8.67e-07 |
| imidacloprid | neonicotinoid | 1.00e-08 | 2.128 | -0.562 | -6.169 | -0.194 |
| ivermectin | avermectin | 1.00e-08 | 0.343 | 0.113 | -2.206 | -0.016 |
| nicotine | alkaloid | 1.00e-08 | 0.036 | 0.266 | 0.243 | 0.093 |
| nitenpyram | neonicotinoid | 1.00e-08 | 0.127 | 0.543 | 0.867 | 0.191 |
| permethrin | type I pyrethroid | 1.00e-08 | 0.024 | 0.395 | 0.165 | 0.500 |
| picrotoxin | plant convulsant | 1.00e-08 | -0.010 | 0.114 | 0.067 | 0.008 |
| spinosad | spinosyn | 1.00e-08 | 0.004 | 0.017 | 0.014 | 0.006 |
| sulfoxaflor | sulfoximine | 1.00e-08 | 6.67e-05 | 5.80e-04 | 4.51e-04 | 2.02e-04 |
| thiamethoxam | neonicotinoid | 1.00e-08 | 1.00e-04 | 8.69e-04 | 6.77e-04 | 3.03e-04 |
| acetamiprid | neonicotinoid | 1.00e-07 | 0.697 | -10.8 | -4.549 | -3.378 |
| acetylcholine | endogenous transmitter | 1.00e-07 | 0.004 | 0.034 | 0.027 | 0.012 |
| caffeine | methylxanthine | 1.00e-07 | 0.000 | 0.000 | 0.000 | 0.000 |
| chlordimeform | formamidine | 1.00e-07 | 0.003 | 0.000 | 0.022 | 0.000 |
| chlorpyrifos_oxon | organophosphate (oxon metabolite) | 1.00e-07 | -0.909 | -10.8 | -5.908 | -4.453 |
| clothianidin | neonicotinoid | 1.00e-07 | 0.809 | -10.8 | -5.269 | -4.774 |
| ddt | organochlorine | 1.00e-07 | 0.002 | 0.025 | 0.010 | 0.031 |
| deltamethrin | type II pyrethroid | 1.00e-07 | 0.699 | 11.3 | 5.173 | 24.0 |
| diazepam | benzodiazepine | 1.00e-07 | 0.000 | 0.000 | 0.000 | 0.000 |
| dieldrin | cyclodiene organochlorine | 1.00e-07 | -0.500 | 5.764 | 3.582 | 0.424 |
| fipronil | phenylpyrazole | 1.00e-07 | -0.809 | 9.328 | 6.108 | 0.738 |
| gaba | endogenous transmitter | 1.00e-07 | 8.94e-05 | -4.12e-04 | -6.05e-04 | -2.74e-05 |
| imidacloprid | neonicotinoid | 1.00e-07 | 2.853 | -10.8 | -6.169 | -5.680 |
| ivermectin | avermectin | 1.00e-07 | 0.924 | -0.021 | -5.423 | -0.059 |
| nicotine | alkaloid | 1.00e-07 | 0.318 | -0.754 | -2.107 | -0.260 |
| nitenpyram | neonicotinoid | 1.00e-07 | 0.697 | -10.8 | -4.549 | -3.378 |
| permethrin | type I pyrethroid | 1.00e-07 | 0.200 | 3.243 | 1.379 | 4.710 |
| picrotoxin | plant convulsant | 1.00e-07 | -0.091 | 1.048 | 0.620 | 0.071 |
| spinosad | spinosyn | 1.00e-07 | 0.039 | 0.157 | 0.133 | 0.055 |
| sulfoxaflor | sulfoximine | 1.00e-07 | 6.66e-04 | 0.006 | 0.005 | 0.002 |
| thiamethoxam | neonicotinoid | 1.00e-07 | 9.99e-04 | 0.009 | 0.007 | 0.003 |
| acetamiprid | neonicotinoid | 1.00e-06 | 0.973 | -10.8 | -6.169 | -6.173 |
| acetylcholine | endogenous transmitter | 1.00e-06 | 0.059 | 0.394 | 0.404 | 0.138 |
| caffeine | methylxanthine | 1.00e-06 | 0.000 | 0.000 | 0.000 | 0.000 |
| chlordimeform | formamidine | 1.00e-06 | 0.032 | 0.000 | 0.219 | 0.000 |
| chlorpyrifos_oxon | organophosphate (oxon metabolite) | 1.00e-06 | -0.990 | -10.8 | -6.169 | -5.555 |
| clothianidin | neonicotinoid | 1.00e-06 | 0.985 | -10.8 | -6.169 | -6.173 |
| ddt | organochlorine | 1.00e-06 | 0.015 | 0.246 | 0.103 | 0.308 |
| deltamethrin | type II pyrethroid | 1.00e-06 | 0.959 | 15.5 | 7.461 | 37.2 |
| diazepam | benzodiazepine | 1.00e-06 | 0.000 | 0.000 | 0.000 | 0.000 |
| dieldrin | cyclodiene organochlorine | 1.00e-06 | -0.909 | 10.5 | 7.000 | 0.851 |
| fipronil | phenylpyrazole | 1.00e-06 | -0.985 | 11.0 | 7.712 | 0.900 |
| gaba | endogenous transmitter | 1.00e-06 | 0.003 | -0.013 | -0.019 | -8.64e-04 |
| imidacloprid | neonicotinoid | 1.00e-06 | 2.989 | -10.8 | -6.169 | -6.173 |
| ivermectin | avermectin | 1.00e-06 | 1.480 | -1.837 | -6.169 | -0.184 |
| nicotine | alkaloid | 1.00e-06 | 0.855 | -10.8 | -5.559 | -5.404 |
| nitenpyram | neonicotinoid | 1.00e-06 | 0.973 | -10.8 | -6.169 | -6.173 |
| permethrin | type I pyrethroid | 1.00e-06 | 0.714 | 11.6 | 5.298 | 24.8 |
| picrotoxin | plant convulsant | 1.00e-06 | -0.500 | 5.764 | 3.582 | 0.424 |
| spinosad | spinosyn | 1.00e-06 | 0.333 | 0.483 | 1.146 | 0.170 |
| sulfoxaflor | sulfoximine | 1.00e-06 | 0.007 | 0.056 | 0.045 | 0.020 |
| thiamethoxam | neonicotinoid | 1.00e-06 | 0.010 | 0.083 | 0.067 | 0.029 |
| acetamiprid | neonicotinoid | 1.00e-05 | 0.998 | -10.8 | -6.169 | -6.173 |
| acetylcholine | endogenous transmitter | 1.00e-05 | 0.500 | -4.348 | -3.285 | -1.438 |
| caffeine | methylxanthine | 1.00e-05 | 0.000 | 0.000 | 0.000 | 0.000 |
| chlordimeform | formamidine | 1.00e-05 | 0.250 | 0.000 | 1.734 | 0.000 |
| chlorpyrifos_oxon | organophosphate (oxon metabolite) | 1.00e-05 | -0.999 | -10.8 | -6.169 | -5.550 |
| clothianidin | neonicotinoid | 1.00e-05 | 0.999 | -10.8 | -6.169 | -6.173 |
| ddt | organochlorine | 1.00e-05 | 0.133 | 2.162 | 0.913 | 2.982 |
| deltamethrin | type II pyrethroid | 1.00e-05 | 0.996 | 16.1 | 7.811 | 38.2 |
| diazepam | benzodiazepine | 1.00e-05 | 0.000 | 0.000 | 0.000 | 0.000 |
| dieldrin | cyclodiene organochlorine | 1.00e-05 | -0.990 | 11.0 | 7.758 | 0.900 |
| fipronil | phenylpyrazole | 1.00e-05 | -0.999 | 11.0 | 7.843 | 0.900 |
| gaba | endogenous transmitter | 1.00e-05 | 0.082 | -0.379 | -0.552 | -0.025 |
| imidacloprid | neonicotinoid | 1.00e-05 | 2.999 | -10.8 | -6.169 | -6.173 |
| ivermectin | avermectin | 1.00e-05 | 1.907 | -3.715 | -6.169 | -0.297 |
| nicotine | alkaloid | 1.00e-05 | 0.987 | -10.8 | -6.169 | -6.173 |
| nitenpyram | neonicotinoid | 1.00e-05 | 0.998 | -10.8 | -6.169 | -6.173 |
| permethrin | type I pyrethroid | 1.00e-05 | 0.962 | 15.6 | 7.487 | 37.2 |
| picrotoxin | plant convulsant | 1.00e-05 | -0.909 | 10.5 | 7.000 | 0.851 |
| spinosad | spinosyn | 1.00e-05 | 1.333 | -9.663 | -4.356 | -3.041 |
| sulfoxaflor | sulfoximine | 1.00e-05 | 0.063 | 0.408 | 0.425 | 0.143 |
| thiamethoxam | neonicotinoid | 1.00e-05 | 0.091 | 0.503 | 0.620 | 0.177 |
