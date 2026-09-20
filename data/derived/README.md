# Derived MaleCNS products (small, committed)

All files here are cut from the public **MaleCNS v1.0** release
(HHMI Janelia FlyEM, CC-BY 4.0). Raw feathers are never committed.

| File | What | Built by |
|---|---|---|
| `malecns_named_neighborhood.json` | 1-hop neighborhood of `MN9` + `DNp01`, edges ≥ 5 synapses (1126 nodes / 1360 edges) | `extract-malecns-subgraph` Action |
| `malecns_taste_motor_neighborhood.json` | 1-hop neighborhood of `MN9`, `DNp01` and labellar GRN types `LB1a-d`, `LB3b/c`, plus induced partner–partner edges ≥ 10 synapses (1841 nodes / 19066 edges) | same Action, `closure_min_weight=10` |
| `malecns_scale_1k.json` | Smallest rung of the **scale ladder**: induced cut on the first 1000 cells of a deterministic ranked-BFS growth from the taste-motor seed set, edges ≥ 5 synapses (1000 nodes / 22857 edges, mean degree 22.9). Carries its own `recipe`, `census` and `composition` blocks | `extract-malecns-subgraph` Action, ladder step |
| `malecns_census_v1.json` | Traced-neuron census: transmitter, superclass and class counts, named cell IDs | `flylab.maps.malecns.load_census` on the atlas |
| `malecns_gustatory_seeds.json` | Labellar GRN body IDs by MaleCNS `type`; sweet/bitter working hypothesis | atlas `type` column |

Do not edit the neighborhood JSONs by hand. Re-run the Action (or
`flylab extract-subgraph`) if the map or thresholds change. Keep each file
≤ ~1 MB; raise `min_weight` / `closure_min_weight` rather than adding hops.

## The scale ladder

`flylab.maps.extract.extract_ladder` builds nested cuts at 1k / 5k / 10k / 25k
/ 50k cells by one fixed rule (`LADDER_RECIPE`), so that a connectome-
dependence verdict can be compared across sizes without the construction
changing underneath it. **Only the 1k rung is small enough for git** (1.1 MB);
5k is 13 MB and 50k is 105 MB, so the rest are uploaded by the Action as the
`malecns-scale-ladder` artifact and are not committed. Drop them anywhere on
`flylab.circuit.rate.DATA_DIRS` (e.g. `~/.flylab/`) and
`flylab.analysis.scale.available_cuts()` will find them.
