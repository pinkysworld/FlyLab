# Derived MaleCNS products (small, committed)

All files here are cut from the public **MaleCNS v1.0** release
(HHMI Janelia FlyEM, CC-BY 4.0). Raw feathers are never committed.

| File | What | Built by |
|---|---|---|
| `malecns_named_neighborhood.json` | 1-hop neighborhood of `MN9` + `DNp01`, edges ≥ 5 synapses (1126 nodes / 1360 edges) | `extract-malecns-subgraph` Action |
| `malecns_taste_motor_neighborhood.json` | 1-hop neighborhood of `MN9`, `DNp01` and labellar GRN types `LB1a-d`, `LB3b/c`, plus induced partner–partner edges ≥ 10 synapses (1841 nodes / 19066 edges) | same Action, `closure_min_weight=10` |
| `malecns_census_v1.json` | Traced-neuron census: transmitter, superclass and class counts, named cell IDs | `flylab.maps.malecns.load_census` on the atlas |
| `malecns_gustatory_seeds.json` | Labellar GRN body IDs by MaleCNS `type`; sweet/bitter working hypothesis | atlas `type` column |

Do not edit the neighborhood JSONs by hand. Re-run the Action (or
`flylab extract-subgraph`) if the map or thresholds change. Keep each file
≤ ~1 MB; raise `min_weight` / `closure_min_weight` rather than adding hops.
