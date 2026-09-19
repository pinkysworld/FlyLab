# FlyLab

Virtual pharmacology bench on **MaleCNS v1.0** (adult male *Drosophila* brain + ventral nerve cord).

https://github.com/pinkysworld/FlyLab

One compound, two scorecards (insect circuit / census + vertebrate receptors), one notebook JSON.

## Install

```bash
git clone https://github.com/pinkysworld/FlyLab
cd FlyLab
python -m pip install -e ".[dev]"
pytest
flylab serve    # http://127.0.0.1:8765
```

## Commands that work today

```bash
flylab occupancy imidacloprid --conc 1e-6
flylab assay --compound imidacloprid --conc 1e-6          # reduced taste→MN9
flylab assay-cns --compound imidacloprid --conc 1e-6      # MaleCNS census (~165k traced)
flylab assay-subgraph --compound imidacloprid --conc 1e-6 # real 1-hop MN9+DNp01 graph
flylab list-drugs
```

Atlas (optional, for census refresh):

```bash
flylab download-malecns           # ~55 MB annotations + transmitters
flylab download-malecns --full    # + 1.1 GB weights
flylab extract-subgraph           # rebuild data/derived/malecns_named_neighborhood.json
```

Or let GitHub Actions do the 1.1 GB cut: **Actions → extract-malecns-subgraph**.

## What is real vs modelled

| Object | Status |
|---|---|
| MaleCNS named cells MN9, DNp01 | real body IDs |
| 1-hop neighborhood (1126 cells, 1360 edges) | real synapse counts from the public matrix |
| Traced-neuron NT census | real map counts |
| Hill occupancy / vertebrate panel | literature-order teaching EC50s |
| Rate dynamics on the subgraph | model, predicted-transmitter signs |
| Full 166k LIF | not shipped |
| Live fly pairing | protocol drafts only |

## Paper

`papers/IJRC_FlyLab_draft.md` — architecture draft. Do not submit as empirical biology until a live table exists.
