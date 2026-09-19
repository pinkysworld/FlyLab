# FlyLab

Virtual pharmacology bench on **MaleCNS v1.0** (adult male *Drosophila* brain + ventral nerve cord).

https://github.com/pinkysworld/FlyLab

**Agents / next model:** read [`HANDOFF.md`](HANDOFF.md) then [`docs/RESEARCH_MAP.md`](docs/RESEARCH_MAP.md).

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
flylab assay --compound imidacloprid --conc 1e-6
flylab assay-cns --compound imidacloprid --conc 1e-6
flylab assay-subgraph --compound imidacloprid --conc 1e-6
flylab list-drugs
```

Atlas / weights (optional):

```bash
flylab download-malecns
flylab download-malecns --full
flylab extract-subgraph --types MN9,DNp01
```

Or **Actions → extract-malecns-subgraph** for the 1.1 GB cut. Only the neighborhood JSON is committed.

## Paper

`papers/IJRC_FlyLab_draft.md` — architecture draft. No live-animal dataset.
