# FlyLab handoff (2026-09-19)

Read this first if you are Claude, GPT, or another agent continuing the repo.
Owner: https://github.com/pinkysworld/FlyLab  
HEAD at write time: `7799126` on `main`.  
Package version: `0.4.0` in `pyproject.toml` (server banner says 0.4.1).

## What the project is

A **virtual pharmacology bench** on the public **MaleCNS v1.0** adult male *Drosophila* CNS (brain + VNC). One compound + concentration → insect occupancy + vertebrate occupancy → notebook JSON.

It is **not** a new connectome, not GLP tox, not a 166k-cell LIF product, not live-animal results.

Canonical plan: `docs/RESEARCH_MAP.md`  
Paper draft: `papers/IJRC_FlyLab_draft.md` (methods / software note only)  
Target venue intent: ijrcom.org / IJRC.

## Do this first

```bash
git pull
python -m pip install -e ".[dev]"
pytest -q
flylab occupancy imidacloprid --conc 1e-6
flylab assay-subgraph --compound imidacloprid --conc 1e-6
flylab serve   # http://127.0.0.1:8765
```

CI tests: `.github/workflows/tests.yml` (was green on earlier pushes).

## What already works

| Piece | Where |
|---|---|
| Hill occupancy + YAML library | `flylab/pharm/` — imidacloprid, nitenpyram, nicotine, diazepam, fipronil |
| Reduced taste→MN9 (Shiu-style bitter veto) | `flylab/circuit/reduced_taste.py`, `flylab assays` |
| MaleCNS census (~165k traced, NT counts, named MN9/DNp01 IDs) | `flylab/assays/wholens.py` needs local atlas **or** falls back depending on code path — check `test_wholens.py` |
| Real 1-hop neighborhood | `data/derived/malecns_named_neighborhood.json` — **1126 nodes, 1360 edges**, seeds MN9 `[10331, 16949]`, DNp01 `[10001, 10010]` |
| Rate model + `g_ACh` / `g_GABA` patches | `flylab/assays/subgraph.py` |
| Bench UI (taste / whole CNS / MN9-DNp01) | `flylab/static/index.html` + `flylab/server.py` |
| Extractor (Arrow batch, no full pandas load of 1.1 GB) | `flylab/maps/extract.py` |
| Actions cut of the 1.1 GB file | `.github/workflows/malecns-subgraph.yml` |
| Protocols (no data) | `protocols/per.md`, `protocols/climbing.md` |

Teaching EC50s live in `flylab/pharm/library.yaml`. They are literature-order, not fitted.

### Occupancy patch rules (do not silently change)

- insect nAChR agonist: `g_ACh = max(0.05, 1 + 0.4θ − 1.6θ2)`
- insect nAChR antagonist: `g_ACh = max(0.05, 1 − θ)`
- insect RDL antagonist: `g_GABA = max(0.05, 1 − θ)`

Imidacloprid at 1e-6 M should show high insect nAChR occupancy and low vertebrate a4b2. Fipronil should drop `g_GABA`, not `g_ACh`.

## Big files (never commit)

From `https://storage.googleapis.com/flyem-male-cns/v1.0/connectome-data/flat-connectome`:

- `body-annotations-male-cns-v1.0-minconf-0.5.feather`
- `body-neurotransmitters-male-cns-v1.0.feather`
- `connectome-weights-male-cns-v1.0-minconf-0.5.feather` (~1.1 GB)

Local dir: `$FLYLAB_MALECNS` or `~/.flylab/malecns_v1`.  
`.gitignore` already drops `*.feather`.  
GitHub Action downloads weights on the runner, writes JSON, deletes the feather.

Weight columns used by extractor: `body_pre`, `body_post`, `weight`. Stream with `pyarrow.ipc.open_file` batches — a full pandas read OOMs.

## Issues

- #1 Gate 0 occupancy CLI — **closed**
- #2 Gate 1 map-extracted taste GRNs — **OPEN**. Reduced circuit veto works. LB* seeds are named but not proven in MaleCNS `type` strings.
- #3 Gate 2 five-drug panel — **closed** (teaching EC50s)

## Next work (in order)

1. **Prove LB type names on MaleCNS**  
   After `flylab download-malecns` (atlas is enough):
   ```python
   import pandas as pd
   ann = pd.read_feather("~/.flylab/malecns_v1/body-annotations-male-cns-v1.0-minconf-0.5.feather")
   print(ann[ann.type.fillna("").str.contains("LB1|LB3|GRN|MN9", case=False)]["type"].value_counts())
   ```
   Cell 2026 gustatory paper names labellar types **LB1a–d** (bitter / Gr33a-like) and **LB3b/c** (sweet / Gr64f-like). MaleCNS `type` strings may differ (`LB1a` vs `lbGRN_LB1a`). `seed_ids` already does exact then prefix match.

2. **Re-run extract with those types** if counts are non-zero:
   ```bash
   flylab extract-subgraph --types MN9,DNp01,LB1a,LB1b,LB1c,LB1d,LB3b,LB3c
   ```
   Or Actions → `extract-malecns-subgraph`. Do not commit the feather. If JSON grows past ~1 MB, raise `min_weight` or drop hops; do not git-lfs the matrix.

3. **Wire sugar/bitter drive onto extracted LB cells** in `subgraph.py` instead of driving only MN9/DNp01. Keep `reduced_taste_v0` as the directional control until the map path also vetoes MN9.

4. **Do not** start a 166k LIF. Out of scope for v0 (`docs/RESEARCH_MAP.md`).

5. **Do not** invent live PER/climbing numbers. Write only into notebook field `live_lab` when a real table exists. Protocols are drafts.

6. Paper: expand `papers/IJRC_FlyLab_draft.md` only as a software/methods note. Do not claim H2–H5 as biology without live n.

## Hard constraints for the next model

- Organism is **fly** (*Drosophila*), not flea. Map is **MaleCNS v1.0**, not FlyWire-only (MN9/DNp01 need the cord).
- Novelty vs RatCVS / PharmVR: those have no named CNS cells. Fly sims have no occupancy. Dual scorecard is the point.
- At most **one** fitted parameter later (`g_ACh` or one EC50). Do not refit MaleCNS weights.
- Notebooks must keep warnings when the circuit is reduced or hops-limited.
- Windows + macOS matter; current UI is FastAPI + static HTML, not Electron.

## Useful commands

```bash
flylab list-drugs
flylab occupancy fipronil --conc 1e-6
flylab assay --compound imidacloprid --conc 1e-6
flylab assay-cns --compound imidacloprid --conc 1e-6
flylab assay-subgraph --compound fipronil --conc 1e-6
flylab download-malecns          # ~55 MB atlas
flylab download-malecns --full   # + 1.1 GB weights
flylab extract-subgraph --types MN9,DNp01 --hops 1 --min-weight 5
```

APIs: `POST /api/assay/taste`, `/api/assay/cns`, `/api/assay/subgraph`; `GET /api/assay/subgraph/dose-response`.

## Files not to rewrite from scratch

`flylab/maps/extract.py` (streaming), `data/derived/malecns_named_neighborhood.json` (CI product), `flylab/pharm/library.yaml` (cite sources if you change EC50s).
