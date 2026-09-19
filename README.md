# FlyLab

Virtual pharmacology bench: a public *Drosophila* connectome **plus** a ligand-receptor layer **plus** a vertebrate comparison panel.

https://github.com/pinkysworld/FlyLab

## Run the bench (Windows / macOS / Linux)

```bash
python -m pip install -e ".[dev]"
pytest
flylab serve
```

Open http://127.0.0.1:8765

Without the UI:

```bash
flylab occupancy imidacloprid --conc 1e-6
flylab assay --compound imidacloprid --conc 1e-6
```

## What v0.2 actually is

| Component | Status |
|---|---|
| Hill occupancy, cited teaching EC50s | working |
| Insect vs vertebrate scorecard | working |
| Reduced taste to MN9 circuit (Shiu-style control) | working |
| Local two-pane bench + notebook JSON export | working |
| Full FlyWire / MaleCNS weight matrix | **not loaded yet** |
| Live fly colony pairing | protocol draft only |

The taste circuit reproduces the published *direction*: sugar drives MN9; bitter vetoes it. Weights are **not** synapse counts from Janelia. Every notebook file says that.

## Paper draft

See `papers/IJRC_FlyLab_draft.md`.
