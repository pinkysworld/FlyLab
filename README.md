# FlyLab

Virtual pharmacology bench on the **MaleCNS v1.0** adult male *Drosophila* nervous system (brain + ventral nerve cord).

https://github.com/pinkysworld/FlyLab

## Map choice

MaleCNS, not FlyWire-only: motor neurons and descending neurons (MN9, DNp01) live in the cord as well as the brain.

## Install and run

```bash
python -m pip install -e ".[dev]"
flylab download-malecns          # ~55 MB atlas (annotations + transmitters)
# flylab download-malecns --full # also the 1.1 GB synapse-weight matrix
pytest
flylab serve                     # http://127.0.0.1:8765
flylab assay-cns --compound imidacloprid --conc 1e-6
```

## What the numbers are

- **Taste assay:** reduced MN9 circuit (directional Shiu control).
- **Whole CNS:** real traced-neuron census from MaleCNS (~165k traced cells; ACh / GABA / Glu counts) + occupancy-patched excitation index + vertebrate panel.
- Notebook JSON export for both.

The whole-CNS layer is **not** yet a 166k-cell LIF on the weight matrix. The notebook says so. Named cells MN9 and DNp01 are looked up from the real map.

## Paper draft

`papers/IJRC_FlyLab_draft.md`
