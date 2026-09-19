# FlyLab

A virtual pharmacology laboratory on the public *Drosophila* connectome.

FlyLab is not another Doom-fly. It does not treat the MaleCNS / FlyWire wiring diagram as a finished animal. It adds a **ligand–receptor layer**, a **vertebrate comparison panel**, and a path to **live fly assays** so a compound can be scored the same way a teaching lab scores a dose–response — except the readout is a named circuit, not a pithed-rat blood-pressure trace.

Repository: https://github.com/pinkysworld/FlyLab

## What works today (v0)

After install you can compute receptor occupancy for a small built-in library and print an insect vs vertebrate comparison. That is the pharmacology object. The connectome runtime is Phase 1.

```bash
python -m pip install -e .
flylab occupancy imidacloprid --conc 1e-6
flylab compare imidacloprid nitenpyram nicotine diazepam --conc 1e-6
```

## What “practically usable” means

A lab is usable when a student or researcher can:

1. pick a named compound and a concentration,
2. see occupancy on insect vs vertebrate receptor classes,
3. (Phase 1+) inject that occupancy into a real connectome circuit,
4. read a named motor output (MN9 feeding, giant-fibre escape, climbing proxy),
5. export a notebook JSON they can put next to a vial assay.

If step 3 is missing, this is still a receptor calculator, not a fly lab. If step 5 is missing, it is a demo, not a lab.

Read **[ROADMAP.md](ROADMAP.md)** for gates. Do not skip Gate 0.

## What this is not

- Not a reconstructed flea brain.
- Not a full mouse synaptic connectome.
- Not a regulatory replacement for GLP toxicology.
- Not a claim that wiring + LIF = a living fly.
- Not RatCVS with insect clip-art.

## Maps we use (do not re-map)

| Resource | What it is | License |
|---|---|---|
| FlyWire v783 | Adult female brain ~139k neurons | CC-BY |
| MaleCNS v1.0 | Adult male brain + VNC ~166.7k neurons | CC-BY |
| Predicted transmitters | Per-neuron NT classes shipped with the maps | with the map |

Cite the map papers in every manuscript and in `CITATION.cff`.

## Repo layout

```
flylab/            package
  pharm/           occupancy, drugs, vertebrate panel
  circuit/         connectome loaders (Phase 1)
  assays/          named in-silico assays
  notebook/        export schema
data/              download scripts only — never commit the 1 GB graphs
protocols/         live *Drosophila* bench protocols
docs/              architecture and novelty
```

## Status

| Gate | Name | State |
|---|---|---|
| 0 | Honesty + repo contract | **open / this commit** |
| 1 | Taste circuit reproduces Shiu / flypoke control | not started |
| 2 | Four-drug pharmacology layer on that circuit | not started |
| 3 | Web bench UI | not started |
| 4 | Vertebrate panel beside the fly | not started |
| 5 | Live fly protocol pack + first paired dataset | not started |
| 6 | IJRC manuscript | not started |

## Licence

Code: MIT. Connectome files remain under their own CC-BY terms and must be downloaded by the user.
