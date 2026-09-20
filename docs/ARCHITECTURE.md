# Architecture

The connectome is a **netlist**. Pharmacology is a **patch** on that netlist. They live in separate files so a reviewer can always see which number came from Janelia and which came from a paper's EC50.

## Dataflow

```
                       compound + concentration
                    (CLI · FastAPI · browser bridge)
                                  |
              +-------------------+-------------------+
              |            occupancy engine           |
              |  Hill θ = C^n / (EC50^n + C^n)        |
              |  library.yaml: EC50, n, direction,    |
              |  efficacy, source, evidence_tier      |
              |  modifiers: genotype · mixture ·      |
              |             exposure C(t) · MC jitter |
              +---------+-------------------+---------+
                        |                   |
            insect panel|                   |vertebrate panel
   nAChR RDL GluCl AChE |                   | a4b2 a7 GABA-A GlyR
        Nav OctR        |                   | AChE Nav1.x
                        v                   v
        +-------------------------+   +--------------------------+
        | mechanism table         |   | scored at the same dose. |
        | (pharm/mechanisms.py)   |   | NEVER patches a circuit. |
        | receptor+direction ->   |   | no vertebrate connectome |
        | g_ach g_gaba g_glu      |   | exists at this quality.  |
        | g_oct g_nav ach_tone    |   +--------------------------+
        +------------+------------+                |
                     |                             |
     optional: expression weighting                |
     (per-node gains, sensitivity only)            |
                     v                             |
   +-------------------------------------+         |
   |  MaleCNS v1.0 netlist (CC-BY 4.0)   |         |
   |  named      1126 nodes / 1360 edges |         |
   |  taste_motor 1841 / 19066 (+closure)|         |
   |  census     165122 traced cells     |         |
   +------------------+------------------+         |
                      v                            |
   +-------------------------------------+         |
   |  circuit runtime                    |         |
   |  rate: clipped leaky iteration      |         |
   |        (deterministic, ~50 ms)      |         |
   |  LIF : Shiu-style, dt 0.1 ms,       |         |
   |        + background Poisson drive   |         |
   +------------------+------------------+         |
                      v                            |
   +-------------------------------------+         |
   |  readouts: MN9 Hz, DNp01 Hz,        |         |
   |  mean Hz, bitter veto ratio,        |         |
   |  per-node / per-edge impact         |         |
   +------------------+------------------+         |
                      |                            |
                      v                            v
   +----------------------------------------------------------+
   |  analysis + validation                                   |
   |  null models (4 degradations, z and permutation p)       |
   |  selectivity landscape (receptor SI vs circuit SI)       |
   |  Hill fit + bootstrap CI  ·  sensitivity tornado         |
   |  rank validation vs published orders (Spearman, Kendall) |
   |  mixtures (Bliss / Loewe) vs published verdicts          |
   |  pre-registered predictions H1-H7 + power                |
   +----------------------------+-----------------------------+
                                v
   +----------------------------------------------------------+
   |  notebook JSON (schema 0.3)                              |
   |  occupancy · gains · readouts · uncertainty · warnings   |
   |  live_lab: null   (only a human may fill it)             |
   |  provenance: flylab version · git sha · library sha256   |
   |              map id + citation · rng seed · platform     |
   +----------------------------------------------------------+
```

## Layers and who owns what

| Layer | Package | Rule |
|---|---|---|
| Pharmacology | `flylab/pharm/` | no circuit knowledge; every constant sourced and tiered |
| Netlist | `data/derived/`, `flylab/maps/` | CI products, never hand-edited; feathers never committed |
| Runtime | `flylab/circuit/` | numerics only; gains arrive as a dict, never computed here |
| Assays | `flylab/assays/` | compose pharm + netlist + runtime into a notebook with warnings |
| Analysis | `flylab/analysis/` | tests the instrument: nulls, selectivity, fits, impact, layout, predictions |
| Validation | `flylab/validation/` | compares the model with published orderings; records discrepancies |
| Interfaces | `flylab/cli.py`, `server.py`, `browser/`, `static/` | three transports, one implementation |
| Reproduction | `scripts/reproduce_paper.py` | the only writer of `papers/figures`, `papers/tables`, `papers/results.json` |

Every assay obtains its gains from `flylab.pharm.mechanisms.gains_from_occupancy`. There is no second copy of a patch rule anywhere in the tree — that is what makes the mechanism table (T2) a truthful description of the software.

## Three transports, one program

`flylab/cli.py` (typer), `flylab/server.py` (FastAPI, 34 routes) and `flylab/browser/bridge.py` (Pyodide, same 34 routes) all call the same functions. `tests/test_browser_bridge.py` asserts the HTTP and in-browser transports return equal JSON. See `docs/PAGES.md`.

## What is deliberately absent

- No vertebrate circuit. The vertebrate side is occupancy only, by design.
- No fitted parameter. v0.5 fits nothing to animal data; the only fits are of the model's own curves, labelled `model_derived`.
- No write path to `live_lab` from code.
- No second implementation of the pharmacology in JavaScript.
