# Architecture

The connectome is a **netlist**. Pharmacology is a **patch** on that netlist. They live in separate files so a reviewer can always see which number came from Janelia and which came from a paper.

Two rules shape everything below. **The evidence is typed**: what a source measured (`Kd`/`Ki` vs `EC50`/`IC50`/`Kb` vs nothing) and how far it sits from this compound/receptor/species jointly decide whether FlyLab may compute a binding occupancy, a functional *engagement*, or nothing at all — a row with no sourced value returns `None`, never a small number, and asking anyway raises `EvidenceTypeError`. **Every prediction carries its dependence**: the analysis layer can say, for any prediction, the weakest degraded graph model that still reproduces it.

## Dataflow

```
                       compound + concentration
                    (CLI · FastAPI · browser bridge)
                                  |
              +-------------------+-------------------+
              |     typed engagement engine           |
              |  theta = C^n / (value^n + C^n)        |
              |  pharm/evidence.py decides WHICH:     |
              |    Kd/Ki  -> binding_occupancy        |
              |    EC50.. -> functional_engagement    |
              |    else   -> not_modelled  (N/A)      |
              |  library.yaml v3: param_type, value,  |
              |  n, direction, relation, species,     |
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
   |  dependence.py  permutation p vs 4 degraded graphs ->    |
   |                 class + necessary information level      |
   |  baselines.py   ablation ladder A/B/C/D                  |
   |  robustness.py  25 gain specifications x 7 conclusions,  |
   |                 plus the amplify/buffer threshold grid   |
   |  uncertainty_global.py  Sobol' budget (Ishigami-checked) |
   |  voi.py         VOI_j = S_j Var(Y) -> ranked experiments |
   |  selectivity.py receptor SI vs circuit SI landscape      |
   |  fit/impact/layout  ·  mixtures  ·  prospective H1-H7    |
   |  validation/rank.py  concordance + source-overlap flag   |
   +----------------------------+-----------------------------+
                                v
   +----------------------------------------------------------+
   |  notebook JSON (schema 0.3)                              |
   |  engagement (+ model and reason) · gains · readouts ·    |
   |  uncertainty · warnings                                  |
   |  live_lab: null   (only a human may fill it)             |
   |  provenance: flylab version · git sha · library sha256   |
   |              map id + citation · rng seed · platform     |
   +----------------------------------------------------------+
```

## Layers and who owns what

| Layer | Package | Rule |
|---|---|---|
| Pharmacology | `flylab/pharm/` | no circuit knowledge; every constant sourced, tiered **and typed**; unsupported rows return N/A |
| Netlist | `data/derived/`, `flylab/maps/` | CI products, never hand-edited; feathers never committed |
| Runtime | `flylab/circuit/` | numerics only; gains arrive as a dict, never computed here |
| Assays | `flylab/assays/` | compose pharm + netlist + runtime into a notebook with warnings |
| Analysis | `flylab/analysis/` | **tests the instrument, and is allowed to fail it**: dependence, ablation, robustness, global uncertainty, VOI, nulls, selectivity, fits, impact, layout, predictions |
| Validation | `flylab/validation/` | compares the model with published orderings; records discrepancies; flags which comparisons share a source with the library and therefore are not out-of-sample |
| Interfaces | `flylab/cli.py`, `server.py`, `browser/`, `static/` | three transports, one implementation |
| Reproduction | `scripts/reproduce_paper.py` | the only writer of `papers/figures`, `papers/tables`, `papers/results.json` |

Every assay obtains its gains from `flylab.pharm.mechanisms.gains_from_occupancy`. There is no second copy of a patch rule anywhere in the tree — that is what makes the mechanism table (T2) a truthful description of the software.

## Three transports, one program

`flylab/cli.py` (typer), `flylab/server.py` (FastAPI) and `flylab/browser/bridge.py` (Pyodide, the same routes) all call the same functions. `tests/test_browser_bridge.py` asserts the HTTP and in-browser transports return equal JSON. See `docs/PAGES.md`.

## What is deliberately absent

- No vertebrate circuit. The vertebrate side is a receptor panel only, by design: no vertebrate connectome of comparable completeness exists.
- No fitted parameter. Nothing is fitted to animal data; the only fits are of the model's own curves, labelled `model_derived`. When one parameter is eventually fitted, the uncertainty budget says which: the engagement→gain transformation.
- No write path to `live_lab` from code.
- No second implementation of the pharmacology in JavaScript.
- No number produced from a receptor row that has no sourced value. This is enforced by the type system in `pharm/evidence.py`, not by reviewer vigilance.
