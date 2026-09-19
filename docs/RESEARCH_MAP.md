# FlyLab research map

This is the research plan, not a results paper. Software that exists today is marked **done**. Live-animal work is marked **blocked**.

## Claim the project is allowed to make

Same compound, same concentration, two scorecards:

1. Insect: occupancy on insect nAChR / RDL → gain patch on a *named* MaleCNS object.
2. Vertebrate: occupancy on a4b2 / GABA-A (teaching library).

One notebook JSON. Map version recorded. No invented EC50. No fake live n.

## Map objects (MaleCNS v1.0)

| Object | What it is | Status |
|---|---|---|
| Atlas census | ~165k traced bodies, consensus transmitters | **done** (`assay-cns`) |
| Named motor seeds | type `MN9`, `DNp01` | **done** (real body IDs) |
| 1-hop neighborhood | 1126 nodes / 1360 edges, weight ≥ 5 | **done** (Actions + `data/derived/`) |
| Labellar GRN types | Cell 2026: LB1a–d bitter (`Gr33a`-like), LB3b/c sweet (`Gr64f`-like) | **named in this map; not yet extracted** |
| Full 25M-edge LIF | every traced cell | **out of scope for v0** |

Public files (do not commit):

- `body-annotations-male-cns-v1.0-minconf-0.5.feather`
- `body-neurotransmitters-male-cns-v1.0.feather`
- `connectome-weights-male-cns-v1.0-minconf-0.5.feather` (~1.1 GB)

Cut on GitHub Actions. Keep only JSON ≤ ~200 KB in git.

## Pharmacological hypotheses (testable, not yet tested in vivo)

H1. At 1 µM, imidacloprid occupancy on insect nAChR ≥ 0.9 and on vertebrate a4b2 ≤ 0.2 in the teaching library.  
**Software check: done.** Live check: blocked.

H2. The same dose lowers model drive through ACh edges in the MN9/DNp01 neighborhood (`g_ACh` falls).  
**Software check: done.** Fit of `g_ACh` to PER: blocked.

H3. Fipronil at 1 µM lowers `g_GABA` (RDL antagonist) and barely moves vertebrate GABA-A.  
**Software check: done.** Live climbing/seizure pairing: blocked.

H4. Diazepam occupies vertebrate GABA-A and does not occupy insect RDL in this library.  
**Software check: done.**

H5. Bitter input vetoes sugar-driven MN9 (Shiu 2024 direction).  
**Software check: done on reduced_taste_v0.** Map-extracted LB1→MN9 path: next extract, not a new thesis.

## Assays

| Assay | Command / API | Reads | Writes |
|---|---|---|---|
| Occupancy | `flylab occupancy` | library.yaml | receptor table |
| Taste (reduced) | `flylab assay` | 5 rate units | MN9 Hz, g_ACh |
| Neighborhood | `flylab assay-subgraph` | committed JSON | named-cell Hz, g_ACh, g_GABA |
| Census | `flylab assay-cns` | atlas feathers | NT counts, excitation index |
| Live PER | `protocols/per.md` | flies | `live_lab` field only |
| Live climbing | `protocols/climbing.md` | flies | `live_lab` field only |

## Parameter policy

Allowed to fit later, at most one coefficient per paper: `g_ACh` or one EC50.  
Not allowed: refitting the MaleCNS weights, inventing transmitters, hiding the reduced-circuit warning.

## Paper path (IJRC)

`papers/IJRC_FlyLab_draft.md` is a methods / software note.

Submit when either:

- the journal accepts a software article with no live table, or
- one PER or climbing table is in the notebook.

Do not submit H2–H5 as empirical biology without that table.

## Next extract (Actions, not a rewrite)

```text
flylab extract-subgraph --types MN9,DNp01,LB1a,LB1b,LB1c,LB1d,LB3b,LB3c
```

If MaleCNS `type` strings differ (`LB1a` vs `lbGRN_LB1a`), prefix match in `flylab.maps.extract.seed_ids` still keeps MN9/DNp01. Empty extra types are recorded as `[]`, not a crash.

## What this map is not

Not a new connectome. Not GLP toxicology. Not a 166k-cell product. Not a claim that the teaching EC50s are the fly's EC50s.
