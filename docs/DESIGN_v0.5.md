# FlyLab v0.5 design contract

This file is the interface contract for the v0.5 build. Agents working on
separate packages must implement exactly these names and JSON shapes so the
server, UI, CLI, tests and paper can integrate without renegotiation.

Non-negotiable scientific rules (from `HANDOFF.md` / `docs/RESEARCH_MAP.md`):

- Organism is *Drosophila*; map is MaleCNS v1.0. No new connectome, no 166k LIF.
- Existing patch rules are preserved exactly (see `flylab/pharm/mechanisms.py`).
- No invented EC50: every library value carries `source`, `evidence_tier`.
- Every notebook keeps `warnings` when circuit is reduced / hops-limited.
- No live-animal numbers are ever generated. `live_lab` stays `null` unless a
  user imports a real table.
- At most one fitted parameter later. v0.5 fits **nothing** to animal data.
  The only fits are Hill fits of the *model's own* dose-response curves
  (these are descriptive summaries, labelled `model_derived`).

## Package layout (owner in brackets)

```
flylab/
  pharm/
    occupancy.py        (existing; keep hill_occupancy, compare_compound)
    mechanisms.py       [pharm]  mechanism -> gain patch rules (single source of truth)
    binding.py          [pharm]  competitive binding, Schild shift, operational model
    exposure.py         [pharm]  1-compartment exposure -> C(t) -> occupancy(t)
    uncertainty.py      [pharm]  Monte-Carlo on log-EC50 -> occupancy CI
    library.yaml        [pharm]  expanded, sourced library (schema below)
  circuit/
    reduced_taste.py    (existing; untouched)
    rate.py             [circuit] rate network extracted from assays/subgraph.py
    lif.py              [circuit] LIF spiking network on the neighborhood graph
  assays/
    taste.py            (existing)
    wholens.py          (existing; add census fallback JSON)
    subgraph.py         (existing API kept; delegates to circuit/rate.py)
    spiking.py          [circuit] LIF assay on neighborhood graph
    ensemble.py         [circuit] replicates, CIs, sensitivity, circuit IC50
    experiment.py       [circuit] batch design: compounds x concs x replicates
  analysis/
    impact.py           [circuit] per-edge / per-node drug impact, path ranking
    fit.py              [circuit] Hill fit of model dose-response with bootstrap CI
  notebook/
    schema.py           (bump to 0.3: add provenance block)
  server.py             [ui]  FastAPI endpoints listed below
  static/               [ui]  bench UI (no build step, CDN libs only)
  cli.py                [ui]  new subcommands
data/derived/
  malecns_named_neighborhood.json   (CI product; never rewritten by hand)
  malecns_census_v1.json            [lead] small census fallback
```

## library.yaml schema (v2)

```yaml
schema_version: 2
receptors:                       # receptor classes, with organism and default transmitter
  insect_nAChR:        {organism: insect,      transmitter: acetylcholine, family: pLGIC-cation}
  insect_RDL:          {organism: insect,      transmitter: gaba,          family: pLGIC-anion}
  insect_GluCl:        {organism: insect,      transmitter: glutamate,     family: pLGIC-anion}
  insect_AChE:         {organism: insect,      transmitter: acetylcholine, family: enzyme}
  insect_Nav:          {organism: insect,      transmitter: null,          family: channel}
  insect_OctR:         {organism: insect,      transmitter: octopamine,    family: GPCR}
  vertebrate_nAChR_a4b2: {organism: vertebrate, transmitter: acetylcholine, family: pLGIC-cation}
  vertebrate_nAChR_a7:   {organism: vertebrate, transmitter: acetylcholine, family: pLGIC-cation}
  vertebrate_GABA_A:     {organism: vertebrate, transmitter: gaba,          family: pLGIC-anion}
  vertebrate_GlyR:       {organism: vertebrate, transmitter: glycine,       family: pLGIC-anion}
  vertebrate_AChE:       {organism: vertebrate, transmitter: acetylcholine, family: enzyme}
  vertebrate_Nav1_x:     {organism: vertebrate, transmitter: null,          family: channel}
compounds:
  imidacloprid:
    name: Imidacloprid
    class: neonicotinoid
    cas: 138261-41-3
    receptors:
      insect_nAChR:
        ec50_M: 2.0e-08
        n: 1.2
        direction: agonist          # agonist | partial_agonist | antagonist | positive_modulator | negative_modulator | inhibitor | none
        efficacy: 1.0               # operational-model tau-like relative efficacy (0..1), optional
        source: "Dederer et al. 2011, Insect Biochem Mol Biol 41:51 (hybrid Dalpha2/beta2 EC50 order)"
        evidence_tier: literature_order   # literature_order | class_placeholder | measured_fit
      ...
```

Rules for the pharm agent:
- Keep the five existing compounds and their existing insect_nAChR / insect_RDL /
  vertebrate_nAChR_a4b2 / vertebrate_GABA_A values (tests depend on them).
- Add compounds only with a named source string. Tier `class_placeholder`
  for receptors where no number exists (ec50 0.01, direction none), as today.
- `compare_compound` output stays backward compatible (same keys), plus
  `selectivity` block (insect_over_vertebrate ratio per receptor pair) and
  `evidence_tier` per row.

## Mechanism -> gain patch (`flylab/pharm/mechanisms.py`)

```python
def gains_from_occupancy(rows: list[dict]) -> dict[str, float]
    # returns {"g_ach": .., "g_gaba": .., "g_glu": .., "g_oct": .., "g_nav": .., "ach_tone": ..}
    # existing rules, verbatim:
    #   insect nAChR agonist:     g_ach  = max(0.05, 1 + 0.4*th - 1.6*th*th)
    #   insect nAChR antagonist:  g_ach  = max(0.05, 1 - th)
    #   insect RDL antagonist:    g_gaba = max(0.05, 1 - th)
    #   insect RDL agonist/PAM:   g_gaba = max(0.05, 1 + 0.4*th)   (kept from subgraph.py)
    # new rules (document each in docstring, keep monotone and bounded):
    #   insect GluCl agonist (avermectins): g_glu = max(0.05, 1 + 0.8*th)  (GluCl is inhibitory; more Cl- current)
    #   insect AChE inhibitor:  ach_tone = 1 + 2.0*th ; g_ach = ach_tone * agonist-desensitisation curve at th_eff
    #   insect Nav modulator (pyrethroid/DDT): g_nav = 1 + 1.5*th  (applied as global excitability multiplier)
    #   insect OctR agonist: g_oct = 1 + 0.5*th (applied to octopamine edges)
```
`assays/taste.py`, `assays/wholens.py`, `assays/subgraph.py` must all call this
function instead of carrying their own copy.

## Circuit runtime contracts

```python
# flylab/circuit/rate.py
class RateNetwork:
    def __init__(self, graph: dict, sign: dict[str, float] = SIGN): ...
    def run(self, drive: dict[int, float], gains: dict[str, float], steps=80, alpha=0.3) -> np.ndarray  # rates per node

# flylab/circuit/lif.py
class LIFNetwork:
    """Shiu-style LIF: tau_m=20ms, V_rest=-52, V_th=-45, V_reset=-52, refractory 2.2ms,
    exponential synapse tau_s=5ms, weight = w_scale * synapse_count * sign(nt) * gain(nt).
    Input: Poisson drive (Hz) per seed node. dt = 0.1 ms."""
    def __init__(self, graph, seed=0, w_scale=..., dt_ms=0.1): ...
    def run(self, drive_hz: dict[int, float], gains: dict[str, float], t_ms=500.0) -> SpikeResult
@dataclass
class SpikeResult:
    spikes: list[tuple[float, int]]      # (t_ms, node_index)
    rates_hz: np.ndarray                 # per node, over the last analysis window
    t_ms: float
```

Assay functions return notebooks (see schema 0.3). Function names:

```python
run_subgraph_assay(compound, conc_M, graph_path_arg=None, drive_hz=40., steps=80)  # existing, keep
run_spiking_assay(compound, conc_M, graph_path_arg=None, drive_hz=40., t_ms=500., seed=0)
run_ensemble(assay: str, compound, conc_M, n_rep=8, ec50_sd_log10=0.3, seed=0)   -> notebook with CI
sensitivity(assay, compound, conc_M, params=("g_ach_coef","ec50","drive_hz"))     -> tornado table
circuit_ic50(assay, compound, readout="mn9_hz", concs=..., n_boot=200)             -> fit + CI
run_experiment(design: dict) -> {"design":..., "rows":[...], "summary":..., "notebooks":[...]}
```

Experiment design JSON:

```json
{"assay": "subgraph", "compounds": ["imidacloprid","nicotine"],
 "concs_M": [1e-9,1e-8,1e-7,1e-6,1e-5], "replicates": 4, "seed": 0,
 "readouts": ["mn9_hz","dnp01_hz","g_ach","g_gaba"]}
```

## Notebook schema 0.3

```json
{
 "flylab_notebook_version": "0.3",
 "assay": "...",
 "created_utc": "...",
 "provenance": {"flylab_version": "...", "git_sha": "...|null", "library_sha256": "...",
                "map": {"name":..., "version":..., "citation":...}, "rng_seed": 0,
                "platform": "..."},
 "map": {...},                      # kept for backward compat
 "compound": ..., "concentration_M": ...,
 "occupancy": [...],                # rows from compare_compound
 "gains": {...},                    # from gains_from_occupancy
 "readouts": {...},
 "uncertainty": null | {"n_rep":..., "ci": {"readout": [lo, hi]}, "method": "..."},
 "live_lab": null,
 "warnings": [...]
}
```

## Server endpoints (v0.5)

```
GET  /api/meta                      version, map, library hash, compounds, receptors, mechanisms
GET  /api/drugs                     (existing)
GET  /api/drugs/{key}               full library entry
GET  /api/occupancy                 (existing) + selectivity
GET  /api/occupancy/curve?compound= per-receptor occupancy vs conc (1e-11..1e-3)
POST /api/assay/taste               (existing)
POST /api/assay/cns                 (existing; works with census fallback)
POST /api/assay/subgraph            (existing)
POST /api/assay/spiking             {compound, conc_M, drive_hz, t_ms, seed}
POST /api/assay/ensemble            {assay, compound, conc_M, n_rep, seed}
POST /api/analysis/sensitivity      {assay, compound, conc_M}
POST /api/analysis/ic50             {assay, compound, readout}
POST /api/experiment                design JSON -> table
GET  /api/graph                     neighborhood nodes/edges (+ layout positions) for the viewer
POST /api/graph/impact              {compound, conc_M} -> per-edge effective weight change, top nodes
POST /api/exposure                  {compound, dose, route, t_h} -> C(t), occupancy(t)
GET  /api/assay/{name}/dose-response (existing two kept)
```

## UI panels (`flylab/static/`)

Single-page bench, tabs: **Dose** (compound, conc, ladder), **Scorecard**
(insect vs vertebrate occupancy + selectivity), **Curves** (occupancy vs conc,
circuit readout vs conc with CI band, model IC50), **Circuit** (interactive
neighborhood graph, nodes coloured by transmitter, sized by rate, edges by
effective weight; before/after slider), **Spikes** (raster + PSTH for named
cells), **Exposure** (C(t) and occupancy(t)), **Experiment** (design table,
run, results grid, CSV export), **Notebook** (JSON, warnings, export, import
of a `live_lab` CSV). Plotly + cytoscape.js from CDN only. Works offline once
loaded except CDN scripts; no build step; Windows/macOS.
