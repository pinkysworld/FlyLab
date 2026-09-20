# FlyLab tutorial

Seven short lessons. Each one gives you a command you can paste, the output it
actually produced on a clean checkout, and a note on **what that output tells
you and what it does not**.

Every command and every block of output on this page was run. Nothing is
illustrative.

You can also run the lessons:

```bash
flylab tutorial --list        # the index below
flylab tutorial 3             # run lesson 3's real commands, with its caveats
```

---

> ## If you only read one thing
>
> FlyLab produces four kinds of number and they look identical on screen:
>
> 1. **a measurement of the system being simulated** — the MaleCNS synaptic
>    wiring, and nothing else;
> 2. **a literature parameter**, which may have been measured in an aphid, on a
>    chimeric receptor, or as a statement about a chemical class;
> 3. **a transformation FlyLab chose** — the Hill engagement, a selectivity
>    index, the row normalisation of the rate engine;
> 4. **an assumption FlyLab asserted and never fitted** — above all the rule
>    that turns receptor engagement into synaptic gain.
>
> A claim card (Lesson 6) says *"exactly one link of the nine-link chain behind
> this number is a measurement of the system being simulated"*. That is the
> sentence to keep in mind while reading everything else.
>
> When you are ready to quote a number, read
> [INTERPRETATION.md](INTERPRETATION.md) first. It is the more important of
> these two documents.

---

## Contents

| # | Lesson | Teaches |
|---|---|---|
| 1 | [Two minutes in the browser](#lesson-1) | what each headline tile is |
| 2 | [The same run locally](#lesson-2) | that CLI, server and browser are one implementation |
| 3 | [A dose-response, and the parameter that stops mattering](#lesson-3) | why the EC50 is irrelevant at a saturating dose |
| 4 | [Two compounds whose selectivity orderings disagree](#lesson-4) | that a receptor ratio does not buy a circuit window |
| 5 | [Did the connectome matter?](#lesson-5) | the three verdicts, and that a verdict belongs to a cut |
| 6 | [A batch run from a spec file](#lesson-6) | a run directory and the claim card it writes |
| 7 | [Bringing your own compound](#lesson-7) | typed evidence, and being refused |
| — | [Troubleshooting](#troubleshooting) | slow first load, the atlas, the slow tests |

Install, if you want the local lessons:

```bash
python -m pip install -e ".[dev,viz]"
```

Lesson 1 needs nothing at all.

---

<a id="lesson-1"></a>
## Lesson 1 — Two minutes in the browser

Open **<https://minh.systems/FlyLab/>**. Nothing is installed and nothing is
uploaded: the page downloads a Python runtime and the FlyLab wheel and runs the
same scientific core in your tab.

The first load takes about ten seconds on a warm cache and a good deal longer on
a cold one — see [Troubleshooting](#troubleshooting). When the status bar at the
bottom says `dashboard ready`, you are looking at **imidacloprid at 1 µM on the
`named` MaleCNS cut**, which is the default.

### The header

```
Imidacloprid · 1.00 µM
class neonicotinoid   target nAChR_beta1 (agonist)   CAS 138261-41-3   graph named
Evidence coverage: 4 of 6 receptor rows are literature-sourced, 2 are not
modelled (no value exists — they are not zero).
```

Three standing warnings sit underneath it, on every run:

```
! Every number on this page is simulated. FlyLab has never dosed a fly.
! Engagement is not occupancy: only a Kd/Ki row reports fractional receptor
  occupancy, an EC50/IC50 row reports normalised functional engagement.
! A missing value means 'not modelled', never zero.
```

### The four headline tiles

| Tile | Value | Sub-label |
|---|---|---|
| Insect engagement | **1.000** | nAChR_beta1 · Kd |
| Vertebrate engagement | **0.091** | nAChR_a4b2 · EC50 |
| Receptor selectivity | **500×** | nAChR pair · potency ratio |
| Evidence tier | literature order | 4 sourced / 2 not modelled |

- **Insect engagement** — the Hill response of the strongest sourced insect
  receptor row at this dose.
- **Vertebrate engagement** — the same for the vertebrate counterpart. The
  vertebrate side is *scored only*; it never patches a circuit.
- **Receptor selectivity** — the ratio of the two potencies. It is labelled
  `MODEL-DERIVED` and captioned *"ratio of potencies — not a safety margin"*.
- **Evidence tier** — a census, not a score: how many rows carry a sourced value.

### Move the concentration slider

Click **10 nM** on the dose ladder. Every card recomputes. The insect tile goes
from `1.000` to `0.992`; the vertebrate tile falls from `0.091` to `0.001`.
Notice that the insect tile has barely moved over two decades — hold on to
that, it is Lesson 3.

### Open a provenance drawer

Click the green `LITERATURE` chip on the **Insect engagement** tile. A drawer
opens on the right with the whole provenance of that one number:

```
insect_nAChR_beta1                                          LITERATURE
PARAMETER TYPE          Kd
VALUE                   8.300e-11 M
ENGAGEMENT AT THIS DOSE 0.9999
ENGAGEMENT MODEL        binding_engagement_proxy
WHAT THAT MEANS         engagement proxy derived from a binding constant (Kd/Ki)
                        measured at evidence distance E1/E2 -- another species,
                        or a related preparation: the same Hill expression, but
                        it is NOT fractional receptor occupancy of the modelled
                        target
EVIDENCE DISTANCE       E1 (cross-species)
PROVENANCE WARNING      EVIDENCE DISTANCE E1 (cross-species): the binding
                        constant was measured in Myzus persicae (native
                        beta1-containing nAChR, susceptible clone 4106A), not in
                        Drosophila melanogaster. ...
SPECIES / PREPARATION   Myzus persicae (native beta1-containing nAChR,
                        susceptible clone 4106A)
SOURCE                  Bass et al. 2011, BMC Neuroscience 12:51 ...
```

Chips that carry no provenance record — the `NOT MODELLED` chip, the legend
chips — are still clickable and explain **the label itself**.

There is a **Guided tour** button in the top bar and a **Take the guided tour**
button beside the dashboard heading. It runs once on a first visit, remembers
that you dismissed it, never blocks the page, and states one interpretation
caveat per panel. Beside it, **How to read this →** links to
[INTERPRETATION.md](INTERPRETATION.md).

### What this tells you, and what it does not

**Tells you:** what compound, at what free concentration, on which cut; how much
of the receptor row is sourced and how much is missing; and, for every number,
the paper it came from.

**Does not tell you:**

- The insect tile reading `1.000` is **not** "the receptor is 100% occupied in a
  fly". Its own sub-label says `Kd` and its drawer says
  `binding_engagement_proxy, E1 (cross-species)`. It is the Hill expression on
  an aphid binding constant. See [INTERPRETATION §1](INTERPRETATION.md#1-engagement-is-not-occupancy).
- "4 of 6 rows sourced" is coverage, not quality. Two of the four are E2/E3
  extrapolations.
- **`500×` is not a safety margin.** It is a ratio of two literature-order
  teaching numbers, and the vertebrate one is a class-level extrapolation.

> **One thing on this page contradicts the rest of it.** The generated *"Why
> this happened"* paragraph ends with *"…the exact MaleCNS wiring, which every
> structure-preserving shuffle **reproduces**"*. The verdict table directly above
> it is more careful: it says `equivalent within the prespecified margin for
> rewire_degree_preserving, weight_permute; merely indeterminate (not
> distinguishable, not shown equivalent) for sign_permute`. The table is right
> and the sentence is not: "reproduces" is licensed only by an
> `equivalent_within_tolerance` verdict, and `sign_permute` did not earn it.
> This is [referee item B5](../papers/reviews/round2_referee_report.md) in the
> wild. Read the table, not the sentence, and see
> [INTERPRETATION §3](INTERPRETATION.md#3-the-three-dependence-verdicts).

---

<a id="lesson-2"></a>
## Lesson 2 — The same run locally

There are three front ends and one implementation. Prove it.

### The CLI

```bash
flylab assay-subgraph --compound imidacloprid --conc 1e-6
```

```
malecns_neighborhood  compound=imidacloprid  conc=1e-06
  n_nodes                        1126
  n_edges                        1360
  mean_hz                       0.458
  max_hz                       39.849
  mn9_hz                       39.406
  dnp01_hz                     39.844
  gains g_ach=0.050  g_gaba=1.000  g_glu=1.000  g_oct=1.000  g_nav=1.000  ach_tone=1.000
  ! Hops-limited MaleCNS neighborhood, not the full 25M-edge CNS.
  ! Signs from predicted transmitters. ACh gain and RDL gain are teaching patches.
```

### The local server

```bash
flylab serve --port 8799
```

```
FlyLab bench -> http://127.0.0.1:8799
INFO:     Started server process [25764]
INFO:     Waiting for application startup.
INFO:     Application startup complete.
INFO:     Uvicorn running on http://127.0.0.1:8799 (Press CTRL+C to quit)
```

That serves the bench UI *and* the HTTP API. From another terminal — the assay
routes take POST, so a bare GET returns `{"detail":"Method Not Allowed"}`:

```bash
curl -s -X POST -H 'content-type: application/json' \
  -d '{"compound":"imidacloprid","conc_M":1e-6}' \
  http://127.0.0.1:8799/api/assay/subgraph \
| python3 -c "import json,sys; d=json.load(sys.stdin); print(json.dumps({k:d['readouts'][k] for k in ('n_nodes','n_edges','mean_hz','mn9_hz','dnp01_hz')},indent=1))"
```

```json
{
 "n_nodes": 1126,
 "n_edges": 1360,
 "mean_hz": 0.45846211979249263,
 "mn9_hz": 39.40567922340328,
 "dnp01_hz": 39.844466847897166
}
```

Those are the CLI's `0.458`, `39.406` and `39.844` before rounding. The
dashboard route, which is what the bench UI calls, agrees as well:

```bash
curl -s "http://127.0.0.1:8799/api/dashboard?compound=imidacloprid&conc_M=1e-6&graph=named&include_dependence=false" \
| python3 -c "
import json,sys; d=json.load(sys.stdin); h=d['headline']
print('insect', round(h['insect_engagement']['value'],3), h['insect_engagement']['engagement_model'], h['insect_engagement']['evidence_distance_label'])
print('circuit', {k:(round(v['vehicle'],3), round(v['treated'],3), round(v['percent'],1)) for k,v in d['circuit']['readouts'].items()})
"
```

```
insect 1.0 binding_engagement_proxy E1 (cross-species)
circuit {'mean_hz': (6.631, 0.458, -93.1), 'mn9_hz': (38.391, 39.406, 2.6), 'dnp01_hz': (41.616, 39.844, -4.3)}
```

The browser build's tiles in Lesson 1 read `1.000` and `0.091` from the same
call, computed by the same Python inside Pyodide.

### The provenance you need to make two runs comparable

```bash
flylab meta
```

```
FlyLab 0.5.0  notebook schema 0.3
map   male-cns:v1.0
library 0.6  sha256 594b8d32c5e61ebd...  21 compounds
graph named          1126 nodes    1360 edges  seeds: MN9=2, DNp01=2
graph taste_motor    1841 nodes   19066 edges  seeds: MN9=2, DNp01=2, LB1a=11, LB1b=6, LB1c=16, LB1d=5, LB3b=11, LB3c=23
mechanism rules: 15
```

### What this tells you, and what it does not

**Tells you:** the CLI, the FastAPI server and the Pyodide build call one
scientific core; a number obtained at the terminal is the number on the
dashboard. The library SHA and the map id are what make two runs comparable —
quote them with any number you take out of the bench.

**Does not tell you:**

- Identical numbers prove the *transports* agree. They say nothing about whether
  the model is right.
- Parity is tested on representative route families plus an exact route-set
  contract, not on every parameter combination
  ([PAGES.md](PAGES.md)). The browser does **not** regenerate the paper.
- The browser's LIF default window is 200 ms against 500 ms natively (same
  0.1 ms step). A shorter run of the same model is not the same experiment; the
  bench says so in the panel.

---

<a id="lesson-3"></a>
## Lesson 3 — A dose-response, and the parameter that stops mattering

This is the cleanest lesson in the tool, because it shows a parameter that
obviously *should* matter and demonstrably does not.

```bash
for c in 1e-9 1e-8 1e-7 1e-6 1e-5; do
  flylab assay-subgraph --compound imidacloprid --conc $c
done
```

Collecting the interesting lines:

| conc (M) | `g_ach` | `mean_hz` | `mn9_hz` | `dnp01_hz` |
|---|---|---|---|---|
| 1e-9 | 1.010 | 6.704 | 38.475 | 41.684 |
| 1e-8 | 0.974 | 6.437 | 38.175 | 41.439 |
| 1e-7 | 0.129 | 0.951 | 38.605 | 39.651 |
| 1e-6 | **0.050** | **0.458** | **39.406** | **39.844** |
| 1e-5 | **0.050** | **0.458** | **39.406** | **39.844** |

1 µM and 10 µM agree to every printed digit. The receptor is full and the gain
patch has reached its floor (`_FLOOR = 0.05` in `flylab/pharm/mechanisms.py` —
"gains never reach zero: a dead network is not informative" — and respected by
every member of the robustness family). Everything between 10 nM and 100 nM, by
contrast, is a cliff.

Now perturb the parameters one at a time at the saturating dose:

```bash
flylab sensitivity --compound imidacloprid --conc 1e-6
```

```
subgraph  imidacloprid @ 1.00e-06 M  readout=mean_hz  x/2.0
parameter                   low       base       high       span
gain_coef                 2.698      0.458      0.458      2.239
drive_hz                  0.229      0.458      0.917      0.688
weight_threshold          0.458      0.458      0.260      0.199
ec50                      0.458      0.458      0.458      0.000
hill_n                    0.458      0.458      0.458      0.000
```

**A two-fold error in the cited potency moves the readout by 0.000 Hz.** The
asserted gain coefficient moves it by 2.24 Hz. The published pharmacology — the
thing a newcomer assumes is the most important input to a pharmacology model —
is at this dose entirely inert, and a number FlyLab made up dominates.

That ordering is exactly what the global uncertainty budget finds. From
`flylab voi imidacloprid --conc 1e-6 --n-base 64`:

```
  1. gain_transform       variance removed   0.516  calibration of receptor engagement against synaptic gain
  ...
  6. lif_seed             variance removed   0.075  none: increase the number of simulated replicates
  7. expression           variance removed   0.075  FISH / scRNA-seq for receptor expression in the named cell class
  8. potency              variance removed   0.069  electrophysiological dose-response on the stated receptor
```

`lif_seed` is a **deliberately null factor**: it has no effect at all on the
deterministic rate engine, so its index is a direct read of the estimator's own
error. `potency` ranks *below* it. Only `gain_transform` has a confidence
interval that excludes zero; everything else is reported as *unresolved at this
sample size*.

### What this tells you, and what it does not

**Tells you:** at a saturating dose, the accuracy of the cited EC50/Kd is not
what limits this prediction. The shape of the engagement-to-gain rule is. That
is a good reason to read [INTERPRETATION §6](INTERPRETATION.md#6-specification-stability)
before quoting any suppression result.

**Does not tell you:**

- The EC50 is irrelevant **at this dose**. Between 10 nM and 100 nM it is the
  steepest thing in the model. "Insensitive to potency" is a statement about one
  point on one curve.
- `span = 0.000` is a property of the readout `mean_hz`. Another readout, or
  another compound whose dose sits on the cliff, will look different.
- The VOI numbers above are variance of a **simulated** readout under **assumed
  input ranges**, not biological variability, and the run prints
  `! not converged` at this sample size. See
  [INTERPRETATION §7](INTERPRETATION.md#7-the-uncertainty-budget-and-value-of-information).

---

<a id="lesson-4"></a>
## Lesson 4 — Two compounds whose selectivity orderings disagree

```bash
flylab compare imidacloprid deltamethrin --conc 1e-6 --no-occupancy
```

```
compare @ 1.00e-06 M  (estimated 2.26 s)
compound         target                  insect    vert   recSI   cirSI     gap  circuit d%  topology
Imidacloprid     insect_nAChR_beta1       1.000   0.091    2.70    1.79   -0.91       -93.1  composition-dominated
Deltamethrin     insect_Nav               0.959   0.167    2.07    2.28    0.22       560.4  composition-dominated
  * The highest receptor selectivity (Imidacloprid, 2.70 log10) is NOT the highest circuit selectivity (Deltamethrin, 2.28 log10).
  ! Receptor selectivity and circuit selectivity are different quantities; a large receptor ratio does not buy a wide circuit window.
  ! Eight compounds in the library have no circuit selectivity index at all: RDL / GluCl block cannot reach the effect threshold on these cuts.
  ! Topology dependence is a documented label from permutation nulls at the stated n, not a hypothesis test with multiplicity control.
```

Imidacloprid wins on receptor selectivity (2.70 against 2.07) and **loses** on
circuit selectivity (1.79 against 2.28). The `gap` column is the difference, and
it has the opposite sign for the two compounds.

The two indices are different quantities. `recSI` is a ratio of potencies at two
receptors. `cirSI` asks how far *this simulation's circuit* moves before the
vertebrate receptor in the teaching library fills. Nothing guarantees they agree,
and here they do not.

### What this tells you, and what it does not

**Tells you:** a receptor-level selectivity ratio does not predict a circuit-level
window. If you are triaging compounds, you have to say which one you are ranking
on.

**Does not tell you:**

- Neither index is a safety margin, a therapeutic index, or a claim about
  vertebrates. The vertebrate scorecard is scored, never simulated.
- Deltamethrin's `+560.4%` is a change in a simulated rate on a 1126-cell cut,
  not a behavioural effect.
- Eight library compounds have **no** `cirSI` at all — a property of these cuts,
  not evidence about those compounds.
- The `topology` column here runs at the CLI default of 20 shuffles, whose
  permutation resolution is 1/21 ≈ 0.048. It is a label, not a test. Lesson 5
  does it properly.

---

<a id="lesson-5"></a>
## Lesson 5 — Did the connectome matter?

The instrument asks: *would a degraded version of this graph have given the same
drug effect?* It destroys one kind of structure at a time and compares.

```bash
flylab dependence imidacloprid --conc 1e-6 --n 200 --graph named
```

```
imidacloprid @ 1.00e-06 M  -  specific wiring evidence: weak (composition-dominated)
  real effect     -6.173 (mean_hz)
  sign_permute               p=0.2587 (resolution 0.0050)  does not beat null   z=-1.23
  weight_permute             p=0.5572 (resolution 0.0050)  does not beat null   z=-0.49
  rewire_degree_preserving   p=0.4179 (resolution 0.0050)  does not beat null   z=-0.74
  erdos_renyi                p=0.0050 (resolution 0.0050)  beats null   z=-47.23
  necessary information level: degree_sequence - a graph model that keeps only every
  node's in/out degree and its transmitter is equivalent to the real cut on this
  readout within the prespecified margin delta = 0.3316 (gap 0.1293), so nothing more
  detailed is necessary for this prediction.
```

### Reading the ladder

The four null models are ordered by how much real network information they keep.
`erdos_renyi` keeps almost nothing and is trivially distinguishable — and its
`p = 0.0050` is sitting exactly on the **resolution floor** `1/(n+1)`, so report
it as `p ≤ floor`, never as an exact value. The three structure-preserving nulls
are not distinguishable.

That is **not** the same as saying they gave the same effect. The CLI prints
`p` and `z`; the payload carries the three-way verdict per mode. Run it at the
paper's `n = 1000` and read the verdicts:

```bash
flylab dependence imidacloprid --conc 1e-6 --n 1000 --graph named --json
```

| mode | p | gap from null median | verdict |
|---|---|---|---|
| `sign_permute` | 0.2747 | 2.652 | **indeterminate** |
| `weight_permute` | 0.5864 | 0.126 | equivalent_within_tolerance |
| `rewire_degree_preserving` | 0.4715 | 0.123 | equivalent_within_tolerance |
| `sign_permute_weight_matched` | 0.3886 | 0.236 | equivalent_within_tolerance |
| `erdos_renyi` | 0.0010 | 6.015 | distinguishable |

(`delta = 0.3316` — five per cent of the vehicle `mean_hz` of 6.631.)
`sign_permute` has a `p` in the middle of the pack and is *still* not
equivalent: its gap from the null median is eight times the margin. The
permutation probability alone could not have told you which of the two
non-rejection answers you had — which is exactly why there are three verdicts
and not two.

### Now run the same compound on the other cut

```bash
flylab dependence imidacloprid --conc 1e-6 --n 200 --graph taste_motor
```

```
imidacloprid @ 1.00e-06 M  -  specific wiring evidence: present (topology-dependent)
  real effect    -12.231 (mean_hz)
  sign_permute               p=0.0796 (resolution 0.0050)  does not beat null   z=-1.28
  weight_permute             p=0.0746 (resolution 0.0050)  does not beat null   z=1.75
  rewire_degree_preserving   p=0.0050 (resolution 0.0050)  beats null   z=-24.27
  erdos_renyi                p=0.0050 (resolution 0.0050)  beats null   z=-25.43
  necessary information level: topology_without_weight_pairing - ...
```

Same compound. Same dose. Same readout. **Opposite verdict.**

Why: the `named` cut is an *in-star*.

```bash
python -c "from flylab.analysis.dependence import cut_census; print(cut_census('named'))"
```

It has **mean degree 1.21**: 84.9% of its 1360 edges terminate on four seed
cells, only 184 of its 1126 nodes have any input at all, and only 24.9% of its
edges can be part of a path longer than one hop. A degree-preserving rewire of
such a graph is close to the identity — there is almost nothing for the shuffle
to destroy. `taste_motor` has 19 066 edges, mean degree 10.4, and 89% of nodes
with input.

### And it moves again with size

The repository also carries a scale ladder — the same connectome cut at 1000,
5000, 10 000, 25 000 and 50 000 cells:

```bash
python -c "from flylab.analysis.scale import available_cuts; [print(c['name'], c['n_nodes'], c['n_edges']) for c in available_cuts()]"
```

```
scale_1k 1000 22857
named 1126 1360
taste_motor 1841 19066
scale_5k 5000 279845
scale_10k 10000 621600
scale_25k 25000 1364375
scale_50k 50000 2216881
```

Walking a dependence profile up that ladder moves the verdict a second time: at
about 1100 cells imidacloprid is composition-dominated, and at about 5000 cells
the same compound at the same concentration is distinguishable from *every*
null, including the weight-matched transmitter null. That measurement is in
flight at the time of writing — `papers/results.json` still records
`scale_study_status: "not yet run"` — so take the number from the shipped record
rather than from this page.

### What this tells you, and what it does not

**Tells you:** whether *this* readout, for *this* compound, at *this* dose, on
*this* cut, of *this size*, can be told apart from each degraded graph model at
the permutation effort you paid for.

**Does not tell you:**

- **Neither verdict is the truth about imidacloprid.** They are two true
  statements about a readout on two different substrates, and a third substrate
  gives a third answer. The interesting question is not which cut is right; it
  is *where the verdict settles as the cut grows* — and that is an open question
  currently being measured, not something either of these runs answers.
- Nothing about the connectome, unless you name the cut *and* its size and mean
  degree. A cut where most edges land on a handful of cells cannot support a
  negative topology result at all.
- A large `p` licenses "not distinguishable from this null ensemble at n
  shuffles" and nothing stronger. Equivalence needs the prespecified margin, and
  the result says when it was met.
- `z` is not a probability. Null distributions on degraded graphs are routinely
  non-normal.
- A run at `n = 200` is **exploratory**; the output says so, and names the two
  prespecified confirmatory compounds. See
  [INTERPRETATION §5](INTERPRETATION.md#5-confirmatory-versus-exploratory).
- Watch for `! the permutation p had not stabilised at this n`. The verdict at
  alpha can be settled while the value of `p` is not yet quotable.

**Never quote a dependence verdict without naming the compound, the
concentration, the readout and the cut.**

---

<a id="lesson-6"></a>
## Lesson 6 — A batch run from a spec file

A spec is one YAML file that fully describes a run, so reusing FlyLab does not
mean learning its Python API.

```bash
flylab spec-schema --example > tutorial.yaml
```

Edit it down to something quick:

```yaml
name: tutorial_fipronil
compound: fipronil
concentrations: [1.0e-8, 1.0e-7, 1.0e-6]
graph: named
readouts: [mean_hz, mn9_hz]
engines: [rate]
analyses:
  dependence: {permutations: 200}
figures: false
```

Validate before you spend anything:

```bash
flylab run tutorial.yaml --dry-run
```

```
resolved spec (tutorial.yaml)
flylab_spec_version: '1.0'
name: tutorial_fipronil
compounds:
- fipronil
concentrations:
- 1.0e-08
- 1.0e-07
- 1.0e-06
graph: named
engines:
- rate
readouts:
- mean_hz
- mn9_hz
replicates: 1
seed: 0
include_vehicle: true
jitter_log10: 0.0
analyses:
  dependence:
    permutations: 200
figures: false
cards: true

estimate: ~10 s (4 assay runs + dependence)
  dependence   ~9.4 s
would write -> runs/tutorial_fipronil/  (cards/, dependence.json, manifest.json, notebooks/, results.csv, spec.yaml)
spec_sha256 aa93e612aeb728780552cabcd3655eaee69249d780d3699962f57e3bd277dff0
  ! Nothing was executed: --dry-run validates the spec and prints the plan.
```

Then run it:

```bash
flylab run tutorial.yaml --outdir runs/tutorial_fipronil
```

```
  assays (rate)
    rate: 3 rows + 1 vehicle rows
  dependence
    dependence: fipronil @ 1.00e-06 M, n=200
  -> runs/tutorial_fipronil  (10 files, digest 4280c4fd2215)
run -> runs/tutorial_fipronil
  rows      4
  notebooks 4
  analyses  dependence
  cards     1
  digest    4280c4fd22150c6f688cbf7adf2c3858a2560225a42357e48af026b01762f6d2
  total     21.821 s
  ! Every value is simulated: teaching EC50 library, hops-limited map, gain patches. No row is a measurement from a living fly.
  ! no figures were written: figures were switched off in the spec
```

The run directory is self-describing:

```
runs/tutorial_fipronil/
  cards/fipronil__rate.json
  cards/fipronil__rate.md
  dependence.json
  manifest.json
  notebooks/rate__fipronil__1.00e-6__r0.json
  notebooks/rate__fipronil__1.00e-7__r0.json
  notebooks/rate__fipronil__1.00e-8__r0.json
  notebooks/rate__vehicle__0__r0.json
  results.csv
  spec.yaml
```

```csv
engine,condition,compound,conc_M,replicate,mean_hz,mn9_hz
rate,drug,fipronil,1e-08,0,6.7996293607879315,40.69386434297126
rate,drug,fipronil,1e-07,0,7.369266774566889,49.4074001737001
rate,drug,fipronil,1e-06,0,7.530851909481498,52.139965246425675
rate,vehicle,,0.0,0,6.631244285926711,38.39073690731564
```

### The claim card

`cards/fipronil__rate.md` has nine sections, always in the same order, always all
present. Its first three:

```markdown
# Claim card -- fipronil at 1.00e-06 M

## Claim
> On the MaleCNS v1.0 `named` cut (1126 cells, 1360 edges, >= 5 synapses),
> FlyLab's rate engine predicts that fipronil at 1.00e-06 M raises mean_hz from
> 6.631 under vehicle to 7.531 (+13.6%).

## Status
This is a prediction of a simulation, not an observation. ... the notebook's
live_lab field is null.

## Evidence inputs
2 of 4 receptor rows carry a sourced, typed parameter. ...

| receptor | parameter type | value (M) | evidence distance | engagement |
|---|---|---|---|---|
| insect_nAChR | unknown | not modelled | unsupported | N/A |
| vertebrate_nAChR_a4b2 | unknown | not modelled | unsupported | N/A |
| vertebrate_GABA_A | IC50 | 1.10e-06 | exact_compound_exact_receptor_other_species | 0.476 |
| insect_RDL | IC50 | 3.00e-08 | exact_compound_exact_receptor_exact_species | 0.985 |
```

and section 7, because this spec did not ask for a robustness analysis:

```markdown
## Specification stability
Specification stability was not assessed in this run. Nothing here says whether
this claim survives the alternative engagement-to-gain rules.

*How to assess it: add `robustness: {specifications: default_family}` to the
spec's `analyses:` block.*
```

A full section-by-section reading of this card is in
[INTERPRETATION §13](INTERPRETATION.md#13-a-claim-card-read-section-by-section).

You can also print a card for any compound without a run:

```bash
flylab card fipronil --conc 1e-6 --markdown
```

### What this tells you, and what it does not

**Tells you:** the spec, the resolved defaults, the seeds, the code version, the
library hash, the map ids and a content hash of every artifact. The same spec
and seed produce a byte-identical `results.csv` and the same
`manifest.determinism.digest`; `tests/test_spec.py` proves it by running the
same spec twice.

**Does not tell you:**

- Sections 6 and 7 read *"not assessed in this run"* when the spec did not ask
  for them. That is an unanswered question, not an absent concern — and it is
  why the sections are never dropped.
- Determinism is per environment. Hashes are not portable across NumPy or
  matplotlib versions; `requirements-lock.txt` records which versions produced
  the committed numbers.
- The LIF engine consumes the RNG, so a different `seed` gives different spike
  trains. The rate engine does not, which is why a rate run's readouts do not
  move with the seed.

---

<a id="lesson-7"></a>
## Lesson 7 — Bringing your own compound

A library row is admitted on two facts and two facts only: what kind of
parameter the source reported (`param_type`) and how far that source sits from
"this compound, at this receptor, in *Drosophila melanogaster*" (`relation`).
Everything else — including the prose in `source` — is documentation.

Run `flylab tutorial 7` for this lesson end to end. Here it is written out.

### An honest row

```python
from flylab.pharm.occupancy import compare_compound, load_library

lib = load_library()
lib["compounds"]["my_neonic"] = {
    "name": "My neonicotinoid",
    "class": "neonicotinoid",
    "receptors": {
        "insect_nAChR": {
            "param_type": "EC50",                 # what the source measured
            "value_M": 4.0e-8,
            "n": 1.1,
            "direction": "agonist",
            "relation": "exact_compound_exact_receptor_exact_species",
            "species": "Drosophila melanogaster",
            "source": "my lab, two-electrode voltage clamp on Dalpha1/Dbeta1",
            "evidence_tier": "literature_order",
        },
        "vertebrate_nAChR_a4b2": {               # nobody measured this
            "param_type": "unknown",
            "value_M": None,
            "direction": "none",
            "relation": "unsupported",
            "source": "no source measured this compound at this receptor",
            "evidence_tier": "class_placeholder",
        },
    },
}
result = compare_compound("my_neonic", 1e-6, library=lib)
for row in result["receptors"]:
    print(row["receptor"], row["engagement"], row["engagement_model"], row["evidence_distance"])
```

```
receptor                 engagement  model                        distance
insect_nAChR                  0.972  functional_engagement        E0
vertebrate_nAChR_a4b2           n/a  not_modelled                 E4
```

Your EC50, measured on target, gives `functional_engagement` at `E0` — the
strongest a *functional* potency can ever be. It is a normalised functional
response, **not** an occupancy. The receptor you left unmeasured reports `n/a`
and is excluded from every number downstream; it does not become a small
plausible effect.

### The same row, in a library file

For anything you intend to keep, copy `flylab/pharm/library.yaml`, add the
compound under `compounds:`, and load your copy by path. The YAML is the same
four facts:

```yaml
  my_neonic:
    name: My neonicotinoid
    class: neonicotinoid
    receptors:
      insect_nAChR:
        param_type: EC50                                  # what was measured
        value_M: 4.0e-8                                   # in molar
        n: 1.1
        direction: agonist
        relation: exact_compound_exact_receptor_exact_species   # -> distance E0
        species: Drosophila melanogaster
        source: >-
          my lab, unpublished: two-electrode voltage clamp on Dalpha1/Dbeta1,
          agonist EC50 40 nM, n = 4 oocytes.
        evidence_tier: literature_order
      vertebrate_nAChR_a4b2:
        param_type: unknown                               # no usable number
        value_M: null                                     # NOT 0, NOT a guess
        direction: none
        relation: unsupported                             # -> distance E4
        source: no source measured this compound at this receptor
        evidence_tier: class_placeholder
```

```python
from pathlib import Path
lib = load_library(Path("my_library.yaml"))
```

Three fields carry all the weight, and it is worth being pedantic about each:

- **`param_type`** — what the instrument produced, not what you wish it had.
  A displacement IC50 is an `IC50`, even if the paper's abstract calls it an
  affinity.
- **`relation`** — how far the preparation is from *this* receptor in
  *Drosophila*. A hybrid or chimeric receptor is
  `exact_compound_related_receptor` (E2), not `..._exact_species`.
- **`source`** — enough for a reader to find the number and disagree with it:
  citation, preparation, what was measured, and anything you had to choose. The
  shipped rows read like short methods paragraphs, and that is the standard.

`evidence_tier` records provenance class; `measured_fit` is reserved for values
fitted to data you imported and nothing in the shipped library uses it.

### Now type it wrongly

`flylab/pharm/evidence.py::check_transformation(param_type, model, relation)`
is the gate. Ask it for more than the evidence supports and it refuses, with the
reason:

```python
from flylab.pharm import evidence as ev

ev.check_transformation("EC50", "binding_occupancy",
                        "exact_compound_exact_receptor_exact_species")
```

**An EC50 relabelled as an occupancy:**

```
EvidenceTypeError: EC50 is not a binding constant: half-maximal effective
concentration of a functional response; the Hill curve is a normalised response,
NOT occupancy. A Hill curve built from a EC50 is a normalised functional
response (engagement_model=functional_engagement), not physical receptor
occupancy and not a binding-derived proxy. Use a Kd or Ki row if you need a
binding-based model.
```

**A binding constant measured in another species, claimed as occupancy:**

```
EvidenceTypeError: Kd at evidence distance E1 (cross-species) cannot drive
engagement_model=binding_occupancy: cross-species: the right receptor, another
organism. Physical fractional occupancy of the modelled receptor requires
evidence distance E0 (on-target), i.e.
relation=exact_compound_exact_receptor_exact_species; this row supports at most
engagement_model=binding_engagement_proxy, which is a proxy, not an occupancy of
this receptor.
```

**A placeholder asked for a number:**

```
EvidenceTypeError: unknown at evidence distance E4 (unsupported) cannot drive
engagement_model=functional_engagement: unsupported: no source supports any
value here. A row of this kind reports N/A (not modelled), never a number.
```

The vocabulary is closed, too — you cannot invent a parameter type or a relation
to get around the table:

```
EvidenceTypeError: unknown parameter type 'Ki50'; allowed: Kd, Ki, EC50, IC50,
Kb, relative_potency, class_order, unknown

EvidenceTypeError: unknown source relation 'measured_in_our_lab'; allowed:
exact_compound_exact_receptor_exact_species,
exact_compound_exact_receptor_other_species, exact_compound_related_receptor,
class_extrapolation, unsupported
```

A misspelt **compound** in a spec is caught the same way, and the message names
every valid option:

```
$ flylab run bad.yaml
bad.yaml: unknown compound 'fibronil'. Valid options: acetamiprid, acetylcholine,
caffeine, chlordimeform, chlorpyrifos_oxon, clothianidin, ddt, deltamethrin,
diazepam, dieldrin, fipronil, gaba, imidacloprid, ivermectin, nicotine,
nitenpyram, permethrin, picrotoxin, spinosad, sulfoxaflor, thiamethoxam.
```

### What this tells you, and what it does not

**Tells you:** exactly which transformation your evidence licenses, and a
refusal with a reason when you ask for more. The full table is
`TRANSFORMATION_TABLE` in `flylab/pharm/evidence.py`.

**Does not tell you:**

- The type system checks that you *described* your measurement honestly. It
  cannot check that the measurement is right, or that the receptor key you chose
  is the one your preparation actually contained.
- Declining to claim as much as your evidence permits is always allowed. Only
  claiming *more* is an error.
- Adding a compound does not add a mechanism. The gain rules in
  `flylab/pharm/mechanisms.py` are keyed on receptors, and a row at a receptor
  with no rule changes no circuit.
- A row you add is not in `flylab/pharm/library.yaml`, so it is not in the
  library SHA a notebook records. Keep your own library file, and cite it.

---

## Troubleshooting

**The browser bench takes a long time to load the first time.**
The page fetches a Python runtime (Pyodide), NumPy, PyYAML and the FlyLab wheel
— around 20 MB — and then installs them in your tab. A cold cache is measured in
tens of seconds; a warm one is about ten. The status bar reports the boot time
when it finishes. Performance depends on browser, CPU and cache state, so do not
treat any single timing as a specification ([PAGES.md](PAGES.md)).

**Does the browser need the network once it has loaded?**
No. The published build vendors Pyodide and the chart libraries, and
`manifest.json` records both decisions as `pyodide_vendored` and `js_vendored`.
A build served locally with every off-origin request blocked boots, runs assays
and draws charts.

**"Plotly did not load from either CDN; charts fall back to tables."**
Harmless. Every chart in the bench has a table beside it carrying the same
values. It means the two CDNs were unreachable, which is normal behind a
restrictive proxy; the served build fetches them from a CDN while the static
build vendors them.

**Do I need to download the MaleCNS atlas?**
No. The committed cuts under `data/derived/` are what every routine command
uses, and they ship with the repository. `flylab download-malecns` and
`flylab extract-subgraph` exist for rebuilding a cut from the full upstream
reconstruction, which is not committed. Skip both unless you are making a new
cut.

**`python -m pytest -q` skips things.**
`pyproject.toml` sets `addopts = "-m 'not slow'"`. The slow marks are the wheel
build, the static (Pyodide) site and the browser boot:

```bash
python -m pytest -q            # the fast contract suite
python -m pytest -q -m slow    # the wheel, the static build, the browser boot
```

Run the slow set before a tagged release; CI runs it weekly and on release tags.

**An analysis is slower than its estimate.**
The estimates printed by `--estimate` / `--dry-run` are order-of-magnitude, in
process, on an unloaded machine, and they say so. Permutation cost is linear in
`n` and in the number of cells: `flylab/analysis/dependence.py` records the
84-cell landscape on the `named` cut at about 11 s for `n = 20`, 55 s for
`n = 100` and 9 minutes for `n = 1000`, or roughly a quarter of that with four
worker processes (`analyses: {dependence: {n_jobs: 4}}` in a spec). Start with
the small `n` and raise it only for the number you intend to quote. If another
process on the machine is busy, every estimate on this page is optimistic.

**`flylab manifest --check` fails after I edited something.**
That is the point. Any edit under `flylab/`, `data/` or `papers/` invalidates
the artifact manifest. Run `flylab manifest --write` and commit the result with
the change that caused it ([WORKFLOW.md](WORKFLOW.md)).

---

## Where to go next

- **[INTERPRETATION.md](INTERPRETATION.md)** — what each output licenses you to
  say. Read this before quoting anything.
- [NOVELTY.md](NOVELTY.md) — the claim guardrail, and the claims this project
  has withdrawn.
- [ARCHITECTURE.md](ARCHITECTURE.md) — the dataflow, the layer boundaries and the
  two cuts.
- [WORKFLOW.md](WORKFLOW.md) — specs, claim cards, the manifest and the release
  gate.
- [PAGES.md](PAGES.md) — the browser build and what transport parity does and
  does not cover.
- [README.md](README.md) — the documentation index and its authority map.
