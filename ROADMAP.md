# FlyLab roadmap

Target: a virtual pharmacology lab that a person can actually run, not a vision deck.

Dates are elapsed time from this commit (2026-09-19). Slip a gate rather than fake it.

---

## Product definition

**FlyLab v1 (usable teaching bench)**  
Install → choose drug + dose → occupancy on insect vs vertebrate receptors → run one validated fly circuit → plot named cell rates → export notebook JSON.

**FlyLab v2 (research bench)**  
v1 + three circuits + live-fly protocol pack + one paired wet-lab dataset.

**FlyLab v3 (paper)**  
v2 written for IJRC as architecture + evaluation, with code and data availability statements.

---

## Gate 0 — contract (this week)

**Done when**

- README states what the maps are and are not.
- Occupancy CLI runs without a connectome download.
- Drug library entries have a citation field. No invented EC50.
- `.gitignore` excludes `*.feather`, `*.parquet`, connectome dumps.

**Exit test**

```text
flylab occupancy imidacloprid --conc 1e-6
# prints occupancy table, source citations, and a warning that circuit injection is Phase 1
```

---

## Gate 1 — the taste lock (weeks 1–3)

Do not touch the full 166k graph until this passes.

**Build**

- Loader for a **taste subgraph**: labellar sweet GRNs, bitter GRNs, second-order cells (Usnea / Rattle / Phantom / G2N-1 if present in the chosen map), MN9.
- LIF or rate model with synapse-count weights and predicted-transmitter signs (Shiu 2024 rule).
- Two conditions: sugar only; sugar + bitter.

**Pass / fail (published control, not our opinion)**

On FlyWire-style models, bitter vetoes feeding: MN9 high with sugar, collapsed with sugar+bitter (order-of-magnitude drop). Reproduce the *direction* and a documented numeric window. Publish the exact map version and seed in the notebook JSON.

If this gate fails, stop. Pharmacology on a dead circuit is theatre.

**Hardware:** laptop, 16 GB RAM is enough for the subgraph. No GPU required.

---

## Gate 2 — pharmacology on that circuit (weeks 3–6)

**Drug objects (minimum four)**

| Compound | Why it is in v1 | Insect target | Vertebrate panel |
|---|---|---|---|
| Imidacloprid | selective insect nAChR agonist | cholinergic CNS gain ↑ then desensitization | mammalian nAChR, weak |
| Nitenpyram | fast flea/fly adulticide, teaching contrast | insect nAChR | mammalian nAChR |
| Nicotine | shared ligand, different potency | insect + vertebrate nAChR | human/rat nAChR |
| Picrotoxin or dieldrin-class / fipronil | GABA-Cl block | RDL / GluCl | GABA-A (different) |

Mechanism in v1 is deliberately crude and stated as such:

```text
occupancy  =  conc^n / (EC50^n + conc^n)
effective_weight_ij  *=  (1 + alpha * occupancy * direction)
```

`direction` is + for agonist drive or − for channel block, taken from the presynaptic predicted transmitter. `alpha` is a fitted scalar, not a molecular dynamics result. Every notebook prints that sentence.

**Pass / fail**

- Imidacloprid at a high occupancy collapses or saturates cholinergic feeding drive relative to vehicle.
- Nicotine moves both insect and vertebrate panels; imidacloprid moves insect ≥ vertebrate.
- Changing EC50 in the YAML changes the curve. No hidden constants.

---

## Gate 3 — bench UI (weeks 5–8)

A single local web page (FastAPI + static, or Streamlit if faster):

- left: fly circuit activity (named cells, not 166k dots on day one)
- right: vertebrate receptor occupancy (“next to brain”)
- top: compound, log-concentration slider, Run, Export notebook

**Pass / fail:** a person who is not you can run the sugar/bitter + imidacloprid demo in ten minutes from the README.

---

## Gate 4 — second circuit + vertebrate panel (weeks 8–10)

Add **loom / giant fibre** (LPLC2 / LC4 → DNp01) as assay two.  
Keep the vertebrate panel as a **reduced receptor set**, not a fake mouse connectome:

- muscle nAChR (α1βγδ)
- neuronal α4β2
- GABA-A α1β2γ2
- (optional) 5-HT3

This is the honest “next to brain.” Same clock, two scorecards.

---

## Gate 5 — live fly pack (weeks 10–16, needs a fly room)

Protocols in `/protocols`, each one page:

1. Negative geotaxis (climbing) after vehicle vs compound.
2. Proboscis extension / feeding suppression (pairs with Gate 1).
3. Optional: bang-sensitive seizure line vs one AED if a lab already keeps those stocks.

**Pass / fail:** one paired table, n stated, ethics/institution named. Fit a single `alpha` or EC50. Report residual. Do not claim the sim predicted the fly until that table exists.

If there is no fly room, FlyLab v1 still ships. Gate 5 is marked `blocked: no colony` and the paper stays computational.

---

## Gate 6 — IJRC manuscript

Only after Gate 2 is green.

Structure: related-work table (RatCVS / PharmVR / flypoke / MaleCNS meme sims) → architecture → taste control → four-drug assays → limits → code availability.

Word budget 6–8k. IEEE refs. AI-use statement. Do not title it “we uploaded a fly.”

---

## What we will not do in v1

- Train the connectome to play games.
- Simulate 86 billion human neurons.
- Fit every synapse as a kinetic scheme.
- Promise 3Rs replacement of mammalian GLP studies.
- Commit 1 GB Feather files to git.

---

## Default stack

- Python 3.11+
- numpy / pydantic / typer
- optional: flypoke or a thin MaleCNS CSR loader in Phase 1
- pytest for occupancy math and assay gates
- FastAPI later

Reuse Shiu / flypoke for the spike kernel. Write the pharmacology layer and the lab notebook ourselves. That split *is* the novelty.
