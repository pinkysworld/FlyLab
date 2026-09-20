# Workflows: specs, claim cards and releases

Three pieces of machinery sit on top of the science layers, and they exist for
one reason: a result should arrive with everything needed to disbelieve it.

| Piece | Module | What it produces |
|---|---|---|
| Experiment spec | `flylab/spec.py` | a self-describing run directory from one YAML file |
| Claim card | `flylab/report/card.py` | nine fixed sections per headline result, JSON and Markdown |
| Artifact manifest | `scripts/manifest.py` | hashes of every input the committed results depend on |

None of them computes pharmacology or circuit maths. They compose the layers
that do, record what was used, and refuse to leave a section out.

---

## 1. The experiment spec

A spec is a YAML file that fully describes a run, so reusing FlyLab does not
mean learning its Python API.

```yaml
name: fipronil_connectome_dependence
compound: fipronil                 # or: compounds: [fipronil, imidacloprid]
concentrations: [1e-8, 1e-7, 1e-6, 1e-5]
graph: named                       # named | taste_motor
readouts: [mean_hz, mn9_hz]
engines: [rate]                    # rate | lif | taste | taste_map
analyses:
  dependence: {permutations: 1000, correction: benjamini-hochberg}
  robustness: {specifications: default_family}
  uncertainty: {samples: 2048}
```

```bash
flylab spec-schema --example > spec.yaml     # a starter file
flylab run spec.yaml --dry-run               # validate + estimate, run nothing
flylab run spec.yaml --outdir runs/fipronil  # the real thing
flylab spec-schema --json                    # the machine-readable schema
```

`flylab experiment run` is a thin alias of `flylab run`, and a v0.5 design file
is a valid spec: `assay`, `concs_M` and `compound` are accepted as aliases of
`engines`, `concentrations` and `compounds`.

### Fields

Every field has a default, and the run directory records the spec *after* those
defaults were filled in. `flylab spec-schema --json` is the authoritative list;
the summary:

| Field | Default | Notes |
|---|---|---|
| `name` | `flylab_run` | names the run directory and the manifest |
| `compounds` | `[imidacloprid]` | must exist in the typed library |
| `concentrations` | five decades | free molar; sorted and de-duplicated |
| `graph` | `named` | `named` or `taste_motor` |
| `engines` | `[rate]` | `rate`, `lif`, `taste`, `taste_map` |
| `readouts` | all five | `mn9_hz`, `dnp01_hz`, `mean_hz`, `g_ach`, `g_gaba` |
| `replicates` | `1` | replicate *r* runs with `seed + r` |
| `seed` | `0` | the rate engine is deterministic; the LIF engine is not |
| `include_vehicle` | `true` | a `compound: null` cell per engine |
| `jitter_log10` | `0.0` | > 0 turns on the library Monte-Carlo |
| `analyses` | `{}` | see below; omit to run the assays only |
| `figures` / `cards` | `true` | figures need matplotlib; cards never do |

Anything else is rejected, and the message names the valid options:

```
$ flylab run bad.yaml
bad.yaml: unknown compound 'fibronil'. Valid options: acetamiprid, acetylcholine,
caffeine, ... spinosad, sulfoxaflor, thiamethoxam.  Did you mean 'fipronil'?
```

The same applies to unknown keys, graphs, engines, readouts, analyses and
analysis options.

### Analyses

| Analysis | Options | Artifact |
|---|---|---|
| `dependence` | `permutations`, `correction`, `modes`, `readout`, `concentrations`, `alpha`, `n_jobs` | `dependence.json` |
| `robustness` | `specifications` (`default_family`/`fast_family`), `shuffles`, `n_jobs` | `robustness.json` |
| `uncertainty` | `samples`, `readout`, `bootstrap`, `concentrations` | `uncertainty.json` |

`concentrations` defaults to `headline` — the highest dose in the spec — because
a permutation analysis at every dose multiplies the cost by the length of the
ladder. Pass `all`, or an explicit list, when you want the whole ladder.

`correction: benjamini-hochberg` applies BH across every (cell, null model)
permutation *p* in the run and stores the result next to the uncorrected
per-mode values; it never overwrites them.

The analysis entry points are looked up by name at call time, never imported at
module scope. If one is renamed, the run says so instead of crashing:

```
analysis 'dependence': none of dependence_profile is exported by
flylab.analysis.dependence. It currently exports: classify, dependence_landscape, ...
Update flylab.spec.ANALYSES if the entry point was renamed.
```

### The run directory

```
runs/fipronil_connectome_dependence/
  manifest.json      spec hash, code version, git sha, library sha, map ids, seeds,
                     timings, and a hash of every artifact below
  spec.yaml          the spec as resolved, defaults included, with its sha256
  notebooks/         one notebook JSON per assay run, warnings and live_lab: null intact
  results.csv        engine, condition, compound, conc_M, replicate, readouts
  dependence.json    only if requested
  uncertainty.json   only if requested
  robustness.json    only if requested
  cards/             one claim card per headline result, as .json and .md
  figures/           only when matplotlib is importable
```

Without matplotlib the run still completes; `manifest.figures.available` is
`false`, the reason is recorded, and a warning says so.

### Determinism

The same spec and the same seed produce a byte-identical `results.csv` and
identical artifact hashes. Two things make that non-trivial, and both are
handled explicitly rather than hidden:

* Notebooks carry `created_utc` and a `platform` string; analysis payloads carry
  `runtime_s`. Those keys are **volatile** (`flylab.spec.VOLATILE_KEYS`).
* Every artifact therefore has two hashes: `sha256` over the bytes on disk and
  `content_sha256` over the same payload with the volatile keys stripped. The
  manifest lists which keys were removed from each file.

`manifest.determinism.digest` is a single hash over the `content_sha256` values.
Two runs of one spec agree on it; `tests/test_spec.py` proves this by running
the same spec twice and comparing.

What is **not** deterministic, and is not claimed to be:

* `manifest.json` itself — it records timings and a timestamp.
* Runs with `jitter_log10 > 0` are reproducible *from the seed*, but the readouts
  are not the unjittered ones; the run warns about this.
* The LIF engine consumes the RNG, so a different `seed` gives different spike
  trains. The rate engine does not, which is why the readouts of a rate run do
  not move with the seed.
* Hashes are not portable across NumPy or matplotlib versions. The digest proves
  reproduction on one environment; `requirements-lock.txt` says which.

---

## 2. The claim card

A card is a per-result artifact with exactly nine sections, always in this
order, always all present:

1. **Claim** — one sentence: compound, dose, graph, engine, readout, vehicle and change.
2. **Status** — model-derived, and what that means.
3. **Evidence inputs** — each receptor row with its parameter type (`Kd`/`IC50`/…) and its evidence distance (the source relation), plus the engagement it produced. An unsupported row shows N/A, never a small number.
4. **Measured, not assumed** — the MaleCNS edges: the one link in the chain that is a measurement of the system being simulated.
5. **Assumptions introduced** — transmitter sign, gain mapping, free concentration, uniform expression, plus everything else the audit chain declares.
6. **Connectome dependence** — which null models the effect is distinguishable from, each with its permutation *p*, the dependence class and the necessary information level.
7. **Specification stability** — how many of the gain-rule family retain each conclusion.
8. **Independent biological validation** — none, and what would change that.
9. **Named unknowns** — the gaps nothing in the chain can fill.

Sections 6 and 7 need analyses a run may not have requested. They are then
filled in with `assessed: false` and the words *"not assessed in this run"* —
never dropped, because a missing section reads as an absent concern rather than
an unanswered question.

```bash
flylab card fipronil --conc 1e-6                 # the summary
flylab card fipronil --conc 1e-6 --markdown      # the full card
flylab card --run runs/fipronil --json           # the cards a run already wrote
```

```python
from flylab.report.card import claim_card, cards_for_run, to_markdown

card = claim_card(notebook, compound="fipronil", conc_M=1e-6, vehicle=vehicle_notebook,
                  dependence=dependence_payload, requested_analyses=["dependence"])
print(to_markdown(card))
cards = cards_for_run("runs/fipronil")   # reads cards/, or rebuilds from notebooks/
```

`claim_card` builds on `flylab.analysis.claims.claim_audit`, which walks the
nine-link dependency chain; the card is the reader-facing projection of that
chain, and `flylab.analysis.claims.evidence_rows` is the typed-evidence
accessor both share.

---

## 3. Manifest, lock and release

### The artifact manifest

`scripts/manifest.py` hashes every class of input the committed results depend
on into `artifact-manifest.json`:

| Section | Covers | Fatal on mismatch |
|---|---|---|
| `code` | every source file git accounts for under `flylab/` and `scripts/`, plus `pyproject.toml` | yes |
| `library` | `flylab/pharm/library.yaml` | yes |
| `graphs` | `data/derived/*.json` | yes |
| `literature` | `data/literature/*.yaml` | yes |
| `paper` | `papers/results.json`, `papers/*.md`, `papers/figures/*`, `papers/tables/*` | yes |
| `dependencies` | resolved versions of the dependency closure | **no** |

The file set comes from `git ls-files --cached --others --exclude-standard`:
everything committed plus everything new that is not ignored. A build artifact,
a `__pycache__` entry or a `*.feather` can therefore never drift into the code
hash, while a module added in a branch is covered before it is committed.

`dependencies` is deliberately not fatal. A manifest is verified on machines
other than the one that wrote it, and a different patch release of NumPy is a
fact to report, not a reason to fail a release. It lands in `notes`.

```bash
flylab manifest --write        # (re)write artifact-manifest.json
flylab manifest --check        # structured diff; exits 2 on a fatal mismatch
flylab manifest --lock         # write requirements-lock.txt
python scripts/manifest.py --check --json
```

`verify_manifest(path)` returns the diff as data: per section, `ok`, `added`,
`removed`, `changed` (each with the committed and the current hash), and a
one-line `summary`. `--check` prints it and exits non-zero when a fatal section
differs.

**Regenerate the manifest with the change that caused it.** Any edit under
`flylab/`, `data/` or `papers/` invalidates it; that is the point. `--check`
fails until `--write` is run and the result committed.

### The dependency lock

`requirements-lock.txt` records the versions the committed results were produced
with: `pip`'s installed set, filtered to the project's actual dependency closure
(the names in `[project] dependencies` and `optional-dependencies`, plus
everything those require transitively).

It is **not** a full environment specification. There are no artifact hashes, no
platform or Python-version pins, no build dependencies, and no guarantee about
packages outside the closure. Use it to reproduce the numbers, not to certify a
supply chain. The file says so in its own header.

### The release gate

`.github/workflows/release.yml` runs on a `v*` tag push and on
`workflow_dispatch`, and fails the release if any of it fails:

1. `pytest -q` — the fast contract suite.
2. `pytest -q -m slow` — the wheel build, the static (Pyodide) site and the browser boot, which ordinary CI skips.
3. `python scripts/reproduce_paper.py --fast --outdir dist/fastrun` — the whole pipeline, into a scratch directory so the committed `papers/` artifacts stay exactly as tagged.
4. `python scripts/build_pages.py --out dist/pages --vendor-pyodide` — the static bench.
5. `flylab manifest --check` — the committed results still follow from the committed inputs.
6. A spec run end to end, to prove the workflow machinery itself still works.

It also diffs the installed versions against `requirements-lock.txt` and prints
the drift, and uploads the fast run, the static manifest, the spec run, the
artifact manifest and the lock as a build artifact.

It creates no tag and no GitHub release. That stays a human act.

`.github/workflows/tests.yml` gained a weekly `slow` job (Mondays, plus
`workflow_dispatch`) for the same reason: `pyproject.toml` sets
`addopts = "-m 'not slow'"`, so between releases the static build and the
browser checks would otherwise never run. Its manifest check is
`continue-on-error`, because `main` moves between releases and a stale manifest
there is a reminder, not a breakage.

### Before tagging

```bash
python scripts/reproduce_paper.py       # the real run, not --fast
python -m pytest -q && python -m pytest -q -m slow
flylab manifest --lock
flylab manifest --write
git add artifact-manifest.json requirements-lock.txt papers/
```

Then tag. `papers/SUBMISSION_CHECKLIST.md` covers what the manuscript needs on
top of this.
