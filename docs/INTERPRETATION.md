# How to read a FlyLab result

FlyLab is easy to drive and easy to misread. Five minutes gets you a dashboard
full of numbers; nothing on that dashboard tells you which of them you are
allowed to repeat out loud.

This page is the answer to that. For each output the bench produces it says
**what it means**, **what it does not mean**, and **the specific misreading to
avoid** — and names the function you can read to check. Every entry is grounded
in shipped behaviour, not in the manuscript.

The companion pages are [TUTORIAL.md](TUTORIAL.md), which walks you through
producing these outputs, and [NOVELTY.md](NOVELTY.md), which is the project's
claim guardrail. Where this page states a rule about what you may write, the
rule comes from NOVELTY.md and that document wins.

> **The reason this page exists.** The project's own second-round referee report
> (`papers/reviews/round2_referee_report.md`) found that the authors had
> misread several of their own results: a non-rejection written up as
> "reproduces", a near-algebraic identity written up as a finding about the
> connectome, a baseline that could only fail written up as a property of the
> substrate, and a noise floor stated wrongly. None of those were dishonest.
> They were all cases of reading a number as a stronger quantity than it is.

---

## Contents

1. [Engagement is not occupancy](#1-engagement-is-not-occupancy)
2. [Not modelled is not zero](#2-not-modelled-is-not-zero)
3. [The three dependence verdicts](#3-the-three-dependence-verdicts)
4. [The verdict belongs to the substrate](#4-the-verdict-belongs-to-the-substrate)
5. [Confirmatory versus exploratory](#5-confirmatory-versus-exploratory)
6. [Specification stability](#6-specification-stability)
7. [The uncertainty budget and value of information](#7-the-uncertainty-budget-and-value-of-information)
8. [What the engine's row normalisation does to composition versus topology](#8-what-the-engines-row-normalisation-does-to-composition-versus-topology)
9. [The ablation ladder](#9-the-ablation-ladder)
10. [Selectivity indices](#10-selectivity-indices)
11. [What no number here can tell you](#11-what-no-number-here-can-tell-you)
12. [Claims you may make, and claims you may not](#12-claims-you-may-make-and-claims-you-may-not)
13. [A claim card, read section by section](#13-a-claim-card-read-section-by-section)

---

## 1. Engagement is not occupancy

**Where you see it.** The `insect engagement` / `vertebrate engagement` tiles on
the dashboard; the `engagement` column of `flylab occupancy`; the
`engagement_model` field on every receptor row.

**What it means.** FlyLab applies the Hill expression `C^n / (V^n + C^n)` to a
sourced parameter `V`. What that number *is* depends entirely on what the
source measured and how far the source sits from "this compound, at this
receptor, in *Drosophila melanogaster*". The pair (parameter type, evidence
distance) decides the permitted transformation, and the name of that
transformation travels with the number:

| Parameter type | Evidence distance | `engagement_model` | What the number is |
|---|---|---|---|
| Kd / Ki | E0 on-target | `binding_occupancy` | fractional occupancy of the modelled binding site |
| Kd / Ki | E1 other species, E2 related preparation | `binding_engagement_proxy` | a labelled extrapolation from a binding constant |
| Kd / Ki | E3 class statement | `functional_engagement_proxy` | weaker still |
| EC50 / IC50 / Kb | E0 on-target | `functional_engagement` | normalised functional response |
| EC50 / IC50 / Kb | E1–E3 | `functional_engagement_proxy` | a transferred functional response |
| anything | E4 unsupported | `not_modelled` | nothing; no number is produced |

The distances are `E0` exact compound + receptor + species, `E1` right receptor
in another organism, `E2` a hybrid / chimeric / native-mixed preparation, `E3` a
class-order statement, `E4` nothing supports a value. Distance can only *lower*
a claim, never raise it.

Read `flylab/pharm/evidence.py`: `TRANSFORMATION_TABLE` is the whole rule,
`model_for()` applies it, and `check_transformation()` raises
`EvidenceTypeError` rather than silently performing a forbidden transform.

**What it does not mean.** A proxy is not a measurement. A
`binding_engagement_proxy` of 0.97 is not "97% of the *Drosophila* receptors are
occupied"; it is "the Hill expression on a constant measured in another animal
gives 0.97, and FlyLab is telling you so".

**The misreading to avoid.**

> Imidacloprid's insect engagement reads **1.000**, therefore the receptor is
> saturated.

The dashboard prints, on the same line, `binding_engagement_proxy,
E1 (cross-species)`. That 1.000 comes from a Kd of 8.3 × 10⁻¹¹ M measured in
*Myzus persicae* (Bass et al. 2011), on a β1-containing receptor. It is a full
bar on an aphid constant carried across to a fly. The same compound's
*Drosophila*-measured row (`insect_nAChR_native_dmel`, Kd 2 × 10⁻⁹ M, Tomizawa
1996) reads 0.998 — and that one *is* a `binding_occupancy` — while the aggregate
`insect_nAChR` key, which is what actually drives the gain rule, is an EC50
measured on a cat-flea/chicken hybrid receptor and reads 0.991 as a
`functional_engagement_proxy`.

Three numbers, three different quantities, one receptor family. Say which one
you mean.

**One more thing worth knowing.** Of the shipped library's 102 receptor rows,
exactly **one** reaches `binding_occupancy` (imidacloprid at
`insect_nAChR_native_dmel`). 56 are `not_modelled`, 35 are
`functional_engagement_proxy`, 9 are `functional_engagement`, and 1 is
`binding_engagement_proxy`. If you are about to write the word "occupancy",
there is exactly one row in the whole library that entitles you to it.

That imbalance is itself a caveat about the type system, and the referee raised
it (item **I14**): a three-way distinction with one branch taken once is close
to a two-way distinction with extra vocabulary. It is still the right
distinction to record — but do not read the vocabulary as evidence that the
library is full of binding data.

Use the umbrella word **engagement** in prose, and let the detailed output
expose the specific model.

---

## 2. Not modelled is not zero

**Where you see it.** `n/a` and `not modelled` in the CLI tables; a hatched
`NOT MODELLED` chip in the bench; `engagement: null` in a notebook.

**What it means.** A row whose relation is `unsupported` (E4), or whose
parameter type carries no concentration scale (`relative_potency`,
`class_order`, `unknown`), produces **no number at all**. It is excluded from
every numeric aggregation path, from the selectivity indices, from the
engagement maxima that feed the mechanism layer, and from the ablation
orderings. The ablation output says so out loud:

```
$ flylab ablation imidacloprid --conc 1e-6
...
! 1 insect receptor row(s) for imidacloprid carry no modellable value
  (engagement is None); they were skipped, never read as an engagement of 0.
```

This is a deliberate repair of a schema-v2 bug in which missing evidence was
stored as `ec50_M: 0.01` and therefore produced a small, plausible, entirely
invented response.

**What it does not mean.** It does not mean the compound has no effect at that
receptor. It means nobody in the library's sources measured it.

**The misreading to avoid.**

> Imidacloprid's RDL row is empty, so imidacloprid does not act at RDL.

The row says `class-order placeholder: no neonicotinoid action at RDL assumed`.
That is an assumption recorded as an absence, not a measurement of zero. A
compound with two sourced rows and four placeholders has a *narrower* evidence
base than one with six sourced rows, not a cleaner selectivity profile.

Check the census with `flylab dashboard <compound> --conc <M>`, which prints
`(n sourced / m not modelled)` beside the evidence tier.

---

## 3. The three dependence verdicts

**Where you see it.** `flylab dependence`; the "Is the effect wiring-dependent?"
card in the bench; section 6 of a claim card.

**What it means.** `flylab/analysis/dependence.py::mode_verdict()` returns
exactly one of three verdicts per null mode, and it is the only claim this
design supports:

- **`distinguishable`** — the empirical two-sided permutation probability is at
  or below alpha. The real graph and that graph model give different drug
  effects.
- **`equivalent_within_tolerance`** — the test did not reject **and** the gap
  between the real effect and the null ensemble's median is below a
  **prespecified** margin `delta`. This is the only verdict that licenses the
  word *reproduces*.
- **`indeterminate`** — the test did not reject and the gap is not below the
  margin. Consistent with a difference this study cannot resolve.

`delta` is prespecified, not fitted: `equivalence_margin()` sets it to
`DEFAULT_DELTA_FRAC` (5%) of the **vehicle readout of the real graph** — the one
scale in the cell that does not depend on the drug, on the null or on the
sample size. Where the vehicle readout is undefined or zero, no margin exists,
equivalence cannot be claimed, and every non-rejected mode is `indeterminate`.

The headline statistic is the empirical permutation `p = (k+1)/(n+1)`. Its
resolution — and its floor — is `1/(n+1)`. A `p` at the floor means "not
resolvable with this n", not "p = 0". The `z` printed beside it is **not** a
probability: null distributions on degraded graphs are routinely non-normal,
and a large `|z|` can simply mean the shuffles collapsed onto one value.

**What it does not mean.** A large `p` is not evidence of sameness. The three
verdicts exist precisely so that "we could not tell them apart" and "they are
the same within a stated margin" cannot be written with the same words.

**The misreading to avoid.** This one is not hypothetical. Earlier versions of
this project's own manuscript wrote the imidacloprid non-rejections as though
the shuffle had given the same effect. Referee item **B5** records that the
abstract asserted "reproduces" twice, in a place the shipped code contains an
explicit comment refusing to say it.

Here is the real profile, at the paper's `n = 1000`, on the `named` cut — all
three verdicts in one cell:

| mode | p (two-sided) | gap from null median | verdict |
|---|---|---|---|
| `sign_permute` | 0.2747 | 2.652 | **indeterminate** |
| `weight_permute` | 0.5864 | 0.126 | equivalent_within_tolerance |
| `rewire_degree_preserving` | 0.4715 | 0.123 | equivalent_within_tolerance |
| `sign_permute_weight_matched` | 0.3886 | 0.236 | equivalent_within_tolerance |
| `erdos_renyi` | 0.0010 | 6.015 | distinguishable |

(`delta = 0.3316`, five per cent of the vehicle `mean_hz` of 6.631.)

Look at `sign_permute`. Its permutation probability, 0.2747, sits comfortably
between the two modes that *did* reach equivalence — and it is the one
structure-preserving mode that is **not** equivalent, because its gap from the
null median is eight times the margin. Ordering the modes by `p` would have put
it in the middle of the equivalent ones. The permutation probability, on its
own, could not tell you which of the two non-rejection answers you had. That is
the whole point of the three-way verdict.

> **A note about this document and NOVELTY.md.** NOVELTY.md states that
> imidacloprid's sign, weight and degree modes at a thousand permutations are
> indeterminate. As shipped today, only `sign_permute` is; the weight and degree
> modes now fall inside the prespecified margin and report
> `equivalent_within_tolerance`. The *rule* NOVELTY.md states is the shipped
> rule and is unchanged — a non-rejection alone licenses only `indeterminate`.
> When quoting a verdict, take it from a run, not from either document.

Two further things the output tells you and you should not drop:

- **Stability of p.** The result warns when `p` had not stabilised at the chosen
  `n` ("it moved by 0.057 between checkpoints, tolerance 0.020"). The verdict at
  alpha may be unchanged while the value of `p` is not yet quotable.
- **The ladder is not always monotone.** `necessary_information_level()` emits a
  warning when a richer null model is distinguishable while a poorer one is not.
  Referee item **I6** is that 8 of 20 positive landscape cells did this and the
  manuscript never mentioned it. If the warning is there, quote it.

---

## 4. The verdict belongs to the substrate

**Where you see it.** The `--graph` flag; the `graph` field of every dependence
result; the Graph selector in the bench sidebar.

**What it means.** FlyLab ships two committed MaleCNS-derived cuts and they are
**not interchangeable**. `flylab/analysis/dependence.py::cut_census()` reports
why:

| | `named` | `taste_motor` |
|---|---|---|
| nodes | 1126 | 1841 |
| edges | 1360 | 19066 |
| mean degree | 1.21 | 10.36 |
| seed cells | 4 | 76 |
| share of edges landing on seeds | **84.9%** | 9.6% |
| nodes with any input at all | **184 (16.3%)** | 1643 (89.2%) |
| `in_star` | **true** | false |

The `named` cut is an in-star. Eighty-five per cent of its edges terminate on
four cells, and five sixths of its nodes have no input whatsoever. A
degree-preserving rewire of such a graph is close to the identity, so "not
distinguishable from the rewiring null" there is an almost structural result.

And the answer changes. Same compound, same dose, same readout:

```
$ flylab dependence imidacloprid --conc 1e-6 --n 200 --graph named
imidacloprid @ 1.00e-06 M  -  specific wiring evidence: weak (composition-dominated)

$ flylab dependence imidacloprid --conc 1e-6 --n 200 --graph taste_motor
imidacloprid @ 1.00e-06 M  -  specific wiring evidence: present (topology-dependent)
```

On `taste_motor`, `rewire_degree_preserving` reaches the resolution floor
(p = 0.0050 at n = 200) and the cell is topology-dependent.

**What it does not mean.** Neither verdict is "the truth about imidacloprid".
Both are true statements about a readout on a substrate.

**The misreading to avoid.** Referee item **B3**: the paper generalised from the
`named` cut anyway. Never write "imidacloprid's effect does not need the
connectome". Write "on the `named` MaleCNS cut, at 1 µM, on `mean_hz`, the real
graph effect was not distinguishable from the degree-preserving null ensemble at
n = 1000" — and if you are making a negative claim, repeat it on `taste_motor`
first.

NOVELTY.md now treats this reversal as the result rather than a caveat ("The
verdict is substrate-dependent, and that is now the result"), and withdraws
"most predictions do not need the connectome" along with it. The instrument is
not what is in doubt: `synthetic_cut`, `ladder_recovery` and `ladder_power` in
`dependence.py` plant a recurrent loop of known strength in a synthetic graph of
matched size, density and composition, check that the ladder recovers it and
that the unplanted control is not called, and map detection rate against effect
size and permutation count. The reversal is about the graphs.

**Always name four things with a dependence verdict**: the compound, the
concentration, the readout, and the cut.

---

## 5. Confirmatory versus exploratory

**Where you see it.** The `design` and `confirmatory` fields on every dependence
result; the `exploratory` warning in the CLI output; `p_adjusted`, `q_value` and
`fdr_alpha` in a landscape.

**What it means.** `dependence_profile()` marks a run **confirmatory** only when
it is one of `CONFIRMATORY_COMPOUNDS` (imidacloprid, fipronil) at `n >= PAPER_N`
(1000). That is a single prespecified test and needs no multiplicity
correction. Everything else is **exploratory** and says so:

```
! exploratory: this profile is not one of the prespecified confirmatory tests
  (imidacloprid, fipronil at n >= 1000) and carries no multiplicity correction
  of its own. Counting verdicts over many such profiles needs the
  FDR-controlled dependence_landscape instead.
```

`dependence_landscape()` runs 21 compounds × 4 concentrations = **84 cells**,
each against the null modes, applies Benjamini–Hochberg across the run's
*structural* tests (`weight_permute` and `rewire_degree_preserving`, so 168 of
them), classifies on the adjusted probabilities and keeps the raw ones beside
them. It reports both counts, so the change in the headline number is visible
rather than silent.

**What it does not mean.** Failing correction inside the landscape does not
retract the prespecified test. They are different objects answering different
questions:

- The **confirmatory profile** asks: *does this one prespecified effect differ
  from this null ensemble?* One test, no correction owed.
- The **landscape** asks: *how many of 84 screened cells are topology-dependent?*
  168 structural tests; an uncorrected count is inflated, so BH is applied and
  the count drops.

**The real fipronil example.** Run both and compare the *same cell*.

As a prespecified profile, fipronil at 1 µM on the `named` cut is
distinguishable from `weight_permute` (p = 0.0100) and
`rewire_degree_preserving` (p = 0.0149) at n = 200, and its class is
`topology-dependent`. It is one of the two `CONFIRMATORY_COMPOUNDS`; at
n = 1000 that profile is a single prespecified test and owes no correction.

Inside the 84-cell landscape at n = 100, the same cell reads:

| mode | raw p | q after BH | |
|---|---|---|---|
| `weight_permute` | 0.0099 (at the floor `1/101`) | 0.0792 | not rejected |
| `rewire_degree_preserving` | 0.0297 | 0.2079 | not rejected |

`class_raw: topology-dependent` → `class: composition-dominated`. Across the
whole landscape the count moves from **14 of 84** raw to **0 of 84** corrected,
over 168 structural tests.

Both statements are correct, and neither cancels the other. What you may not do
is take the landscape's corrected count and report it as though the prespecified
test had failed, or take the prespecified result and add it to a landscape
count.

There is a third trap, and the landscape guards against it explicitly. When the
permutation resolution `1/(n+1)` is too coarse for the adjusted threshold
`fdr_alpha / m` to be reachable by any single test, no cell *can* be rejected.
That is exactly what happened above, and the run says so rather than silently
reporting zero:

```
! the FDR-adjusted threshold is resolution-limited: with 168 structural tests at
  resolution 0.00990, at least 34 tests must sit at the floor together before any
  is rejected, so an isolated strong cell cannot survive the correction at this n.
! EXPLORATORY: this landscape's permutation resolution (0.00990) does not reach
  the adjusted threshold fdr_alpha/m = 0.000298, so its counts are a screen to be
  confirmed, not a confirmatory result.
```

So the "0 of 84" above is a **resolution** statement, not a biological one. The
same failure mode appears in the fast robustness family — see the next section.

Any count taken from a landscape must be quoted with its `design` label.

---

## 6. Specification stability

**Where you see it.** `flylab stability`; section 7 of a claim card;
`flylab/analysis/robustness.py`.

**What it means.** The engagement-to-gain rules in `flylab/pharm/mechanisms.py`
(`1 + 0.4t - 1.6t²`, `1 + 1.5t`, `1 + 2t`) are **asserted**, never fitted. The
robustness layer re-derives each of seven qualitative conclusions under a
prespecified family of alternative rules and reports the fraction of the family
that retains it, here with the 9-member fast family at 25 shuffles:

```
$ flylab stability --fast --conc 1e-6
  C4_fipronil_topology_exceeds       retained 0/9 FRAGILE
  C1_nicotinic_suppression           retained 3/9 FRAGILE
  C5_nicotinic_buffering             retained 7/9 FRAGILE
  C6_nav_ache_amplification          retained 8/9 FRAGILE
  C2_rdl_disinhibition               retained 9/9 stable
  C3_imidacloprid_topology_not_distinguishable retained 9/9 stable
  C7_map_bitter_veto                 retained 9/9 stable
```

**The nicotinic suppression reversal.** `C1` is "a nicotinic agonist
(imidacloprid, 1 µM) suppresses network activity" — the bench's most visible
single result, a 93% drop in `mean_hz`. It is retained by **3 of 9** members of
the fast family. Under monotone alternatives to the default biphasic rule, the
sign flips.

So the **sign** of the bench's most visible result is a property of the
engagement-to-gain rule at least as much as of the pharmacology or the
connectome. That is one of the most useful things this project has found out
about itself, and it has to travel with the result: quoting the 93% suppression
without it turns a rule-dependent conclusion into a prediction. Referee item
**I15** is exactly that failure — the paper disowned this result in §3.4 and
then relied on it, unflagged, in the prospective predictions.

**What it does not mean.** A stability fraction is **not** a posterior
probability over models. It is the fraction of one finite, prespecified family
— which shares a connectome cut, a library and a runtime with every other
member. Those assumptions are held fixed here and attacked elsewhere.

**The misreading to avoid — and one this page has to soften.** `C4` reads
`retained 0/9`, and it would be very easy to write "the fipronil topology result
does not survive the gain-rule family". That is wrong, and the output says so:

```
! C4_fipronil_topology_exceeds: multiplicity changes the fraction - 9/9 retained
  on the raw permutation p, 0/9 after Benjamini-Hochberg across the run.
! the topology conclusions were evaluated at n_shuffles = 25 per specification
  (permutation resolution 1/(n+1) = 0.0385).
```

At 25 shuffles the finest attainable p is 0.0385; after BH across 36 structural
tests the adjusted threshold is far below that, so **no** specification could
have retained `C4`. The 0/9 is a resolution artefact of the `--fast` family, not
a finding about the gain rules. Use `--full` (25 specifications, 100 shuffles)
before quoting a topology row of this table, and quote both fractions.

Finally, note that the *same shuffle stream* (seed 0) is used for every
specification's topology cell, so the structural tests of one run are positively
dependent rather than independent. BH controls the FDR under positive
dependence, but the corrected fractions are a property of that shuffle stream.

---

## 7. The uncertainty budget and value of information

**Where you see it.** `flylab uncertainty`, `flylab voi`, the Uncertainty panel
in the bench; `flylab/analysis/uncertainty_global.py` and
`flylab/analysis/voi.py`.

**What it means.** A Sobol' decomposition of the variance of a **simulated
readout** across nine declared factors, under **assumed input ranges**. Each
factor is then mapped to a concrete experiment, and `VOI_j = max(0, S_j) · Var(Y)`
is the variance that would disappear if that factor were resolved exactly.

Three labels matter, and the output prints all three:

- **resolved** — the factor's first-order bootstrap confidence interval excludes
  zero. Only these get a number you may quote.
- **unresolved at this sample size** — the CI includes zero. The sample cannot
  separate the effect from nothing. This is a statement about the sample, not
  about the factor, and it carries no ranking claim. Note that the CLI still
  *prints* a number for such a factor: `voi_fraction` in `voi.py` is
  `max(0, S_j)` for every row, and it is the `status` field — not the number —
  that tells you whether the number means anything. Read the status line before
  the table.
- **null control** — `lif_seed` has no effect at all on the deterministic rate
  engine. Its index is therefore a *direct read of the estimator's own error* on
  a factor that is exactly zero.

Here is the real ranking:

```
$ flylab voi imidacloprid --conc 1e-6 --n-base 64
  1. gain_transform       variance removed   0.516  calibration of receptor engagement against synaptic gain
  2. weight_threshold     variance removed   0.474  no experiment: a reconstruction-confidence analysis
  3. drive                variance removed   0.238  in vivo baseline firing rates of the driven cell classes
  4. transmitter          variance removed   0.087  immunostaining or ground-truth transmitter labels
  5. hill_n               variance removed   0.078  the same dose-response, read for slope rather than midpoint
  6. lif_seed             variance removed   0.075  none: increase the number of simulated replicates
  7. expression           variance removed   0.075  FISH / scRNA-seq for receptor expression in the named cell class
  8. potency              variance removed   0.069  electrophysiological dose-response on the stated receptor
  9. gain_coef            variance removed   0.055  the same calibration, read for amplitude rather than shape
```

Read the ranking against row 6. `lif_seed`, the factor that is *exactly zero by
construction*, scores 0.075 — and `potency`, the published parameter that a
newcomer assumes is the most important thing in a pharmacology model, scores
0.069, **below the null control**. Only `gain_transform` has a confidence
interval excluding zero. The asserted engagement-to-gain rule is the single
resolvable source of variance in this model; the cited potency is at the noise
floor.

That is the same fact the tornado in [Lesson 3](TUTORIAL.md#lesson-3) shows
directly: at a saturating dose a two-fold error in the EC50 moves `mean_hz` by
0.000 Hz.

This run also prints `! not converged: doubling n_base moved a first-order index
by 0.23`, so the *ordering* below `gain_transform` is not yet stable. The shape
of the finding — one resolvable model assumption, everything else at the
estimator's noise floor — survives the convergence warning; the individual
shares do not.

**What it does not mean.**

- It is **not** biological variability. Nothing here describes variation between
  flies, between cells or across time.
- It is **not** the expected value of an experiment in the world. It is variance
  of a model output under ranges declared in
  `flylab.analysis.uncertainty_global.FACTORS`. Widening a factor's declared
  range raises its index almost mechanically.
- A real measurement narrows a factor's range; it does not collapse it to a
  point. Every VOI is an upper bound on what its experiment would actually buy.

**The misreadings to avoid.**

> Factor X has a negative index, so resolving it would make things worse.

No. Jansen first-order estimates are kept unclipped on purpose, and at finite
`N` a truly-zero index is estimated as a small signed number. Referee item **I5**
records that v0.6 of `voi.py` passed those signs straight through into the table
(expression −0.018, potency −0.019, gain_coef −0.144 Hz²) as though they were
negative information. `voi.py` now separates `voi_fraction_raw` (unclipped,
a convergence diagnostic) from the decision VOI (never negative).

> Factor Y beat the null control, so its effect is real.

Also no, and the output says so in as many words: *"an index is NOT resolved
merely by exceeding them"*. The null control is **one draw** of the estimator's
error, not a symmetric tolerance band. Resolution is decided per factor by its
own confidence interval.

> These shares are the answer.

Check for `! not converged` first. `flylab uncertainty imidacloprid --conc 1e-6`
at the browser's default `n_base = 32` ranks `weight_threshold` first with a
share of 0.734, puts `potency` at 0.481 — and the null-control `lif_seed`, whose
true index is zero, at **0.480**. Seven of the nine factors have confidence
intervals excluding zero at that sample size, which is the estimator being
confident about its own noise. At `n_base = 128` the same call ranks
`gain_transform` first at 0.362, puts `potency` at 0.054 and `lif_seed` at
0.052, and resolves only `gain_transform`. Both runs print `not converged`
(a first-order index moved by 0.33 and 0.19 respectively when `n_base` was
doubled). The paper uses 128 or more; raise `n_base` before quoting a share, and
read the null control every time.

---

## 8. What the engine's row normalisation does to composition versus topology

**Where you see it.** Nowhere on the dashboard — which was referee item **B1**.
It is in `flylab/circuit/rate.py`: `NORMALISATIONS`, `DEFAULT_NORMALISATION`,
and the `NORMALISATION_NOTES` table.

**What it means.** The rate engine divides each row of its weight matrix by that
row's absolute weight sum (`normalise="row_abs"`, the shipped default and the
one behind every published FlyLab number). Gains are then applied as a
presynaptic column scale. The consequence, in the module's own words:

> each cell's recurrent input is a composition-weighted average of its
> presynaptic gains

That is a modelling choice that bounds rates and keeps the iteration stable. It
is not a measurement of the circuit. And it **partly builds
"composition-dominated" into the readout that the dependence analysis then
tests**.

**What it does not mean.** It does not make the composition-dominated verdicts
wrong. It makes them conditional.

**The misreading to avoid.** Referee item **I3**, and a claim NOVELTY.md has
formally withdrawn: a bare rank correlation of about 0.99 between the
composition-only ablation level and the full model, quoted as a finding about
the connectome. Two reasons it is not:

1. **It has a structural floor.** `composition_reference_distribution()` draws
   pharmacology-free pseudo-compounds as arbitrary gain vectors — compounds that
   do not exist and that no pharmacology connects — and NOVELTY.md records that
   they already reach a median rank correlation around 0.8. The observed value
   is that floor plus a small excess. The same function also shows that
   permuting the gain vectors across compound labels leaves the correlation
   *provably unchanged*, which is the sharpest form of the point: the number
   carries no compound-level information at all.
2. **It does not survive a different normalisation.**
   `composition_dominance_under_normalisations()` re-runs the same comparison
   under each of the three modes:

   | normalisation | Spearman ρ (composition vs full) | Pearson r | reproduces ordering |
   |---|---|---|---|
   | `row_abs` (default) | 0.988 | 0.756 | yes |
   | `none` | 0.906 | 0.935 | yes |
   | `degree` | **0.666** | 0.827 | **no** |

   Under a degree-corrected engine the agreement falls away, and at lower
   concentrations it inverts.

**The rule.** Never quote a composition-versus-topology comparison without both
its matched reference distribution and the normalisation it was computed under.

---

## 9. The ablation ladder

**Where you see it.** `flylab ablation`; `flylab/analysis/baselines.py`.

**What it means.** Four levels — `A_receptor_only`, `B_composition_only`,
`C_topology_only`, `D_full_flylab` — are scored across the library, and each is
correlated against the full model's ordering. The question is which layer of
information carries the ordering.

**What it does not mean.** The four levels **do not share units**. Only their
*ordering* of the compound library is comparable, which is what the correlations
measure. The output says so.

**Two misreadings to avoid, both of them the project's own.**

1. *"The real cut with a generic, mechanism-free multiplier is the worst level
   of the four."* Withdrawn (referee **I1**, **I2**). That reading rested on a
   topology-only rule that applied a **depression-only** multiplier to every
   transmitter alike. A baseline that can only reduce activity cannot order a
   library containing compounds whose full-model effect is positive; its poor
   score was a property of the rule. Under the direction-aware rule now shipped
   as the default (`GENERIC_RULES` in `baselines.py`), which keeps a generic
   magnitude and takes one sign bit from the mechanism table, that level's rank
   correlation with the full model rises from about 0.24 to about 0.96. The old
   rule survives as a reported floor, never as a competitor.
2. *"Level B reproduces the full model."* See section 8. `reproduces_full` in
   this output is computed on `|ρ|`, so a **perfectly inverted** ordering would
   count as reproduction (referee **I4**). Read the sign.

Also note the warning that levels B and D use opposite signs for glutamate (one
uses the `wholens` convention, the rate engine signs glutamate −0.4);
`detail.effect_other_convention` gives the same compound under the other
convention, and `glutamate_sign_reconciliation()` exists for exactly this.

---

## 10. Selectivity indices

**Where you see it.** The `Receptor selectivity` tile; `recSI` and `cirSI` in
`flylab compare`.

**What it means.** `receptor_si_log10` is a ratio of potencies between an insect
and a vertebrate receptor at the same free concentration. `circuit_si_log10` is
a different quantity: how far this simulation's circuit moves before the
vertebrate receptor in the teaching library fills.

They rank compounds differently, and the bench will tell you when:

```
$ flylab compare imidacloprid deltamethrin --conc 1e-6 --no-occupancy
compound         target                  insect    vert   recSI   cirSI     gap  circuit d%  topology
Imidacloprid     insect_nAChR_beta1       1.000   0.091    2.70    1.79   -0.91       -93.1  composition-dominated
Deltamethrin     insect_Nav               0.959   0.167    2.07    2.28    0.22       560.4  composition-dominated
  * The highest receptor selectivity (Imidacloprid, 2.70 log10) is NOT the highest circuit selectivity (Deltamethrin, 2.28 log10).
  ! Receptor selectivity and circuit selectivity are different quantities; a large receptor ratio does not buy a wide circuit window.
```

**What it does not mean.** Neither index is a safety margin, a therapeutic
index, or a statement about a vertebrate. The vertebrate side is a teaching
scorecard that is *scored only* and never patches a circuit
(`docs/ARCHITECTURE.md`). Eight library compounds have no circuit selectivity
index at all, because RDL/GluCl block cannot reach the effect threshold on these
cuts — that is a property of the cut, not evidence of safety.

**The misreading to avoid.** "500× selective for the insect receptor" is a
ratio of two literature-order teaching numbers, one of which
(`vertebrate_nAChR_a4b2` for imidacloprid) is an `E3` class extrapolation. Never
write "safe" from this scorecard.

---

## 11. What no number here can tell you

FlyLab has never dosed a fly. Its own warning stream says so on every run. In
particular:

- **Nothing about a living animal.** There is no independent, out-of-sample
  biological validation of any circuit-level claim. Comparisons against
  published compound rankings are **literature concordance**, and some of them
  share a source with the library, so they are not even out-of-sample. Agreement
  between the rate and LIF engines is implementation consistency; the two agree
  on *direction* only, and absolute rate agreement should not be assumed.
- **No dose.** The concentration you type is a **free concentration at the
  receptor**, applied directly. No model maps an administered dose to it.
  Cuticular penetration, haemolymph binding, the blood–brain barrier and
  metabolism are all outside the model. The exposure panel's PK constants are
  labelled `placeholder` on the dashboard's own trust table.
- **No exposure.** `flylab exposure` is a one-compartment teaching model, listed
  in NOVELTY.md under "supporting capabilities, not primary novelty".
- **No biological variance.** Between flies, between cells, across time: not
  represented. The rate engine is deterministic; running it with a different
  seed changes nothing. The ensemble panel's spread is Monte-Carlo over declared
  library uncertainty, not measured variation.
- **No whole fly.** The cuts are hops-limited neighbourhoods of MaleCNS v1.0
  with a synapse-count floor: 1126 or 1841 cells, against a reconstruction with
  tens of millions of edges.
- **No measured transmitter identities.** Every edge sign comes from MaleCNS
  *predicted* consensus transmitters. A wrong prediction flips a synapse, and
  the drug patch follows the label rather than the biology.
- **No per-cell receptor expression.** The gain is applied uniformly to every
  cell of a transmitter class; about 20% of cells are mapped. A cell that does
  not express the receptor is patched exactly like one that does.

And there is deliberately **no single confidence score** anywhere in the bench.
The layers fail independently, and one blended number would hide which of them
is weak.

---

## 12. Claims you may make, and claims you may not

This table is the operational form of [NOVELTY.md](NOVELTY.md). Where the two
differ, NOVELTY.md wins.

| You may say | You may not say |
|---|---|
| "not distinguishable from the degree-preserving null ensemble at n = 1000" | "the degree-preserving graph reproduces the effect" |
| "equivalent within the prespecified margin (δ = 5% of vehicle) for `weight_permute`" | "equivalent" for any non-rejection that did not meet the margin |
| "on the `named` cut, at 1 µM, on `mean_hz`" | "the connectome matters for imidacloprid" |
| "functional engagement, from an EC50 transferred from a hybrid receptor" | "97% receptor occupancy" |
| "binding-derived occupancy" only when parameter type is Kd/Ki **and** distance is E0 | "occupancy" for any proxy row |
| "no sourced value; the row is not modelled" | "no effect at this receptor" |
| "literature concordance with published orderings" | "validated" |
| "model-derived" / "prospective software prediction" | "pre-registered", while the prediction lives in a mutable repository |
| "this conclusion is retained by 3 of 9 prespecified gain rules" | "our model predicts network suppression" |
| "the model's variance under declared input ranges is dominated by the gain rule" | "the gain rule explains 52% of biological variability" |
| "topology-sensitive classifications concentrate among perturbations affecting inhibitory signalling, with exceptions" | "exactly the chloride-channel blockers are topology-dependent" |
| "composition-dominated under `row_abs` normalisation, against a pseudo-compound reference median of ≈0.8" | a bare ρ ≈ 0.99 between composition-only and the full model |
| "the vertebrate receptor scorecard is scored, never simulated" | "safe for vertebrates" |
| "a simulation on a public connectome" | "a digital twin of a fly" |

---

## 13. A claim card, read section by section

A claim card is the reader-facing projection of the nine-link audit chain in
`flylab/analysis/claims.py`. It has exactly nine sections, always in the same
order, always all present — so a claim can never quietly omit the section that
would have undermined it (`flylab/report/card.py`, `SECTIONS`).

This is the card written by the spec run in
[Lesson 6](TUTORIAL.md#lesson-6), read one section at a time.

**1. Claim.**

> On the MaleCNS v1.0 `named` cut (1126 cells, 1360 edges, >= 5 synapses),
> FlyLab's rate engine predicts that fipronil at 1.00e-06 M raises mean_hz from
> 6.631 under vehicle to 7.531 (+13.6%).

*Read:* every qualifier you need in order to quote it is in the sentence — cut,
size, synapse floor, engine, compound, dose, readout, vehicle, direction. If you
find yourself shortening it, you are dropping one of them.

**2. Status.** "This is a prediction of a simulation, not an observation... the
notebook's `live_lab` field is null."

*Read:* `live_lab` is the only slot in a FlyLab notebook that can hold a measured
animal number, and nothing in FlyLab ever writes to it. A human has to import a
CSV. `null` means no measurement has ever been attached to this claim.

**3. Evidence inputs.** Four rows, two of them `not modelled`:

| receptor | type | value (M) | evidence distance | engagement |
|---|---|---|---|---|
| `insect_nAChR` | unknown | not modelled | unsupported | N/A |
| `vertebrate_nAChR_a4b2` | unknown | not modelled | unsupported | N/A |
| `vertebrate_GABA_A` | IC50 | 1.10e-06 | `exact_compound_exact_receptor_other_species` | 0.476 |
| `insect_RDL` | IC50 | 3.00e-08 | `exact_compound_exact_receptor_exact_species` | 0.985 |

*Read:* two rows are sourced and only **one** of them can move a circuit.
`insect_RDL` is an IC50 measured in *Drosophila* (Hosie et al. 1995) at E0 — so
it is `functional_engagement`, a normalised functional response, and **not** an
occupancy even though it is on-target and reads 0.985. `vertebrate_GABA_A` is
the same parameter type at E1 (human recombinant receptors), so it is a proxy,
and it is scored rather than simulated: the vertebrate side never patches a
circuit. The two N/A rows are absences, not zeros. See sections
[1](#1-engagement-is-not-occupancy) and [2](#2-not-modelled-is-not-zero).

**4. Measured, not assumed.** "Exactly one link of the nine-link chain behind
this number is a measurement of the system being simulated: the synaptic
wiring."

*Read:* this is the section that sets the ceiling on the whole card. One
observation; everything downstream is assertion or computation.

**5. Assumptions introduced.** Transmitter sign, gain mapping, free
concentration and uniform expression, followed by ten more drawn from the audit
chain (the Hill relation itself, the discarded assay conditions behind an
EC50/IC50, the mechanism rule and its uniform application, the weight floor on
the cut, the predicted transmitter labels, the hand-set drive) — each with a
sentence on how it could change the sign or the size of the claim.

*Read:* `gain_mapping` is flagged here as "the single assumption worth an
experiment", which is the VOI ranking in section
[7](#7-the-uncertainty-budget-and-value-of-information) arriving in the card.
`free_concentration` is where the absence of pharmacokinetics is recorded.

**6. Connectome dependence.**

> Distinguishable from `weight_permute` (p = 0.0100), `rewire_degree_preserving`
> (p = 0.0149), `erdos_renyi` (p = 0.0050) at n = 200. Not distinguishable from
> `sign_permute` (p = 0.1194), `sign_permute_weight_matched` (p = 0.1841).
> Dependence class: **topology-dependent**. Necessary information level:
> `wiring_with_weighted_transmitter_balance`.

*Read:* "not distinguishable" — not "the same". This run is at n = 200, so it is
**exploratory**; the confirmatory version is the same compound at n = 1000.
`erdos_renyi` at p = 0.0050 is sitting exactly on the resolution floor
`1/(n+1)`, so report it as `p ≤ floor`, not as an exact value.

If the spec did not request a dependence analysis, this section reads *"not
assessed in this run"* — which is an unanswered question, not an absent concern.

**7. Specification stability.** In this run: *"not assessed... add
`robustness: {specifications: default_family}` to the spec's `analyses:`
block."*

*Read:* unanswered. The card tells you the exact spec line that would answer it.
See section [6](#6-specification-stability) for what the answer looks like and
how to read the topology rows.

**8. Independent biological validation.** "None." Plus what would change it: "an
electrophysiological or behavioural measurement in adult *Drosophila*, obtained
after this prediction was recorded, at a dose and readout the model actually
names."

*Read:* this section is always "None" for a circuit-level claim in the current
release. It is present anyway, because a dropped section reads as a satisfied
one.

**9. Named unknowns.** CNS free concentration, MN9's receptor complement,
biological variance, and one entry per receptor with no sourced value — each
with the measurement that would resolve it.

*Read:* these are named rather than parameterised on purpose. FlyLab does not
supply a plausible number in their place.

And then the footer, which is the honest summary of all nine:

> FlyLab is a simulation on a public connectome. Nothing on this card is a
> measured drug effect in a living fly.

---

## Related documentation

- [TUTORIAL.md](TUTORIAL.md) — a worked walkthrough that produces every output on this page
- [NOVELTY.md](NOVELTY.md) — the claim guardrail; authoritative where this page and it disagree
- [ARCHITECTURE.md](ARCHITECTURE.md) — the dataflow and the layer boundaries
- [WORKFLOW.md](WORKFLOW.md) — specs, claim cards and the artifact manifest
- [README.md](README.md) — the documentation index and its authority map
