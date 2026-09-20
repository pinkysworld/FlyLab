"""Scale-ladder constants and naming: the part of the extractor with no heavy deps.

:mod:`flylab.maps.extract` needs pandas and pyarrow to *build* a cut.  The
analysis layer only needs to know which rungs exist, what they are called and
by what rule they were made, and :mod:`flylab.analysis.scale` is reachable from
:mod:`flylab.browser.bridge`, which must import with pandas, pyarrow, pydantic,
fastapi and typer all hidden (that is what makes the Pyodide bench possible and
``tests/test_browser_bridge.py`` enforces it).

So the pure data and naming live here, with a stdlib-only import list, and
``extract.py`` re-exports them.  Nothing in this module may import pandas,
pyarrow or numpy.
"""
from __future__ import annotations

DEFAULT_TYPES: tuple[str, ...] = ("MN9", "DNp01")
GUSTATORY_TYPES: tuple[str, ...] = ("LB1a", "LB1b", "LB1c", "LB1d", "LB3b", "LB3c")
TASTE_MOTOR_TYPES: tuple[str, ...] = DEFAULT_TYPES + GUSTATORY_TYPES

#: Node budgets the ladder targets.  The rungs are **nested**: the node set of
#: a smaller rung is a prefix of the next one's, so the ladder is a sequence of
#: views of one object rather than five unrelated graphs.
LADDER_SIZES: tuple[int, ...] = (1000, 5000, 10000, 25000, 50000)

#: Seeds the ladder grows from.  Deliberately the ``taste_motor`` seed set:
#: MN9 and DNp01 keep ``mn9_hz`` / ``dnp01_hz`` defined at every rung, and the
#: labellar GRNs keep the taste drive defined, so the *readout* is the same
#: quantity all the way up and only the surrounding graph changes.
LADDER_SEED_TYPES: tuple[str, ...] = TASTE_MOTOR_TYPES

#: Synapse-count floor held fixed across the whole ladder.  5 is the floor the
#: committed cuts use.
LADDER_MIN_WEIGHT = 5

#: A cut above this many bytes of JSON is a CI artifact, not a git object.
COMMIT_BYTE_BUDGET = 1_400_000

LADDER_RECIPE = """\
How a cut is grown, and why this way
------------------------------------
A rung is the **induced subgraph on the first K cells of one fixed,
deterministic growth order** out of the seed set.  The order is:

1. rank 0: the seed cells themselves (all traced cells whose MaleCNS type
   matches ``seed_types``), in ascending bodyId;
2. rank r+1: every traced cell not yet selected that shares at least one
   edge of weight >= ``min_weight`` with the set selected after rank r,
   sorted by the **total synaptic weight joining it to that set**
   (descending, summed over both directions), ties broken by ascending
   bodyId;
3. stop when K cells have been taken.

The cut is then every edge of weight >= ``min_weight`` whose two ends are
both selected (an induced subgraph, the same closure rule the committed
``taste_motor`` cut uses).

Why not the alternatives:

* **More hops.**  The traced-to-traced MaleCNS graph has 25 563 197 edges
  at weight 1 (mean degree 155) and 6 235 682 at weight 5 (mean degree 38).
  One hop from the 76 seeds already reaches 1 304 cells and two hops reach
  55 127.  Hop count is not a dial that can be set to 5 000 or 25 000, and
  the jump from one hop to two changes the cut by a factor of 42.
* **A lower weight floor.**  Growing by admitting weaker synapses changes
  the *edge* population as well as the node population, so mean degree and
  the topology/composition balance would move for two reasons at once and
  a verdict change could not be attributed to scale.  The floor is
  therefore frozen at ``min_weight`` for every rung.
* **A whole neuropil.**  Well defined, but the resulting cuts are not
  nested, do not contain the same seed cells, and have no size dial.

What the rule does *not* fix: mean degree is not constant along the ladder
(the growth order takes the densely connected core first), and each cut's
census is recorded in its metadata precisely so that a reader can see the
structural change alongside the verdict change.

Determinism: the order is a function of the weight matrix, the seed types
and the floor only.  No RNG, no iteration over Python sets, no
platform-dependent sort -- ``np.lexsort`` on (bodyId asc, -score) is a
stable total order.  The one residual is float summation order inside
``np.bincount``, which can in principle tie two cells whose scores differ
by <1 ulp; scores are integer synapse counts summed in float64 and the
largest is ~1e5, so this cannot happen below 2**53.
"""


def ladder_filename(n_target: int) -> str:
    """Committed/artifact filename for a rung (``malecns_scale_5k.json``)."""
    k = int(n_target)
    label = f"{k // 1000}k" if k >= 1000 and k % 1000 == 0 else str(k)
    return f"malecns_scale_{label}.json"


__all__ = [
    "DEFAULT_TYPES",
    "GUSTATORY_TYPES",
    "TASTE_MOTOR_TYPES",
    "LADDER_SIZES",
    "LADDER_SEED_TYPES",
    "LADDER_MIN_WEIGHT",
    "LADDER_RECIPE",
    "COMMIT_BYTE_BUDGET",
    "ladder_filename",
]
