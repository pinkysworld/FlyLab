"""Mechanism -> circuit gain patches (single source of truth).

Every assay (`taste`, `wholens`, `subgraph`, spiking) must obtain its gain
dictionary from :func:`gains_from_occupancy` so that one receptor rule change
propagates everywhere. The four rules shipped in v0.4 are reproduced verbatim
(`HANDOFF.md` forbids changing them silently).

Gains are dimensionless multipliers on the synaptic weight of a transmitter
class (or, for ``g_nav``, on node excitability). 1.0 means "vehicle".
"""

from __future__ import annotations

from typing import Any, Iterable, Mapping

__all__ = [
    "GAIN_KEYS",
    "MECHANISM_TABLE",
    "default_gains",
    "gains_from_occupancy",
    "mechanism_table_rows",
    "describe_gain",
]

GAIN_KEYS = ("g_ach", "g_gaba", "g_glu", "g_oct", "g_nav", "ach_tone")

#: Which transmitter each gain multiplies, for the UI legend.
GAIN_TARGET = {
    "g_ach": "acetylcholine synapses",
    "g_gaba": "GABA synapses",
    "g_glu": "glutamate synapses",
    "g_oct": "octopamine synapses",
    "g_nav": "global node excitability",
    "ach_tone": "ambient acetylcholine level (not a synaptic weight)",
}

_FLOOR = 0.05  # gains never reach zero: a dead network is not informative.

# Direction synonyms. Insect receptors are excited by agonists and allosteric
# activators alike, so both map onto the same rule family.
_ACTIVATING = {"agonist", "partial_agonist", "positive_modulator"}
_BLOCKING = {"antagonist", "negative_modulator", "inhibitor"}


def _agonist_curve(th: float) -> float:
    """v0.4 nAChR agonist patch: gain rises, then desensitises/blocks."""
    return 1.0 + 0.4 * th - 1.6 * th * th


def default_gains() -> dict[str, float]:
    """All gains at vehicle (1.0)."""
    return {k: 1.0 for k in GAIN_KEYS}


def gains_from_occupancy(rows: Iterable[Mapping[str, Any]]) -> dict[str, float]:
    """Map occupancy rows (from :func:`flylab.pharm.occupancy.compare_compound`)
    onto circuit gains.

    Rules, one line of rationale each:

    * ``insect_nAChR`` activating -> ``g_ach = max(0.05, 1 + 0.4*th - 1.6*th^2)``
      -- low agonist occupancy excites cholinergic synapses, high occupancy
      desensitises / depolarisation-blocks them (v0.4 rule, verbatim).
    * ``insect_nAChR`` blocking -> ``g_ach = max(0.05, 1 - th)``
      -- a competitive blocker removes that fraction of cholinergic drive
      (v0.4 rule, verbatim).
    * ``insect_RDL`` blocking -> ``g_gaba = max(0.05, 1 - th)``
      -- blocking the GABA-gated chloride channel removes inhibition
      (v0.4 rule, verbatim).
    * ``insect_RDL`` activating -> ``g_gaba = max(0.05, 1 + 0.4*th)``
      -- direct or allosteric opening of RDL adds inhibitory chloride current
      (v0.4 rule from ``assays/subgraph.py``, verbatim).
    * ``insect_GluCl`` activating -> ``g_glu = max(0.05, 1 + 0.8*th)``
      -- avermectins open an inhibitory glutamate-gated chloride channel, so
      glutamatergic (inhibitory-signed) edges get stronger; the coefficient is
      larger than RDL's because GluCl opening is effectively irreversible.
    * ``insect_GluCl`` blocking -> ``g_glu = max(0.05, 1 - th)``
      -- mirror of the activating rule, for completeness.
    * ``insect_AChE`` inhibitor -> ``ach_tone = 1 + 2.0*th`` and
      ``g_ach *= ach_tone * (1 + 0.4*th - 1.6*th^2)``
      -- unhydrolysed acetylcholine accumulates (tone up), and that same excess
      first excites then desensitises the nicotinic receptors, which is the
      classic organophosphate excitation-then-block picture.
    * ``insect_Nav`` activating -> ``g_nav = 1 + 1.5*th``
      -- pyrethroids and DDT hold sodium channels open, raising excitability of
      every node regardless of transmitter.
    * ``insect_Nav`` blocking -> ``g_nav = max(0.05, 1 - th)``
      -- a channel blocker (TTX-like) lowers excitability instead.
    * ``insect_OctR`` activating -> ``g_oct = 1 + 0.5*th``
      -- formamidines mimic octopamine, a modulatory amine, so the coefficient
      is deliberately smaller than the fast-transmitter ones.
    * ``insect_OctR`` blocking -> ``g_oct = max(0.05, 1 - th)``.

    Vertebrate receptor rows never change a fly gain; they exist for the
    selectivity scorecard only. Unknown receptors and ``direction: none`` rows
    are ignored.

    Args:
        rows: occupancy rows with ``receptor``, ``occupancy`` and ``direction``.

    Returns:
        dict with keys ``g_ach``, ``g_gaba``, ``g_glu``, ``g_oct``, ``g_nav``,
        ``ach_tone``; all default to 1.0.
    """
    g = default_gains()
    if not rows:
        return g

    ache_th = 0.0
    for row in rows:
        receptor = str(row.get("receptor", ""))
        direction = str(row.get("direction") or "none")
        th = float(row.get("occupancy") or 0.0)
        th = min(max(th, 0.0), 1.0)

        if receptor == "insect_nAChR":
            if direction in _ACTIVATING:
                g["g_ach"] = max(_FLOOR, _agonist_curve(th))
            elif direction in _BLOCKING:
                g["g_ach"] = max(_FLOOR, 1.0 - th)
        elif receptor == "insect_RDL":
            if direction in _BLOCKING:
                g["g_gaba"] = max(_FLOOR, 1.0 - th)
            elif direction in _ACTIVATING:
                g["g_gaba"] = max(_FLOOR, 1.0 + 0.4 * th)
        elif receptor == "insect_GluCl":
            if direction in _ACTIVATING:
                g["g_glu"] = max(_FLOOR, 1.0 + 0.8 * th)
            elif direction in _BLOCKING:
                g["g_glu"] = max(_FLOOR, 1.0 - th)
        elif receptor == "insect_Nav":
            if direction in _ACTIVATING:
                g["g_nav"] = 1.0 + 1.5 * th
            elif direction in _BLOCKING:
                g["g_nav"] = max(_FLOOR, 1.0 - th)
        elif receptor == "insect_OctR":
            if direction in _ACTIVATING:
                g["g_oct"] = 1.0 + 0.5 * th
            elif direction in _BLOCKING:
                g["g_oct"] = max(_FLOOR, 1.0 - th)
        elif receptor == "insect_AChE":
            if direction in _BLOCKING:
                ache_th = max(ache_th, th)

    if ache_th > 0.0:
        g["ach_tone"] = 1.0 + 2.0 * ache_th
        g["g_ach"] = max(_FLOOR, g["g_ach"] * g["ach_tone"] * _agonist_curve(ache_th))
    return g


def _rule(receptor: str, direction: str, gain: str, formula: str, rationale: str) -> dict[str, str]:
    return {
        "key": f"{receptor}:{direction}:{gain}",
        "receptor": receptor,
        "direction": direction,
        "gain": gain,
        "formula": formula,
        "rationale": rationale,
        "applies_to": GAIN_TARGET[gain],
    }


#: ``"<receptor>:<direction>" -> rule`` description of every patch rule.
#: JSON-serialisable so the server, UI and paper can render the same table.
MECHANISM_TABLE: dict[str, dict[str, str]] = {
    r["key"]: r
    for r in [
        _rule("insect_nAChR", "agonist", "g_ach", "max(0.05, 1 + 0.4*th - 1.6*th^2)",
              "Cholinergic excitation at low occupancy, desensitisation/block at high occupancy (v0.4 rule)."),
        _rule("insect_nAChR", "partial_agonist", "g_ach", "max(0.05, 1 + 0.4*th - 1.6*th^2)",
              "Partial agonists open the same cation channel; same curve as full agonists."),
        _rule("insect_nAChR", "positive_modulator", "g_ach", "max(0.05, 1 + 0.4*th - 1.6*th^2)",
              "Allosteric activators (spinosyns) open the same channel without using the ACh site."),
        _rule("insect_nAChR", "antagonist", "g_ach", "max(0.05, 1 - th)",
              "A blocker removes that fraction of cholinergic drive (v0.4 rule)."),
        _rule("insect_RDL", "antagonist", "g_gaba", "max(0.05, 1 - th)",
              "Blocking the GABA-gated chloride channel removes inhibition (v0.4 rule)."),
        _rule("insect_RDL", "agonist", "g_gaba", "max(0.05, 1 + 0.4*th)",
              "Opening RDL adds inhibitory chloride current (v0.4 rule from subgraph.py)."),
        _rule("insect_RDL", "positive_modulator", "g_gaba", "max(0.05, 1 + 0.4*th)",
              "Allosteric potentiation of RDL is treated like partial agonism."),
        _rule("insect_GluCl", "agonist", "g_glu", "max(0.05, 1 + 0.8*th)",
              "GluCl is inhibitory; avermectin opening is near-irreversible, hence the larger coefficient."),
        _rule("insect_GluCl", "antagonist", "g_glu", "max(0.05, 1 - th)",
              "Mirror rule: blocking GluCl removes glutamatergic chloride current."),
        _rule("insect_AChE", "inhibitor", "ach_tone", "1 + 2.0*th",
              "Unhydrolysed acetylcholine accumulates, raising ambient cholinergic tone."),
        _rule("insect_AChE", "inhibitor", "g_ach", "g_ach * ach_tone * (1 + 0.4*th - 1.6*th^2)",
              "The excess ACh first excites and then desensitises nicotinic receptors: OP excitation then block."),
        _rule("insect_Nav", "positive_modulator", "g_nav", "1 + 1.5*th",
              "Pyrethroids/DDT hold sodium channels open, raising excitability of every node."),
        _rule("insect_Nav", "antagonist", "g_nav", "max(0.05, 1 - th)",
              "Sodium-channel block lowers excitability of every node."),
        _rule("insect_OctR", "agonist", "g_oct", "1 + 0.5*th",
              "Formamidines mimic octopamine; modulatory amine, so a smaller coefficient than fast transmitters."),
        _rule("insect_OctR", "antagonist", "g_oct", "max(0.05, 1 - th)",
              "Octopamine-receptor block removes modulatory drive."),
    ]
}


def mechanism_table_rows() -> list[dict[str, str]]:
    """MECHANISM_TABLE as a list, for tables in the UI and the paper."""
    return list(MECHANISM_TABLE.values())


def describe_gain(key: str) -> str:
    """One-line description of what a gain multiplies."""
    return GAIN_TARGET.get(key, "unknown gain")
