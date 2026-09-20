"""Retrospective validation against published potency orderings.

Data source: ``data/literature/published_rank_orders.yaml`` (12 entries, each
with its assay, species, readout, ``use_in_flylab`` flag and citation).  Those
entries are **targets**, never inputs: nothing here is fitted to them, and a
mismatch is reported rather than repaired.

For every entry whose compounds are all in the library, the model's ordering
(receptor occupancy, or circuit displacement when a circuit assay is named) is
compared with the published ordering by Spearman's rho and Kendall's tau-b,
both implemented here without scipy, plus the full pairwise concordance table.

Two known discrepancies are surfaced explicitly by :func:`validate_all` and are
**not** to be tuned away (see :data:`KNOWN_DISCREPANCIES`):

1. Nitenpyram sits near the top of the model's receptor-potency ordering and at
   the bottom of the *Drosophila* whole-animal bioassay ordering (Perry et al.
   2012).
2. Clothianidin versus imidacloprid is inverted: the library makes imidacloprid
   marginally the more potent, both the acute bioassay and the chronic-lifespan
   ordering put clothianidin first.

The explanation is not a bug: a receptor EC50 and a whole-animal LC50 are
different quantities.  The whole-animal number contains uptake, distribution,
metabolism (nitenpyram has poor persistence) and excretion; the receptor number
does not.  A model that reproduces receptor pharmacology is not obliged to
reproduce a bioassay ordering, and saying so is part of the result.
"""

from __future__ import annotations

import copy
import math
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import yaml

from flylab.pharm.occupancy import compare_compound, engagement, load_library, spec_value_M

__all__ = [
    "RANK_ORDERS_PATH",
    "KNOWN_DISCREPANCIES",
    "ASSAYS",
    "spearman_rho",
    "kendall_tau",
    "load_rank_orders",
    "validate_entry",
    "validate_all",
    "validation_report_markdown",
]

_DATA_DIRS = [
    Path(__file__).resolve().parents[2] / "data" / "literature",
    Path("data/literature"),
]
RANK_ORDERS_PATH = next(
    (d / "published_rank_orders.yaml" for d in _DATA_DIRS if (d / "published_rank_orders.yaml").exists()),
    _DATA_DIRS[0] / "published_rank_orders.yaml",
)

ASSAYS = ("occupancy", "subgraph")

QUANTITY_NOTE = (
    "A receptor EC50 and a whole-animal LC50/EC50 are DIFFERENT QUANTITIES. "
    "The whole-animal number contains uptake, cuticular penetration, "
    "distribution, metabolism and excretion; the receptor number does not. "
    "FlyLab models receptor occupancy at a stated free concentration, so it is "
    "not obliged to reproduce a bioassay ordering - and where it does not, the "
    "difference is reported, not corrected."
)

#: Vertebrate-preparation markers in an entry's ``species`` string.
_VERTEBRATE_WORDS = ("rattus", "rat ", "human", "homo", "mouse", "mus ", "mammal", "vertebrate", "chicken")
_INSECT_WORDS = ("drosophila", "musca", "myzus", "insect", "nilaparvata", "apis", "plutella", "house fly")


# ---------------------------------------------------------------------------
# rank statistics (no scipy)
# ---------------------------------------------------------------------------
def _ranks(values: Sequence[float]) -> list[float]:
    """Average ranks (1-based), ties shared - the tie handling rho needs."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            out[order[k]] = avg
        i = j + 1
    return out


def spearman_rho(x: Sequence[float], y: Sequence[float]) -> float:
    """Spearman's rank correlation: Pearson's r on average ranks.

    Implemented directly (FlyLab does not depend on scipy).  Returns ``nan``
    when either ranking has no variance (e.g. every published rank is tied).
    """
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    if len(x) < 2:
        return float("nan")
    rx, ry = np.array(_ranks(list(x))), np.array(_ranks(list(y)))
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denom = math.sqrt(float((rx * rx).sum()) * float((ry * ry).sum()))
    if denom <= 0:
        return float("nan")
    return float((rx * ry).sum() / denom)


def kendall_tau(x: Sequence[float], y: Sequence[float]) -> float:
    """Kendall's tau-b (ties corrected), implemented without scipy."""
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    n = len(x)
    if n < 2:
        return float("nan")
    concordant = discordant = tx = ty = 0
    for i in range(n):
        for j in range(i + 1, n):
            dx = np.sign(x[i] - x[j])
            dy = np.sign(y[i] - y[j])
            if dx == 0 and dy == 0:
                tx += 1
                ty += 1
            elif dx == 0:
                tx += 1
            elif dy == 0:
                ty += 1
            elif dx * dy > 0:
                concordant += 1
            else:
                discordant += 1
    n0 = concordant + discordant
    denom = math.sqrt((n0 + tx) * (n0 + ty))
    if denom <= 0:
        return float("nan")
    return float((concordant - discordant) / denom)


# ---------------------------------------------------------------------------
# dataset
# ---------------------------------------------------------------------------
@lru_cache(maxsize=4)
def _read_rank_orders(path_str: str) -> dict[str, Any]:
    p = Path(path_str)
    if not p.exists():  # pragma: no cover
        return {"entries": [], "_missing": str(p)}
    data = yaml.safe_load(p.read_text()) or {}
    data.setdefault("entries", [])
    return data


def load_rank_orders(path: Path | None = None) -> dict[str, Any]:
    """Deep copy of the published-rank-order dataset."""
    return copy.deepcopy(_read_rank_orders(str(path or RANK_ORDERS_PATH)))


def _entry_by_id(entry_id: str, path: Path | None = None) -> dict[str, Any]:
    for entry in load_rank_orders(path).get("entries", []):
        if entry.get("id") == entry_id:
            return entry
    raise KeyError(f"unknown rank-order entry {entry_id!r}")


def _norm(name: Any) -> str | None:
    if not isinstance(name, str):
        return None
    key = name.strip().lower()
    if not key or "(" in key or " " in key:
        return None
    return key.replace("-", "_")


def _numeric(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    return None


# ---------------------------------------------------------------------------
# published ordering
# ---------------------------------------------------------------------------
def _published_ranks(entry: Mapping[str, Any]) -> tuple[list[float] | None, str, str]:
    """``(ranks, order_source, note)`` - rank 1 is the most potent."""
    items = list(entry.get("ordered_compounds") or [])
    if len(items) < 2:
        return None, "none", "entry lists fewer than two compounds"

    ranks = [_numeric(it.get("rank")) for it in items]
    values = [_numeric(it.get("value")) for it in items]
    readout = str(entry.get("readout") or "").lower()
    descending = any(word in readout for word in ("resistance", "selectivity", "potentiation"))

    if all(r is not None for r in ranks):
        if len(set(ranks)) < 2:
            return None, "tied", "every published rank is identical, so the ordering has no variance"
        return [float(r) for r in ranks], "published_rank", "ranks taken verbatim from the dataset"

    if sum(v is not None for v in values) >= 2:
        inactive = [
            v is None
            and any(
                word in str(it.get("efficacy") or it.get("direction") or "").lower()
                for word in ("no agonist", "no effect", "inactive")
            )
            for v, it in zip(values, items)
        ]
        if all(v is not None or inactive[i] for i, v in enumerate(values)):
            numeric = [v for v in values if v is not None]
            worst = max(numeric) if not descending else min(numeric)
            filled = [
                (worst * 10.0 if not descending else worst / 10.0) if v is None else v
                for v in values
            ]
            signed = [(-v if descending else v) for v in filled]
            note = (
                f"ranked by the published {entry.get('readout')} "
                f"({'higher' if descending else 'lower'} = more potent)"
            )
            if any(inactive):
                note += (
                    "; compounds the paper reports as having NO effect are placed last "
                    "(their value is null because there is no curve, not because it was unavailable)"
                )
            return _ranks(signed), "published_value", note

    directions = [str(it.get("direction") or "") for it in items]
    if all(directions) and len({("vertebrate" in d.lower()) for d in directions}) > 1:
        return (
            [2.0 if "vertebrate" in d.lower() else 1.0 for d in directions],
            "direction_selectivity",
            "the dataset records no numbers for this entry; the usable claim is the "
            "direction of selectivity, so insect-selective compounds are ranked above "
            "vertebrate-selective ones",
        )

    return None, "none", (
        "the dataset deliberately records no ranks or values for this entry "
        f"({entry.get('notes', '')[:120].strip()}...), so no ordering is asserted"
    )


# ---------------------------------------------------------------------------
# model ordering
# ---------------------------------------------------------------------------
def _organism_filter(entry: Mapping[str, Any]) -> str | None:
    species = str(entry.get("species") or "").lower() + " " + str(entry.get("assay") or "").lower()
    if any(w in species for w in _VERTEBRATE_WORDS) and not any(w in species for w in _INSECT_WORDS):
        return "vertebrate_"
    if any(w in species for w in _INSECT_WORDS) and not any(w in species for w in _VERTEBRATE_WORDS):
        return "insect_"
    return None


#: text in an entry's assay/readout that names a receptor directly
_RECEPTOR_HINTS: tuple[tuple[str, str], ...] = (
    ("alpha7", "vertebrate_nAChR_a7"),
    ("alpha 7", "vertebrate_nAChR_a7"),
    ("a7 ", "vertebrate_nAChR_a7"),
    ("alpha4beta2", "vertebrate_nAChR_a4b2"),
    ("glucl", "insect_GluCl"),
    ("rdl", "insect_RDL"),
    ("glycine", "vertebrate_GlyR"),
    ("sodium channel", "insect_Nav"),
)


def _hinted_receptor(entry: Mapping[str, Any], keys: Sequence[str], lib: dict[str, Any]) -> str | None:
    text = " ".join(str(entry.get(k) or "") for k in ("assay", "species", "readout")).lower()
    for token, receptor in _RECEPTOR_HINTS:
        if token in text and sum(
            1 for key in keys if receptor in (lib["compounds"][key].get("receptors") or {})
        ) >= 2:
            return receptor
    return None


def _infer_receptor(
    keys: Sequence[str], entry: Mapping[str, Any], library: dict[str, Any] | None
) -> tuple[str | None, str]:
    lib = library or load_library()
    hinted = _hinted_receptor(entry, keys, lib)
    if hinted:
        return hinted, f"{hinted} named by the published assay description"
    prefix = _organism_filter(entry)
    counts: dict[str, list[float]] = {}
    for key in keys:
        for receptor, spec in (lib["compounds"][key].get("receptors") or {}).items():
            if prefix and not receptor.startswith(prefix):
                continue
            if str(spec.get("direction", "none")) == "none":
                continue
            value = spec_value_M(spec)
            if value is None:
                continue  # schema v3 placeholder: no number to rank on
            counts.setdefault(receptor, []).append(value)
    if not counts:
        return None, "no receptor is active for these compounds in the library"
    best = max(counts, key=lambda r: (len(counts[r]), -float(np.mean(np.log10(counts[r])))))
    return best, (
        f"{best} covers {len(counts[best])}/{len(keys)} compounds"
        + (f" and matches the {prefix.rstrip('_')} preparation of the published assay" if prefix else "")
    )


_CIRCUIT_CACHE: dict[tuple, float] = {}
_VEHICLE_RATES: dict[tuple, Any] = {}


def _circuit_score(compound: str, conc_M: float, graph: str | None = None, **kw: Any) -> float:
    """Fraction of vehicle network activity displaced by the compound.

    Direction-free on purpose: the nAChR patch is non-monotone, so a signed
    readout cannot rank compounds by strength of effect.
    """
    from flylab.circuit.rate import DEFAULT_GAINS, compute_gains, rate_network

    key = (compound, conc_M, str(graph))
    if key in _CIRCUIT_CACHE:
        return _CIRCUIT_CACHE[key]
    net = rate_network(graph)
    drive = net.drive_vector(net.seed_drive(float(kw.get("drive_hz", 40.0))))
    steps = int(kw.get("steps", 80))
    veh_key = (str(graph), float(kw.get("drive_hz", 40.0)), steps)
    if veh_key not in _VEHICLE_RATES:
        _VEHICLE_RATES[veh_key] = net.run(drive, DEFAULT_GAINS, steps=steps)
    veh = _VEHICLE_RATES[veh_key]
    gains, _ = compute_gains(compound, conc_M)
    r = net.run(drive, gains, steps=steps)
    score = float(np.mean(np.abs(r - veh)) / max(float(np.mean(veh)), 1e-12))
    _CIRCUIT_CACHE[key] = score
    return score


def _model_scores(
    items: Sequence[Mapping[str, Any]],
    keys: Sequence[str],
    receptor: str | None,
    entry: Mapping[str, Any],
    assay: str,
    conc_M: float,
    order_source: str,
    library: dict[str, Any] | None,
    **kw: Any,
) -> tuple[list[float], list[str | None], str]:
    """``(scores, per-item receptor, description)``; a higher score is more potent."""
    lib = library or load_library()
    scores: list[float] = []
    used: list[str | None] = []

    if assay == "subgraph":
        for key in keys:
            scores.append(_circuit_score(key, conc_M, **kw))
            used.append(None)
        return scores, used, (
            "circuit displacement mean|r_treated - r_vehicle| / mean(r_vehicle) on the "
            "MaleCNS neighborhood rate model"
        )

    if order_source == "direction_selectivity":
        for key in keys:
            sel = compare_compound(key, conc_M, library=lib)["selectivity"]
            pair = next(
                (p for p in sel.values() if p["insect_receptor"] == receptor),
                sel.get("nAChR"),
            )
            ratio = pair.get("ec50_ratio_vert_over_insect") if pair else None
            # ratio is None when either side of the pair is a placeholder row:
            # the entry is then skipped, not scored off missing evidence.
            scores.append(float(np.log10(ratio)) if ratio else float("nan"))
            used.append(pair["insect_receptor"] if pair else None)
        return scores, used, (
            "log10(vertebrate EC50 / insect EC50) from the library's selectivity block; "
            "positive means insect-selective"
        )

    for item, key in zip(items, keys):
        target = item.get("target") or receptor
        spec = (lib["compounds"][key].get("receptors") or {}).get(target)
        value = spec_value_M(spec) if spec else None
        if not spec or value is None:
            # no sourced row (or a placeholder): NaN -> the entry is skipped
            scores.append(float("nan"))
            used.append(target)
            continue
        scores.append(
            engagement(conc_M, value, float(spec.get("n", 1.0)),
                       param_type=spec.get("param_type"), relation=spec.get("relation"))
        )
        used.append(target)
    return scores, used, f"receptor occupancy at {receptor} and {conc_M:g} M"


# ---------------------------------------------------------------------------
# per-entry validation
# ---------------------------------------------------------------------------
def validate_entry(
    entry: Mapping[str, Any] | str,
    assay: str = "occupancy",
    conc_M: float = 1e-6,
    receptor: str | None = None,
    library: dict[str, Any] | None = None,
    path: Path | None = None,
    **kw: Any,
) -> dict[str, Any]:
    """Compare the model's ordering with one published ordering.

    ``entry`` is a dataset row or its ``id``.  Returns Spearman's rho, Kendall's
    tau-b, an ``exact_match`` flag and the pairwise table; entries the library
    cannot cover come back with ``skipped: True`` and ``missing_compounds``.
    """
    if assay not in ASSAYS:
        raise ValueError(f"unknown assay {assay!r}; expected one of {ASSAYS}")
    if isinstance(entry, str):
        entry = _entry_by_id(entry, path)
    lib = library or load_library()
    known = set(lib["compounds"])

    items = list(entry.get("ordered_compounds") or [])
    base = {
        "id": entry.get("id"),
        "assay_published": entry.get("assay"),
        "species": entry.get("species"),
        "readout": entry.get("readout"),
        "use_in_flylab": entry.get("use_in_flylab"),
        "evidence_tier": entry.get("evidence_tier"),
        "source": entry.get("source"),
        "model_assay": assay,
        "conc_M": conc_M,
    }

    labels = [str(it.get("compound")) for it in items]
    keys = [_norm(it.get("compound")) for it in items]
    missing = sorted({lab for lab, key in zip(labels, keys) if key is None or key not in known})
    if missing:
        return {
            **base,
            "skipped": True,
            "reason": "library does not contain every compound in this entry",
            "missing_compounds": missing,
            "n_compounds": len(items),
        }

    published, order_source, order_note = _published_ranks(entry)
    if published is None:
        return {
            **base,
            "skipped": True,
            "reason": order_note,
            "missing_compounds": [],
            "n_compounds": len(items),
        }

    receptor_used = receptor
    inference_note = "receptor given by the caller"
    if receptor_used is None and assay == "occupancy":
        receptor_used, inference_note = _infer_receptor([k for k in keys if k], entry, lib)

    scores, per_item_receptor, score_note = _model_scores(
        items, [k for k in keys if k], receptor_used, entry, assay, conc_M, order_source, lib, **kw
    )
    if any(not np.isfinite(s) for s in scores):
        return {
            **base,
            "skipped": True,
            "reason": "the model produced no finite score for at least one compound "
            f"(receptor {receptor_used})",
            "missing_compounds": [],
            "n_compounds": len(items),
        }

    if len({round(float(s), 12) for s in scores}) < 2:
        return {
            **base,
            "skipped": True,
            "reason": (
                "the model gives every compound in this entry the same score at "
                f"{conc_M:g} M (saturated or identical), so no ordering can be compared"
            ),
            "missing_compounds": [],
            "n_compounds": len(items),
            "model_scores": [float(s) for s in scores],
        }

    # rank 1 = most potent = highest model score
    model_ranks = _ranks([-s for s in scores])
    rho = spearman_rho(published, model_ranks)
    tau = kendall_tau(published, model_ranks)
    exact = [round(r, 6) for r in published] == [round(r, 6) for r in model_ranks]

    pairwise = []
    n_conc = n_disc = 0
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            pub = np.sign(published[j] - published[i])  # +1 -> i is more potent
            mod = np.sign(model_ranks[j] - model_ranks[i])
            agree = bool(pub == mod) if pub and mod else None
            if agree is True:
                n_conc += 1
            elif agree is False:
                n_disc += 1
            pairwise.append(
                {
                    "a": labels[i],
                    "b": labels[j],
                    "published": "a>b" if pub > 0 else ("b>a" if pub < 0 else "tie"),
                    "model": "a>b" if mod > 0 else ("b>a" if mod < 0 else "tie"),
                    "agree": agree,
                }
            )

    # Secondary ordering by EC50 alone. Occupancy far below EC50 is sensitive to
    # the Hill coefficient, so the two can disagree for reasons that are an
    # artefact of n rather than a pharmacological claim; when they do, say so.
    ec50_block: dict[str, Any] = {}
    if assay == "occupancy" and order_source != "direction_selectivity":
        ec50s = []
        for item, key in zip(items, [k for k in keys if k]):
            target = item.get("target") or receptor_used
            spec = (lib["compounds"][key].get("receptors") or {}).get(target)
            value = spec_value_M(spec) if spec else None
            ec50s.append(value if value is not None else float("inf"))
        if all(math.isfinite(e) for e in ec50s) and len(set(ec50s)) > 1:
            ec50_ranks = _ranks(ec50s)  # lower EC50 = rank 1 = more potent
            ec50_block = {
                "ec50_M": ec50s,
                "model_rank_by_ec50": ec50_ranks,
                "spearman_rho_by_ec50": spearman_rho(published, ec50_ranks),
                "agrees_with_occupancy_ordering": [round(r, 6) for r in ec50_ranks]
                == [round(r, 6) for r in model_ranks],
            }
            if not ec50_block["agrees_with_occupancy_ordering"]:
                ec50_block["note"] = (
                    f"At {conc_M:g} M the occupancy ordering differs from the EC50 ordering. "
                    "Occupancy well below EC50 is dominated by the Hill coefficient "
                    f"(n = {[round(float((lib['compounds'][k].get('receptors') or {}).get(item.get('target') or receptor_used, {}).get('n', 1.0)), 2) for k, item in zip([x for x in keys if x], items)]}), "
                    "so this disagreement is an artefact of n, not a potency claim. "
                    "Compare spearman_rho_by_ec50 as well."
                )

    direction_conflicts = []
    if order_source == "direction_selectivity":
        for i, item in enumerate(items):
            published_vertebrate = "vertebrate" in str(item.get("direction") or "").lower()
            model_vertebrate = scores[i] < 0
            if published_vertebrate != model_vertebrate:
                direction_conflicts.append(
                    {
                        "compound": labels[i],
                        "published_direction": item.get("direction"),
                        "model_log10_ratio_vert_over_insect": scores[i],
                        "note": "the model's selectivity SIGN disagrees with the published direction "
                        "even where the ordering agrees",
                    }
                )

    return {
        **base,
        "skipped": False,
        "n_compounds": len(items),
        "ec50_ordering": ec50_block,
        "direction_conflicts": direction_conflicts,
        "receptor": receptor_used,
        "receptor_note": inference_note,
        "order_source": order_source,
        "order_note": order_note,
        "model_score_note": score_note,
        "compounds": [
            {
                "compound": labels[i],
                "key": keys[i],
                "target": per_item_receptor[i],
                "published_rank": published[i],
                "published_value": items[i].get("value"),
                "model_score": scores[i],
                "model_rank": model_ranks[i],
            }
            for i in range(len(items))
        ],
        "spearman_rho": rho,
        "kendall_tau": tau,
        "exact_match": exact,
        "pairwise": pairwise,
        "n_concordant": n_conc,
        "n_discordant": n_disc,
        "warnings": [QUANTITY_NOTE] if "bioassay" in str(entry.get("assay", "")).lower() else [],
    }


# ---------------------------------------------------------------------------
# known discrepancies
# ---------------------------------------------------------------------------
KNOWN_DISCREPANCIES = (
    {
        "id": "nitenpyram_rank_inverted",
        "compounds": ["nitenpyram", "clothianidin", "imidacloprid", "acetamiprid", "thiamethoxam"],
        "claim": (
            "Nitenpyram is near the top of the model's insect_nAChR potency ordering but is "
            "the LEAST potent neonicotinoid in Drosophila whole-animal bioassays "
            "(8x weaker than clothianidin in larvae, 25x weaker in adults)."
        ),
        "published_source": "Perry T, Heckel DG, McKenzie JA, Batterham P (2012) Pestic Biochem Physiol 102(1):56-60",
        "resolution": "reported, not corrected",
    },
    {
        "id": "clothianidin_vs_imidacloprid_inverted",
        "compounds": ["clothianidin", "imidacloprid"],
        "claim": (
            "The library makes imidacloprid (2.0e-8 M) marginally more potent than "
            "clothianidin (3.0e-8 M) at insect_nAChR; both the acute Drosophila bioassay "
            "and the chronic-lifespan ordering put clothianidin first. The gap is inside "
            "the uncertainty of either number but the sign is opposite."
        ),
        "published_source": (
            "Perry et al. (2012) Pestic Biochem Physiol 102:56; "
            "Tasman K, Rands SA, Hodge JJL (2021) Front Physiol 12:659440"
        ),
        "resolution": "reported, not corrected",
    },
)


def _discrepancy_report(conc_M: float, library: dict[str, Any] | None = None) -> dict[str, Any]:
    lib = library or load_library()
    known = set(lib["compounds"])
    out = []
    for spec in KNOWN_DISCREPANCIES:
        keys = [c for c in spec["compounds"] if c in known]
        rows = []
        for key in keys:
            entry = (lib["compounds"][key].get("receptors") or {}).get("insect_nAChR")
            if not entry or str(entry.get("direction", "none")) == "none":
                continue
            value = spec_value_M(entry)
            if value is None:
                continue  # placeholder row: excluded from the ordering
            eng = engagement(conc_M, value, float(entry.get("n", 1.0)),
                             param_type=entry.get("param_type"), relation=entry.get("relation"))
            rows.append(
                {
                    "compound": key,
                    "param_type": entry.get("param_type"),
                    "param_value_M": value,
                    "ec50_M": value,
                    "engagement": eng,
                    "occupancy": eng,
                }
            )
        rows.sort(key=lambda r: r["ec50_M"])
        model_order = [r["compound"] for r in rows]
        block = {
            **{k: v for k, v in spec.items()},
            "model_order_by_insect_nAChR_ec50": model_order,
            "model_rows": rows,
            "confirmed_present": False,
        }
        if spec["id"] == "nitenpyram_rank_inverted" and "nitenpyram" in model_order:
            position = model_order.index("nitenpyram") + 1
            block["model_rank_of_nitenpyram"] = position
            block["published_rank_of_nitenpyram"] = len(model_order)
            block["confirmed_present"] = position < len(model_order)
            block["detail"] = (
                f"model ranks nitenpyram {position} of {len(model_order)} by insect_nAChR EC50; "
                f"Perry et al. 2012 rank it last ({len(model_order)} of {len(model_order)})."
            )
        if spec["id"] == "clothianidin_vs_imidacloprid_inverted" and {"clothianidin", "imidacloprid"} <= set(model_order):
            inverted = model_order.index("imidacloprid") < model_order.index("clothianidin")
            block["confirmed_present"] = bool(inverted)
            block["detail"] = (
                "model order " + " < ".join(model_order) + " (by EC50); published order puts "
                "clothianidin before imidacloprid. Inverted: " + str(inverted)
            )
        out.append(block)
    return {
        "discrepancies": out,
        "explanation": QUANTITY_NOTE,
        "policy": (
            "The library is NOT tuned to remove these. Both are consequences of comparing a "
            "receptor-level EC50 with a whole-animal potency, which is the explanation as "
            "well as the caveat."
        ),
    }


# ---------------------------------------------------------------------------
# whole dataset
# ---------------------------------------------------------------------------
def validate_all(
    conc_M: float = 1e-6,
    assays: Iterable[str] = ("occupancy", "subgraph"),
    library: dict[str, Any] | None = None,
    path: Path | None = None,
    **kw: Any,
) -> dict[str, Any]:
    """Validate every dataset entry under each named assay.

    Returns ``{"assays": {name: {...}}, "summary": {...},
    "known_discrepancies": {...}}``.  ``summary`` carries ``n_entries``,
    ``n_evaluated``, ``n_skipped`` and ``mean_rho`` (over the assays' evaluated
    entries with a finite rho).
    """
    data = load_rank_orders(path)
    entries = data.get("entries", [])
    assays = list(assays)
    per_assay: dict[str, Any] = {}
    all_rho: list[float] = []
    n_eval = n_skip = 0
    for assay in assays:
        rows = [validate_entry(e, assay=assay, conc_M=conc_M, library=library, path=path, **kw) for e in entries]
        rhos = [r["spearman_rho"] for r in rows if not r["skipped"] and np.isfinite(r["spearman_rho"])]
        evaluated = [r for r in rows if not r["skipped"]]
        per_assay[assay] = {
            "rows": rows,
            "n_entries": len(rows),
            "n_evaluated": len(evaluated),
            "n_skipped": len(rows) - len(evaluated),
            "n_exact": sum(1 for r in evaluated if r.get("exact_match")),
            "mean_rho": float(np.mean(rhos)) if rhos else None,
            "rho_by_entry": {r["id"]: r["spearman_rho"] for r in evaluated},
        }
        all_rho.extend(rhos)
        n_eval += len(evaluated)
        n_skip += len(rows) - len(evaluated)

    assay_names = list(assays)
    return {
        "conc_M": conc_M,
        "assays": per_assay,
        "summary": {
            "n_entries": len(entries),
            "n_assays": len(assay_names),
            "n_runs": len(entries) * len(assay_names),
            "n_evaluated": n_eval,
            "n_skipped": n_skip,
            "mean_rho": float(np.mean(all_rho)) if all_rho else None,
        },
        "known_discrepancies": _discrepancy_report(conc_M, library),
        "not_sourced": data.get("not_sourced", []),
        "warnings": [
            QUANTITY_NOTE,
            "These entries are validation TARGETS. Nothing in flylab/pharm/library.yaml "
            "is fitted to them, and mismatches are reported rather than removed.",
        ],
    }


def validation_report_markdown(
    conc_M: float = 1e-6,
    assays: Iterable[str] = ("occupancy", "subgraph"),
    library: dict[str, Any] | None = None,
    path: Path | None = None,
    **kw: Any,
) -> str:
    """The retrospective-validation table, for the paper."""
    result = validate_all(conc_M, assays, library, path, **kw)
    lines = [
        "# FlyLab retrospective validation",
        "",
        f"Model concentration: {conc_M:g} M. Source dataset: "
        "`data/literature/published_rank_orders.yaml`.",
        "",
        "| entry | published assay | species | n | model assay | rho | tau | exact | note |",
        "|---|---|---|---:|---|---:|---:|:-:|---|",
    ]
    for assay, block in result["assays"].items():
        for row in block["rows"]:
            if row["skipped"]:
                note = f"skipped: {row['reason']}"
                if row.get("missing_compounds"):
                    note += " (" + ", ".join(row["missing_compounds"]) + ")"
                lines.append(
                    f"| {row['id']} | {row.get('assay_published', '')} | {row.get('species', '')} | "
                    f"{row.get('n_compounds', 0)} | {assay} | - | - | - | {note} |"
                )
            else:
                lines.append(
                    f"| {row['id']} | {row.get('assay_published', '')} | {row.get('species', '')} | "
                    f"{row['n_compounds']} | {assay} | {row['spearman_rho']:.3f} | "
                    f"{row['kendall_tau']:.3f} | {'yes' if row['exact_match'] else 'no'} | "
                    f"{row['order_source']}; {row['model_score_note']} |"
                )
    summary = result["summary"]
    lines += [
        "",
        f"**Summary.** {summary['n_evaluated']} evaluated, {summary['n_skipped']} skipped, "
        f"mean Spearman rho = "
        + (f"{summary['mean_rho']:.3f}" if summary["mean_rho"] is not None else "n/a")
        + ".",
        "",
        "## Known discrepancies (reported, not corrected)",
        "",
    ]
    for block in result["known_discrepancies"]["discrepancies"]:
        lines += [
            f"### {block['id']}",
            "",
            block["claim"],
            "",
            f"- Model: {block.get('detail', 'n/a')}",
            f"- Published: {block['published_source']}",
            f"- Present in this build: {'yes' if block['confirmed_present'] else 'no'}",
            "",
        ]
    lines += ["", result["known_discrepancies"]["explanation"], "", result["known_discrepancies"]["policy"], ""]
    return "\n".join(lines)
