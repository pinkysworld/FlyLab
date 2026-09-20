from __future__ import annotations

from typing import Any

from flylab.maps.malecns import MAP_CITATION, MAP_ID, load_census
from flylab.notebook.schema import empty_notebook, set_map
from flylab.pharm.mechanisms import default_gains, gains_from_occupancy
from flylab.pharm.occupancy import compare_compound

#: v0.4 set the GABA gain from the vertebrate GABA-A row. v0.5 uses insect RDL.
GABA_SOURCE_CHANGE = (
    "GABA gain now comes from insect RDL occupancy; FlyLab v0.4 used the "
    "vertebrate GABA-A row for this insect assay. Vertebrate-only GABA ligands "
    "(e.g. diazepam) therefore no longer move the fly GABA gain."
)


def _gains(compound: str | None, conc_M: float):
    """Transmitter gains for the census model, from pharm.mechanisms."""
    mech = default_gains()
    if not compound:
        return _transmitter_gains(mech), None, mech
    occ = compare_compound(compound, conc_M)
    mech = gains_from_occupancy(occ["receptors"])
    return _transmitter_gains(mech), occ, mech


def _transmitter_gains(mech: dict[str, float]) -> dict[str, float]:
    return {
        "acetylcholine": mech["g_ach"],
        "gaba": mech["g_gaba"],
        "glutamate": mech["g_glu"],
        "octopamine": mech["g_oct"],
        "other": 1.0,
    }


def _settle(n, gains, g_nav: float = 1.0):
    tot = max(sum(n.values()), 1)
    drive = {k: 40.0 * g_nav * (v / tot) for k, v in n.items()}
    ach = drive.get("acetylcholine", 0.0) * gains["acetylcholine"]
    gaba = drive.get("gaba", 0.0) * gains["gaba"]
    glu = drive.get("glutamate", 0.0) * gains.get("glutamate", 1.0)
    excitation = max(0.0, ach + 0.3 * glu - 1.1 * gaba)
    inhibition = max(0.0, gaba + 0.7 * glu)
    return {
        "acetylcholine_drive": ach, "gaba_drive": gaba, "glutamate_drive": glu,
        "cns_excitation_index": excitation, "cns_inhibition_index": inhibition,
        "excitation_inhibition_ratio": excitation / inhibition if inhibition > 1e-9 else None,
    }


def run_wholens_assay(compound=None, conc_M=0.0, dest=None):
    census = load_census(dest)
    n = {k: int(v) for k, v in census["neurotransmitter_counts"].items()}
    gains, occ, mech = _gains(compound, conc_M)
    treated = _settle(n, gains, mech["g_nav"])
    vehicle = _settle(n, {k: 1.0 for k in gains}, 1.0)
    nb = empty_notebook("whole_cns_census")
    set_map(nb, MAP_ID, "v1.0", MAP_CITATION)
    nb["compound"] = compound
    nb["concentration_M"] = conc_M
    nb["occupancy"] = occ["receptors"] if occ else []
    nb["gains"] = mech
    nb["readouts"] = {
        "n_traced": census["n_traced"], "neurotransmitter_counts": n,
        "superclass_counts": census["superclass_counts"], "named_cells": census["named_cells"],
        "gains": gains, "mechanism_gains": mech, "vehicle": vehicle, "treated": treated,
        "delta_excitation": treated["cns_excitation_index"] - vehicle["cns_excitation_index"],
        "weights_present": census["weights_present"],
    }
    nb["warnings"] = [
        "Whole-CNS layer uses MaleCNS traced-neuron census + predicted transmitters.",
        "It is not a synapse-resolved 166k-cell LIF simulation.",
        "Insect GABA-Cl (RDL) and GluCl gains are phenomenological patches.",
        "Do not treat excitation_index as a measured firing rate.",
    ]
    if not census["weights_present"]:
        nb["warnings"].append("Weight matrix not downloaded. Run: flylab download-malecns --full")
    if occ:
        nb["warnings"].append(GABA_SOURCE_CHANGE)
        by = {r["receptor"]: r for r in occ["receptors"]}
        vert = by.get("vertebrate_GABA_A")
        insect = by.get("insect_RDL")
        vert_engagement = vert.get("engagement", vert.get("occupancy")) if vert else None
        vert_active = bool(
            vert and vert["direction"] != "none"
            and vert_engagement is not None and vert_engagement > 0.01
        )
        insect_active = bool(insect and insect["direction"] != "none")
        if vert_active and not insect_active:
            nb["warnings"].append(
                f"{occ['compound']} acts on vertebrate GABA-A in this library but has no insect RDL row, "
                "so its GABA gain is 1.0 in this fly assay."
            )
        nb["warnings"].append(occ["disclaimer"])
    return nb
