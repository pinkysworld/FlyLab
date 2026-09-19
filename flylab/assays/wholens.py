from __future__ import annotations
from typing import Any
from flylab.maps.malecns import MAP_CITATION, MAP_ID, load_census
from flylab.notebook.schema import empty_notebook
from flylab.pharm.occupancy import compare_compound

def _gains(compound: str | None, conc_M: float):
    g = {k: 1.0 for k in ("acetylcholine", "gaba", "glutamate", "other")}
    if not compound:
        return g, None
    occ = compare_compound(compound, conc_M)
    by = {r["receptor"]: r for r in occ["receptors"]}
    ach = by.get("insect_nAChR")
    if ach:
        th, d = ach["occupancy"], ach["direction"]
        if d == "agonist":
            g["acetylcholine"] = max(0.05, 1.0 + 0.4 * th - 1.6 * th * th)
        elif d == "antagonist":
            g["acetylcholine"] = max(0.05, 1.0 - th)
    gaba = by.get("vertebrate_GABA_A")
    if gaba and gaba["direction"] in {"antagonist", "positive_modulator"}:
        th = gaba["occupancy"]
        g["gaba"] = max(0.05, 1.0 - th) if gaba["direction"] == "antagonist" else 1.0 + 0.8 * th
    return g, occ

def _settle(n, gains):
    tot = max(sum(n.values()), 1)
    drive = {k: 40.0 * (v / tot) for k, v in n.items()}
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
    gains, occ = _gains(compound, conc_M)
    treated = _settle(n, gains)
    vehicle = _settle(n, {k: 1.0 for k in gains})
    nb = empty_notebook("whole_cns_census")
    nb["map"] = {"name": MAP_ID, "version": "v1.0", "citation": MAP_CITATION}
    nb["compound"] = compound
    nb["concentration_M"] = conc_M
    nb["occupancy"] = occ["receptors"] if occ else []
    nb["readouts"] = {
        "n_traced": census["n_traced"], "neurotransmitter_counts": n,
        "superclass_counts": census["superclass_counts"], "named_cells": census["named_cells"],
        "gains": gains, "vehicle": vehicle, "treated": treated,
        "delta_excitation": treated["cns_excitation_index"] - vehicle["cns_excitation_index"],
        "weights_present": census["weights_present"],
    }
    nb["warnings"] = [
        "Whole-CNS layer uses MaleCNS traced-neuron census + predicted transmitters.",
        "It is not a synapse-resolved 166k-cell LIF simulation.",
        "Insect GABA-Cl (RDL/GluCl) is only approximated; library GABA-A is vertebrate.",
        "Do not treat excitation_index as a measured firing rate.",
    ]
    if not census["weights_present"]:
        nb["warnings"].append("Weight matrix not downloaded. Run: flylab download-malecns --full")
    if occ:
        nb["warnings"].append(occ["disclaimer"])
    return nb
