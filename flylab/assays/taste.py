from __future__ import annotations
from typing import Any
from flylab.circuit.reduced_taste import CIRCUIT_CITATION, CIRCUIT_ID, run_taste_circuit
from flylab.notebook.schema import empty_notebook
from flylab.pharm.occupancy import compare_compound, load_library

def _ach_gain(compound: str | None, conc_M: float) -> tuple[float, dict[str, Any] | None]:
    if not compound:
        return 1.0, None
    occ = compare_compound(compound, conc_M)
    insect = next(r for r in occ["receptors"] if r["receptor"] == "insect_nAChR")
    direction = insect["direction"]
    occupancy = insect["occupancy"]
    if direction == "agonist":
        gain = max(0.05, 1.0 + 0.4 * occupancy - 1.6 * occupancy**2)
    elif direction == "antagonist":
        gain = max(0.05, 1.0 - occupancy)
    else:
        gain = 1.0
    return gain, occ

def run_taste_assay(compound: str | None = None, conc_M: float = 0.0, sugar_hz: float = 150.0, bitter_hz: float = 0.0) -> dict[str, Any]:
    g_ach, occ = _ach_gain(compound, conc_M)
    sugar = run_taste_circuit(sugar_drive_hz=sugar_hz, bitter_drive_hz=0.0, g_ach=g_ach)
    both = run_taste_circuit(sugar_drive_hz=sugar_hz, bitter_drive_hz=150.0, g_ach=g_ach)
    vehicle_sugar = run_taste_circuit(sugar_drive_hz=sugar_hz, bitter_drive_hz=0.0, g_ach=1.0)
    notebook = empty_notebook("taste_mn9")
    notebook["map"] = {"name": CIRCUIT_ID, "version": "v0", "citation": CIRCUIT_CITATION}
    notebook["compound"] = compound
    notebook["concentration_M"] = conc_M
    notebook["occupancy"] = occ["receptors"] if occ else []
    notebook["readouts"] = {
        "g_ach": g_ach,
        "mn9_sugar_hz": sugar.mn9_hz,
        "mn9_sugar_bitter_hz": both.mn9_hz,
        "mn9_vehicle_sugar_hz": vehicle_sugar.mn9_hz,
        "bitter_veto_ratio": (both.mn9_hz / sugar.mn9_hz if sugar.mn9_hz > 1e-6 else None),
        "cells": {
            "sweet_grn": sugar.sweet_grn,
            "bitter_grn": both.bitter_grn,
            "second_order_exc": sugar.second_order_exc,
            "second_order_inh": both.second_order_inh,
            "mn9": sugar.mn9_hz,
        },
    }
    notebook["warnings"] = [
        "Circuit is reduced_taste_v0, not a FlyWire/MaleCNS synapse extract.",
        "g_ach is a phenomenological patch from insect nAChR occupancy.",
        "Do not treat MN9 Hz as a measured spike rate from a living fly.",
    ]
    if occ:
        notebook["warnings"].append(occ["disclaimer"])
    return notebook

def dose_response(compound: str, concs: list[float], sugar_hz: float = 150.0) -> list[dict[str, Any]]:
    rows = []
    for c in concs:
        nb = run_taste_assay(compound=compound, conc_M=c, sugar_hz=sugar_hz, bitter_hz=0.0)
        rows.append({"concentration_M": c, "mn9_hz": nb["readouts"]["mn9_sugar_hz"], "g_ach": nb["readouts"]["g_ach"]})
    return rows

def library_keys() -> list[str]:
    return sorted(load_library()["compounds"])
