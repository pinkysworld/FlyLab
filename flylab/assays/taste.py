from __future__ import annotations

from typing import Any

from flylab.circuit.reduced_taste import CIRCUIT_CITATION, CIRCUIT_ID, run_taste_circuit
from flylab.notebook.schema import empty_notebook, set_map
from flylab.pharm.mechanisms import default_gains, gains_from_occupancy
from flylab.pharm.occupancy import compare_compound, load_library

#: Gains the reduced taste circuit can actually inject.
_SUPPORTED_GAINS = ("g_ach",)


def _gains(compound: str | None, conc_M: float) -> tuple[dict[str, float], dict[str, Any] | None]:
    """Mechanism gains for one dose (single source of truth: pharm.mechanisms)."""
    if not compound:
        return default_gains(), None
    occ = compare_compound(compound, conc_M)
    return gains_from_occupancy(occ["receptors"]), occ


def _ach_gain(compound: str | None, conc_M: float) -> tuple[float, dict[str, Any] | None]:
    """Backward-compatible helper: cholinergic gain only."""
    gains, occ = _gains(compound, conc_M)
    return gains["g_ach"], occ


def run_taste_assay(compound: str | None = None, conc_M: float = 0.0, sugar_hz: float = 150.0, bitter_hz: float = 0.0) -> dict[str, Any]:
    gains, occ = _gains(compound, conc_M)
    g_ach = gains["g_ach"]
    sugar = run_taste_circuit(sugar_drive_hz=sugar_hz, bitter_drive_hz=0.0, g_ach=g_ach)
    both = run_taste_circuit(sugar_drive_hz=sugar_hz, bitter_drive_hz=150.0, g_ach=g_ach)
    vehicle_sugar = run_taste_circuit(sugar_drive_hz=sugar_hz, bitter_drive_hz=0.0, g_ach=1.0)
    notebook = empty_notebook("taste_mn9")
    set_map(notebook, CIRCUIT_ID, "v0", CIRCUIT_CITATION)
    notebook["compound"] = compound
    notebook["concentration_M"] = conc_M
    notebook["occupancy"] = occ["receptors"] if occ else []
    notebook["gains"] = gains
    notebook["readouts"] = {
        "g_ach": g_ach,
        "mn9_sugar_hz": sugar.mn9_hz,
        "mn9_sugar_bitter_hz": both.mn9_hz,
        "mn9_vehicle_sugar_hz": vehicle_sugar.mn9_hz,
        "bitter_veto_ratio": (both.mn9_hz / sugar.mn9_hz if sugar.mn9_hz > 1e-6 else None),
        "gains": gains,
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
    ignored = [k for k, v in gains.items() if k not in _SUPPORTED_GAINS and abs(v - 1.0) > 1e-9]
    if ignored:
        notebook["warnings"].append(
            "reduced_taste_v0 injects only g_ach; "
            + ", ".join(sorted(ignored))
            + " are reported in `gains` but not applied to this circuit."
        )
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
