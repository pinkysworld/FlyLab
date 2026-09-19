from __future__ import annotations

from typing import Any

NOTEBOOK_VERSION = "0.1"


def empty_notebook(assay: str) -> dict[str, Any]:
    return {
        "flylab_notebook_version": NOTEBOOK_VERSION,
        "assay": assay,
        "map": {"name": None, "version": None, "citation": None},
        "compound": None,
        "concentration_M": None,
        "occupancy": [],
        "readouts": {},
        "live_lab": None,
        "warnings": [
            "Circuit runtime not attached.",
            "Do not cite occupancy-only output as a connectome result.",
        ],
    }
