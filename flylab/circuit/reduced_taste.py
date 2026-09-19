"""Reduced taste-to-MN9 circuit until a map extract lands.

Qualitative control after Shiu et al., Nature 2024. Not a connectome extract.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

CIRCUIT_ID = "reduced_taste_v0"
CIRCUIT_CITATION = (
    "Qualitative structure after Shiu et al., Nature 634, 2024 "
    "(sugar GRN drive of MN9; bitter veto). Weights are not a connectome extract."
)

@dataclass
class TasteState:
    sweet_grn: float
    bitter_grn: float
    second_order_exc: float
    second_order_inh: float
    mn9_hz: float

def _relu(x: float) -> float:
    return float(x) if x > 0.0 else 0.0

def run_taste_circuit(
    sugar_drive_hz: float = 150.0,
    bitter_drive_hz: float = 0.0,
    g_ach: float = 1.0,
    g_gaba: float = 1.0,
    dt_s: float = 0.002,
    t_s: float = 0.4,
) -> TasteState:
    if g_ach < 0 or g_gaba < 0:
        raise ValueError("gains must be >= 0")
    tau = 0.03
    w_s, w_b, w_e, w_i, scale_mn9 = 0.012, 0.014, 0.85, 1.15, 220.0
    s = b = e = i = m = 0.0
    for _ in range(max(int(t_s / dt_s), 1)):
        s += dt_s * (-s + sugar_drive_hz) / tau
        b += dt_s * (-b + bitter_drive_hz) / tau
        e += dt_s * (-e + w_s * s * g_ach) / tau
        i += dt_s * (-i + w_b * b * g_gaba) / tau
        target = scale_mn9 * _relu(w_e * e * g_ach - w_i * i)
        m += dt_s * (-m + target) / tau
    return TasteState(float(s), float(b), float(e), float(i), float(np.clip(m, 0.0, 300.0)))
