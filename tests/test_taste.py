from flylab.assays.taste import run_taste_assay
from flylab.circuit.reduced_taste import run_taste_circuit


def test_bitter_vetoes_mn9():
    sugar = run_taste_circuit(sugar_drive_hz=150, bitter_drive_hz=0)
    both = run_taste_circuit(sugar_drive_hz=150, bitter_drive_hz=150)
    assert sugar.mn9_hz > 40
    assert both.mn9_hz < 0.25 * sugar.mn9_hz


def test_high_imidacloprid_suppresses_mn9_vs_vehicle():
    veh = run_taste_assay(compound=None, conc_M=0.0)
    treated = run_taste_assay(compound="imidacloprid", conc_M=1e-5)
    assert treated["readouts"]["mn9_sugar_hz"] < veh["readouts"]["mn9_sugar_hz"]
    assert treated["map"]["name"] == "reduced_taste_v0"


def test_notebook_warns_not_connectome():
    nb = run_taste_assay("nicotine", 1e-6)
    assert any("not a FlyWire" in w for w in nb["warnings"])
