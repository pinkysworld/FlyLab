"""Pre-registered predictions, power and the taste library hook."""
from __future__ import annotations

import json
import math
import time

import pytest

from flylab.analysis.predictions import (
    HYPOTHESES,
    behavioral_baseline,
    hypothesis,
    power_for_continuous,
    power_for_per,
    predict,
    prediction_table,
    to_markdown,
)


def test_seven_hypotheses_are_declared():
    assert list(HYPOTHESES) == ["H1", "H2", "H3", "H4", "H5", "H6", "H7"]
    for h_id, h in HYPOTHESES.items():
        assert h["statement"] and h["assay"] and h["readout"]
        assert h["direction"] in ("increase", "decrease", "none")
        assert h["h_id"] == h_id
    assert HYPOTHESES["H6"]["compound"] == "fipronil"
    assert HYPOTHESES["H6"]["readout"] == "bitter_veto_ratio"
    assert HYPOTHESES["H7"]["compound"] == "picrotoxin"
    with pytest.raises(KeyError):
        hypothesis("H99")


def test_power_for_per_is_a_sensible_sample_size():
    res = power_for_per(0.3, baseline_p=0.8)
    assert 20 <= res["n_per_group"] <= 80
    assert res["p1"] == pytest.approx(0.8)
    assert res["p2"] == pytest.approx(0.56)
    assert res["label"] == "model_derived" and res["warnings"]
    # a smaller effect needs more flies
    assert power_for_per(0.15, baseline_p=0.8)["n_per_group"] > res["n_per_group"]
    assert power_for_per(0.0, baseline_p=0.8)["n_per_group"] is None


def test_power_for_continuous_matches_the_closed_form():
    res = power_for_continuous(1.0)
    assert res["n_per_group"] == pytest.approx(math.ceil(2 * (1.959963985 + 0.8416212) ** 2), abs=1)
    assert power_for_continuous(0.5)["n_per_group"] > res["n_per_group"]
    assert power_for_continuous(0.0)["n_per_group"] is None
    # model-internal d is capped before the calculation
    capped = power_for_continuous(50.0)
    assert capped["d_used"] == 1.0 and capped["n_per_group"] == res["n_per_group"]


def test_behavioral_baseline_falls_back_and_says_so(tmp_path):
    base = behavioral_baseline()
    assert 0.0 < base["baseline_per"] <= 1.0
    assert base["sd"] > 0
    assert base["warnings"]
    assert base["source"]
    junk = tmp_path / "behavioral_assays.yaml"
    junk.write_text("this: [is, not, a, per, table]\n")
    fallback = behavioral_baseline(junk)
    assert fallback["baseline_per"] == 0.8 and fallback["sd"] == 0.15
    assert any("default" in w for w in fallback["warnings"])


def test_behavioral_baseline_reads_a_plausible_yaml(tmp_path):
    p = tmp_path / "behavioral_assays.yaml"
    p.write_text("assays:\n  per:\n    baseline_proportion: 0.72\n    sd: 0.11\n")
    base = behavioral_baseline(p)
    assert base["loaded"] is True
    assert base["baseline_per"] == pytest.approx(0.72)
    assert base["sd"] == pytest.approx(0.11)


@pytest.mark.parametrize("h_id", ["H1", "H2", "H6"])
def test_predict_returns_an_effect_a_ci_and_a_labelled_d(h_id):
    res = predict(h_id, n_rep=3, seed=0)
    assert res["predicted_effect"] is not None
    assert math.isfinite(res["predicted_effect"])
    assert len(res["ci"]) == 2 and res["ci"][0] <= res["ci"][1]
    assert res["d_meaning"] == "model_internal_variability"
    assert res["status"] == "software_prediction"
    assert res["live_result"] is None
    assert res["direction_matches"], f"{h_id} went the wrong way: {res['observed_direction']}"
    assert res["label"] == "model_derived" and res["warnings"]
    json.dumps(res)


def test_prediction_table_has_seven_rows_with_finite_effects():
    t0 = time.perf_counter()
    table = prediction_table(n_rep=4, seed=0)
    # Wall-clock budget: this guards against an accidentally quadratic
    # implementation, it is not a benchmark. CI runners and parallel test
    # runs are slow and highly variable, so the bound is deliberately loose.
    assert time.perf_counter() - t0 < 300.0
    rows = table["rows"]
    assert len(rows) == 7
    assert [r["h_id"] for r in rows] == ["H1", "H2", "H3", "H4", "H5", "H6", "H7"]
    for r in rows:
        assert r["predicted_effect"] is not None and math.isfinite(r["predicted_effect"])
        assert r["status"] == "software_prediction"
        assert r["live_result"] is None
        if r["h_id"] == "H4":
            # predicted null: no effect, so no sample size
            assert r["predicted_effect"] == pytest.approx(0.0, abs=1e-9)
            assert r["suggested_n_per_group"] is None
        else:
            assert r["suggested_n_per_group"] is not None
            assert 1 <= r["suggested_n_per_group"] <= 10000
    md = to_markdown(table)
    assert md.count("\n|") >= 8
    assert "H7" in md and "software_prediction" in md
    assert to_markdown(rows) == md
    json.dumps(table)


def test_taste_assay_accepts_a_library_override():
    """`run_taste_assay(library=...)` must thread into compare_compound."""
    import copy

    from flylab.assays.taste import run_taste_assay
    from flylab.pharm.occupancy import load_library

    lib = copy.deepcopy(load_library())
    spec = lib["compounds"]["imidacloprid"]["receptors"]["insect_nAChR"]
    spec["ec50_M"] = float(spec["ec50_M"]) * 1e4  # far less potent

    base = run_taste_assay("imidacloprid", 1e-6)
    weak = run_taste_assay("imidacloprid", 1e-6, library=lib)
    assert weak["gains"]["g_ach"] > base["gains"]["g_ach"]
    assert weak["readouts"]["mn9_sugar_hz"] != base["readouts"]["mn9_sugar_hz"]
    # unchanged behaviour when no library is passed
    assert run_taste_assay("imidacloprid", 1e-6)["readouts"] == base["readouts"]
