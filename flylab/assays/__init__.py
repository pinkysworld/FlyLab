"""FlyLab assays: each returns a notebook (schema 0.3)."""
from flylab.assays.ensemble import circuit_ic50, run_ensemble, sensitivity
from flylab.assays.fullcns import estimate_fullcns_cost, load_full_cns, run_fullcns_assay
from flylab.assays.experiment import ExperimentDesign, design_from_yaml, run_experiment
from flylab.assays.spiking import run_spiking_assay
from flylab.assays.subgraph import dose_response_subgraph, run_subgraph_assay
from flylab.assays.taste import dose_response, run_taste_assay
from flylab.assays.taste_map import run_taste_map_assay

__all__ = [
    "run_fullcns_assay",
    "load_full_cns",
    "estimate_fullcns_cost",
    "run_taste_assay",
    "dose_response",
    "run_subgraph_assay",
    "dose_response_subgraph",
    "run_spiking_assay",
    "run_taste_map_assay",
    "run_ensemble",
    "sensitivity",
    "circuit_ic50",
    "run_experiment",
    "design_from_yaml",
    "ExperimentDesign",
]
