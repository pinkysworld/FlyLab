"""Circuit runtimes: reduced taste control, rate network, LIF network."""
from flylab.circuit.lif import LIFNetwork, SpikeResult, lif_network
from flylab.circuit.rate import (
    DEFAULT_GAINS,
    GRAPHS,
    SIGN,
    RateNetwork,
    compute_gains,
    load_graph,
    rate_network,
    resolve_graph,
)
from flylab.circuit.reduced_taste import CIRCUIT_ID, run_taste_circuit

__all__ = [
    "run_taste_circuit",
    "CIRCUIT_ID",
    "RateNetwork",
    "rate_network",
    "LIFNetwork",
    "SpikeResult",
    "lif_network",
    "SIGN",
    "GRAPHS",
    "DEFAULT_GAINS",
    "compute_gains",
    "load_graph",
    "resolve_graph",
]
