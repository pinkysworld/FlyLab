"""Analysis layer: dose-response fits, drug impact maps, viewer layout,
connectome null models, Connectome-Dependence Analysis, the model ablation
ladder, circuit selectivity and prospective predictions."""
from flylab.analysis.baselines import (
    LEVELS,
    ablation,
    ablation_table,
    level_predictions,
)
from flylab.analysis.dependence import (
    CLASSES,
    DEFAULT_MODES,
    INFORMATION_LADDER,
    VERDICTS,
    balance_report,
    benjamini_hochberg,
    cut_census,
    dependence_landscape,
    dependence_profile,
    equivalence_margin,
    ladder_power,
    ladder_recovery,
    mode_verdict,
    necessary_information_level,
    synthetic_cut,
)
from flylab.analysis.fit import bootstrap_fit, hill4, hill4_fit
from flylab.analysis.impact import (
    edge_impact,
    grn_to_mn9_paths,
    node_impact,
    path_impact,
    summarize_impact,
)
from flylab.analysis.layout import graph_for_viewer, graph_layout
from flylab.analysis.nullmodels import (
    MODES,
    connectome_information_score,
    convergence_report,
    empirical_p,
    null_distribution,
    null_panel,
    p_resolution,
    shuffle_graph,
)
from flylab.analysis.predictions import (
    HYPOTHESES,
    power_for_continuous,
    power_for_per,
    predict,
    prediction_table,
    to_markdown,
)
from flylab.analysis.selectivity import (
    circuit_selectivity_index,
    circuit_threshold_conc,
    graph_composition,
    receptor_selectivity_table,
    selectivity_landscape,
)
from flylab.analysis.robustness import (
    CONCLUSIONS,
    MECHANISM_FAMILY,
    SHAPES,
    conclusion_stability,
    family_table,
    make_spec,
    mechanism_spec,
    spec_gains,
    threshold_sensitivity,
)
from flylab.analysis.uncertainty_global import (
    FACTORS,
    ISHIGAMI_REFERENCE,
    RESOLUTION_STATES,
    RateSurrogate,
    ishigami,
    jansen_indices,
    saltelli_matrices,
    resolution_summary,
    sobol_analysis,
    uncertainty_budget,
)
from flylab.analysis.voi import EXPERIMENTS, VOI_WARNING, value_of_information

__all__ = [
    "hill4",
    "hill4_fit",
    "bootstrap_fit",
    "edge_impact",
    "node_impact",
    "path_impact",
    "grn_to_mn9_paths",
    "summarize_impact",
    "graph_layout",
    "graph_for_viewer",
    # null models
    "MODES",
    "shuffle_graph",
    "null_distribution",
    "null_panel",
    "connectome_information_score",
    "empirical_p",
    "p_resolution",
    "convergence_report",
    # connectome-dependence analysis
    "dependence_profile",
    "dependence_landscape",
    "necessary_information_level",
    "INFORMATION_LADDER",
    "DEFAULT_MODES",
    "CLASSES",
    "VERDICTS",
    "mode_verdict",
    "equivalence_margin",
    "benjamini_hochberg",
    "cut_census",
    "balance_report",
    "synthetic_cut",
    "ladder_recovery",
    "ladder_power",
    # ablation ladder
    "LEVELS",
    "level_predictions",
    "ablation",
    "ablation_table",
    # selectivity
    "receptor_selectivity_table",
    "circuit_threshold_conc",
    "circuit_selectivity_index",
    "selectivity_landscape",
    "graph_composition",
    # predictions
    "HYPOTHESES",
    "predict",
    "power_for_per",
    "power_for_continuous",
    "prediction_table",
    "to_markdown",
    # mechanism-rule robustness
    "SHAPES",
    "MECHANISM_FAMILY",
    "CONCLUSIONS",
    "make_spec",
    "spec_gains",
    "mechanism_spec",
    "family_table",
    "conclusion_stability",
    "threshold_sensitivity",
    # global uncertainty attribution
    "FACTORS",
    "RateSurrogate",
    "saltelli_matrices",
    "jansen_indices",
    "sobol_analysis",
    "uncertainty_budget",
    "resolution_summary",
    "RESOLUTION_STATES",
    "ishigami",
    "ISHIGAMI_REFERENCE",
    # value of information (exported under a name that does not shadow the
    # flylab.analysis.voi submodule)
    "value_of_information",
    "EXPERIMENTS",
    "VOI_WARNING",
]
