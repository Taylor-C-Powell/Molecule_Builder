"""Process engineering: reactors, solvents, purification, costing, safety.

This package provides modules for reactor selection, solvent scoring,
purification strategy, reaction condition optimisation, cost estimation,
safety assessment, and scale-up analysis.
"""

# -- reactor --
from molbuilder.process.reactor import (
    ReactorType,
    ReactorSpec,
    select_reactor,
)

# -- solvent systems --
from molbuilder.process.solvent_systems import (
    SolventRecommendation,
    select_solvent,
)

# -- purification --
from molbuilder.process.purification import (
    PurificationMethod,
    PurificationStep,
    recommend_purification,
)

# -- conditions --
from molbuilder.process.conditions import (
    ReactionConditions,
    optimize_conditions,
)

# -- costing --
from molbuilder.process.costing import (
    CostBreakdown,
    CostEstimate,
    estimate_cost,
)

# -- safety --
from molbuilder.process.safety import (
    GHS_PICTOGRAMS,
    GHS_HAZARD_STATEMENTS,
    HazardInfo,
    SafetyAssessment,
    assess_safety,
)

# -- scale-up --
from molbuilder.process.scale_up import (
    ScaleUpAnalysis,
    analyze_scale_up,
)

__all__ = [
    # reactor
    "ReactorType",
    "ReactorSpec",
    "select_reactor",
    # solvent_systems
    "SolventRecommendation",
    "select_solvent",
    # purification
    "PurificationMethod",
    "PurificationStep",
    "recommend_purification",
    # conditions
    "ReactionConditions",
    "optimize_conditions",
    # costing
    "CostBreakdown",
    "CostEstimate",
    "estimate_cost",
    # safety
    "GHS_PICTOGRAMS",
    "GHS_HAZARD_STATEMENTS",
    "HazardInfo",
    "SafetyAssessment",
    "assess_safety",
    # scale_up
    "ScaleUpAnalysis",
    "analyze_scale_up",
]
