"""Report generation: synthesis, safety, costing, molecule summaries.

Public API
----------
Text formatting utilities:
    section_header, subsection_header, ascii_table, word_wrap,
    bullet_list, key_value_block, horizontal_bar,
    format_currency, format_percent

Report generators:
    generate_molecule_report  -- comprehensive molecule analysis
    generate_synthesis_report -- multi-step synthesis route
    generate_safety_report    -- hazard and PPE assessment
    generate_cost_report      -- cost breakdown with charts

Protocol types:
    SynthesisStepLike   -- structural type for synthesis step objects
    SynthesisRouteLike  -- structural type for synthesis route objects
    CostEstimateLike    -- structural type for cost estimation objects
    SafetyAssessmentLike -- structural type for safety assessment objects
    MoleculeLike        -- structural type for molecule-like objects
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable

from molbuilder.reports.text_formatter import (
    section_header,
    subsection_header,
    ascii_table,
    word_wrap,
    bullet_list,
    key_value_block,
    horizontal_bar,
    format_currency,
    format_percent,
)

from molbuilder.reports.molecule_report import generate_molecule_report
from molbuilder.reports.synthesis_report import generate_synthesis_report
from molbuilder.reports.safety_report import generate_safety_report
from molbuilder.reports.cost_report import generate_cost_report


# =====================================================================
#  Protocol types for duck-typed interfaces
# =====================================================================

@runtime_checkable
class SynthesisStepLike(Protocol):
    """Structural type for synthesis step objects.

    Used by reports/ and process/ modules that accept step objects
    with a .template attribute (e.g. SynthesisStep from synthesis_route).
    """
    step_number: int
    template: object   # ReactionTemplate (avoid circular import)
    precursors: list
    product_smiles: str
    product_name: str
    expected_yield: float


@runtime_checkable
class SynthesisRouteLike(Protocol):
    """Structural type for synthesis route objects.

    Used by generate_synthesis_report() and process engineering modules.
    """
    target_smiles: str
    target_name: str
    steps: list
    overall_yield: float
    starting_materials: list
    total_steps: int
    longest_linear_sequence: int


@runtime_checkable
class CostBreakdownLike(Protocol):
    """Structural type for the cost breakdown sub-object."""
    raw_materials_usd: float
    labor_usd: float
    equipment_usd: float
    energy_usd: float
    waste_disposal_usd: float
    overhead_usd: float


@runtime_checkable
class CostEstimateLike(Protocol):
    """Structural type for cost estimate objects.

    Used by generate_cost_report().
    """
    total_usd: float
    per_kg_usd: float
    scale_kg: float
    breakdown: CostBreakdownLike
    notes: list


@runtime_checkable
class SafetyAssessmentLike(Protocol):
    """Structural type for per-step safety assessment objects.

    Used by generate_safety_report().
    """
    step_number: int
    step_name: str
    hazards: list
    ppe_required: list
    engineering_controls: list
    emergency_procedures: list
    incompatible_materials: list
    waste_classification: str
    risk_level: str


@runtime_checkable
class MoleculeLike(Protocol):
    """Structural type for molecule-like objects.

    Used by generate_molecule_report() and functional group detectors.
    """
    name: str
    atoms: list
    bonds: list

    def neighbors(self, idx: int) -> list[int]: ...
    def get_bond(self, i: int, j: int) -> object | None: ...


__all__ = [
    # Protocols
    "SynthesisStepLike",
    "SynthesisRouteLike",
    "CostBreakdownLike",
    "CostEstimateLike",
    "SafetyAssessmentLike",
    "MoleculeLike",
    # Formatters
    "section_header",
    "subsection_header",
    "ascii_table",
    "word_wrap",
    "bullet_list",
    "key_value_block",
    "horizontal_bar",
    "format_currency",
    "format_percent",
    # Report generators
    "generate_molecule_report",
    "generate_synthesis_report",
    "generate_safety_report",
    "generate_cost_report",
]
