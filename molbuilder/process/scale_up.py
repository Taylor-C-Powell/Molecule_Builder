"""Scale-up analysis for transitioning from lab to production.

Evaluates batch vs. continuous manufacturing, estimates cycle times,
annual capacity, capital costs, and identifies scale-up risks.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from molbuilder.reactions.reaction_types import ReactionCategory, ReactionTemplate
from molbuilder.process import DEFAULT_SOLVENT_L_PER_KG


# =====================================================================
#  Data class
# =====================================================================

@dataclass
class ScaleUpAnalysis:
    """Complete scale-up analysis result."""

    target_annual_kg: float
    recommended_mode: str       # "batch", "continuous", "semi-continuous"
    batch_size_kg: float | None
    batches_per_year: int | None
    cycle_time_hours: float
    annual_capacity_kg: float
    capital_cost_usd: float
    operating_cost_annual_usd: float
    scale_up_risks: list[str]
    recommendations: list[str]


# =====================================================================
#  Configuration constants
# =====================================================================

# Annual operating hours (24/7 with ~15% downtime for maintenance)
_ANNUAL_OPERATING_HOURS = 7500.0

# Maximum practical batch size for glass-lined reactors (litres)
_MAX_BATCH_VOLUME_L = 16_000.0

# Typical product density assumption (kg/L)
_PRODUCT_DENSITY = 1.0

# Solvent-to-product volume ratio (shared constant)
_SOLVENT_RATIO = DEFAULT_SOLVENT_L_PER_KG

# Batch turnaround time (cleaning, charging, discharging) in hours
_TURNAROUND_HOURS = 2.0

# Labour rate for production operators (USD/h)
_OPERATOR_RATE_USD_H = 55.0

# Number of operators per shift
_OPERATORS_PER_SHIFT = 2

# Shifts per day for continuous operation
_SHIFTS_PER_DAY = 3

# Utility cost per operating hour (steam, cooling water, electricity)
_UTILITY_COST_PER_HOUR = 45.0

# Maintenance as fraction of capital cost per year
_MAINTENANCE_FRACTION = 0.05

# Raw material cost placeholder (USD per kg product)
_RAW_MATERIAL_COST_PER_KG = 200.0


# =====================================================================
#  Internal helpers
# =====================================================================

def _estimate_cycle_time(steps: list[Any]) -> float:
    """Estimate total cycle time in hours for one batch through all steps.

    Includes reaction time, workup, and turnaround per step.
    """
    total = 0.0
    for step in steps:
        template: ReactionTemplate = step.template
        mean_t = (template.temperature_range[0] + template.temperature_range[1]) / 2.0

        # Base reaction time estimate from yield and temperature
        _, yield_hi = template.typical_yield
        rxn_hours = 1.0 + (100.0 - yield_hi) * 0.04

        # Temperature correction: cryogenic requires cool-down, slow
        # addition, and warm-up (Arrhenius: lower T -> slower kinetics)
        if mean_t < -20:
            rxn_hours *= 1.8
        elif mean_t > 120:
            rxn_hours *= 0.6

        # Workup time estimate
        workup_hours = 1.5
        if template.category in {ReactionCategory.COUPLING, ReactionCategory.OXIDATION}:
            workup_hours = 2.0
        if template.category == ReactionCategory.POLYMERIZATION:
            workup_hours = 3.0

        total += rxn_hours + workup_hours + _TURNAROUND_HOURS

    return round(max(total, 2.0), 1)


def _max_batch_kg() -> float:
    """Maximum practical batch size in kg of product."""
    # Max vessel volume / solvent ratio gives product mass
    return _MAX_BATCH_VOLUME_L / _SOLVENT_RATIO / _PRODUCT_DENSITY


def _is_continuous_candidate(steps: list[Any], target_annual_kg: float) -> bool:
    """Decide if continuous processing is feasible.

    Continuous is preferred when:
    - Annual volume > 50,000 kg
    - All steps are relatively fast (no step > 8 h)
    - No solid-handling steps that complicate continuous flow
    """
    if target_annual_kg < 50_000:
        return False

    for step in steps:
        template: ReactionTemplate = step.template
        _, yield_hi = template.typical_yield
        rxn_hours = 1.0 + (100.0 - yield_hi) * 0.04
        if rxn_hours > 8.0:
            return False
        # Solid handling makes continuous harder
        if template.category in {ReactionCategory.POLYMERIZATION}:
            return False

    return True


def _capital_cost_batch(batch_size_kg: float, num_steps: int) -> float:
    """Estimate capital cost for a batch plant.

    Uses the six-tenths rule scaled from a reference 100 kg plant
    costing $500k per step.
    """
    ref_cost_per_step = 500_000.0
    scale_factor = (max(batch_size_kg, 1.0) / 100.0) ** 0.6
    return round(ref_cost_per_step * max(scale_factor, 0.2) * num_steps, -3)


def _capital_cost_continuous(target_annual_kg: float, num_steps: int) -> float:
    """Estimate capital cost for a continuous plant.

    Continuous plants are more expensive per step but have lower
    operating costs at high throughput.
    """
    ref_cost_per_step = 800_000.0
    throughput_factor = (max(target_annual_kg, 1000.0) / 100_000.0) ** 0.5
    return round(ref_cost_per_step * max(throughput_factor, 0.3) * num_steps, -3)


def _identify_risks(steps: list[Any], target_annual_kg: float, mode: str) -> list[str]:
    """Identify scale-up risks for the route."""
    risks: list[str] = []

    for idx, step in enumerate(steps):
        template: ReactionTemplate = step.template
        mean_t = (template.temperature_range[0] + template.temperature_range[1]) / 2.0

        # Heat transfer risk
        if mean_t < -40 or mean_t > 150:
            risks.append(
                f"Step {idx+1} ({template.name}): Extreme temperature ({mean_t:.0f} C) "
                "-- heat transfer may be limiting at scale"
            )

        # Mixing sensitivity
        if template.category in {ReactionCategory.COUPLING, ReactionCategory.RADICAL}:
            risks.append(
                f"Step {idx+1} ({template.name}): Mixing-sensitive reaction "
                "-- verify mass transfer at larger vessel dimensions"
            )

        # Exothermic risk
        if template.category in {ReactionCategory.ADDITION, ReactionCategory.POLYMERIZATION}:
            risks.append(
                f"Step {idx+1} ({template.name}): Potentially exothermic "
                "-- calorimetry data required; ensure adequate cooling capacity"
            )

        # Chromatography at scale
        if template.category in {ReactionCategory.REARRANGEMENT, ReactionCategory.MISC}:
            if target_annual_kg > 1000:
                risks.append(
                    f"Step {idx+1} ({template.name}): Complex mixture expected "
                    "-- chromatographic purification impractical at this scale"
                )

        # Safety
        if template.safety_notes:
            risks.append(
                f"Step {idx+1} ({template.name}): Safety concern -- {template.safety_notes}"
            )

    # General risks
    if target_annual_kg > 10_000:
        risks.append(
            "Regulatory: Production at >10 t/year may require environmental "
            "impact assessment and process safety review (PSM/COMAH)"
        )
    if mode == "continuous":
        risks.append(
            "Continuous operation: Requires robust process analytical technology "
            "(PAT) and automated control systems"
        )
    if len(steps) > 5:
        risks.append(
            f"Long linear sequence ({len(steps)} steps): cumulative yield loss "
            "is significant; consider convergent synthesis"
        )

    return risks


def _generate_recommendations(
    steps: list[Any],
    target_annual_kg: float,
    mode: str,
    batch_size_kg: float | None,
) -> list[str]:
    """Generate actionable scale-up recommendations."""
    recs: list[str] = []

    if mode == "batch":
        recs.append(
            "Perform reaction calorimetry (RC1) for each step to characterise "
            "thermal hazard and design cooling systems"
        )
        recs.append(
            "Conduct solvent-swap studies to identify process-friendly solvents "
            "that simplify workup and waste treatment"
        )
        if batch_size_kg and batch_size_kg > 500:
            recs.append(
                "At >500 kg batch size, consider in-line analytics (FTIR, Raman) "
                "for real-time reaction monitoring"
            )

    if mode == "continuous":
        recs.append(
            "Develop residence time distribution (RTD) model for each "
            "continuous unit operation"
        )
        recs.append(
            "Implement process analytical technology (PAT) for real-time "
            "quality control"
        )
        recs.append(
            "Design start-up and shutdown procedures; define off-spec "
            "material handling protocols"
        )

    if mode == "semi-continuous":
        recs.append(
            "Identify which steps benefit from continuous operation "
            "(e.g. fast reactions, dangerous intermediates) and which "
            "are better run batchwise"
        )

    # Universal recommendations
    recs.append(
        "Conduct HAZOP study before commissioning"
    )
    recs.append(
        "Validate analytical methods for in-process control and "
        "final product release"
    )
    if target_annual_kg > 1000:
        recs.append(
            "Evaluate catalyst recycling and solvent recovery to "
            "reduce operating costs and environmental impact"
        )

    return recs


# =====================================================================
#  Public API
# =====================================================================

def analyze_scale_up(
    steps: list[Any],
    target_annual_kg: float,
) -> ScaleUpAnalysis:
    """Analyse scale-up feasibility for *steps* at *target_annual_kg*.

    Parameters
    ----------
    steps : list
        Each element must have a ``.template`` attribute
        (:class:`ReactionTemplate`) and a ``.precursors`` attribute.
        Duck typing is used.
    target_annual_kg : float
        Desired annual production in kilograms.

    Returns
    -------
    ScaleUpAnalysis
        Comprehensive scale-up analysis including mode recommendation,
        batch sizing, capital and operating cost estimates, risks, and
        recommendations.
    """
    if not steps:
        return ScaleUpAnalysis(
            target_annual_kg=target_annual_kg,
            recommended_mode="batch",
            batch_size_kg=None,
            batches_per_year=None,
            cycle_time_hours=0.0,
            annual_capacity_kg=0.0,
            capital_cost_usd=0.0,
            operating_cost_annual_usd=0.0,
            scale_up_risks=[],
            recommendations=[],
        )
    for i, step in enumerate(steps):
        if not hasattr(step, 'template'):
            raise TypeError(
                f"Step {i} must have a 'template' attribute, "
                f"got {type(step).__name__}"
            )

    if target_annual_kg <= 0:
        target_annual_kg = 1.0

    num_steps = len(steps)
    cycle_time = _estimate_cycle_time(steps)

    # --- Mode selection ---
    if _is_continuous_candidate(steps, target_annual_kg):
        mode = "continuous"
    elif target_annual_kg > 10_000 and num_steps <= 3:
        mode = "semi-continuous"
    else:
        mode = "batch"

    # --- Batch sizing ---
    batch_size_kg: float | None = None
    batches_per_year: int | None = None

    if mode in ("batch", "semi-continuous"):
        max_kg = _max_batch_kg()
        # Size batch to produce target with reasonable number of batches
        batches_needed_min = max(1, math.ceil(target_annual_kg / max_kg))
        available_batches = int(_ANNUAL_OPERATING_HOURS / cycle_time)
        if available_batches < 1:
            available_batches = 1

        if batches_needed_min <= available_batches:
            batches_per_year = batches_needed_min
            batch_size_kg = round(target_annual_kg / batches_per_year, 1)
        else:
            # Cannot meet target with single train; use max batch size
            batches_per_year = available_batches
            batch_size_kg = round(max_kg, 1)

    # --- Annual capacity ---
    if mode == "continuous":
        # Continuous throughput based on cycle time per step and parallel capacity
        hourly_throughput_kg = target_annual_kg / _ANNUAL_OPERATING_HOURS
        annual_capacity_kg = round(hourly_throughput_kg * _ANNUAL_OPERATING_HOURS * 1.1, 0)
    else:
        if batches_per_year and batch_size_kg:
            annual_capacity_kg = round(batches_per_year * batch_size_kg, 0)
        else:
            annual_capacity_kg = target_annual_kg

    # --- Capital cost ---
    if mode == "continuous":
        capital_cost = _capital_cost_continuous(target_annual_kg, num_steps)
    else:
        capital_cost = _capital_cost_batch(batch_size_kg or 100.0, num_steps)

    # --- Operating cost ---
    # Labour
    if mode == "continuous":
        annual_labor_hours = _ANNUAL_OPERATING_HOURS * _OPERATORS_PER_SHIFT
    else:
        batches = batches_per_year or 1
        annual_labor_hours = batches * cycle_time * _OPERATORS_PER_SHIFT

    labor_cost = annual_labor_hours * _OPERATOR_RATE_USD_H

    # Utilities
    if mode == "continuous":
        utility_cost = _ANNUAL_OPERATING_HOURS * _UTILITY_COST_PER_HOUR
    else:
        batches = batches_per_year or 1
        utility_cost = batches * cycle_time * _UTILITY_COST_PER_HOUR

    # Raw materials
    raw_material_cost = target_annual_kg * _RAW_MATERIAL_COST_PER_KG

    # Maintenance
    maintenance_cost = capital_cost * _MAINTENANCE_FRACTION

    operating_cost_annual = round(
        labor_cost + utility_cost + raw_material_cost + maintenance_cost, -2
    )

    # --- Risks and recommendations ---
    risks = _identify_risks(steps, target_annual_kg, mode)
    recommendations = _generate_recommendations(steps, target_annual_kg, mode, batch_size_kg)

    return ScaleUpAnalysis(
        target_annual_kg=target_annual_kg,
        recommended_mode=mode,
        batch_size_kg=batch_size_kg,
        batches_per_year=batches_per_year,
        cycle_time_hours=cycle_time,
        annual_capacity_kg=annual_capacity_kg,
        capital_cost_usd=capital_cost,
        operating_cost_annual_usd=operating_cost_annual,
        scale_up_risks=risks,
        recommendations=recommendations,
    )
