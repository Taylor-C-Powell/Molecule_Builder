"""Process costing for multi-step synthesis routes.

Estimates total manufacturing cost by summing raw-material, labour,
equipment, energy, waste-disposal, and overhead contributions across
all synthesis steps.  Reagent pricing comes from ``REAGENT_DB``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, List

from molbuilder.reactions.reaction_types import ReactionCategory, ReactionTemplate
from molbuilder.reactions.reagent_data import REAGENT_DB, SOLVENT_DB, get_reagent, get_solvent


# =====================================================================
#  Data classes
# =====================================================================

@dataclass
class CostBreakdown:
    """Itemised cost breakdown in USD."""

    raw_materials_usd: float
    labor_usd: float
    equipment_usd: float
    energy_usd: float
    waste_disposal_usd: float
    overhead_usd: float


@dataclass
class CostEstimate:
    """Complete cost estimate for a synthesis route."""

    total_usd: float
    per_kg_usd: float
    breakdown: CostBreakdown
    scale_kg: float
    notes: list[str]


# =====================================================================
#  Configuration constants (realistic 2024 values)
# =====================================================================

# Labour rate USD/h (loaded, including benefits)
_LABOR_RATE_USD_H = 75.0

# Hours of labour per synthesis step (depends on complexity)
_BASE_LABOR_HOURS_PER_STEP = 3.0

# Energy cost per kWh
_ENERGY_COST_KWH = 0.10

# Waste disposal cost per kg
_WASTE_DISPOSAL_PER_KG = 2.50

# Overhead multiplier on direct costs (insurance, QC, facility)
_OVERHEAD_FRACTION = 0.25

# Equipment depreciation rate (fraction of capital cost per batch)
_EQUIPMENT_DEPRECIATION_PER_BATCH = 0.002

# Default reagent equivalents if not otherwise specified
_DEFAULT_REAGENT_EQUIV = 1.2

# Typical solvent volume: litres per kg of product
_SOLVENT_L_PER_KG = 8.0


# =====================================================================
#  Internal helpers
# =====================================================================

def _normalise_key(name: str) -> str:
    """Normalise a reagent/solvent name to a REAGENT_DB lookup key."""
    return name.lower().replace(" ", "_").replace("-", "_")


def _reagent_cost_for_step(template: ReactionTemplate, scale_kg: float) -> float:
    """Estimate raw-material cost for one step.

    Looks up each reagent in REAGENT_DB.  Falls back to a conservative
    default of $100/kg if the reagent is not found.
    """
    total = 0.0
    for rname in template.reagents:
        reagent = get_reagent(rname)
        price_per_kg = reagent.cost_per_kg if reagent and reagent.cost_per_kg > 0 else 100.0
        # Assume reagent mass ~ DEFAULT_EQUIV * product mass (rough stoichiometry)
        reagent_kg = scale_kg * _DEFAULT_REAGENT_EQUIV
        total += price_per_kg * reagent_kg

    # Catalyst cost (used in smaller amounts: ~0.05 equiv)
    for cname in template.catalysts:
        cat_reagent = get_reagent(cname)
        cat_price = cat_reagent.cost_per_kg if cat_reagent and cat_reagent.cost_per_kg > 0 else 500.0
        catalyst_kg = scale_kg * 0.05
        total += cat_price * catalyst_kg

    # Solvent cost
    for sname in template.solvents[:1]:  # primary solvent only
        solvent = get_solvent(sname)
        cost_per_L = solvent.cost_per_L if solvent and solvent.cost_per_L > 0 else 15.0
        solvent_L = scale_kg * _SOLVENT_L_PER_KG
        total += cost_per_L * solvent_L

    return total


def _labor_cost_for_step(template: ReactionTemplate, scale_kg: float) -> float:
    """Estimate labour cost for one step.

    Base 3 h per step.  Add time for:
    - Cryogenic work (+1 h)
    - Chromatography purification (+2 h)
    - Scale > 50 kg (+1 h for logistics)
    """
    hours = _BASE_LABOR_HOURS_PER_STEP

    mean_t = (template.temperature_range[0] + template.temperature_range[1]) / 2.0
    if mean_t < -20:
        hours += 1.0
    if template.category in {ReactionCategory.REARRANGEMENT, ReactionCategory.MISC}:
        hours += 1.0   # complexity
    if scale_kg > 50:
        hours += 1.0   # material handling
    if scale_kg > 200:
        hours += 1.0   # additional logistics

    return hours * _LABOR_RATE_USD_H


def _equipment_cost_for_step(template: ReactionTemplate, scale_kg: float) -> float:
    """Depreciation-based equipment cost per batch.

    Uses six-tenths rule for vessel capital cost then applies per-batch
    depreciation fraction.
    """
    volume_L = scale_kg * 6.0  # rough vessel sizing
    base_capital = 15_000.0    # reference 100 L vessel
    scale_factor = (max(volume_L, 1.0) / 100.0) ** 0.6
    capital = base_capital * max(scale_factor, 0.3)
    return capital * _EQUIPMENT_DEPRECIATION_PER_BATCH


def _energy_cost_for_step(template: ReactionTemplate, scale_kg: float) -> float:
    """Estimate energy consumption in USD.

    Heating/cooling energy based on temperature delta, agitation power,
    and distillation energy for solvent removal.
    """
    mean_t = (template.temperature_range[0] + template.temperature_range[1]) / 2.0
    delta_t = abs(mean_t - 20.0)  # from ambient

    # Heating or cooling energy (kWh) ~ mass * Cp * deltaT / 3600
    # Approximate solution Cp ~ 2 kJ/(kg*K), plus safety factor
    mass_kg = scale_kg * _SOLVENT_L_PER_KG  # total charge mass approx
    energy_kwh = mass_kg * 2.0 * delta_t / 3600.0

    # Agitation: ~0.5 kW per 100 L for 2 h
    volume_L = scale_kg * 6.0
    agitation_kwh = (volume_L / 100.0) * 0.5 * 2.0

    # Solvent removal (rotovap / distillation): ~0.3 kWh per L
    solvent_removal_kwh = scale_kg * _SOLVENT_L_PER_KG * 0.3

    total_kwh = energy_kwh + agitation_kwh + solvent_removal_kwh
    return total_kwh * _ENERGY_COST_KWH


def _waste_cost_for_step(template: ReactionTemplate, scale_kg: float) -> float:
    """Estimate waste disposal cost.

    Waste includes spent solvent, aqueous washes, and by-products.
    Typically 5-10x product mass.
    """
    # Hazardous waste multiplier
    hazardous = False
    for rname in template.reagents:
        reagent = get_reagent(rname)
        if reagent:
            dangerous_codes = {"H300", "H310", "H330", "H340", "H350", "H360"}
            if set(reagent.ghs_hazards) & dangerous_codes:
                hazardous = True
                break

    waste_multiplier = 8.0 if hazardous else 5.0
    waste_kg = scale_kg * waste_multiplier
    disposal_rate = _WASTE_DISPOSAL_PER_KG * (2.0 if hazardous else 1.0)
    return waste_kg * disposal_rate


# =====================================================================
#  Public API
# =====================================================================

def estimate_cost(
    steps: list[Any],
    scale_kg: float,
) -> CostEstimate:
    """Estimate total manufacturing cost for a synthesis route.

    Parameters
    ----------
    steps : list
        Each element is expected to have a ``.template`` attribute
        (:class:`ReactionTemplate`) and a ``.precursors`` attribute
        (used for context but not directly priced -- precursor costs are
        captured through the reagent lookup on each template).
        Duck typing is used; no specific class is required.
    scale_kg : float
        Target production scale in kilograms of final product.

    Returns
    -------
    CostEstimate
        Complete cost estimate with per-kg pricing and itemised breakdown.
    """
    if scale_kg <= 0:
        scale_kg = 0.001  # avoid division by zero

    total_materials = 0.0
    total_labor = 0.0
    total_equipment = 0.0
    total_energy = 0.0
    total_waste = 0.0
    notes: list[str] = []

    num_steps = len(steps)
    if num_steps == 0:
        notes.append("No synthesis steps provided; returning zero cost.")
        breakdown = CostBreakdown(0, 0, 0, 0, 0, 0)
        return CostEstimate(
            total_usd=0.0, per_kg_usd=0.0,
            breakdown=breakdown, scale_kg=scale_kg, notes=notes,
        )

    # Cumulative yield loss: each step reduces effective scale
    cumulative_scale = scale_kg
    for step in reversed(steps):
        template = step.template
        # Average yield for this step
        avg_yield = (template.typical_yield[0] + template.typical_yield[1]) / 2.0 / 100.0
        if avg_yield <= 0:
            avg_yield = 0.5
        # Scale needed at this step to deliver cumulative_scale
        step_scale = cumulative_scale / avg_yield

        total_materials += _reagent_cost_for_step(template, step_scale)
        total_labor += _labor_cost_for_step(template, step_scale)
        total_equipment += _equipment_cost_for_step(template, step_scale)
        total_energy += _energy_cost_for_step(template, step_scale)
        total_waste += _waste_cost_for_step(template, step_scale)

        cumulative_scale = step_scale

    # Overhead
    direct_costs = total_materials + total_labor + total_equipment + total_energy + total_waste
    total_overhead = direct_costs * _OVERHEAD_FRACTION

    total_usd = direct_costs + total_overhead

    # Notes
    notes.append(f"Costing based on {num_steps} synthesis step(s) at {scale_kg:.2f} kg scale.")
    if scale_kg < 0.1:
        notes.append(
            "Lab-scale pricing; per-kg cost will decrease significantly at larger scale."
        )
    if scale_kg > 100:
        notes.append(
            "Bulk pricing may apply for reagents; actual costs could be 20-40% lower."
        )
    if any(
        get_reagent(r) and get_reagent(r).cost_per_kg > 5000  # type: ignore[union-attr]
        for step in steps for r in step.template.reagents
    ):
        notes.append(
            "One or more expensive reagents detected (>$5000/kg); "
            "consider catalyst recycling or alternative chemistry."
        )

    breakdown = CostBreakdown(
        raw_materials_usd=round(total_materials, 2),
        labor_usd=round(total_labor, 2),
        equipment_usd=round(total_equipment, 2),
        energy_usd=round(total_energy, 2),
        waste_disposal_usd=round(total_waste, 2),
        overhead_usd=round(total_overhead, 2),
    )

    return CostEstimate(
        total_usd=round(total_usd, 2),
        per_kg_usd=round(total_usd / scale_kg, 2),
        breakdown=breakdown,
        scale_kg=scale_kg,
        notes=notes,
    )
