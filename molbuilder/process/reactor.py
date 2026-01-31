"""Reactor selection and specification for process engineering.

Maps reaction characteristics (exothermicity, phase, scale) to an appropriate
reactor type and provides a fully populated :class:`ReactorSpec` dataclass.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Optional

from molbuilder.reactions.reaction_types import ReactionCategory, ReactionTemplate


# =====================================================================
#  Enums
# =====================================================================

class ReactorType(Enum):
    BATCH = auto()
    SEMI_BATCH = auto()
    CSTR = auto()          # Continuous stirred-tank reactor
    PFR = auto()           # Plug flow reactor
    MICROREACTOR = auto()
    FIXED_BED = auto()


# =====================================================================
#  Reactor specification
# =====================================================================

@dataclass
class ReactorSpec:
    """Complete specification of a process reactor."""

    reactor_type: ReactorType
    volume_L: float
    temperature_C: float
    pressure_atm: float
    residence_time_min: float
    mixing_type: str        # "mechanical", "static", "none"
    heat_transfer: str      # "jacketed", "coil", "adiabatic"
    material: str           # "glass", "stainless steel", "hastelloy"
    estimated_cost_usd: float
    notes: str


# =====================================================================
#  Internal helpers
# =====================================================================

# Categories that are typically highly exothermic or fast
_EXOTHERMIC_CATEGORIES = {
    ReactionCategory.ADDITION,
    ReactionCategory.RADICAL,
    ReactionCategory.POLYMERIZATION,
}

# Categories that are often multiphase (solid catalyst, heterogeneous)
_MULTIPHASE_CATEGORIES = {
    ReactionCategory.COUPLING,
    ReactionCategory.REDUCTION,
}

# Categories commonly run under inert or pressurised conditions
_HIGH_PRESSURE_CATEGORIES = {
    ReactionCategory.REDUCTION,
    ReactionCategory.POLYMERIZATION,
}


def _mean_temp(template: ReactionTemplate) -> float:
    """Return the midpoint of the template temperature range."""
    lo, hi = template.temperature_range
    return (lo + hi) / 2.0


def _is_cryogenic(template: ReactionTemplate) -> bool:
    return _mean_temp(template) < -20.0


def _is_high_temp(template: ReactionTemplate) -> bool:
    return _mean_temp(template) > 150.0


def _select_material(template: ReactionTemplate) -> str:
    """Choose vessel material based on temperature and corrosive reagents."""
    corrosive_keywords = {"hcl", "h2so4", "hno3", "hf", "tfa", "socl2", "ticl4"}
    reagent_keys = {r.lower().replace(" ", "_").replace("-", "_")
                    for r in template.reagents}
    if reagent_keys & corrosive_keywords:
        return "hastelloy"
    if _is_high_temp(template):
        return "stainless steel"
    return "glass"


def _estimate_reactor_cost(reactor_type: ReactorType, volume_L: float) -> float:
    """Rough capital cost estimate in USD.

    Based on typical 2024 equipment pricing for chemical process vessels.
    """
    base_costs = {
        ReactorType.BATCH: 8_000,
        ReactorType.SEMI_BATCH: 12_000,
        ReactorType.CSTR: 25_000,
        ReactorType.PFR: 30_000,
        ReactorType.MICROREACTOR: 50_000,
        ReactorType.FIXED_BED: 35_000,
    }
    base = base_costs.get(reactor_type, 10_000)
    # Scale by volume using the six-tenths rule (cost ~ volume^0.6)
    reference_volume = 100.0  # litres
    if volume_L <= 0:
        volume_L = 1.0
    scale_factor = (volume_L / reference_volume) ** 0.6
    return round(base * max(scale_factor, 0.3), -2)


def _volume_for_scale(scale_kg: float, concentration_factor: float = 5.0) -> float:
    """Estimate vessel volume in litres.

    Assumes ~5 L of solvent+reagent per kg of product (tuneable via
    *concentration_factor*) and a 75 % fill level.
    """
    raw = scale_kg * concentration_factor
    return round(raw / 0.75, 1)


# =====================================================================
#  Public API
# =====================================================================

def select_reactor(
    template: ReactionTemplate,
    scale_kg: float,
) -> ReactorSpec:
    """Select an appropriate reactor for *template* at *scale_kg*.

    Decision tree
    -------------
    1. Fast / exothermic at small scale  -> MICROREACTOR
    2. Fast / exothermic at large scale  -> CSTR with jacketed cooling
    3. Slow multiphase or catalytic       -> BATCH (small) or FIXED_BED (large)
    4. High-volume commodity (>500 kg)    -> PFR
    5. Moderate scale, needs controlled   -> SEMI_BATCH
       addition
    6. Default                            -> BATCH
    """
    mean_t = _mean_temp(template)
    volume = _volume_for_scale(scale_kg)
    material = _select_material(template)

    is_exothermic = template.category in _EXOTHERMIC_CATEGORIES
    is_multiphase = template.category in _MULTIPHASE_CATEGORIES
    has_catalyst = len(template.catalysts) > 0

    # --- 1. Fast exothermic, small scale -> microreactor ---
    if is_exothermic and scale_kg < 1.0:
        rt = ReactorType.MICROREACTOR
        vol = max(0.05, scale_kg * 0.5)
        return ReactorSpec(
            reactor_type=rt,
            volume_L=vol,
            temperature_C=mean_t,
            pressure_atm=1.0 if mean_t < 100 else 2.0,
            residence_time_min=2.0,
            mixing_type="static",
            heat_transfer="coil",
            material="stainless steel",
            estimated_cost_usd=_estimate_reactor_cost(rt, vol),
            notes=(
                "Microreactor recommended for fast exothermic reaction at "
                "sub-kilogram scale.  Excellent heat removal and mixing."
            ),
        )

    # --- 2. Fast exothermic, large scale -> CSTR ---
    if is_exothermic and scale_kg >= 1.0:
        rt = ReactorType.CSTR
        return ReactorSpec(
            reactor_type=rt,
            volume_L=volume,
            temperature_C=mean_t,
            pressure_atm=1.0 if mean_t < 100 else 3.0,
            residence_time_min=30.0,
            mixing_type="mechanical",
            heat_transfer="jacketed",
            material=material,
            estimated_cost_usd=_estimate_reactor_cost(rt, volume),
            notes=(
                "CSTR selected for exothermic reaction at production scale.  "
                "Jacket cooling essential; consider cascade of 2-3 CSTRs for "
                "improved conversion."
            ),
        )

    # --- 3. Multiphase / catalytic ---
    if is_multiphase or has_catalyst:
        if scale_kg > 100.0 and has_catalyst:
            rt = ReactorType.FIXED_BED
            return ReactorSpec(
                reactor_type=rt,
                volume_L=volume,
                temperature_C=mean_t,
                pressure_atm=3.0,
                residence_time_min=15.0,
                mixing_type="none",
                heat_transfer="coil",
                material=material,
                estimated_cost_usd=_estimate_reactor_cost(rt, volume),
                notes=(
                    "Fixed-bed reactor for heterogeneous catalytic process "
                    "at >100 kg scale.  Catalyst lifetime and regeneration "
                    "strategy must be defined."
                ),
            )
        rt = ReactorType.BATCH
        return ReactorSpec(
            reactor_type=rt,
            volume_L=volume,
            temperature_C=mean_t,
            pressure_atm=1.0,
            residence_time_min=120.0,
            mixing_type="mechanical",
            heat_transfer="jacketed",
            material=material,
            estimated_cost_usd=_estimate_reactor_cost(rt, volume),
            notes=(
                "Batch reactor for multiphase or catalytic reaction.  "
                "Ensure adequate agitation for mass transfer."
            ),
        )

    # --- 4. High-volume commodity -> PFR ---
    if scale_kg > 500.0:
        rt = ReactorType.PFR
        return ReactorSpec(
            reactor_type=rt,
            volume_L=volume,
            temperature_C=mean_t,
            pressure_atm=5.0,
            residence_time_min=20.0,
            mixing_type="static",
            heat_transfer="coil",
            material="stainless steel",
            estimated_cost_usd=_estimate_reactor_cost(rt, volume),
            notes=(
                "Plug-flow reactor for high-volume continuous production.  "
                "Back-mixing minimised; good for high conversion targets."
            ),
        )

    # --- 5. Controlled addition needed -> semi-batch ---
    needs_slow_addition = (
        template.category in {ReactionCategory.CARBONYL, ReactionCategory.SUBSTITUTION}
        and scale_kg > 10.0
    )
    if needs_slow_addition:
        rt = ReactorType.SEMI_BATCH
        return ReactorSpec(
            reactor_type=rt,
            volume_L=volume,
            temperature_C=mean_t,
            pressure_atm=1.0,
            residence_time_min=90.0,
            mixing_type="mechanical",
            heat_transfer="jacketed",
            material=material,
            estimated_cost_usd=_estimate_reactor_cost(rt, volume),
            notes=(
                "Semi-batch reactor allows controlled reagent addition to "
                "manage selectivity and heat release."
            ),
        )

    # --- 6. Default -> batch ---
    rt = ReactorType.BATCH
    return ReactorSpec(
        reactor_type=rt,
        volume_L=volume,
        temperature_C=mean_t,
        pressure_atm=1.0,
        residence_time_min=60.0,
        mixing_type="mechanical",
        heat_transfer="jacketed" if volume > 20 else "adiabatic",
        material=material,
        estimated_cost_usd=_estimate_reactor_cost(rt, volume),
        notes="Standard batch reactor; suitable for most laboratory and pilot-scale work.",
    )
