"""Reaction condition optimisation with scale-dependent adjustments.

Provides :func:`optimize_conditions` which returns a fully populated
:class:`ReactionConditions` dataclass tuned for the target production scale.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from molbuilder.reactions.reaction_types import ReactionCategory, ReactionTemplate
from molbuilder.reactions.reagent_data import normalize_reagent_name


# =====================================================================
#  Data class
# =====================================================================

@dataclass
class ReactionConditions:
    """Optimised reaction conditions for a given template and scale."""

    temperature_C: float
    pressure_atm: float
    solvent: str
    concentration_M: float
    addition_rate: str          # "all at once", "dropwise over X min", etc.
    reaction_time_hours: float
    atmosphere: str             # "air", "N2", "Ar"
    workup_procedure: str
    notes: str


# =====================================================================
#  Internal helpers
# =====================================================================

# Categories that typically need inert atmosphere
_INERT_CATEGORIES = {
    ReactionCategory.COUPLING,
    ReactionCategory.REDUCTION,
    ReactionCategory.RADICAL,
}

# Categories requiring slow / controlled addition
_CONTROLLED_ADDITION_CATEGORIES = {
    ReactionCategory.ADDITION,
    ReactionCategory.CARBONYL,
    ReactionCategory.POLYMERIZATION,
}

# Default concentration (mol/L) by category
_DEFAULT_CONCENTRATIONS: dict[ReactionCategory, float] = {
    ReactionCategory.SUBSTITUTION: 0.5,
    ReactionCategory.ELIMINATION: 0.3,
    ReactionCategory.ADDITION: 0.5,
    ReactionCategory.OXIDATION: 0.2,
    ReactionCategory.REDUCTION: 0.3,
    ReactionCategory.COUPLING: 0.1,
    ReactionCategory.CARBONYL: 0.5,
    ReactionCategory.PROTECTION: 0.5,
    ReactionCategory.DEPROTECTION: 0.3,
    ReactionCategory.REARRANGEMENT: 0.2,
    ReactionCategory.RADICAL: 0.2,
    ReactionCategory.PERICYCLIC: 1.0,
    ReactionCategory.POLYMERIZATION: 2.0,
    ReactionCategory.MISC: 0.3,
}


def _select_atmosphere(template: ReactionTemplate) -> str:
    """Choose atmosphere based on reaction category and reagent sensitivity."""
    if template.category in _INERT_CATEGORIES:
        return "N2"
    # Check for air-sensitive reagents
    sensitive_keywords = {"lialh4", "nabh4", "n_buli", "dibal", "grignard",
                          "memgbr", "etmgbr", "phmgbr", "lda", "lhmds",
                          "nahmds", "khmds", "nah", "red_al"}
    reagent_keys = {normalize_reagent_name(r) for r in template.reagents}
    if reagent_keys & sensitive_keywords:
        return "Ar"
    return "air"


def _addition_rate(template: ReactionTemplate, scale_kg: float) -> str:
    """Determine reagent addition rate.

    Larger scales require slower, controlled addition for thermal management.
    """
    needs_control = (
        template.category in _CONTROLLED_ADDITION_CATEGORIES
        or scale_kg > 10.0
    )
    if not needs_control:
        if scale_kg < 0.1:
            return "all at once"
        return "portion-wise over 10 min"

    # Scale-dependent drip rate
    if scale_kg < 1.0:
        return "dropwise over 15-30 min"
    if scale_kg < 10.0:
        return "dropwise over 30-60 min via addition funnel"
    if scale_kg < 100.0:
        return "metered addition over 1-2 h via peristaltic pump"
    return "metered addition over 2-4 h via mass-flow-controlled pump"


def _estimate_reaction_time(template: ReactionTemplate, scale_kg: float) -> float:
    """Estimate reaction time in hours.

    At larger scale, heat/mass transfer limitations extend reaction time.
    """
    # Base time from yield midpoint heuristic: higher yield takes longer
    _, yield_hi = template.typical_yield
    base_hours = 0.5 + (yield_hi / 100.0) * 2.0  # higher yield -> longer time

    # Temperature effect: cryogenic reactions are faster but need hold time
    mean_t = (template.temperature_range[0] + template.temperature_range[1]) / 2.0
    if mean_t < -40:
        base_hours *= 0.5   # fast at cryo, but add hold time
        base_hours += 0.5
    elif mean_t > 120:
        base_hours *= 0.7   # faster at high temp

    # Scale factor: roughly +20% per decade of scale
    if scale_kg > 1.0:
        import math
        decades = math.log10(scale_kg)
        base_hours *= (1.0 + 0.2 * decades)

    return round(max(0.5, base_hours), 1)


def _workup_procedure(template: ReactionTemplate, scale_kg: float) -> str:
    """Generate a workup procedure string."""
    cat = template.category
    parts: list[str] = []

    if cat in {ReactionCategory.OXIDATION, ReactionCategory.REDUCTION}:
        parts.append("Quench reaction carefully (ice bath).")
        parts.append("Dilute with water or saturated NH4Cl.")
        parts.append("Extract 3x with ethyl acetate or DCM.")
        parts.append("Wash combined organics with brine, dry over Na2SO4.")
    elif cat == ReactionCategory.COUPLING:
        parts.append("Filter through Celite to remove catalyst residues.")
        parts.append("Wash filtrate with water and brine.")
        parts.append("Dry organic layer over MgSO4, concentrate in vacuo.")
    elif cat in {ReactionCategory.SUBSTITUTION, ReactionCategory.ELIMINATION}:
        parts.append("Pour into ice water to quench.")
        parts.append("Extract with DCM or ethyl acetate (3x).")
        parts.append("Wash with saturated NaHCO3, then brine.")
        parts.append("Dry over Na2SO4, filter, concentrate.")
    elif cat == ReactionCategory.CARBONYL:
        parts.append("Quench with saturated NH4Cl solution.")
        parts.append("Extract with ethyl acetate (3x).")
        parts.append("Wash organic layer with brine, dry over MgSO4.")
        parts.append("Filter and concentrate under reduced pressure.")
    elif cat in {ReactionCategory.PROTECTION, ReactionCategory.DEPROTECTION}:
        parts.append("Dilute with ethyl acetate.")
        parts.append("Wash with 1M HCl, saturated NaHCO3, then brine.")
        parts.append("Dry over Na2SO4, concentrate.")
    elif cat == ReactionCategory.POLYMERIZATION:
        parts.append("Precipitate polymer into cold anti-solvent (methanol or hexanes).")
        parts.append("Filter, wash precipitate, dry under vacuum.")
    else:
        parts.append("Standard aqueous workup: dilute, extract, wash, dry, concentrate.")

    if scale_kg > 50:
        parts.append("At this scale, use mixer-settler for extraction.")

    return "  ".join(parts)


def _select_solvent(template: ReactionTemplate) -> str:
    """Pick the first listed solvent from the template, or a sensible default."""
    if template.solvents:
        return template.solvents[0]
    defaults = {
        ReactionCategory.COUPLING: "THF",
        ReactionCategory.REDUCTION: "MeOH",
        ReactionCategory.OXIDATION: "DCM",
        ReactionCategory.CARBONYL: "THF",
        ReactionCategory.SUBSTITUTION: "DMF",
        ReactionCategory.ELIMINATION: "THF",
        ReactionCategory.PROTECTION: "DCM",
        ReactionCategory.DEPROTECTION: "DCM",
        ReactionCategory.POLYMERIZATION: "toluene",
    }
    return defaults.get(template.category, "THF")


# =====================================================================
#  Public API
# =====================================================================

def optimize_conditions(
    template: ReactionTemplate,
    scale_kg: float,
) -> ReactionConditions:
    """Return optimised :class:`ReactionConditions` for *template* at *scale_kg*.

    Scale-dependent adjustments include:
    * Slower addition at larger scale for thermal control
    * Lower concentration at large scale to aid mixing
    * Longer reaction time accounting for heat/mass transfer
    * Appropriate workup procedures
    """
    if not hasattr(template, 'temperature_range'):
        raise TypeError(
            f"template must have a 'temperature_range' attribute, "
            f"got {type(template).__name__}"
        )

    mean_t = (template.temperature_range[0] + template.temperature_range[1]) / 2.0

    # Concentration adjustment: dilute slightly at large scale for mixing
    base_conc = _DEFAULT_CONCENTRATIONS.get(template.category, 0.3)
    if scale_kg > 100.0:
        concentration = base_conc * 0.7
    elif scale_kg > 10.0:
        concentration = base_conc * 0.85
    else:
        concentration = base_conc

    # Pressure: normally atmospheric unless high-temp or hydrogenation
    pressure = 1.0
    if mean_t > 150.0:
        pressure = 3.0
    if template.category == ReactionCategory.REDUCTION:
        # Check for hydrogenation reagents
        h2_reagents = {"h2_pd_c", "pd_c_10", "lindlar"}
        reagent_keys = {normalize_reagent_name(r) for r in template.reagents}
        if reagent_keys & h2_reagents:
            pressure = 4.0  # typical balloon to autoclave

    notes_parts: list[str] = []
    if scale_kg > 100:
        notes_parts.append(
            "Large-scale operation: verify heat removal capacity before start."
        )
    if mean_t < -40:
        notes_parts.append(
            "Cryogenic conditions: use dry ice/acetone or liquid N2 cooling."
        )
    if template.safety_notes:
        notes_parts.append(f"Safety: {template.safety_notes}")

    return ReactionConditions(
        temperature_C=round(mean_t, 1),
        pressure_atm=pressure,
        solvent=_select_solvent(template),
        concentration_M=round(concentration, 2),
        addition_rate=_addition_rate(template, scale_kg),
        reaction_time_hours=_estimate_reaction_time(template, scale_kg),
        atmosphere=_select_atmosphere(template),
        workup_procedure=_workup_procedure(template, scale_kg),
        notes="  ".join(notes_parts) if notes_parts else "No special notes.",
    )
