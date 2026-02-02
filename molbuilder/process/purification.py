"""Purification strategy recommendation for synthesis products.

Selects one or more :class:`PurificationStep` instances based on the reaction
type, product characteristics, and production scale.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import List

from molbuilder.reactions.reaction_types import ReactionCategory, ReactionTemplate


# =====================================================================
#  Enums and data classes
# =====================================================================

class PurificationMethod(Enum):
    DISTILLATION = auto()
    RECRYSTALLIZATION = auto()
    COLUMN_CHROMATOGRAPHY = auto()
    FLASH_CHROMATOGRAPHY = auto()
    EXTRACTION = auto()
    FILTRATION = auto()
    PRECIPITATION = auto()
    SUBLIMATION = auto()


@dataclass
class PurificationStep:
    """A single purification operation with expected performance."""

    method: PurificationMethod
    description: str
    estimated_recovery: float   # percent (0-100)
    estimated_purity: float     # percent (0-100)
    scale_appropriate: bool
    notes: str


# =====================================================================
#  Internal helpers
# =====================================================================

# Categories whose products are typically liquids at room temperature
_LIQUID_PRODUCT_CATEGORIES = {
    ReactionCategory.ELIMINATION,
    ReactionCategory.RADICAL,
}

# Categories whose products are typically solids
_SOLID_PRODUCT_CATEGORIES = {
    ReactionCategory.COUPLING,
    ReactionCategory.PROTECTION,
    ReactionCategory.CARBONYL,
    ReactionCategory.PERICYCLIC,
}

# Categories that often produce complex mixtures requiring chromatography
_COMPLEX_MIXTURE_CATEGORIES = {
    ReactionCategory.REARRANGEMENT,
    ReactionCategory.RADICAL,
    ReactionCategory.MISC,
}


def _chromatography_appropriate(scale_kg: float) -> bool:
    """Chromatography is generally impractical above ~5 kg scale."""
    return scale_kg <= 5.0


def _extraction_step(scale_kg: float) -> PurificationStep:
    return PurificationStep(
        method=PurificationMethod.EXTRACTION,
        description=(
            "Liquid-liquid extraction with aqueous wash (saturated NaHCO3, "
            "then brine) to remove polar impurities and inorganic salts."
        ),
        estimated_recovery=92.0,
        estimated_purity=70.0,
        scale_appropriate=True,
        notes="Use separating funnel (lab) or mixer-settler (plant).",
    )


def _distillation_step(scale_kg: float) -> PurificationStep:
    return PurificationStep(
        method=PurificationMethod.DISTILLATION,
        description=(
            "Simple or fractional distillation under reduced pressure "
            "if product bp is below 200 degC.  Use short-path distillation "
            "for heat-sensitive materials."
        ),
        estimated_recovery=85.0,
        estimated_purity=95.0,
        scale_appropriate=True,
        notes=(
            "Highly scalable.  Ensure delta-bp between product and "
            "impurities is >15 degC for simple distillation."
        ),
    )


def _recrystallization_step(scale_kg: float) -> PurificationStep:
    return PurificationStep(
        method=PurificationMethod.RECRYSTALLIZATION,
        description=(
            "Dissolve crude in minimum hot solvent (e.g. ethanol, ethyl "
            "acetate, or toluene), filter hot, cool slowly to crystallise.  "
            "Collect crystals by vacuum filtration."
        ),
        estimated_recovery=75.0,
        estimated_purity=97.0,
        scale_appropriate=True,
        notes="Solvent screening recommended.  May need 2 crops to maximise yield.",
    )


def _flash_chromatography_step(scale_kg: float) -> PurificationStep:
    return PurificationStep(
        method=PurificationMethod.FLASH_CHROMATOGRAPHY,
        description=(
            "Flash column chromatography on silica gel (40-63 um) with "
            "gradient elution (e.g. hexanes/ethyl acetate)."
        ),
        estimated_recovery=80.0,
        estimated_purity=95.0,
        scale_appropriate=_chromatography_appropriate(scale_kg),
        notes=(
            "Practical up to ~5 kg.  Above that, consider preparative HPLC "
            "or alternative purification strategies."
        ),
    )


def _column_chromatography_step(scale_kg: float) -> PurificationStep:
    return PurificationStep(
        method=PurificationMethod.COLUMN_CHROMATOGRAPHY,
        description=(
            "Gravity column chromatography on silica gel.  Suitable when "
            "flash equipment is unavailable; slower but gentler."
        ),
        estimated_recovery=75.0,
        estimated_purity=93.0,
        scale_appropriate=_chromatography_appropriate(scale_kg),
        notes="Load ratio: ~30:1 silica-to-crude by weight.",
    )


def _filtration_step(scale_kg: float) -> PurificationStep:
    return PurificationStep(
        method=PurificationMethod.FILTRATION,
        description=(
            "Vacuum filtration through Celite or sintered-glass funnel "
            "to remove catalyst residues and insoluble by-products."
        ),
        estimated_recovery=95.0,
        estimated_purity=60.0,
        scale_appropriate=True,
        notes="Often the first purification step for heterogeneous reactions.",
    )


def _precipitation_step(scale_kg: float) -> PurificationStep:
    return PurificationStep(
        method=PurificationMethod.PRECIPITATION,
        description=(
            "Add anti-solvent (water, hexanes, or diethyl ether) to "
            "precipitate product from solution.  Collect by filtration."
        ),
        estimated_recovery=80.0,
        estimated_purity=88.0,
        scale_appropriate=True,
        notes="Works best when product is much less soluble than impurities in the anti-solvent.",
    )


def _sublimation_step(scale_kg: float) -> PurificationStep:
    return PurificationStep(
        method=PurificationMethod.SUBLIMATION,
        description=(
            "Vacuum sublimation at reduced pressure.  Best for low-MW "
            "solids with high vapour pressure (e.g. naphthalene, ferrocene)."
        ),
        estimated_recovery=70.0,
        estimated_purity=99.0,
        scale_appropriate=scale_kg <= 1.0,
        notes="Limited throughput; primarily a lab-scale technique.",
    )


# =====================================================================
#  Public API
# =====================================================================

def recommend_purification(
    template: ReactionTemplate,
    scale_kg: float,
) -> list[PurificationStep]:
    """Return an ordered list of purification steps for *template* at *scale_kg*.

    Strategy
    --------
    * Liquid products -> extraction then distillation
    * Solid products  -> extraction then recrystallization
    * Complex mixtures -> extraction + chromatography (small scale) or
      extraction + precipitation + recrystallization (large scale)
    * Catalytic reactions always start with filtration
    * Large scale avoids chromatography
    """
    if not hasattr(template, 'category'):
        raise TypeError(
            f"template must have a 'category' attribute, "
            f"got {type(template).__name__}"
        )

    steps: list[PurificationStep] = []
    cat = template.category
    has_catalyst = len(template.catalysts) > 0

    # --- Step 0: Filtration for catalytic reactions ---
    if has_catalyst:
        steps.append(_filtration_step(scale_kg))

    # --- Liquid products ---
    if cat in _LIQUID_PRODUCT_CATEGORIES:
        steps.append(_extraction_step(scale_kg))
        steps.append(_distillation_step(scale_kg))
        return steps

    # --- Complex mixtures ---
    if cat in _COMPLEX_MIXTURE_CATEGORIES:
        steps.append(_extraction_step(scale_kg))
        if _chromatography_appropriate(scale_kg):
            steps.append(_flash_chromatography_step(scale_kg))
        else:
            steps.append(_precipitation_step(scale_kg))
            steps.append(_recrystallization_step(scale_kg))
        return steps

    # --- Solid products ---
    if cat in _SOLID_PRODUCT_CATEGORIES:
        steps.append(_extraction_step(scale_kg))
        steps.append(_recrystallization_step(scale_kg))
        return steps

    # --- Oxidation / Reduction: often aqueous workup + extraction ---
    if cat in {ReactionCategory.OXIDATION, ReactionCategory.REDUCTION}:
        steps.append(_extraction_step(scale_kg))
        if _chromatography_appropriate(scale_kg) and scale_kg < 0.5:
            steps.append(_flash_chromatography_step(scale_kg))
        else:
            steps.append(_distillation_step(scale_kg))
        return steps

    # --- Substitution / Addition: general workflow ---
    if cat in {ReactionCategory.SUBSTITUTION, ReactionCategory.ADDITION}:
        steps.append(_extraction_step(scale_kg))
        if _chromatography_appropriate(scale_kg) and scale_kg < 1.0:
            steps.append(_flash_chromatography_step(scale_kg))
        else:
            steps.append(_distillation_step(scale_kg))
        return steps

    # --- Deprotection ---
    if cat == ReactionCategory.DEPROTECTION:
        steps.append(_extraction_step(scale_kg))
        steps.append(_precipitation_step(scale_kg))
        return steps

    # --- Polymerization ---
    if cat == ReactionCategory.POLYMERIZATION:
        steps.append(_precipitation_step(scale_kg))
        steps.append(_filtration_step(scale_kg))
        return steps

    # --- Default fallback ---
    steps.append(_extraction_step(scale_kg))
    if _chromatography_appropriate(scale_kg):
        steps.append(_flash_chromatography_step(scale_kg))
    else:
        steps.append(_recrystallization_step(scale_kg))
    return steps
