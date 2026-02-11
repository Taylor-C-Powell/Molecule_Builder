"""Solvent selection and scoring for reaction optimisation.

Scores every solvent in ``SOLVENT_DB`` against the requirements of a
:class:`ReactionTemplate` and returns a ranked list of
:class:`SolventRecommendation` instances.
"""

from __future__ import annotations

from dataclasses import dataclass

from molbuilder.reactions.reaction_types import ReactionCategory, ReactionTemplate
from molbuilder.reactions.reagent_data import SOLVENT_DB, Solvent


# =====================================================================
#  Data class
# =====================================================================

@dataclass
class SolventRecommendation:
    """A scored solvent recommendation."""

    solvent_name: str
    score: float              # 0-100 composite score
    reasons: list[str]
    green_score: int          # 1-10 (1 = greenest)
    estimated_cost_per_L: float


# =====================================================================
#  Internal scoring helpers
# =====================================================================

# Desired polarity index ranges per reaction category
_POLARITY_PREFERENCES: dict[ReactionCategory, tuple[float, float]] = {
    ReactionCategory.SUBSTITUTION: (3.0, 7.0),
    ReactionCategory.ELIMINATION: (2.0, 5.0),
    ReactionCategory.ADDITION: (2.0, 6.0),
    ReactionCategory.OXIDATION: (4.0, 8.0),
    ReactionCategory.REDUCTION: (3.0, 7.0),
    ReactionCategory.COUPLING: (3.0, 7.0),
    ReactionCategory.CARBONYL: (3.0, 7.0),
    ReactionCategory.PROTECTION: (2.0, 5.0),
    ReactionCategory.DEPROTECTION: (3.0, 7.0),
    ReactionCategory.REARRANGEMENT: (2.0, 5.0),
    ReactionCategory.RADICAL: (1.0, 4.0),
    ReactionCategory.PERICYCLIC: (2.0, 5.0),
    ReactionCategory.POLYMERIZATION: (1.0, 5.0),
    ReactionCategory.MISC: (2.0, 7.0),
}


def _polarity_score(solvent: Solvent, category: ReactionCategory) -> tuple[float, str]:
    """Score 0-30 based on polarity index match."""
    lo, hi = _POLARITY_PREFERENCES.get(category, (2.0, 7.0))
    pi = solvent.polarity_index
    if lo <= pi <= hi:
        return 30.0, "Polarity index within ideal range for this reaction type"
    distance = min(abs(pi - lo), abs(pi - hi))
    score = max(0.0, 30.0 - distance * 6.0)
    reason = "Polarity index outside preferred range" if score < 15 else "Polarity acceptable"
    return round(score, 1), reason


def _bp_score(solvent: Solvent, template: ReactionTemplate) -> tuple[float, str]:
    """Score 0-25 based on boiling-point suitability.

    The solvent bp should be at least 10 degC above the reaction temperature
    (for reflux margin) but not excessively high (hard to remove).
    """
    mean_t = (template.temperature_range[0] + template.temperature_range[1]) / 2.0
    bp = solvent.bp

    if bp < mean_t:
        return 0.0, "Boiling point below reaction temperature"
    margin = bp - mean_t
    if 10.0 <= margin <= 60.0:
        return 25.0, "Ideal bp margin above reaction temperature"
    if margin < 10.0:
        return 10.0, "bp margin tight; risk of solvent loss"
    # margin > 60 -- harder to remove
    penalty = min(15.0, (margin - 60.0) * 0.3)
    return round(25.0 - penalty, 1), "High bp; may be difficult to remove"


def _green_score(solvent: Solvent) -> tuple[float, str]:
    """Score 0-20 based on green chemistry rating (lower = greener)."""
    gs = solvent.green_score
    score = max(0.0, 20.0 - (gs - 1) * 2.2)
    if gs <= 3:
        reason = "Excellent green chemistry profile"
    elif gs <= 5:
        reason = "Acceptable green chemistry profile"
    elif gs <= 7:
        reason = "Moderate environmental concern"
    else:
        reason = "Poor green chemistry profile; consider alternatives"
    return round(score, 1), reason


def _cost_score(solvent: Solvent, scale_kg: float) -> tuple[float, str]:
    """Score 0-15 based on cost (more important at larger scale)."""
    cost = solvent.cost_per_L
    if cost <= 0:
        return 15.0, "Negligible solvent cost"
    # At large scale, cost matters more
    threshold = 20.0 if scale_kg < 10 else 12.0
    if cost <= threshold:
        return 15.0, "Solvent cost acceptable for this scale"
    excess = cost - threshold
    score = max(0.0, 15.0 - excess * 0.5)
    return round(score, 1), "Solvent cost may be significant at scale"


def _separation_score(solvent: Solvent) -> tuple[float, str]:
    """Score 0-10 based on ease of separation (low bp, immiscible with water)."""
    score = 5.0
    reasons = []
    if not solvent.miscible_with_water:
        score += 3.0
        reasons.append("immiscible with water (easy extraction)")
    if solvent.bp < 80.0:
        score += 2.0
        reasons.append("low bp (easy rotovap removal)")
    elif solvent.bp > 150.0:
        score -= 2.0
        reasons.append("high bp (difficult to remove)")
    reason = "; ".join(reasons) if reasons else "Average separation characteristics"
    return round(min(score, 10.0), 1), reason


def _listed_solvent_bonus(solvent: Solvent, template: ReactionTemplate) -> float:
    """Return 10 bonus points if this solvent appears in the template's solvent list."""
    template_names = {s.lower() for s in template.solvents}
    if solvent.name.lower() in template_names:
        return 10.0
    # Also check common abbreviation-style keys
    for key, db_solvent in SOLVENT_DB.items():
        if db_solvent is solvent and key.lower() in template_names:
            return 10.0
    return 0.0


# =====================================================================
#  Public API
# =====================================================================

def select_solvent(
    template: ReactionTemplate,
    scale_kg: float = 1.0,
) -> list[SolventRecommendation]:
    """Score and rank solvents from ``SOLVENT_DB`` for *template*.

    Returns a list of :class:`SolventRecommendation` sorted best-first.
    Only solvents scoring above 30 are included.
    """
    recommendations: list[SolventRecommendation] = []

    for _key, solvent in SOLVENT_DB.items():
        total = 0.0
        reasons: list[str] = []

        ps, pr = _polarity_score(solvent, template.category)
        total += ps
        reasons.append(pr)

        bs, br = _bp_score(solvent, template)
        total += bs
        reasons.append(br)

        gs, gr = _green_score(solvent)
        total += gs
        reasons.append(gr)

        cs, cr = _cost_score(solvent, scale_kg)
        total += cs
        reasons.append(cr)

        ss, sr = _separation_score(solvent)
        total += ss
        reasons.append(sr)

        bonus = _listed_solvent_bonus(solvent, template)
        if bonus:
            total += bonus
            reasons.append("Listed as suitable solvent in reaction template")

        total = min(total, 100.0)

        if total >= 30.0:
            recommendations.append(
                SolventRecommendation(
                    solvent_name=solvent.name,
                    score=round(total, 1),
                    reasons=reasons,
                    green_score=solvent.green_score,
                    estimated_cost_per_L=solvent.cost_per_L,
                )
            )

    recommendations.sort(key=lambda r: r.score, reverse=True)
    return recommendations
