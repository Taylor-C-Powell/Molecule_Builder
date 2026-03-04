"""Feature extraction for ML-based retrosynthetic disconnection scoring.

Extracts a 90-feature fixed-length vector for any disconnection, regardless
of precursor count.  Features cover:
  A. Target molecule descriptors (13)
  B. Precursor aggregates (10)
  C. Relationship features (8)
  D. Template features (19)
  E. Target FG one-hot (36)
  F. Strategic flags (4)

This module is pure Python (no scikit-learn dependency) and can be used for
both offline training and online inference.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

from molbuilder.molecule.graph import Molecule
from molbuilder.core.elements import atomic_weight
from molbuilder.molecule.properties import (
    heavy_atom_count,
    rotatable_bond_count,
    crippen_logp,
    hydrogen_bond_donors,
    hydrogen_bond_acceptors,
    topological_polar_surface_area,
    predict_pka,
)
from molbuilder.molecule.sa_score import sa_score
from molbuilder.reactions.functional_group_detect import (
    detect_functional_groups,
    FunctionalGroup,
)
from molbuilder.reactions.reaction_types import ReactionTemplate, ReactionCategory
from molbuilder.process.ml_features import (
    _FG_NAMES,
    _CATEGORY_NAMES,
    _molecular_weight,
    _count_rings,
)

if TYPE_CHECKING:
    from molbuilder.reactions.retrosynthesis import Precursor


# =====================================================================
#  Helpers
# =====================================================================

def _safe_heavy_atom_count(smiles: str) -> int:
    """Count heavy atoms from SMILES, returning 0 on failure."""
    try:
        from molbuilder.smiles.parser import parse
        mol = parse(smiles)
        return heavy_atom_count(mol)
    except Exception:
        return 0


def _safe_mw(smiles: str) -> float:
    """Compute molecular weight from SMILES, returning 0.0 on failure."""
    try:
        from molbuilder.smiles.parser import parse
        mol = parse(smiles)
        return _molecular_weight(mol)
    except Exception:
        return 0.0


def _safe_sa_score(mol: Molecule) -> float:
    """Compute SA score, returning 5.0 (midpoint) on failure."""
    try:
        result = sa_score(mol)
        return result.sa_score
    except Exception:
        return 5.0


def _is_purchasable(smiles: str) -> bool:
    """Check purchasability without circular import at module level."""
    from molbuilder.reactions.retrosynthesis import is_purchasable
    return is_purchasable(smiles)


def _precursor_cost(precursor: "Precursor") -> float:
    """Get cost_per_kg for a precursor, defaulting to 500.0 for unknown."""
    if precursor.cost_per_kg > 0:
        return precursor.cost_per_kg
    return 500.0


# =====================================================================
#  CC coupling keywords (shared with score_disconnection heuristic)
# =====================================================================

_CC_KEYWORDS = frozenset({
    "coupling", "grignard", "aldol", "wittig", "suzuki",
    "heck", "sonogashira", "stille", "negishi",
    "horner", "claisen condensation", "michael",
    "robinson", "traube", "methylation",
})


# =====================================================================
#  Feature name constants
# =====================================================================

# Block A: Target molecule descriptors (13)
_TARGET_DESCRIPTOR_NAMES = [
    "target_mw", "target_heavy_atoms", "target_num_bonds",
    "target_num_rings", "target_rotatable_bonds", "target_logp",
    "target_hbd", "target_hba", "target_tpsa",
    "target_chiral_centers", "target_ionizable_groups",
    "target_fg_count", "target_sa_score",
]

# Block B: Precursor aggregates (10)
_PRECURSOR_AGGREGATE_NAMES = [
    "precursor_count", "precursor_mean_heavy", "precursor_max_heavy",
    "precursor_total_heavy", "precursor_mean_mw",
    "precursor_purchasable_frac", "precursor_mean_cost",
    "precursor_max_cost", "precursor_mean_fg_count",
    "precursor_any_purchasable",
]

# Block C: Relationship features (8)
_RELATIONSHIP_NAMES = [
    "heavy_atom_ratio", "mass_balance_ratio", "complexity_reduction",
    "sa_score_reduction", "cost_estimate", "best_precursor_purchasable",
    "all_precursors_purchasable", "depth_in_tree",
]

# Block D: Template features (5 numeric + 14 category one-hot = 19)
_TEMPLATE_NUMERIC_NAMES = [
    "template_yield_lo", "template_yield_hi", "template_yield_mid",
    "template_num_reagents", "template_num_solvents",
]
_TEMPLATE_CATEGORY_NAMES = [f"tcat_{c.lower()}" for c in _CATEGORY_NAMES]

# Block E: Target FG one-hot (36)
_TARGET_FG_NAMES = [f"tfg_{n}" for n in _FG_NAMES]

# Block F: Strategic flags (4)
_STRATEGIC_FLAG_NAMES = [
    "is_cc_coupling", "is_named_reaction",
    "has_incompatible_fgs", "num_fg_required",
]

ALL_RETRO_FEATURE_NAMES: list[str] = (
    _TARGET_DESCRIPTOR_NAMES
    + _PRECURSOR_AGGREGATE_NAMES
    + _RELATIONSHIP_NAMES
    + _TEMPLATE_NUMERIC_NAMES
    + _TEMPLATE_CATEGORY_NAMES
    + _TARGET_FG_NAMES
    + _STRATEGIC_FLAG_NAMES
)


# =====================================================================
#  Main extraction function
# =====================================================================

def extract_retro_features(
    target_mol: Molecule,
    target_smiles: str,
    template: ReactionTemplate,
    precursors: list["Precursor"],
    depth: int = 0,
    target_fgs: list[FunctionalGroup] | None = None,
) -> dict[str, float]:
    """Extract 90-feature vector for a retrosynthetic disconnection.

    Parameters
    ----------
    target_mol : Molecule
        The parsed target molecule.
    target_smiles : str
        SMILES string for the target.
    template : ReactionTemplate
        The reaction template applied in reverse.
    precursors : list[Precursor]
        Precursor molecules from this disconnection.
    depth : int
        Depth in the retro tree (root = 0).
    target_fgs : list[FunctionalGroup] | None
        Pre-computed functional groups for the target (avoids recomputation).

    Returns
    -------
    dict[str, float]
        Feature dictionary with exactly 90 entries.
    """
    features: dict[str, float] = {}

    # Pre-compute target FGs if not provided
    if target_fgs is None:
        target_fgs = detect_functional_groups(target_mol)
    fg_names = {fg.name for fg in target_fgs}

    # === Block A: Target molecule descriptors (13) ===
    target_heavy = heavy_atom_count(target_mol)
    target_mw = _molecular_weight(target_mol)
    target_sa = _safe_sa_score(target_mol)

    features["target_mw"] = target_mw
    features["target_heavy_atoms"] = float(target_heavy)
    features["target_num_bonds"] = float(len(target_mol.bonds))
    features["target_num_rings"] = float(_count_rings(target_mol))
    features["target_rotatable_bonds"] = float(rotatable_bond_count(target_mol))
    features["target_logp"] = crippen_logp(target_mol)
    features["target_hbd"] = float(hydrogen_bond_donors(target_mol))
    features["target_hba"] = float(hydrogen_bond_acceptors(target_mol))
    features["target_tpsa"] = topological_polar_surface_area(target_mol)
    features["target_chiral_centers"] = float(
        sum(1 for a in target_mol.atoms if a.chirality is not None)
    )
    features["target_ionizable_groups"] = float(len(predict_pka(target_mol)))
    features["target_fg_count"] = float(len(target_fgs))
    features["target_sa_score"] = target_sa

    # === Block B: Precursor aggregates (10) ===
    n_prec = len(precursors)
    features["precursor_count"] = float(n_prec)

    if n_prec > 0:
        prec_heavies = [_safe_heavy_atom_count(p.smiles) for p in precursors]
        prec_mws = [_safe_mw(p.smiles) for p in precursors]
        prec_purchasable = [_is_purchasable(p.smiles) for p in precursors]
        prec_costs = [_precursor_cost(p) for p in precursors]

        # Count FGs per precursor
        prec_fg_counts: list[int] = []
        for p in precursors:
            try:
                from molbuilder.smiles.parser import parse
                pmol = parse(p.smiles)
                prec_fg_counts.append(len(detect_functional_groups(pmol)))
            except Exception:
                prec_fg_counts.append(0)

        features["precursor_mean_heavy"] = sum(prec_heavies) / n_prec
        features["precursor_max_heavy"] = float(max(prec_heavies))
        features["precursor_total_heavy"] = float(sum(prec_heavies))
        features["precursor_mean_mw"] = sum(prec_mws) / n_prec
        features["precursor_purchasable_frac"] = (
            sum(1 for p in prec_purchasable if p) / n_prec
        )
        features["precursor_mean_cost"] = sum(prec_costs) / n_prec
        features["precursor_max_cost"] = float(max(prec_costs))
        features["precursor_mean_fg_count"] = sum(prec_fg_counts) / n_prec
        features["precursor_any_purchasable"] = (
            1.0 if any(prec_purchasable) else 0.0
        )
    else:
        features["precursor_mean_heavy"] = 0.0
        features["precursor_max_heavy"] = 0.0
        features["precursor_total_heavy"] = 0.0
        features["precursor_mean_mw"] = 0.0
        features["precursor_purchasable_frac"] = 0.0
        features["precursor_mean_cost"] = 0.0
        features["precursor_max_cost"] = 0.0
        features["precursor_mean_fg_count"] = 0.0
        features["precursor_any_purchasable"] = 0.0

    # === Block C: Relationship features (8) ===
    total_prec_heavy = features["precursor_total_heavy"]
    max_prec_heavy = features["precursor_max_heavy"]

    if target_heavy > 0 and n_prec > 0:
        features["heavy_atom_ratio"] = max_prec_heavy / target_heavy
        features["mass_balance_ratio"] = total_prec_heavy / target_heavy
        reduction = (target_heavy - max_prec_heavy) / target_heavy
        features["complexity_reduction"] = max(0.0, reduction)
    else:
        features["heavy_atom_ratio"] = 0.0
        features["mass_balance_ratio"] = 0.0
        features["complexity_reduction"] = 0.0

    # SA score reduction: compute only for the largest precursor
    if n_prec > 0:
        largest_idx = 0
        largest_heavy = 0
        for i, p in enumerate(precursors):
            h = _safe_heavy_atom_count(p.smiles)
            if h > largest_heavy:
                largest_heavy = h
                largest_idx = i
        try:
            from molbuilder.smiles.parser import parse
            largest_mol = parse(precursors[largest_idx].smiles)
            largest_sa = _safe_sa_score(largest_mol)
            features["sa_score_reduction"] = target_sa - largest_sa
        except Exception:
            features["sa_score_reduction"] = 0.0
    else:
        features["sa_score_reduction"] = 0.0

    # Cost estimate (sum of precursor costs per kg)
    if n_prec > 0:
        features["cost_estimate"] = sum(
            _precursor_cost(p) for p in precursors
        )
    else:
        features["cost_estimate"] = 0.0

    # Best/all precursor purchasable flags
    if n_prec > 0:
        # "Best" = largest precursor
        features["best_precursor_purchasable"] = (
            1.0 if _is_purchasable(precursors[largest_idx].smiles) else 0.0
        )
        features["all_precursors_purchasable"] = (
            1.0 if all(_is_purchasable(p.smiles) for p in precursors) else 0.0
        )
    else:
        features["best_precursor_purchasable"] = 0.0
        features["all_precursors_purchasable"] = 0.0

    features["depth_in_tree"] = float(depth)

    # === Block D: Template features (19) ===
    lo, hi = template.typical_yield
    features["template_yield_lo"] = float(lo)
    features["template_yield_hi"] = float(hi)
    features["template_yield_mid"] = (lo + hi) / 2.0
    features["template_num_reagents"] = float(len(template.reagents))
    features["template_num_solvents"] = float(len(template.solvents))

    # Category one-hot (14 categories)
    cat_name = template.category.name
    for cn in _CATEGORY_NAMES:
        features[f"tcat_{cn.lower()}"] = 1.0 if cn == cat_name else 0.0

    # === Block E: Target FG one-hot (36) ===
    for fn in _FG_NAMES:
        features[f"tfg_{fn}"] = 1.0 if fn in fg_names else 0.0

    # === Block F: Strategic flags (4) ===
    name_lower = template.name.lower()
    named_lower = (template.named_reaction or "").lower()
    features["is_cc_coupling"] = (
        1.0 if any(
            kw in name_lower or kw in named_lower for kw in _CC_KEYWORDS
        ) else 0.0
    )
    features["is_named_reaction"] = (
        1.0 if template.named_reaction else 0.0
    )

    # Incompatible FG check
    incompatible = set(template.functional_group_incompatible or [])
    features["has_incompatible_fgs"] = (
        1.0 if incompatible & fg_names else 0.0
    )
    features["num_fg_required"] = float(
        len(template.functional_group_required or [])
    )

    return features
