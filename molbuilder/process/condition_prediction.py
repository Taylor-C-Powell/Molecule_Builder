"""Substrate-aware reaction condition prediction.

Given a SMILES string, predict optimal reaction conditions by:
1. Parsing the substrate and detecting functional groups
2. Analysing steric and electronic properties
3. Matching against the reaction template knowledge base
4. Scoring and ranking candidate templates
5. Computing optimised conditions for top matches

Public API:
    predict_conditions(smiles, reaction_name=None, scale_kg=1.0, max_candidates=5)
        -> ConditionPrediction
"""

from __future__ import annotations

from dataclasses import dataclass, field

from molbuilder.smiles.parser import parse
from molbuilder.molecule.graph import Molecule
from molbuilder.molecule.properties import crippen_logp, heavy_atom_count
from molbuilder.reactions.functional_group_detect import (
    detect_functional_groups,
    FunctionalGroup,
)
from molbuilder.reactions.knowledge_base import (
    lookup_by_functional_group,
    lookup_by_name,
)
from molbuilder.reactions.reaction_types import ReactionCategory, ReactionTemplate
from molbuilder.process.conditions import optimize_conditions, ReactionConditions
from molbuilder.process.solvent_systems import select_solvent


# =====================================================================
#  Return types
# =====================================================================

@dataclass
class SubstrateAnalysis:
    """Analysis of a substrate molecule's properties."""
    detected_functional_groups: list[str]
    molecular_weight: float
    heavy_atom_count: int
    steric_class: str           # "unhindered" | "moderately_hindered" | "hindered"
    electronic_character: str   # "electron_rich" | "neutral" | "electron_poor"
    sensitive_groups: list[str]


@dataclass
class SolventScore:
    """A scored solvent for a particular template-substrate combination."""
    solvent_name: str
    composite_score: float      # 0-100
    green_score: int            # 1-10
    reasons: list[str]


@dataclass
class TemplateMatch:
    """A ranked template match with computed conditions."""
    template_name: str
    named_reaction: str | None
    category: str
    match_score: float          # 0-100
    match_reasons: list[str]
    conditions: ReactionConditions
    recommended_solvents: list[SolventScore]
    adjusted_yield_range: tuple[float, float]
    warnings: list[str]
    data_source: str = "heuristic"  # "heuristic" or "ORD (n=...)"


@dataclass
class ConditionPrediction:
    """Full prediction result with substrate analysis and ranked candidates."""
    smiles: str
    substrate_analysis: SubstrateAnalysis
    candidates: list[TemplateMatch]
    best_match: TemplateMatch | None
    overall_confidence: str     # "high" | "medium" | "low"


# =====================================================================
#  Substrate analysis helpers
# =====================================================================

# Functional groups that are sensitive and could be damaged by harsh conditions
_SENSITIVE_GROUPS: set[str] = {
    "aldehyde", "alkene", "alkyne", "epoxide", "thiol",
    "boronic_acid", "acid_chloride", "anhydride",
    "primary_amine", "secondary_amine",
}

# Groups that donate electron density (EDG)
_EDG_NAMES: set[str] = {
    "alcohol", "ether", "primary_amine", "secondary_amine",
    "tertiary_amine", "thiol", "alkene",
}

# Groups that withdraw electron density (EWG)
_EWG_NAMES: set[str] = {
    "nitrile", "nitro", "ketone", "aldehyde", "carboxylic_acid",
    "ester", "amide", "sulfonamide", "sulfone", "acid_chloride",
}


def _classify_steric_environment(
    mol: Molecule, fgs: list[FunctionalGroup],
) -> str:
    """Classify steric hindrance near functional group centers.

    Counts heavy-atom branching (degree >= 3) near each FG center using
    ``mol.neighbors()``.  Returns one of:
    ``"unhindered"``, ``"moderately_hindered"``, ``"hindered"``.
    """
    if not fgs:
        return "unhindered"

    max_branching = 0
    for fg in fgs:
        center = fg.center
        if center < 0 or center >= len(mol.atoms):
            continue
        # Count heavy (non-H) neighbours of the FG center
        heavy_nbrs = [
            n for n in mol.neighbors(center)
            if mol.atoms[n].symbol != "H"
        ]
        degree = len(heavy_nbrs)
        # Also check alpha neighbours for branching
        alpha_branching = 0
        for n in heavy_nbrs:
            n_heavy = sum(
                1 for nn in mol.neighbors(n)
                if mol.atoms[nn].symbol != "H" and nn != center
            )
            alpha_branching = max(alpha_branching, n_heavy)
        total = degree + alpha_branching
        max_branching = max(max_branching, total)

    if max_branching >= 6:
        return "hindered"
    if max_branching >= 4:
        return "moderately_hindered"
    return "unhindered"


def _classify_electronic_character(
    fg_names: list[str], mol: Molecule,
) -> str:
    """Classify overall electronic character of the substrate.

    Counts EDG vs EWG from detected functional groups and aromatic ring
    presence.  Returns ``"electron_rich"``, ``"neutral"``, or
    ``"electron_poor"``.
    """
    has_aromatic = "aromatic_ring" in fg_names
    edg_count = 0.0
    for name in fg_names:
        if name in _EDG_NAMES:
            # Alkene detections inside aromatic rings are artefacts of the
            # FG detector; do not count them as separate EDGs.
            if name == "alkene" and has_aromatic:
                continue
            edg_count += 1
    ewg_count = sum(1 for name in fg_names if name in _EWG_NAMES)

    # Aromatic rings are mildly electron-donating (pi-system)
    aromatic_count = sum(1 for name in fg_names if name == "aromatic_ring")
    edg_count += aromatic_count * 0.5

    diff = edg_count - ewg_count
    if diff >= 1.5:
        return "electron_rich"
    if diff <= -1.5:
        return "electron_poor"
    return "neutral"


def _identify_sensitive_groups(fg_names: list[str]) -> list[str]:
    """Return names of functional groups that are sensitive to damage."""
    return sorted(set(name for name in fg_names if name in _SENSITIVE_GROUPS))


def _analyze_substrate(
    mol: Molecule, fgs: list[FunctionalGroup],
) -> SubstrateAnalysis:
    """Build a SubstrateAnalysis from a parsed molecule and its FGs."""
    from molbuilder.core.elements import atomic_weight

    fg_names = [fg.name for fg in fgs]
    mw = sum(atomic_weight(a.symbol) for a in mol.atoms)

    return SubstrateAnalysis(
        detected_functional_groups=sorted(set(fg_names)),
        molecular_weight=round(mw, 2),
        heavy_atom_count=heavy_atom_count(mol),
        steric_class=_classify_steric_environment(mol, fgs),
        electronic_character=_classify_electronic_character(fg_names, mol),
        sensitive_groups=_identify_sensitive_groups(fg_names),
    )


# =====================================================================
#  Template matching & scoring
# =====================================================================

# Categories that benefit from electron-rich substrates
_ELECTROPHILIC_CATEGORIES: set[ReactionCategory] = {
    ReactionCategory.OXIDATION,
    ReactionCategory.ADDITION,
    ReactionCategory.PERICYCLIC,
}

# Categories that benefit from electron-poor substrates
_NUCLEOPHILIC_CATEGORIES: set[ReactionCategory] = {
    ReactionCategory.SUBSTITUTION,
    ReactionCategory.REDUCTION,
}


def _gather_candidate_templates(
    fg_names: list[str],
    reaction_name: str | None,
) -> list[ReactionTemplate]:
    """Collect unique candidate templates by FG lookup and optional name hint."""
    seen_names: set[str] = set()
    candidates: list[ReactionTemplate] = []

    # Lookup by each detected FG
    for fg in set(fg_names):
        for tmpl in lookup_by_functional_group(fg):
            if tmpl.name not in seen_names:
                seen_names.add(tmpl.name)
                candidates.append(tmpl)

    # Lookup by name hint
    if reaction_name:
        for tmpl in lookup_by_name(reaction_name):
            if tmpl.name not in seen_names:
                seen_names.add(tmpl.name)
                candidates.append(tmpl)

    return candidates


def _score_template(
    template: ReactionTemplate,
    fg_names: list[str],
    substrate: SubstrateAnalysis,
    reaction_name: str | None,
) -> tuple[float, list[str]]:
    """Score a template against the substrate (0-100). Returns (score, reasons)."""
    score = 0.0
    reasons: list[str] = []

    fg_name_set = set(fg_names)
    required = [f.lower() for f in template.functional_group_required]

    # --- FG match completeness (max 30) ---
    if required:
        matched = sum(1 for r in required if r in fg_name_set)
        fg_score = 30.0 * (matched / len(required))
        score += fg_score
        if matched == len(required):
            reasons.append("All required functional groups present")
        else:
            reasons.append(
                f"{matched}/{len(required)} required FGs matched"
            )
    else:
        # Templates with no FG requirements get a partial score
        score += 15.0
        reasons.append("Template has no specific FG requirements")

    # --- Filter: if not all required FGs present, cap at low score ---
    if required:
        matched = sum(1 for r in required if r in fg_name_set)
        if matched < len(required):
            # Still score but penalize -- don't outright reject
            score *= 0.5

    # --- Yield midpoint (max 20) ---
    mid_yield = (template.typical_yield[0] + template.typical_yield[1]) / 2.0
    yield_score = 20.0 * (mid_yield / 100.0)
    score += yield_score
    reasons.append(f"Expected yield {template.typical_yield[0]:.0f}-{template.typical_yield[1]:.0f}%")

    # --- Steric fit (max 10) ---
    steric_map = {"unhindered": 10.0, "moderately_hindered": 6.0, "hindered": 2.0}
    steric_score = steric_map.get(substrate.steric_class, 6.0)
    # SN2 reactions are extra-penalized for hindered substrates
    if template.category == ReactionCategory.SUBSTITUTION and "SN2" in template.name:
        if substrate.steric_class == "hindered":
            steric_score = 0.0
        elif substrate.steric_class == "moderately_hindered":
            steric_score = 3.0
    score += steric_score
    if steric_score < 5.0:
        reasons.append(f"Steric hindrance may reduce efficiency ({substrate.steric_class})")
    else:
        reasons.append(f"Steric environment compatible ({substrate.steric_class})")

    # --- Electronic match (max 10) ---
    electronic_score = 5.0  # neutral default
    if template.category in _ELECTROPHILIC_CATEGORIES:
        if substrate.electronic_character == "electron_rich":
            electronic_score = 10.0
            reasons.append("Electron-rich substrate favors this reaction type")
        elif substrate.electronic_character == "electron_poor":
            electronic_score = 2.0
            reasons.append("Electron-poor substrate may slow this reaction")
    elif template.category in _NUCLEOPHILIC_CATEGORIES:
        if substrate.electronic_character == "electron_poor":
            electronic_score = 10.0
            reasons.append("Electron-poor substrate favors nucleophilic attack")
        elif substrate.electronic_character == "electron_rich":
            electronic_score = 3.0
            reasons.append("Electron-rich substrate may reduce electrophilicity")
    else:
        electronic_score = 5.0
        reasons.append("Electronic character neutral for this reaction type")
    score += electronic_score

    # --- Reaction name hint bonus (max 15) ---
    if reaction_name:
        q = reaction_name.lower()
        name_match = (
            q in template.name.lower()
            or (template.named_reaction and q in template.named_reaction.lower())
            or template.category.name.lower() == q
        )
        if name_match:
            score += 15.0
            reasons.append(f"Matches reaction hint '{reaction_name}'")

    # --- Template specificity bonus (max 15) ---
    if required:
        specificity = min(15.0, len(required) * 5.0)
        score += specificity
        reasons.append(f"Template requires {len(required)} specific FG(s)")
    else:
        score += 3.0

    # --- Compatibility penalty (up to -40) ---
    incompatible = {f.lower() for f in template.functional_group_incompatible}
    clashes = incompatible & fg_name_set
    if clashes:
        penalty = min(40.0, len(clashes) * 20.0)
        score -= penalty
        reasons.append(
            f"Incompatible FGs present: {', '.join(sorted(clashes))}"
        )

    score = max(0.0, min(100.0, score))
    return round(score, 1), reasons


# =====================================================================
#  Solvent scoring with substrate awareness
# =====================================================================

def _score_solvents(
    template: ReactionTemplate,
    mol: Molecule,
    substrate: SubstrateAnalysis,
    scale_kg: float,
) -> list[SolventScore]:
    """Extend base solvent scoring with logP-based solubility heuristic."""
    base_recs = select_solvent(template, scale_kg)
    logp = crippen_logp(mol)

    scored: list[SolventScore] = []
    for rec in base_recs[:8]:  # take top 8 from base scorer
        adj_score = rec.score
        reasons = list(rec.reasons)

        # LogP-based solubility heuristic: polar substrates dissolve better
        # in polar solvents; nonpolar substrates in nonpolar solvents
        if logp > 3.0 and rec.green_score <= 3:
            # Very nonpolar substrate in a green (often polar) solvent
            adj_score -= 5.0
            reasons.append("Substrate may have limited solubility (high logP)")
        elif logp < 0.0 and rec.green_score >= 7:
            adj_score -= 5.0
            reasons.append("Polar substrate in nonpolar solvent")

        # Check FG-solvent compatibility for sensitive groups
        solvent_lower = rec.solvent_name.lower()
        if "thiol" in substrate.sensitive_groups and "dcm" in solvent_lower:
            adj_score -= 3.0
            reasons.append("Thiols may react with chlorinated solvents")
        if "primary_amine" in substrate.sensitive_groups and solvent_lower in ("dmf", "dmso"):
            # Amines are generally fine in DMF/DMSO, no penalty
            pass

        adj_score = max(0.0, min(100.0, adj_score))
        scored.append(SolventScore(
            solvent_name=rec.solvent_name,
            composite_score=round(adj_score, 1),
            green_score=rec.green_score,
            reasons=reasons,
        ))

    scored.sort(key=lambda s: s.composite_score, reverse=True)
    return scored[:5]


# =====================================================================
#  Yield adjustment
# =====================================================================

def _adjust_yield_range(
    template: ReactionTemplate,
    substrate: SubstrateAnalysis,
) -> tuple[float, float]:
    """Adjust expected yield range based on substrate properties."""
    lo, hi = template.typical_yield

    # Steric penalty
    if substrate.steric_class == "hindered":
        lo *= 0.75
        hi *= 0.85
    elif substrate.steric_class == "moderately_hindered":
        lo *= 0.90
        hi *= 0.95

    # Electronic mismatch penalty
    if template.category in _ELECTROPHILIC_CATEGORIES:
        if substrate.electronic_character == "electron_poor":
            lo *= 0.85
            hi *= 0.90
    elif template.category in _NUCLEOPHILIC_CATEGORIES:
        if substrate.electronic_character == "electron_rich":
            lo *= 0.90
            hi *= 0.95

    return (round(max(0.0, lo), 1), round(min(100.0, hi), 1))


# =====================================================================
#  Build warnings
# =====================================================================

def _build_warnings(
    template: ReactionTemplate,
    substrate: SubstrateAnalysis,
) -> list[str]:
    """Generate warnings for a template-substrate combination."""
    warnings: list[str] = []

    # Warn about sensitive groups that might be affected
    incompatible = {f.lower() for f in template.functional_group_incompatible}
    for sg in substrate.sensitive_groups:
        if sg in incompatible:
            warnings.append(
                f"Sensitive group '{sg}' is incompatible with this reaction"
            )
        elif sg not in [f.lower() for f in template.functional_group_required]:
            warnings.append(
                f"Sensitive group '{sg}' present -- verify it survives reaction conditions"
            )

    if substrate.steric_class == "hindered":
        warnings.append("Hindered substrate may require extended reaction time or forcing conditions")

    if substrate.molecular_weight > 500:
        warnings.append("High MW substrate -- solubility and mixing may be challenging")

    return warnings


# =====================================================================
#  Public API
# =====================================================================

def predict_conditions(
    smiles: str,
    reaction_name: str | None = None,
    scale_kg: float = 1.0,
    max_candidates: int = 5,
) -> ConditionPrediction:
    """Predict optimal reaction conditions for a substrate.

    Parameters
    ----------
    smiles : str
        Substrate SMILES string.
    reaction_name : str | None
        Optional hint to prefer a reaction type (e.g. ``"oxidation"``).
    scale_kg : float
        Target production scale in kilograms.
    max_candidates : int
        Maximum number of ranked template matches to return.

    Returns
    -------
    ConditionPrediction
        Full prediction with substrate analysis, ranked candidates, and
        overall confidence assessment.

    Raises
    ------
    ValueError
        If the SMILES string cannot be parsed.
    """
    # Check ML predictor first (returns None if no model loaded)
    from molbuilder.process.ml_predict import get_predictor
    ml_conditions = get_predictor().predict(smiles, reaction_name, scale_kg)

    mol = parse(smiles)
    fgs = detect_functional_groups(mol)
    substrate = _analyze_substrate(mol, fgs)

    fg_names = [fg.name for fg in fgs]
    candidates = _gather_candidate_templates(fg_names, reaction_name)

    # Score and rank
    scored: list[tuple[float, list[str], ReactionTemplate]] = []
    for tmpl in candidates:
        sc, reasons = _score_template(tmpl, fg_names, substrate, reaction_name)
        scored.append((sc, reasons, tmpl))

    scored.sort(key=lambda x: x[0], reverse=True)
    top = scored[:max_candidates]

    # Build TemplateMatch for each top candidate
    matches: list[TemplateMatch] = []
    for idx, (sc, reasons, tmpl) in enumerate(top):
        # Use ML-predicted conditions for best match if available
        if idx == 0 and ml_conditions is not None:
            conditions = ml_conditions
            data_source = "ML model"
        else:
            conditions = optimize_conditions(tmpl, scale_kg)
            data_source = conditions.data_source
        solvents = _score_solvents(tmpl, mol, substrate, scale_kg)
        adj_yield = _adjust_yield_range(tmpl, substrate)
        warnings = _build_warnings(tmpl, substrate)

        matches.append(TemplateMatch(
            template_name=tmpl.name,
            named_reaction=tmpl.named_reaction,
            category=tmpl.category.name,
            match_score=sc,
            match_reasons=reasons,
            conditions=conditions,
            recommended_solvents=solvents,
            adjusted_yield_range=adj_yield,
            warnings=warnings,
            data_source=data_source,
        ))

    best = matches[0] if matches else None

    # Overall confidence
    if best and best.match_score >= 70:
        confidence = "high"
    elif best and best.match_score >= 40:
        confidence = "medium"
    else:
        confidence = "low"

    return ConditionPrediction(
        smiles=smiles,
        substrate_analysis=substrate,
        candidates=matches,
        best_match=best,
        overall_confidence=confidence,
    )
