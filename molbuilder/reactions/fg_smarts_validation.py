"""Cross-validation of heuristic FG detection against SMARTS patterns.

Provides canonical SMARTS patterns for each of the 24 core functional groups,
a cross-validation function that compares heuristic vs. SMARTS detection,
and a validated detection function that unions both methods.

This module does NOT replace the heuristic detector (which is fast and
well-tested).  It adds an optional validation layer for higher confidence.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from molbuilder.molecule.graph import Molecule
from molbuilder.reactions.functional_group_detect import (
    FunctionalGroup,
    detect_functional_groups,
)
from molbuilder.smarts import detect_by_smarts, has_match, parse_smarts


# =====================================================================
#  Canonical SMARTS patterns for 24 core FGs
# =====================================================================

FG_SMARTS_PATTERNS: dict[str, str] = {
    "alcohol":          "[OX2H1]",
    "aldehyde":         "[CX3H1](=O)",
    "ketone":           "[#6][CX3](=O)[#6]",
    "carboxylic_acid":  "[CX3](=O)[OX2H1]",
    "ester":            "[CX3](=O)[OX2][#6]",
    "amide":            "[CX3](=O)[NX3]",
    "primary_amine":    "[NX3H2]",
    "secondary_amine":  "[NX3H1]([#6])[#6]",
    "tertiary_amine":   "[NX3]([#6])([#6])[#6]",
    "alkene":           "[CX3]=[CX3]",
    "alkyne":           "[CX2]#[CX2]",
    "ether":            "[OX2]([#6])[#6]",
    "thiol":            "[SX2H1]",
    "nitrile":          "[CX2]#[NX1]",
    "nitro":            "[NX3](=O)=O",
    "epoxide":          "[OX2]1[CX4][CX4]1",
    "acid_chloride":    "[CX3](=O)[Cl]",
    "anhydride":        "[CX3](=O)[OX2][CX3](=O)",
    "sulfoxide":        "[SX3](=O)([#6])[#6]",
    "sulfone":          "[SX4](=O)(=O)([#6])[#6]",
    "imine":            "[CX3]=[NX2]",
    "alkyl_halide":     "[CX4][F,Cl,Br,I]",
    "boronic_acid":     "[BX3](O)(O)",
    "sulfonamide":      "[SX4](=O)(=O)[NX3]",
}


# =====================================================================
#  Data classes
# =====================================================================

@dataclass
class FGValidationResult:
    """Result of cross-validating heuristic vs. SMARTS FG detection."""

    agreed: list[str] = field(default_factory=list)
    heuristic_only: list[str] = field(default_factory=list)
    smarts_only: list[str] = field(default_factory=list)
    confidence: float = 0.0


# =====================================================================
#  Core functions
# =====================================================================

def _smarts_detect_all(mol: Molecule) -> set[str]:
    """Run all SMARTS patterns and return the set of matched FG names."""
    found: set[str] = set()
    for fg_name, smarts_str in FG_SMARTS_PATTERNS.items():
        try:
            pattern = parse_smarts(smarts_str)
            if has_match(mol, pattern):
                found.add(fg_name)
        except Exception:
            # Skip patterns that fail to parse or match
            continue
    return found


def _normalise_fg_name(name: str) -> str:
    """Map variant names to their canonical forms for comparison."""
    n = name.lower().strip()
    # Map specific halide variants to the generic
    if n.startswith("alkyl_halide_"):
        return "alkyl_halide"
    # Map amine to primary_amine (the heuristic often reports just "amine")
    if n == "amine":
        return "primary_amine"
    if n == "lactam":
        return "amide"
    if n == "cyclic_imine":
        return "imine"
    return n


def cross_validate_fg(mol: Molecule) -> FGValidationResult:
    """Run both heuristic and SMARTS detection and report agreement.

    Returns an :class:`FGValidationResult` with:
    - ``agreed`` -- FG names detected by both methods
    - ``heuristic_only`` -- detected only by the heuristic
    - ``smarts_only`` -- detected only by SMARTS
    - ``confidence`` -- agreement ratio (0.0 to 1.0)
    """
    heuristic_fgs = detect_functional_groups(mol)
    heuristic_names = {_normalise_fg_name(fg.name) for fg in heuristic_fgs}

    smarts_names = {_normalise_fg_name(n) for n in _smarts_detect_all(mol)}

    agreed = sorted(heuristic_names & smarts_names)
    heuristic_only = sorted(heuristic_names - smarts_names)
    smarts_only = sorted(smarts_names - heuristic_names)

    total = len(heuristic_names | smarts_names)
    confidence = len(agreed) / total if total > 0 else 1.0

    return FGValidationResult(
        agreed=agreed,
        heuristic_only=heuristic_only,
        smarts_only=smarts_only,
        confidence=confidence,
    )


def detect_functional_groups_validated(mol: Molecule) -> list[FunctionalGroup]:
    """Detect functional groups using both heuristic and SMARTS methods.

    Returns the union of both methods.  Each FunctionalGroup from the
    heuristic detector is returned as-is.  Any additional groups found
    only by SMARTS are appended with ``smarts_like="SMARTS-only"``.
    """
    heuristic_fgs = detect_functional_groups(mol)
    heuristic_names = {_normalise_fg_name(fg.name) for fg in heuristic_fgs}

    smarts_names = _smarts_detect_all(mol)

    result = list(heuristic_fgs)

    for fg_name in smarts_names:
        normalised = _normalise_fg_name(fg_name)
        if normalised not in heuristic_names:
            # Found by SMARTS but not by heuristic -- add it
            smarts_str = FG_SMARTS_PATTERNS.get(fg_name, "")
            smarts_groups = detect_by_smarts(mol, smarts_str, fg_name)
            result.extend(smarts_groups)

    return result
