"""Synthetic Accessibility (SA) scoring for molecules.

Simplified Ertl-style SA score (1 = easy to synthesize, 10 = hard)
using features MolBuilder already computes.  No RDKit required.

Score components
----------------
- Ring complexity   (0-4):  cyclomatic number + fused/bridged penalty
- Stereo penalty    (0-2):  chiral centres and E/Z bonds
- Fragment penalty  (0-2):  rare / complex functional groups
- Size penalty      (0-1):  heavy-atom count beyond 35
- sp3 bonus       (0 to -1): high sp3 fraction reduces score

Final score = 3 (base) + components, clamped to [1, 10].
"""

from __future__ import annotations

from dataclasses import dataclass

from molbuilder.molecule.graph import Molecule, Hybridization


# =====================================================================
#  Result dataclass
# =====================================================================

@dataclass
class SAScoreResult:
    """Breakdown of synthetic accessibility score for a molecule."""
    sa_score: float
    ring_complexity: float
    stereo_penalty: float
    fragment_penalty: float
    size_penalty: float
    sp3_bonus: float
    heavy_atom_count: int
    ring_count: int
    stereo_count: int


# =====================================================================
#  Private helpers
# =====================================================================

def _heavy_atoms(mol: Molecule) -> list[int]:
    """Return indices of non-hydrogen atoms."""
    return [i for i, a in enumerate(mol.atoms) if a.symbol != "H"]


def _ring_complexity(mol: Molecule) -> float:
    """Score ring complexity from cyclomatic number (0-4).

    cyclomatic = |edges| - |vertices| + connected_components  (for heavy graph)
    Fused ring systems (cyclomatic >= 3) get extra penalty.
    """
    heavy = set(_heavy_atoms(mol))
    if not heavy:
        return 0.0

    # Build heavy-atom subgraph
    edges = 0
    for bond in mol.bonds:
        if bond.atom_i in heavy and bond.atom_j in heavy:
            edges += 1

    # Connected components via union-find on heavy atoms
    parent: dict[int, int] = {i: i for i in heavy}

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for bond in mol.bonds:
        if bond.atom_i in heavy and bond.atom_j in heavy:
            union(bond.atom_i, bond.atom_j)

    components = len({find(i) for i in heavy})
    cyclomatic = edges - len(heavy) + components

    if cyclomatic == 0:
        return 0.0
    if cyclomatic == 1:
        return 1.0
    if cyclomatic == 2:
        return 2.0
    # Fused / bridged systems
    return min(1.5 + cyclomatic * 0.7, 4.0)


def _stereo_penalty(mol: Molecule) -> float:
    """Penalty for stereochemistry (0-2).

    Each chiral centre or E/Z constraint adds 0.5, capped at 2.
    """
    count = 0
    for atom in mol.atoms:
        if atom.chirality and atom.chirality not in ("", "NONE", "None"):
            count += 1
    return min(count * 0.5, 2.0)


def _stereo_count(mol: Molecule) -> int:
    """Count stereocentres (chiral atoms)."""
    count = 0
    for atom in mol.atoms:
        if atom.chirality and atom.chirality not in ("", "NONE", "None"):
            count += 1
    return count


def _fragment_penalty(mol: Molecule) -> float:
    """Penalty for complex or rare functional groups (0-2).

    Heuristic: count distinct FG types.  Having many different types
    increases synthetic difficulty.
    """
    try:
        from molbuilder.reactions.functional_group_detect import detect_functional_groups
        fgs = detect_functional_groups(mol)
        unique = {fg.name for fg in fgs}
        # Rare / complex patterns
        complex_fgs = {
            "epoxide", "aziridine", "nitrile_oxide", "diazo",
            "isocyanate", "isothiocyanate", "acyl_azide",
            "acyl_halide", "anhydride", "peroxide",
        }
        complex_count = len(unique & complex_fgs)
        diversity = len(unique)
        score = 0.0
        if diversity > 4:
            score += 0.5 * min(diversity - 4, 3)
        score += complex_count * 0.5
        return min(score, 2.0)
    except Exception:
        return 0.0


def _size_penalty(mol: Molecule) -> float:
    """Penalty for large molecules (0-1).  Kicks in above 35 heavy atoms."""
    n = len(_heavy_atoms(mol))
    if n <= 35:
        return 0.0
    return min((n - 35) * 0.1, 1.0)


def _sp3_fraction(mol: Molecule) -> float:
    """Fraction of heavy atoms with sp3 hybridisation."""
    heavy = _heavy_atoms(mol)
    if not heavy:
        return 0.0
    sp3 = sum(
        1 for i in heavy
        if mol.atoms[i].hybridization == Hybridization.SP3
    )
    return sp3 / len(heavy)


def _sp3_bonus(mol: Molecule) -> float:
    """Bonus (negative contribution) for high sp3 fraction (0 to -1).

    High sp3 fraction correlates with simpler synthesis in many contexts.
    """
    frac = _sp3_fraction(mol)
    if frac >= 0.5:
        return -min((frac - 0.5) * 2.0, 1.0)
    return 0.0


# =====================================================================
#  Public API
# =====================================================================

def sa_score(mol: Molecule) -> SAScoreResult:
    """Compute synthetic accessibility score for a molecule.

    Returns an SAScoreResult with overall score (1-10, lower = easier)
    and component breakdown.
    """
    ring = _ring_complexity(mol)
    stereo = _stereo_penalty(mol)
    fragment = _fragment_penalty(mol)
    size = _size_penalty(mol)
    sp3 = _sp3_bonus(mol)

    raw = 3.0 + ring + stereo + fragment + size + sp3
    clamped = max(1.0, min(10.0, round(raw, 2)))

    heavy = _heavy_atoms(mol)
    # Count rings: cyclomatic number
    heavy_set = set(heavy)
    edges = sum(
        1 for b in mol.bonds
        if b.atom_i in heavy_set and b.atom_j in heavy_set
    )
    parent: dict[int, int] = {i: i for i in heavy_set}

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for b in mol.bonds:
        if b.atom_i in heavy_set and b.atom_j in heavy_set:
            union(b.atom_i, b.atom_j)

    components = len({find(i) for i in heavy_set}) if heavy_set else 1
    ring_count = edges - len(heavy_set) + components

    return SAScoreResult(
        sa_score=clamped,
        ring_complexity=round(ring, 2),
        stereo_penalty=round(stereo, 2),
        fragment_penalty=round(fragment, 2),
        size_penalty=round(size, 2),
        sp3_bonus=round(sp3, 2),
        heavy_atom_count=len(heavy),
        ring_count=max(ring_count, 0),
        stereo_count=_stereo_count(mol),
    )
