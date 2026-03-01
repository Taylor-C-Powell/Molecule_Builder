"""Aqueous solubility and crystallization risk prediction.

Provides ESOL (Delaney 2004) and GSE (Yalkowsky-Valvani) solubility
estimates, melting point heuristic, crystallization risk, and polymorph
risk assessment.  All calculations use Lipinski properties already
computed by MolBuilder -- no external dependencies.

References
----------
- Delaney, J.S. (2004). ESOL: Estimating Aqueous Solubility Directly
  from Molecular Structure. J. Chem. Inf. Comput. Sci. 44, 1000-1005.
- Yalkowsky, S.H. & Valvani, S.C. (1980). Solubility and partitioning
  I: Solubility of nonelectrolytes in water. J. Pharm. Sci. 69, 912-922.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from molbuilder.molecule.graph import Molecule
from molbuilder.molecule.properties import lipinski_properties


# =====================================================================
#  Result dataclass
# =====================================================================

@dataclass
class SolubilityResult:
    """Aqueous solubility and crystallization prediction for a molecule."""
    log_s_esol: float
    log_s_gse: float
    solubility_mg_ml: float
    solubility_class: str
    estimated_melting_point_c: float
    crystallization_risk: str
    polymorph_risk: str


# =====================================================================
#  Private helpers
# =====================================================================

def _aromatic_proportion(mol: Molecule) -> float:
    """Fraction of heavy atoms that are aromatic (have aromatic bonds)."""
    heavy = [a for a in mol.atoms if a.symbol != "H"]
    if not heavy:
        return 0.0

    aromatic_indices: set[int] = set()
    for bond in mol.bonds:
        # Bond order 4 or 1.5 typically marks aromatic in our graph
        if bond.order == 4:
            aromatic_indices.add(bond.atom_i)
            aromatic_indices.add(bond.atom_j)

    # Filter to heavy atoms only
    heavy_indices = {a.index for a in heavy}
    aromatic_heavy = aromatic_indices & heavy_indices
    return len(aromatic_heavy) / len(heavy) if heavy else 0.0


def _ring_count(mol: Molecule) -> int:
    """Count rings via cyclomatic number on heavy-atom subgraph."""
    heavy = {a.index for a in mol.atoms if a.symbol != "H"}
    if not heavy:
        return 0
    edges = sum(
        1 for b in mol.bonds
        if b.atom_i in heavy and b.atom_j in heavy
    )
    # Union-find for connected components
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

    for b in mol.bonds:
        if b.atom_i in heavy and b.atom_j in heavy:
            union(b.atom_i, b.atom_j)

    components = len({find(i) for i in heavy})
    return max(edges - len(heavy) + components, 0)


def _aromatic_ring_count(mol: Molecule) -> int:
    """Estimate aromatic ring count from aromatic atom proportion."""
    ap = _aromatic_proportion(mol)
    heavy = sum(1 for a in mol.atoms if a.symbol != "H")
    aromatic_atoms = int(ap * heavy)
    # Rough estimate: 6 aromatic atoms per ring, but at least count if any
    if aromatic_atoms >= 5:
        return max(aromatic_atoms // 5, 1)
    return 0


def _estimate_melting_point(mol: Molecule, props: object) -> float:
    """Heuristic melting point estimate (Celsius).

    base 80 + 20 per ring + 10 per HBD + 0.01 * MW
    """
    rings = _ring_count(mol)
    mp = 80.0 + 20.0 * rings + 10.0 * props.hbd + 0.01 * props.molecular_weight
    return round(mp, 1)


def _esol(props: object, ap: float) -> float:
    """ESOL model (Delaney 2004): logS = 0.16 - 0.63*logP - 0.0062*MW + 0.066*RB - 0.74*AP."""
    log_s = (
        0.16
        - 0.63 * props.logp
        - 0.0062 * props.molecular_weight
        + 0.066 * props.rotatable_bonds
        - 0.74 * ap
    )
    return round(log_s, 3)


def _gse(props: object, mp: float) -> float:
    """GSE model (Yalkowsky-Valvani): logS = 0.5 - logP - 0.01*(MP - 25)."""
    log_s = 0.5 - props.logp - 0.01 * (mp - 25.0)
    return round(log_s, 3)


def _solubility_class(log_s: float) -> str:
    """Map logS (mol/L) to a human-readable solubility class."""
    if log_s > 0:
        return "highly soluble"
    if log_s > -2:
        return "soluble"
    if log_s > -4:
        return "moderately soluble"
    if log_s > -6:
        return "poorly soluble"
    return "insoluble"


def _crystallization_risk(mol: Molecule, props: object, mp: float) -> str:
    """Assess crystallization risk from MW, aromatic rings, melting point."""
    score = 0
    if props.molecular_weight > 400:
        score += 1
    if _aromatic_ring_count(mol) >= 3:
        score += 1
    if mp > 200:
        score += 1
    if score >= 2:
        return "high"
    if score == 1:
        return "moderate"
    return "low"


def _polymorph_risk(mol: Molecule, props: object) -> str:
    """Assess polymorph risk from H-bond donor/acceptor count."""
    hb_total = props.hbd + props.hba
    if hb_total > 8:
        return "high"
    if hb_total > 4:
        return "moderate"
    return "low"


# =====================================================================
#  Public API
# =====================================================================

def predict_solubility(mol: Molecule) -> SolubilityResult:
    """Predict aqueous solubility and crystallization properties.

    Returns a SolubilityResult with ESOL and GSE predictions,
    solubility class, melting point estimate, and risk assessments.
    """
    props = lipinski_properties(mol)
    ap = _aromatic_proportion(mol)
    mp = _estimate_melting_point(mol, props)

    log_s_esol = _esol(props, ap)
    log_s_gse = _gse(props, mp)

    # Use average of both models for classification and mg/mL conversion
    avg_log_s = (log_s_esol + log_s_gse) / 2.0

    # Convert logS (log mol/L) to mg/mL:  mg/mL = 10^logS * MW / 1000
    if avg_log_s > -10:
        sol_mg_ml = math.pow(10, avg_log_s) * props.molecular_weight / 1000.0
    else:
        sol_mg_ml = 0.0

    return SolubilityResult(
        log_s_esol=log_s_esol,
        log_s_gse=log_s_gse,
        solubility_mg_ml=round(max(sol_mg_ml, 0.0), 6),
        solubility_class=_solubility_class(avg_log_s),
        estimated_melting_point_c=mp,
        crystallization_risk=_crystallization_risk(mol, props, mp),
        polymorph_risk=_polymorph_risk(mol, props),
    )
