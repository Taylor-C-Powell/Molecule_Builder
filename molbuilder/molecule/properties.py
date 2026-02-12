"""Lipinski Rule-of-5 and related molecular property calculations.

Provides atom-additive logP (Wildman-Crippen), topological polar surface
area (Ertl TPSA), hydrogen-bond donor/acceptor counts, rotatable bonds,
and Lipinski Ro5 violation assessment.

All calculations work on ``Molecule`` objects from
``molbuilder.molecule.graph``.  SMILES-parsed molecules have explicit
hydrogen atoms, so N-H / O-H bonds are directly visible.
"""

from __future__ import annotations

from dataclasses import dataclass

from molbuilder.molecule.graph import Molecule, Hybridization
from molbuilder.molecule.property_data import (
    CRIPPEN_LOGP,
    TPSA_CONTRIBUTIONS,
    STANDARD_VALENCE,
)


# =====================================================================
#  Result dataclass
# =====================================================================

@dataclass
class LipinskiProperties:
    """Lipinski Rule-of-5 property set for a molecule."""
    molecular_weight: float
    logp: float
    hbd: int
    hba: int
    rotatable_bonds: int
    tpsa: float
    heavy_atom_count: int
    lipinski_violations: int
    lipinski_pass: bool


# =====================================================================
#  Private helpers
# =====================================================================

def _element(mol: Molecule, idx: int) -> str:
    return mol.atoms[idx].symbol


def _neighbors(mol: Molecule, idx: int) -> list[int]:
    return mol.neighbors(idx)


def _neighbor_symbols(mol: Molecule, idx: int) -> list[str]:
    return [_element(mol, n) for n in _neighbors(mol, idx)]


def _bond_order(mol: Molecule, i: int, j: int) -> int:
    bond = mol.get_bond(i, j)
    return bond.order if bond is not None else 0


def _sum_bond_orders(mol: Molecule, idx: int) -> int:
    total = 0
    for n in _neighbors(mol, idx):
        total += _bond_order(mol, idx, n)
    return total


def _h_count(mol: Molecule, idx: int) -> int:
    """Explicit + implicit H count on atom *idx*."""
    explicit_h = sum(1 for s in _neighbor_symbols(mol, idx) if s == "H")
    if explicit_h > 0:
        return explicit_h
    sym = _element(mol, idx)
    valences = STANDARD_VALENCE.get(sym)
    if valences is None:
        return 0
    bo = _sum_bond_orders(mol, idx)
    for v in valences:
        implicit = v - bo
        if implicit >= 0:
            return implicit
    return 0


def _heavy_neighbors(mol: Molecule, idx: int) -> list[int]:
    """Return indices of non-H neighbors."""
    return [n for n in _neighbors(mol, idx) if _element(mol, n) != "H"]


def _has_double_bond_to(mol: Molecule, idx: int, elem: str) -> bool:
    for n in _neighbors(mol, idx):
        if _element(mol, n) == elem and _bond_order(mol, idx, n) == 2:
            return True
    return False


def _is_in_ring(mol: Molecule, i: int, j: int) -> bool:
    return mol.is_in_ring(i, j)


def _is_aromatic_atom(mol: Molecule, idx: int) -> bool:
    """Heuristic: SP2-hybridized atom that participates in a ring."""
    atom = mol.atoms[idx]
    if atom.hybridization != Hybridization.SP2:
        return False
    for n in _neighbors(mol, idx):
        if _is_in_ring(mol, idx, n):
            return True
    return False


def _count_double_bonds(mol: Molecule, idx: int) -> int:
    """Count number of double bonds on atom."""
    count = 0
    for n in _neighbors(mol, idx):
        if _bond_order(mol, idx, n) == 2:
            count += 1
    return count


def _count_triple_bonds(mol: Molecule, idx: int) -> int:
    count = 0
    for n in _neighbors(mol, idx):
        if _bond_order(mol, idx, n) == 3:
            count += 1
    return count


# =====================================================================
#  LogP atom classification
# =====================================================================

def _classify_atom_for_logp(mol: Molecule, idx: int) -> str:
    """Return the Crippen logP atom-type key for atom *idx*."""
    sym = _element(mol, idx)
    nbr_syms = _neighbor_symbols(mol, idx)
    heavy_nbrs = [s for s in nbr_syms if s != "H"]
    aromatic = _is_aromatic_atom(mol, idx)

    if sym == "H":
        # H on carbon vs heteroatom
        for n in _neighbors(mol, idx):
            if _element(mol, n) != "H":
                if _element(mol, n) == "C":
                    return "H_on_C"
                return "H_on_heteroatom"
        return "H_on_C"

    if sym == "C":
        hyb = mol.atoms[idx].hybridization
        if hyb == Hybridization.SP:
            return "C_sp"
        if hyb == Hybridization.SP2 or aromatic:
            if aromatic:
                hetero_nbrs = [s for s in heavy_nbrs if s not in ("C", "H")]
                if hetero_nbrs:
                    return "C_sp2_aromatic_het"
                return "C_sp2_aromatic"
            if _has_double_bond_to(mol, idx, "O"):
                return "C_sp2_carbonyl"
            return "C_sp2_alkene"
        # SP3
        hetero_nbrs = [s for s in heavy_nbrs if s not in ("C", "H")]
        if hetero_nbrs:
            return "C_sp3_hetero"
        return "C_sp3_all_C_H"

    if sym == "N":
        if aromatic:
            return "N_aromatic"
        # Check nitrile
        if _count_triple_bonds(mol, idx) > 0:
            return "N_nitrile"
        # Check nitro (N bonded to 2+ O)
        o_count = sum(1 for s in heavy_nbrs if s == "O")
        if o_count >= 2:
            return "N_nitro"
        # Check amide: N bonded to a C that has C=O
        for n in _heavy_neighbors(mol, idx):
            if _element(mol, n) == "C" and _has_double_bond_to(mol, n, "O"):
                return "N_amide"
        h = _h_count(mol, idx)
        if h >= 2:
            return "N_primary_amine"
        if h == 1:
            return "N_secondary_amine"
        return "N_tertiary_amine"

    if sym == "O":
        h = _h_count(mol, idx)
        if h >= 1:
            return "O_alcohol"
        if _count_double_bonds(mol, idx) > 0:
            # Check if attached to N (nitro)
            for n in _neighbors(mol, idx):
                if _element(mol, n) == "N" and _bond_order(mol, idx, n) == 2:
                    return "O_nitro"
            return "O_carbonyl"
        # Ether or ester
        for n in _heavy_neighbors(mol, idx):
            if _element(mol, n) == "C" and _has_double_bond_to(mol, n, "O"):
                return "O_ester"
        return "O_ether"

    if sym == "S":
        if aromatic:
            return "S_aromatic"
        h = _h_count(mol, idx)
        if h >= 1:
            return "S_thiol"
        return "S_thioether"

    if sym in ("F", "Cl", "Br", "I"):
        return sym

    if sym == "P":
        return "P_any"

    return "other"


# =====================================================================
#  TPSA atom classification
# =====================================================================

def _classify_atom_for_tpsa(mol: Molecule, idx: int) -> tuple | None:
    """Return lookup key for TPSA_CONTRIBUTIONS, or None if non-polar."""
    sym = _element(mol, idx)
    if sym not in ("N", "O", "S", "P"):
        return None

    n_heavy = len(_heavy_neighbors(mol, idx))
    n_h = _h_count(mol, idx)
    n_dbl = _count_double_bonds(mol, idx)
    # Count triple bonds as 2 double bonds for TPSA classification
    n_dbl += 2 * _count_triple_bonds(mol, idx)
    aromatic = _is_aromatic_atom(mol, idx)
    charge = mol.atoms[idx].formal_charge

    return (sym, n_heavy, n_h, n_dbl, aromatic, charge)


# =====================================================================
#  Public API
# =====================================================================

def heavy_atom_count(mol: Molecule) -> int:
    """Count non-hydrogen atoms."""
    return sum(1 for a in mol.atoms if a.symbol != "H")


def hydrogen_bond_donors(mol: Molecule) -> int:
    """Count H-bond donors: N-H and O-H bonds.

    Counts the number of N and O atoms that have at least one hydrogen
    attached (explicit or implicit).
    """
    count = 0
    for idx, atom in enumerate(mol.atoms):
        if atom.symbol in ("N", "O"):
            if _h_count(mol, idx) > 0:
                count += 1
    return count


def hydrogen_bond_acceptors(mol: Molecule) -> int:
    """Count H-bond acceptors: N and O atoms."""
    return sum(1 for a in mol.atoms if a.symbol in ("N", "O"))


def rotatable_bond_count(mol: Molecule) -> int:
    """Count rotatable bonds.

    A bond is rotatable if:
    - It is marked rotatable (single bond, order=1)
    - Neither endpoint is a terminal atom (degree 1 considering only heavy atoms)
    - It is not in a ring
    - It is not an amide C-N bond (C(=O)-N)
    """
    count = 0
    for bond in mol.bonds:
        if not bond.rotatable:
            continue
        if bond.order != 1:
            continue

        i, j = bond.atom_i, bond.atom_j
        sym_i = _element(mol, i)
        sym_j = _element(mol, j)

        # Skip bonds involving H
        if sym_i == "H" or sym_j == "H":
            continue

        # Skip terminal heavy atoms (only 1 heavy neighbor)
        heavy_i = len(_heavy_neighbors(mol, i))
        heavy_j = len(_heavy_neighbors(mol, j))
        if heavy_i <= 1 or heavy_j <= 1:
            continue

        # Skip ring bonds
        if _is_in_ring(mol, i, j):
            continue

        # Skip amide C-N bonds
        if sym_i == "C" and sym_j == "N" and _has_double_bond_to(mol, i, "O"):
            continue
        if sym_j == "C" and sym_i == "N" and _has_double_bond_to(mol, j, "O"):
            continue

        count += 1
    return count


def crippen_logp(mol: Molecule) -> float:
    """Calculate Wildman-Crippen atom-additive logP estimate."""
    total = 0.0
    for idx in range(len(mol.atoms)):
        key = _classify_atom_for_logp(mol, idx)
        total += CRIPPEN_LOGP.get(key, 0.0)
    return round(total, 2)


def topological_polar_surface_area(mol: Molecule) -> float:
    """Calculate Ertl topological polar surface area (TPSA) in Angstroms^2."""
    total = 0.0
    for idx in range(len(mol.atoms)):
        key = _classify_atom_for_tpsa(mol, idx)
        if key is None:
            continue
        contribution = TPSA_CONTRIBUTIONS.get(key)
        if contribution is not None:
            total += contribution
        # If no exact match found, try without charge then without aromatic
        elif key[5] != 0:
            fallback = (key[0], key[1], key[2], key[3], key[4], 0)
            total += TPSA_CONTRIBUTIONS.get(fallback, 0.0)
        elif key[4]:
            fallback = (key[0], key[1], key[2], key[3], False, key[5])
            total += TPSA_CONTRIBUTIONS.get(fallback, 0.0)
    return round(total, 2)


def lipinski_properties(mol: Molecule) -> LipinskiProperties:
    """Calculate Lipinski Rule-of-5 properties for a molecule.

    Returns a ``LipinskiProperties`` dataclass with all fields populated.
    A molecule passes Ro5 if it has at most 1 violation of:
    - MW <= 500
    - logP <= 5
    - HBD <= 5
    - HBA <= 10
    """
    from molbuilder.core.elements import atomic_weight

    mw = sum(atomic_weight(a.symbol) for a in mol.atoms)
    mw = round(mw, 4)
    lp = crippen_logp(mol)
    hbd = hydrogen_bond_donors(mol)
    hba = hydrogen_bond_acceptors(mol)
    rot = rotatable_bond_count(mol)
    tpsa = topological_polar_surface_area(mol)
    heavy = heavy_atom_count(mol)

    violations = 0
    if mw > 500:
        violations += 1
    if lp > 5:
        violations += 1
    if hbd > 5:
        violations += 1
    if hba > 10:
        violations += 1

    return LipinskiProperties(
        molecular_weight=mw,
        logp=lp,
        hbd=hbd,
        hba=hba,
        rotatable_bonds=rot,
        tpsa=tpsa,
        heavy_atom_count=heavy,
        lipinski_violations=violations,
        lipinski_pass=(violations <= 1),
    )
