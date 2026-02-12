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


# =====================================================================
#  pKa prediction
# =====================================================================

@dataclass
class pKaPrediction:
    """Predicted pKa for a single ionizable group in a molecule."""
    group_name: str
    atom_index: int
    pka_value: float
    acidic: bool


# Base pKa values for ionizable functional groups
_PKA_BASE: dict[str, float] = {
    "carboxylic acid": 4.76,
    "phenol": 10.0,
    "alcohol": 16.0,
    "thiol": 10.5,
    "primary amine": 10.6,
    "secondary amine": 11.0,
    "amide": 15.0,
    "sulfonamide": 10.0,
}


def _ewg_correction(mol: Molecule, idx: int, group_name: str = "") -> float:
    """Additive correction for electron-withdrawing groups on neighbors.

    Scans the alpha heavy-atom neighbors of atom *idx* and applies
    stabilisation-based pKa shifts:
      * halogen (F, Cl, Br) on an alpha carbon:  -1.0 each
      * carbonyl (=O) on an alpha carbon:         -1.5 each

    For carboxylic acids the C=O that defines the group is excluded so
    that only *additional* EWGs contribute.
    """
    correction = 0.0
    for nbr in _heavy_neighbors(mol, idx):
        # For carboxylic acids, the direct C neighbor bearing the
        # defining C=O is part of the functional group -- skip its
        # carbonyl but still count halogens on that carbon.
        skip_carbonyl_on_nbr = (
            group_name == "carboxylic acid"
            and _element(mol, nbr) == "C"
            and _has_double_bond_to(mol, nbr, "O")
        )
        for sub in _neighbors(mol, nbr):
            if sub == idx:
                continue
            sub_sym = _element(mol, sub)
            if sub_sym in ("F", "Cl", "Br"):
                correction -= 1.0
            elif sub_sym == "O" and _bond_order(mol, nbr, sub) == 2:
                if skip_carbonyl_on_nbr:
                    # Skip the first (defining) carbonyl, then reset
                    skip_carbonyl_on_nbr = False
                    continue
                correction -= 1.5
    return correction


def _classify_oxygen(mol: Molecule, o_idx: int) -> tuple[str, int, bool] | None:
    """Classify an O-H oxygen as carboxylic acid, phenol, or alcohol.

    Returns (group_name, ionizable_atom_index, acidic) or None.
    """
    if _h_count(mol, o_idx) < 1:
        return None

    # Look at the heavy neighbors of this oxygen
    heavy = _heavy_neighbors(mol, o_idx)
    if not heavy:
        # Bare water-like O-H, treat as alcohol
        return ("alcohol", o_idx, True)

    for nbr in heavy:
        sym = _element(mol, nbr)
        if sym == "C":
            # Carboxylic acid: O-H where the C also has a C=O
            if _has_double_bond_to(mol, nbr, "O"):
                return ("carboxylic acid", o_idx, True)
            # Phenol: O-H where C is aromatic
            if _is_aromatic_atom(mol, nbr):
                return ("phenol", o_idx, True)
    # Default: ordinary alcohol
    return ("alcohol", o_idx, True)


def _classify_nitrogen(mol: Molecule, n_idx: int) -> tuple[str, int, bool] | None:
    """Classify an N-H nitrogen as amine, amide, or sulfonamide.

    Amines are classified as conjugate-acid donors (acidic=False means the
    protonated form is the acid, so they are base sites).
    """
    h = _h_count(mol, n_idx)
    if h < 1:
        return None

    heavy = _heavy_neighbors(mol, n_idx)

    # Sulfonamide: N-H bonded to S
    for nbr in heavy:
        if _element(mol, nbr) == "S":
            return ("sulfonamide", n_idx, True)

    # Amide: N-H bonded to C that bears C=O
    for nbr in heavy:
        if _element(mol, nbr) == "C" and _has_double_bond_to(mol, nbr, "O"):
            return ("amide", n_idx, True)

    # Primary or secondary amine
    n_heavy = len(heavy)
    if h >= 2 and n_heavy <= 1:
        return ("primary amine", n_idx, False)
    if h >= 1 and n_heavy == 2:
        return ("secondary amine", n_idx, False)

    # Fallback: treat as primary amine if it has H
    return ("primary amine", n_idx, False)


def _classify_sulfur(mol: Molecule, s_idx: int) -> tuple[str, int, bool] | None:
    """Classify an S-H sulfur as a thiol."""
    if _h_count(mol, s_idx) < 1:
        return None
    return ("thiol", s_idx, True)


def _aromatic_ring_correction(mol: Molecule, idx: int, group: str) -> float:
    """Correction when the ionizable atom is near an aromatic ring.

    For most groups the check is on the direct heavy neighbors.  For
    carboxylic acids the ionizable O-H is two bonds from the ring
    (O-H -> C(=O) -> Ar), so we also check the beta (next-neighbor)
    heavy atoms.
    """
    # Direct neighbors
    for nbr in _heavy_neighbors(mol, idx):
        if _is_aromatic_atom(mol, nbr):
            if group in ("carboxylic acid", "phenol"):
                return -1.0
            if group in ("primary amine", "secondary amine"):
                return -2.0
    # Beta neighbors (one bond further) -- relevant for carboxylic acids
    # and amines attached via a linker carbon
    if group in ("carboxylic acid",):
        for nbr in _heavy_neighbors(mol, idx):
            for nbr2 in _heavy_neighbors(mol, nbr):
                if nbr2 == idx:
                    continue
                if _is_aromatic_atom(mol, nbr2):
                    return -1.0
    return 0.0


def _nitro_on_ring_correction(mol: Molecule, idx: int, group: str) -> float:
    """Extra correction for nitro groups on a ring bearing a phenol."""
    if group != "phenol":
        return 0.0
    # Walk: O -> aromatic C -> other aromatic C -> N(=O)(=O)
    for c_nbr in _heavy_neighbors(mol, idx):
        if not _is_aromatic_atom(mol, c_nbr):
            continue
        for ring_nbr in _heavy_neighbors(mol, c_nbr):
            if ring_nbr == idx:
                continue
            if _is_aromatic_atom(mol, ring_nbr):
                for sub in _heavy_neighbors(mol, ring_nbr):
                    if _element(mol, sub) == "N":
                        o_count = sum(
                            1 for nn in _neighbors(mol, sub)
                            if _element(mol, nn) == "O"
                            and _bond_order(mol, sub, nn) == 2
                        )
                        if o_count >= 2:
                            return -3.0
    return 0.0


def predict_pka(mol: Molecule) -> list[pKaPrediction]:
    """Predict pKa values for all ionizable groups in *mol*.

    Uses a Hammett-style additive approach with base pKa values and
    corrections for electron-withdrawing substituents, aromatic rings,
    and nitro groups.  Results are sorted from most acidic (lowest pKa)
    to least acidic.
    """
    results: list[pKaPrediction] = []
    seen: set[int] = set()

    for idx in range(len(mol.atoms)):
        sym = _element(mol, idx)
        classification = None

        if sym == "O":
            classification = _classify_oxygen(mol, idx)
        elif sym == "N":
            classification = _classify_nitrogen(mol, idx)
        elif sym == "S":
            classification = _classify_sulfur(mol, idx)

        if classification is None:
            continue

        group_name, atom_idx, acidic = classification

        # Avoid double-counting the same atom
        if atom_idx in seen:
            continue
        seen.add(atom_idx)

        base = _PKA_BASE[group_name]
        correction = 0.0
        correction += _ewg_correction(mol, atom_idx, group_name)
        correction += _aromatic_ring_correction(mol, atom_idx, group_name)
        correction += _nitro_on_ring_correction(mol, atom_idx, group_name)

        pka = round(base + correction, 1)

        results.append(pKaPrediction(
            group_name=group_name,
            atom_index=atom_idx,
            pka_value=pka,
            acidic=acidic,
        ))

    return sorted(results, key=lambda p: p.pka_value)
