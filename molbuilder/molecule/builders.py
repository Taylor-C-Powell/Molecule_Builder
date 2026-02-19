"""
Molecule builder functions.

Convenience functions that construct common molecules with correct 3D
geometry: ethane, butane, cyclohexane, 2-butene, and a generic chiral
centre.
"""

from __future__ import annotations

import math

import numpy as np

from molbuilder.molecule.graph import Molecule, Hybridization, RingConformation
from molbuilder.core.bond_data import bond_length, SP3_ANGLE, SP2_ANGLE
from molbuilder.core.geometry import normalize, available_tetrahedral_dirs, add_sp3_hydrogens


def build_ethane(dihedral_deg: float = 60.0) -> Molecule:
    """Build ethane (C2H6) with a specified H-C-C-H dihedral.

    60 deg = staggered (default, lowest energy).
    0 deg = eclipsed (highest energy).
    """
    mol = Molecule("ethane")
    CC = bond_length("C", "C", 1)
    CH = bond_length("C", "H", 1)

    c0 = mol.add_atom("C", [0.0, 0.0, 0.0], Hybridization.SP3)
    c1 = mol.add_atom("C", [0.0, 0.0, CC], Hybridization.SP3)
    mol.add_bond(c0, c1)

    # Place 3 H on C0 at tetrahedral positions pointing away from C1
    c0_bond_dir = normalize(mol.atoms[c1].position - mol.atoms[c0].position)
    c0_h_dirs = available_tetrahedral_dirs([c0_bond_dir], 3)
    for d in c0_h_dirs:
        h_pos = mol.atoms[c0].position + CH * d
        h_idx = mol.add_atom("H", h_pos)
        mol.add_bond(c0, h_idx, order=1, rotatable=False)

    # Place 3 H on C1 at tetrahedral positions pointing away from C0,
    # then rotate by dihedral_deg
    c1_bond_dir = normalize(mol.atoms[c0].position - mol.atoms[c1].position)
    c1_h_dirs = available_tetrahedral_dirs([c1_bond_dir], 3)
    for d in c1_h_dirs:
        h_pos = mol.atoms[c1].position + CH * d
        h_idx = mol.add_atom("H", h_pos)
        mol.add_bond(c1, h_idx, order=1, rotatable=False)

    # Set the desired dihedral angle between the first H on each carbon
    h_on_c0 = [i for i in mol.neighbors(c0) if mol.atoms[i].symbol == "H"]
    h_on_c1 = [i for i in mol.neighbors(c1) if mol.atoms[i].symbol == "H"]
    if h_on_c0 and h_on_c1:
        mol.set_dihedral(h_on_c0[0], c0, c1, h_on_c1[0], dihedral_deg)

    if abs(dihedral_deg - 60.0) < 1:
        label = "staggered"
    elif abs(dihedral_deg) < 1 or abs(dihedral_deg - 360) < 1:
        label = "eclipsed"
    else:
        label = f"dihedral={dihedral_deg:.0f}"
    mol.name = f"ethane ({label})"
    return mol


def build_butane(central_dihedral_deg: float = 180.0) -> Molecule:
    """Build n-butane (C4H10) with a specified C-C-C-C dihedral.

    180 = anti (default, lowest energy).
    60  = gauche.
    0   = eclipsed (syn-periplanar).
    120 = eclipsed (anti-clinal).
    """
    mol = Molecule("butane")
    CC = bond_length("C", "C", 1)

    # Carbon backbone
    c0 = mol.add_atom("C", [0.0, 0.0, 0.0], Hybridization.SP3)
    c1 = mol.add_atom_bonded("C", c0, bond_length=CC,
                              hybridization=Hybridization.SP3)
    c2 = mol.add_atom_bonded("C", c1, angle_ref=c0,
                              bond_angle_deg=SP3_ANGLE,
                              dihedral_deg=180.0,
                              hybridization=Hybridization.SP3)
    c3 = mol.add_atom_bonded("C", c2, angle_ref=c1,
                              dihedral_ref=c0,
                              bond_angle_deg=SP3_ANGLE,
                              dihedral_deg=central_dihedral_deg,
                              hybridization=Hybridization.SP3)

    # Hydrogens: C0 gets 3, C1 gets 2, C2 gets 2, C3 gets 3
    add_sp3_hydrogens(mol, c0, 3)
    add_sp3_hydrogens(mol, c1, 2)
    add_sp3_hydrogens(mol, c2, 2)
    add_sp3_hydrogens(mol, c3, 3)

    if abs(central_dihedral_deg - 180) < 1:
        label = "anti"
    elif abs(abs(central_dihedral_deg) - 60) < 1:
        label = "gauche"
    elif abs(central_dihedral_deg) < 1:
        label = "eclipsed (syn)"
    elif abs(central_dihedral_deg - 120) < 1:
        label = "eclipsed (anti-clinal)"
    else:
        label = f"dihedral={central_dihedral_deg:.0f}"
    mol.name = f"butane ({label})"
    return mol


def build_cyclohexane(
        conformation: RingConformation = RingConformation.CHAIR,
) -> Molecule:
    """Build cyclohexane (C6H12) in chair or boat conformation.

    Chair:  alternating +/- puckering from crystallographic geometry.
    Boat:   atoms 0 and 3 displaced to the same side of the ring plane.
    """
    mol = Molecule(f"cyclohexane ({conformation.name.lower()})")
    CC = bond_length("C", "C", 1)
    CH = bond_length("C", "H", 1)

    # Ring radius so that adjacent C-C distance = CC
    # Adjacent carbons are separated by 60 deg in angle and 2*d in z:
    #   CC^2 = r^2 + (2d)^2  =>  r = sqrt(CC^2 - 4*d^2)
    d = 0.253  # puckering amplitude (Angstroms)
    r = math.sqrt(CC ** 2 - (2.0 * d) ** 2)

    if conformation == RingConformation.CHAIR:
        puckering = [d, -d, d, -d, d, -d]
    elif conformation == RingConformation.BOAT:
        puckering = [d, -d, 0, d, -d, 0]
    else:
        raise ValueError(
            f"Unsupported conformation: {conformation}. "
            f"Use CHAIR or BOAT.")

    # Place ring carbons
    for i in range(6):
        angle = math.radians(60.0 * i)
        x = r * math.cos(angle)
        y = r * math.sin(angle)
        z = puckering[i]
        mol.add_atom("C", [x, y, z], Hybridization.SP3)

    # Close ring bonds (non-rotatable)
    for i in range(6):
        mol.add_bond(i, (i + 1) % 6, order=1, rotatable=False)

    # Add axial and equatorial hydrogens
    center = np.mean([mol.atoms[i].position for i in range(6)], axis=0)
    up = np.array([0.0, 0.0, 1.0])

    for c_idx in range(6):
        c_pos = mol.atoms[c_idx].position
        radial = normalize(c_pos - center)

        # Axial direction alternates with puckering
        if puckering[c_idx] > 0:
            ax_dir = up
        else:
            ax_dir = -up

        # Equatorial H: tetrahedral angle (~109.5 deg) from axial direction,
        # projected outward.  Weight chosen so H-C-H ~ 107 deg.
        eq_dir = normalize(radial + ax_dir * (-0.33))

        # Axial hydrogen
        h_ax_pos = c_pos + CH * normalize(ax_dir)
        h_ax = mol.add_atom("H", h_ax_pos)
        mol.add_bond(c_idx, h_ax, order=1, rotatable=False)

        # Equatorial hydrogen
        h_eq_pos = c_pos + CH * normalize(eq_dir)
        h_eq = mol.add_atom("H", h_eq_pos)
        mol.add_bond(c_idx, h_eq, order=1, rotatable=False)

    return mol


def build_2_butene(is_cis: bool = True) -> Molecule:
    """Build 2-butene (CH3-CH=CH-CH3) as cis (Z) or trans (E)."""
    label = "Z/cis" if is_cis else "E/trans"
    mol = Molecule(f"2-butene ({label})")
    CC_s = bond_length("C", "C", 1)
    CC_d = bond_length("C", "C", 2)

    # C0(sp3) - C1(sp2) = C2(sp2) - C3(sp3)
    c0 = mol.add_atom("C", [0.0, 0.0, 0.0], Hybridization.SP3)
    c1 = mol.add_atom_bonded("C", c0, bond_length=CC_s,
                              hybridization=Hybridization.SP2)
    c2 = mol.add_atom_bonded("C", c1, bond_order=2,
                              angle_ref=c0,
                              bond_length=CC_d,
                              bond_angle_deg=SP2_ANGLE,
                              dihedral_deg=180.0,
                              hybridization=Hybridization.SP2)

    # cis: C0 and C3 on same side (dihedral C0-C1-C2-C3 = 0)
    # trans: opposite sides (dihedral = 180)
    dih = 0.0 if is_cis else 180.0
    c3 = mol.add_atom_bonded("C", c2, angle_ref=c1,
                              dihedral_ref=c0,
                              bond_length=CC_s,
                              bond_angle_deg=SP2_ANGLE,
                              dihedral_deg=dih,
                              hybridization=Hybridization.SP3)

    # Vinyl H on C1 and C2
    mol.add_atom_bonded("H", c1, angle_ref=c0,
                         bond_angle_deg=SP2_ANGLE,
                         dihedral_deg=180.0,
                         rotatable=False)
    mol.add_atom_bonded("H", c2, angle_ref=c1,
                         bond_angle_deg=SP2_ANGLE,
                         dihedral_deg=180.0,
                         rotatable=False)

    # Methyl H's
    add_sp3_hydrogens(mol, c0, 3)
    add_sp3_hydrogens(mol, c3, 3)

    return mol


def build_chiral_molecule(
        substituents: list[str] | None = None,
) -> Molecule:
    """Build a tetrahedral carbon with four different substituents.

    Default: CHFClBr (bromochlorofluoromethane).
    """
    if substituents is None:
        substituents = ["H", "F", "Cl", "Br"]

    if len(substituents) != 4:
        raise ValueError("Exactly 4 substituents required")
    if len(set(substituents)) != 4:
        raise ValueError("All 4 substituents must be different for chirality")

    mol = Molecule(f"C({''.join(substituents)}) -- chiral centre")

    c = mol.add_atom("C", [0.0, 0.0, 0.0], Hybridization.SP3)

    tet_dirs = [
        np.array([1.0, 1.0, 1.0]) / math.sqrt(3),
        np.array([1.0, -1.0, -1.0]) / math.sqrt(3),
        np.array([-1.0, 1.0, -1.0]) / math.sqrt(3),
        np.array([-1.0, -1.0, 1.0]) / math.sqrt(3),
    ]

    for i, sym in enumerate(substituents):
        bl = bond_length("C", sym, 1)
        pos = tet_dirs[i] * bl
        idx = mol.add_atom(sym, pos)
        mol.add_bond(c, idx, order=1)

    return mol
