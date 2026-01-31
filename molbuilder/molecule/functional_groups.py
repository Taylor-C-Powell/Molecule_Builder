"""
Functional group builder functions.

Each function attaches a functional group to an existing atom in a
Molecule and returns a dict mapping atom labels to their new indices.
"""

from __future__ import annotations

import math

import numpy as np

from molbuilder.molecule.graph import Molecule, Hybridization
from molbuilder.core.bond_data import bond_length, SP3_ANGLE, SP2_ANGLE, SP_ANGLE
from molbuilder.core.geometry import normalize, available_tetrahedral_dirs, place_atom_zmatrix


# ===================================================================
# Constants
# ===================================================================

# Bond lengths not covered by STANDARD_BOND_LENGTHS (Angstroms)
AMIDE_CN = 1.33       # amide C-N (partial double bond character)
SS_BOND = 2.05        # disulfide S-S
CC_AROMATIC = 1.40    # aromatic C-C (average)
CN_AROMATIC = 1.34    # aromatic C-N (average)
CO_CARBOXYL = 1.25    # carboxylate C-O (resonance average)


def add_hydroxyl(mol: Molecule, attach_to: int,
                 angle_ref: int | None = None,
                 dihedral_ref: int | None = None,
                 dihedral_deg: float = 180.0) -> dict[str, int]:
    """Attach -OH to an atom.  Returns {'O': idx, 'H': idx}."""
    if angle_ref is None:
        nbrs = mol.neighbors(attach_to)
        angle_ref = nbrs[0] if nbrs else 0

    o_idx = mol.add_atom_bonded(
        "O", attach_to, bond_order=1, angle_ref=angle_ref,
        dihedral_ref=dihedral_ref, dihedral_deg=dihedral_deg,
        hybridization=Hybridization.SP3)

    h_idx = mol.add_atom_bonded(
        "H", o_idx, bond_order=1, angle_ref=attach_to,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=180.0,
        rotatable=False)

    return {"O": o_idx, "H": h_idx}


def add_amino(mol: Molecule, attach_to: int,
              angle_ref: int | None = None,
              dihedral_ref: int | None = None,
              dihedral_deg: float = 180.0) -> dict[str, int]:
    """Attach -NH2 to an atom.  Returns {'N': idx, 'H1': idx, 'H2': idx}."""
    if angle_ref is None:
        nbrs = mol.neighbors(attach_to)
        angle_ref = nbrs[0] if nbrs else 0

    n_idx = mol.add_atom_bonded(
        "N", attach_to, bond_order=1, angle_ref=angle_ref,
        dihedral_ref=dihedral_ref, dihedral_deg=dihedral_deg,
        hybridization=Hybridization.SP3)

    h1 = mol.add_atom_bonded(
        "H", n_idx, bond_order=1, angle_ref=attach_to,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=120.0,
        rotatable=False)

    h2 = mol.add_atom_bonded(
        "H", n_idx, bond_order=1, angle_ref=attach_to,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=-120.0,
        rotatable=False)

    return {"N": n_idx, "H1": h1, "H2": h2}


def add_carboxyl(mol: Molecule, attach_to: int,
                 angle_ref: int | None = None,
                 dihedral_ref: int | None = None,
                 dihedral_deg: float = 180.0) -> dict[str, int]:
    """Attach -COOH to an atom.  Returns {'C': idx, 'O1': idx, 'O2': idx, 'H': idx}."""
    if angle_ref is None:
        nbrs = mol.neighbors(attach_to)
        angle_ref = nbrs[0] if nbrs else 0

    c_idx = mol.add_atom_bonded(
        "C", attach_to, bond_order=1, angle_ref=angle_ref,
        dihedral_ref=dihedral_ref, dihedral_deg=dihedral_deg,
        hybridization=Hybridization.SP2)

    # C=O (double bond)
    o1 = mol.add_atom_bonded(
        "O", c_idx, bond_order=2, angle_ref=attach_to,
        bond_angle_deg=SP2_ANGLE, dihedral_deg=0.0,
        rotatable=False)

    # C-OH (single bond)
    o2 = mol.add_atom_bonded(
        "O", c_idx, bond_order=1, angle_ref=attach_to,
        bond_angle_deg=SP2_ANGLE, dihedral_deg=180.0,
        hybridization=Hybridization.SP3)

    h = mol.add_atom_bonded(
        "H", o2, bond_order=1, angle_ref=c_idx,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=0.0,
        rotatable=False)

    return {"C": c_idx, "O1": o1, "O2": o2, "H": h}


def add_carbonyl(mol: Molecule, attach_to: int,
                 angle_ref: int | None = None,
                 dihedral_deg: float = 0.0) -> dict[str, int]:
    """Attach C=O to an atom.  Returns {'O': idx}."""
    if angle_ref is None:
        nbrs = mol.neighbors(attach_to)
        angle_ref = nbrs[0] if nbrs else 0

    o_idx = mol.add_atom_bonded(
        "O", attach_to, bond_order=2, angle_ref=angle_ref,
        bond_angle_deg=SP2_ANGLE, dihedral_deg=dihedral_deg,
        rotatable=False)

    return {"O": o_idx}


def add_amide(mol: Molecule, attach_to: int,
              angle_ref: int | None = None,
              dihedral_ref: int | None = None,
              dihedral_deg: float = 180.0) -> dict[str, int]:
    """Attach -CONH2 to an atom.  Returns dict of created indices."""
    if angle_ref is None:
        nbrs = mol.neighbors(attach_to)
        angle_ref = nbrs[0] if nbrs else 0

    c_idx = mol.add_atom_bonded(
        "C", attach_to, bond_order=1, angle_ref=angle_ref,
        dihedral_ref=dihedral_ref, dihedral_deg=dihedral_deg,
        hybridization=Hybridization.SP2)

    o_idx = mol.add_atom_bonded(
        "O", c_idx, bond_order=2, angle_ref=attach_to,
        bond_angle_deg=SP2_ANGLE, dihedral_deg=0.0,
        rotatable=False)

    n_idx = mol.add_atom_bonded(
        "N", c_idx, bond_order=1, angle_ref=attach_to,
        bond_angle_deg=SP2_ANGLE, dihedral_deg=180.0,
        bond_length=AMIDE_CN,
        hybridization=Hybridization.SP2)

    h1 = mol.add_atom_bonded(
        "H", n_idx, bond_order=1, angle_ref=c_idx,
        bond_angle_deg=SP2_ANGLE, dihedral_deg=0.0,
        rotatable=False)

    h2 = mol.add_atom_bonded(
        "H", n_idx, bond_order=1, angle_ref=c_idx,
        bond_angle_deg=SP2_ANGLE, dihedral_deg=180.0,
        rotatable=False)

    return {"C": c_idx, "O": o_idx, "N": n_idx, "H1": h1, "H2": h2}


def add_thiol(mol: Molecule, attach_to: int,
              angle_ref: int | None = None,
              dihedral_ref: int | None = None,
              dihedral_deg: float = 180.0) -> dict[str, int]:
    """Attach -SH to an atom.  Returns {'S': idx, 'H': idx}."""
    if angle_ref is None:
        nbrs = mol.neighbors(attach_to)
        angle_ref = nbrs[0] if nbrs else 0

    s_idx = mol.add_atom_bonded(
        "S", attach_to, bond_order=1, angle_ref=angle_ref,
        dihedral_ref=dihedral_ref, dihedral_deg=dihedral_deg,
        hybridization=Hybridization.SP3)

    h_idx = mol.add_atom_bonded(
        "H", s_idx, bond_order=1, angle_ref=attach_to,
        bond_angle_deg=96.0,  # H-S-C angle ~96 deg
        dihedral_deg=180.0,
        rotatable=False)

    return {"S": s_idx, "H": h_idx}


def _place_planar_ring(mol: Molecule, attach_to: int,
                       symbols: list[str],
                       bond_length_val: float = CC_AROMATIC,
                       angle_ref: int | None = None,
                       dihedral_deg: float = 90.0) -> list[int]:
    """Place a planar ring of atoms bonded to *attach_to*.

    The ring is built analytically: vertices of a regular polygon with
    the first atom bonded to attach_to.  Ring normal is set perpendicular
    to the attach_to bond direction.

    Returns a list of atom indices for the ring atoms.
    """
    n_ring = len(symbols)
    if angle_ref is None:
        nbrs = mol.neighbors(attach_to)
        angle_ref = nbrs[0] if nbrs else 0

    center_pos = mol.atoms[attach_to].position
    ref_pos = mol.atoms[angle_ref].position

    # Bond direction from attach_to toward angle_ref
    bond_dir = normalize(center_pos - ref_pos)

    # Ring normal: perpendicular to bond direction
    up = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(bond_dir, up)) > 0.9:
        up = np.array([0.0, 1.0, 0.0])
    ring_normal = normalize(np.cross(bond_dir, up))

    # Rotate ring_normal by dihedral_deg around bond_dir
    from molbuilder.core.geometry import rotation_matrix
    R_dih = rotation_matrix(bond_dir, math.radians(dihedral_deg))
    ring_normal = R_dih @ ring_normal

    # In-plane direction perpendicular to ring normal and along bond
    in_plane = normalize(np.cross(ring_normal, bond_dir))

    # Ring radius from bond length
    # For a regular n-gon with edge length L: R = L / (2 * sin(pi/n))
    ring_radius = bond_length_val / (2.0 * math.sin(math.pi / n_ring))

    # Center of ring: along bond_dir from attach_to
    # The first ring atom bonds to attach_to, so the ring center
    # is offset by ring_radius from the first atom position.
    # First ring atom is at distance bond_length from attach_to
    # in the bond_dir direction.
    first_atom_pos = center_pos + bond_dir * bond_length_val
    # Ring center is at first_atom_pos - ring_radius * bond_dir
    # (ring is oriented so first atom is closest to attach_to)
    # Actually, let the first atom be at the "top" of the ring.
    ring_center = first_atom_pos - ring_radius * bond_dir

    ring_indices = []
    for i in range(n_ring):
        angle = 2.0 * math.pi * i / n_ring
        pos = (ring_center
               + ring_radius * math.cos(angle) * bond_dir
               + ring_radius * math.sin(angle) * in_plane)
        hyb = Hybridization.SP2
        idx = mol.add_atom(symbols[i], pos, hyb)
        ring_indices.append(idx)

    # Bond ring atoms to each other
    for i in range(n_ring):
        j = (i + 1) % n_ring
        if i < n_ring - 1:
            mol.add_bond(ring_indices[i], ring_indices[j],
                         order=1, rotatable=False)
        else:
            mol.close_ring(ring_indices[i], ring_indices[j])

    # Bond first ring atom to attach_to
    mol.add_bond(attach_to, ring_indices[0], order=1, rotatable=True)

    return ring_indices


def add_phenyl_ring(mol: Molecule, attach_to: int,
                    angle_ref: int | None = None,
                    dihedral_deg: float = 90.0) -> dict[str, int]:
    """Attach a phenyl ring (-C6H5) to an atom.

    Returns {'ring': [6 carbon indices], 'H': [5 hydrogen indices]}.
    """
    symbols = ["C"] * 6
    ring = _place_planar_ring(mol, attach_to, symbols,
                              bond_length_val=CC_AROMATIC,
                              angle_ref=angle_ref,
                              dihedral_deg=dihedral_deg)

    # Add H to each ring C except the one bonded to attach_to (ring[0])
    h_indices = []
    for i in range(1, 6):
        c_pos = mol.atoms[ring[i]].position
        # H points outward from ring center
        ring_center = np.mean(
            [mol.atoms[ring[j]].position for j in range(6)], axis=0)
        outward = normalize(c_pos - ring_center)
        h_pos = c_pos + outward * bond_length("C", "H", 1)
        h_idx = mol.add_atom("H", h_pos)
        mol.add_bond(ring[i], h_idx, order=1, rotatable=False)
        h_indices.append(h_idx)

    return {"ring": ring, "H": h_indices}


def add_imidazole_ring(mol: Molecule, attach_to: int,
                       angle_ref: int | None = None,
                       dihedral_deg: float = 90.0) -> dict[str, int]:
    """Attach an imidazole ring to an atom (5-membered: C-N=C-NH-C).

    Numbering: CG-ND1=CE1-NE2(H)-CD2, with CG bonded to attach_to.
    Returns dict with atom name -> index.
    """
    # Imidazole: 5-membered ring C3N2
    # CG(0)-ND1(1)=CE1(2)-NE2(3)-CD2(4)-CG(close)
    symbols = ["C", "N", "C", "N", "C"]
    ring = _place_planar_ring(mol, attach_to, symbols,
                              bond_length_val=CN_AROMATIC,
                              angle_ref=angle_ref,
                              dihedral_deg=dihedral_deg)

    ring_center = np.mean(
        [mol.atoms[ring[j]].position for j in range(5)], axis=0)

    # Add H to CE1 (ring[2]) and CD2 (ring[4])
    h_indices = {}
    for label, i in [("HE1", 2), ("HD2", 4)]:
        c_pos = mol.atoms[ring[i]].position
        outward = normalize(c_pos - ring_center)
        h_pos = c_pos + outward * bond_length("C", "H", 1)
        h_idx = mol.add_atom("H", h_pos)
        mol.add_bond(ring[i], h_idx, order=1, rotatable=False)
        h_indices[label] = h_idx

    # Add H to NE2 (ring[3]) -- the NH of imidazole
    n_pos = mol.atoms[ring[3]].position
    outward = normalize(n_pos - ring_center)
    h_pos = n_pos + outward * bond_length("N", "H", 1)
    he2 = mol.add_atom("H", h_pos)
    mol.add_bond(ring[3], he2, order=1, rotatable=False)
    h_indices["HE2"] = he2

    return {
        "CG": ring[0], "ND1": ring[1], "CE1": ring[2],
        "NE2": ring[3], "CD2": ring[4], **h_indices,
    }


def add_indole_ring(mol: Molecule, attach_to: int,
                    angle_ref: int | None = None,
                    dihedral_deg: float = 90.0) -> dict[str, int]:
    """Attach an indole ring system (fused 6+5 ring, as in tryptophan).

    The 5-membered ring is bonded to attach_to via CG.
    Returns dict with atom name -> index.
    """
    # Build the 6-membered ring (benzene part) first, then the 5-membered
    # ring fused to it sharing two atoms (CD2, CE2).
    # Indole numbering:
    #   5-ring: CG-CD1=NE1-CE2-CD2-CG
    #   6-ring: CE2-CZ2-CH2-CZ3-CE3-CD2

    # Place 5-membered ring
    symbols_5 = ["C", "C", "N", "C", "C"]
    ring5 = _place_planar_ring(mol, attach_to, symbols_5,
                               bond_length_val=CN_AROMATIC,
                               angle_ref=angle_ref,
                               dihedral_deg=dihedral_deg)

    # CG=ring5[0], CD1=ring5[1], NE1=ring5[2], CE2=ring5[3], CD2=ring5[4]
    # Now build the 6-membered ring fused at CE2-CD2
    ce2_pos = mol.atoms[ring5[3]].position
    cd2_pos = mol.atoms[ring5[4]].position

    # Ring center of the 5-ring
    ring5_center = np.mean(
        [mol.atoms[ring5[j]].position for j in range(5)], axis=0)

    # The 6-ring extends outward from the CE2-CD2 edge
    midpoint = (ce2_pos + cd2_pos) / 2.0
    outward = normalize(midpoint - ring5_center)

    # Build 4 new atoms for the 6-ring (CZ2, CH2, CZ3, CE3)
    # The edge CE2-CD2 is shared; we need 4 more vertices to close a hexagon.
    edge_vec = normalize(ce2_pos - cd2_pos)
    edge_len = float(np.linalg.norm(ce2_pos - cd2_pos))

    # Use the ring normal from the 5-ring
    ring_normal = normalize(np.cross(
        mol.atoms[ring5[1]].position - mol.atoms[ring5[0]].position,
        mol.atoms[ring5[4]].position - mol.atoms[ring5[0]].position,
    ))

    # For a regular hexagon with edge length = CC_AROMATIC:
    h_offset = CC_AROMATIC * math.sin(math.radians(60))  # height
    half_edge = CC_AROMATIC * math.cos(math.radians(60))  # 0.5 * edge

    # 4 new atoms arranged as a hexagon extension
    cz2_pos = ce2_pos + outward * h_offset + edge_vec * half_edge
    ce3_pos = cd2_pos + outward * h_offset - edge_vec * half_edge
    ch2_pos = cz2_pos + outward * h_offset - edge_vec * half_edge
    cz3_pos = ce3_pos + outward * h_offset + edge_vec * half_edge

    cz2 = mol.add_atom("C", cz2_pos, Hybridization.SP2)
    ch2 = mol.add_atom("C", ch2_pos, Hybridization.SP2)
    cz3 = mol.add_atom("C", cz3_pos, Hybridization.SP2)
    ce3 = mol.add_atom("C", ce3_pos, Hybridization.SP2)

    # Bonds for the 6-ring: CE2-CZ2-CH2-CZ3-CE3-CD2
    mol.add_bond(ring5[3], cz2, order=1, rotatable=False)  # CE2-CZ2
    mol.add_bond(cz2, ch2, order=1, rotatable=False)       # CZ2-CH2
    mol.add_bond(ch2, cz3, order=1, rotatable=False)       # CH2-CZ3
    mol.add_bond(cz3, ce3, order=1, rotatable=False)       # CZ3-CE3
    mol.close_ring(ce3, ring5[4])                           # CE3-CD2

    # Add H to exposed atoms
    all_ring = ring5 + [cz2, ch2, cz3, ce3]
    overall_center = np.mean(
        [mol.atoms[idx].position for idx in all_ring], axis=0)

    h_dict = {}
    # CD1 (ring5[1]), NE1 gets H, CZ2, CH2, CZ3, CE3 get H
    for label, idx in [("HD1", ring5[1]), ("HZ2", cz2),
                       ("HH2", ch2), ("HZ3", cz3), ("HE3", ce3)]:
        a_pos = mol.atoms[idx].position
        out = normalize(a_pos - overall_center)
        h_pos = a_pos + out * bond_length("C", "H", 1)
        h_i = mol.add_atom("H", h_pos)
        mol.add_bond(idx, h_i, order=1, rotatable=False)
        h_dict[label] = h_i

    # NE1 hydrogen
    ne1_pos = mol.atoms[ring5[2]].position
    out = normalize(ne1_pos - overall_center)
    h_pos = ne1_pos + out * bond_length("N", "H", 1)
    he1 = mol.add_atom("H", h_pos)
    mol.add_bond(ring5[2], he1, order=1, rotatable=False)
    h_dict["HE1"] = he1

    return {
        "CG": ring5[0], "CD1": ring5[1], "NE1": ring5[2],
        "CE2": ring5[3], "CD2": ring5[4],
        "CZ2": cz2, "CH2": ch2, "CZ3": cz3, "CE3": ce3,
        **h_dict,
    }


def add_guanidinium(mol: Molecule, attach_to: int,
                    angle_ref: int | None = None,
                    dihedral_ref: int | None = None,
                    dihedral_deg: float = 180.0) -> dict[str, int]:
    """Attach a guanidinium group -C(=NH)(NH2)(NH2) to an atom.

    Returns dict with atom name -> index.
    """
    if angle_ref is None:
        nbrs = mol.neighbors(attach_to)
        angle_ref = nbrs[0] if nbrs else 0

    cz = mol.add_atom_bonded(
        "C", attach_to, bond_order=1, angle_ref=angle_ref,
        dihedral_ref=dihedral_ref, dihedral_deg=dihedral_deg,
        hybridization=Hybridization.SP2)

    # NH1 (=NH)
    nh1 = mol.add_atom_bonded(
        "N", cz, bond_order=2, angle_ref=attach_to,
        bond_angle_deg=SP2_ANGLE, dihedral_deg=0.0,
        bond_length=1.33, hybridization=Hybridization.SP2)
    hh11 = mol.add_atom_bonded(
        "H", nh1, bond_order=1, angle_ref=cz,
        bond_angle_deg=SP2_ANGLE, dihedral_deg=180.0,
        rotatable=False)

    # NH2 (-NH2)
    nh2 = mol.add_atom_bonded(
        "N", cz, bond_order=1, angle_ref=attach_to,
        bond_angle_deg=SP2_ANGLE, dihedral_deg=180.0,
        bond_length=AMIDE_CN, hybridization=Hybridization.SP2)
    hh21 = mol.add_atom_bonded(
        "H", nh2, bond_order=1, angle_ref=cz,
        bond_angle_deg=SP2_ANGLE, dihedral_deg=0.0,
        rotatable=False)
    hh22 = mol.add_atom_bonded(
        "H", nh2, bond_order=1, angle_ref=cz,
        bond_angle_deg=SP2_ANGLE, dihedral_deg=180.0,
        rotatable=False)

    return {
        "CZ": cz, "NH1": nh1, "HH11": hh11,
        "NH2": nh2, "HH21": hh21, "HH22": hh22,
    }
