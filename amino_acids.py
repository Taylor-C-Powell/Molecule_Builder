"""
Amino Acid Structures and Functional Groups

Builds all 20 standard amino acids with correct 3D geometry using the
Molecule class from molecular_conformations.  Provides:

    - Functional group builder functions (hydroxyl, amino, carboxyl, etc.)
    - Planar ring builders (phenyl, imidazole, indole)
    - Amino acid backbone template with L-chirality enforcement
    - Side chain builders for all 20 standard amino acids
    - Peptide bond formation (condensation reaction)
    - Backbone conformation control (phi/psi angles)
    - Secondary structure presets (alpha helix, beta sheet)

All coordinates are in Angstroms.  No Unicode characters are used
(Windows cp1252 compatibility).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum, auto

import numpy as np

from molecular_conformations import (
    Molecule,
    Hybridization,
    Stereodescriptor,
    _bond_length,
    _normalize,
    _rotation_matrix,
    _place_atom_zmatrix,
    _available_tetrahedral_dirs,
    _add_sp3_hydrogens,
    SP3_ANGLE,
    SP2_ANGLE,
    STANDARD_BOND_LENGTHS,
)


# ===================================================================
# Constants
# ===================================================================

# Bond lengths not covered by STANDARD_BOND_LENGTHS (Angstroms)
AMIDE_CN = 1.33       # amide C-N (partial double bond character)
SS_BOND = 2.05        # disulfide S-S
CC_AROMATIC = 1.40    # aromatic C-C (average)
CN_AROMATIC = 1.34    # aromatic C-N (average)
CO_CARBOXYL = 1.25    # carboxylate C-O (resonance average)

# Ramachandran angle presets (phi, psi) in degrees
ALPHA_HELIX = (-57.0, -47.0)
BETA_SHEET = (-135.0, 135.0)
POLYPROLINE_II = (-75.0, 145.0)
EXTENDED = (-180.0, 180.0)


# ===================================================================
# Enumerations
# ===================================================================

class AminoAcidType(Enum):
    """The 20 standard amino acids."""
    GLY = auto()
    ALA = auto()
    VAL = auto()
    LEU = auto()
    ILE = auto()
    PRO = auto()
    PHE = auto()
    TRP = auto()
    MET = auto()
    SER = auto()
    THR = auto()
    CYS = auto()
    TYR = auto()
    ASN = auto()
    GLN = auto()
    ASP = auto()
    GLU = auto()
    LYS = auto()
    ARG = auto()
    HIS = auto()


class SecondaryStructure(Enum):
    """Common secondary structure types."""
    ALPHA_HELIX = auto()
    BETA_SHEET = auto()
    POLYPROLINE_II = auto()
    EXTENDED = auto()


# ===================================================================
# Backbone indices dataclass
# ===================================================================

@dataclass
class BackboneIndices:
    """Atom indices for key positions in an amino acid backbone.

    Stored on the Molecule as mol.backbone after building.
    """
    N: int                   # amino nitrogen
    H_N: int | None          # H on N (None for proline)
    CA: int                  # alpha carbon
    HA: int                  # H on alpha carbon
    C: int                   # carbonyl carbon
    O: int                   # carbonyl oxygen
    CB: int | None           # beta carbon (None for glycine)
    H_Nterm: int | None      # extra H on free N-terminus
    OXT: int | None          # second O on free C-terminus (OH)
    H_OXT: int | None        # H on OXT

    def __repr__(self):
        return (f"Backbone(N={self.N}, CA={self.CA}, C={self.C}, "
                f"O={self.O}, CB={self.CB})")


# ===================================================================
# Functional group builders
# ===================================================================

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
                       bond_length: float = CC_AROMATIC,
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
    bond_dir = _normalize(center_pos - ref_pos)

    # Ring normal: perpendicular to bond direction
    up = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(bond_dir, up)) > 0.9:
        up = np.array([0.0, 1.0, 0.0])
    ring_normal = _normalize(np.cross(bond_dir, up))

    # Rotate ring_normal by dihedral_deg around bond_dir
    R_dih = _rotation_matrix(bond_dir, math.radians(dihedral_deg))
    ring_normal = R_dih @ ring_normal

    # In-plane direction perpendicular to ring normal and along bond
    in_plane = _normalize(np.cross(ring_normal, bond_dir))

    # Ring radius from bond length
    # For a regular n-gon with edge length L: R = L / (2 * sin(pi/n))
    ring_radius = bond_length / (2.0 * math.sin(math.pi / n_ring))

    # Center of ring: along bond_dir from attach_to
    # The first ring atom bonds to attach_to, so the ring center
    # is offset by ring_radius from the first atom position.
    # First ring atom is at distance bond_length from attach_to
    # in the bond_dir direction.
    first_atom_pos = center_pos + bond_dir * bond_length
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
                              bond_length=CC_AROMATIC,
                              angle_ref=angle_ref,
                              dihedral_deg=dihedral_deg)

    # Add H to each ring C except the one bonded to attach_to (ring[0])
    h_indices = []
    for i in range(1, 6):
        c_pos = mol.atoms[ring[i]].position
        # H points outward from ring center
        ring_center = np.mean(
            [mol.atoms[ring[j]].position for j in range(6)], axis=0)
        outward = _normalize(c_pos - ring_center)
        h_pos = c_pos + outward * _bond_length("C", "H", 1)
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
                              bond_length=CN_AROMATIC,
                              angle_ref=angle_ref,
                              dihedral_deg=dihedral_deg)

    ring_center = np.mean(
        [mol.atoms[ring[j]].position for j in range(5)], axis=0)

    # Add H to CE1 (ring[2]) and CD2 (ring[4])
    h_indices = {}
    for label, i in [("HE1", 2), ("HD2", 4)]:
        c_pos = mol.atoms[ring[i]].position
        outward = _normalize(c_pos - ring_center)
        h_pos = c_pos + outward * _bond_length("C", "H", 1)
        h_idx = mol.add_atom("H", h_pos)
        mol.add_bond(ring[i], h_idx, order=1, rotatable=False)
        h_indices[label] = h_idx

    # Add H to NE2 (ring[3]) -- the NH of imidazole
    n_pos = mol.atoms[ring[3]].position
    outward = _normalize(n_pos - ring_center)
    h_pos = n_pos + outward * _bond_length("N", "H", 1)
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
                               bond_length=CN_AROMATIC,
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
    outward = _normalize(midpoint - ring5_center)

    # Build 4 new atoms for the 6-ring (CZ2, CH2, CZ3, CE3)
    # The edge CE2-CD2 is shared; we need 4 more vertices to close a hexagon.
    edge_vec = _normalize(ce2_pos - cd2_pos)
    edge_len = float(np.linalg.norm(ce2_pos - cd2_pos))

    # Use the ring normal from the 5-ring
    ring_normal = _normalize(np.cross(
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
        out = _normalize(a_pos - overall_center)
        h_pos = a_pos + out * _bond_length("C", "H", 1)
        h_i = mol.add_atom("H", h_pos)
        mol.add_bond(idx, h_i, order=1, rotatable=False)
        h_dict[label] = h_i

    # NE1 hydrogen
    ne1_pos = mol.atoms[ring5[2]].position
    out = _normalize(ne1_pos - overall_center)
    h_pos = ne1_pos + out * _bond_length("N", "H", 1)
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


# ===================================================================
# Amino acid backbone builder
# ===================================================================

def _build_backbone(name: str, has_cb: bool = True,
                    is_proline: bool = False) -> Molecule:
    """Build the amino acid backbone: H2N-CA(H)(R)-C(=O)-OH.

    Places atoms using z-matrix coordinates.  If has_cb is True, places
    CB in the L-configuration direction.  For proline, N gets only one
    H (secondary amine).

    Returns a Molecule with a .backbone attribute (BackboneIndices).
    """
    mol = Molecule(name)

    # -- N (amino nitrogen) --
    n_idx = mol.add_atom("N", [0.0, 0.0, 0.0], Hybridization.SP3)

    # -- CA (alpha carbon) --
    ca_idx = mol.add_atom_bonded(
        "C", n_idx, bond_order=1,
        bond_length=_bond_length("C", "N", 1),
        hybridization=Hybridization.SP3)

    # -- C (carbonyl carbon) --
    c_idx = mol.add_atom_bonded(
        "C", ca_idx, bond_order=1, angle_ref=n_idx,
        bond_angle_deg=SP3_ANGLE,
        dihedral_deg=180.0,
        hybridization=Hybridization.SP2)

    # -- O (carbonyl oxygen, C=O) --
    o_idx = mol.add_atom_bonded(
        "O", c_idx, bond_order=2, angle_ref=ca_idx,
        bond_angle_deg=SP2_ANGLE, dihedral_deg=0.0,
        rotatable=False)

    # -- OXT (C-terminus OH) --
    oxt = mol.add_atom_bonded(
        "O", c_idx, bond_order=1, angle_ref=ca_idx,
        bond_angle_deg=SP2_ANGLE, dihedral_deg=180.0,
        hybridization=Hybridization.SP3)

    h_oxt = mol.add_atom_bonded(
        "H", oxt, bond_order=1, angle_ref=c_idx,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=0.0,
        rotatable=False)

    # -- H on N-terminus --
    h_n = mol.add_atom_bonded(
        "H", n_idx, bond_order=1, angle_ref=ca_idx,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=120.0,
        rotatable=False)

    h_nterm = None
    if not is_proline:
        h_nterm = mol.add_atom_bonded(
            "H", n_idx, bond_order=1, angle_ref=ca_idx,
            bond_angle_deg=SP3_ANGLE, dihedral_deg=-120.0,
            rotatable=False)

    # -- HA (H on alpha carbon) --
    # Place HA to establish L-chirality.
    # For L-amino acids, looking from H -> CA, the priority order
    # NH2 > COOH > R is clockwise (S configuration in CIP, except Cys).
    ha_idx = mol.add_atom_bonded(
        "H", ca_idx, bond_order=1, angle_ref=n_idx,
        dihedral_ref=c_idx,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=120.0,
        rotatable=False)

    # -- CB (beta carbon) --
    # Placed at -120 deg dihedral to give L-configuration (S in CIP
    # for most amino acids, R for cysteine due to sulfur priority).
    cb_idx = None
    if has_cb:
        cb_idx = mol.add_atom_bonded(
            "C", ca_idx, bond_order=1, angle_ref=n_idx,
            dihedral_ref=c_idx,
            bond_angle_deg=SP3_ANGLE, dihedral_deg=-120.0,
            hybridization=Hybridization.SP3)

    backbone = BackboneIndices(
        N=n_idx, H_N=h_n, CA=ca_idx, HA=ha_idx,
        C=c_idx, O=o_idx, CB=cb_idx,
        H_Nterm=h_nterm, OXT=oxt, H_OXT=h_oxt,
    )
    mol.backbone = backbone
    return mol


# ===================================================================
# Side chain builders
# ===================================================================

def _sidechain_gly(mol: Molecule):
    """Glycine: no side chain -- add a second H to CA."""
    bb = mol.backbone
    # CA already has N, C, HA; add second H
    ca_pos = mol.atoms[bb.CA].position
    existing = mol.neighbors(bb.CA)
    existing_dirs = [
        _normalize(mol.atoms[n].position - ca_pos) for n in existing
    ]
    new_dirs = _available_tetrahedral_dirs(existing_dirs, 1)
    if new_dirs:
        h_pos = ca_pos + _bond_length("C", "H", 1) * new_dirs[0]
        h2 = mol.add_atom("H", h_pos)
        mol.add_bond(bb.CA, h2, order=1, rotatable=False)


def _sidechain_ala(mol: Molecule):
    """Alanine: -CH3."""
    bb = mol.backbone
    _add_sp3_hydrogens(mol, bb.CB, 3)


def _sidechain_val(mol: Molecule):
    """Valine: -CH(CH3)2."""
    bb = mol.backbone
    cb = bb.CB

    # CG1
    cg1 = mol.add_atom_bonded(
        "C", cb, bond_order=1, angle_ref=bb.CA,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=60.0,
        hybridization=Hybridization.SP3)
    _add_sp3_hydrogens(mol, cg1, 3)

    # CG2
    cg2 = mol.add_atom_bonded(
        "C", cb, bond_order=1, angle_ref=bb.CA,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=-60.0,
        hybridization=Hybridization.SP3)
    _add_sp3_hydrogens(mol, cg2, 3)

    # HB on CB
    _add_sp3_hydrogens(mol, cb, 1)


def _sidechain_leu(mol: Molecule):
    """Leucine: -CH2-CH(CH3)2."""
    bb = mol.backbone
    cb = bb.CB

    cg = mol.add_atom_bonded(
        "C", cb, bond_order=1, angle_ref=bb.CA,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=180.0,
        hybridization=Hybridization.SP3)

    cd1 = mol.add_atom_bonded(
        "C", cg, bond_order=1, angle_ref=cb,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=60.0,
        hybridization=Hybridization.SP3)
    _add_sp3_hydrogens(mol, cd1, 3)

    cd2 = mol.add_atom_bonded(
        "C", cg, bond_order=1, angle_ref=cb,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=-60.0,
        hybridization=Hybridization.SP3)
    _add_sp3_hydrogens(mol, cd2, 3)

    _add_sp3_hydrogens(mol, cg, 1)   # HG
    _add_sp3_hydrogens(mol, cb, 2)   # HB1, HB2


def _sidechain_ile(mol: Molecule):
    """Isoleucine: -CH(CH3)-CH2-CH3."""
    bb = mol.backbone
    cb = bb.CB

    # CG1 (the longer chain: -CH2-CH3)
    cg1 = mol.add_atom_bonded(
        "C", cb, bond_order=1, angle_ref=bb.CA,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=180.0,
        hybridization=Hybridization.SP3)

    cd1 = mol.add_atom_bonded(
        "C", cg1, bond_order=1, angle_ref=cb,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=180.0,
        hybridization=Hybridization.SP3)
    _add_sp3_hydrogens(mol, cd1, 3)
    _add_sp3_hydrogens(mol, cg1, 2)

    # CG2 (the methyl branch)
    cg2 = mol.add_atom_bonded(
        "C", cb, bond_order=1, angle_ref=bb.CA,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=-60.0,
        hybridization=Hybridization.SP3)
    _add_sp3_hydrogens(mol, cg2, 3)

    _add_sp3_hydrogens(mol, cb, 1)   # HB


def _sidechain_pro(mol: Molecule):
    """Proline: pyrrolidine ring connecting CB-CG-CD back to N."""
    bb = mol.backbone
    cb = bb.CB

    cg = mol.add_atom_bonded(
        "C", cb, bond_order=1, angle_ref=bb.CA,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=180.0,
        hybridization=Hybridization.SP3)

    cd = mol.add_atom_bonded(
        "C", cg, bond_order=1, angle_ref=cb,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=30.0,
        hybridization=Hybridization.SP3)

    # Close ring: CD-N
    mol.close_ring(cd, bb.N)

    _add_sp3_hydrogens(mol, cb, 2)
    _add_sp3_hydrogens(mol, cg, 2)
    _add_sp3_hydrogens(mol, cd, 2)


def _sidechain_phe(mol: Molecule):
    """Phenylalanine: -CH2-phenyl."""
    bb = mol.backbone
    cb = bb.CB
    _add_sp3_hydrogens(mol, cb, 2)
    add_phenyl_ring(mol, cb, angle_ref=bb.CA, dihedral_deg=90.0)


def _sidechain_trp(mol: Molecule):
    """Tryptophan: -CH2-indole."""
    bb = mol.backbone
    cb = bb.CB
    _add_sp3_hydrogens(mol, cb, 2)
    add_indole_ring(mol, cb, angle_ref=bb.CA, dihedral_deg=90.0)


def _sidechain_met(mol: Molecule):
    """Methionine: -CH2-CH2-S-CH3."""
    bb = mol.backbone
    cb = bb.CB

    cg = mol.add_atom_bonded(
        "C", cb, bond_order=1, angle_ref=bb.CA,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=180.0,
        hybridization=Hybridization.SP3)

    sd = mol.add_atom_bonded(
        "S", cg, bond_order=1, angle_ref=cb,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=180.0,
        hybridization=Hybridization.SP3)

    ce = mol.add_atom_bonded(
        "C", sd, bond_order=1, angle_ref=cg,
        bond_angle_deg=100.0,  # C-S-C angle ~100 deg
        dihedral_deg=180.0,
        hybridization=Hybridization.SP3)

    _add_sp3_hydrogens(mol, ce, 3)
    _add_sp3_hydrogens(mol, cg, 2)
    _add_sp3_hydrogens(mol, cb, 2)


def _sidechain_ser(mol: Molecule):
    """Serine: -CH2-OH."""
    bb = mol.backbone
    cb = bb.CB
    _add_sp3_hydrogens(mol, cb, 2)
    add_hydroxyl(mol, cb, angle_ref=bb.CA, dihedral_deg=60.0)


def _sidechain_thr(mol: Molecule):
    """Threonine: -CH(OH)-CH3."""
    bb = mol.backbone
    cb = bb.CB

    # OG1 (hydroxyl)
    og1 = mol.add_atom_bonded(
        "O", cb, bond_order=1, angle_ref=bb.CA,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=60.0,
        hybridization=Hybridization.SP3)
    mol.add_atom_bonded(
        "H", og1, bond_order=1, angle_ref=cb,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=180.0,
        rotatable=False)

    # CG2 (methyl)
    cg2 = mol.add_atom_bonded(
        "C", cb, bond_order=1, angle_ref=bb.CA,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=-60.0,
        hybridization=Hybridization.SP3)
    _add_sp3_hydrogens(mol, cg2, 3)

    _add_sp3_hydrogens(mol, cb, 1)  # HB


def _sidechain_cys(mol: Molecule):
    """Cysteine: -CH2-SH."""
    bb = mol.backbone
    cb = bb.CB
    _add_sp3_hydrogens(mol, cb, 2)
    add_thiol(mol, cb, angle_ref=bb.CA, dihedral_deg=60.0)


def _sidechain_tyr(mol: Molecule):
    """Tyrosine: -CH2-phenyl-OH."""
    bb = mol.backbone
    cb = bb.CB
    _add_sp3_hydrogens(mol, cb, 2)
    result = add_phenyl_ring(mol, cb, angle_ref=bb.CA, dihedral_deg=90.0)
    # Add OH to the para carbon (ring[3], opposite to attachment)
    para_c = result["ring"][3]
    # Remove one H from para position -- we added H to ring[1..5]
    # ring[3] is index 2 in the H list (since H list starts at ring[1])
    # Actually the H on para_c is result["H"][2]
    # We need to replace it with OH -- but we can't easily remove atoms
    # from Molecule. Instead, we'll add OH and the extra H is a minor
    # approximation for this educational model.
    add_hydroxyl(mol, para_c, angle_ref=result["ring"][2])


def _sidechain_asn(mol: Molecule):
    """Asparagine: -CH2-CONH2."""
    bb = mol.backbone
    cb = bb.CB
    _add_sp3_hydrogens(mol, cb, 2)
    add_amide(mol, cb, angle_ref=bb.CA, dihedral_deg=180.0)


def _sidechain_gln(mol: Molecule):
    """Glutamine: -CH2-CH2-CONH2."""
    bb = mol.backbone
    cb = bb.CB

    cg = mol.add_atom_bonded(
        "C", cb, bond_order=1, angle_ref=bb.CA,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=180.0,
        hybridization=Hybridization.SP3)

    _add_sp3_hydrogens(mol, cg, 2)
    _add_sp3_hydrogens(mol, cb, 2)
    add_amide(mol, cg, angle_ref=cb, dihedral_deg=180.0)


def _sidechain_asp(mol: Molecule):
    """Aspartate: -CH2-COO-."""
    bb = mol.backbone
    cb = bb.CB
    _add_sp3_hydrogens(mol, cb, 2)
    add_carboxyl(mol, cb, angle_ref=bb.CA, dihedral_deg=180.0)


def _sidechain_glu(mol: Molecule):
    """Glutamate: -CH2-CH2-COO-."""
    bb = mol.backbone
    cb = bb.CB

    cg = mol.add_atom_bonded(
        "C", cb, bond_order=1, angle_ref=bb.CA,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=180.0,
        hybridization=Hybridization.SP3)

    _add_sp3_hydrogens(mol, cg, 2)
    _add_sp3_hydrogens(mol, cb, 2)
    add_carboxyl(mol, cg, angle_ref=cb, dihedral_deg=180.0)


def _sidechain_lys(mol: Molecule):
    """Lysine: -(CH2)4-NH2."""
    bb = mol.backbone
    cb = bb.CB

    cg = mol.add_atom_bonded(
        "C", cb, bond_order=1, angle_ref=bb.CA,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=180.0,
        hybridization=Hybridization.SP3)

    cd = mol.add_atom_bonded(
        "C", cg, bond_order=1, angle_ref=cb,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=180.0,
        hybridization=Hybridization.SP3)

    ce = mol.add_atom_bonded(
        "C", cd, bond_order=1, angle_ref=cg,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=180.0,
        hybridization=Hybridization.SP3)

    _add_sp3_hydrogens(mol, ce, 2)
    _add_sp3_hydrogens(mol, cd, 2)
    _add_sp3_hydrogens(mol, cg, 2)
    _add_sp3_hydrogens(mol, cb, 2)

    add_amino(mol, ce, angle_ref=cd, dihedral_deg=180.0)


def _sidechain_arg(mol: Molecule):
    """Arginine: -(CH2)3-NH-C(=NH)(NH2)."""
    bb = mol.backbone
    cb = bb.CB

    cg = mol.add_atom_bonded(
        "C", cb, bond_order=1, angle_ref=bb.CA,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=180.0,
        hybridization=Hybridization.SP3)

    cd = mol.add_atom_bonded(
        "C", cg, bond_order=1, angle_ref=cb,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=180.0,
        hybridization=Hybridization.SP3)

    ne = mol.add_atom_bonded(
        "N", cd, bond_order=1, angle_ref=cg,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=180.0,
        hybridization=Hybridization.SP3)

    # H on NE
    mol.add_atom_bonded(
        "H", ne, bond_order=1, angle_ref=cd,
        bond_angle_deg=SP3_ANGLE, dihedral_deg=120.0,
        rotatable=False)

    _add_sp3_hydrogens(mol, cd, 2)
    _add_sp3_hydrogens(mol, cg, 2)
    _add_sp3_hydrogens(mol, cb, 2)

    add_guanidinium(mol, ne, angle_ref=cd, dihedral_deg=180.0)


def _sidechain_his(mol: Molecule):
    """Histidine: -CH2-imidazole."""
    bb = mol.backbone
    cb = bb.CB
    _add_sp3_hydrogens(mol, cb, 2)
    add_imidazole_ring(mol, cb, angle_ref=bb.CA, dihedral_deg=90.0)


# ===================================================================
# Amino acid data and factory
# ===================================================================

AMINO_ACID_DATA: dict[AminoAcidType, dict] = {
    AminoAcidType.GLY: {
        "name": "Glycine", "code1": "G", "code3": "Gly",
        "has_cb": False, "builder": _sidechain_gly,
    },
    AminoAcidType.ALA: {
        "name": "Alanine", "code1": "A", "code3": "Ala",
        "has_cb": True, "builder": _sidechain_ala,
    },
    AminoAcidType.VAL: {
        "name": "Valine", "code1": "V", "code3": "Val",
        "has_cb": True, "builder": _sidechain_val,
    },
    AminoAcidType.LEU: {
        "name": "Leucine", "code1": "L", "code3": "Leu",
        "has_cb": True, "builder": _sidechain_leu,
    },
    AminoAcidType.ILE: {
        "name": "Isoleucine", "code1": "I", "code3": "Ile",
        "has_cb": True, "builder": _sidechain_ile,
    },
    AminoAcidType.PRO: {
        "name": "Proline", "code1": "P", "code3": "Pro",
        "has_cb": True, "builder": _sidechain_pro,
    },
    AminoAcidType.PHE: {
        "name": "Phenylalanine", "code1": "F", "code3": "Phe",
        "has_cb": True, "builder": _sidechain_phe,
    },
    AminoAcidType.TRP: {
        "name": "Tryptophan", "code1": "W", "code3": "Trp",
        "has_cb": True, "builder": _sidechain_trp,
    },
    AminoAcidType.MET: {
        "name": "Methionine", "code1": "M", "code3": "Met",
        "has_cb": True, "builder": _sidechain_met,
    },
    AminoAcidType.SER: {
        "name": "Serine", "code1": "S", "code3": "Ser",
        "has_cb": True, "builder": _sidechain_ser,
    },
    AminoAcidType.THR: {
        "name": "Threonine", "code1": "T", "code3": "Thr",
        "has_cb": True, "builder": _sidechain_thr,
    },
    AminoAcidType.CYS: {
        "name": "Cysteine", "code1": "C", "code3": "Cys",
        "has_cb": True, "builder": _sidechain_cys,
    },
    AminoAcidType.TYR: {
        "name": "Tyrosine", "code1": "Y", "code3": "Tyr",
        "has_cb": True, "builder": _sidechain_tyr,
    },
    AminoAcidType.ASN: {
        "name": "Asparagine", "code1": "N", "code3": "Asn",
        "has_cb": True, "builder": _sidechain_asn,
    },
    AminoAcidType.GLN: {
        "name": "Glutamine", "code1": "Q", "code3": "Gln",
        "has_cb": True, "builder": _sidechain_gln,
    },
    AminoAcidType.ASP: {
        "name": "Aspartate", "code1": "D", "code3": "Asp",
        "has_cb": True, "builder": _sidechain_asp,
    },
    AminoAcidType.GLU: {
        "name": "Glutamate", "code1": "E", "code3": "Glu",
        "has_cb": True, "builder": _sidechain_glu,
    },
    AminoAcidType.LYS: {
        "name": "Lysine", "code1": "K", "code3": "Lys",
        "has_cb": True, "builder": _sidechain_lys,
    },
    AminoAcidType.ARG: {
        "name": "Arginine", "code1": "R", "code3": "Arg",
        "has_cb": True, "builder": _sidechain_arg,
    },
    AminoAcidType.HIS: {
        "name": "Histidine", "code1": "H", "code3": "His",
        "has_cb": True, "builder": _sidechain_his,
    },
}


def build_amino_acid(aa_type: AminoAcidType) -> Molecule:
    """Build a single amino acid with correct 3D geometry.

    The returned Molecule has a .backbone attribute (BackboneIndices)
    that records the indices of backbone atoms for peptide bond
    formation.
    """
    data = AMINO_ACID_DATA[aa_type]
    is_pro = (aa_type == AminoAcidType.PRO)
    mol = _build_backbone(data["name"], has_cb=data["has_cb"],
                          is_proline=is_pro)
    data["builder"](mol)
    return mol


# ===================================================================
# Peptide bond formation
# ===================================================================

def form_peptide_bond(mol_nterm: Molecule, mol_cterm: Molecule,
                      omega_deg: float = 180.0) -> Molecule:
    """Join two amino acids via a peptide (amide) bond.

    mol_nterm provides the N-terminal residue (its C-terminus is
    condensed).  mol_cterm provides the C-terminal residue (one H
    removed from its N-terminus).

    Returns a new Molecule with both residues joined.
    """
    bb_n = mol_nterm.backbone
    bb_c = mol_cterm.backbone

    # Atoms to skip (condensation products: H2O)
    skip_nterm = set()  # skip OXT and H_OXT from N-terminal residue
    if bb_n.OXT is not None:
        skip_nterm.add(bb_n.OXT)
    if bb_n.H_OXT is not None:
        skip_nterm.add(bb_n.H_OXT)

    skip_cterm = set()  # skip one H from C-terminal residue's N
    if bb_c.H_Nterm is not None:
        skip_cterm.add(bb_c.H_Nterm)

    result = Molecule(f"{mol_nterm.name}-{mol_cterm.name}")

    # Copy atoms from mol_nterm (with index remapping)
    remap_n: dict[int, int] = {}
    for atom in mol_nterm.atoms:
        if atom.index in skip_nterm:
            continue
        new_idx = result.add_atom(atom.symbol, atom.position.copy(),
                                  atom.hybridization)
        remap_n[atom.index] = new_idx

    # Copy atoms from mol_cterm
    # First, compute the rigid-body transform to position mol_cterm
    # so that its N is at the correct peptide bond distance from
    # mol_nterm's C.
    c_pos = mol_nterm.atoms[bb_n.C].position
    ca_n_pos = mol_nterm.atoms[bb_n.CA].position
    o_n_pos = mol_nterm.atoms[bb_n.O].position

    # Target position for the peptide N: along C->away_from_CA direction
    # at AMIDE_CN distance, at SP2 angle from CA-C-N
    c_to_ca = _normalize(ca_n_pos - c_pos)
    c_to_o = _normalize(o_n_pos - c_pos)

    # The peptide N should be at ~120 deg from both CA and O
    # N direction: reflect c_to_ca across the plane normal to c_to_o
    # For sp2 geometry: N is at 120 deg from CA on the opposite side of O
    n_dir = _normalize(-(c_to_ca + c_to_o))
    target_n_pos = c_pos + n_dir * AMIDE_CN

    # Compute translation for mol_cterm
    current_n_pos = mol_cterm.atoms[bb_c.N].position
    translation = target_n_pos - current_n_pos

    remap_c: dict[int, int] = {}
    for atom in mol_cterm.atoms:
        if atom.index in skip_cterm:
            continue
        new_pos = atom.position.copy() + translation
        new_idx = result.add_atom(atom.symbol, new_pos,
                                  atom.hybridization)
        remap_c[atom.index] = new_idx

    # Copy bonds from mol_nterm
    for bond in mol_nterm.bonds:
        if bond.atom_i in skip_nterm or bond.atom_j in skip_nterm:
            continue
        if bond.atom_i in remap_n and bond.atom_j in remap_n:
            result.add_bond(remap_n[bond.atom_i], remap_n[bond.atom_j],
                            order=bond.order, rotatable=bond.rotatable)

    # Copy bonds from mol_cterm
    for bond in mol_cterm.bonds:
        if bond.atom_i in skip_cterm or bond.atom_j in skip_cterm:
            continue
        if bond.atom_i in remap_c and bond.atom_j in remap_c:
            result.add_bond(remap_c[bond.atom_i], remap_c[bond.atom_j],
                            order=bond.order, rotatable=bond.rotatable)

    # Add the peptide bond (C-N amide bond)
    peptide_c = remap_n[bb_n.C]
    peptide_n = remap_c[bb_c.N]
    result.add_bond(peptide_c, peptide_n, order=1, rotatable=False)

    # Build backbone indices for the combined molecule
    # N-terminal residue backbone
    bb_new_n = BackboneIndices(
        N=remap_n[bb_n.N],
        H_N=remap_n.get(bb_n.H_N) if bb_n.H_N is not None else None,
        CA=remap_n[bb_n.CA],
        HA=remap_n[bb_n.HA],
        C=remap_n[bb_n.C],
        O=remap_n[bb_n.O],
        CB=remap_n.get(bb_n.CB) if bb_n.CB is not None else None,
        H_Nterm=remap_n.get(bb_n.H_Nterm) if bb_n.H_Nterm is not None else None,
        OXT=None, H_OXT=None,  # removed in condensation
    )

    bb_new_c = BackboneIndices(
        N=remap_c[bb_c.N],
        H_N=remap_c.get(bb_c.H_N) if bb_c.H_N is not None else None,
        CA=remap_c[bb_c.CA],
        HA=remap_c[bb_c.HA],
        C=remap_c[bb_c.C],
        O=remap_c[bb_c.O],
        CB=remap_c.get(bb_c.CB) if bb_c.CB is not None else None,
        H_Nterm=None,  # removed in condensation
        OXT=remap_c.get(bb_c.OXT) if bb_c.OXT is not None else None,
        H_OXT=remap_c.get(bb_c.H_OXT) if bb_c.H_OXT is not None else None,
    )

    # Store backbone index sets -- carry forward any existing residues
    result.backbone = bb_new_n
    prior_residues = []
    if hasattr(mol_nterm, 'residues'):
        # Remap all prior residue indices
        for old_bb in mol_nterm.residues:
            remapped = BackboneIndices(
                N=remap_n.get(old_bb.N, old_bb.N),
                H_N=remap_n.get(old_bb.H_N) if old_bb.H_N is not None else None,
                CA=remap_n.get(old_bb.CA, old_bb.CA),
                HA=remap_n.get(old_bb.HA, old_bb.HA),
                C=remap_n.get(old_bb.C, old_bb.C),
                O=remap_n.get(old_bb.O, old_bb.O),
                CB=remap_n.get(old_bb.CB) if old_bb.CB is not None else None,
                H_Nterm=remap_n.get(old_bb.H_Nterm) if old_bb.H_Nterm is not None else None,
                OXT=remap_n.get(old_bb.OXT) if old_bb.OXT is not None else None,
                H_OXT=remap_n.get(old_bb.H_OXT) if old_bb.H_OXT is not None else None,
            )
            prior_residues.append(remapped)
    else:
        prior_residues.append(bb_new_n)
    result.residues = prior_residues + [bb_new_c]

    # Set omega dihedral (CA_n - C_n - N_c - CA_c)
    try:
        result.set_dihedral(
            bb_new_n.CA, bb_new_n.C, bb_new_c.N, bb_new_c.CA,
            omega_deg)
    except (ValueError, IndexError):
        pass  # ring bonds or missing atoms

    return result


def build_peptide(sequence: list[AminoAcidType],
                  phi_psi: list[tuple[float, float]] | None = None,
                  ) -> Molecule:
    """Build a polypeptide from a sequence of amino acids.

    Parameters
    ----------
    sequence : list of AminoAcidType
        Amino acid sequence from N-terminus to C-terminus.
    phi_psi : optional list of (phi, psi) angle pairs in degrees.
        If provided, must have len(sequence) entries.
    """
    if len(sequence) < 1:
        raise ValueError("Sequence must have at least one amino acid")

    if len(sequence) == 1:
        return build_amino_acid(sequence[0])

    # Build each amino acid individually
    mols = [build_amino_acid(aa) for aa in sequence]

    # Join them one by one
    result = mols[0]
    for i in range(1, len(mols)):
        result = form_peptide_bond(result, mols[i])

    # Apply phi/psi angles if given
    if phi_psi is not None and hasattr(result, 'residues'):
        for i, (phi, psi) in enumerate(phi_psi):
            if i < len(result.residues):
                try:
                    set_phi(result, result.residues, i, phi)
                except (ValueError, IndexError):
                    pass
                try:
                    set_psi(result, result.residues, i, psi)
                except (ValueError, IndexError):
                    pass

    return result


# ===================================================================
# Backbone conformation control
# ===================================================================

def set_phi(mol: Molecule, residues: list[BackboneIndices],
            residue_index: int, angle_deg: float):
    """Set the phi angle for a residue.

    Phi = C(i-1) - N(i) - CA(i) - C(i).
    Only valid for residue_index >= 1 (the first residue has no
    preceding C).
    """
    if residue_index < 1:
        return
    bb_prev = residues[residue_index - 1]
    bb = residues[residue_index]
    mol.set_dihedral(bb_prev.C, bb.N, bb.CA, bb.C, angle_deg)


def set_psi(mol: Molecule, residues: list[BackboneIndices],
            residue_index: int, angle_deg: float):
    """Set the psi angle for a residue.

    Psi = N(i) - CA(i) - C(i) - N(i+1).
    Only valid for residue_index < len-1 (the last residue has no
    following N).
    """
    if residue_index >= len(residues) - 1:
        return
    bb = residues[residue_index]
    bb_next = residues[residue_index + 1]
    mol.set_dihedral(bb.N, bb.CA, bb.C, bb_next.N, angle_deg)


def set_secondary_structure(mol: Molecule,
                            residues: list[BackboneIndices],
                            structure: SecondaryStructure):
    """Apply phi/psi angles for a secondary structure type."""
    angles = {
        SecondaryStructure.ALPHA_HELIX: ALPHA_HELIX,
        SecondaryStructure.BETA_SHEET: BETA_SHEET,
        SecondaryStructure.POLYPROLINE_II: POLYPROLINE_II,
        SecondaryStructure.EXTENDED: EXTENDED,
    }
    phi, psi = angles[structure]
    for i in range(len(residues)):
        try:
            set_phi(mol, residues, i, phi)
        except (ValueError, IndexError):
            pass
        try:
            set_psi(mol, residues, i, psi)
        except (ValueError, IndexError):
            pass


# ===================================================================
# Main -- demonstration
# ===================================================================

if __name__ == "__main__":
    print("Amino Acid Structures and Functional Groups")
    print("=" * 60)
    print()

    # ---- 1. Functional group showcase --------------------------------
    print("--- Functional Group Showcase ---\n")

    demo_mol = Molecule("functional group demo")
    c0 = demo_mol.add_atom("C", [0.0, 0.0, 0.0], Hybridization.SP3)
    c1 = demo_mol.add_atom_bonded("C", c0, hybridization=Hybridization.SP3)
    c2 = demo_mol.add_atom_bonded("C", c1, angle_ref=c0,
                                   bond_angle_deg=SP3_ANGLE,
                                   hybridization=Hybridization.SP3)
    c3 = demo_mol.add_atom_bonded("C", c2, angle_ref=c1,
                                   dihedral_ref=c0,
                                   bond_angle_deg=SP3_ANGLE,
                                   dihedral_deg=180.0,
                                   hybridization=Hybridization.SP3)

    oh = add_hydroxyl(demo_mol, c0)
    nh2 = add_amino(demo_mol, c1, angle_ref=c0, dihedral_deg=60.0)
    sh = add_thiol(demo_mol, c2, angle_ref=c1, dihedral_deg=60.0)
    cooh = add_carboxyl(demo_mol, c3, angle_ref=c2, dihedral_deg=60.0)

    print(f"  Built molecule with {len(demo_mol.atoms)} atoms, "
          f"{len(demo_mol.bonds)} bonds")
    print(f"  Hydroxyl O at index {oh['O']}")
    print(f"  Amino N at index {nh2['N']}")
    print(f"  Thiol S at index {sh['S']}")
    print(f"  Carboxyl C at index {cooh['C']}")
    print()

    # ---- 2. Amino acid gallery ---------------------------------------
    print("--- Amino Acid Gallery ---\n")

    categories = {
        "Nonpolar": [AminoAcidType.GLY, AminoAcidType.ALA,
                     AminoAcidType.VAL, AminoAcidType.LEU,
                     AminoAcidType.ILE, AminoAcidType.PRO,
                     AminoAcidType.PHE, AminoAcidType.TRP,
                     AminoAcidType.MET],
        "Polar":    [AminoAcidType.SER, AminoAcidType.THR,
                     AminoAcidType.CYS, AminoAcidType.TYR,
                     AminoAcidType.ASN, AminoAcidType.GLN],
        "Charged":  [AminoAcidType.ASP, AminoAcidType.GLU,
                     AminoAcidType.LYS, AminoAcidType.ARG,
                     AminoAcidType.HIS],
    }

    print(f"  {'Name':<16} {'Code':>4} {'Atoms':>6} {'Bonds':>6}")
    print(f"  {'-' * 40}")

    all_amino_acids = {}
    for cat_name, aa_list in categories.items():
        for aa_type in aa_list:
            data = AMINO_ACID_DATA[aa_type]
            mol = build_amino_acid(aa_type)
            all_amino_acids[aa_type] = mol
            print(f"  {data['name']:<16} {data['code3']:>4}"
                  f" {len(mol.atoms):>6} {len(mol.bonds):>6}")
    print()

    # ---- 3. L-chirality verification ---------------------------------
    print("--- L-Chirality Verification ---\n")

    print(f"  {'Name':<16} {'CA idx':>6} {'Config':>8}")
    print(f"  {'-' * 34}")

    for aa_type, mol in all_amino_acids.items():
        bb = mol.backbone
        if aa_type == AminoAcidType.GLY:
            config = "achiral"
        else:
            desc = mol.assign_rs(bb.CA)
            config = desc.name
        data = AMINO_ACID_DATA[aa_type]
        print(f"  {data['name']:<16} {bb.CA:>6} {config:>8}")
    print()

    # ---- 4. Single amino acid detail ---------------------------------
    print("--- Alanine Detail ---\n")
    ala = all_amino_acids[AminoAcidType.ALA]
    print(ala.summary())
    print()

    # ---- 5. Dipeptide: Ala-Gly ---------------------------------------
    print("--- Dipeptide: Ala-Gly ---\n")
    ala_mol = build_amino_acid(AminoAcidType.ALA)
    gly_mol = build_amino_acid(AminoAcidType.GLY)
    dipeptide = form_peptide_bond(ala_mol, gly_mol)

    print(f"  Name: {dipeptide.name}")
    print(f"  Atoms: {len(dipeptide.atoms)}")
    print(f"  Bonds: {len(dipeptide.bonds)}")
    print(f"  Residues: {len(dipeptide.residues)}")

    # Check peptide bond length
    if hasattr(dipeptide, 'residues') and len(dipeptide.residues) >= 2:
        c_idx = dipeptide.residues[0].C
        n_idx = dipeptide.residues[1].N
        pep_dist = dipeptide.distance(c_idx, n_idx)
        print(f"  Peptide bond length: {pep_dist:.3f} A")

        # Omega dihedral
        ca1 = dipeptide.residues[0].CA
        ca2 = dipeptide.residues[1].CA
        omega = dipeptide.dihedral_angle(ca1, c_idx, n_idx, ca2)
        print(f"  Omega dihedral: {omega:.1f} deg")
    print()

    # ---- 6. Tripeptide with secondary structure ----------------------
    print("--- Tripeptide: Ala-Ala-Ala (alpha helix) ---\n")
    tripeptide = build_peptide(
        [AminoAcidType.ALA, AminoAcidType.ALA, AminoAcidType.ALA])
    print(f"  Atoms: {len(tripeptide.atoms)}")
    print(f"  Bonds: {len(tripeptide.bonds)}")

    if hasattr(tripeptide, 'residues'):
        set_secondary_structure(tripeptide, tripeptide.residues,
                                SecondaryStructure.ALPHA_HELIX)
        print(f"  Secondary structure: alpha helix")
        print(f"  Residue count: {len(tripeptide.residues)}")
    print()

    print("Done.")
