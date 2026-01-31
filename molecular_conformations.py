"""
Molecular Conformations and Stereochemistry

Represents multi-centre molecules with explicit atom connectivity and 3D
coordinates built from internal coordinates (bond lengths, bond angles,
dihedral angles).  Supports:

    - Z-matrix style atom placement for multi-centre molecules
    - Dihedral rotation about single bonds
    - Named conformations (staggered, eclipsed, gauche, anti)
    - Chair / boat cyclohexane
    - E/Z geometric isomerism and R/S chirality
    - Torsional strain energy estimation (Pitzer potential)
    - Newman projection data extraction

Complements the existing VSEPR model (single-central-atom molecules) by
handling chains, branches, and rings such as ethane, butane, and
cyclohexane.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum, auto

import numpy as np

from bohr_model import SYMBOL_TO_Z
from element_data import (
    estimated_bond_length_angstrom,
    cpk_color,
)


# ===================================================================
# Constants -- standard bond geometries
# ===================================================================

# Standard bond lengths in Angstroms (experimental averages).
# Keys are (symbol_a, symbol_b, order) with symbols in alphabetical order.
STANDARD_BOND_LENGTHS: dict[tuple[str, str, int], float] = {
    ("C", "C", 1): 1.54,
    ("C", "C", 2): 1.34,
    ("C", "C", 3): 1.20,
    ("C", "H", 1): 1.09,
    ("C", "O", 1): 1.43,
    ("C", "O", 2): 1.23,
    ("C", "N", 1): 1.47,
    ("C", "N", 2): 1.29,
    ("C", "N", 3): 1.16,
    ("C", "F", 1): 1.35,
    ("C", "Cl", 1): 1.77,
    ("Br", "C", 1): 1.94,
    ("C", "I", 1): 2.14,
    ("C", "S", 1): 1.82,
    ("H", "O", 1): 0.96,
    ("H", "N", 1): 1.01,
    ("H", "S", 1): 1.34,
    ("F", "F", 1): 1.42,
    ("Cl", "Cl", 1): 1.99,
    ("Br", "Br", 1): 2.28,
    ("H", "H", 1): 0.74,
    ("S", "S", 1): 2.05,
    ("N", "O", 1): 1.36,
    ("N", "O", 2): 1.21,
    ("O", "P", 1): 1.63,
    ("O", "P", 2): 1.48,
    ("H", "P", 1): 1.44,
}

# Standard bond angles in degrees
SP3_ANGLE = 109.47   # tetrahedral
SP2_ANGLE = 120.0    # trigonal planar
SP_ANGLE = 180.0     # linear

# Torsion barrier parameters (kJ/mol) for the Pitzer potential.
# The per-pair energy uses OPLS-AA sign conventions:
#   V(phi) = (V1/2)(1 + cos(phi))       -- onefold, max at 0 (eclipsed)
#          + (V2/2)(1 - cos(2*phi))      -- twofold, max at 90
#          + (V3/2)(1 + cos(3*phi))      -- threefold, max at 0/120/240
# Parameters are per dihedral pair; the total is summed over all
# substituent pairs across the bond.
TORSION_BARRIERS: dict[str, dict[str, float]] = {
    "H_sp3_sp3_H": {"V1": 0.0,  "V2": 0.0,   "V3": 1.39},
    "C_sp3_sp3_C": {"V1": 2.73, "V2": -0.53,  "V3": 0.84},
    "C_sp3_sp3_H": {"V1": 0.0,  "V2": 0.0,    "V3": 0.76},
    "H_sp3_sp3_C": {"V1": 0.0,  "V2": 0.0,    "V3": 0.76},
    "default":     {"V1": 0.0,  "V2": 0.0,    "V3": 1.00},
}


# ===================================================================
# 3D geometry helpers
# ===================================================================

def _normalize(v: np.ndarray) -> np.ndarray:
    """Return unit vector, or zero vector if input is near-zero."""
    n = np.linalg.norm(v)
    if n < 1e-12:
        return np.zeros(3)
    return v / n


def _rotation_matrix(axis: np.ndarray, theta: float) -> np.ndarray:
    """Rodrigues rotation matrix: rotate by *theta* radians about *axis*."""
    u = _normalize(axis)
    K = np.array([
        [0.0,   -u[2],  u[1]],
        [u[2],   0.0,  -u[0]],
        [-u[1],  u[0],  0.0],
    ])
    return np.eye(3) + math.sin(theta) * K + (1.0 - math.cos(theta)) * (K @ K)


def _place_atom_zmatrix(pos_ref: np.ndarray,
                         pos_angle_ref: np.ndarray,
                         pos_dihedral_ref: np.ndarray,
                         bond_length: float,
                         bond_angle_deg: float,
                         dihedral_deg: float) -> np.ndarray:
    """Place an atom using internal (z-matrix) coordinates.

    Given three reference positions (j, i, k) place a new atom m so
    that distance(j,m)=bond_length, angle(i-j-m)=bond_angle_deg, and
    dihedral(k-i-j-m)=dihedral_deg.

    Parameters
    ----------
    pos_ref : position of bonded atom j.
    pos_angle_ref : position of atom i (defines bond angle).
    pos_dihedral_ref : position of atom k (defines dihedral).
    bond_length, bond_angle_deg, dihedral_deg : internal coordinates.
    """
    theta = math.radians(bond_angle_deg)
    phi = math.radians(dihedral_deg)

    v_ij = _normalize(pos_ref - pos_angle_ref)
    v_ki = _normalize(pos_angle_ref - pos_dihedral_ref)

    # Plane normal
    n = np.cross(v_ki, v_ij)
    if np.linalg.norm(n) < 1e-10:
        # Collinear fallback
        perp = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(v_ij, perp)) > 0.9:
            perp = np.array([0.0, 1.0, 0.0])
        n = np.cross(v_ij, perp)
    n = _normalize(n)

    # In-plane perpendicular
    d = _normalize(np.cross(n, v_ij))

    new_dir = (
        -v_ij * math.cos(theta)
        + d * math.sin(theta) * math.cos(phi)
        + n * math.sin(theta) * math.sin(phi)
    )
    return pos_ref + bond_length * _normalize(new_dir)


def _bond_length(sym_a: str, sym_b: str, order: int = 1) -> float:
    """Look up standard bond length in Angstroms.

    Checks STANDARD_BOND_LENGTHS first, falls back to the
    element_data covalent-radii estimate.
    """
    a, b = sorted([sym_a, sym_b])
    key = (a, b, order)
    if key in STANDARD_BOND_LENGTHS:
        return STANDARD_BOND_LENGTHS[key]
    return estimated_bond_length_angstrom(sym_a, sym_b, order)


# ===================================================================
# Enumerations
# ===================================================================

class Hybridization(Enum):
    """Hybridisation state of a bonding centre."""
    SP3 = auto()
    SP2 = auto()
    SP = auto()


class ConformationType(Enum):
    """Named conformation types for rotational isomers."""
    ECLIPSED = auto()
    STAGGERED = auto()
    GAUCHE = auto()
    ANTI = auto()
    CUSTOM = auto()


class RingConformation(Enum):
    """Named ring conformations."""
    CHAIR = auto()
    BOAT = auto()
    TWIST_BOAT = auto()
    HALF_CHAIR = auto()
    FLAT = auto()


class Stereodescriptor(Enum):
    """Stereochemical configuration descriptors."""
    R = auto()
    S = auto()
    E = auto()
    Z = auto()
    NONE = auto()


# ===================================================================
# Core data classes
# ===================================================================

@dataclass
class Atom:
    """An atom with a 3D position in a multi-centre molecule."""
    symbol: str
    position: np.ndarray
    index: int
    hybridization: Hybridization | None = None

    def __repr__(self):
        x, y, z = self.position
        return f"Atom({self.symbol}[{self.index}] @ ({x:.3f}, {y:.3f}, {z:.3f}))"


@dataclass
class Bond:
    """A bond between two atoms in a multi-centre molecule."""
    atom_i: int
    atom_j: int
    order: int = 1
    rotatable: bool = False

    def __repr__(self):
        sym = {1: "-", 2: "=", 3: "#"}.get(self.order, "?")
        rot = " (rot)" if self.rotatable else ""
        return f"Bond({self.atom_i}{sym}{self.atom_j}{rot})"


@dataclass
class TorsionAngle:
    """A dihedral / torsion angle defined by four atom indices."""
    atom_i: int
    atom_j: int
    atom_k: int
    atom_l: int
    angle_deg: float = 0.0

    def __repr__(self):
        return (f"Torsion({self.atom_i}-{self.atom_j}-"
                f"{self.atom_k}-{self.atom_l}: {self.angle_deg:.1f} deg)")


@dataclass
class NewmanProjection:
    """Data for a Newman projection along a bond.

    Substituent angles are measured clockwise from 12-o'clock.
    """
    front_atom: int
    back_atom: int
    front_substituents: list[tuple[int, str, float]] = field(default_factory=list)
    back_substituents: list[tuple[int, str, float]] = field(default_factory=list)
    dihedral_deg: float = 0.0

    def summary(self) -> str:
        lines = [
            f"  Newman Projection along bond "
            f"{self.front_atom}-{self.back_atom}",
            f"  Dihedral: {self.dihedral_deg:.1f} deg",
            f"  Front substituents (atom {self.front_atom}):",
        ]
        for idx, sym, ang in self.front_substituents:
            lines.append(f"    {sym}[{idx}] at {ang:.1f} deg")
        lines.append(f"  Back substituents (atom {self.back_atom}):")
        for idx, sym, ang in self.back_substituents:
            lines.append(f"    {sym}[{idx}] at {ang:.1f} deg")
        return "\n".join(lines)


@dataclass
class StrainEnergy:
    """Torsional strain energy breakdown."""
    total_kj_per_mol: float
    contributions: list[tuple[TorsionAngle, float]] = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            f"  Torsional Strain: {self.total_kj_per_mol:.2f} kJ/mol",
        ]
        for torsion, e in self.contributions:
            lines.append(f"    {torsion}: {e:.2f} kJ/mol")
        return "\n".join(lines)


# ===================================================================
# Molecule -- multi-centre molecular graph with 3D coordinates
# ===================================================================

class Molecule:
    """A multi-centre molecule with explicit connectivity and 3D geometry.

    Unlike VSEPRMolecule (single-central-atom systems), this class can
    represent chains, branches, and rings such as ethane, butane, and
    cyclohexane.

    Atoms are placed incrementally using internal coordinates (z-matrix
    style): bond length, bond angle, and dihedral angle.

    Parameters
    ----------
    name : str
        Human-readable name for the molecule.
    """

    def __init__(self, name: str = ""):
        self.name = name
        self.atoms: list[Atom] = []
        self.bonds: list[Bond] = []
        self._adj: dict[int, list[int]] = {}

    # ---- building ----

    def add_atom(self, symbol: str, position,
                 hybridization: Hybridization | None = None) -> int:
        """Add an atom at an explicit 3D position.  Returns new index."""
        idx = len(self.atoms)
        self.atoms.append(Atom(
            symbol=symbol,
            position=np.array(position, dtype=float),
            index=idx,
            hybridization=hybridization,
        ))
        self._adj[idx] = []
        return idx

    def add_bond(self, i: int, j: int, order: int = 1,
                 rotatable: bool | None = None) -> Bond:
        """Add a bond between atoms *i* and *j*.  Returns the Bond."""
        if rotatable is None:
            rotatable = (order == 1)
        bond = Bond(atom_i=i, atom_j=j, order=order, rotatable=rotatable)
        self.bonds.append(bond)
        self._adj[i].append(j)
        self._adj[j].append(i)
        return bond

    def add_atom_bonded(self, symbol: str, bonded_to: int,
                         bond_order: int = 1,
                         angle_ref: int | None = None,
                         dihedral_ref: int | None = None,
                         bond_length: float | None = None,
                         bond_angle_deg: float | None = None,
                         dihedral_deg: float = 0.0,
                         hybridization: Hybridization | None = None,
                         rotatable: bool | None = None) -> int:
        """Add an atom bonded to an existing atom via internal coordinates.

        For the first atom (index 0): placed at origin.
        For the second atom (index 1): placed along +z.
        For the third: placed in the xz-plane using a synthetic dihedral
        reference.  For all subsequent atoms: full z-matrix placement.

        Returns the index of the new atom.
        """
        if bond_length is None:
            parent_sym = self.atoms[bonded_to].symbol
            bond_length = _bond_length(parent_sym, symbol, bond_order)

        n = len(self.atoms)

        # --- first atom ------------------------------------------------
        if n == 0:
            return self.add_atom(symbol, [0.0, 0.0, 0.0], hybridization)

        # --- second atom -----------------------------------------------
        if n == 1:
            pos = (self.atoms[bonded_to].position
                   + np.array([0.0, 0.0, bond_length]))
            idx = self.add_atom(symbol, pos, hybridization)
            self.add_bond(bonded_to, idx, bond_order, rotatable)
            return idx

        # --- default bond angle from parent hybridisation --------------
        if bond_angle_deg is None:
            parent_hyb = self.atoms[bonded_to].hybridization
            if parent_hyb == Hybridization.SP2:
                bond_angle_deg = SP2_ANGLE
            elif parent_hyb == Hybridization.SP:
                bond_angle_deg = SP_ANGLE
            else:
                bond_angle_deg = SP3_ANGLE

        # --- auto-select angle reference -------------------------------
        if angle_ref is None:
            nbrs = self._adj.get(bonded_to, [])
            if nbrs:
                angle_ref = nbrs[0]
            else:
                angle_ref = 0 if bonded_to != 0 else 1

        # --- third atom (synthetic dihedral ref) -----------------------
        if n == 2:
            pos_j = self.atoms[bonded_to].position
            pos_i = self.atoms[angle_ref].position
            synthetic_k = pos_i + np.array([0.0, 1.0, 0.0])
            pos = _place_atom_zmatrix(pos_j, pos_i, synthetic_k,
                                       bond_length, bond_angle_deg,
                                       dihedral_deg)
            idx = self.add_atom(symbol, pos, hybridization)
            self.add_bond(bonded_to, idx, bond_order, rotatable)
            return idx

        # --- general case: full z-matrix -------------------------------
        if dihedral_ref is None:
            ar_nbrs = self._adj.get(angle_ref, [])
            candidates = [x for x in ar_nbrs if x != bonded_to]
            if candidates:
                dihedral_ref = candidates[0]
            else:
                for a in self.atoms:
                    if a.index not in (bonded_to, angle_ref):
                        dihedral_ref = a.index
                        break

        pos = _place_atom_zmatrix(
            self.atoms[bonded_to].position,
            self.atoms[angle_ref].position,
            self.atoms[dihedral_ref].position,
            bond_length, bond_angle_deg, dihedral_deg,
        )
        idx = self.add_atom(symbol, pos, hybridization)
        self.add_bond(bonded_to, idx, bond_order, rotatable)
        return idx

    def close_ring(self, i: int, j: int, order: int = 1):
        """Bond two existing atoms to close a ring (non-rotatable)."""
        self.add_bond(i, j, order=order, rotatable=False)

    # ---- geometry queries ----

    def distance(self, i: int, j: int) -> float:
        """Distance between atoms *i* and *j* in Angstroms."""
        return float(np.linalg.norm(
            self.atoms[i].position - self.atoms[j].position))

    def bond_angle(self, i: int, j: int, k: int) -> float:
        """Bond angle i-j-k in degrees."""
        vi = self.atoms[i].position - self.atoms[j].position
        vk = self.atoms[k].position - self.atoms[j].position
        cos_a = np.clip(
            np.dot(vi, vk) / (np.linalg.norm(vi) * np.linalg.norm(vk)),
            -1.0, 1.0)
        return math.degrees(math.acos(cos_a))

    def dihedral_angle(self, i: int, j: int, k: int, l: int) -> float:
        """Signed dihedral angle i-j-k-l in degrees (-180 to 180)."""
        b1 = self.atoms[j].position - self.atoms[i].position
        b2 = self.atoms[k].position - self.atoms[j].position
        b3 = self.atoms[l].position - self.atoms[k].position

        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)
        n1_n = np.linalg.norm(n1)
        n2_n = np.linalg.norm(n2)
        if n1_n < 1e-12 or n2_n < 1e-12:
            return 0.0
        n1 /= n1_n
        n2 /= n2_n
        b2_hat = b2 / np.linalg.norm(b2)
        x = float(np.dot(n1, n2))
        y = float(np.dot(np.cross(n1, b2_hat), n2))
        return math.degrees(math.atan2(y, x))

    def neighbors(self, idx: int) -> list[int]:
        """Indices of atoms bonded to atom *idx*."""
        return list(self._adj.get(idx, []))

    def get_bond(self, i: int, j: int) -> Bond | None:
        """Return the Bond between *i* and *j*, or None."""
        for b in self.bonds:
            if {b.atom_i, b.atom_j} == {i, j}:
                return b
        return None

    def is_in_ring(self, i: int, j: int) -> bool:
        """True if removing the i-j edge leaves i and j connected."""
        visited: set[int] = set()
        stack = [i]
        while stack:
            cur = stack.pop()
            if cur == j:
                return True
            if cur in visited:
                continue
            visited.add(cur)
            for nb in self._adj.get(cur, []):
                if cur == i and nb == j:
                    continue
                if cur == j and nb == i:
                    continue
                if nb not in visited:
                    stack.append(nb)
        return False

    # ---- dihedral manipulation ----

    def rotate_dihedral(self, j: int, k: int, angle_deg: float):
        """Rotate all atoms on the k-side of bond j-k by *angle_deg*.

        Atoms on the j-side remain fixed.  Raises ValueError if the
        bond is part of a ring.
        """
        if self.is_in_ring(j, k):
            raise ValueError(
                f"Cannot rotate ring bond {j}-{k}.")

        # BFS from k, not crossing back to j
        k_side: set[int] = set()
        stack = [k]
        while stack:
            cur = stack.pop()
            if cur in k_side:
                continue
            k_side.add(cur)
            for nb in self._adj.get(cur, []):
                if cur == k and nb == j:
                    continue
                if nb not in k_side:
                    stack.append(nb)

        axis = self.atoms[k].position - self.atoms[j].position
        pivot = self.atoms[j].position
        R = _rotation_matrix(axis, math.radians(angle_deg))

        for idx in k_side:
            rel = self.atoms[idx].position - pivot
            self.atoms[idx].position = pivot + R @ rel

    def set_dihedral(self, i: int, j: int, k: int, l: int,
                     target_deg: float):
        """Set the dihedral i-j-k-l to *target_deg* by rotating the
        k-side of bond j-k."""
        current = self.dihedral_angle(i, j, k, l)
        self.rotate_dihedral(j, k, target_deg - current)

    # ---- torsional energy ----

    def torsional_energy(self, j: int, k: int) -> StrainEnergy:
        """Estimate torsional strain about bond j-k (Pitzer potential).

        Sums over every (i on j) x (l on k) substituent pair.
        """
        j_subs = [n for n in self.neighbors(j) if n != k]
        k_subs = [n for n in self.neighbors(k) if n != j]

        total = 0.0
        contributions: list[tuple[TorsionAngle, float]] = []

        for i_idx in j_subs:
            for l_idx in k_subs:
                phi = self.dihedral_angle(i_idx, j, k, l_idx)
                phi_rad = math.radians(phi)

                key = self._torsion_key(i_idx, j, k, l_idx)
                params = TORSION_BARRIERS.get(key, TORSION_BARRIERS["default"])
                V1, V2, V3 = params["V1"], params["V2"], params["V3"]

                energy = (
                    (V1 / 2.0) * (1.0 + math.cos(phi_rad))
                    + (V2 / 2.0) * (1.0 - math.cos(2.0 * phi_rad))
                    + (V3 / 2.0) * (1.0 + math.cos(3.0 * phi_rad))
                )
                torsion = TorsionAngle(i_idx, j, k, l_idx, phi)
                contributions.append((torsion, energy))
                total += energy

        return StrainEnergy(total_kj_per_mol=total,
                            contributions=contributions)

    def _torsion_key(self, i: int, j: int, k: int, l: int) -> str:
        sym_i = self.atoms[i].symbol
        sym_l = self.atoms[l].symbol

        def _hyb_str(idx):
            h = self.atoms[idx].hybridization
            if h == Hybridization.SP3:
                return "sp3"
            if h == Hybridization.SP2:
                return "sp2"
            return "sp"

        a, b = sym_i, sym_l
        ha, hb = _hyb_str(j), _hyb_str(k)
        if a > b:
            a, b = b, a
            ha, hb = hb, ha
        return f"{a}_{ha}_{hb}_{b}"

    # ---- Newman projection ----

    def newman_projection(self, j: int, k: int) -> NewmanProjection:
        """Extract Newman projection data looking from j toward k.

        j is the front atom; k is the back atom.  Substituent angles
        are measured clockwise from 12-o'clock in the projection plane.
        """
        pos_j = self.atoms[j].position
        pos_k = self.atoms[k].position
        view = _normalize(pos_k - pos_j)

        # Orthonormal basis for the projection plane
        up = np.array([0.0, 1.0, 0.0])
        if abs(np.dot(view, up)) > 0.9:
            up = np.array([1.0, 0.0, 0.0])
        right = _normalize(np.cross(view, up))
        proj_up = _normalize(np.cross(right, view))

        def _proj_angle(sub_pos, center_pos):
            v = sub_pos - center_pos
            v_proj = v - np.dot(v, view) * view
            x = float(np.dot(v_proj, right))
            y = float(np.dot(v_proj, proj_up))
            ang = math.degrees(math.atan2(x, y))
            return ang % 360.0

        front = []
        for n in self.neighbors(j):
            if n == k:
                continue
            ang = _proj_angle(self.atoms[n].position, pos_j)
            front.append((n, self.atoms[n].symbol, ang))
        front.sort(key=lambda t: t[2])

        back = []
        for n in self.neighbors(k):
            if n == j:
                continue
            ang = _proj_angle(self.atoms[n].position, pos_k)
            back.append((n, self.atoms[n].symbol, ang))
        back.sort(key=lambda t: t[2])

        dih = 0.0
        if front and back:
            dih = self.dihedral_angle(front[0][0], j, k, back[0][0])

        return NewmanProjection(
            front_atom=j, back_atom=k,
            front_substituents=front,
            back_substituents=back,
            dihedral_deg=dih,
        )

    # ---- stereochemistry ----

    def cip_priority(self, center: int, nbrs: list[int],
                     depth: int = 4) -> list[int]:
        """Rank substituents by CIP (Cahn-Ingold-Prelog) priority.

        Returns neighbour indices from highest to lowest priority.
        Uses BFS expansion with phantom-atom duplication for multiple
        bonds.
        """
        def _expand(start, exclude, d):
            levels = []
            layer = [start]
            visited = {exclude, start}
            for _ in range(d):
                z_list: list[int] = []
                nxt: list[int] = []
                for node in layer:
                    z_list.append(SYMBOL_TO_Z.get(self.atoms[node].symbol, 0))
                    for nb in self.neighbors(node):
                        if nb not in visited:
                            visited.add(nb)
                            nxt.append(nb)
                            bond = self.get_bond(node, nb)
                            if bond and bond.order > 1:
                                for _ in range(bond.order - 1):
                                    z_list.append(
                                        SYMBOL_TO_Z.get(
                                            self.atoms[nb].symbol, 0))
                z_list.sort(reverse=True)
                levels.append(tuple(z_list))
                layer = nxt
                if not layer:
                    break
            return levels

        pmap = {nb: _expand(nb, center, depth) for nb in nbrs}
        return sorted(nbrs, key=lambda nb: pmap[nb], reverse=True)

    def assign_rs(self, center: int) -> Stereodescriptor:
        """Assign R or S to a tetrahedral stereocenter.

        Returns NONE if the centre does not have exactly four different
        substituents.
        """
        nbrs = self.neighbors(center)
        if len(nbrs) != 4:
            return Stereodescriptor.NONE

        ranked = self.cip_priority(center, nbrs)
        # ranked[0]=highest ... ranked[3]=lowest

        c = self.atoms[center].position
        p1 = self.atoms[ranked[0]].position - c
        p2 = self.atoms[ranked[1]].position - c
        p3 = self.atoms[ranked[2]].position - c
        p4 = self.atoms[ranked[3]].position - c

        normal = np.cross(p1 - p2, p3 - p2)
        if np.dot(normal, p4) > 0:
            return Stereodescriptor.R
        return Stereodescriptor.S

    def assign_ez(self, j: int, k: int) -> Stereodescriptor:
        """Assign E or Z about a double bond j=k.

        Uses the dihedral angle between the highest-priority
        substituents on each side.  If they are on the same side
        (|dihedral| < 90) the configuration is Z (zusammen); if on
        opposite sides (|dihedral| > 90) it is E (entgegen).

        Returns NONE if either side has fewer than two substituents.
        """
        j_subs = [n for n in self.neighbors(j) if n != k]
        k_subs = [n for n in self.neighbors(k) if n != j]

        if len(j_subs) < 2 or len(k_subs) < 2:
            return Stereodescriptor.NONE

        j_ranked = self.cip_priority(j, j_subs)
        k_ranked = self.cip_priority(k, k_subs)

        high_j = j_ranked[0]
        high_k = k_ranked[0]

        dih = self.dihedral_angle(high_j, j, k, high_k)

        if abs(dih) < 90.0:
            return Stereodescriptor.Z
        return Stereodescriptor.E

    # ---- visualisation compatibility ----

    def to_coordinates_dict(self) -> dict:
        """Return coordinates in the format used by molecule_visualization.

        Produces the same dict structure as
        vsepr_model.generate_3d_coordinates().
        """
        atom_positions = [
            (a.symbol, a.position.copy()) for a in self.atoms
        ]
        bonds = [(b.atom_i, b.atom_j, b.order) for b in self.bonds]
        return {
            "atom_positions": atom_positions,
            "bonds": bonds,
            "lone_pair_positions": [],
            "central_index": 0,
        }

    # ---- display ----

    def __repr__(self):
        return (f"Molecule({self.name!r}, "
                f"{len(self.atoms)} atoms, {len(self.bonds)} bonds)")

    def summary(self) -> str:
        lines = [
            f"{'=' * 60}",
            f"  Molecule: {self.name}",
            f"  Atoms: {len(self.atoms)}    Bonds: {len(self.bonds)}",
            f"{'=' * 60}",
            f"  Coordinates (Angstroms):",
            f"  {'Idx':<5} {'Sym':<4} {'Hyb':<5}"
            f"  {'x':>8}  {'y':>8}  {'z':>8}",
            f"  {'-' * 48}",
        ]
        for atom in self.atoms:
            hyb = atom.hybridization.name if atom.hybridization else "---"
            x, y, z = atom.position
            lines.append(
                f"  {atom.index:<5} {atom.symbol:<4} {hyb:<5}"
                f"  {x:>8.4f}  {y:>8.4f}  {z:>8.4f}")
        lines.append(f"  {'-' * 48}")
        lines.append("  Bonds:")
        for bond in self.bonds:
            d = self.distance(bond.atom_i, bond.atom_j)
            rot = "rotatable" if bond.rotatable else "fixed"
            sa = self.atoms[bond.atom_i].symbol
            sb = self.atoms[bond.atom_j].symbol
            sym = {1: "-", 2: "=", 3: "#"}.get(bond.order, "?")
            lines.append(
                f"    {sa}[{bond.atom_i}]{sym}"
                f"{sb}[{bond.atom_j}]"
                f"  {d:.3f} A  ({rot})")
        lines.append(f"{'=' * 60}")
        return "\n".join(lines)


# ===================================================================
# Molecule builders
# ===================================================================

def _available_tetrahedral_dirs(
        existing_dirs: list[np.ndarray],
        count: int,
) -> list[np.ndarray]:
    """Compute *count* tetrahedral directions that avoid *existing_dirs*.

    Given the unit-direction vectors of bonds already on a tetrahedral
    centre, returns *count* new unit-direction vectors that complete the
    tetrahedral arrangement.

    Parameters
    ----------
    existing_dirs : directions of existing bonds (unit vectors FROM centre
                    TOWARD each bonded neighbour).
    count : how many new directions to generate (max 4 - len(existing)).
    """
    n = len(existing_dirs)
    tet_angle = math.acos(-1.0 / 3.0)  # ~109.47 deg

    if n == 0:
        # Return up to 4 standard tetrahedral directions
        dirs = [
            np.array([1, 1, 1]) / math.sqrt(3),
            np.array([1, -1, -1]) / math.sqrt(3),
            np.array([-1, 1, -1]) / math.sqrt(3),
            np.array([-1, -1, 1]) / math.sqrt(3),
        ]
        return dirs[:count]

    if n == 1:
        # Three remaining positions: cone at 109.47 deg from v0,
        # spaced 120 deg apart
        v0 = _normalize(existing_dirs[0])
        perp = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(v0, perp)) > 0.9:
            perp = np.array([0.0, 1.0, 0.0])
        p = _normalize(np.cross(v0, perp))
        q = _normalize(np.cross(v0, p))

        ct = math.cos(tet_angle)  # -1/3
        st = math.sin(tet_angle)  # sqrt(8/9)

        dirs = []
        for i in range(count):
            phi = math.radians(120.0 * i)
            d = v0 * ct + (p * math.cos(phi) + q * math.sin(phi)) * st
            dirs.append(_normalize(d))
        return dirs

    if n == 2:
        # Two remaining positions: symmetric about the bisector plane
        v0 = _normalize(existing_dirs[0])
        v1 = _normalize(existing_dirs[1])
        bisector = _normalize(v0 + v1)
        perp_in_plane = _normalize(v0 - v1)
        out_of_plane = _normalize(np.cross(v0, v1))

        # The two new bonds are reflections across the v0-v1 plane,
        # at the tetrahedral angle from both v0 and v1.
        # They lie along: -bisector * cos(alpha) +/- out_of_plane * sin(alpha)
        # where alpha satisfies the tetrahedral constraint.
        # For ideal tetrahedral with 2 existing bonds at 109.47:
        ct = -1.0 / 3.0
        half_angle = math.acos(np.clip(np.dot(v0, v1), -1, 1)) / 2.0
        cos_half = math.cos(half_angle)
        # Component along -bisector
        a = (ct - cos_half ** 2 * ct) / (1.0 - cos_half ** 2) \
            if abs(1.0 - cos_half ** 2) > 1e-10 else ct
        # Actually, let me use a simpler approach: rotate a standard
        # tetrahedral set to align two vertices with v0 and v1.
        # For robustness, compute directly.

        # New directions in the plane perpendicular to bisector:
        # d_new = -bisector * cb + out_of_plane * (+-sb)
        # where cb and sb satisfy: dot(d_new, v0) = -1/3
        # dot(-bisector*cb + out*sb, v0) = -cb*dot(bisector,v0) + sb*dot(out,v0)
        # dot(bisector, v0) = dot(normalize(v0+v1), v0) = (1 + dot(v0,v1)) / |v0+v1|
        bv0 = float(np.dot(bisector, v0))
        ov0 = float(np.dot(out_of_plane, v0))  # should be ~0
        # -cb * bv0 + sb * ov0 = -1/3
        # cb^2 + sb^2 = 1
        # Since ov0 ≈ 0: cb ≈ (1/3) / bv0, sb = sqrt(1 - cb^2)
        if abs(bv0) > 1e-10:
            cb = (1.0 / 3.0) / bv0
        else:
            cb = 0.0
        cb = np.clip(cb, -1.0, 1.0)
        sb = math.sqrt(max(0.0, 1.0 - cb * cb))

        dirs = []
        if count >= 1:
            dirs.append(_normalize(-bisector * cb + out_of_plane * sb))
        if count >= 2:
            dirs.append(_normalize(-bisector * cb - out_of_plane * sb))
        return dirs

    if n >= 3:
        # One remaining position: opposite to centroid of existing
        centroid = _normalize(sum(existing_dirs))
        return [_normalize(-centroid)][:count]

    return []


def _add_sp3_hydrogens(mol: Molecule, carbon_idx: int, count: int):
    """Add *count* hydrogens to fill remaining tetrahedral positions.

    Uses direct 3D geometry to find available tetrahedral directions
    relative to the existing bonds on the carbon.
    """
    c_pos = mol.atoms[carbon_idx].position
    existing = mol.neighbors(carbon_idx)
    if not existing:
        return

    existing_dirs = [
        _normalize(mol.atoms[n].position - c_pos) for n in existing
    ]

    CH = _bond_length("C", "H", 1)
    new_dirs = _available_tetrahedral_dirs(existing_dirs, count)

    for d in new_dirs:
        h_pos = c_pos + CH * d
        h_idx = mol.add_atom("H", h_pos)
        mol.add_bond(carbon_idx, h_idx, order=1, rotatable=False)


def build_ethane(dihedral_deg: float = 60.0) -> Molecule:
    """Build ethane (C2H6) with a specified H-C-C-H dihedral.

    60 deg = staggered (default, lowest energy).
    0 deg = eclipsed (highest energy).
    """
    mol = Molecule("ethane")
    CC = _bond_length("C", "C", 1)
    CH = _bond_length("C", "H", 1)

    c0 = mol.add_atom("C", [0.0, 0.0, 0.0], Hybridization.SP3)
    c1 = mol.add_atom("C", [0.0, 0.0, CC], Hybridization.SP3)
    mol.add_bond(c0, c1)

    # Place 3 H on C0 at tetrahedral positions pointing away from C1
    c0_bond_dir = _normalize(mol.atoms[c1].position - mol.atoms[c0].position)
    c0_h_dirs = _available_tetrahedral_dirs([c0_bond_dir], 3)
    for d in c0_h_dirs:
        h_pos = mol.atoms[c0].position + CH * d
        h_idx = mol.add_atom("H", h_pos)
        mol.add_bond(c0, h_idx, order=1, rotatable=False)

    # Place 3 H on C1 at tetrahedral positions pointing away from C0,
    # then rotate by dihedral_deg
    c1_bond_dir = _normalize(mol.atoms[c0].position - mol.atoms[c1].position)
    c1_h_dirs = _available_tetrahedral_dirs([c1_bond_dir], 3)
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
    CC = _bond_length("C", "C", 1)

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
    _add_sp3_hydrogens(mol, c0, 3)
    _add_sp3_hydrogens(mol, c1, 2)
    _add_sp3_hydrogens(mol, c2, 2)
    _add_sp3_hydrogens(mol, c3, 3)

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
    CC = _bond_length("C", "C", 1)
    CH = _bond_length("C", "H", 1)

    # Ring radius so that adjacent C-C distance = CC
    # Adjacent carbons are separated by 60 deg in angle and 2*d in z:
    #   CC^2 = r^2 + (2d)^2  =>  r = sqrt(CC^2 - 4*d^2)
    d = 0.253  # puckering amplitude (Angstroms)
    r = math.sqrt(CC ** 2 - (2.0 * d) ** 2)

    if conformation == RingConformation.CHAIR:
        puckering = [d, -d, d, -d, d, -d]
    elif conformation == RingConformation.BOAT:
        puckering = [d, -d * 0.5, -d, d, -d, -d * 0.5]
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
        radial = _normalize(c_pos - center)

        # Axial direction alternates with puckering
        if puckering[c_idx] > 0:
            ax_dir = up
        else:
            ax_dir = -up

        eq_dir = _normalize(radial + ax_dir * 0.2)

        # Axial hydrogen
        h_ax_pos = c_pos + CH * _normalize(ax_dir)
        h_ax = mol.add_atom("H", h_ax_pos)
        mol.add_bond(c_idx, h_ax, order=1, rotatable=False)

        # Equatorial hydrogen
        h_eq_pos = c_pos + CH * _normalize(eq_dir)
        h_eq = mol.add_atom("H", h_eq_pos)
        mol.add_bond(c_idx, h_eq, order=1, rotatable=False)

    return mol


def build_2_butene(is_cis: bool = True) -> Molecule:
    """Build 2-butene (CH3-CH=CH-CH3) as cis (Z) or trans (E)."""
    label = "Z/cis" if is_cis else "E/trans"
    mol = Molecule(f"2-butene ({label})")
    CC_s = _bond_length("C", "C", 1)
    CC_d = _bond_length("C", "C", 2)

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
    _add_sp3_hydrogens(mol, c0, 3)
    _add_sp3_hydrogens(mol, c3, 3)

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
        bl = _bond_length("C", sym, 1)
        pos = tet_dirs[i] * bl
        idx = mol.add_atom(sym, pos)
        mol.add_bond(c, idx, order=1)

    return mol


# ===================================================================
# Conformational analysis utilities
# ===================================================================

def classify_conformation(dihedral_deg: float) -> ConformationType:
    """Classify a dihedral angle into a named conformation.

    Ranges (normalised to -180..180):
        |d| < 10          -> ECLIPSED  (syn-periplanar)
        |d| ~ 60          -> GAUCHE
        |d| ~ 120         -> ECLIPSED  (anti-clinal)
        |d| > 170         -> ANTI      (antiperiplanar)
        otherwise         -> CUSTOM
    """
    d = dihedral_deg % 360
    if d > 180:
        d -= 360

    if abs(d) < 10:
        return ConformationType.ECLIPSED
    if abs(abs(d) - 60) < 15:
        return ConformationType.GAUCHE
    if abs(abs(d) - 180) < 15:
        return ConformationType.ANTI
    if abs(abs(d) - 120) < 15:
        return ConformationType.ECLIPSED
    return ConformationType.CUSTOM


def scan_torsion(mol: Molecule, j: int, k: int,
                 ref_i: int, ref_l: int,
                 steps: int = 36) -> list[tuple[float, float]]:
    """Scan torsional energy as a function of dihedral angle.

    Rotates the k-side of bond j-k through 360 degrees, computing
    strain at each step.  Restores original geometry between each
    step to avoid cumulative numerical drift.

    Returns list of (angle_deg, energy_kJ_per_mol).
    """
    originals = [a.position.copy() for a in mol.atoms]

    results: list[tuple[float, float]] = []
    step_size = 360.0 / steps

    for s in range(steps):
        # Reset to original before each step
        for i, pos in enumerate(originals):
            mol.atoms[i].position = pos.copy()

        target = -180.0 + s * step_size
        mol.set_dihedral(ref_i, j, k, ref_l, target)
        energy = mol.torsional_energy(j, k)
        results.append((target, energy.total_kj_per_mol))

    # Restore original
    for i, pos in enumerate(originals):
        mol.atoms[i].position = pos

    return results


# ===================================================================
# Main -- demonstration
# ===================================================================

if __name__ == "__main__":
    print("Molecular Conformations and Stereochemistry\n")

    # ---- 1. Ethane conformations ----------------------------------
    print("--- Ethane: Staggered vs Eclipsed ---\n")
    staggered = build_ethane(dihedral_deg=60.0)
    eclipsed = build_ethane(dihedral_deg=0.0)

    print(staggered.summary())
    e_s = staggered.torsional_energy(0, 1)
    print(e_s.summary())
    print()

    print(eclipsed.summary())
    e_e = eclipsed.torsional_energy(0, 1)
    print(e_e.summary())
    print()

    # Newman projections
    print("Newman (staggered):")
    print(staggered.newman_projection(0, 1).summary())
    print()
    print("Newman (eclipsed):")
    print(eclipsed.newman_projection(0, 1).summary())
    print()

    # ---- 2. Butane conformations ----------------------------------
    print("--- Butane: Anti / Gauche / Eclipsed ---\n")
    print(f"  {'Conf':<14} {'Dihedral':>8}  {'Energy':>12}  {'Type'}")
    print(f"  {'-' * 50}")
    for dih in [180, 60, 0, 120]:
        mol = build_butane(dih)
        energy = mol.torsional_energy(1, 2)
        conf = classify_conformation(dih)
        print(f"  {mol.name:<14} {dih:>6} deg"
              f"  {energy.total_kj_per_mol:>8.2f} kJ/mol"
              f"  {conf.name}")
    print()

    # ---- 3. Torsion energy scan -----------------------------------
    print("--- Butane: Torsion Energy Scan (C0-C1-C2-C3) ---\n")
    butane = build_butane(180.0)
    scan = scan_torsion(butane, 1, 2, 0, 3, steps=36)
    print(f"  {'Angle':>7}  {'kJ/mol':>8}  Bar")
    print(f"  {'-' * 40}")
    for angle, energy in scan:
        bar = "#" * max(0, int(energy / 2))
        print(f"  {angle:>7.1f}  {energy:>8.2f}  {bar}")
    print()

    # ---- 4. Cyclohexane -------------------------------------------
    print("--- Cyclohexane: Chair vs Boat ---\n")
    chair = build_cyclohexane(RingConformation.CHAIR)
    boat = build_cyclohexane(RingConformation.BOAT)
    print(chair.summary())
    print()
    print(boat.summary())
    print()

    # Verify C-C distances in chair
    print("  Chair C-C distances:")
    for i in range(6):
        j = (i + 1) % 6
        print(f"    C{i}-C{j}: {chair.distance(i, j):.3f} A")
    print()

    # ---- 5. E/Z Isomerism ----------------------------------------
    print("--- 2-Butene: E/Z Isomers ---\n")
    cis = build_2_butene(is_cis=True)
    trans = build_2_butene(is_cis=False)
    ez_cis = cis.assign_ez(1, 2)
    ez_trans = trans.assign_ez(1, 2)
    print(f"  cis-2-butene:   {ez_cis.name}")
    print(f"  trans-2-butene: {ez_trans.name}")
    print()

    # ---- 6. R/S Chirality ----------------------------------------
    print("--- Chirality: R/S Assignment ---\n")
    chiral = build_chiral_molecule(["H", "F", "Cl", "Br"])
    descriptor = chiral.assign_rs(0)
    print(f"  CHFClBr configuration: {descriptor.name}")
    print(chiral.summary())
    print()

    print("Done.")
