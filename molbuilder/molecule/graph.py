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

from molbuilder.core.elements import SYMBOL_TO_Z
from molbuilder.core.element_properties import cpk_color
from molbuilder.core.bond_data import bond_length, SP3_ANGLE, SP2_ANGLE, SP_ANGLE, TORSION_BARRIERS
from molbuilder.core.geometry import normalize, rotation_matrix, place_atom_zmatrix


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
                         bond_length_val: float | None = None,
                         bond_angle_deg: float | None = None,
                         dihedral_deg: float = 0.0,
                         hybridization: Hybridization | None = None,
                         rotatable: bool | None = None,
                         **kwargs) -> int:
        """Add an atom bonded to an existing atom via internal coordinates.

        For the first atom (index 0): placed at origin.
        For the second atom (index 1): placed along +z.
        For the third: placed in the xz-plane using a synthetic dihedral
        reference.  For all subsequent atoms: full z-matrix placement.

        Returns the index of the new atom.
        """
        # Support both 'bond_length' and 'bond_length_val' parameter names
        bl = bond_length_val if bond_length_val is not None else kwargs.get('bond_length', None)
        if bl is None:
            parent_sym = self.atoms[bonded_to].symbol
            bl = bond_length(parent_sym, symbol, bond_order)

        n = len(self.atoms)

        # --- first atom ------------------------------------------------
        if n == 0:
            return self.add_atom(symbol, [0.0, 0.0, 0.0], hybridization)

        # --- second atom -----------------------------------------------
        if n == 1:
            pos = (self.atoms[bonded_to].position
                   + np.array([0.0, 0.0, bl]))
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
            pos = place_atom_zmatrix(pos_j, pos_i, synthetic_k,
                                       bl, bond_angle_deg,
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

        pos = place_atom_zmatrix(
            self.atoms[bonded_to].position,
            self.atoms[angle_ref].position,
            self.atoms[dihedral_ref].position,
            bl, bond_angle_deg, dihedral_deg,
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
        R = rotation_matrix(axis, math.radians(angle_deg))

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
        view = normalize(pos_k - pos_j)

        # Orthonormal basis for the projection plane
        up = np.array([0.0, 1.0, 0.0])
        if abs(np.dot(view, up)) > 0.9:
            up = np.array([1.0, 0.0, 0.0])
        right = normalize(np.cross(view, up))
        proj_up = normalize(np.cross(right, view))

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
