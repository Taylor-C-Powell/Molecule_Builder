"""Force field engine for molecular dynamics.

Implements a classical mechanics force field with five energy terms:

    E_total = E_bond + E_angle + E_torsion + E_LJ + E_coulomb

Scientific basis:
    - OPLS-AA (Jorgensen et al., J. Am. Chem. Soc. 1996) for torsion
    - UFF (Rappe et al., J. Am. Chem. Soc. 1992) for Lennard-Jones
    - Harmonic approximation for bond stretching and angle bending
    - Electronegativity-based partial charges (Gasteiger-like)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from molbuilder.core.bond_data import (
    bond_length as ref_bond_length,
    BDE_TABLE,
    TORSION_BARRIERS,
)
from molbuilder.core.element_properties import (
    electronegativity,
    PAULING_ELECTRONEGATIVITY,
)
from molbuilder.core.elements import ELEMENTS, SYMBOL_TO_Z

if TYPE_CHECKING:
    from molbuilder.molecule.graph import Molecule


# ===================================================================
# UFF Lennard-Jones parameters  (sigma in Angstroms, epsilon in kJ/mol)
# Source: Rappe et al., J. Am. Chem. Soc. 1992, 114, 10024-10035
# ===================================================================

UFF_LJ: dict[str, tuple[float, float]] = {
    # symbol: (sigma_A, epsilon_kJ_mol)
    "H":  (2.886, 0.184),
    "He": (2.362, 0.056),
    "Li": (2.451, 0.025),
    "Be": (2.745, 0.085),
    "B":  (4.083, 0.389),
    "C":  (3.851, 0.439),
    "N":  (3.660, 0.447),
    "O":  (3.500, 0.460),
    "F":  (3.364, 0.050),
    "Ne": (3.243, 0.042),
    "Na": (2.983, 0.030),
    "Mg": (3.021, 0.111),
    "Al": (4.499, 0.505),
    "Si": (4.295, 0.402),
    "P":  (4.147, 0.305),
    "S":  (4.035, 1.046),
    "Cl": (3.947, 0.950),
    "Ar": (3.868, 0.849),
    "K":  (3.812, 0.035),
    "Ca": (3.399, 0.119),
    "Br": (4.189, 1.220),
    "I":  (4.500, 1.581),
    "Fe": (2.912, 0.013),
    "Cu": (3.495, 0.005),
    "Zn": (2.763, 0.124),
    "Se": (4.205, 0.291),
}

# Coulomb constant in kJ*A / (mol * e^2)
# k_e = 1389.354578 kJ*A/mol for charges in elementary charge units
_COULOMB_KJ_A = 1389.354578


# ===================================================================
# Helper: estimate partial charges from electronegativity
# ===================================================================

def _estimate_partial_charges(symbols: list[str],
                               bond_pairs: list[tuple[int, int, int]],
                               ) -> np.ndarray:
    """Gasteiger-like partial charge estimation from electronegativity.

    Iteratively equalizes electronegativity across bonds.  A simplified
    single-pass version is used here for speed.

    Parameters
    ----------
    symbols : list[str]
        Element symbols for each atom.
    bond_pairs : list[tuple[int, int, int]]
        (atom_i, atom_j, bond_order) for each bond.

    Returns
    -------
    charges : ndarray of shape (n_atoms,)
        Partial charges in elementary charge units.
    """
    n = len(symbols)
    charges = np.zeros(n)
    en = np.array([electronegativity(s) for s in symbols])

    for i, j, order in bond_pairs:
        if en[i] == 0.0 or en[j] == 0.0:
            continue
        delta = (en[j] - en[i]) / (en[i] + en[j])
        transfer = 0.15 * order * delta
        charges[i] += transfer
        charges[j] -= transfer

    return charges


# ===================================================================
# Force field spring constant estimation
# ===================================================================

def _estimate_bond_k(sym_a: str, sym_b: str, order: int) -> float:
    """Estimate harmonic bond force constant in kJ/(mol*A^2).

    Uses BDE as a rough guide: k ~ 2 * BDE / r0^2 (very approximate).
    Falls back to a generic value if BDE is unavailable.
    """
    a, b = sorted([sym_a, sym_b])
    bde = BDE_TABLE.get((a, b, order))
    r0 = ref_bond_length(sym_a, sym_b, order)
    if bde is not None and r0 > 0:
        return 2.0 * bde / (r0 * r0)
    # Generic fallback: ~500 kJ/(mol*A^2) for single bond
    return 500.0 * order


def _estimate_angle_k() -> float:
    """Generic harmonic angle force constant in kJ/(mol*rad^2).

    Typical values range from 300-600 kJ/(mol*rad^2).
    """
    return 400.0


# ===================================================================
# Data classes
# ===================================================================

@dataclass
class ForceFieldParams:
    """Per-atom and bonded-term parameters for the force field.

    Attributes
    ----------
    n_atoms : int
        Number of atoms.
    masses : ndarray of shape (n_atoms,)
        Atomic masses in AMU.
    symbols : list[str]
        Element symbol per atom.
    sigma : ndarray of shape (n_atoms,)
        LJ sigma in Angstroms.
    epsilon : ndarray of shape (n_atoms,)
        LJ epsilon in kJ/mol.
    charges : ndarray of shape (n_atoms,)
        Partial charges in elementary charge units.
    bond_indices : ndarray of shape (n_bonds, 2)
        Atom index pairs for each bond.
    bond_r0 : ndarray of shape (n_bonds,)
        Equilibrium bond length in Angstroms.
    bond_k : ndarray of shape (n_bonds,)
        Bond force constant in kJ/(mol*A^2).
    angle_indices : ndarray of shape (n_angles, 3)
        Atom index triples (i, j, k) where j is the central atom.
    angle_theta0 : ndarray of shape (n_angles,)
        Equilibrium angle in radians.
    angle_k : ndarray of shape (n_angles,)
        Angle force constant in kJ/(mol*rad^2).
    torsion_indices : ndarray of shape (n_torsions, 4)
        Atom index quads (i, j, k, l).
    torsion_V : ndarray of shape (n_torsions, 3)
        V1, V2, V3 torsion parameters in kJ/mol.
    exclusion_14 : set[frozenset[int]]
        Atom pairs connected by <= 3 bonds (excluded from LJ/Coulomb).
    """
    n_atoms: int
    masses: np.ndarray
    symbols: list[str]
    sigma: np.ndarray
    epsilon: np.ndarray
    charges: np.ndarray
    bond_indices: np.ndarray
    bond_r0: np.ndarray
    bond_k: np.ndarray
    angle_indices: np.ndarray
    angle_theta0: np.ndarray
    angle_k: np.ndarray
    torsion_indices: np.ndarray
    torsion_V: np.ndarray
    exclusion_14: set = field(default_factory=set)


@dataclass
class ForceResult:
    """Result of a force computation.

    Attributes
    ----------
    forces : ndarray of shape (n_atoms, 3)
        Force on each atom in kJ/(mol*A).
    energy_bond : float
    energy_angle : float
    energy_torsion : float
    energy_lj : float
    energy_coulomb : float
    energy_total : float
        All energies in kJ/mol.
    """
    forces: np.ndarray
    energy_bond: float = 0.0
    energy_angle: float = 0.0
    energy_torsion: float = 0.0
    energy_lj: float = 0.0
    energy_coulomb: float = 0.0

    @property
    def energy_total(self) -> float:
        return (self.energy_bond + self.energy_angle + self.energy_torsion
                + self.energy_lj + self.energy_coulomb)


# ===================================================================
# ForceField
# ===================================================================

class ForceField:
    """Vectorized classical force field.

    Computes forces and energies for bond stretching (harmonic),
    angle bending (harmonic), torsional rotation (OPLS-AA cosine),
    Lennard-Jones van der Waals, and Coulomb electrostatics.

    Parameters
    ----------
    params : ForceFieldParams
        Pre-built parameter set.
    """

    def __init__(self, params: ForceFieldParams):
        self.params = params

    @classmethod
    def from_molecule(cls, mol: Molecule) -> ForceField:
        """Auto-parameterize a ForceField from a Molecule instance.

        Reads atom types, bond connectivity, and geometry from the
        Molecule's atoms, bonds, and adjacency information, then builds
        angle and torsion lists by traversing the molecular graph.

        Parameters
        ----------
        mol : Molecule
            A molbuilder Molecule with atoms and bonds defined.

        Returns
        -------
        ForceField
            Parameterized force field ready for ``compute()``.
        """
        n = len(mol.atoms)
        symbols = [a.symbol for a in mol.atoms]

        # Masses
        masses = np.array([
            ELEMENTS[SYMBOL_TO_Z.get(s, 1)][2] for s in symbols
        ])

        # LJ parameters
        default_lj = (3.5, 0.3)
        sigma = np.array([UFF_LJ.get(s, default_lj)[0] for s in symbols])
        epsilon = np.array([UFF_LJ.get(s, default_lj)[1] for s in symbols])

        # Bonds
        bond_list = [(b.atom_i, b.atom_j, b.order) for b in mol.bonds]
        bond_indices = np.array(
            [(b.atom_i, b.atom_j) for b in mol.bonds],
            dtype=int,
        ).reshape(-1, 2)
        bond_r0 = np.array([
            ref_bond_length(symbols[b.atom_i], symbols[b.atom_j], b.order)
            for b in mol.bonds
        ])
        bond_k = np.array([
            _estimate_bond_k(symbols[b.atom_i], symbols[b.atom_j], b.order)
            for b in mol.bonds
        ])

        # Charges
        charges = _estimate_partial_charges(symbols, bond_list)

        # Build adjacency for angle/torsion enumeration
        adj: dict[int, list[int]] = {i: [] for i in range(n)}
        for b in mol.bonds:
            adj[b.atom_i].append(b.atom_j)
            adj[b.atom_j].append(b.atom_i)

        # Angles: all i-j-k where j is the central atom
        angles: list[tuple[int, int, int]] = []
        for j in range(n):
            nbrs = adj[j]
            for ii in range(len(nbrs)):
                for kk in range(ii + 1, len(nbrs)):
                    angles.append((nbrs[ii], j, nbrs[kk]))

        angle_indices = np.array(angles, dtype=int).reshape(-1, 3)
        angle_theta0_list = []
        for i_a, j_a, k_a in angles:
            vi = mol.atoms[i_a].position - mol.atoms[j_a].position
            vk = mol.atoms[k_a].position - mol.atoms[j_a].position
            ni = np.linalg.norm(vi)
            nk = np.linalg.norm(vk)
            if ni > 1e-12 and nk > 1e-12:
                cos_a = np.clip(np.dot(vi, vk) / (ni * nk), -1.0, 1.0)
                angle_theta0_list.append(math.acos(cos_a))
            else:
                angle_theta0_list.append(math.radians(109.47))
        angle_theta0 = np.array(angle_theta0_list)
        angle_k = np.full(len(angles), _estimate_angle_k())

        # Torsions: all i-j-k-l where j-k is a bond
        torsions: list[tuple[int, int, int, int]] = []
        torsion_params: list[tuple[float, float, float]] = []

        def _hyb_str(idx: int) -> str:
            h = mol.atoms[idx].hybridization
            if h is not None:
                return h.name.lower()
            return "sp3"

        for b in mol.bonds:
            j, k = b.atom_i, b.atom_j
            j_nbrs = [x for x in adj[j] if x != k]
            k_nbrs = [x for x in adj[k] if x != j]
            for i_t in j_nbrs:
                for l_t in k_nbrs:
                    torsions.append((i_t, j, k, l_t))
                    # Look up torsion parameters
                    si = symbols[i_t]
                    sl = symbols[l_t]
                    a, b_s = si, sl
                    ha, hb = _hyb_str(j), _hyb_str(k)
                    if a > b_s:
                        a, b_s = b_s, a
                        ha, hb = hb, ha
                    key = f"{a}_{ha}_{hb}_{b_s}"
                    p = TORSION_BARRIERS.get(key, TORSION_BARRIERS["default"])
                    torsion_params.append((p["V1"], p["V2"], p["V3"]))

        torsion_indices = np.array(torsions, dtype=int).reshape(-1, 4)
        torsion_V = np.array(torsion_params).reshape(-1, 3)

        # Build 1-2, 1-3, 1-4 exclusion set
        exclusions: set[frozenset[int]] = set()
        for b in mol.bonds:
            exclusions.add(frozenset((b.atom_i, b.atom_j)))
        for _, j_e, k_e in angles:
            for nb_j in adj[j_e]:
                exclusions.add(frozenset((nb_j, k_e)))
                for nb_k in adj[k_e]:
                    exclusions.add(frozenset((nb_j, nb_k)))
        # Also exclude 1-2 pairs found in angles
        for i_e, j_e, k_e in angles:
            exclusions.add(frozenset((i_e, k_e)))

        params = ForceFieldParams(
            n_atoms=n,
            masses=masses,
            symbols=symbols,
            sigma=sigma,
            epsilon=epsilon,
            charges=charges,
            bond_indices=bond_indices,
            bond_r0=bond_r0,
            bond_k=bond_k,
            angle_indices=angle_indices,
            angle_theta0=angle_theta0,
            angle_k=angle_k,
            torsion_indices=torsion_indices,
            torsion_V=torsion_V,
            exclusion_14=exclusions,
        )
        return cls(params)

    def compute(self, positions: np.ndarray) -> ForceResult:
        """Compute forces and energies at the given atomic positions.

        Parameters
        ----------
        positions : ndarray of shape (n_atoms, 3)
            Atomic positions in Angstroms.

        Returns
        -------
        ForceResult
            Forces in kJ/(mol*A) and energies in kJ/mol.
        """
        p = self.params
        forces = np.zeros_like(positions)

        e_bond = self._compute_bonds(positions, forces)
        e_angle = self._compute_angles(positions, forces)
        e_torsion = self._compute_torsions(positions, forces)
        e_lj, e_coul = self._compute_nonbonded(positions, forces)

        return ForceResult(
            forces=forces,
            energy_bond=e_bond,
            energy_angle=e_angle,
            energy_torsion=e_torsion,
            energy_lj=e_lj,
            energy_coulomb=e_coul,
        )

    def _compute_bonds(self, pos: np.ndarray, forces: np.ndarray) -> float:
        """Harmonic bond stretching: E = 0.5 * k * (r - r0)^2."""
        p = self.params
        if p.bond_indices.size == 0:
            return 0.0

        energy = 0.0
        for b_idx in range(len(p.bond_r0)):
            i, j = p.bond_indices[b_idx]
            rij = pos[j] - pos[i]
            r = np.linalg.norm(rij)
            if r < 1e-12:
                continue
            dr = r - p.bond_r0[b_idx]
            e = 0.5 * p.bond_k[b_idx] * dr * dr
            energy += e
            # Force: -dE/dr * r_hat
            f_mag = -p.bond_k[b_idx] * dr
            f_vec = f_mag * (rij / r)
            forces[i] -= f_vec
            forces[j] += f_vec

        return energy

    def _compute_angles(self, pos: np.ndarray, forces: np.ndarray) -> float:
        """Harmonic angle bending: E = 0.5 * k * (theta - theta0)^2."""
        p = self.params
        if p.angle_indices.size == 0:
            return 0.0

        energy = 0.0
        for a_idx in range(len(p.angle_theta0)):
            i, j, k = p.angle_indices[a_idx]
            rji = pos[i] - pos[j]
            rjk = pos[k] - pos[j]
            nji = np.linalg.norm(rji)
            njk = np.linalg.norm(rjk)
            if nji < 1e-12 or njk < 1e-12:
                continue

            cos_theta = np.clip(np.dot(rji, rjk) / (nji * njk), -1.0, 1.0)
            theta = math.acos(cos_theta)
            d_theta = theta - p.angle_theta0[a_idx]
            e = 0.5 * p.angle_k[a_idx] * d_theta * d_theta
            energy += e

            # Gradient of angle w.r.t. positions
            sin_theta = math.sin(theta)
            if abs(sin_theta) < 1e-12:
                continue

            dE_dtheta = p.angle_k[a_idx] * d_theta

            # Force on atom i
            rji_hat = rji / nji
            rjk_hat = rjk / njk
            fi = (dE_dtheta / (nji * sin_theta)) * (
                cos_theta * rji_hat - rjk_hat)
            fk = (dE_dtheta / (njk * sin_theta)) * (
                cos_theta * rjk_hat - rji_hat)

            forces[i] -= fi
            forces[k] -= fk
            forces[j] += fi + fk  # Newton's third law

        return energy

    def _compute_torsions(self, pos: np.ndarray,
                           forces: np.ndarray) -> float:
        """OPLS-AA torsional: E = sum V_n/2 * (1 + cos(n*phi - delta))."""
        p = self.params
        if p.torsion_indices.size == 0:
            return 0.0

        energy = 0.0
        for t_idx in range(len(p.torsion_V)):
            i, j, k, l = p.torsion_indices[t_idx]
            V1, V2, V3 = p.torsion_V[t_idx]

            b1 = pos[j] - pos[i]
            b2 = pos[k] - pos[j]
            b3 = pos[l] - pos[k]

            n1 = np.cross(b1, b2)
            n2 = np.cross(b2, b3)
            n1_norm = np.linalg.norm(n1)
            n2_norm = np.linalg.norm(n2)
            if n1_norm < 1e-12 or n2_norm < 1e-12:
                continue
            n1 /= n1_norm
            n2 /= n2_norm

            b2_hat = b2 / np.linalg.norm(b2)
            x = float(np.dot(n1, n2))
            y = float(np.dot(np.cross(n1, b2_hat), n2))
            phi = math.atan2(y, x)

            e = ((V1 / 2.0) * (1.0 + math.cos(phi))
                 + (V2 / 2.0) * (1.0 - math.cos(2.0 * phi))
                 + (V3 / 2.0) * (1.0 + math.cos(3.0 * phi)))
            energy += e

            # Numerical gradient for torsional forces (analytical is complex)
            eps = 1e-5
            for atom_idx in (i, j, k, l):
                for dim in range(3):
                    pos_p = pos.copy()
                    pos_p[atom_idx, dim] += eps
                    phi_p = self._dihedral(pos_p, i, j, k, l)
                    e_p = ((V1 / 2.0) * (1.0 + math.cos(phi_p))
                           + (V2 / 2.0) * (1.0 - math.cos(2.0 * phi_p))
                           + (V3 / 2.0) * (1.0 + math.cos(3.0 * phi_p)))

                    pos_m = pos.copy()
                    pos_m[atom_idx, dim] -= eps
                    phi_m = self._dihedral(pos_m, i, j, k, l)
                    e_m = ((V1 / 2.0) * (1.0 + math.cos(phi_m))
                           + (V2 / 2.0) * (1.0 - math.cos(2.0 * phi_m))
                           + (V3 / 2.0) * (1.0 + math.cos(3.0 * phi_m)))

                    forces[atom_idx, dim] -= (e_p - e_m) / (2.0 * eps)

        return energy

    @staticmethod
    def _dihedral(pos: np.ndarray, i: int, j: int, k: int, l: int) -> float:
        """Compute dihedral angle for four atom indices."""
        b1 = pos[j] - pos[i]
        b2 = pos[k] - pos[j]
        b3 = pos[l] - pos[k]
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)
        n1n = np.linalg.norm(n1)
        n2n = np.linalg.norm(n2)
        if n1n < 1e-12 or n2n < 1e-12:
            return 0.0
        n1 /= n1n
        n2 /= n2n
        b2h = b2 / np.linalg.norm(b2)
        x = float(np.dot(n1, n2))
        y = float(np.dot(np.cross(n1, b2h), n2))
        return math.atan2(y, x)

    def _compute_nonbonded(self, pos: np.ndarray,
                            forces: np.ndarray) -> tuple[float, float]:
        """Lennard-Jones and Coulomb non-bonded interactions."""
        p = self.params
        e_lj = 0.0
        e_coul = 0.0

        for i in range(p.n_atoms):
            for j in range(i + 1, p.n_atoms):
                if frozenset((i, j)) in p.exclusion_14:
                    continue

                rij = pos[j] - pos[i]
                r = np.linalg.norm(rij)
                if r < 0.5:  # Prevent singularity
                    r = 0.5

                rij_hat = rij / r

                # Lennard-Jones (Lorentz-Berthelot combining rules)
                sig_ij = 0.5 * (p.sigma[i] + p.sigma[j])
                eps_ij = math.sqrt(p.epsilon[i] * p.epsilon[j])

                if eps_ij > 0:
                    sr6 = (sig_ij / r) ** 6
                    sr12 = sr6 * sr6
                    e_lj += 4.0 * eps_ij * (sr12 - sr6)
                    f_lj = 4.0 * eps_ij * (12.0 * sr12 - 6.0 * sr6) / r
                    f_vec = f_lj * rij_hat
                    forces[i] -= f_vec
                    forces[j] += f_vec

                # Coulomb
                qi, qj = p.charges[i], p.charges[j]
                if abs(qi) > 1e-10 and abs(qj) > 1e-10:
                    e_c = _COULOMB_KJ_A * qi * qj / r
                    e_coul += e_c
                    f_c = _COULOMB_KJ_A * qi * qj / (r * r)
                    f_vec_c = f_c * rij_hat
                    forces[i] -= f_vec_c
                    forces[j] += f_vec_c

        return e_lj, e_coul
