"""Geometry optimization via L-BFGS-B using the existing ForceField.

Builds a ForceField with ideal equilibrium values (not measured from current
positions) and minimizes energy to refine DG-embedded coordinates.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import minimize

from molbuilder.core.bond_data import bond_length as ref_bond_length
from molbuilder.dynamics.forcefield import (
    ForceField, ForceFieldParams, UFF_LJ,
    _estimate_partial_charges, _estimate_bond_k, _estimate_angle_k,
)
from molbuilder.core.elements import ELEMENTS, SYMBOL_TO_Z
from molbuilder.molecule.graph import Hybridization

if TYPE_CHECKING:
    from molbuilder.molecule.graph import Molecule


# Ideal angles by hybridization in radians
_IDEAL_ANGLES = {
    Hybridization.SP3: math.radians(109.47),
    Hybridization.SP2: math.radians(120.0),
    Hybridization.SP: math.radians(180.0),
}

# Aromatic bond length for C-C (Angstroms)
_AROMATIC_CC = 1.40


def optimize_geometry(mol: Molecule, max_iter: int = 200) -> None:
    """Optimize molecular geometry using L-BFGS-B with ideal equilibrium values.

    Modifies atom positions in-place.

    Parameters
    ----------
    mol : Molecule
        Molecule with DG-embedded coordinates.
    max_iter : int
        Maximum L-BFGS-B iterations.
    """
    n = len(mol.atoms)
    if n < 2:
        return

    ff = _build_ideal_forcefield(mol)
    sp2_centers = _find_sp2_centers(mol)
    positions = np.array([a.position for a in mol.atoms])

    def objective(x):
        pos = x.reshape(n, 3)
        result = ff.compute(pos)
        energy = result.energy_total
        grad = -result.forces

        # Add out-of-plane penalty for SP2 centers
        e_oop, g_oop = _out_of_plane_penalty(pos, sp2_centers, mol)
        energy += e_oop
        grad += g_oop

        return energy, grad.ravel()

    result = minimize(
        objective, positions.ravel(), method='L-BFGS-B', jac=True,
        options={'maxiter': max_iter, 'ftol': 1e-4, 'gtol': 1e-3},
    )

    new_pos = result.x.reshape(n, 3)
    for i, atom in enumerate(mol.atoms):
        atom.position = new_pos[i]


def _find_sp2_centers(mol: Molecule) -> list[tuple[int, list[int]]]:
    """Find SP2 atoms with 3+ neighbors (candidates for planarity enforcement)."""
    centers = []
    for atom in mol.atoms:
        if atom.hybridization == Hybridization.SP2:
            nbrs = mol.neighbors(atom.index)
            if len(nbrs) >= 3:
                centers.append((atom.index, nbrs[:3]))
    return centers


def _out_of_plane_penalty(
    pos: np.ndarray,
    sp2_centers: list[tuple[int, list[int]]],
    mol: Molecule,
) -> tuple[float, np.ndarray]:
    """Compute out-of-plane penalty for SP2 centers.

    For each SP2 center with 3 neighbors, measures the distance of the center
    from the plane defined by its neighbors and penalizes deviation.
    """
    k_oop = 200.0  # kJ/(mol*A^2)
    energy = 0.0
    grad = np.zeros_like(pos)

    for center_idx, nbrs in sp2_centers:
        if len(nbrs) < 3:
            continue
        # Plane defined by center and two neighbors; third should be in-plane
        p0 = pos[center_idx]
        p1 = pos[nbrs[0]]
        p2 = pos[nbrs[1]]
        p3 = pos[nbrs[2]]

        # Normal to the plane formed by (p1, p2, p3)
        v1 = p2 - p1
        v2 = p3 - p1
        normal = np.cross(v1, v2)
        norm_len = np.linalg.norm(normal)
        if norm_len < 1e-10:
            continue
        normal = normal / norm_len

        # Distance of center from the plane
        d = np.dot(p0 - p1, normal)
        energy += 0.5 * k_oop * d * d

        # Gradient on center
        grad[center_idx] += k_oop * d * normal

        # Approximate gradient on plane atoms (distribute equally)
        f_plane = -k_oop * d * normal / 3.0
        grad[nbrs[0]] += f_plane
        grad[nbrs[1]] += f_plane
        grad[nbrs[2]] += f_plane

    return energy, grad


def _is_aromatic_bond(mol: Molecule, b) -> bool:
    """Check if a bond is aromatic (both atoms SP2, non-H, and in a ring)."""
    a_i = mol.atoms[b.atom_i]
    a_j = mol.atoms[b.atom_j]
    if a_i.hybridization != Hybridization.SP2:
        return False
    if a_j.hybridization != Hybridization.SP2:
        return False
    if a_i.symbol == "H" or a_j.symbol == "H":
        return False
    return mol.is_in_ring(b.atom_i, b.atom_j)


def _build_ideal_forcefield(mol: Molecule) -> ForceField:
    """Build ForceField with ideal equilibrium angles from hybridization.

    Unlike ForceField.from_molecule() which reads angle_theta0 from current
    (possibly wrong) geometry, this uses ideal angles from hybridization.

    Also uses aromatic bond lengths (1.40 A) for bonds between SP2 atoms
    in rings, instead of kekulized single/double bond lengths.
    """
    from molbuilder.core.bond_data import TORSION_BARRIERS

    n = len(mol.atoms)
    symbols = [a.symbol for a in mol.atoms]

    # Masses
    masses = np.array([
        ELEMENTS[SYMBOL_TO_Z.get(s, 1)][2] for s in symbols
    ])

    # LJ parameters -- scale down epsilon for geometry optimization.
    # Full LJ repulsion overpowers bonded terms during optimization.
    _LJ_EPS_SCALE = 0.1
    default_lj = (3.5, 0.3)
    sigma = np.array([UFF_LJ.get(s, default_lj)[0] for s in symbols])
    epsilon = np.array([
        UFF_LJ.get(s, default_lj)[1] * _LJ_EPS_SCALE for s in symbols
    ])

    # Bonds -- with aromatic bond length detection
    bond_list = [(b.atom_i, b.atom_j, b.order) for b in mol.bonds]
    bond_indices = np.array(
        [(b.atom_i, b.atom_j) for b in mol.bonds], dtype=int,
    ).reshape(-1, 2)

    # For geometry optimization, bond force constants must be much stiffer
    # than for MD. The BDE-based estimates (~300 kJ/mol/A^2 for C-C) are too
    # weak relative to LJ repulsion. Real force fields (MMFF94, UFF) use
    # ~2000-4000 kJ/mol/A^2. We scale up by 10x.
    _BOND_K_SCALE = 10.0

    bond_r0_list = []
    bond_k_list = []
    for b in mol.bonds:
        si, sj = symbols[b.atom_i], symbols[b.atom_j]
        if _is_aromatic_bond(mol, b) and si == "C" and sj == "C":
            bond_r0_list.append(_AROMATIC_CC)
            bond_k_list.append(_estimate_bond_k(si, sj, 1) * _BOND_K_SCALE * 1.5)
        else:
            bond_r0_list.append(ref_bond_length(si, sj, b.order))
            bond_k_list.append(_estimate_bond_k(si, sj, b.order) * _BOND_K_SCALE)

    bond_r0 = np.array(bond_r0_list)
    bond_k = np.array(bond_k_list)

    # Charges
    charges = _estimate_partial_charges(symbols, bond_list)

    # Adjacency
    adj: dict[int, list[int]] = {i: [] for i in range(n)}
    for b in mol.bonds:
        adj[b.atom_i].append(b.atom_j)
        adj[b.atom_j].append(b.atom_i)

    # Angles with IDEAL theta0 from hybridization
    angles: list[tuple[int, int, int]] = []
    for j in range(n):
        nbrs = adj[j]
        for ii in range(len(nbrs)):
            for kk in range(ii + 1, len(nbrs)):
                angles.append((nbrs[ii], j, nbrs[kk]))

    angle_indices = np.array(angles, dtype=int).reshape(-1, 3)
    angle_theta0 = np.array([
        _IDEAL_ANGLES.get(mol.atoms[j].hybridization, math.radians(109.47))
        for _, j, _ in angles
    ])
    angle_k_base = _estimate_angle_k()

    # Boost angle force constant for SP2 centers (planarity enforcement)
    angle_k_arr = np.array([
        angle_k_base * 3.0 if mol.atoms[j].hybridization == Hybridization.SP2
        else angle_k_base
        for _, j, _ in angles
    ])

    # Torsions
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

    # Exclusions (1-2, 1-3, 1-4)
    exclusions: set[frozenset[int]] = set()
    for b in mol.bonds:
        exclusions.add(frozenset((b.atom_i, b.atom_j)))
    for i_a, j_a, k_a in angles:
        exclusions.add(frozenset((i_a, k_a)))
    for _, j_e, k_e in angles:
        for nb_j in adj[j_e]:
            exclusions.add(frozenset((nb_j, k_e)))
            for nb_k in adj[k_e]:
                exclusions.add(frozenset((nb_j, nb_k)))

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
        angle_k=angle_k_arr,
        torsion_indices=torsion_indices,
        torsion_V=torsion_V,
        exclusion_14=exclusions,
    )
    return ForceField(params)
