"""Distance bounds matrix for distance geometry embedding.

Builds an N x N lower/upper bounds matrix from molecular topology:
    - 1-2 pairs (bonded): bond_length +/- 0.01 A
    - 1-3 pairs (angle): law-of-cosines from bond lengths + ideal angles
    - 1-4+ pairs: lower = sum VDW radii * 0.5; upper = sum of path bond lengths

Triangle inequality smoothing tightens bounds iteratively.
"""

from __future__ import annotations

import math
from collections import deque
from typing import TYPE_CHECKING

import numpy as np

from molbuilder.core.bond_data import bond_length as ref_bond_length
from molbuilder.dynamics.forcefield import UFF_LJ
from molbuilder.molecule.graph import Hybridization

if TYPE_CHECKING:
    from molbuilder.molecule.graph import Molecule


def _ideal_angle_rad(hyb: Hybridization | None) -> float:
    """Return the ideal bond angle in radians for a hybridization."""
    if hyb == Hybridization.SP:
        return math.radians(180.0)
    if hyb == Hybridization.SP2:
        return math.radians(120.0)
    return math.radians(109.47)


def _vdw_radius(symbol: str) -> float:
    """VDW radius in Angstroms (half of UFF sigma, floored at 1.2)."""
    sigma = UFF_LJ.get(symbol, (3.5, 0.3))[0]
    return max(sigma * 0.5, 1.2)


def build_bounds(mol: Molecule) -> tuple[np.ndarray, np.ndarray]:
    """Build lower and upper distance bounds matrices.

    Parameters
    ----------
    mol : Molecule
        Molecule with atoms and bonds defined (positions not needed).

    Returns
    -------
    lower, upper : ndarray of shape (N, N)
        Lower and upper distance bounds in Angstroms.
    """
    n = len(mol.atoms)
    lower = np.zeros((n, n))
    upper = np.full((n, n), 1e6)

    # Diagonal
    for i in range(n):
        upper[i, i] = 0.0

    # Build adjacency with bond orders
    adj: dict[int, list[tuple[int, int]]] = {i: [] for i in range(n)}
    bond_map: dict[tuple[int, int], float] = {}
    for b in mol.bonds:
        bl = ref_bond_length(mol.atoms[b.atom_i].symbol,
                             mol.atoms[b.atom_j].symbol, b.order)
        adj[b.atom_i].append((b.atom_j, b.order))
        adj[b.atom_j].append((b.atom_i, b.order))
        bond_map[(b.atom_i, b.atom_j)] = bl
        bond_map[(b.atom_j, b.atom_i)] = bl

    # 1-2 pairs (bonded): tight bounds around ideal bond length
    for b in mol.bonds:
        i, j = b.atom_i, b.atom_j
        bl = bond_map[(i, j)]
        lower[i, j] = lower[j, i] = bl - 0.01
        upper[i, j] = upper[j, i] = bl + 0.01

    # 1-3 pairs (angle): law of cosines
    for j in range(n):
        nbrs = [nb for nb, _ in adj[j]]
        hyb = mol.atoms[j].hybridization
        theta = _ideal_angle_rad(hyb)
        for ii in range(len(nbrs)):
            for kk in range(ii + 1, len(nbrs)):
                i_idx, k_idx = nbrs[ii], nbrs[kk]
                r_ij = bond_map.get((i_idx, j), 1.5)
                r_jk = bond_map.get((j, k_idx), 1.5)
                # Law of cosines: d_ik^2 = r_ij^2 + r_jk^2 - 2*r_ij*r_jk*cos(theta)
                d_sq = r_ij**2 + r_jk**2 - 2.0 * r_ij * r_jk * math.cos(theta)
                d = math.sqrt(max(d_sq, 0.01))
                tol = 0.1
                new_lower = d - tol
                new_upper = d + tol
                if new_lower > lower[i_idx, k_idx]:
                    lower[i_idx, k_idx] = lower[k_idx, i_idx] = new_lower
                if new_upper < upper[i_idx, k_idx]:
                    upper[i_idx, k_idx] = upper[k_idx, i_idx] = new_upper

    # 1-4+ pairs: VDW lower bound, shortest path upper bound
    # Compute shortest path distances via BFS from each atom
    for start in range(n):
        visited = {start: 0.0}
        queue = deque([(start, 0.0)])
        while queue:
            node, dist = queue.popleft()
            for nb, _ in adj[node]:
                new_dist = dist + bond_map.get((node, nb), 1.5)
                if nb not in visited or new_dist < visited[nb]:
                    visited[nb] = new_dist
                    queue.append((nb, new_dist))
        for end, path_dist in visited.items():
            if end <= start:
                continue
            # Only set if not already constrained by 1-2 or 1-3
            vdw_low = (_vdw_radius(mol.atoms[start].symbol)
                       + _vdw_radius(mol.atoms[end].symbol)) * 0.5
            if vdw_low > lower[start, end]:
                lower[start, end] = lower[end, start] = vdw_low
            if path_dist < upper[start, end]:
                upper[start, end] = upper[end, start] = path_dist

    # Ensure lower <= upper
    mask = lower > upper
    avg = (lower + upper) / 2.0
    lower[mask] = avg[mask] - 0.01
    upper[mask] = avg[mask] + 0.01

    # Triangle inequality smoothing (limited iterations for speed)
    _triangle_smooth(lower, upper, n)

    return lower, upper


def _triangle_smooth(lower: np.ndarray, upper: np.ndarray, n: int,
                     max_iters: int = 3) -> None:
    """Apply triangle inequality smoothing in-place.

    For each triple (i, j, k):
        upper[i,k] <= upper[i,j] + upper[j,k]
        lower[i,k] >= lower[i,j] - upper[j,k]
    """
    for _ in range(max_iters):
        changed = False
        for k in range(n):
            for i in range(n):
                if i == k:
                    continue
                for j in range(i + 1, n):
                    if j == k:
                        continue
                    # Tighten upper bound
                    u_new = upper[i, k] + upper[k, j]
                    if u_new < upper[i, j]:
                        upper[i, j] = upper[j, i] = u_new
                        changed = True
                    # Tighten lower bound
                    l_new = lower[i, k] - upper[k, j]
                    if l_new > lower[i, j]:
                        lower[i, j] = lower[j, i] = l_new
                        changed = True
                    l_new2 = lower[j, k] - upper[k, i]
                    if l_new2 > lower[i, j]:
                        lower[i, j] = lower[j, i] = l_new2
                        changed = True
        if not changed:
            break
