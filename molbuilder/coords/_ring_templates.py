"""Ideal ring coordinate templates and Kabsch alignment.

Provides pre-built ideal coordinates for common ring types and
aligns them onto embedded coordinates via SVD (Kabsch algorithm).
Absorbed from smiles/ring_geometry.py for use in the DG pipeline.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from molbuilder.core.bond_data import bond_length as ref_bond_length
from molbuilder.core.geometry import normalize
from molbuilder.molecule.graph import Hybridization

if TYPE_CHECKING:
    from molbuilder.molecule.graph import Molecule


def apply_ring_templates(mol: Molecule) -> None:
    """Detect rings and apply ideal geometry templates via Kabsch alignment.

    Modifies atom positions in-place. Should be called after DG embedding
    and before FF optimization.
    """
    rings = _find_rings(mol)
    if not rings:
        return

    rings.sort(key=len)
    corrected: set[int] = set()

    for ring in rings:
        if all(idx in corrected for idx in ring):
            continue

        old_positions = {i: mol.atoms[i].position.copy() for i in ring}
        ideal = _ideal_coords(mol, ring)
        shared = {idx for idx in ring if idx in corrected}

        if len(shared) >= 2:
            aligned = _align_fused(ideal, ring, shared, mol)
        else:
            target = np.array([mol.atoms[i].position for i in ring])
            aligned = kabsch_align(ideal, target)

        for k, idx in enumerate(ring):
            mol.atoms[idx].position = aligned[k]

        _reposition_substituents(mol, ring, old_positions, corrected)
        corrected.update(ring)


def _find_rings(mol: Molecule) -> list[list[int]]:
    """Find all simple rings up to size 8 (heavy atoms only)."""
    heavy = {a.index for a in mol.atoms if a.symbol != "H"}
    adj: dict[int, list[int]] = {i: [] for i in heavy}
    for b in mol.bonds:
        if b.atom_i in heavy and b.atom_j in heavy:
            adj[b.atom_i].append(b.atom_j)
            adj[b.atom_j].append(b.atom_i)

    max_ring = 8
    rings: list[list[int]] = []
    ring_set: set[frozenset[int]] = set()

    for start in sorted(adj.keys()):
        stack: list[tuple[int, list[int]]] = [(start, [start])]
        while stack:
            node, path = stack.pop()
            if len(path) > max_ring:
                continue
            for neighbor in adj[node]:
                if neighbor == start and len(path) >= 3:
                    ring_key = frozenset(path)
                    if ring_key not in ring_set:
                        ring_set.add(ring_key)
                        rings.append(list(path))
                    continue
                if neighbor in path:
                    continue
                if neighbor < start:
                    continue
                if len(path) < max_ring:
                    stack.append((neighbor, path + [neighbor]))

    return rings


def _is_planar(mol: Molecule, ring: list[int]) -> bool:
    """True if the ring should be planar (all SP2 or size <= 3)."""
    if len(ring) <= 3:
        return True
    return all(
        mol.atoms[i].hybridization == Hybridization.SP2
        for i in ring
    )


def _avg_bond_len(mol: Molecule, ring: list[int]) -> float:
    """Average ideal bond length for the ring."""
    n = len(ring)
    all_carbon = all(mol.atoms[ring[k]].symbol == "C" for k in range(n))
    if all_carbon and n in (5, 6, 7):
        orders = []
        for k in range(n):
            b = mol.get_bond(ring[k], ring[(k + 1) % n])
            orders.append(b.order if b else 1)
        if 1 in orders and 2 in orders:
            return 1.40  # aromatic C-C

    total = 0.0
    for k in range(n):
        i = ring[k]
        j = ring[(k + 1) % n]
        b = mol.get_bond(i, j)
        order = b.order if b else 1
        total += ref_bond_length(mol.atoms[i].symbol, mol.atoms[j].symbol, order)
    return total / n


def _ideal_coords(mol: Molecule, ring: list[int]) -> np.ndarray:
    """Generate ideal ring coordinates."""
    if _is_planar(mol, ring):
        return _ideal_planar(mol, ring)
    return _ideal_puckered(mol, ring)


def _ideal_planar(mol: Molecule, ring: list[int]) -> np.ndarray:
    """Regular polygon in XY plane."""
    n = len(ring)
    avg_bl = _avg_bond_len(mol, ring)
    r = avg_bl / (2.0 * math.sin(math.pi / n))
    coords = np.zeros((n, 3))
    for k in range(n):
        angle = 2.0 * math.pi * k / n
        coords[k, 0] = r * math.cos(angle)
        coords[k, 1] = r * math.sin(angle)
    return coords


def _ideal_puckered(mol: Molecule, ring: list[int]) -> np.ndarray:
    """Puckered ring: chair for 6-membered, envelope otherwise."""
    n = len(ring)
    avg_bl = _avg_bond_len(mol, ring)
    coords = np.zeros((n, 3))

    if n == 6:
        pucker = 0.25
        r = math.sqrt(max(avg_bl**2 - (2.0 * pucker)**2, 0.01))
        for k in range(n):
            angle = 2.0 * math.pi * k / n
            coords[k, 0] = r * math.cos(angle)
            coords[k, 1] = r * math.sin(angle)
            coords[k, 2] = pucker if k % 2 == 0 else -pucker
    else:
        r = avg_bl / (2.0 * math.sin(math.pi / n))
        pucker = 0.3
        for k in range(n):
            angle = 2.0 * math.pi * k / n
            coords[k, 0] = r * math.cos(angle)
            coords[k, 1] = r * math.sin(angle)
            if k == n - 1:
                coords[k, 2] = pucker

    return coords


def kabsch_align(ideal: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Kabsch SVD alignment of ideal coordinates to target positions."""
    ideal_center = ideal.mean(axis=0)
    target_center = target.mean(axis=0)
    ideal_c = ideal - ideal_center
    target_c = target - target_center

    H = ideal_c.T @ target_c
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1.0, 1.0, np.sign(d)])
    R = Vt.T @ sign_matrix @ U.T

    return ideal_c @ R.T + target_center


def _align_fused(ideal: np.ndarray, ring: list[int],
                 shared: set[int], mol: Molecule) -> np.ndarray:
    """Align ideal ring using shared atoms as anchors for fused rings."""
    n = len(ring)
    ka, kb = None, None
    for k in range(n):
        k_next = (k + 1) % n
        if ring[k] in shared and ring[k_next] in shared:
            ka, kb = k, k_next
            break

    if ka is None:
        target = np.array([mol.atoms[i].position for i in ring])
        return kabsch_align(ideal, target)

    idx_a, idx_b = ring[ka], ring[kb]
    pa = mol.atoms[idx_a].position.copy()
    pb = mol.atoms[idx_b].position.copy()
    qa = ideal[ka].copy()
    qb = ideal[kb].copy()

    ideal_edge_mid = (qa + qb) / 2.0
    ideal_edge = qb - qa
    ideal_edge_hat = normalize(ideal_edge)
    ideal_center = ideal.mean(axis=0)
    ideal_to_center = ideal_center - ideal_edge_mid
    ideal_to_center -= np.dot(ideal_to_center, ideal_edge_hat) * ideal_edge_hat
    ideal_perp_hat = normalize(ideal_to_center)
    ideal_normal = normalize(np.cross(ideal_edge_hat, ideal_perp_hat))

    actual_edge_mid = (pa + pb) / 2.0
    actual_edge = pb - pa
    actual_edge_hat = normalize(actual_edge)

    non_shared_k = [k for k in range(n) if ring[k] not in shared]
    if non_shared_k:
        bfs_pos = mol.atoms[ring[non_shared_k[0]]].position
        hint = bfs_pos - actual_edge_mid
        hint -= np.dot(hint, actual_edge_hat) * actual_edge_hat
        if np.linalg.norm(hint) > 1e-10:
            actual_perp_hat = normalize(hint)
        else:
            perp = np.array([1.0, 0.0, 0.0])
            if abs(np.dot(actual_edge_hat, perp)) > 0.9:
                perp = np.array([0.0, 1.0, 0.0])
            actual_perp_hat = normalize(
                perp - np.dot(perp, actual_edge_hat) * actual_edge_hat)
    else:
        perp = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(actual_edge_hat, perp)) > 0.9:
            perp = np.array([0.0, 1.0, 0.0])
        actual_perp_hat = normalize(
            perp - np.dot(perp, actual_edge_hat) * actual_edge_hat)

    actual_normal = normalize(np.cross(actual_edge_hat, actual_perp_hat))

    ideal_frame = np.column_stack([ideal_edge_hat, ideal_perp_hat, ideal_normal])
    actual_frame = np.column_stack([actual_edge_hat, actual_perp_hat, actual_normal])
    R = actual_frame @ ideal_frame.T

    ideal_edge_len = np.linalg.norm(ideal_edge)
    actual_edge_len = np.linalg.norm(actual_edge)
    scale = actual_edge_len / ideal_edge_len if ideal_edge_len > 1e-10 else 1.0

    aligned = np.zeros_like(ideal)
    for k in range(n):
        rel = ideal[k] - ideal_edge_mid
        aligned[k] = actual_edge_mid + scale * (R @ rel)

    return aligned


def _reposition_substituents(mol: Molecule, ring: list[int],
                             old_positions: dict[int, np.ndarray],
                             protected: set[int]) -> None:
    """Reposition non-ring atoms after ring correction."""
    from collections import deque

    ring_set = set(ring) | protected

    for ring_idx in ring:
        old_pos = old_positions[ring_idx]
        new_pos = mol.atoms[ring_idx].position

        non_ring_nbs = [nb for nb in mol.neighbors(ring_idx) if nb not in ring_set]
        if not non_ring_nbs:
            continue

        ring_nbs = [nb for nb in mol.neighbors(ring_idx) if nb in ring_set]

        if len(ring_nbs) < 2:
            delta = new_pos - old_pos
            for nb in non_ring_nbs:
                visited: set[int] = set()
                queue = deque([nb])
                while queue:
                    idx = queue.popleft()
                    if idx in visited or idx in ring_set:
                        continue
                    visited.add(idx)
                    mol.atoms[idx].position = mol.atoms[idx].position + delta
                    for nb2 in mol.neighbors(idx):
                        if nb2 not in visited and nb2 not in ring_set:
                            queue.append(nb2)
            continue

        rn0, rn1 = ring_nbs[0], ring_nbs[1]
        old_v0 = old_positions.get(rn0, mol.atoms[rn0].position) - old_pos
        old_v1 = old_positions.get(rn1, mol.atoms[rn1].position) - old_pos
        old_normal = np.cross(old_v0, old_v1)
        old_norm_len = np.linalg.norm(old_normal)

        new_v0 = mol.atoms[rn0].position - new_pos
        new_v1 = mol.atoms[rn1].position - new_pos
        new_normal = np.cross(new_v0, new_v1)
        new_norm_len = np.linalg.norm(new_normal)

        if old_norm_len < 1e-10 or new_norm_len < 1e-10:
            delta = new_pos - old_pos
            for nb in non_ring_nbs:
                visited = set()
                queue = deque([nb])
                while queue:
                    idx = queue.popleft()
                    if idx in visited or idx in ring_set:
                        continue
                    visited.add(idx)
                    mol.atoms[idx].position = mol.atoms[idx].position + delta
                    for nb2 in mol.neighbors(idx):
                        if nb2 not in visited and nb2 not in ring_set:
                            queue.append(nb2)
            continue

        old_x = old_v0 / np.linalg.norm(old_v0)
        old_z = old_normal / old_norm_len
        old_y = np.cross(old_z, old_x)

        new_x = new_v0 / np.linalg.norm(new_v0)
        new_z = new_normal / new_norm_len
        new_y = np.cross(new_z, new_x)

        old_frame = np.column_stack([old_x, old_y, old_z])
        new_frame = np.column_stack([new_x, new_y, new_z])
        R = new_frame @ old_frame.T

        for nb in non_ring_nbs:
            visited = set()
            queue = deque([nb])
            while queue:
                idx = queue.popleft()
                if idx in visited or idx in ring_set:
                    continue
                visited.add(idx)
                rel = mol.atoms[idx].position - old_pos
                mol.atoms[idx].position = new_pos + R @ rel
                for nb2 in mol.neighbors(idx):
                    if nb2 not in visited and nb2 not in ring_set:
                        queue.append(nb2)
