"""Post-processing ring geometry correction for SMILES-parsed molecules.

The BFS-based coordinate assignment in parser.py places atoms along a tree
traversal.  Ring closure bonds simply connect two already-placed atoms without
adjusting positions, producing distorted rings.  This module detects rings,
generates ideal geometry (regular polygons for planar rings, chair for SP3
6-membered rings), and Kabsch-aligns the ideal coordinates onto the BFS
positions to preserve molecular context.

Public API
----------
correct_ring_geometry(mol) -> None
    Detect rings, replace their coordinates with correct geometry, and
    reposition substituents.  Modifies the Molecule in-place.
"""

from __future__ import annotations

import math
from collections import deque

import numpy as np

from molbuilder.core.bond_data import bond_length
from molbuilder.core.geometry import normalize
from molbuilder.molecule.graph import Molecule, Hybridization


# ===================================================================
# Ring detection (DFS-based, all atoms, max size 8)
# ===================================================================

def find_all_rings(mol: Molecule) -> list[list[int]]:
    """Find all simple rings up to size 8 in the molecular graph.

    Uses bounded DFS from each atom, deduplicating by frozenset.
    Only considers heavy atoms (non-hydrogen) for ring membership.

    Returns
    -------
    list[list[int]]
        Each element is an ordered list of atom indices forming a ring.
    """
    # Build adjacency for heavy atoms only
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


# ===================================================================
# Ring classification
# ===================================================================

def _is_planar_ring(mol: Molecule, ring: list[int]) -> bool:
    """True if the ring should be planar.

    A ring is planar if all ring atoms are SP2-hybridized, or if the
    ring is 3-membered (always planar by geometry).
    """
    if len(ring) <= 3:
        return True
    return all(
        mol.atoms[i].hybridization == Hybridization.SP2
        for i in ring
    )


# ===================================================================
# Ideal geometry generation
# ===================================================================

def _avg_ring_bond_length(mol: Molecule, ring: list[int]) -> float:
    """Average ideal bond length for consecutive atoms in the ring.

    For aromatic rings (alternating single/double in Kekule form), uses
    the known aromatic bond length (1.40 A for C-C) instead of averaging
    single and double.
    """
    n = len(ring)
    # Check if this is an aromatic ring (all sp2 carbons with alternating
    # single/double bonds -- Kekule representation of aromatic system)
    all_carbon = all(mol.atoms[ring[k]].symbol == "C" for k in range(n))
    if all_carbon and n in (5, 6, 7):
        orders = []
        for k in range(n):
            b = mol.get_bond(ring[k], ring[(k + 1) % n])
            orders.append(b.order if b else 1)
        has_single = 1 in orders
        has_double = 2 in orders
        if has_single and has_double:
            # Aromatic: use experimental aromatic C-C bond length
            return 1.40

    total = 0.0
    for k in range(n):
        i = ring[k]
        j = ring[(k + 1) % n]
        b = mol.get_bond(i, j)
        order = b.order if b else 1
        total += bond_length(mol.atoms[i].symbol, mol.atoms[j].symbol, order)
    return total / n


def _ideal_planar_ring(mol: Molecule, ring: list[int]) -> np.ndarray:
    """Generate ideal planar ring coordinates as a regular polygon in XY.

    The polygon radius is chosen so that edge lengths match the average
    expected bond length for the ring.

    Returns
    -------
    ndarray of shape (N, 3)
        Ideal positions with z=0.
    """
    n = len(ring)
    avg_bl = _avg_ring_bond_length(mol, ring)
    # edge = 2 * r * sin(pi/n)  =>  r = edge / (2 * sin(pi/n))
    r = avg_bl / (2.0 * math.sin(math.pi / n))

    coords = np.zeros((n, 3))
    for k in range(n):
        angle = 2.0 * math.pi * k / n
        coords[k, 0] = r * math.cos(angle)
        coords[k, 1] = r * math.sin(angle)
    return coords


def _ideal_puckered_ring(mol: Molecule, ring: list[int]) -> np.ndarray:
    """Generate ideal puckered ring coordinates.

    For 6-membered SP3 rings: chair conformation (alternating +/- 0.25 A).
    For other sizes: envelope approximation (one atom displaced).

    Returns
    -------
    ndarray of shape (N, 3)
        Ideal positions.
    """
    n = len(ring)
    avg_bl = _avg_ring_bond_length(mol, ring)

    coords = np.zeros((n, 3))

    if n == 6:
        # Chair conformation: alternating z displacement.
        # Adjacent atoms differ by 2*pucker in z, so in-plane radius must
        # be reduced so that 3D C-C distance equals avg_bl:
        #   avg_bl^2 = edge_planar^2 + (2*pucker)^2
        #   edge_planar = 2 * r * sin(pi/6) = r
        #   => r = sqrt(avg_bl^2 - (2*pucker)^2)
        pucker = 0.25  # Angstroms
        r = math.sqrt(max(avg_bl ** 2 - (2.0 * pucker) ** 2, 0.01))
        for k in range(n):
            angle = 2.0 * math.pi * k / n
            coords[k, 0] = r * math.cos(angle)
            coords[k, 1] = r * math.sin(angle)
            coords[k, 2] = pucker if k % 2 == 0 else -pucker
    else:
        # Envelope: last atom displaced, rest roughly planar
        r = avg_bl / (2.0 * math.sin(math.pi / n))
        pucker = 0.3
        for k in range(n):
            angle = 2.0 * math.pi * k / n
            coords[k, 0] = r * math.cos(angle)
            coords[k, 1] = r * math.sin(angle)
            if k == n - 1:
                coords[k, 2] = pucker

    return coords


# ===================================================================
# Kabsch alignment (SVD-based)
# ===================================================================

def _kabsch_align(ideal: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Optimally align ideal coordinates to target positions via SVD.

    Uses the Kabsch algorithm:
    1. Center both sets
    2. Compute cross-covariance H = ideal_centered.T @ target_centered
    3. SVD: U, S, Vt = svd(H)
    4. Optimal rotation R = Vt.T @ U.T (with reflection fix)
    5. Apply rotation and translate to target centroid

    Parameters
    ----------
    ideal : ndarray (N, 3)
        Ideal ring coordinates to be aligned.
    target : ndarray (N, 3)
        BFS-placed coordinates to align to.

    Returns
    -------
    ndarray (N, 3)
        Aligned ideal coordinates.
    """
    ideal_center = ideal.mean(axis=0)
    target_center = target.mean(axis=0)

    ideal_c = ideal - ideal_center
    target_c = target - target_center

    H = ideal_c.T @ target_c
    U, S, Vt = np.linalg.svd(H)

    # Reflection fix: ensure proper rotation (det = +1)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1.0, 1.0, np.sign(d)])
    R = Vt.T @ sign_matrix @ U.T

    aligned = ideal_c @ R.T + target_center
    return aligned


# ===================================================================
# Substituent repositioning
# ===================================================================

def _reposition_substituents(
    mol: Molecule,
    ring_indices: list[int],
    old_positions: dict[int, np.ndarray],
    protected: set[int] | None = None,
) -> None:
    """Reposition non-ring substituents after ring geometry correction.

    For each ring atom, compute the displacement from old to new position
    and the rotation from old local frame to new local frame.  BFS through
    non-ring subtrees from each ring atom and apply the rigid transform.

    Parameters
    ----------
    protected : set[int] | None
        Additional atom indices to exclude from repositioning (e.g.
        atoms already corrected by other rings in fused systems).
    """
    ring_set = set(ring_indices)
    if protected:
        ring_set = ring_set | protected

    for ring_idx in ring_indices:
        old_pos = old_positions[ring_idx]
        new_pos = mol.atoms[ring_idx].position

        # Find non-ring neighbors
        non_ring_neighbors = [
            nb for nb in mol.neighbors(ring_idx) if nb not in ring_set
        ]
        if not non_ring_neighbors:
            continue

        # Compute old and new local frame from ring neighbors
        ring_neighbors = [
            nb for nb in mol.neighbors(ring_idx) if nb in ring_set
        ]

        if len(ring_neighbors) < 2:
            # Simple translation for atoms with < 2 ring neighbors
            delta = new_pos - old_pos
            for nb in non_ring_neighbors:
                _translate_subtree(mol, nb, ring_set, delta)
            continue

        # Build local frames from the two ring neighbors
        rn0, rn1 = ring_neighbors[0], ring_neighbors[1]

        old_v0 = old_positions[rn0] - old_pos
        old_v1 = old_positions[rn1] - old_pos
        old_normal = np.cross(old_v0, old_v1)
        old_norm_len = np.linalg.norm(old_normal)

        new_v0 = mol.atoms[rn0].position - new_pos
        new_v1 = mol.atoms[rn1].position - new_pos
        new_normal = np.cross(new_v0, new_v1)
        new_norm_len = np.linalg.norm(new_normal)

        if old_norm_len < 1e-10 or new_norm_len < 1e-10:
            # Degenerate: just translate
            delta = new_pos - old_pos
            for nb in non_ring_neighbors:
                _translate_subtree(mol, nb, ring_set, delta)
            continue

        # Build orthonormal bases
        old_x = old_v0 / np.linalg.norm(old_v0)
        old_z = old_normal / old_norm_len
        old_y = np.cross(old_z, old_x)

        new_x = new_v0 / np.linalg.norm(new_v0)
        new_z = new_normal / new_norm_len
        new_y = np.cross(new_z, new_x)

        # Rotation matrix from old frame to new frame
        old_frame = np.column_stack([old_x, old_y, old_z])
        new_frame = np.column_stack([new_x, new_y, new_z])
        R = new_frame @ old_frame.T

        # Apply rigid transform to substituent subtrees
        for nb in non_ring_neighbors:
            _transform_subtree(mol, nb, ring_set, old_pos, new_pos, R)


def _translate_subtree(
    mol: Molecule,
    start: int,
    exclude: set[int],
    delta: np.ndarray,
) -> None:
    """BFS from start, translating all atoms not in exclude by delta."""
    visited: set[int] = set()
    queue = deque([start])

    while queue:
        idx = queue.popleft()
        if idx in visited or idx in exclude:
            continue
        visited.add(idx)
        mol.atoms[idx].position = mol.atoms[idx].position + delta
        for nb in mol.neighbors(idx):
            if nb not in visited and nb not in exclude:
                queue.append(nb)


def _transform_subtree(
    mol: Molecule,
    start: int,
    exclude: set[int],
    old_pivot: np.ndarray,
    new_pivot: np.ndarray,
    R: np.ndarray,
) -> None:
    """BFS from start, applying rotation R (around old_pivot) + translation."""
    visited: set[int] = set()
    queue = deque([start])

    while queue:
        idx = queue.popleft()
        if idx in visited or idx in exclude:
            continue
        visited.add(idx)
        rel = mol.atoms[idx].position - old_pivot
        mol.atoms[idx].position = new_pivot + R @ rel
        for nb in mol.neighbors(idx):
            if nb not in visited and nb not in exclude:
                queue.append(nb)


# ===================================================================
# Fused ring alignment (edge-anchored)
# ===================================================================

def _align_fused_ring(
    ideal: np.ndarray,
    ring: list[int],
    shared_atom_set: set[int],
    mol: Molecule,
) -> np.ndarray:
    """Align ideal ring coordinates using shared atoms as exact anchors.

    For fused rings, the shared edge is already at the correct position
    from the first ring's correction.  This function places the ideal
    polygon so the shared edge matches exactly, and the remaining atoms
    are on the correct side.

    Parameters
    ----------
    ideal : ndarray (N, 3)
        Ideal ring coordinates (regular polygon).
    ring : list[int]
        Atom indices in traversal order.
    shared_atom_set : set[int]
        Atom indices already corrected by previous rings.
    mol : Molecule
        The molecule (shared atoms already at corrected positions).

    Returns
    -------
    ndarray (N, 3)
        Aligned coordinates.
    """
    n = len(ring)

    # Find two shared atoms that are ADJACENT in the ring (the shared edge).
    # In fused rings, the shared bond connects atoms consecutive in the ring.
    ka, kb = None, None
    for k in range(n):
        k_next = (k + 1) % n
        if ring[k] in shared_atom_set and ring[k_next] in shared_atom_set:
            ka, kb = k, k_next
            break

    if ka is None:
        # No adjacent shared pair found, fall back to Kabsch
        target = np.array([mol.atoms[i].position for i in ring])
        return _kabsch_align(ideal, target)

    idx_a, idx_b = ring[ka], ring[kb]

    pa = mol.atoms[idx_a].position.copy()
    pb = mol.atoms[idx_b].position.copy()

    qa = ideal[ka].copy()
    qb = ideal[kb].copy()

    # -- Build ideal frame from the shared edge --
    ideal_edge_mid = (qa + qb) / 2.0
    ideal_edge = qb - qa
    ideal_edge_hat = normalize(ideal_edge)

    # Direction from shared edge midpoint to polygon center (in-plane perp)
    ideal_center = ideal.mean(axis=0)
    ideal_to_center = ideal_center - ideal_edge_mid
    ideal_to_center -= np.dot(ideal_to_center, ideal_edge_hat) * ideal_edge_hat
    ideal_perp_hat = normalize(ideal_to_center)

    # Normal to the ideal polygon plane
    ideal_normal = normalize(np.cross(ideal_edge_hat, ideal_perp_hat))

    # -- Build actual frame from the shared edge --
    actual_edge_mid = (pa + pb) / 2.0
    actual_edge = pb - pa
    actual_edge_hat = normalize(actual_edge)

    # Determine which side of the shared edge the new ring atoms should go.
    # Use BFS positions of non-shared atoms as a directional hint.
    non_shared_k = [k for k in range(n) if ring[k] not in shared_atom_set]

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

    # -- Compute rotation from ideal frame to actual frame --
    ideal_frame = np.column_stack([ideal_edge_hat, ideal_perp_hat, ideal_normal])
    actual_frame = np.column_stack([actual_edge_hat, actual_perp_hat, actual_normal])
    R = actual_frame @ ideal_frame.T

    # Scale to match shared edge length
    ideal_edge_len = np.linalg.norm(ideal_edge)
    actual_edge_len = np.linalg.norm(actual_edge)
    scale = actual_edge_len / ideal_edge_len if ideal_edge_len > 1e-10 else 1.0

    # Apply transform: rotate + scale around shared edge midpoint
    aligned = np.zeros_like(ideal)
    for k in range(n):
        rel = ideal[k] - ideal_edge_mid
        aligned[k] = actual_edge_mid + scale * (R @ rel)

    return aligned


# ===================================================================
# Public entry point
# ===================================================================

def correct_ring_geometry(mol: Molecule) -> None:
    """Detect rings and correct their 3D geometry in-place.

    Algorithm
    ---------
    1. Find all rings (up to size 8) via DFS.
    2. Sort smallest-first so small rings are corrected before larger ones.
    3. For each ring:
       a. Skip if all atoms already corrected (redundant macrocycle)
       b. Save old positions
       c. Generate ideal geometry (planar or puckered)
       d. For fused rings (>=2 shared atoms): edge-anchored alignment
          For standalone rings: Kabsch alignment to BFS positions
       e. Update atom positions
       f. Reposition substituent subtrees
    """
    rings = find_all_rings(mol)
    if not rings:
        return

    # Sort smallest-first
    rings.sort(key=len)

    # Track which atoms have been corrected (for fused ring handling)
    corrected: set[int] = set()

    for ring in rings:
        # Skip redundant rings where all atoms are already corrected
        # (e.g. the 10-membered perimeter of naphthalene)
        if all(idx in corrected for idx in ring):
            continue

        # Save old positions
        old_positions = {i: mol.atoms[i].position.copy() for i in ring}

        # Generate ideal geometry
        if _is_planar_ring(mol, ring):
            ideal = _ideal_planar_ring(mol, ring)
        else:
            ideal = _ideal_puckered_ring(mol, ring)

        # Check for shared atoms with already-corrected rings
        shared = {idx for idx in ring if idx in corrected}

        if len(shared) >= 2:
            # Fused ring: use edge-anchored alignment
            aligned = _align_fused_ring(ideal, ring, shared, mol)
        else:
            # Standard Kabsch alignment to BFS positions
            target = np.array([mol.atoms[i].position for i in ring])
            aligned = _kabsch_align(ideal, target)

        # Update atom positions
        for k, idx in enumerate(ring):
            mol.atoms[idx].position = aligned[k]

        # Reposition substituents (protect already-corrected ring atoms)
        _reposition_substituents(mol, ring, old_positions, corrected)

        corrected.update(ring)
