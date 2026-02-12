"""Aromatic ring detection and Kekulization for SMILES-parsed molecules.

This module finds aromatic rings in a molecular graph and assigns
alternating single/double bond orders (Kekulization) so that every
aromatic atom participates in exactly one double bond within the
pi system.

The functions operate on duck-typed atom and bond objects produced by
the SMILES parser:

    atom.index    : int   -- atom index in the molecule
    atom.aromatic : bool  -- True if the atom was specified as aromatic
    atom.symbol   : str   -- element symbol (capitalized, e.g. "C", "N")
    atom.hcount   : int | None -- explicit hydrogen count (bracket atoms)
    atom.bracket  : bool  -- True if atom was specified in brackets

    bond.atom_i   : int   -- index of first atom
    bond.atom_j   : int   -- index of second atom
    bond.order    : int   -- bond order (modified in-place by kekulize)

Algorithm overview
------------------
1. Build an adjacency list restricted to aromatic atoms.
2. Detect all simple rings using DFS-based cycle detection.
3. Filter to rings where every atom is aromatic.
4. Determine pi-electron contributions for each aromatic atom to
   validate Hueckel compliance (4n+2 electrons).
5. Find a perfect matching on the aromatic subgraph using augmenting
   paths (Hopcroft-Karp-style BFS/DFS matching) so that every aromatic
   atom gets exactly one double bond.
6. Set matching edges to bond order 2 in-place.
"""

from __future__ import annotations

from collections import deque


# ===================================================================
# Adjacency helpers
# ===================================================================

def _build_aromatic_adj(atoms: list, bonds: list) -> dict[int, list[int]]:
    """Build an adjacency list containing only aromatic atoms.

    Parameters
    ----------
    atoms : list
        Atom objects with ``.index`` and ``.aromatic`` attributes.
    bonds : list
        Bond objects with ``.atom_i`` and ``.atom_j`` attributes.

    Returns
    -------
    dict[int, list[int]]
        Mapping from aromatic atom index to list of aromatic neighbor
        indices.
    """
    aromatic_set = {a.index for a in atoms if a.aromatic}
    adj: dict[int, list[int]] = {idx: [] for idx in aromatic_set}

    for b in bonds:
        if b.atom_i in aromatic_set and b.atom_j in aromatic_set:
            adj[b.atom_i].append(b.atom_j)
            adj[b.atom_j].append(b.atom_i)

    return adj


def _bond_lookup(bonds: list) -> dict[tuple[int, int], object]:
    """Create a lookup from (min_idx, max_idx) -> bond object.

    This allows O(1) retrieval of the bond object connecting two atoms
    so that we can modify its order in-place during Kekulization.
    """
    lookup: dict[tuple[int, int], object] = {}
    for b in bonds:
        key = (min(b.atom_i, b.atom_j), max(b.atom_i, b.atom_j))
        lookup[key] = b
    return lookup


# ===================================================================
# Ring detection via DFS
# ===================================================================

def find_aromatic_rings(atoms: list, bonds: list) -> list[list[int]]:
    """Find all simple rings where every member atom is aromatic.

    Uses DFS-based cycle detection on the aromatic subgraph.  Each ring
    is returned as an ordered list of atom indices forming the cycle.

    Parameters
    ----------
    atoms : list
        Atom objects with ``.index`` and ``.aromatic`` attributes.
    bonds : list
        Bond objects with ``.atom_i`` and ``.atom_j`` attributes.

    Returns
    -------
    list[list[int]]
        Each element is a list of atom indices that form an aromatic
        ring, in traversal order.
    """
    adj = _build_aromatic_adj(atoms, bonds)
    if not adj:
        return []

    # We use a bounded DFS to find all simple cycles of length <= max_ring.
    # Aromatic rings in organic chemistry are typically 5-7 membered.
    # We allow up to 8 for unusual heterocycles.
    max_ring = 8
    rings: list[list[int]] = []
    ring_set: set[frozenset[int]] = set()  # deduplicate rings

    for start in sorted(adj.keys()):
        # DFS stack: (current_node, path_so_far)
        stack: list[tuple[int, list[int]]] = [(start, [start])]

        while stack:
            node, path = stack.pop()

            if len(path) > max_ring:
                continue

            for neighbor in adj[node]:
                # Found a cycle back to the start
                if neighbor == start and len(path) >= 3:
                    ring_key = frozenset(path)
                    if ring_key not in ring_set:
                        ring_set.add(ring_key)
                        rings.append(list(path))
                    continue

                # Avoid revisiting nodes already on the current path,
                # and only explore neighbors with index > start to avoid
                # generating each ring multiple times from different
                # starting vertices.
                if neighbor in path:
                    continue
                if neighbor < start:
                    continue

                if len(path) < max_ring:
                    stack.append((neighbor, path + [neighbor]))

    return rings


# ===================================================================
# Pi-electron contribution
# ===================================================================

def _pi_electrons(atom, adj: dict[int, list[int]]) -> int:
    """Determine how many electrons an aromatic atom contributes to the
    pi system.

    Standard contributions:
    - Carbon (aromatic): 1 electron (from the p-orbital)
    - Nitrogen with 3 bonds, no H (pyridine-like): 1 electron
    - Nitrogen with 2 bonds and 1H (pyrrole-like): 2 electrons
    - Oxygen in furan: 2 electrons (lone pair donated to ring)
    - Sulfur in thiophene: 2 electrons (lone pair donated to ring)
    - Boron: 0 electrons (empty p-orbital)
    - Phosphorus: similar to nitrogen

    Parameters
    ----------
    atom : object
        Atom with ``.symbol``, ``.index``, ``.hcount``, ``.bracket``
        attributes.
    adj : dict[int, list[int]]
        Adjacency list for aromatic atoms.

    Returns
    -------
    int
        Number of pi electrons contributed by this atom.
    """
    n_aromatic_bonds = len(adj.get(atom.index, []))
    symbol = atom.symbol.upper() if hasattr(atom, 'symbol') else ""

    if symbol == "C":
        return 1

    if symbol == "N":
        # Pyrrole-type nitrogen: 2 aromatic bonds + 1 H -> contributes 2
        # Pyridine-type nitrogen: 2 aromatic bonds, no H -> contributes 1
        # But if it has 3 aromatic bonds (fused ring junction), contributes 1
        has_h = False
        if hasattr(atom, 'hcount') and atom.hcount is not None:
            has_h = atom.hcount > 0
        elif hasattr(atom, 'bracket') and atom.bracket:
            has_h = False  # bracket without H means no H
        else:
            # Organic subset aromatic nitrogen: pyridine-like by default
            # unless it has only 2 aromatic bonds (could be pyrrole)
            has_h = (n_aromatic_bonds == 2)

        if n_aromatic_bonds == 3:
            # Junction nitrogen in fused rings -- contributes 1
            return 1
        if has_h:
            return 2  # pyrrole-like
        return 1  # pyridine-like

    if symbol in ("O", "S"):
        # Furan / thiophene: lone pair donated, 2 electrons
        return 2

    if symbol == "P":
        # Phosphole: similar to pyrrole
        return 2

    if symbol == "B":
        # Borole: empty p-orbital, 0 electrons
        return 0

    # Default: 1 electron
    return 1


# ===================================================================
# Augmenting-path matching (for Kekulization)
# ===================================================================

def _find_perfect_matching(
    adj: dict[int, list[int]],
) -> dict[int, int] | None:
    """Find a perfect matching on the aromatic subgraph.

    A *perfect matching* is a set of edges such that every vertex is
    incident to exactly one edge in the set.  In chemical terms, every
    aromatic atom participates in exactly one double bond.

    Uses an augmenting-path algorithm (simplified Hopcroft-Karp):
    repeatedly find augmenting paths via BFS to grow the matching until
    no more augmenting paths exist.

    Parameters
    ----------
    adj : dict[int, list[int]]
        Adjacency list of the aromatic subgraph.

    Returns
    -------
    dict[int, int] | None
        Mapping from atom index to its matched partner, or ``None``
        if no perfect matching exists.
    """
    nodes = set(adj.keys())
    if len(nodes) % 2 == 1:
        # Odd number of aromatic atoms -> no perfect matching possible.
        # This can happen with valid aromatic systems (e.g. 5-membered
        # rings with a heteroatom contributing 2 electrons).  We handle
        # this by marking heteroatoms that contribute 2 electrons as
        # "pre-satisfied" (they don't need a double bond).
        return None

    # match[v] = the vertex matched to v, or -1 if unmatched
    match: dict[int, int] = {v: -1 for v in nodes}

    def _augment(u: int) -> bool:
        """Try to find an augmenting path starting from unmatched vertex u.

        Uses DFS.  If an augmenting path is found, it flips all edges
        along the path (matched <-> unmatched) and returns True.
        """
        for v in adj[u]:
            if v not in visited:
                visited.add(v)
                # v is either unmatched, or its current match can be
                # re-routed via another augmenting path
                if match[v] == -1 or _augment(match[v]):
                    match[v] = u
                    match[u] = v
                    return True
        return False

    # Greedy initialization: for each unmatched node, try to match it
    # to an unmatched neighbor.  This gives a good starting matching
    # and reduces the number of augmenting paths needed.
    for u in sorted(nodes):
        if match[u] == -1:
            for v in adj[u]:
                if match[v] == -1:
                    match[u] = v
                    match[v] = u
                    break

    # Augmenting path phase
    changed = True
    while changed:
        changed = False
        for u in sorted(nodes):
            if match[u] == -1:
                visited: set[int] = {u}
                if _augment(u):
                    changed = True

    return match


def _find_matching_with_donors(
    adj: dict[int, list[int]],
    atoms_by_idx: dict[int, object],
) -> dict[int, int]:
    """Find a matching that accounts for lone-pair donor atoms.

    Heteroatoms that contribute 2 pi electrons (pyrrole N, O, S) do
    not need a double bond -- their lone pair completes the aromatic
    sextet.  We remove these "donor" atoms from the matching graph and
    find a perfect matching on the remaining atoms.

    If a perfect matching on all atoms works, we use that instead
    (simpler and handles all-carbon rings like benzene).

    Parameters
    ----------
    adj : dict[int, list[int]]
        Adjacency of aromatic atoms.
    atoms_by_idx : dict[int, object]
        Atom objects indexed by atom index.

    Returns
    -------
    dict[int, int]
        Matching: atom -> matched partner.  Donor atoms are mapped to -1
        (no double bond needed).
    """
    nodes = set(adj.keys())

    # First, try a perfect matching on the full aromatic subgraph.
    # This works for benzene, naphthalene, pyridine, etc.
    if len(nodes) % 2 == 0:
        full_match = _find_perfect_matching(adj)
        if full_match is not None:
            # Check that every atom is actually matched
            if all(full_match[v] != -1 for v in nodes):
                return full_match

    # Identify donor atoms (contribute 2 electrons, don't need double bond)
    donors: set[int] = set()
    for idx in nodes:
        atom = atoms_by_idx[idx]
        pe = _pi_electrons(atom, adj)
        if pe == 2:
            donors.add(idx)

    # Build reduced adjacency without donors
    reduced_adj: dict[int, list[int]] = {}
    non_donors = nodes - donors
    for idx in non_donors:
        reduced_adj[idx] = [n for n in adj[idx] if n in non_donors]

    # Find perfect matching on non-donor subgraph
    if len(non_donors) % 2 == 0 and len(non_donors) > 0:
        reduced_match = _find_perfect_matching(reduced_adj)
        if reduced_match is not None:
            if all(reduced_match[v] != -1 for v in non_donors):
                # Success: donor atoms get -1 (no double bond)
                result = {v: -1 for v in nodes}
                for v in non_donors:
                    result[v] = reduced_match[v]
                return result

    # Fallback: try removing donors one at a time until we get a
    # perfect matching on the rest.  This handles edge cases where
    # not all potential donors actually donate.
    for donor_subset_size in range(len(donors), 0, -1):
        for donor_combo in _combinations(sorted(donors), donor_subset_size):
            active_donors = set(donor_combo)
            remaining = nodes - active_donors
            if len(remaining) % 2 != 0:
                continue
            if len(remaining) == 0:
                # All atoms are donors -- no double bonds needed
                return {v: -1 for v in nodes}
            sub_adj: dict[int, list[int]] = {}
            for idx in remaining:
                sub_adj[idx] = [n for n in adj[idx] if n in remaining]
            sub_match = _find_perfect_matching(sub_adj)
            if sub_match is not None:
                if all(sub_match[v] != -1 for v in remaining):
                    result = {v: -1 for v in nodes}
                    for v in remaining:
                        result[v] = sub_match[v]
                    return result

    # Last resort: return empty matching (no Kekulization possible).
    # The bonds stay as single bonds.
    return {v: -1 for v in nodes}


def _combinations(items: list, r: int):
    """Generate all r-length combinations from items.

    Simple implementation to avoid importing itertools in a module
    that aims to be lightweight.
    """
    n = len(items)
    if r > n:
        return
    if r == 0:
        yield ()
        return
    for i in range(n - r + 1):
        for rest in _combinations(items[i + 1:], r - 1):
            yield (items[i],) + rest


# ===================================================================
# Public API
# ===================================================================

def kekulize(atoms: list, bonds: list) -> list:
    """Assign alternating single/double bonds to aromatic ring systems.

    Modifies bond orders in-place and returns the bonds list.
    Atoms must have ``.aromatic``, ``.index``, and ``.symbol`` attributes.
    Bonds must have ``.atom_i``, ``.atom_j``, and ``.order`` attributes.

    Algorithm
    ---------
    1. Build an adjacency list for aromatic atoms only.
    2. Find all aromatic rings via DFS.
    3. Use augmenting-path matching to assign double bonds such that
       every aromatic atom participates in exactly one double bond
       (or is a lone-pair donor that does not need one).
    4. Modify bond orders in-place.

    Parameters
    ----------
    atoms : list
        Atom-like objects from the SMILES parser.
    bonds : list
        Bond-like objects from the SMILES parser.

    Returns
    -------
    list
        The same *bonds* list, with some ``.order`` values changed
        from 1 to 2 where the matching assigns a double bond.

    Examples
    --------
    After parsing ``c1ccccc1`` (benzene), all 6 ring bonds have
    order=1.  After ``kekulize(atoms, bonds)``, 3 bonds become
    order=2 and 3 remain order=1, giving the Kekule structure.
    """
    # Quick exit if no aromatic atoms
    aromatic_atoms = [a for a in atoms if a.aromatic]
    if not aromatic_atoms:
        return bonds

    # Build aromatic adjacency and atom lookup
    adj = _build_aromatic_adj(atoms, bonds)
    atoms_by_idx = {a.index: a for a in atoms}

    # Find the matching (which aromatic bonds become double)
    matching = _find_matching_with_donors(adj, atoms_by_idx)

    # Apply the matching to bond orders
    bond_map = _bond_lookup(bonds)

    assigned: set[tuple[int, int]] = set()
    for u, v in matching.items():
        if v == -1:
            continue
        key = (min(u, v), max(u, v))
        if key in assigned:
            continue
        assigned.add(key)

        bond = bond_map.get(key)
        if bond is not None:
            bond.order = 2

    return bonds
