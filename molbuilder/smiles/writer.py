"""SMILES writer: Molecule -> canonical SMILES string.

Uses Morgan's algorithm for canonical atom numbering, then DFS traversal
to generate the SMILES string.
"""

from __future__ import annotations

from collections import deque

from molbuilder.molecule.graph import Molecule
from molbuilder.core.elements import SYMBOL_TO_Z
from molbuilder.smiles.tokenizer import ORGANIC_SUBSET, DEFAULT_VALENCE


# ===================================================================
# Morgan canonical ordering
# ===================================================================

def _morgan_canonical_order(mol: Molecule) -> dict[int, int]:
    """Compute canonical atom ranks using Morgan's extended connectivity.

    Algorithm
    ---------
    1. Initialize each atom's invariant to its degree (number of
       neighbours).
    2. Iteratively replace each atom's invariant with the sum of its
       neighbours' invariants until the number of distinct values
       stabilises.
    3. Break ties using atomic number (higher Z = higher rank).
    4. Return a mapping from atom index to rank (0 = lowest priority,
       first in canonical order).

    Parameters
    ----------
    mol : Molecule
        The molecule whose atoms should be ranked.

    Returns
    -------
    dict[int, int]
        Mapping from atom index to canonical rank.
    """
    n = len(mol.atoms)
    if n == 0:
        return {}

    # Initial connectivity value = degree
    ec = [len(mol.neighbors(i)) for i in range(n)]

    prev_classes = 0
    for _iteration in range(100):
        new_ec = [0] * n
        for i in range(n):
            new_ec[i] = sum(ec[j] for j in mol.neighbors(i))

        # Count distinct classes
        classes = len(set(new_ec))
        if classes == prev_classes:
            break
        prev_classes = classes
        ec = new_ec

    # Build (ec_value, atomic_number, original_index) for stable sort
    sort_keys = []
    for i in range(n):
        z = SYMBOL_TO_Z.get(mol.atoms[i].symbol, 0)
        sort_keys.append((ec[i], z, i))

    sorted_atoms = sorted(sort_keys)
    order: dict[int, int] = {}
    for rank, (_, _, idx) in enumerate(sorted_atoms):
        order[idx] = rank

    return order


# ===================================================================
# Connected components
# ===================================================================

def _connected_components(mol: Molecule) -> list[list[int]]:
    """Return lists of atom indices for each connected fragment.

    Considers only heavy (non-hydrogen) atoms.
    """
    n = len(mol.atoms)
    heavy = {i for i in range(n) if mol.atoms[i].symbol != "H"}
    visited: set[int] = set()
    components: list[list[int]] = []

    for start in heavy:
        if start in visited:
            continue
        comp: list[int] = []
        stack = [start]
        while stack:
            cur = stack.pop()
            if cur in visited:
                continue
            visited.add(cur)
            if cur not in heavy:
                continue
            comp.append(cur)
            for nb in mol.neighbors(cur):
                if nb not in visited and nb in heavy:
                    stack.append(nb)
        if comp:
            components.append(comp)

    return components


# ===================================================================
# Bond order symbol
# ===================================================================

def _bond_symbol(mol: Molecule, i: int, j: int) -> str:
    """Return the SMILES bond symbol between atoms *i* and *j*.

    Single bonds return ``""`` (implicit), double ``"="``, triple ``"#"``.
    """
    bond = mol.get_bond(i, j)
    if bond is None:
        return ""
    if bond.order == 2:
        return "="
    if bond.order == 3:
        return "#"
    return ""


# ===================================================================
# DFS SMILES generation
# ===================================================================

def _dfs_smiles(mol: Molecule, order: dict[int, int],
                component: list[int]) -> str:
    """Generate a SMILES string for one connected component via DFS.

    Parameters
    ----------
    mol : Molecule
        The full molecule.
    order : dict[int, int]
        Canonical ranking from ``_morgan_canonical_order``.
    component : list[int]
        Atom indices belonging to this fragment (heavy atoms only).

    Returns
    -------
    str
        SMILES string for the fragment.
    """
    if not component:
        return ""

    comp_set = set(component)

    # Start DFS from the atom with the lowest canonical rank in this component
    start = min(component, key=lambda i: order.get(i, 0))

    visited: set[int] = set()
    # Track ring closure pairs: when DFS finds a back-edge, both
    # endpoints get a digit.
    ring_bonds: dict[tuple[int, int], int] = {}
    next_ring_digit = 1

    parts: list[str] = []

    def _heavy_neighbors(idx: int) -> list[int]:
        """Return non-hydrogen neighbours in this component, sorted by
        canonical rank (highest rank = most preferred child = last in
        the sort so it gets visited first on the main chain)."""
        nbs = []
        for nb in mol.neighbors(idx):
            if mol.atoms[nb].symbol == "H":
                continue
            if nb not in comp_set:
                continue
            nbs.append(nb)
        # Sort: the *last* element becomes the "main chain" child (no
        # parentheses).  Lower rank = visited first = branch.
        nbs.sort(key=lambda x: order.get(x, 0))
        return nbs

    def _atom_str(idx: int) -> str:
        """Return the SMILES atom token for atom *idx*.

        Outputs bracket notation ``[<isotope><symbol><chirality><hcount><charge>]``
        when the atom has non-default properties (chirality, isotope, charge,
        or is not in the organic subset).  Organic subset atoms without
        special properties are written without brackets.
        """
        atom = mol.atoms[idx]
        sym = atom.symbol
        has_chirality = atom.chirality is not None
        has_isotope = atom.isotope is not None
        has_charge = atom.formal_charge != 0
        needs_bracket = has_chirality or has_isotope or has_charge or sym not in ORGANIC_SUBSET

        if not needs_bracket:
            return sym

        # Build bracket atom string: [<isotope><symbol><chirality><Hn><charge>]
        parts: list[str] = []
        if has_isotope:
            parts.append(str(atom.isotope))
        parts.append(sym)
        if has_chirality:
            parts.append(atom.chirality)

        # Compute implicit H count: count explicit H neighbours
        h_count = sum(1 for nb in mol.neighbors(idx)
                      if mol.atoms[nb].symbol == "H")
        if h_count == 1:
            parts.append("H")
        elif h_count > 1:
            parts.append(f"H{h_count}")

        if has_charge:
            ch = atom.formal_charge
            if ch == 1:
                parts.append("+")
            elif ch == -1:
                parts.append("-")
            elif ch > 0:
                parts.append(f"+{ch}")
            else:
                parts.append(str(ch))

        return f"[{''.join(parts)}]"

    def _emit_ring_closure(idx: int, nb: int) -> None:
        """Register and emit a ring closure digit between idx and nb."""
        nonlocal next_ring_digit
        edge = (min(idx, nb), max(idx, nb))
        if edge not in ring_bonds:
            ring_bonds[edge] = next_ring_digit
            next_ring_digit += 1
        digit = ring_bonds[edge]
        bsym = _bond_symbol(mol, idx, nb)
        if digit < 10:
            parts.append(f"{bsym}{digit}")
        else:
            parts.append(f"{bsym}%{digit:02d}")

    def _dfs(idx: int, parent: int | None = None) -> None:
        visited.add(idx)
        parts.append(_atom_str(idx))

        # Classify neighbours into ring closures and tree children
        nbs = _heavy_neighbors(idx)
        parent_consumed = False
        ring_nbs: list[int] = []
        unvisited_nbs: list[int] = []

        for nb in nbs:
            if nb == parent and not parent_consumed:
                # Skip the tree edge we arrived on (consume once)
                parent_consumed = True
                continue
            if nb in visited:
                ring_nbs.append(nb)
            else:
                unvisited_nbs.append(nb)

        # Emit ring closure digits at this atom
        for nb in ring_nbs:
            _emit_ring_closure(idx, nb)

        if not unvisited_nbs:
            return

        # Process children.  All except the last become branches
        # (wrapped in parentheses); the last is the main chain.
        # After each branch, re-check whether later children have
        # been visited (they may have been reached through a ring
        # inside the branch -- in that case emit the ring closure
        # digit at this atom's position).
        remaining = list(unvisited_nbs)
        while remaining:
            # Promote children visited during a sibling branch to ring
            # closures: their partner already emitted the digit, so we
            # must emit the matching digit here.
            still_unvisited: list[int] = []
            for c in remaining:
                if c in visited:
                    edge = (min(idx, c), max(idx, c))
                    if edge in ring_bonds:
                        _emit_ring_closure(idx, c)
                    # else: already handled or no ring bond -- skip
                else:
                    still_unvisited.append(c)
            remaining = still_unvisited

            if not remaining:
                return

            if len(remaining) == 1:
                # Last remaining child: main chain, no parentheses
                child = remaining[0]
                bsym = _bond_symbol(mol, idx, child)
                parts.append(bsym)
                _dfs(child, parent=idx)
                return

            # Branch child (not the last one)
            child = remaining.pop(0)
            if child in visited:
                # Reached through another path; ring closure
                _emit_ring_closure(idx, child)
            else:
                bsym = _bond_symbol(mol, idx, child)
                parts.append(f"({bsym}")
                _dfs(child, parent=idx)
                parts.append(")")

    _dfs(start)

    # Also emit ring closure digits for atoms that have back-edge
    # partners not yet annotated.  (The partner side is handled when
    # the DFS visits that atom -- it adds the digit there.)
    # In the standard algorithm above, both sides are already handled.

    return "".join(parts)


# ===================================================================
# Public API
# ===================================================================

def to_smiles(mol: Molecule) -> str:
    """Convert a Molecule to a SMILES string.

    Hydrogen atoms are omitted (they are implicit in SMILES notation).
    Multi-fragment molecules are joined with ``.`` separators.

    Parameters
    ----------
    mol : Molecule
        A molecule built with the ``molbuilder`` framework.

    Returns
    -------
    str
        A SMILES string representing the molecule.

    Examples
    --------
    >>> from molbuilder.smiles.parser import parse
    >>> mol = parse("CCO")
    >>> to_smiles(mol)  # may return "CCO" or "OCC" depending on canonicalization
    '...'
    """
    if not mol.atoms:
        return ""

    order = _morgan_canonical_order(mol)
    components = _connected_components(mol)

    if not components:
        return ""

    # Sort components for deterministic output (largest first,
    # then by lowest canonical rank of their start atom)
    components.sort(key=lambda c: (-len(c), min(order.get(i, 0) for i in c)))

    fragments = []
    for comp in components:
        smi = _dfs_smiles(mol, order, comp)
        if smi:
            fragments.append(smi)

    return ".".join(fragments)
