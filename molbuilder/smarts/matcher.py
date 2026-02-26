"""VF2-style subgraph isomorphism matcher for SMARTS patterns.

Matches a ``SmartsPattern`` against a ``Molecule`` by backtracking
over candidate mappings.  Atom and bond constraints from the SMARTS
specification are enforced at each step.

Public API
----------
find_matches(mol, pattern) -> list[dict[int, int]]
    Find all subgraph matches (pattern_idx -> mol_idx).
has_match(mol, pattern) -> bool
    Return True on first match found.
detect_by_smarts(mol, smarts_str, group_name) -> list[FunctionalGroup]
    Convenience wrapper returning FunctionalGroup instances.
"""

from __future__ import annotations

from molbuilder.molecule.graph import Molecule, Hybridization
from molbuilder.smarts.parser import SmartsAtom, SmartsBond, SmartsPattern, parse_smarts
from molbuilder.core.elements import SYMBOL_TO_Z


# ===================================================================
# Helpers for querying molecule atom properties
# ===================================================================

# Standard valences for implicit hydrogen inference (mirrors the data
# in functional_group_detect.py so that the SMARTS matcher does not
# depend on the SMILES subpackage).
_STANDARD_VALENCE: dict[str, list[int]] = {
    "B": [3], "C": [4], "N": [3, 5], "O": [2], "P": [3, 5],
    "S": [2, 4, 6], "F": [1], "Cl": [1], "Br": [1], "I": [1],
}


def _h_count(mol: Molecule, idx: int) -> int:
    """Total hydrogen count on atom *idx* (explicit + implicit).

    Explicit H atoms are counted from the neighbour list.  Implicit H
    atoms are inferred from standard valence rules when the molecule
    has fewer explicit neighbours than expected.
    """
    explicit_h = sum(
        1 for nb in mol.neighbors(idx) if mol.atoms[nb].symbol == "H"
    )
    if explicit_h > 0:
        return explicit_h

    sym = mol.atoms[idx].symbol
    valences = _STANDARD_VALENCE.get(sym)
    if valences is None:
        return 0
    bond_order_sum = _sum_bond_orders(mol, idx)
    for v in valences:
        implicit = v - bond_order_sum
        if implicit >= 0:
            return implicit
    return 0


def _sum_bond_orders(mol: Molecule, idx: int) -> int:
    """Sum of bond orders for all bonds on atom *idx*."""
    total = 0
    for nb in mol.neighbors(idx):
        bond = mol.get_bond(idx, nb)
        if bond is not None:
            total += bond.order
    return total


def _atom_is_aromatic(mol: Molecule, idx: int) -> bool:
    """Return True if atom *idx* is considered aromatic.

    An atom is aromatic if it is SP2-hybridized and participates in at
    least one ring bond.
    """
    atom = mol.atoms[idx]
    if atom.hybridization != Hybridization.SP2:
        return False
    # Check if any bond from this atom is a ring bond
    for nb in mol.neighbors(idx):
        if mol.is_in_ring(idx, nb):
            return True
    return False


def _count_ring_memberships(mol: Molecule, idx: int) -> int:
    """Count the number of distinct rings containing atom *idx*.

    Uses a simple approach: count ring bonds incident on the atom.
    Each ring bond indicates membership in at least one ring.  For a
    more accurate count we would need full SSSR computation, but this
    gives a reasonable approximation (0 = not in ring, >0 = in ring).
    """
    count = 0
    for nb in mol.neighbors(idx):
        if mol.is_in_ring(idx, nb):
            count += 1
    # Each ring contributes 2 ring bonds to any member atom, so divide
    # by 2 and round up to approximate ring count.
    return (count + 1) // 2 if count > 0 else 0


def _smallest_ring_size(mol: Molecule, idx: int) -> int:
    """Find the size of the smallest ring containing atom *idx*.

    Uses BFS from *idx* through ring bonds.  Returns 0 if the atom
    is not in any ring.
    """
    # Quick check: atom must be in at least one ring
    ring_nbs = [nb for nb in mol.neighbors(idx) if mol.is_in_ring(idx, nb)]
    if not ring_nbs:
        return 0

    # BFS to find shortest cycle through idx
    best = len(mol.atoms) + 1
    for start_nb in ring_nbs:
        # BFS from start_nb, trying to reach idx without going back directly
        visited: dict[int, int] = {start_nb: 1}
        frontier = [start_nb]
        found = False
        while frontier and not found:
            next_frontier: list[int] = []
            for node in frontier:
                for nb in mol.neighbors(node):
                    if nb == idx and visited[node] >= 2:
                        # Found a cycle of length visited[node] + 1
                        ring_len = visited[node] + 1
                        if ring_len < best:
                            best = ring_len
                        found = True
                        break
                    if nb == idx:
                        continue  # don't go back to start too early
                    if nb not in visited:
                        visited[nb] = visited[node] + 1
                        next_frontier.append(nb)
                if found:
                    break
            frontier = next_frontier

    return best if best <= len(mol.atoms) else 0


# ===================================================================
# Atom matching
# ===================================================================

def _atom_matches(
    mol_atom_idx: int,
    mol: Molecule,
    smarts_atom: SmartsAtom,
    pattern: SmartsPattern,
) -> bool:
    """Check if molecule atom at *mol_atom_idx* satisfies *smarts_atom*.

    Parameters
    ----------
    mol_atom_idx : int
        Index of the candidate atom in the molecule.
    mol : Molecule
        The molecule being searched.
    smarts_atom : SmartsAtom
        The SMARTS atom specification to match against.
    pattern : SmartsPattern
        The full pattern (needed for recursive SMARTS context).

    Returns
    -------
    bool
    """
    # Handle NOT
    if smarts_atom.is_not:
        # Create a non-negated copy to test, then invert
        non_neg = SmartsAtom(
            symbol=smarts_atom.symbol,
            atomic_number=smarts_atom.atomic_number,
            degree=smarts_atom.degree,
            valence=smarts_atom.valence,
            hcount=smarts_atom.hcount,
            charge=smarts_atom.charge,
            aromatic=smarts_atom.aromatic,
            ring_membership=smarts_atom.ring_membership,
            ring_size=smarts_atom.ring_size,
            is_not=False,
            recursive_smarts=smarts_atom.recursive_smarts,
            or_atoms=smarts_atom.or_atoms,
            and_atoms=smarts_atom.and_atoms,
        )
        return not _atom_matches(mol_atom_idx, mol, non_neg, pattern)

    # Handle OR: any alternative must match
    if smarts_atom.or_atoms is not None:
        return any(
            _atom_matches(mol_atom_idx, mol, alt, pattern)
            for alt in smarts_atom.or_atoms
        )

    # Handle AND: all specs must match
    if smarts_atom.and_atoms is not None:
        return all(
            _atom_matches(mol_atom_idx, mol, spec, pattern)
            for spec in smarts_atom.and_atoms
        )

    mol_atom = mol.atoms[mol_atom_idx]

    # Symbol check (None = wildcard, matches any)
    if smarts_atom.symbol is not None:
        if mol_atom.symbol != smarts_atom.symbol:
            return False

    # Atomic number check
    if smarts_atom.atomic_number is not None:
        mol_z = SYMBOL_TO_Z.get(mol_atom.symbol, 0)
        if mol_z != smarts_atom.atomic_number:
            return False

    # Degree check (number of explicit heavy-atom neighbors)
    if smarts_atom.degree is not None:
        mol_degree = len(mol.neighbors(mol_atom_idx))
        if mol_degree != smarts_atom.degree:
            return False

    # Valence check (sum of bond orders)
    if smarts_atom.valence is not None:
        mol_valence = _sum_bond_orders(mol, mol_atom_idx)
        if mol_valence != smarts_atom.valence:
            return False

    # H-count check
    if smarts_atom.hcount is not None:
        mol_hcount = _h_count(mol, mol_atom_idx)
        if mol_hcount != smarts_atom.hcount:
            return False

    # Charge check
    if smarts_atom.charge is not None:
        if mol_atom.formal_charge != smarts_atom.charge:
            return False

    # Aromatic check
    if smarts_atom.aromatic is not None:
        is_arom = _atom_is_aromatic(mol, mol_atom_idx)
        if smarts_atom.aromatic and not is_arom:
            return False
        if not smarts_atom.aromatic and is_arom:
            return False

    # Ring membership check
    if smarts_atom.ring_membership is not None:
        ring_count = _count_ring_memberships(mol, mol_atom_idx)
        if smarts_atom.ring_membership == -1:
            # Bare R: must be in at least one ring
            if ring_count == 0:
                return False
        elif smarts_atom.ring_membership == 0:
            # R0: must NOT be in a ring
            if ring_count != 0:
                return False
        else:
            if ring_count != smarts_atom.ring_membership:
                return False

    # Ring size check
    if smarts_atom.ring_size is not None:
        srs = _smallest_ring_size(mol, mol_atom_idx)
        if smarts_atom.ring_size == -1:
            # Bare r: must be in at least one ring
            if srs == 0:
                return False
        else:
            if srs != smarts_atom.ring_size:
                return False

    # Recursive SMARTS check
    if smarts_atom.recursive_smarts is not None:
        rec_pattern = parse_smarts(smarts_atom.recursive_smarts)
        # The recursive pattern must match with mol_atom_idx as one of
        # the matched atoms (specifically the first atom of the pattern).
        rec_matches = find_matches(mol, rec_pattern)
        if not any(m.get(0) == mol_atom_idx for m in rec_matches):
            return False

    return True


# ===================================================================
# Bond matching
# ===================================================================

def _bond_matches(
    mol: Molecule,
    mol_i: int,
    mol_j: int,
    smarts_bond: SmartsBond,
) -> bool:
    """Check if the molecule bond between *mol_i* and *mol_j* matches *smarts_bond*.

    Parameters
    ----------
    mol : Molecule
        The molecule being searched.
    mol_i, mol_j : int
        Atom indices in the molecule.
    smarts_bond : SmartsBond
        The SMARTS bond specification.

    Returns
    -------
    bool
    """
    mol_bond = mol.get_bond(mol_i, mol_j)
    if mol_bond is None:
        return False

    # NOT bond
    if smarts_bond.is_not:
        non_neg = SmartsBond(
            order=smarts_bond.order,
            is_any=smarts_bond.is_any,
            is_ring=smarts_bond.is_ring,
            is_not=False,
            atom_i=smarts_bond.atom_i,
            atom_j=smarts_bond.atom_j,
        )
        return not _bond_matches(mol, mol_i, mol_j, non_neg)

    # Any bond (~)
    if smarts_bond.is_any:
        return True

    # Ring bond (@)
    if smarts_bond.is_ring:
        return mol.is_in_ring(mol_i, mol_j)

    # Specific bond order
    if smarts_bond.order is not None:
        mol_order = float(mol_bond.order)
        if abs(mol_order - smarts_bond.order) > 0.1:
            return False
        return True

    # Default (order is None, is_any False, is_ring False): matches any bond
    return True


# ===================================================================
# VF2-style subgraph matching
# ===================================================================

def find_matches(
    mol: Molecule,
    pattern: SmartsPattern,
) -> list[dict[int, int]]:
    """Find all subgraph matches of *pattern* in *mol*.

    Parameters
    ----------
    mol : Molecule
        The molecule to search within.
    pattern : SmartsPattern
        The SMARTS pattern to match.

    Returns
    -------
    list[dict[int, int]]
        A list of mappings, each mapping pattern atom index to
        molecule atom index.  Empty list if no matches.
    """
    n_pattern = len(pattern.atoms)
    if n_pattern == 0:
        return []

    n_mol = len(mol.atoms)
    if n_mol == 0:
        return []

    # Pre-build adjacency for the pattern
    pattern_adj: dict[int, list[int]] = {i: [] for i in range(n_pattern)}
    for b in pattern.bonds:
        pattern_adj[b.atom_i].append(b.atom_j)
        pattern_adj[b.atom_j].append(b.atom_i)

    # Pre-filter: for each pattern atom, compute candidate mol atoms
    candidates: list[list[int]] = []
    for pi in range(n_pattern):
        p_atom = pattern.atoms[pi]
        cands: list[int] = []
        for mi in range(n_mol):
            if _atom_matches(mi, mol, p_atom, pattern):
                cands.append(mi)
        candidates.append(cands)
        # Early termination: if any pattern atom has no candidates, no match
        if not cands:
            return []

    results: list[dict[int, int]] = []
    mapping: dict[int, int] = {}     # pattern_idx -> mol_idx
    used_mol: set[int] = set()       # mol atoms already mapped

    def _get_bond_between(pi_a: int, pi_b: int) -> SmartsBond | None:
        """Return the SmartsBond between pattern atoms pi_a and pi_b."""
        for b in pattern.bonds:
            if (b.atom_i == pi_a and b.atom_j == pi_b) or \
               (b.atom_i == pi_b and b.atom_j == pi_a):
                return b
        return None

    def _backtrack(depth: int) -> None:
        """Recursively extend the mapping one pattern atom at a time."""
        if depth == n_pattern:
            # Complete mapping found
            results.append(dict(mapping))
            return

        p_idx = depth  # map pattern atoms in order 0, 1, 2, ...

        for m_idx in candidates[p_idx]:
            if m_idx in used_mol:
                continue

            # Check bond constraints with already-mapped neighbors
            ok = True
            for p_nb in pattern_adj[p_idx]:
                if p_nb not in mapping:
                    continue  # neighbor not yet mapped
                m_nb = mapping[p_nb]
                bond_spec = _get_bond_between(p_idx, p_nb)
                if bond_spec is not None:
                    if not _bond_matches(mol, m_idx, m_nb, bond_spec):
                        ok = False
                        break
                else:
                    # Implicit bond (should not happen but check connectivity)
                    if mol.get_bond(m_idx, m_nb) is None:
                        ok = False
                        break

            if not ok:
                continue

            # Extend mapping
            mapping[p_idx] = m_idx
            used_mol.add(m_idx)

            _backtrack(depth + 1)

            # Undo
            del mapping[p_idx]
            used_mol.discard(m_idx)

    _backtrack(0)
    return results


def has_match(mol: Molecule, pattern: SmartsPattern) -> bool:
    """Return True if *pattern* matches anywhere in *mol*.

    This is an optimized version of ``find_matches`` that stops on
    the first match found.

    Parameters
    ----------
    mol : Molecule
        The molecule to search within.
    pattern : SmartsPattern
        The SMARTS pattern to match.

    Returns
    -------
    bool
    """
    n_pattern = len(pattern.atoms)
    if n_pattern == 0:
        return False

    n_mol = len(mol.atoms)
    if n_mol == 0:
        return False

    # Pre-build adjacency for the pattern
    pattern_adj: dict[int, list[int]] = {i: [] for i in range(n_pattern)}
    for b in pattern.bonds:
        pattern_adj[b.atom_i].append(b.atom_j)
        pattern_adj[b.atom_j].append(b.atom_i)

    # Pre-filter candidates
    candidates: list[list[int]] = []
    for pi in range(n_pattern):
        p_atom = pattern.atoms[pi]
        cands = [
            mi for mi in range(n_mol)
            if _atom_matches(mi, mol, p_atom, pattern)
        ]
        candidates.append(cands)
        if not cands:
            return False

    mapping: dict[int, int] = {}
    used_mol: set[int] = set()

    def _get_bond_between(pi_a: int, pi_b: int) -> SmartsBond | None:
        for b in pattern.bonds:
            if (b.atom_i == pi_a and b.atom_j == pi_b) or \
               (b.atom_i == pi_b and b.atom_j == pi_a):
                return b
        return None

    def _backtrack(depth: int) -> bool:
        if depth == n_pattern:
            return True

        p_idx = depth

        for m_idx in candidates[p_idx]:
            if m_idx in used_mol:
                continue

            ok = True
            for p_nb in pattern_adj[p_idx]:
                if p_nb not in mapping:
                    continue
                m_nb = mapping[p_nb]
                bond_spec = _get_bond_between(p_idx, p_nb)
                if bond_spec is not None:
                    if not _bond_matches(mol, m_idx, m_nb, bond_spec):
                        ok = False
                        break
                else:
                    if mol.get_bond(m_idx, m_nb) is None:
                        ok = False
                        break

            if not ok:
                continue

            mapping[p_idx] = m_idx
            used_mol.add(m_idx)

            if _backtrack(depth + 1):
                return True

            del mapping[p_idx]
            used_mol.discard(m_idx)

        return False

    return _backtrack(0)


# ===================================================================
# Convenience function
# ===================================================================

def detect_by_smarts(
    mol: Molecule,
    smarts_str: str,
    group_name: str = "custom",
) -> list:
    """Detect functional groups by SMARTS pattern.

    This is a convenience wrapper that parses a SMARTS string, matches
    it against *mol*, and returns a list of ``FunctionalGroup`` objects
    with the matched atom indices.

    Parameters
    ----------
    mol : Molecule
        The molecule to search.
    smarts_str : str
        A SMARTS pattern string.
    group_name : str
        Name to assign to detected groups (default ``"custom"``).

    Returns
    -------
    list[FunctionalGroup]
        One FunctionalGroup per match found.
    """
    from molbuilder.reactions.functional_group_detect import FunctionalGroup

    pattern = parse_smarts(smarts_str)
    matches = find_matches(mol, pattern)

    results: list = []
    for match in matches:
        atom_indices = sorted(match.values())
        center = atom_indices[0] if atom_indices else -1
        results.append(FunctionalGroup(
            name=group_name,
            smarts_like=smarts_str,
            atoms=atom_indices,
            center=center,
        ))
    return results
