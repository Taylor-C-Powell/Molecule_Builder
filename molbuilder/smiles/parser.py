"""SMILES parser: tokens -> Molecule with 3D coordinates.

Algorithm:
1. Parse tokens using a stack for branch tracking and a dict for ring closures
2. Build connectivity graph (atoms + bonds)
3. Add implicit hydrogens based on standard valence
4. Assign 3D coordinates via BFS z-matrix placement
"""

from __future__ import annotations

from collections import deque

import numpy as np

from molbuilder.smiles.tokenizer import (
    tokenize, Token, TokenType, DEFAULT_VALENCE, ORGANIC_SUBSET, AROMATIC_ATOMS,
)
from molbuilder.molecule.graph import Molecule, Hybridization
from molbuilder.core.bond_data import bond_length, SP3_ANGLE, SP2_ANGLE, SP_ANGLE
from molbuilder.core.geometry import (
    normalize, place_atom_zmatrix,
)


# ===================================================================
# Bond order from SMILES bond symbol
# ===================================================================

_BOND_ORDER = {
    "-": 1,
    "=": 2,
    "#": 3,
    ":": 1,   # aromatic bond treated as order 1 for connectivity
    "/": 1,   # E/Z bond direction indicator (single bond connectivity)
    "\\": 1,  # E/Z bond direction indicator (single bond connectivity)
}


# ===================================================================
# Internal data structures used during graph construction
# ===================================================================

class _AtomInfo:
    """Lightweight bookkeeping record for an atom during parsing."""

    __slots__ = ("index", "symbol", "aromatic", "bracket",
                 "isotope", "hcount", "charge", "chirality")

    def __init__(self, index: int, symbol: str, aromatic: bool = False,
                 bracket: bool = False, isotope: int | None = None,
                 hcount: int | None = None, charge: int = 0,
                 chirality: str | None = None):
        self.index = index
        self.symbol = symbol
        self.aromatic = aromatic
        self.bracket = bracket
        self.isotope = isotope
        self.hcount = hcount
        self.charge = charge
        self.chirality = chirality


class _BondInfo:
    """Lightweight bookkeeping record for a bond during parsing."""

    __slots__ = ("atom_i", "atom_j", "order")

    def __init__(self, atom_i: int, atom_j: int, order: int = 1):
        self.atom_i = atom_i
        self.atom_j = atom_j
        self.order = order


# ===================================================================
# Graph construction from tokens
# ===================================================================

def _build_graph(tokens: list[Token]) -> tuple[list[_AtomInfo], list[_BondInfo]]:
    """Walk the token list and build atom / bond lists.

    Uses a stack for branch handling and a dictionary for ring closures.

    Returns
    -------
    atoms : list[_AtomInfo]
    bonds : list[_BondInfo]
    """
    atoms: list[_AtomInfo] = []
    bonds: list[_BondInfo] = []

    stack: list[int] = []          # branch stack of atom indices
    ring_closures: dict[str, tuple[int, int]] = {}  # digit -> (atom_idx, bond_order)
    prev: int | None = None        # index of the most recent atom
    pending_bond_order: int | None = None  # explicit bond symbol waiting

    for tok in tokens:
        # ---- atom ----
        if tok.type == TokenType.ATOM:
            idx = len(atoms)
            is_bracket = (tok.hcount is not None or tok.charge != 0
                          or tok.isotope is not None
                          or tok.chirality is not None
                          or (tok.value not in ORGANIC_SUBSET
                              and tok.value not in AROMATIC_ATOMS))
            # Canonical symbol: aromatic lowercase -> titlecase for storage
            symbol = tok.value
            if tok.aromatic and symbol.islower():
                symbol = symbol.capitalize()

            atoms.append(_AtomInfo(
                index=idx,
                symbol=symbol,
                aromatic=tok.aromatic,
                bracket=is_bracket,
                isotope=tok.isotope,
                hcount=tok.hcount,
                charge=tok.charge,
                chirality=tok.chirality,
            ))

            # Bond to previous atom
            if prev is not None:
                order = pending_bond_order if pending_bond_order else 1
                # Aromatic bond default: if both atoms are aromatic and no
                # explicit bond, use order 1 (aromatic bonds are kekulized
                # later or left as single for coordinate purposes).
                bonds.append(_BondInfo(prev, idx, order))
            pending_bond_order = None
            prev = idx
            continue

        # ---- bond symbol ----
        if tok.type == TokenType.BOND:
            pending_bond_order = _BOND_ORDER[tok.value]
            continue

        # ---- branch open ----
        if tok.type == TokenType.BRANCH_OPEN:
            if prev is not None:
                stack.append(prev)
            continue

        # ---- branch close ----
        if tok.type == TokenType.BRANCH_CLOSE:
            if stack:
                prev = stack.pop()
            pending_bond_order = None
            continue

        # ---- ring closure ----
        if tok.type == TokenType.RING_DIGIT:
            digit = tok.value
            if digit in ring_closures:
                # Close the ring
                other_idx, ring_order = ring_closures.pop(digit)
                order = pending_bond_order if pending_bond_order else ring_order
                bonds.append(_BondInfo(prev, other_idx, order))
                pending_bond_order = None
            else:
                # Open a ring
                order = pending_bond_order if pending_bond_order else 1
                ring_closures[digit] = (prev, order)
                pending_bond_order = None
            continue

        # ---- dot (disconnection) ----
        if tok.type == TokenType.DOT:
            prev = None
            pending_bond_order = None
            continue

    if ring_closures:
        open_digits = ", ".join(ring_closures.keys())
        raise ValueError(
            f"Unclosed ring closure(s) for digit(s): {open_digits}")

    return atoms, bonds


# ===================================================================
# Implicit hydrogen addition
# ===================================================================

def _explicit_valence(atom_idx: int, bonds: list[_BondInfo]) -> int:
    """Sum of bond orders touching *atom_idx*."""
    total = 0
    for b in bonds:
        if b.atom_i == atom_idx or b.atom_j == atom_idx:
            total += b.order
    return total


def _add_implicit_hydrogens(
        atoms: list[_AtomInfo],
        bonds: list[_BondInfo],
) -> tuple[list[_AtomInfo], list[_BondInfo]]:
    """Add implicit H atoms to organic-subset atoms.

    Bracket atoms with an explicit ``hcount`` use that count.  Organic-
    subset atoms use ``DEFAULT_VALENCE`` to infer the missing hydrogens.

    Returns the (possibly extended) atoms and bonds lists.
    """
    heavy_count = len(atoms)

    for ai in range(heavy_count):
        atom = atoms[ai]
        ev = _explicit_valence(ai, bonds)

        # Bracket atom: use explicit H count, or 0 if none specified
        # (per OpenSMILES spec, bracket atoms with no H get zero implicit H)
        if atom.bracket:
            n_h = atom.hcount if atom.hcount is not None else 0
        else:
            # Look up default valence for organic subset atoms.
            # Aromatic atoms use the same valence table after
            # kekulization assigns proper bond orders.
            lookup_sym = atom.symbol.lower() if atom.aromatic else atom.symbol
            if lookup_sym not in DEFAULT_VALENCE and atom.symbol not in DEFAULT_VALENCE:
                continue  # unknown atom -- no implicit H
            valences = DEFAULT_VALENCE.get(
                lookup_sym, DEFAULT_VALENCE.get(atom.symbol, []))
            if not valences:
                continue

            # Pick the smallest default valence >= explicit valence
            target = None
            for v in sorted(valences):
                if v >= ev:
                    target = v
                    break
            if target is None:
                target = max(valences)

            # Kekulization has already assigned double bonds to aromatic
            # atoms, so their explicit valence already reflects the pi
            # electron.  No additional reduction needed.
            n_h = max(0, target - ev)

        # Add H atoms
        for _ in range(n_h):
            h_idx = len(atoms)
            atoms.append(_AtomInfo(
                index=h_idx, symbol="H", aromatic=False, bracket=False))
            bonds.append(_BondInfo(ai, h_idx, 1))

    return atoms, bonds


# ===================================================================
# Hybridization determination
# ===================================================================

def _determine_hybridization(
        atom_idx: int,
        atoms: list[_AtomInfo],
        bonds: list[_BondInfo],
) -> Hybridization:
    """Infer hybridization from bond orders around an atom.

    Rules
    -----
    - Any triple bond -> SP
    - Any double bond -> SP2
    - All single bonds -> SP3
    - Aromatic atoms -> SP2
    """
    if atoms[atom_idx].aromatic:
        return Hybridization.SP2

    has_double = False
    for b in bonds:
        if b.atom_i == atom_idx or b.atom_j == atom_idx:
            if b.order == 3:
                return Hybridization.SP
            if b.order == 2:
                has_double = True

    if has_double:
        return Hybridization.SP2
    return Hybridization.SP3


# ===================================================================
# 3D coordinate assignment via BFS
# ===================================================================

def _angle_for_hyb(hyb: Hybridization) -> float:
    """Return the ideal bond angle in degrees for a hybridization."""
    if hyb == Hybridization.SP:
        return SP_ANGLE
    if hyb == Hybridization.SP2:
        return SP2_ANGLE
    return SP3_ANGLE


def _assign_3d_coordinates(mol: Molecule) -> None:
    """Place atoms in 3D using BFS from atom 0.

    Algorithm
    ---------
    1. Place atom 0 at the origin.
    2. BFS outward; for each newly visited atom, use its parent (and
       grandparent if available) for z-matrix placement.
    3. Distribute multiple substituents around each centre at regular
       dihedral intervals based on hybridization.
    """
    n = len(mol.atoms)
    if n == 0:
        return

    # Place first atom at origin
    mol.atoms[0].position = np.array([0.0, 0.0, 0.0])
    placed = {0}

    if n == 1:
        return

    # Build adjacency from molecule bonds
    adj: dict[int, list[tuple[int, int]]] = {i: [] for i in range(n)}
    for b in mol.bonds:
        adj[b.atom_i].append((b.atom_j, b.order))
        adj[b.atom_j].append((b.atom_i, b.order))

    # BFS queue: (atom_index, parent_index, grandparent_index_or_None)
    queue: deque[tuple[int, int | None, int | None]] = deque()

    # Track children scheduled per parent to assign dihedral offsets
    child_counter: dict[int, int] = {}

    # Seed the BFS from atom 0: schedule all neighbours
    for nb_idx, nb_order in adj[0]:
        queue.append((nb_idx, 0, None))

    while queue:
        atom_idx, parent_idx, grandparent_idx = queue.popleft()

        if atom_idx in placed:
            continue

        parent_pos = mol.atoms[parent_idx].position
        parent_hyb = mol.atoms[parent_idx].hybridization
        angle = _angle_for_hyb(parent_hyb) if parent_hyb else SP3_ANGLE

        # Determine bond order for bond length
        b_order = 1
        for nb, bo in adj[parent_idx]:
            if nb == atom_idx:
                b_order = bo
                break

        bl = bond_length(
            mol.atoms[parent_idx].symbol,
            mol.atoms[atom_idx].symbol,
            b_order,
        )

        # Child counter for dihedral offset
        child_num = child_counter.get(parent_idx, 0)
        child_counter[parent_idx] = child_num + 1

        # Dihedral step depends on parent hybridization
        if parent_hyb == Hybridization.SP2:
            dihedral_step = 120.0
        elif parent_hyb == Hybridization.SP:
            dihedral_step = 180.0
        else:
            dihedral_step = 120.0  # tetrahedral uses ~120 deg between projections

        dihedral = dihedral_step * child_num

        if grandparent_idx is not None and grandparent_idx in placed:
            # Normal z-matrix placement
            gp_pos = mol.atoms[grandparent_idx].position
            pos = place_atom_zmatrix(
                parent_pos, gp_pos,
                _dihedral_ref_pos(mol, parent_idx, grandparent_idx, placed),
                bl, angle, dihedral,
            )
        elif len(placed) == 1:
            # Second atom: place along +z
            pos = parent_pos + np.array([0.0, 0.0, bl])
        else:
            # No grandparent yet -- use a synthetic reference
            ref_pos = parent_pos + np.array([0.0, 0.0, -1.0])
            synth_k = ref_pos + np.array([0.0, 1.0, 0.0])
            pos = place_atom_zmatrix(
                parent_pos, ref_pos, synth_k,
                bl, angle, dihedral,
            )

        mol.atoms[atom_idx].position = pos
        placed.add(atom_idx)

        # Enqueue unvisited neighbours
        for nb_idx, nb_order in adj[atom_idx]:
            if nb_idx not in placed:
                queue.append((nb_idx, atom_idx, parent_idx))


def _dihedral_ref_pos(
        mol: Molecule,
        parent_idx: int,
        grandparent_idx: int,
        placed: set[int],
) -> np.ndarray:
    """Find a third reference position for z-matrix dihedral.

    Looks for a placed neighbour of the grandparent that is not the
    parent.  Falls back to a synthetic offset if none is found.
    """
    gp_pos = mol.atoms[grandparent_idx].position
    for nb in mol.neighbors(grandparent_idx):
        if nb != parent_idx and nb in placed:
            return mol.atoms[nb].position

    # Synthetic fallback: offset perpendicular to parent-grandparent axis
    axis = mol.atoms[parent_idx].position - gp_pos
    perp = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(normalize(axis), perp)) > 0.9:
        perp = np.array([0.0, 1.0, 0.0])
    return gp_pos + np.cross(axis, perp) * 0.5


# ===================================================================
# Public API
# ===================================================================

def parse(smiles: str) -> Molecule:
    """Parse a SMILES string and return a Molecule with 3D coordinates.

    Parameters
    ----------
    smiles : str
        A SMILES string, e.g. ``"CCO"`` (ethanol), ``"c1ccccc1"``
        (benzene), ``"CC(=O)O"`` (acetic acid).

    Returns
    -------
    Molecule
        A fully constructed molecule with atoms, bonds, hybridization,
        and approximate 3D coordinates.

    Raises
    ------
    ValueError
        If the SMILES string is invalid.

    Examples
    --------
    >>> mol = parse("C")
    >>> len([a for a in mol.atoms if a.symbol == "C"])
    1
    >>> len([a for a in mol.atoms if a.symbol == "H"])
    4
    """
    tokens = tokenize(smiles)
    atoms, bonds = _build_graph(tokens)

    # Kekulize: assign alternating single/double bonds to aromatic systems
    from molbuilder.smiles.aromatic import kekulize
    kekulize(atoms, bonds)

    atoms, bonds = _add_implicit_hydrogens(atoms, bonds)

    # Create molecule
    mol = Molecule(name=smiles)

    # First pass: add atoms with placeholder positions and hybridization
    for ai in atoms:
        hyb = _determine_hybridization(ai.index, atoms, bonds)
        mol.add_atom(
            symbol=ai.symbol,
            position=[0.0, 0.0, 0.0],
            hybridization=hyb,
            chirality=ai.chirality,
            isotope=ai.isotope,
            formal_charge=ai.charge,
        )

    # Add bonds
    for bi in bonds:
        # Determine if rotatable (single bonds between heavy atoms)
        rot = (bi.order == 1
               and atoms[bi.atom_i].symbol != "H"
               and atoms[bi.atom_j].symbol != "H")
        mol.add_bond(bi.atom_i, bi.atom_j, order=bi.order, rotatable=rot)

    # Assign 3D coordinates via BFS
    _assign_3d_coordinates(mol)

    # Correct ring geometry (planar aromatic rings, chair cyclohexane, etc.)
    from molbuilder.smiles.ring_geometry import correct_ring_geometry
    correct_ring_geometry(mol)

    return mol
