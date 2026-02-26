"""SMARTS parser: token stream -> pattern graph.

Parses the output of ``tokenize_smarts()`` into a ``SmartsPattern``
graph consisting of ``SmartsAtom`` nodes and ``SmartsBond`` edges.
The graph can then be matched against ``Molecule`` objects using the
matcher module.

The parser follows the same stack-based branch / ring-closure logic
as the SMILES parser but with SMARTS-specific semantics:

- Default bond between consecutive atoms is *any* bond (not single).
- Atom specifications may include logical combinations (OR/AND).
- Wildcard atoms (``*``) match any element.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from molbuilder.smarts.tokenizer import (
    tokenize_smarts,
    SmartsToken,
    SmartsTokenType,
)


# ===================================================================
# Pattern data classes
# ===================================================================

@dataclass
class SmartsAtom:
    """A single atom specification in a SMARTS pattern.

    Attributes
    ----------
    symbol : str | None
        Element symbol, or None for wildcard (matches any element).
    atomic_number : int | None
        Atomic number constraint (from ``#N``).
    degree : int | None
        Degree (number of explicit connections) constraint.
    valence : int | None
        Total valence (sum of bond orders) constraint.
    hcount : int | None
        Hydrogen count constraint.
    charge : int | None
        Formal charge constraint.  None means unconstrained.
    aromatic : bool | None
        True = must be aromatic, False = must be aliphatic, None = either.
    ring_membership : int | None
        Ring membership count.  -1 means "in at least one ring".
    ring_size : int | None
        Smallest ring size containing this atom.
    is_not : bool
        If True, the entire atom specification is negated.
    recursive_smarts : str | None
        Recursive SMARTS pattern string from ``$(...)`` notation.
    or_atoms : list[SmartsAtom] | None
        List of alternative SmartsAtom specs (OR combination).
    and_atoms : list[SmartsAtom] | None
        List of conjunctive SmartsAtom specs (AND combination).
    """
    symbol: str | None = None
    atomic_number: int | None = None
    degree: int | None = None
    valence: int | None = None
    hcount: int | None = None
    charge: int | None = None
    aromatic: bool | None = None
    ring_membership: int | None = None
    ring_size: int | None = None
    is_not: bool = False
    recursive_smarts: str | None = None
    or_atoms: list | None = None
    and_atoms: list | None = None


@dataclass
class SmartsBond:
    """A bond specification in a SMARTS pattern.

    Attributes
    ----------
    order : float | None
        Required bond order.  None means any bond (default in SMARTS).
        Values: 1.0 (single), 2.0 (double), 3.0 (triple), 1.5 (aromatic).
    is_any : bool
        True for ``~`` (matches any bond).
    is_ring : bool
        True for ``@`` (bond must be in a ring).
    is_not : bool
        True if preceded by ``!`` (negate bond match).
    atom_i : int
        Index of the first atom in the pattern.
    atom_j : int
        Index of the second atom in the pattern.
    """
    order: float | None = None
    is_any: bool = False
    is_ring: bool = False
    is_not: bool = False
    atom_i: int = -1
    atom_j: int = -1


@dataclass
class SmartsPattern:
    """A parsed SMARTS pattern graph.

    Attributes
    ----------
    atoms : list[SmartsAtom]
        The atom specifications in the pattern.
    bonds : list[SmartsBond]
        The bond specifications connecting pattern atoms.
    """
    atoms: list[SmartsAtom] = field(default_factory=list)
    bonds: list[SmartsBond] = field(default_factory=list)

    def neighbor_bonds(self, idx: int) -> list[SmartsBond]:
        """Return all bonds involving pattern atom *idx*."""
        return [b for b in self.bonds if b.atom_i == idx or b.atom_j == idx]

    def neighbors(self, idx: int) -> list[int]:
        """Return indices of atoms bonded to pattern atom *idx*."""
        result: list[int] = []
        for b in self.bonds:
            if b.atom_i == idx:
                result.append(b.atom_j)
            elif b.atom_j == idx:
                result.append(b.atom_i)
        return result


# ===================================================================
# Token -> SmartsAtom / SmartsBond conversion helpers
# ===================================================================

def _token_to_atom(tok: SmartsToken) -> SmartsAtom:
    """Convert a SmartsToken into a SmartsAtom."""
    # Handle OR combinations
    if tok.or_specs is not None:
        or_atoms = [_token_to_atom(sub) for sub in tok.or_specs]
        return SmartsAtom(or_atoms=or_atoms)

    # Handle AND combinations
    if tok.and_specs is not None:
        and_atoms = [_token_to_atom(sub) for sub in tok.and_specs]
        return SmartsAtom(and_atoms=and_atoms)

    # Determine charge: use None (unconstrained) unless explicitly set
    charge_val: int | None = None
    if tok.charge != 0:
        charge_val = tok.charge

    return SmartsAtom(
        symbol=tok.symbol,
        atomic_number=tok.atomic_number,
        degree=tok.degree,
        valence=tok.valence,
        hcount=tok.hcount,
        charge=charge_val,
        aromatic=tok.aromatic,
        ring_membership=tok.ring_membership,
        ring_size=tok.ring_size,
        is_not=tok.is_not,
        recursive_smarts=tok.recursive,
    )


def _bond_token_to_smarts_bond(tok: SmartsToken) -> SmartsBond:
    """Convert a BOND token into a SmartsBond specification."""
    val = tok.value
    if val == "~":
        return SmartsBond(is_any=True)
    if val == "@":
        return SmartsBond(is_ring=True)
    if val == "-":
        return SmartsBond(order=1.0)
    if val == "=":
        return SmartsBond(order=2.0)
    if val == "#":
        return SmartsBond(order=3.0)
    if val == ":":
        return SmartsBond(order=1.5)
    # Fallback: any bond
    return SmartsBond()


# ===================================================================
# Main parser
# ===================================================================

def parse_smarts(smarts: str) -> SmartsPattern:
    """Parse a SMARTS string into a SmartsPattern graph.

    Parameters
    ----------
    smarts : str
        A SMARTS pattern string, e.g. ``"[OH]"``, ``"[CD3](=O)[OD2]"``,
        ``"c1ccccc1"``.

    Returns
    -------
    SmartsPattern
        A pattern graph with atoms and bonds suitable for subgraph
        matching against Molecule objects.

    Raises
    ------
    ValueError
        If the SMARTS string is malformed.
    """
    tokens = tokenize_smarts(smarts)

    atoms: list[SmartsAtom] = []
    bonds: list[SmartsBond] = []

    stack: list[int] = []                       # branch stack
    ring_closures: dict[str, tuple[int, SmartsBond | None]] = {}
    prev: int | None = None                     # most recent atom index
    pending_bond: SmartsBond | None = None       # explicit bond waiting

    for tok in tokens:
        # ---- atom or wildcard ----
        if tok.type in (SmartsTokenType.ATOM, SmartsTokenType.WILDCARD):
            idx = len(atoms)

            if tok.type == SmartsTokenType.WILDCARD:
                atoms.append(SmartsAtom(symbol=None))
            else:
                atoms.append(_token_to_atom(tok))

            # Bond to previous atom
            if prev is not None:
                if pending_bond is not None:
                    bond = pending_bond
                else:
                    # Default SMARTS bond: any bond (not single!)
                    bond = SmartsBond()
                bond.atom_i = prev
                bond.atom_j = idx
                bonds.append(bond)

            pending_bond = None
            prev = idx
            continue

        # ---- bond symbol ----
        if tok.type == SmartsTokenType.BOND:
            pending_bond = _bond_token_to_smarts_bond(tok)
            continue

        # ---- branch open ----
        if tok.type == SmartsTokenType.BRANCH_OPEN:
            if prev is not None:
                stack.append(prev)
            continue

        # ---- branch close ----
        if tok.type == SmartsTokenType.BRANCH_CLOSE:
            if stack:
                prev = stack.pop()
            pending_bond = None
            continue

        # ---- ring closure ----
        if tok.type == SmartsTokenType.RING_DIGIT:
            digit = tok.value
            if digit in ring_closures:
                other_idx, ring_bond = ring_closures.pop(digit)
                if pending_bond is not None:
                    bond = pending_bond
                elif ring_bond is not None:
                    bond = ring_bond
                else:
                    bond = SmartsBond()  # default: any bond
                bond.atom_i = prev if prev is not None else 0
                bond.atom_j = other_idx
                bonds.append(bond)
                pending_bond = None
            else:
                ring_closures[digit] = (
                    prev if prev is not None else 0,
                    pending_bond,
                )
                pending_bond = None
            continue

        # ---- dot (disconnection) ----
        if tok.type == SmartsTokenType.DOT:
            prev = None
            pending_bond = None
            continue

    if ring_closures:
        open_digits = ", ".join(ring_closures.keys())
        raise ValueError(
            f"Unclosed ring closure(s) for digit(s): {open_digits}")

    return SmartsPattern(atoms=atoms, bonds=bonds)
