"""
Lewis Structure Determination

Parses a molecular formula, identifies the central atom, distributes
bonding pairs and lone pairs, and checks octet satisfaction.

Handles:
    - Single, double, and triple bonds
    - Expanded octets (period 3+ central atoms)
    - Incomplete octets (H, Be, B)
    - Polyatomic ions via charge parameter
"""

import re
from dataclasses import dataclass
from bohr_model import SYMBOL_TO_Z
from element_data import (
    electronegativity,
    target_electrons,
    can_expand_octet,
    PAULING_ELECTRONEGATIVITY,
)


# ===================================================================
# Formula parsing
# ===================================================================

def parse_formula(formula: str) -> list[str]:
    """Parse a molecular formula string into a list of element symbols.

    Examples
    --------
    >>> parse_formula('H2O')
    ['H', 'H', 'O']
    >>> parse_formula('CH4')
    ['C', 'H', 'H', 'H', 'H']
    >>> parse_formula('SF6')
    ['S', 'F', 'F', 'F', 'F', 'F', 'F']
    """
    tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    atoms = []
    for symbol, count_str in tokens:
        if not symbol:
            continue
        if symbol not in SYMBOL_TO_Z:
            raise ValueError(f"Unknown element symbol: {symbol}")
        count = int(count_str) if count_str else 1
        atoms.extend([symbol] * count)
    return atoms


# ===================================================================
# Valence electron counting
# ===================================================================

# Main-group valence electron lookup by group number.
# VSEPR is primarily for main-group elements.
_GROUP_VALENCE = {
    1: 1, 2: 2,
    13: 3, 14: 4, 15: 5, 16: 6, 17: 7, 18: 8,
}

# Map atomic number -> main-group number for common elements
_Z_TO_GROUP = {
    1: 1, 2: 18,
    3: 1, 4: 2, 5: 13, 6: 14, 7: 15, 8: 16, 9: 17, 10: 18,
    11: 1, 12: 2, 13: 13, 14: 14, 15: 15, 16: 16, 17: 17, 18: 18,
    19: 1, 20: 2,
    31: 13, 32: 14, 33: 15, 34: 16, 35: 17, 36: 18,
    37: 1, 38: 2,
    49: 13, 50: 14, 51: 15, 52: 16, 53: 17, 54: 18,
    55: 1, 56: 2,
    81: 13, 82: 14, 83: 15, 84: 16, 85: 17, 86: 18,
}


def get_valence_electrons(symbol: str) -> int:
    """Get the number of VSEPR-relevant valence electrons for an element.

    Uses group-number lookup for main-group elements.
    Falls back to quantum_model.QuantumAtom for transition metals.
    """
    z = SYMBOL_TO_Z.get(symbol)
    if z is None:
        raise ValueError(f"Unknown element: {symbol}")

    group = _Z_TO_GROUP.get(z)
    if group is not None:
        ve = _GROUP_VALENCE.get(group, 0)
        # Helium special case
        if symbol == "He":
            return 2
        return ve

    # Fallback: use QuantumAtom
    from quantum_model import QuantumAtom
    atom = QuantumAtom(z)
    return atom.valence_electrons


# ===================================================================
# Central atom identification
# ===================================================================

def identify_central_atom(atoms: list[str]) -> int:
    """Identify the central atom index using VSEPR conventions.

    Rules (in priority order):
        1. Hydrogen is never the central atom.
        2. The least electronegative non-H atom is central.
        3. If tied, the atom appearing fewest times is central
           (the unique atom tends to be central).
        4. For diatomics with no non-H atoms, pick index 0.

    Returns the index into the atoms list.
    """
    if len(atoms) <= 2:
        # Diatomic: pick the first non-H, or index 0
        for i, sym in enumerate(atoms):
            if sym != "H":
                return i
        return 0

    # Collect non-H candidates
    candidates = [(i, sym) for i, sym in enumerate(atoms) if sym != "H"]

    if not candidates:
        return 0  # all hydrogen (unusual)

    if len(candidates) == 1:
        return candidates[0][0]

    # Sort by electronegativity (ascending), then by frequency (ascending)
    from collections import Counter
    freq = Counter(atoms)
    candidates.sort(key=lambda t: (electronegativity(t[1]), freq[t[1]]))
    return candidates[0][0]


# ===================================================================
# Data classes
# ===================================================================

@dataclass
class Bond:
    """A bond between two atom indices."""
    atom_a: int
    atom_b: int
    order: int = 1  # 1=single, 2=double, 3=triple

    def __repr__(self):
        order_str = {1: "-", 2: "=", 3: "#"}.get(self.order, "?")
        return f"Bond({self.atom_a}{order_str}{self.atom_b})"


@dataclass
class LonePair:
    """A lone pair on a specific atom."""
    atom_index: int

    def __repr__(self):
        return f"LP({self.atom_index})"


# ===================================================================
# Lewis Structure
# ===================================================================

class LewisStructure:
    """Determines and stores the Lewis structure of a molecule.

    Parameters
    ----------
    formula : str
        Molecular formula (e.g., 'H2O', 'CH4', 'CO2').
    charge : int
        Net charge (0 for neutral, +1 for cation, -1 for anion, etc.).
    """

    def __init__(self, formula: str, charge: int = 0):
        self.formula = formula
        self.charge = charge
        self.atoms = parse_formula(formula)

        if len(self.atoms) < 1:
            raise ValueError(f"Formula '{formula}' produced no atoms")

        self.central_index = identify_central_atom(self.atoms)
        self.central_symbol = self.atoms[self.central_index]
        self.terminal_indices = [
            i for i in range(len(self.atoms)) if i != self.central_index
        ]

        self.bonds: list[Bond] = []
        self.lone_pairs: list[LonePair] = []
        self.total_valence_electrons = 0

        self._solve()

    # ------ solver ------

    def _solve(self):
        """Main Lewis structure algorithm."""
        # Step 1: count total valence electrons
        self.total_valence_electrons = (
            sum(get_valence_electrons(sym) for sym in self.atoms) - self.charge
        )

        remaining = self.total_valence_electrons

        # Step 2: place single bonds from central to each terminal
        for ti in self.terminal_indices:
            self.bonds.append(Bond(self.central_index, ti, order=1))
            remaining -= 2

        if remaining < 0:
            # Electron-deficient molecule (e.g., H2 with charge)
            remaining = 0

        # Step 3: distribute lone pairs to terminal atoms first
        for ti in self.terminal_indices:
            sym = self.atoms[ti]
            target = target_electrons(sym)
            # Terminal already has 2 electrons from its bond
            needed = target - 2
            pairs_to_add = needed // 2
            for _ in range(pairs_to_add):
                if remaining >= 2:
                    self.lone_pairs.append(LonePair(ti))
                    remaining -= 2

        # Step 4: place remaining electrons on central atom as lone pairs
        while remaining >= 2:
            self.lone_pairs.append(LonePair(self.central_index))
            remaining -= 2

        # Step 5: check if central atom needs more electrons (form multiple bonds)
        central_target = target_electrons(self.central_symbol)
        central_electrons = self._electrons_around(self.central_index)

        # Only promote to multiple bonds if central atom is NOT allowed
        # to have an expanded octet, or if it actually needs more electrons
        if not can_expand_octet(self.central_symbol):
            while central_electrons < central_target:
                promoted = False
                # Sort bonds by order (ascending) so we promote the
                # lowest-order bond first -- distributes evenly
                for bond in sorted(self.bonds, key=lambda b: b.order):
                    if bond.order >= 3:
                        continue
                    # Find the terminal atom of this bond
                    ti = bond.atom_b if bond.atom_a == self.central_index else bond.atom_a
                    # Check if terminal has a lone pair we can convert
                    lp_idx = self._find_lone_pair_on(ti)
                    if lp_idx is not None:
                        self.lone_pairs.pop(lp_idx)
                        bond.order += 1
                        central_electrons += 2
                        promoted = True
                        break
                if not promoted:
                    break  # no more promotions possible

    def _electrons_around(self, atom_index: int) -> int:
        """Count electrons around a given atom (bonding + lone pairs)."""
        count = 0
        for bond in self.bonds:
            if bond.atom_a == atom_index or bond.atom_b == atom_index:
                count += bond.order * 2
        for lp in self.lone_pairs:
            if lp.atom_index == atom_index:
                count += 2
        return count

    def _find_lone_pair_on(self, atom_index: int):
        """Find the index of a lone pair on atom_index, or None."""
        for i, lp in enumerate(self.lone_pairs):
            if lp.atom_index == atom_index:
                return i
        return None

    # ------ query methods ------

    def bonding_pairs_on_central(self) -> int:
        """Number of bonding groups (sigma bonds) on the central atom."""
        count = 0
        for bond in self.bonds:
            if bond.atom_a == self.central_index or bond.atom_b == self.central_index:
                count += 1
        return count

    def lone_pairs_on_central(self) -> int:
        """Number of lone pairs on the central atom."""
        return sum(1 for lp in self.lone_pairs if lp.atom_index == self.central_index)

    def steric_number(self) -> int:
        """Steric number = bonding groups + lone pairs on central."""
        return self.bonding_pairs_on_central() + self.lone_pairs_on_central()

    def bond_order_to(self, terminal_index: int) -> int:
        """Get bond order between central atom and a terminal atom."""
        for bond in self.bonds:
            a, b = bond.atom_a, bond.atom_b
            if (a == self.central_index and b == terminal_index) or \
               (b == self.central_index and a == terminal_index):
                return bond.order
        return 0

    # ------ display ------

    def __repr__(self):
        bp = self.bonding_pairs_on_central()
        lp = self.lone_pairs_on_central()
        return f"LewisStructure({self.formula}, central={self.central_symbol}, X={bp}, E={lp})"

    def summary(self) -> str:
        """Return an ASCII summary of the Lewis structure."""
        lines = [
            f"{'='*55}",
            f"  Lewis Structure: {self.formula}",
            f"  Charge: {self.charge:+d}",
            f"  Total valence electrons: {self.total_valence_electrons}",
            f"  Central atom: {self.central_symbol} (index {self.central_index})",
            f"{'='*55}",
        ]

        # Bonds
        lines.append("  Bonds:")
        for bond in self.bonds:
            sym_a = self.atoms[bond.atom_a]
            sym_b = self.atoms[bond.atom_b]
            order_label = {1: "single", 2: "double", 3: "triple"}.get(bond.order, "?")
            order_sym = {1: "-", 2: "=", 3: "#"}.get(bond.order, "?")
            lines.append(f"    {sym_a}{order_sym}{sym_b}  ({order_label})")

        # Lone pairs per atom
        lines.append("  Lone pairs:")
        for i, sym in enumerate(self.atoms):
            count = sum(1 for lp in self.lone_pairs if lp.atom_index == i)
            if count > 0:
                lines.append(f"    {sym} (index {i}): {count} lone pair(s)")

        lines.append(f"  Steric number: {self.steric_number()}")
        lines.append(f"{'='*55}")
        return "\n".join(lines)


# ===================================================================
# Main -- demo
# ===================================================================

if __name__ == "__main__":
    print("Lewis Structure Determination\n")

    test_cases = [
        ("H2",   0),
        ("H2O",  0),
        ("NH3",  0),
        ("CH4",  0),
        ("CO2",  0),
        ("BF3",  0),
        ("PCl5", 0),
        ("SF6",  0),
        ("SF4",  0),
        ("ClF3", 0),
        ("XeF2", 0),
        ("XeF4", 0),
        ("IF5",  0),
        ("NH4",  1),   # ammonium ion
    ]

    for formula, charge in test_cases:
        lewis = LewisStructure(formula, charge)
        print(lewis.summary())
        print()
