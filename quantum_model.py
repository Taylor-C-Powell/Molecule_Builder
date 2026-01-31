"""
Quantum Mechanical Model of the Atom

Represents atoms using the quantum mechanical framework:
    - Four quantum numbers (n, l, m_l, m_s) per electron
    - Electron configuration via Aufbau principle, Hund's rule, Pauli exclusion
    - Slater's rules for effective nuclear charge (Z_eff)
    - Approximate orbital energies for multi-electron atoms
    - Known configuration exceptions for d- and f-block elements
"""

import math
from dataclasses import dataclass, field
from quantum_wavefunctions import (
    SUBSHELL_LETTER,
    SUBSHELL_NUMBER,
    RYDBERG_ENERGY_EV,
    EV_TO_JOULES,
    BOHR_RADIUS,
    orbital_angular_momentum,
    angular_momentum_z,
)

# Import element data from the existing bohr model module
from bohr_model import ELEMENTS, SYMBOL_TO_Z


# ===================================================================
# Quantum number container
# ===================================================================

@dataclass(frozen=True, order=True)
class QuantumState:
    """A single-electron quantum state defined by four quantum numbers.

    n  : principal quantum number          (1, 2, 3, ...)
    l  : orbital angular momentum number   (0 .. n-1)
    ml : magnetic quantum number           (-l .. +l)
    ms : spin magnetic quantum number      (+0.5 or -0.5)
    """
    n: int
    l: int
    ml: int
    ms: float

    def __post_init__(self):
        if self.n < 1:
            raise ValueError(f"n must be >= 1, got {self.n}")
        if not (0 <= self.l < self.n):
            raise ValueError(f"l must be 0..{self.n-1}, got {self.l}")
        if not (-self.l <= self.ml <= self.l):
            raise ValueError(f"ml must be {-self.l}..{self.l}, got {self.ml}")
        if self.ms not in (0.5, -0.5):
            raise ValueError(f"ms must be +0.5 or -0.5, got {self.ms}")

    @property
    def subshell_label(self) -> str:
        return f"{self.n}{SUBSHELL_LETTER.get(self.l, '?')}"

    def __repr__(self):
        sign = "+" if self.ms > 0 else "-"
        return f"({self.n},{self.l},{self.ml},{sign}1/2)"


# ===================================================================
# Subshell
# ===================================================================

@dataclass
class Subshell:
    """A subshell (n, l) that can hold up to 2*(2l+1) electrons."""
    n: int
    l: int
    electron_count: int = 0

    @property
    def capacity(self) -> int:
        return 2 * (2 * self.l + 1)

    @property
    def is_full(self) -> bool:
        return self.electron_count >= self.capacity

    @property
    def is_half_filled(self) -> bool:
        return self.electron_count == (2 * self.l + 1)

    @property
    def label(self) -> str:
        return f"{self.n}{SUBSHELL_LETTER.get(self.l, '?')}"

    def quantum_states(self) -> list[QuantumState]:
        """Generate individual QuantumState objects following Hund's rules.

        Fill order: first all m_l values with spin +1/2 (maximise spin),
        then fill with spin -1/2 starting from most negative m_l.
        """
        states = []
        remaining = self.electron_count
        ml_values = list(range(-self.l, self.l + 1))

        # First pass: spin-up (+1/2) across all m_l
        for ml in ml_values:
            if remaining <= 0:
                break
            states.append(QuantumState(self.n, self.l, ml, +0.5))
            remaining -= 1

        # Second pass: spin-down (-1/2)
        for ml in ml_values:
            if remaining <= 0:
                break
            states.append(QuantumState(self.n, self.l, ml, -0.5))
            remaining -= 1

        return states

    def __repr__(self):
        return f"{self.label}^{self.electron_count}"


# ===================================================================
# Aufbau filling order
# ===================================================================

def aufbau_order(max_n: int = 7) -> list[tuple[int, int]]:
    """Generate subshell filling order using the (n+l, n) rule.

    Returns list of (n, l) tuples in Aufbau filling order:
    1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s, 4d, 5p, 6s, 4f, 5d, 6p, ...
    """
    subshells = []
    for n in range(1, max_n + 1):
        for l in range(n):
            subshells.append((n, l))
    # Sort by (n+l, n) -- Madelung rule
    subshells.sort(key=lambda nl: (nl[0] + nl[1], nl[0]))
    return subshells


# ===================================================================
# Known electron configuration exceptions
# ===================================================================

# Exceptions encoded as {Z: [(n, l, count), ...]} -- full ground state config.
# Only elements whose config deviates from Aufbau prediction are listed.
AUFBAU_EXCEPTIONS: dict[int, list[tuple[int, int, int]]] = {
    24: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,5),(4,0,1)],         # Cr
    29: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,1)],        # Cu
    41: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),
         (4,2,4),(5,0,1)],                                                  # Nb
    42: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),
         (4,2,5),(5,0,1)],                                                  # Mo
    44: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),
         (4,2,7),(5,0,1)],                                                  # Ru
    45: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),
         (4,2,8),(5,0,1)],                                                  # Rh
    46: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),
         (4,2,10)],                                                          # Pd
    47: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),
         (4,2,10),(5,0,1)],                                                  # Ag
    57: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),
         (4,2,10),(5,0,2),(5,1,6),(5,2,1),(6,0,2)],                         # La
    58: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),
         (4,2,10),(4,3,1),(5,0,2),(5,1,6),(5,2,1),(6,0,2)],                 # Ce
    64: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),
         (4,2,10),(4,3,7),(5,0,2),(5,1,6),(5,2,1),(6,0,2)],                 # Gd
    78: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),
         (4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,9),(6,0,1)],                # Pt
    79: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),
         (4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,1)],               # Au
}

# Noble gas cores for shorthand notation
NOBLE_GASES = {
    2:  "He",  10: "Ne",  18: "Ar",
    36: "Kr",  54: "Xe",  86: "Rn",
    118: "Og",
}


# ===================================================================
# Slater's rules for effective nuclear charge
# ===================================================================

def slater_zeff(Z: int, n: int, l: int,
                config: list[tuple[int, int, int]]) -> float:
    """Compute Slater's effective nuclear charge Z_eff for an electron
    in subshell (n, l) of an atom with atomic number Z.

    Slater grouping: (1s)(2s,2p)(3s,3p)(3d)(4s,4p)(4d)(4f)(5s,5p)...

    Shielding rules:
        Same group:     0.35 each  (0.30 for 1s)
        (n-1) group:    0.85 each  [for s,p electrons]
        All lower:      1.00 each  [for s,p electrons]
        For d,f:  same group 0.35, all lower groups 1.00
    """
    def slater_group(ni, li):
        """Assign Slater group index for sorting."""
        if li <= 1:
            return (ni, 0)  # s and p are in the same group
        else:
            return (ni, li)  # d, f each in their own group

    target_group = slater_group(n, l)
    is_sp = (l <= 1)  # s or p electron

    # Build list of all occupied groups with their electron counts
    groups = {}
    for (ni, li, count) in config:
        g = slater_group(ni, li)
        groups[g] = groups.get(g, 0) + count

    sigma = 0.0

    for g, count in groups.items():
        if g == target_group:
            # Electrons in the same group (exclude self)
            same_count = count - 1
            if same_count > 0:
                s_each = 0.30 if (n == 1 and l == 0) else 0.35
                sigma += same_count * s_each
        elif g < target_group:
            if is_sp:
                # For s/p electron: (n-1) group shields 0.85, lower shields 1.0
                g_n = g[0] if g[1] == 0 else g[0]  # principal quantum number of group
                target_n = n
                if g_n == target_n - 1 and g[1] <= 1:
                    sigma += count * 0.85
                elif g_n == target_n and g[1] > 1:
                    # Same n, but d/f group (below s/p in Slater ordering)
                    sigma += count * 0.85
                else:
                    sigma += count * 1.00
            else:
                # For d/f electron: all lower groups shield 1.0
                sigma += count * 1.00

    return Z - sigma


# ===================================================================
# The Quantum Atom
# ===================================================================

class QuantumAtom:
    """Represents an atom using the quantum mechanical model.

    Parameters
    ----------
    atomic_number : int
        Number of protons (Z). Determines the element.
    mass_number : int or None
        Protons + neutrons. Defaults to rounded standard atomic weight.
    charge : int
        Net ionic charge (0 = neutral).
    """

    def __init__(self, atomic_number: int, mass_number: int = None, charge: int = 0):
        if atomic_number < 1 or atomic_number > 118:
            raise ValueError(f"Atomic number must be 1-118, got {atomic_number}")

        self.atomic_number = atomic_number
        self.protons = atomic_number
        symbol, name, weight = ELEMENTS[atomic_number]
        self.symbol = symbol
        self.name = name
        self.standard_atomic_weight = weight

        self.mass_number = mass_number if mass_number is not None else round(weight)
        self.neutrons = self.mass_number - self.protons
        self.charge = charge
        self.num_electrons = self.protons - charge

        if self.num_electrons < 0:
            raise ValueError("Charge exceeds proton count")

        # Build ground-state electron configuration
        self.subshells = self._build_configuration()
        self.quantum_states = self._enumerate_quantum_states()

    # ------ configuration ------

    def _build_configuration(self) -> list[Subshell]:
        """Build the ground-state electron configuration."""
        # Check for known exceptions (neutral atoms only)
        if self.charge == 0 and self.atomic_number in AUFBAU_EXCEPTIONS:
            return [
                Subshell(n, l, count)
                for n, l, count in AUFBAU_EXCEPTIONS[self.atomic_number]
            ]

        # Standard Aufbau filling
        order = aufbau_order()
        remaining = self.num_electrons
        subshells = []

        for n, l in order:
            if remaining <= 0:
                break
            capacity = 2 * (2 * l + 1)
            count = min(remaining, capacity)
            subshells.append(Subshell(n, l, count))
            remaining -= count

        return subshells

    def _enumerate_quantum_states(self) -> list[QuantumState]:
        """List all individual quantum states following Hund's rules."""
        states = []
        for ss in self.subshells:
            states.extend(ss.quantum_states())
        return states

    # ------ configuration as tuples (for Slater's rules) ------

    @property
    def config_tuples(self) -> list[tuple[int, int, int]]:
        """Configuration as [(n, l, count), ...]."""
        return [(ss.n, ss.l, ss.electron_count) for ss in self.subshells]

    # ------ string representations ------

    def electron_configuration_string(self) -> str:
        """Full configuration string, e.g. '1s2 2s2 2p6 3s2 3p2'."""
        return " ".join(
            f"{ss.label}{ss.electron_count}" for ss in self.subshells
        )

    def noble_gas_notation(self) -> str:
        """Shorthand using noble gas core, e.g. '[Ne] 3s2 3p2'."""
        core_z = 0
        core_symbol = ""
        electrons_accounted = 0

        for z, sym in sorted(NOBLE_GASES.items()):
            if z < self.atomic_number and z <= self.num_electrons:
                core_z = z
                core_symbol = sym

        if core_z == 0:
            return self.electron_configuration_string()

        # Find how many subshells form the core
        electrons_accounted = 0
        core_end_idx = 0
        for i, ss in enumerate(self.subshells):
            electrons_accounted += ss.electron_count
            core_end_idx = i + 1
            if electrons_accounted == core_z:
                break

        valence_str = " ".join(
            f"{ss.label}{ss.electron_count}"
            for ss in self.subshells[core_end_idx:]
        )
        return f"[{core_symbol}] {valence_str}".strip()

    # ------ valence / core ------

    @property
    def valence_shell(self) -> int:
        """Highest principal quantum number with electrons."""
        if not self.subshells:
            return 0
        return max(ss.n for ss in self.subshells)

    @property
    def valence_electrons(self) -> int:
        """Number of electrons in the valence shell."""
        vs = self.valence_shell
        return sum(ss.electron_count for ss in self.subshells if ss.n == vs)

    @property
    def core_electrons(self) -> int:
        """Number of core (non-valence) electrons."""
        return self.num_electrons - self.valence_electrons

    # ------ effective nuclear charge ------

    def effective_nuclear_charge(self, n: int, l: int) -> float:
        """Slater's Z_eff for an electron in subshell (n, l)."""
        return slater_zeff(self.atomic_number, n, l, self.config_tuples)

    # ------ energy approximations ------

    def orbital_energy_eV(self, n: int, l: int) -> float:
        """Approximate orbital energy using Slater's Z_eff.

        E_{n,l} ~ -13.6 eV * Z_eff^2 / n^2
        """
        z_eff = self.effective_nuclear_charge(n, l)
        return -RYDBERG_ENERGY_EV * z_eff**2 / n**2

    def ionization_energy_eV(self) -> float:
        """Approximate first ionization energy (Slater model), in eV.

        Energy required to remove the outermost electron.
        """
        if not self.subshells:
            return 0.0
        outermost = self.subshells[-1]
        return -self.orbital_energy_eV(outermost.n, outermost.l)

    # ------ angular momentum ------

    def total_spin_quantum_number(self) -> float:
        """Total spin S = sum of m_s values (maximised by Hund's rule)."""
        return abs(sum(qs.ms for qs in self.quantum_states))

    def spin_multiplicity(self) -> int:
        """2S + 1 -- spin multiplicity of the ground state."""
        S = self.total_spin_quantum_number()
        return int(2 * S + 1)

    # ------ display ------

    def __repr__(self):
        charge_str = ""
        if self.charge > 0:
            charge_str = f"(+{self.charge})"
        elif self.charge < 0:
            charge_str = f"({self.charge})"
        return (f"QuantumAtom({self.symbol}, Z={self.atomic_number}, "
                f"e={self.num_electrons}{charge_str})")

    def summary(self) -> str:
        """Print a detailed summary of the atom's quantum description."""
        lines = [
            f"{'='*60}",
            f"  {self.name} ({self.symbol})",
            f"  Atomic number:  {self.atomic_number}",
            f"  Mass number:    {self.mass_number}",
            f"  Protons: {self.protons}  Neutrons: {self.neutrons}  "
            f"Electrons: {self.num_electrons}",
            f"  Charge: {self.charge:+d}",
            f"{'='*60}",
            f"  Electron configuration:",
            f"    {self.electron_configuration_string()}",
            f"    {self.noble_gas_notation()}",
            f"  Valence electrons: {self.valence_electrons}",
            f"  Spin multiplicity: {self.spin_multiplicity()}",
            f"{'='*60}",
            f"  Subshell         e-   Z_eff    Energy (eV)",
            f"  {'-'*48}",
        ]
        for ss in self.subshells:
            z_eff = self.effective_nuclear_charge(ss.n, ss.l)
            e_eV = self.orbital_energy_eV(ss.n, ss.l)
            lines.append(
                f"  {ss.label:<6}  {ss.electron_count:>4}   "
                f"{z_eff:>7.3f}    {e_eV:>10.4f}"
            )
        lines.append(f"  {'-'*48}")
        lines.append(
            f"  Ionization energy (Slater): "
            f"{self.ionization_energy_eV():.4f} eV"
        )
        lines.append("")

        # Show all quantum states for small atoms
        if self.num_electrons <= 18:
            lines.append(f"  All quantum states (n, l, ml, ms):")
            for qs in self.quantum_states:
                lines.append(f"    {qs}")

        return "\n".join(lines)


# ===================================================================
# Convenience constructors
# ===================================================================

def from_symbol(symbol: str, mass_number: int = None, charge: int = 0):
    """Create a QuantumAtom from an element symbol."""
    sym = symbol.strip()
    if len(sym) > 1:
        sym = sym[0].upper() + sym[1:].lower()
    else:
        sym = sym.upper()
    if sym not in SYMBOL_TO_Z:
        raise ValueError(f"Unknown element symbol: {symbol}")
    return QuantumAtom(SYMBOL_TO_Z[sym], mass_number, charge)


def from_name(name: str, mass_number: int = None, charge: int = 0):
    """Create a QuantumAtom from an element name."""
    name_lower = name.strip().lower()
    for z, (sym, elem_name, _) in ELEMENTS.items():
        if elem_name.lower() == name_lower:
            return QuantumAtom(z, mass_number, charge)
    raise ValueError(f"Unknown element name: {name}")


# ===================================================================
# Main -- demo
# ===================================================================

if __name__ == "__main__":
    print("Quantum Mechanical Model of the Atom\n")

    hydrogen = QuantumAtom(1)
    carbon = QuantumAtom(6)
    iron = from_symbol("Fe")
    chromium = from_symbol("Cr")  # Aufbau exception
    copper = from_symbol("Cu")    # Aufbau exception

    for atom in [hydrogen, carbon, iron, chromium, copper]:
        print(atom.summary())
        print()

    # Show Aufbau filling order
    print("Aufbau filling order:")
    order = aufbau_order()
    labels = [f"{n}{SUBSHELL_LETTER[l]}" for n, l in order]
    print("  " + " -> ".join(labels))
