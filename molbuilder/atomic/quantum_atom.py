"""
Quantum Mechanical Model of the Atom

Represents atoms using the quantum mechanical framework:
    - Four quantum numbers (n, l, m_l, m_s) per electron
    - Electron configuration via Aufbau principle, Hund's rule, Pauli exclusion
    - Slater's rules for effective nuclear charge (Z_eff)
    - Approximate orbital energies for multi-electron atoms
    - Known configuration exceptions for d- and f-block elements
"""

from molbuilder.core.constants import RYDBERG_ENERGY_EV
from molbuilder.core.elements import ELEMENTS, SYMBOL_TO_Z, NOBLE_GASES
from molbuilder.atomic.quantum_numbers import (
    QuantumState,
    Subshell,
    aufbau_order,
    AUFBAU_EXCEPTIONS,
)


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
                # For s/p electron: (n-1) shell shields 0.85, lower shields 1.0
                # The (n-1) shell includes ALL groups whose principal quantum
                # number equals n-1, regardless of angular momentum (s/p/d/f).
                g_n = g[0]  # principal quantum number of the inner group
                if g_n == n - 1:
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

def from_symbol(symbol: str, mass_number: int = None, charge: int = 0) -> QuantumAtom:
    """Create a QuantumAtom from an element symbol."""
    sym = symbol.strip()
    if len(sym) > 1:
        sym = sym[0].upper() + sym[1:].lower()
    else:
        sym = sym.upper()
    if sym not in SYMBOL_TO_Z:
        raise ValueError(f"Unknown element symbol: {symbol}")
    return QuantumAtom(SYMBOL_TO_Z[sym], mass_number, charge)


def from_name(name: str, mass_number: int = None, charge: int = 0) -> QuantumAtom:
    """Create a QuantumAtom from an element name."""
    name_lower = name.strip().lower()
    for z, (sym, elem_name, _) in ELEMENTS.items():
        if elem_name.lower() == name_lower:
            return QuantumAtom(z, mass_number, charge)
    raise ValueError(f"Unknown element name: {name}")
