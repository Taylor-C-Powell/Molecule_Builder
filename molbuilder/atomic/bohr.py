"""
Bohr Model of the Atom

Represents atoms using the Bohr model with quantized energy levels,
electron shell configurations, and orbital mechanics. Provides both
computational physics and animated visualization.

Physics reference:
    - Orbital radius:  r_n = n^2 * a_0 / Z
    - Energy level:    E_n = -13.6 eV * Z^2 / n^2
    - Photon energy:   delta_E = E_final - E_initial
    - Wavelength:      lambda = h*c / |delta_E|
    - Orbital velocity: v_n = (Z * e^2) / (n * hbar)
"""

import math
import warnings
from molbuilder.core.constants import (
    BOHR_RADIUS_M,
    PLANCK_CONSTANT,
    HBAR,
    SPEED_OF_LIGHT,
    ELECTRON_CHARGE,
    COULOMB_CONSTANT,
    EV_TO_JOULES,
    RYDBERG_ENERGY_EV,
    VACUUM_PERMITTIVITY,
)
from molbuilder.core.elements import ELEMENTS, SYMBOL_TO_Z


# ---------------------------------------------------------------------------
# Core Bohr model class
# ---------------------------------------------------------------------------
class BohrAtom:
    """Represents an atom using the Bohr model.

    The Bohr model is only physically valid for hydrogen-like (one-electron)
    atoms and ions (H, He+, Li2+, etc.). For multi-electron atoms, results
    are qualitative approximations. Use QuantumAtom from
    molbuilder.atomic.quantum_atom for multi-electron systems.

    Parameters
    ----------
    atomic_number : int
        Number of protons (Z). Determines the element.
    mass_number : int or None
        Protons + neutrons. Defaults to rounded standard atomic weight.
    charge : int
        Net ionic charge (0 = neutral, +1 = lost one electron, etc.).
    """

    def __init__(self, atomic_number: int, mass_number: int = None, charge: int = 0):
        if atomic_number < 1 or atomic_number > 118:
            raise ValueError(f"Atomic number must be 1-118, got {atomic_number}")

        self.atomic_number = atomic_number  # Z
        self.protons = atomic_number
        symbol, name, weight = ELEMENTS[atomic_number]
        self.symbol = symbol
        self.name = name
        self.standard_atomic_weight = weight

        if mass_number is None:
            self.mass_number = round(weight)
        else:
            self.mass_number = mass_number

        self.neutrons = self.mass_number - self.protons
        self.charge = charge
        self.num_electrons = self.protons - charge

        if self.num_electrons < 0:
            raise ValueError("Charge exceeds number of protons (no electrons left)")

        if self.num_electrons > 1:
            warnings.warn(
                f"BohrAtom(Z={atomic_number}): Bohr model is only exact for hydrogen-like "
                f"(one-electron) systems. Results for Z={atomic_number} are qualitative. "
                f"Consider using QuantumAtom for multi-electron atoms.",
                stacklevel=2,
            )

        self.shell_config = self._compute_shell_config()
        self.num_shells = len(self.shell_config)

    # ----- electron configuration -----

    def _compute_shell_config(self) -> list[int]:
        """Fill electron shells following the 2n^2 rule."""
        remaining = self.num_electrons
        shells = []
        n = 1
        while remaining > 0 and n <= 7:
            capacity = 2 * n**2
            electrons_in_shell = min(remaining, capacity)
            shells.append(electrons_in_shell)
            remaining -= electrons_in_shell
            n += 1
        return shells

    # ----- Bohr radius & energy -----

    def orbital_radius(self, n: int) -> float:
        """Radius of the n-th Bohr orbit in metres.

        r_n = n^2 * a_0 / Z
        """
        if n < 1:
            raise ValueError("Principal quantum number n must be >= 1")
        return (n**2 * BOHR_RADIUS_M) / self.atomic_number

    def orbital_radius_pm(self, n: int) -> float:
        """Radius of the n-th Bohr orbit in picometres."""
        return self.orbital_radius(n) * 1e12

    def energy_level(self, n: int) -> float:
        """Energy of the n-th level in eV.

        E_n = -13.6 eV * Z^2 / n^2
        """
        if n < 1:
            raise ValueError("Principal quantum number n must be >= 1")
        return -RYDBERG_ENERGY_EV * (self.atomic_number**2) / (n**2)

    def energy_level_joules(self, n: int) -> float:
        """Energy of the n-th level in joules."""
        return self.energy_level(n) * EV_TO_JOULES

    def orbital_velocity(self, n: int) -> float:
        """Velocity of the electron in the n-th orbit in m/s.

        v_n = Z * e^2 / (n * 4*pi*eps_0 * hbar)
            = Z * alpha * c / n
        where alpha ~ 1/137 is the fine-structure constant.
        """
        alpha = ELECTRON_CHARGE**2 / (4 * math.pi * VACUUM_PERMITTIVITY * HBAR * SPEED_OF_LIGHT)
        return self.atomic_number * alpha * SPEED_OF_LIGHT / n

    def orbital_period(self, n: int) -> float:
        """Orbital period of the electron in shell n, in seconds."""
        r = self.orbital_radius(n)
        v = self.orbital_velocity(n)
        return 2 * math.pi * r / v

    # ----- transitions -----

    def transition_energy(self, n_initial: int, n_final: int) -> float:
        """Energy of a photon emitted (positive) or absorbed (negative)
        during a transition, in eV."""
        return self.energy_level(n_initial) - self.energy_level(n_final)

    def transition_wavelength(self, n_initial: int, n_final: int) -> float:
        """Wavelength in metres of a photon from a transition between levels."""
        delta_e = abs(self.transition_energy(n_initial, n_final)) * EV_TO_JOULES
        if delta_e == 0:
            return float('inf')
        return PLANCK_CONSTANT * SPEED_OF_LIGHT / delta_e

    def transition_wavelength_nm(self, n_initial: int, n_final: int) -> float:
        """Wavelength in nanometres."""
        return self.transition_wavelength(n_initial, n_final) * 1e9

    def ionization_energy(self) -> float:
        """Minimum energy (eV) to remove the outermost electron (Bohr approx)."""
        if self.num_shells == 0:
            return 0.0
        return -self.energy_level(self.num_shells)

    # ----- electrostatic force on electron -----

    def coulomb_force_on_electron(self, n: int) -> float:
        """Attractive Coulomb force on an electron in shell n, in newtons."""
        r = self.orbital_radius(n)
        return COULOMB_CONSTANT * (self.protons * ELECTRON_CHARGE**2) / r**2

    # ----- representation -----

    def __repr__(self):
        iso = f"-{self.mass_number}" if self.mass_number != round(self.standard_atomic_weight) else ""
        charge_str = ""
        if self.charge > 0:
            charge_str = f" (+{self.charge})"
        elif self.charge < 0:
            charge_str = f" ({self.charge})"
        return (f"BohrAtom({self.symbol}{iso}, Z={self.atomic_number}, "
                f"e={self.num_electrons}{charge_str}, "
                f"shells={self.shell_config})")

    def summary(self) -> str:
        """Return a human-readable summary of the atom."""
        lines = [
            f"{'='*50}",
            f"  {self.name} ({self.symbol})",
            f"  Atomic number (Z): {self.atomic_number}",
            f"  Mass number (A):   {self.mass_number}",
            f"  Protons:  {self.protons}   Neutrons: {self.neutrons}",
            f"  Electrons: {self.num_electrons}  Charge: {self.charge:+d}",
            f"  Shell configuration: {self.shell_config}",
            f"{'='*50}",
            f"  Shell   Electrons   Radius (pm)     Energy (eV)",
            f"  {'-'*46}",
        ]
        for i, count in enumerate(self.shell_config, start=1):
            r = self.orbital_radius_pm(i)
            e = self.energy_level(i)
            lines.append(f"  n={i:<4}  {count:<10}  {r:>10.2f}      {e:>10.4f}")
        lines.append(f"  {'-'*46}")
        lines.append(f"  Ionization energy (Bohr): {self.ionization_energy():.4f} eV")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Convenience constructors
# ---------------------------------------------------------------------------
def from_symbol(symbol: str, mass_number: int = None, charge: int = 0) -> BohrAtom:
    """Create a BohrAtom from an element symbol string."""
    symbol = symbol.capitalize()
    if len(symbol) > 1:
        symbol = symbol[0].upper() + symbol[1:].lower()
    if symbol not in SYMBOL_TO_Z:
        raise ValueError(f"Unknown element symbol: {symbol}")
    return BohrAtom(SYMBOL_TO_Z[symbol], mass_number, charge)


def from_name(name: str, mass_number: int = None, charge: int = 0) -> BohrAtom:
    """Create a BohrAtom from an element name string."""
    name_lower = name.strip().lower()
    for z, (sym, elem_name, _) in ELEMENTS.items():
        if elem_name.lower() == name_lower:
            return BohrAtom(z, mass_number, charge)
    raise ValueError(f"Unknown element name: {name}")
