"""
Covalent Bond Modelling

Models the different types of covalent bonds and their physical properties:

    - Single, double, and triple bonds (sigma/pi decomposition)
    - Polar vs nonpolar classification
    - Coordinate (dative) bonds
    - Bond length, dissociation energy, dipole moment estimation
    - Molecular polarity analysis from bond vectors

Integrates with the existing Lewis structure and element data modules.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum, auto

from molbuilder.core.elements import SYMBOL_TO_Z
from molbuilder.core.element_properties import (
    electronegativity,
    estimated_bond_length_pm,
    covalent_radius_pm,
    PAULING_ELECTRONEGATIVITY,
)
from molbuilder.core.bond_data import BDE_TABLE, NONPOLAR_THRESHOLD, POLAR_COVALENT_MAX
from molbuilder.core.constants import DEBYE_PER_E_ANGSTROM


# ===================================================================
# Enumerations
# ===================================================================

class BondPolarity(Enum):
    """Classification of a bond by electronegativity difference."""
    NONPOLAR_COVALENT = auto()
    POLAR_COVALENT = auto()
    IONIC = auto()


class OrbitalType(Enum):
    """Orbital overlap type contributing to a covalent bond."""
    SIGMA = auto()   # head-on overlap along the bond axis
    PI = auto()       # lateral overlap above/below the bond axis


# ===================================================================
# Sigma / Pi orbital composition
# ===================================================================

@dataclass(frozen=True)
class OrbitalContribution:
    """Describes one orbital overlap contributing to a bond."""
    orbital_type: OrbitalType
    description: str

    def __repr__(self):
        label = "sigma" if self.orbital_type is OrbitalType.SIGMA else "pi"
        return f"{label}({self.description})"


def sigma_pi_composition(bond_order: int) -> list[OrbitalContribution]:
    """Return the sigma/pi orbital makeup for a given bond order.

    Single bond:  1 sigma
    Double bond:  1 sigma + 1 pi
    Triple bond:  1 sigma + 2 pi
    """
    if bond_order < 1 or bond_order > 3:
        raise ValueError(f"Bond order must be 1, 2, or 3; got {bond_order}")

    orbitals = [OrbitalContribution(OrbitalType.SIGMA,
                                    "head-on overlap along bond axis")]
    if bond_order >= 2:
        orbitals.append(OrbitalContribution(OrbitalType.PI,
                                            "lateral overlap, perpendicular plane 1"))
    if bond_order == 3:
        orbitals.append(OrbitalContribution(OrbitalType.PI,
                                            "lateral overlap, perpendicular plane 2"))
    return orbitals


# ===================================================================
# BDE key helper
# ===================================================================

def _bde_key(sym_a: str, sym_b: str, order: int) -> tuple[str, str, int]:
    """Normalise to alphabetical order for BDE table lookup."""
    a, b = sorted([sym_a, sym_b])
    return (a, b, order)


# ===================================================================
# CovalentBond
# ===================================================================

@dataclass
class CovalentBond:
    """Model of a covalent bond between two atoms.

    Parameters
    ----------
    symbol_a : str
        Element symbol of first atom.
    symbol_b : str
        Element symbol of second atom.
    bond_order : int
        1 = single, 2 = double, 3 = triple.
    is_coordinate : bool
        True if both bonding electrons were donated by one atom
        (coordinate / dative bond).
    donor : str or None
        For coordinate bonds, the symbol of the electron-pair donor.
    """
    symbol_a: str
    symbol_b: str
    bond_order: int = 1
    is_coordinate: bool = False
    donor: str | None = None

    def __post_init__(self):
        if self.symbol_a not in SYMBOL_TO_Z:
            raise ValueError(f"Unknown element: {self.symbol_a}")
        if self.symbol_b not in SYMBOL_TO_Z:
            raise ValueError(f"Unknown element: {self.symbol_b}")
        if self.bond_order not in (1, 2, 3):
            raise ValueError(f"Bond order must be 1, 2, or 3; got {self.bond_order}")
        if self.is_coordinate and self.donor is None:
            raise ValueError("Coordinate bond requires a donor symbol")

    # --- Electronegativity & polarity ---

    @property
    def en_a(self) -> float:
        return electronegativity(self.symbol_a)

    @property
    def en_b(self) -> float:
        return electronegativity(self.symbol_b)

    @property
    def delta_en(self) -> float:
        """Absolute electronegativity difference (Pauling scale)."""
        return abs(self.en_a - self.en_b)

    @property
    def polarity(self) -> BondPolarity:
        """Classify the bond by its electronegativity difference."""
        d = self.delta_en
        if d < NONPOLAR_THRESHOLD:
            return BondPolarity.NONPOLAR_COVALENT
        elif d < POLAR_COVALENT_MAX:
            return BondPolarity.POLAR_COVALENT
        else:
            return BondPolarity.IONIC

    @property
    def partial_positive(self) -> str:
        """Symbol of the atom bearing partial positive charge (delta+)."""
        if self.en_a <= self.en_b:
            return self.symbol_a
        return self.symbol_b

    @property
    def partial_negative(self) -> str:
        """Symbol of the atom bearing partial negative charge (delta-)."""
        if self.en_a > self.en_b:
            return self.symbol_a
        return self.symbol_b

    @property
    def percent_ionic_character(self) -> float:
        """Estimate ionic character using Pauling's equation.

        %ionic = 100 * (1 - exp(-0.25 * delta_EN^2))
        """
        return 100.0 * (1.0 - math.exp(-0.25 * self.delta_en ** 2))

    # --- Bond length ---

    @property
    def bond_length_pm(self) -> float:
        """Estimated bond length in picometres."""
        return estimated_bond_length_pm(self.symbol_a, self.symbol_b,
                                        self.bond_order)

    @property
    def bond_length_angstrom(self) -> float:
        """Estimated bond length in angstroms."""
        return self.bond_length_pm / 100.0

    # --- Bond dissociation energy ---

    @property
    def dissociation_energy_kj(self) -> float | None:
        """Mean bond dissociation energy in kJ/mol, or None if unknown."""
        return BDE_TABLE.get(_bde_key(self.symbol_a, self.symbol_b, self.bond_order))

    @property
    def dissociation_energy_ev(self) -> float | None:
        """Mean bond dissociation energy in eV per bond, or None if unknown."""
        e = self.dissociation_energy_kj
        if e is None:
            return None
        return e / 96.485  # 1 eV = 96.485 kJ/mol

    # --- Orbital composition ---

    @property
    def orbital_contributions(self) -> list[OrbitalContribution]:
        """Sigma and pi orbital contributions making up this bond."""
        return sigma_pi_composition(self.bond_order)

    @property
    def sigma_bonds(self) -> int:
        """Number of sigma bonds (always 1 for a covalent bond)."""
        return 1

    @property
    def pi_bonds(self) -> int:
        """Number of pi bonds (0, 1, or 2)."""
        return self.bond_order - 1

    # --- Dipole moment estimate ---

    @property
    def dipole_moment_debye(self) -> float:
        """Estimate bond dipole moment in Debye using the point-charge model.

        mu = q * d * 4.8032

        where q is the partial charge in electron units (percent_ionic / 100)
        and d is the bond length in Angstroms.  The constant 4.8032 converts
        electron-Angstroms to Debye (1 D = 3.336e-30 C*m).

        This is a simplified model.  Real dipole moments depend on the full
        3D electron distribution.
        """
        if self.delta_en == 0.0:
            return 0.0
        partial_charge_e = self.percent_ionic_character / 100.0
        return partial_charge_e * self.bond_length_angstrom * DEBYE_PER_E_ANGSTROM

    # --- Display ---

    @property
    def order_symbol(self) -> str:
        return {1: "-", 2: "=", 3: "#"}.get(self.bond_order, "?")

    @property
    def order_label(self) -> str:
        return {1: "single", 2: "double", 3: "triple"}.get(self.bond_order, "?")

    def __repr__(self):
        tag = " (coordinate)" if self.is_coordinate else ""
        return f"CovalentBond({self.symbol_a}{self.order_symbol}{self.symbol_b}, {self.order_label}{tag})"

    def summary(self) -> str:
        """Return a multi-line human-readable summary of this bond."""
        lines = [
            f"  {self.symbol_a}{self.order_symbol}{self.symbol_b}  "
            f"({self.order_label} bond)",
        ]

        # Coordinate info
        if self.is_coordinate:
            lines.append(f"    Type          : coordinate (dative) bond")
            lines.append(f"    Donor         : {self.donor}")
        else:
            lines.append(f"    Type          : standard covalent bond")

        # Polarity
        pol = self.polarity
        pol_label = {
            BondPolarity.NONPOLAR_COVALENT: "nonpolar covalent",
            BondPolarity.POLAR_COVALENT:    "polar covalent",
            BondPolarity.IONIC:             "ionic",
        }[pol]
        lines.append(f"    Polarity      : {pol_label}  "
                      f"(delta EN = {self.delta_en:.2f})")
        if pol == BondPolarity.POLAR_COVALENT:
            lines.append(f"    Dipole        : "
                          f"delta+ on {self.partial_positive}, "
                          f"delta- on {self.partial_negative}")
            lines.append(f"    % ionic char  : {self.percent_ionic_character:.1f}%")
            lines.append(f"    Dipole moment : ~{self.dipole_moment_debye:.2f} D")

        # Geometry
        lines.append(f"    Bond length   : {self.bond_length_pm:.0f} pm  "
                      f"({self.bond_length_angstrom:.2f} A)")

        # Energy
        bde = self.dissociation_energy_kj
        if bde is not None:
            lines.append(f"    BDE           : {bde:.0f} kJ/mol  "
                          f"({self.dissociation_energy_ev:.2f} eV)")
        else:
            lines.append(f"    BDE           : no reference data")

        # Orbital decomposition
        orb_strs = [repr(o) for o in self.orbital_contributions]
        lines.append(f"    Orbitals      : {' + '.join(orb_strs)}")

        return "\n".join(lines)


# ===================================================================
# Molecular polarity analysis
# ===================================================================

@dataclass
class MolecularBondAnalysis:
    """Analyse all covalent bonds in a molecule from its Lewis structure.

    Constructs CovalentBond objects for each bond and determines whether
    the molecule as a whole is likely polar or nonpolar based on geometry
    symmetry and individual bond polarities.

    Parameters
    ----------
    formula : str
        Molecular formula (e.g. 'H2O', 'CO2', 'CH4').
    charge : int
        Net charge on the species.
    """
    formula: str
    charge: int = 0
    bonds: list[CovalentBond] = field(default_factory=list, init=False)
    _lewis: object = field(default=None, init=False, repr=False)

    def __post_init__(self):
        from molbuilder.bonding.lewis import LewisStructure
        self._lewis = LewisStructure(self.formula, self.charge)
        self._build_bonds()

    def _build_bonds(self):
        """Create CovalentBond objects from the Lewis structure."""
        lew = self._lewis
        for bond in lew.bonds:
            sym_a = lew.atoms[bond.atom_a]
            sym_b = lew.atoms[bond.atom_b]
            self.bonds.append(CovalentBond(sym_a, sym_b, bond.order))

    @property
    def total_sigma_bonds(self) -> int:
        return sum(b.sigma_bonds for b in self.bonds)

    @property
    def total_pi_bonds(self) -> int:
        return sum(b.pi_bonds for b in self.bonds)

    @property
    def all_bonds_nonpolar(self) -> bool:
        return all(b.polarity == BondPolarity.NONPOLAR_COVALENT
                   for b in self.bonds)

    @property
    def has_lone_pairs_on_central(self) -> bool:
        return self._lewis.lone_pairs_on_central() > 0

    @property
    def is_symmetric(self) -> bool:
        """Heuristic: molecule is symmetric if all terminal atoms are
        identical and the central atom has no lone pairs.

        Diatomic molecules (A-B) are symmetric only if both atoms are
        the same element (homonuclear).
        """
        lew = self._lewis
        if len(lew.atoms) <= 2:
            # Homonuclear diatomic (H2, O2, N2) -> symmetric
            # Heteronuclear diatomic (HCl, CO) -> not symmetric
            return len(set(lew.atoms)) == 1
        terminal_syms = {lew.atoms[i] for i in lew.terminal_indices}
        orders = {b.bond_order for b in self.bonds}
        return (len(terminal_syms) == 1
                and len(orders) == 1
                and not self.has_lone_pairs_on_central)

    @property
    def molecular_polarity(self) -> str:
        """Predict whether the molecule is polar or nonpolar.

        A molecule is nonpolar if:
            - all bonds are nonpolar, OR
            - it has a symmetric geometry that cancels dipoles
        Otherwise it is polar.
        """
        if self.all_bonds_nonpolar:
            return "nonpolar"
        if self.is_symmetric:
            return "nonpolar (symmetric -- dipoles cancel)"
        return "polar"

    def summary(self) -> str:
        lines = [
            f"{'=' * 60}",
            f"  Covalent Bond Analysis: {self.formula}"
            + (f" (charge {self.charge:+d})" if self.charge else ""),
            f"{'=' * 60}",
        ]
        for b in self.bonds:
            lines.append(b.summary())
            lines.append("")

        lines.append(f"  Total sigma bonds : {self.total_sigma_bonds}")
        lines.append(f"  Total pi bonds    : {self.total_pi_bonds}")
        lines.append(f"  Molecular polarity: {self.molecular_polarity}")
        lines.append(f"{'=' * 60}")
        return "\n".join(lines)


# ===================================================================
# Convenience constructors
# ===================================================================

def single_bond(sym_a: str, sym_b: str) -> CovalentBond:
    """Create a single covalent bond."""
    return CovalentBond(sym_a, sym_b, bond_order=1)


def double_bond(sym_a: str, sym_b: str) -> CovalentBond:
    """Create a double covalent bond."""
    return CovalentBond(sym_a, sym_b, bond_order=2)


def triple_bond(sym_a: str, sym_b: str) -> CovalentBond:
    """Create a triple covalent bond."""
    return CovalentBond(sym_a, sym_b, bond_order=3)


def coordinate_bond(donor: str, acceptor: str,
                     bond_order: int = 1) -> CovalentBond:
    """Create a coordinate (dative) covalent bond.

    In a coordinate bond both shared electrons originate from the donor.
    Example: NH3 -> BF3  (N donates a lone pair to B).
    """
    return CovalentBond(donor, acceptor, bond_order=bond_order,
                        is_coordinate=True, donor=donor)
