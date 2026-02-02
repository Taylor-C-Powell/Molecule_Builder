"""
Physical constants (SI units) used throughout molbuilder.

Consolidated from bohr_model.py, quantum_wavefunctions.py, and
atomic_particle_physical_laws.py.
"""

import math

# ---------------------------------------------------------------------------
# Fundamental constants
# ---------------------------------------------------------------------------
BOHR_RADIUS_M = 5.29177210903e-11          # metres
BOHR_RADIUS_PM = 52.9177210903             # picometres
PLANCK_CONSTANT = 6.62607015e-34           # J s
HBAR = PLANCK_CONSTANT / (2 * math.pi)    # reduced Planck constant (J s)
SPEED_OF_LIGHT = 2.99792458e8             # m/s
ELECTRON_CHARGE = 1.602176634e-19         # C
ELECTRON_MASS = 9.1093837015e-31          # kg
COULOMB_CONSTANT = 8.9875517873681764e9   # N m^2 / C^2
VACUUM_PERMITTIVITY = 8.8541878128e-12    # F/m

# ---------------------------------------------------------------------------
# Unit conversions
# ---------------------------------------------------------------------------
EV_TO_JOULES = 1.602176634e-19
RYDBERG_ENERGY_EV = 13.605693122994       # eV
DEBYE_PER_E_ANGSTROM = 4.8032            # 1 e*A = 4.8032 D

# ---------------------------------------------------------------------------
# Fine structure constant
# ---------------------------------------------------------------------------
FINE_STRUCTURE_ALPHA = (
    ELECTRON_CHARGE**2
    / (4 * math.pi * VACUUM_PERMITTIVITY * HBAR * SPEED_OF_LIGHT)
)

# ---------------------------------------------------------------------------
# Shell capacity
# ---------------------------------------------------------------------------
MAX_ELECTRONS_PER_SHELL = [2 * n**2 for n in range(1, 8)]  # shells 1-7


# ---------------------------------------------------------------------------
# Electrostatics
# ---------------------------------------------------------------------------
def coulombs_law(q1: float, q2: float, r: float) -> float:
    """Coulomb force between two charges q1, q2 (Coulombs) separated by r (metres)."""
    if r <= 0:
        raise ValueError("Distance must be positive")
    return COULOMB_CONSTANT * (abs(q1) * abs(q2)) / r**2
