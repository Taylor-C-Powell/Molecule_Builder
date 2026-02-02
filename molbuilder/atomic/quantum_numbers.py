"""
Quantum Numbers and Electron Configuration

Defines quantum number containers, subshell representations,
Aufbau filling order, and known configuration exceptions.
"""

import math
from dataclasses import dataclass, field
from molbuilder.core.constants import HBAR
from molbuilder.core.elements import NOBLE_GASES


# ---------------------------------------------------------------------------
# Orbital labels
# ---------------------------------------------------------------------------
SUBSHELL_LETTER = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g", 5: "h", 6: "i"}
SUBSHELL_NUMBER = {v: k for k, v in SUBSHELL_LETTER.items()}

ORBITAL_NAMES = {
    (0, 0): "s",
    (1, -1): "p_y", (1, 0): "p_z", (1, 1): "p_x",
    (2, -2): "d_xy", (2, -1): "d_yz", (2, 0): "d_z2",
    (2, 1): "d_xz", (2, 2): "d_x2-y2",
    (3, -3): "f_y(3x2-y2)", (3, -2): "f_xyz", (3, -1): "f_yz2",
    (3, 0): "f_z3", (3, 1): "f_xz2", (3, 2): "f_z(x2-y2)",
    (3, 3): "f_x(x2-3y2)",
}


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
    89: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),
         (4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2),(6,1,6),
         (6,2,1),(7,0,2)],                                                   # Ac
    90: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),
         (4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(6,0,2),(6,1,6),
         (6,2,2),(7,0,2)],                                                   # Th
    91: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),
         (4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(5,3,2),(6,0,2),(6,1,6),
         (6,2,1),(7,0,2)],                                                   # Pa
    92: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),
         (4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(5,3,3),(6,0,2),(6,1,6),
         (6,2,1),(7,0,2)],                                                   # U
    93: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),
         (4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(5,3,4),(6,0,2),(6,1,6),
         (6,2,1),(7,0,2)],                                                   # Np
    96: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),
         (4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(5,3,7),(6,0,2),(6,1,6),
         (6,2,1),(7,0,2)],                                                   # Cm
    103: [(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,10),(4,0,2),(4,1,6),
          (4,2,10),(4,3,14),(5,0,2),(5,1,6),(5,2,10),(5,3,14),(6,0,2),(6,1,6),
          (7,0,2),(7,1,1)],                                                   # Lr
}
