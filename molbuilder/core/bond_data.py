"""
Bond reference data: standard bond lengths, dissociation energies,
torsion barrier parameters.

Consolidated from molecular_conformations.py (STANDARD_BOND_LENGTHS,
TORSION_BARRIERS) and covalent_bonds.py (BDE_TABLE).
"""

from molbuilder.core.element_properties import estimated_bond_length_angstrom

# ===================================================================
# Standard bond lengths in Angstroms (experimental averages)
# ===================================================================
# Keys: (symbol_a, symbol_b, order) with symbols in alphabetical order.

STANDARD_BOND_LENGTHS: dict[tuple[str, str, int], float] = {
    ("C", "C", 1): 1.54,
    ("C", "C", 2): 1.34,
    ("C", "C", 3): 1.20,
    ("C", "H", 1): 1.09,
    ("C", "O", 1): 1.43,
    ("C", "O", 2): 1.23,
    ("C", "N", 1): 1.47,
    ("C", "N", 2): 1.29,
    ("C", "N", 3): 1.16,
    ("C", "F", 1): 1.35,
    ("C", "Cl", 1): 1.77,
    ("Br", "C", 1): 1.94,
    ("C", "I", 1): 2.14,
    ("C", "S", 1): 1.82,
    ("H", "O", 1): 0.96,
    ("H", "N", 1): 1.01,
    ("H", "S", 1): 1.34,
    ("F", "F", 1): 1.42,
    ("Cl", "Cl", 1): 1.99,
    ("Br", "Br", 1): 2.28,
    ("H", "H", 1): 0.74,
    ("S", "S", 1): 2.05,
    ("N", "O", 1): 1.36,
    ("N", "O", 2): 1.21,
    ("O", "P", 1): 1.63,
    ("O", "P", 2): 1.48,
    ("H", "P", 1): 1.44,
}

# Standard bond angles in degrees
SP3_ANGLE = 109.47   # tetrahedral
SP2_ANGLE = 120.0    # trigonal planar
SP_ANGLE = 180.0     # linear


# ===================================================================
# Bond dissociation energies (kJ/mol)
# ===================================================================
# Sources: CRC Handbook, Darwent (NSRDS-NBS 31), Kerr & Stocker.
# Keys: (symbol_a, symbol_b, bond_order) alphabetical order.

BDE_TABLE: dict[tuple[str, str, int], float] = {
    # Single bonds
    ("H", "H", 1):  436,
    ("C", "H", 1):  413,
    ("C", "C", 1):  348,
    ("C", "N", 1):  293,
    ("C", "O", 1):  358,
    ("C", "F", 1):  485,
    ("C", "Cl", 1): 328,
    ("Br", "C", 1): 276,
    ("C", "I", 1):  240,
    ("C", "S", 1):  272,
    ("H", "N", 1):  391,
    ("N", "N", 1):  163,
    ("N", "O", 1):  201,
    ("H", "O", 1):  463,
    ("O", "O", 1):  146,
    ("F", "F", 1):  155,
    ("F", "H", 1):  567,
    ("Cl", "Cl", 1): 242,
    ("Cl", "H", 1): 431,
    ("Br", "Br", 1): 193,
    ("Br", "H", 1): 366,
    ("H", "I", 1):  297,
    ("I", "I", 1):  151,
    ("H", "S", 1):  363,
    ("S", "S", 1):  266,
    ("H", "P", 1):  322,
    ("O", "Si", 1): 452,
    ("H", "Si", 1): 318,
    # Double bonds
    ("C", "C", 2):  614,
    ("C", "N", 2):  615,
    ("C", "O", 2):  799,
    ("N", "N", 2):  418,
    ("O", "O", 2):  498,
    ("N", "O", 2):  607,
    ("C", "S", 2):  577,
    ("O", "S", 2):  522,
    # Additional single bonds (P2-DATA-2)
    ("C", "P", 1):  264,   # Luo 2007
    ("C", "Si", 1): 318,   # CRC Handbook
    ("N", "S", 1):  272,   # Benson estimates
    ("Cl", "N", 1): 200,   # chloramine / N-Cl
    ("F", "O", 1):  190,   # OF2
    ("O", "S", 1):  265,   # dimethyl sulfoxide (single)
    ("F", "S", 1):  327,   # SF6 type
    ("N", "P", 1):  230,   # P-N amine
    ("F", "Si", 1): 565,   # SiF4
    # Triple bonds
    ("C", "C", 3):  839,
    ("C", "N", 3):  891,
    ("N", "N", 3):  941,
    ("C", "O", 3):  1072,
}


# ===================================================================
# Torsion barrier parameters (kJ/mol) -- OPLS-AA
# ===================================================================

TORSION_BARRIERS: dict[str, dict[str, float]] = {
    # sp3-sp3 torsions (alkanes)
    "H_sp3_sp3_H": {"V1": 0.0,  "V2": 0.0,   "V3": 1.39},
    "C_sp3_sp3_C": {"V1": 2.73, "V2": -0.53,  "V3": 0.84},
    "C_sp3_sp3_H": {"V1": 0.0,  "V2": 0.0,    "V3": 0.76},
    "H_sp3_sp3_C": {"V1": 0.0,  "V2": 0.0,    "V3": 0.76},
    # sp2-sp3 torsions (P2-DATA-3, OPLS-AA Jorgensen et al. JACS 1996)
    "C_sp2_sp3_C": {"V1": 0.0,  "V2": 0.0,    "V3": 0.50},
    "C_sp2_sp3_H": {"V1": 0.0,  "V2": 0.0,    "V3": 0.50},
    "H_sp2_sp3_C": {"V1": 0.0,  "V2": 0.0,    "V3": 0.50},
    "H_sp2_sp3_H": {"V1": 0.0,  "V2": 0.0,    "V3": 0.50},
    "H_sp3_sp2_C": {"V1": 0.0,  "V2": 0.0,    "V3": 0.50},
    "C_sp3_sp2_C": {"V1": 0.0,  "V2": 0.0,    "V3": 0.50},
    "C_sp3_sp2_H": {"V1": 0.0,  "V2": 0.0,    "V3": 0.50},
    # sp2-sp2 torsions (conjugated systems)
    "C_sp2_sp2_C": {"V1": 0.0,  "V2": 12.0,   "V3": 0.0},
    "C_sp2_sp2_H": {"V1": 0.0,  "V2": 12.0,   "V3": 0.0},
    "H_sp2_sp2_C": {"V1": 0.0,  "V2": 12.0,   "V3": 0.0},
    "H_sp2_sp2_H": {"V1": 0.0,  "V2": 12.0,   "V3": 0.0},
    "default":     {"V1": 0.0,  "V2": 0.0,    "V3": 1.00},
}


# ===================================================================
# Electronegativity thresholds for bond polarity
# ===================================================================

NONPOLAR_THRESHOLD = 0.4
POLAR_COVALENT_MAX = 1.7


# ===================================================================
# Lookup functions
# ===================================================================

def bond_length(sym_a: str, sym_b: str, order: int = 1) -> float:
    """Look up standard bond length in Angstroms.

    Checks STANDARD_BOND_LENGTHS first, falls back to
    element_properties covalent-radii estimate.
    """
    a, b = sorted([sym_a, sym_b])
    key = (a, b, order)
    if key in STANDARD_BOND_LENGTHS:
        return STANDARD_BOND_LENGTHS[key]
    return estimated_bond_length_angstrom(sym_a, sym_b, order)


def bde_lookup(sym_a: str, sym_b: str, order: int = 1) -> float | None:
    """Look up mean bond dissociation energy in kJ/mol, or None."""
    a, b = sorted([sym_a, sym_b])
    return BDE_TABLE.get((a, b, order))
