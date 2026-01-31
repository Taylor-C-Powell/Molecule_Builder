"""
Element Chemical Property Data

Centralised reference data for chemical bonding:
    - Pauling electronegativity
    - Single-bond covalent radii (Cordero 2008)
    - CPK (Corey-Pauling-Koltun) atom colours
"""

from bohr_model import ELEMENTS, SYMBOL_TO_Z

# ---------------------------------------------------------------------------
# Pauling electronegativity  (symbol -> float or None)
# Sources: IUPAC recommended values / Allred-Rochow / CRC Handbook
# None indicates no established value (noble gases, some superheavies)
# ---------------------------------------------------------------------------
PAULING_ELECTRONEGATIVITY = {
    "H":  2.20, "He": None,
    "Li": 0.98, "Be": 1.57, "B":  2.04, "C":  2.55, "N":  3.04,
    "O":  3.44, "F":  3.98, "Ne": None,
    "Na": 0.93, "Mg": 1.31, "Al": 1.61, "Si": 1.90, "P":  2.19,
    "S":  2.58, "Cl": 3.16, "Ar": None,
    "K":  0.82, "Ca": 1.00, "Sc": 1.36, "Ti": 1.54, "V":  1.63,
    "Cr": 1.66, "Mn": 1.55, "Fe": 1.83, "Co": 1.88, "Ni": 1.91,
    "Cu": 1.90, "Zn": 1.65, "Ga": 1.81, "Ge": 2.01, "As": 2.18,
    "Se": 2.55, "Br": 2.96, "Kr": 3.00,
    "Rb": 0.82, "Sr": 0.95, "Y":  1.22, "Zr": 1.33, "Nb": 1.60,
    "Mo": 2.16, "Tc": 1.90, "Ru": 2.20, "Rh": 2.28, "Pd": 2.20,
    "Ag": 1.93, "Cd": 1.69, "In": 1.78, "Sn": 1.96, "Sb": 2.05,
    "Te": 2.10, "I":  2.66, "Xe": 2.60,
    "Cs": 0.79, "Ba": 0.89, "La": 1.10, "Ce": 1.12, "Pr": 1.13,
    "Nd": 1.14, "Pm": 1.13, "Sm": 1.17, "Eu": 1.20, "Gd": 1.20,
    "Tb": 1.10, "Dy": 1.22, "Ho": 1.23, "Er": 1.24, "Tm": 1.25,
    "Yb": 1.10, "Lu": 1.27, "Hf": 1.30, "Ta": 1.50, "W":  2.36,
    "Re": 1.90, "Os": 2.20, "Ir": 2.20, "Pt": 2.28, "Au": 2.54,
    "Hg": 2.00, "Tl": 1.62, "Pb": 1.87, "Bi": 2.02, "Po": 2.00,
    "At": 2.20, "Rn": 2.20,
    "Fr": 0.70, "Ra": 0.90, "Ac": 1.10, "Th": 1.30, "Pa": 1.50,
    "U":  1.38, "Np": 1.36, "Pu": 1.28, "Am": 1.30, "Cm": 1.30,
    "Bk": 1.30, "Cf": 1.30, "Es": 1.30, "Fm": 1.30, "Md": 1.30,
    "No": 1.30, "Lr": 1.30,
    "Rf": None, "Db": None, "Sg": None, "Bh": None, "Hs": None,
    "Mt": None, "Ds": None, "Rg": None, "Cn": None, "Nh": None,
    "Fl": None, "Mc": None, "Lv": None, "Ts": None, "Og": None,
}

# ---------------------------------------------------------------------------
# Single-bond covalent radii in picometres (Cordero et al. 2008)
# ---------------------------------------------------------------------------
COVALENT_RADII_PM = {
    "H":  31,  "He": 28,
    "Li": 128, "Be": 96,  "B":  84,  "C":  76,  "N":  71,
    "O":  66,  "F":  57,  "Ne": 58,
    "Na": 166, "Mg": 141, "Al": 121, "Si": 111, "P":  107,
    "S":  105, "Cl": 102, "Ar": 106,
    "K":  203, "Ca": 176, "Sc": 170, "Ti": 160, "V":  153,
    "Cr": 139, "Mn": 139, "Fe": 132, "Co": 126, "Ni": 124,
    "Cu": 132, "Zn": 122, "Ga": 122, "Ge": 120, "As": 119,
    "Se": 120, "Br": 120, "Kr": 116,
    "Rb": 220, "Sr": 195, "Y":  190, "Zr": 175, "Nb": 164,
    "Mo": 154, "Tc": 147, "Ru": 146, "Rh": 142, "Pd": 139,
    "Ag": 145, "Cd": 144, "In": 142, "Sn": 139, "Sb": 139,
    "Te": 138, "I":  139, "Xe": 140,
    "Cs": 244, "Ba": 215, "La": 207, "Ce": 204, "Pr": 203,
    "Nd": 201, "Pm": 199, "Sm": 198, "Eu": 198, "Gd": 196,
    "Tb": 194, "Dy": 192, "Ho": 192, "Er": 189, "Tm": 190,
    "Yb": 187, "Lu": 187, "Hf": 175, "Ta": 170, "W":  162,
    "Re": 151, "Os": 144, "Ir": 141, "Pt": 136, "Au": 136,
    "Hg": 132, "Tl": 145, "Pb": 146, "Bi": 148, "Po": 140,
    "At": 150, "Rn": 150,
    "Fr": 260, "Ra": 221, "Ac": 215, "Th": 206, "Pa": 200,
    "U":  196, "Np": 190, "Pu": 187, "Am": 180, "Cm": 169,
}

# ---------------------------------------------------------------------------
# CPK colours (hex strings) for atom visualisation
# Based on Corey-Pauling-Koltun convention
# ---------------------------------------------------------------------------
CPK_COLORS = {
    "H":  "#FFFFFF",  "He": "#D9FFFF",
    "Li": "#CC80FF",  "Be": "#C2FF00",  "B":  "#FFB5B5",  "C":  "#909090",
    "N":  "#3050F8",  "O":  "#FF0D0D",  "F":  "#90E050",  "Ne": "#B3E3F5",
    "Na": "#AB5CF2",  "Mg": "#8AFF00",  "Al": "#BFA6A6",  "Si": "#F0C8A0",
    "P":  "#FF8000",  "S":  "#FFFF30",  "Cl": "#1FF01F",  "Ar": "#80D1E3",
    "K":  "#8F40D4",  "Ca": "#3DFF00",  "Sc": "#E6E6E6",  "Ti": "#BFC2C7",
    "V":  "#A6A6AB",  "Cr": "#8A99C7",  "Mn": "#9C7AC7",  "Fe": "#E06633",
    "Co": "#F090A0",  "Ni": "#50D050",  "Cu": "#C88033",  "Zn": "#7D80B0",
    "Ga": "#C28F8F",  "Ge": "#668F8F",  "As": "#BD80E3",  "Se": "#FFA100",
    "Br": "#A62929",  "Kr": "#5CB8D1",
    "Rb": "#702EB0",  "Sr": "#00FF00",  "Y":  "#94FFFF",  "Zr": "#94E0E0",
    "Nb": "#73C2C9",  "Mo": "#54B5B5",  "Tc": "#3B9E9E",  "Ru": "#248F8F",
    "Rh": "#0A7D8C",  "Pd": "#006985",  "Ag": "#C0C0C0",  "Cd": "#FFD98F",
    "In": "#A67573",  "Sn": "#668080",  "Sb": "#9E63B5",  "Te": "#D47A00",
    "I":  "#940094",  "Xe": "#429EB0",
    "default": "#FF69B4",
}

# ---------------------------------------------------------------------------
# Typical target electron counts (for Lewis structure octet rules)
# Elements with known incomplete octets
# ---------------------------------------------------------------------------
TARGET_ELECTRONS = {
    "H": 2, "He": 2, "Li": 2, "Be": 4, "B": 6,
}
# All others default to 8 (or expandable for period 3+)


# ===================================================================
# Helper functions
# ===================================================================

def electronegativity(symbol: str) -> float:
    """Return Pauling electronegativity for an element symbol.

    Returns 0.0 for elements without an established value (noble gases, etc.).
    """
    val = PAULING_ELECTRONEGATIVITY.get(symbol)
    if val is None:
        return 0.0
    return val


def covalent_radius_pm(symbol: str) -> float:
    """Return single-bond covalent radius in picometres.

    Falls back to 150 pm for elements not in the table.
    """
    return COVALENT_RADII_PM.get(symbol, 150.0)


def estimated_bond_length_pm(symbol_a: str, symbol_b: str,
                              bond_order: int = 1) -> float:
    """Estimate bond length in picometres from covalent radii.

    Single bond:  r(A) + r(B)
    Double bond:  ~0.87 * single
    Triple bond:  ~0.78 * single
    """
    single = covalent_radius_pm(symbol_a) + covalent_radius_pm(symbol_b)
    if bond_order == 2:
        return single * 0.87
    elif bond_order == 3:
        return single * 0.78
    return single


def estimated_bond_length_angstrom(symbol_a: str, symbol_b: str,
                                    bond_order: int = 1) -> float:
    """Estimate bond length in Angstroms (pm / 100)."""
    return estimated_bond_length_pm(symbol_a, symbol_b, bond_order) / 100.0


def cpk_color(symbol: str) -> str:
    """Return CPK hex colour for an element, with fallback."""
    return CPK_COLORS.get(symbol, CPK_COLORS["default"])


def target_electrons(symbol: str) -> int:
    """Target electron count for octet/duet rule.

    Returns 2 for H/He, 4 for Be, 6 for B, 8 for everything else.
    """
    return TARGET_ELECTRONS.get(symbol, 8)


def period(symbol: str) -> int:
    """Return the period (row) of an element."""
    z = SYMBOL_TO_Z.get(symbol, 0)
    if z <= 2:
        return 1
    if z <= 10:
        return 2
    if z <= 18:
        return 3
    if z <= 36:
        return 4
    if z <= 54:
        return 5
    if z <= 86:
        return 6
    return 7


def can_expand_octet(symbol: str) -> bool:
    """Whether an element can have an expanded octet (period 3+)."""
    return period(symbol) >= 3


# ===================================================================
# Main -- spot checks
# ===================================================================
if __name__ == "__main__":
    print("Element Data -- Spot Checks\n")

    test_elements = ["H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I", "Xe"]
    print(f"  {'Sym':<4} {'EN':>5}  {'Radius':>7}  {'Period':>6}  {'Expand':>7}  Color")
    print(f"  {'-'*50}")
    for sym in test_elements:
        en = electronegativity(sym)
        r = covalent_radius_pm(sym)
        p = period(sym)
        exp = can_expand_octet(sym)
        col = cpk_color(sym)
        print(f"  {sym:<4} {en:>5.2f}  {r:>5.0f} pm  {p:>5}  {'yes' if exp else 'no':>7}  {col}")

    print()
    print(f"  Bond length C-H single:  {estimated_bond_length_pm('C','H',1):.0f} pm")
    print(f"  Bond length C=O double:  {estimated_bond_length_pm('C','O',2):.0f} pm")
    print(f"  Bond length N#N triple:  {estimated_bond_length_pm('N','N',3):.0f} pm")
