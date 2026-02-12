"""Data tables for molecular property calculations.

Wildman-Crippen logP atom contributions (simplified) and Ertl TPSA
fragment contributions.  Used by ``properties.py``.

References
----------
- Wildman & Crippen, J. Chem. Inf. Comput. Sci. 1999, 39, 868-873.
- Ertl, Rohde & Selzer, J. Med. Chem. 2000, 43, 3714-3717.
"""

from __future__ import annotations

# =====================================================================
#  Wildman-Crippen logP atom-type contributions
# =====================================================================
#
# Keys are (element, description) tuples used by the classifier in
# properties.py.  Values are the logP contribution for that atom type.
# This is a simplified scheme (~30 types) covering common drug-like atoms.

CRIPPEN_LOGP: dict[str, float] = {
    # --- Carbon types ---
    "C_sp3_all_C_H":       0.1441,   # CH3, CH2, CH in alkyl chains
    "C_sp3_hetero":       -0.0516,   # C bonded to N, O, S, or halogen
    "C_sp2_aromatic":      0.1305,   # aromatic C (no heteroatom neighbor)
    "C_sp2_aromatic_het":  0.0800,   # aromatic C next to heteroatom
    "C_sp2_alkene":        0.0800,   # non-aromatic C=C
    "C_sp2_carbonyl":     -0.1032,   # C=O (aldehyde, ketone, acid, ester, amide)
    "C_sp":               -0.0170,   # triple-bonded C (nitrile, alkyne)

    # --- Nitrogen types ---
    "N_primary_amine":    -0.7827,   # -NH2
    "N_secondary_amine":  -0.4458,   # -NH-
    "N_tertiary_amine":   -0.2840,   # -N<
    "N_aromatic":         -0.3239,   # aromatic N (pyridine)
    "N_amide":            -0.4806,   # -N(H)-C(=O)- or -N<-C(=O)-
    "N_nitrile":          -0.2893,   # C#N
    "N_nitro":            -0.1820,   # -NO2

    # --- Oxygen types ---
    "O_alcohol":          -0.3567,   # -OH (alcohol, phenol)
    "O_ether":            -0.2893,   # -O- (ether, furan)
    "O_carbonyl":         -0.3339,   # =O on C
    "O_carboxyl_oh":      -0.3567,   # -OH of -COOH
    "O_ester":            -0.2893,   # -O- of ester
    "O_nitro":            -0.0684,   # =O of nitro

    # --- Sulfur types ---
    "S_thiol":            -0.0024,   # -SH
    "S_thioether":         0.6237,   # -S-
    "S_aromatic":          0.4000,   # aromatic S (thiophene)

    # --- Halogen types ---
    "F":                   0.4118,   # fluorine
    "Cl":                  0.6895,   # chlorine
    "Br":                  0.8813,   # bromine
    "I":                   1.0500,   # iodine

    # --- Hydrogen types ---
    "H_on_C":              0.1230,   # H on carbon
    "H_on_heteroatom":    -0.2677,   # H on N, O, S

    # --- Phosphorus ---
    "P_any":               0.2836,   # any P

    # --- Fallback ---
    "other":               0.0000,
}


# =====================================================================
#  Ertl TPSA fragment contributions (Angstroms^2)
# =====================================================================
#
# Keys are tuples: (element, n_heavy_neighbors, n_hydrogens, n_double_bonds,
#                    is_aromatic, formal_charge)
# Special value None means "any" for that position.

TPSA_CONTRIBUTIONS: dict[tuple, float] = {
    # --- Nitrogen ---
    # Primary amine: -NH2 (1 heavy neighbor, 2 H)
    ("N", 1, 2, 0, False, 0):   26.02,
    # Secondary amine: -NH- (2 heavy neighbors, 1 H)
    ("N", 2, 1, 0, False, 0):   12.36,  # corrected from original: was used as 12.03 in some refs
    # Tertiary amine: -N< (3 heavy neighbors, 0 H)
    ("N", 3, 0, 0, False, 0):    3.24,
    # =N- (imine, 1 heavy neighbor)
    ("N", 1, 0, 1, False, 0):   23.79,
    # =N- (imine, 2 heavy neighbors)
    ("N", 2, 0, 1, False, 0):   12.89,
    # Aromatic N: pyridine-like (2 heavy, 0 H)
    ("N", 2, 0, 0, True, 0):    12.89,
    # Aromatic NH: pyrrole-like (2 heavy, 1 H)
    ("N", 2, 1, 0, True, 0):    15.79,
    # Nitrile: C#N (1 heavy, 0 H, effectively 1 triple bond counted as 2 doubles)
    ("N", 1, 0, 2, False, 0):   23.79,
    # Amide: N-C(=O) treated as secondary amine above; context handles it
    # Nitro: =N(=O)- handled by sum of N + O contributions
    # Charged N: -NH3+
    ("N", 1, 3, 0, False, 1):   27.64,
    # Charged N: -NH2+- (secondary)
    ("N", 2, 2, 0, False, 1):   25.59,
    # Charged N: =NH+- (iminium)
    ("N", 2, 1, 1, False, 1):   23.47,

    # --- Oxygen ---
    # Hydroxyl: -OH (1 heavy neighbor, 1 H)
    ("O", 1, 1, 0, False, 0):   20.23,
    # Ether / ester oxygen: -O- (2 heavy neighbors, 0 H)
    ("O", 2, 0, 0, False, 0):    9.23,
    # Carbonyl: =O (1 heavy neighbor, 0 H)
    ("O", 1, 0, 1, False, 0):   17.07,
    # Aromatic O (furan)
    ("O", 2, 0, 0, True, 0):     9.23,
    # Charged O: -O- (carboxylate)
    ("O", 1, 0, 0, False, -1):  23.06,

    # --- Sulfur ---
    # Thiol: -SH (1 heavy, 1 H)
    ("S", 1, 1, 0, False, 0):   38.80,
    # Thioether: -S- (2 heavy, 0 H)
    ("S", 2, 0, 0, False, 0):   25.30,
    # Aromatic S (thiophene)
    ("S", 2, 0, 0, True, 0):    28.24,
    # Sulfoxide: S=O counted as S + O contributions
    # Sulfonyl: S(=O)2 similarly additive

    # --- Phosphorus ---
    # Phosphate-like: =P or -P
    ("P", 2, 0, 0, False, 0):   13.59,
    ("P", 3, 0, 0, False, 0):    8.81,
    ("P", 2, 0, 1, False, 0):   34.14,
}


# =====================================================================
#  Standard valences for implicit H inference
# =====================================================================

STANDARD_VALENCE: dict[str, list[int]] = {
    "B": [3], "C": [4], "N": [3, 5], "O": [2], "P": [3, 5],
    "S": [2, 4, 6], "F": [1], "Cl": [1], "Br": [1], "I": [1],
}
