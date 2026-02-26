"""Feature extraction for ML-based condition prediction.

Extracts numerical feature vectors from molecules and reaction context
that can be used to train ML models for predicting optimal reaction
conditions (temperature, solvent, catalyst, etc.).

The feature set includes:
- Molecular descriptors (MW, heavy atoms, rings, rotatable bonds, logP)
- Functional group one-hot encoding (24 binary features)
- Reaction template category one-hot encoding
- Scale features (log of production scale)

This module is fully functional and produces real feature dictionaries.
It can be used for offline training scripts and is integrated into
the ConditionPredictor stub in ml_predict.py.
"""

from __future__ import annotations

import math

from molbuilder.smiles.parser import parse
from molbuilder.molecule.graph import Molecule
from molbuilder.core.elements import atomic_weight
from molbuilder.molecule.properties import (
    heavy_atom_count,
    rotatable_bond_count,
    crippen_logp,
    hydrogen_bond_donors,
    hydrogen_bond_acceptors,
)
from molbuilder.reactions.functional_group_detect import detect_functional_groups
from molbuilder.reactions.reaction_types import ReactionCategory


def _molecular_weight(mol: Molecule) -> float:
    """Compute molecular weight from atom symbols."""
    return round(sum(atomic_weight(a.symbol) for a in mol.atoms), 4)


# The 24 FG detector names (stable ordering for one-hot encoding)
_FG_NAMES: list[str] = [
    "alcohol",
    "aldehyde",
    "alkene",
    "alkyl_halide_br",
    "alkyl_halide_cl",
    "alkyl_halide_f",
    "alkyl_halide_i",
    "alkyne",
    "amide",
    "anhydride",
    "aromatic_ring",
    "boronic_acid",
    "carboxylic_acid",
    "epoxide",
    "ester",
    "ether",
    "imine",
    "ketone",
    "nitrile",
    "nitro",
    "phosphonate",
    "primary_amine",
    "secondary_amine",
    "sulfonamide",
    "sulfone",
    "sulfoxide",
    "tertiary_amine",
    "thiol",
]

# Reaction category names (stable ordering for one-hot encoding)
_CATEGORY_NAMES: list[str] = [c.name for c in ReactionCategory]


def _count_rings(mol: Molecule) -> int:
    """Count the number of ring bonds (proxy for ring count)."""
    ring_bonds = 0
    for bond in mol.bonds:
        if mol.is_in_ring(bond.atom_i, bond.atom_j):
            ring_bonds += 1
    # Approximate ring count via cyclomatic number: E - V + 1
    # where E = ring_bonds, V = atoms in rings
    ring_atoms = set()
    for bond in mol.bonds:
        if mol.is_in_ring(bond.atom_i, bond.atom_j):
            ring_atoms.add(bond.atom_i)
            ring_atoms.add(bond.atom_j)
    if not ring_atoms:
        return 0
    return ring_bonds - len(ring_atoms) + 1


def extract_features(
    smiles: str,
    template_name: str | None = None,
    scale_kg: float = 1.0,
) -> dict[str, float]:
    """Extract ML-ready features from a molecule and reaction context.

    Parameters
    ----------
    smiles : str
        Substrate SMILES string.
    template_name : str | None
        Optional reaction template category name (e.g. "OXIDATION").
    scale_kg : float
        Production scale in kilograms.

    Returns
    -------
    dict[str, float]
        Feature dictionary with string keys and float values.
        All values are numeric (int/float). Boolean features use 0.0/1.0.
    """
    mol = parse(smiles)
    fgs = detect_functional_groups(mol)
    fg_names = {fg.name for fg in fgs}

    features: dict[str, float] = {}

    # --- Molecular descriptors ---
    features["mw"] = _molecular_weight(mol)
    features["heavy_atom_count"] = float(heavy_atom_count(mol))
    features["num_bonds"] = float(len(mol.bonds))
    features["num_rings"] = float(_count_rings(mol))
    features["rotatable_bonds"] = float(rotatable_bond_count(mol))
    features["logp"] = crippen_logp(mol)
    features["hbd"] = float(hydrogen_bond_donors(mol))
    features["hba"] = float(hydrogen_bond_acceptors(mol))

    # --- Functional group one-hot ---
    for fg_name in _FG_NAMES:
        features[f"fg_{fg_name}"] = 1.0 if fg_name in fg_names else 0.0

    # --- FG count (total number of detected groups) ---
    features["fg_count"] = float(len(fgs))

    # --- Category one-hot ---
    for cat_name in _CATEGORY_NAMES:
        if template_name and template_name.upper() == cat_name:
            features[f"cat_{cat_name.lower()}"] = 1.0
        else:
            features[f"cat_{cat_name.lower()}"] = 0.0

    # --- Scale feature ---
    features["log_scale_kg"] = math.log1p(scale_kg)

    return features


# Expose FG and category names for external use (e.g. training scripts)
FG_FEATURE_NAMES = [f"fg_{n}" for n in _FG_NAMES]
CATEGORY_FEATURE_NAMES = [f"cat_{n.lower()}" for n in _CATEGORY_NAMES]
DESCRIPTOR_NAMES = [
    "mw", "heavy_atom_count", "num_bonds", "num_rings",
    "rotatable_bonds", "logp", "hbd", "hba",
]
ALL_FEATURE_NAMES = DESCRIPTOR_NAMES + FG_FEATURE_NAMES + ["fg_count"] + CATEGORY_FEATURE_NAMES + ["log_scale_kg"]
