"""Peptide construction: bond formation, phi/psi, secondary structure.

The implementation lives in amino_acids.py. This module provides
convenient direct imports.
"""
from molbuilder.molecule.amino_acids import (
    form_peptide_bond,
    build_peptide,
    set_phi_psi,
    apply_secondary_structure,
    SecondaryStructure,
    BackboneIndices,
)
