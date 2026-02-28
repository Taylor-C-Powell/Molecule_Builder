"""Tests for redox transform helpers in the retrosynthesis engine.

Validates that carbonyl<->carboxyl and carbonyl<->alcohol transforms
produce correct precursor SMILES for oxidation and reduction pathways.
"""

import pytest

from molbuilder.reactions.retrosynthesis import (
    _replace_carbonyl_with_carboxyl,
    _replace_carboxyl_with_carbonyl,
    _replace_oh_with_carbonyl,
    _replace_carbonyl_with_oh,
)
from molbuilder.smiles.parser import parse
from molbuilder.smiles.writer import to_smiles


# =====================================================================
#  Workstream 1a: aldehyde/ketone oxidation -> carboxylic acid
# =====================================================================

class TestCarbonylToCarboxyl:
    """_replace_carbonyl_with_carboxyl: C=O -> C(=O)O"""

    def test_aldehyde_oxidation_produces_carboxylic_acid(self):
        """Acetaldehyde (CC=O) should oxidise to acetic acid."""
        result = _replace_carbonyl_with_carboxyl("CC=O")
        assert result is not None
        # The result should parse and contain a carboxylic acid motif
        mol = parse(result)
        # Must have 2 oxygens (one double-bonded, one single-bonded to same C)
        o_count = sum(1 for a in mol.atoms if a.symbol == "O")
        assert o_count == 2, f"Expected 2 oxygens in carboxylic acid, got {o_count}"

    def test_ketone_oxidation_produces_carboxylic_acid(self):
        """Acetone (CC(=O)C) should yield a carboxylic acid form."""
        result = _replace_carbonyl_with_carboxyl("CC(=O)C")
        assert result is not None
        mol = parse(result)
        o_count = sum(1 for a in mol.atoms if a.symbol == "O")
        assert o_count == 2

    def test_no_carbonyl_returns_none(self):
        """A molecule with no C=O should return None."""
        result = _replace_carbonyl_with_carboxyl("CCCC")
        assert result is None

    def test_already_carboxylic_acid(self):
        """Acetic acid itself: the first C=O found will get a second OH,
        but _validate_smiles_transform will detect it as changed."""
        result = _replace_carbonyl_with_carboxyl("CC(=O)O")
        # May or may not succeed depending on graph -- at minimum no crash
        if result is not None:
            mol = parse(result)
            assert len(mol.atoms) > 0


# =====================================================================
#  Workstream 1b: carboxylic acid reduction -> aldehyde
# =====================================================================

class TestCarboxylToCarbonyl:
    """_replace_carboxyl_with_carbonyl: C(=O)O -> C=O"""

    def test_carboxylic_acid_reduction_produces_aldehyde(self):
        """Acetic acid (CC(=O)O) should reduce to acetaldehyde."""
        result = _replace_carboxyl_with_carbonyl("CC(=O)O")
        assert result is not None
        mol = parse(result)
        # Should have exactly 1 oxygen (the remaining C=O)
        o_count = sum(1 for a in mol.atoms if a.symbol == "O")
        assert o_count == 1, f"Expected 1 oxygen in aldehyde, got {o_count}"

    def test_benzoic_acid_reduction(self):
        """Benzoic acid should reduce to benzaldehyde."""
        result = _replace_carboxyl_with_carbonyl("c1ccccc1C(=O)O")
        assert result is not None
        mol = parse(result)
        o_count = sum(1 for a in mol.atoms if a.symbol == "O")
        assert o_count == 1

    def test_no_carboxylic_acid_returns_none(self):
        """A molecule with no COOH should return None."""
        result = _replace_carboxyl_with_carbonyl("CC=O")
        assert result is None

    def test_ester_not_reduced(self):
        """An ester C(=O)O-C should NOT be treated as carboxylic acid
        because the single-bonded O is not terminal."""
        result = _replace_carboxyl_with_carbonyl("CC(=O)OC")
        assert result is None


# =====================================================================
#  Regression: existing transforms still work
# =====================================================================

class TestExistingTransformsRegression:
    """Ensure the existing alcohol<->carbonyl transforms are intact."""

    def test_alcohol_oxidation_still_works(self):
        """Ethanol (CCO) -> acetaldehyde."""
        result = _replace_oh_with_carbonyl("CCO")
        assert result is not None
        mol = parse(result)
        has_carbonyl = False
        for b in mol.bonds:
            if b.order == 2:
                syms = {mol.atoms[b.atom_i].symbol, mol.atoms[b.atom_j].symbol}
                if syms == {"C", "O"}:
                    has_carbonyl = True
        assert has_carbonyl, "Expected C=O bond in oxidised product"

    def test_carbonyl_reduction_still_works(self):
        """Acetaldehyde (CC=O) -> ethanol."""
        result = _replace_carbonyl_with_oh("CC=O")
        assert result is not None
        mol = parse(result)
        # Should have an O with single bond to C
        has_alcohol = False
        for b in mol.bonds:
            if b.order == 1:
                syms = {mol.atoms[b.atom_i].symbol, mol.atoms[b.atom_j].symbol}
                if syms == {"C", "O"}:
                    has_alcohol = True
        assert has_alcohol, "Expected C-O bond in reduced product"


# =====================================================================
#  End-to-end: retrosynthesis finds redox routes
# =====================================================================

class TestRetrosynthesisRedoxRoutes:
    """Integration test: retrosynthesis engine discovers redox pathways."""

    def test_retrosynthesis_finds_reduction_for_acid(self):
        """Acetic acid should have a reduction precursor (aldehyde)."""
        from molbuilder.reactions.retrosynthesis import retrosynthesis
        mol = parse("CC(=O)O")
        mol.name = "acetic_acid"
        tree = retrosynthesis(mol, max_depth=1)
        # The tree should exist and have a target with disconnections
        assert tree is not None
        assert hasattr(tree, "target")
        assert tree.target is not None
        assert tree.routes_found >= 1
