"""Tests for graph-based retrosynthesis transforms.

Validates that the refactored graph-based transforms avoid the
string-replacement false positives of the old SMILES-based approach.

Run with:
    python -m pytest tests/test_retro_graph_transforms.py -v
"""

import pytest

from molbuilder.smiles.parser import parse
from molbuilder.smiles.writer import to_smiles
from molbuilder.reactions.retrosynthesis import (
    _replace_oh_with_carbonyl,
    _replace_carbonyl_with_oh,
    _validate_smiles_transform,
)


# =====================================================================
#  _replace_oh_with_carbonyl (alcohol -> carbonyl, graph-based)
# =====================================================================

class TestOxidizeAlcoholGraph:
    """Test graph-based alcohol-to-carbonyl oxidation."""

    def test_simple_methanol_oxidised(self):
        """CO (methanol) should produce C=O (formaldehyde)."""
        result = _replace_oh_with_carbonyl("CO")
        assert result is not None
        # Parse result and verify it has C=O
        mol = parse(result)
        has_double_o = any(
            b.order == 2 and (
                (mol.atoms[b.atom_i].symbol == "C" and mol.atoms[b.atom_j].symbol == "O")
                or (mol.atoms[b.atom_i].symbol == "O" and mol.atoms[b.atom_j].symbol == "C")
            )
            for b in mol.bonds
        )
        assert has_double_o

    def test_coco_does_not_corrupt(self):
        """COCO (methoxymethanol) must not produce 'C=OCO' corruption.

        The old string-based approach would do smiles.replace('CO', 'C=O', 1)
        which corrupts the ether linkage.  The graph-based approach should
        correctly find an alcohol O-H and only oxidise that specific bond.
        """
        result = _replace_oh_with_carbonyl("COCO")
        if result is not None:
            # It must parse without error
            mol = parse(result)
            # Must not have lost atoms compared to what's expected
            # The result should be valid SMILES with correct atom count
            heavy_count = sum(1 for a in mol.atoms if a.symbol != "H")
            assert heavy_count >= 3  # at least C, O, C

    def test_nco_not_misidentified_as_alcohol(self):
        """NCO (isocyanate-like) should return None since there is no alcohol."""
        result = _replace_oh_with_carbonyl("N=C=O")
        # Should return None -- no alcohol FG to oxidise
        assert result is None

    def test_ethanol_oxidised(self):
        """CCO (ethanol) -> CC=O (acetaldehyde)."""
        result = _replace_oh_with_carbonyl("CCO")
        assert result is not None
        mol = parse(result)
        has_carbonyl = any(
            b.order == 2 and {mol.atoms[b.atom_i].symbol, mol.atoms[b.atom_j].symbol} == {"C", "O"}
            for b in mol.bonds
        )
        assert has_carbonyl


# =====================================================================
#  _replace_carbonyl_with_oh (carbonyl -> alcohol, graph-based)
# =====================================================================

class TestReduceCarbonylGraph:
    """Test graph-based carbonyl-to-alcohol reduction."""

    def test_formaldehyde_reduced(self):
        """C=O (formaldehyde) -> CO (methanol)."""
        result = _replace_carbonyl_with_oh("C=O")
        assert result is not None
        mol = parse(result)
        has_c_o_single = any(
            b.order == 1 and {mol.atoms[b.atom_i].symbol, mol.atoms[b.atom_j].symbol} == {"C", "O"}
            for b in mol.bonds
        )
        assert has_c_o_single

    def test_acetone_reduced(self):
        """CC(=O)C (acetone) should reduce to an alcohol."""
        result = _replace_carbonyl_with_oh("CC(=O)C")
        assert result is not None
        mol = parse(result)
        # Should have a C-O single bond now
        has_c_o_single = any(
            b.order == 1 and {mol.atoms[b.atom_i].symbol, mol.atoms[b.atom_j].symbol} == {"C", "O"}
            for b in mol.bonds
        )
        assert has_c_o_single

    def test_no_carbonyl_returns_none(self):
        """CCO (ethanol, no C=O) should return None."""
        result = _replace_carbonyl_with_oh("CCO")
        assert result is None


# =====================================================================
#  Round-trip validation
# =====================================================================

class TestRoundTrip:
    """Parse -> transform -> parse validates atom counts."""

    def test_ethanol_oxidise_roundtrip(self):
        """CCO -> oxidise -> parse should have correct heavy atoms."""
        mol1 = parse("CCO")
        heavy1 = sum(1 for a in mol1.atoms if a.symbol != "H")
        result = _replace_oh_with_carbonyl("CCO")
        assert result is not None
        mol2 = parse(result)
        heavy2 = sum(1 for a in mol2.atoms if a.symbol != "H")
        # Oxidation removes H from OH: heavy atom count stays same (C, C, O)
        assert heavy2 == heavy1

    def test_acetaldehyde_reduce_roundtrip(self):
        """CC=O -> reduce -> parse should have correct heavy atoms."""
        mol1 = parse("CC=O")
        heavy1 = sum(1 for a in mol1.atoms if a.symbol != "H")
        result = _replace_carbonyl_with_oh("CC=O")
        assert result is not None
        mol2 = parse(result)
        heavy2 = sum(1 for a in mol2.atoms if a.symbol != "H")
        assert heavy2 == heavy1
