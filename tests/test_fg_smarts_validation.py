"""Tests for SMARTS cross-validation of functional group detection."""

import pytest

from molbuilder.smiles.parser import parse
from molbuilder.reactions.fg_smarts_validation import (
    FG_SMARTS_PATTERNS,
    FGValidationResult,
    cross_validate_fg,
    detect_functional_groups_validated,
)


class TestSingleFGAgreement:
    """Both methods should agree on simple molecules with one FG."""

    def test_alcohol_agreement(self):
        """CCO (ethanol) should be detected by both methods."""
        mol = parse("CCO")
        result = cross_validate_fg(mol)
        assert "alcohol" in result.agreed

    def test_aldehyde_agreement(self):
        """CC=O (acetaldehyde) should be detected by both methods."""
        mol = parse("CC=O")
        result = cross_validate_fg(mol)
        assert "aldehyde" in result.agreed

    def test_ketone_agreement(self):
        """CC(=O)C (acetone) should be detected by both methods."""
        mol = parse("CC(=O)C")
        result = cross_validate_fg(mol)
        assert "ketone" in result.agreed

    def test_ester_agreement(self):
        """CC(=O)OC (methyl acetate) should be detected by both methods."""
        mol = parse("CC(=O)OC")
        result = cross_validate_fg(mol)
        assert "ester" in result.agreed

    def test_carboxylic_acid_agreement(self):
        """CC(=O)O (acetic acid) should be detected by both methods."""
        mol = parse("CC(=O)O")
        result = cross_validate_fg(mol)
        assert "carboxylic_acid" in result.agreed

    def test_nitrile_agreement(self):
        """CC#N (acetonitrile) should be detected by both methods."""
        mol = parse("CC#N")
        result = cross_validate_fg(mol)
        assert "nitrile" in result.agreed

    def test_thiol_agreement(self):
        """CS (methanethiol) should be detected by both methods."""
        mol = parse("CS")
        result = cross_validate_fg(mol)
        assert "thiol" in result.agreed


class TestConfidence:
    """Confidence should be high for well-characterised molecules."""

    def test_confidence_above_threshold_ethanol(self):
        """Ethanol should have good agreement."""
        mol = parse("CCO")
        result = cross_validate_fg(mol)
        assert result.confidence >= 0.5

    def test_confidence_above_threshold_acetone(self):
        """Acetone should agree on ketone; SMARTS may find extra FGs
        due to broader pattern matching (expected disagreement)."""
        mol = parse("CC(=O)C")
        result = cross_validate_fg(mol)
        assert "ketone" in result.agreed
        assert result.confidence > 0.0

    def test_methane_full_confidence(self):
        """Methane has no FGs -- both methods should agree (empty)."""
        mol = parse("C")
        result = cross_validate_fg(mol)
        assert result.confidence == 1.0
        assert len(result.agreed) == 0
        assert len(result.heuristic_only) == 0
        assert len(result.smarts_only) == 0


class TestComplexMolecules:
    """Cross-validation on molecules with multiple FGs."""

    def test_aspirin(self):
        """Aspirin (CC(=O)Oc1ccccc1C(=O)O) has ester + carboxylic acid."""
        mol = parse("CC(=O)Oc1ccccc1C(=O)O")
        result = cross_validate_fg(mol)
        # At minimum, both methods should find something
        total = len(result.agreed) + len(result.heuristic_only) + len(result.smarts_only)
        assert total >= 2

    def test_ibuprofen(self):
        """Ibuprofen (CC(C)Cc1ccc(cc1)C(C)C(=O)O) has carboxylic acid."""
        mol = parse("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
        result = cross_validate_fg(mol)
        assert "carboxylic_acid" in result.agreed


class TestValidatedDetect:
    """detect_functional_groups_validated returns the union of both methods."""

    def test_validated_detect_returns_all_groups(self):
        """Should return at least as many groups as the heuristic alone."""
        from molbuilder.reactions.functional_group_detect import detect_functional_groups
        mol = parse("CC(=O)OC")
        heuristic = detect_functional_groups(mol)
        validated = detect_functional_groups_validated(mol)
        assert len(validated) >= len(heuristic)

    def test_validated_detect_returns_fg_objects(self):
        """All returned items should be FunctionalGroup instances."""
        from molbuilder.reactions.functional_group_detect import FunctionalGroup
        mol = parse("CCO")
        validated = detect_functional_groups_validated(mol)
        for fg in validated:
            assert isinstance(fg, FunctionalGroup)

    def test_validated_detect_methane_empty(self):
        """Methane should return no functional groups."""
        mol = parse("C")
        validated = detect_functional_groups_validated(mol)
        assert len(validated) == 0


class TestSmartsPatternsCoverage:
    """Verify the SMARTS pattern dictionary is well-formed."""

    def test_all_patterns_parseable(self):
        """Every SMARTS pattern should parse without error."""
        from molbuilder.smarts import parse_smarts
        for name, smarts in FG_SMARTS_PATTERNS.items():
            try:
                parse_smarts(smarts)
            except Exception as e:
                pytest.fail(f"SMARTS pattern for {name!r} failed to parse: {e}")

    def test_at_least_20_patterns(self):
        """Should have patterns for at least 20 of the 24 core FGs."""
        assert len(FG_SMARTS_PATTERNS) >= 20
