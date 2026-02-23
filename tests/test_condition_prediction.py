"""Tests for substrate-aware reaction condition prediction."""

import pytest

from molbuilder.process.condition_prediction import (
    predict_conditions,
    ConditionPrediction,
    SubstrateAnalysis,
    TemplateMatch,
    SolventScore,
    _analyze_substrate,
    _classify_steric_environment,
    _classify_electronic_character,
    _identify_sensitive_groups,
    _score_template,
    _adjust_yield_range,
    _build_warnings,
    _gather_candidate_templates,
)
from molbuilder.smiles.parser import parse
from molbuilder.reactions.functional_group_detect import detect_functional_groups


# =====================================================================
#  Helper
# =====================================================================

def _parse_and_analyze(smiles: str) -> tuple:
    """Parse SMILES, detect FGs, return (mol, fgs, substrate)."""
    mol = parse(smiles)
    fgs = detect_functional_groups(mol)
    substrate = _analyze_substrate(mol, fgs)
    return mol, fgs, substrate


# =====================================================================
#  predict_conditions() end-to-end tests
# =====================================================================

class TestPredictConditions:
    """End-to-end tests for the public predict_conditions() function."""

    def test_ethanol_returns_prediction(self):
        result = predict_conditions("CCO")
        assert isinstance(result, ConditionPrediction)
        assert result.smiles == "CCO"
        assert isinstance(result.substrate_analysis, SubstrateAnalysis)
        assert result.overall_confidence in ("high", "medium", "low")

    def test_ethanol_finds_candidates(self):
        result = predict_conditions("CCO")
        assert len(result.candidates) > 0
        assert result.best_match is not None

    def test_ethanol_oxidation_hint(self):
        result = predict_conditions("CCO", reaction_name="oxidation")
        assert result.best_match is not None
        # Should find oxidation-related templates
        categories = [c.category for c in result.candidates]
        assert "OXIDATION" in categories

    def test_alkyl_halide_substrate(self):
        # Chloroethane -- primary alkyl halide
        result = predict_conditions("CCCl")
        assert result.best_match is not None
        # Should find substitution templates
        categories = [c.category for c in result.candidates]
        assert "SUBSTITUTION" in categories

    def test_ketone_substrate(self):
        # Acetone
        result = predict_conditions("CC(=O)C")
        assert result.best_match is not None
        assert "ketone" in result.substrate_analysis.detected_functional_groups

    def test_aromatic_substrate(self):
        # Toluene
        result = predict_conditions("Cc1ccccc1")
        assert "aromatic_ring" in result.substrate_analysis.detected_functional_groups

    def test_max_candidates_respected(self):
        result = predict_conditions("CCO", max_candidates=2)
        assert len(result.candidates) <= 2

    def test_max_candidates_one(self):
        result = predict_conditions("CCO", max_candidates=1)
        assert len(result.candidates) <= 1

    def test_scale_affects_conditions(self):
        small = predict_conditions("CCO", reaction_name="oxidation", scale_kg=0.1)
        large = predict_conditions("CCO", reaction_name="oxidation", scale_kg=100.0)
        if small.best_match and large.best_match:
            # Larger scale should have different addition rate or time
            assert (
                small.best_match.conditions.addition_rate
                != large.best_match.conditions.addition_rate
                or small.best_match.conditions.reaction_time_hours
                != large.best_match.conditions.reaction_time_hours
            )

    def test_invalid_smiles_raises(self):
        with pytest.raises((ValueError, KeyError)):
            predict_conditions("INVALID!!!")

    def test_ibuprofen_substrate(self):
        # Ibuprofen SMILES
        result = predict_conditions("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
        assert result.substrate_analysis.molecular_weight > 100
        assert result.substrate_analysis.heavy_atom_count > 10
        assert "carboxylic_acid" in result.substrate_analysis.detected_functional_groups

    def test_best_match_has_conditions(self):
        result = predict_conditions("CCO", reaction_name="oxidation")
        if result.best_match:
            cond = result.best_match.conditions
            assert cond.temperature_C is not None
            assert cond.solvent is not None
            assert cond.atmosphere in ("air", "N2", "Ar")

    def test_best_match_has_solvents(self):
        result = predict_conditions("CCO", reaction_name="oxidation")
        if result.best_match:
            solvents = result.best_match.recommended_solvents
            assert isinstance(solvents, list)
            if solvents:
                s = solvents[0]
                assert isinstance(s, SolventScore)
                assert 0 <= s.composite_score <= 100
                assert 1 <= s.green_score <= 10

    def test_overall_confidence_values(self):
        result = predict_conditions("CCO", reaction_name="oxidation")
        assert result.overall_confidence in ("high", "medium", "low")

    def test_adjusted_yield_range(self):
        result = predict_conditions("CCO", reaction_name="oxidation")
        if result.best_match:
            lo, hi = result.best_match.adjusted_yield_range
            assert 0 <= lo <= hi <= 100


# =====================================================================
#  Substrate analysis tests
# =====================================================================

class TestSubstrateAnalysis:
    """Tests for _analyze_substrate and its helpers."""

    def test_ethanol_is_unhindered(self):
        _, fgs, substrate = _parse_and_analyze("CCO")
        assert substrate.steric_class == "unhindered"

    def test_tert_butanol_is_hindered(self):
        # 2-methyl-2-propanol (tert-butanol)
        _, fgs, substrate = _parse_and_analyze("CC(C)(C)O")
        assert substrate.steric_class in ("moderately_hindered", "hindered")

    def test_electronic_character_alcohol(self):
        _, fgs, substrate = _parse_and_analyze("CCO")
        # Alcohol is EDG, should be electron_rich or neutral
        assert substrate.electronic_character in ("electron_rich", "neutral")

    def test_electronic_character_dinitro(self):
        # 1,3-dinitrobenzene -- strongly electron-poor
        _, fgs, substrate = _parse_and_analyze("c1cc([N+](=O)[O-])cc([N+](=O)[O-])c1")
        assert substrate.electronic_character in ("electron_poor", "neutral")

    def test_sensitive_groups_aldehyde(self):
        # Acetaldehyde
        _, fgs, substrate = _parse_and_analyze("CC=O")
        assert "aldehyde" in substrate.sensitive_groups

    def test_no_sensitive_groups_ketone(self):
        # Acetone -- ketone is not in the sensitive set
        _, fgs, substrate = _parse_and_analyze("CC(=O)C")
        assert "ketone" not in substrate.sensitive_groups

    def test_heavy_atom_count(self):
        _, _, substrate = _parse_and_analyze("CCO")
        assert substrate.heavy_atom_count == 3  # C, C, O

    def test_molecular_weight_ethanol(self):
        _, _, substrate = _parse_and_analyze("CCO")
        # Ethanol: MW ~ 46 (with implicit H it includes H atoms)
        assert substrate.molecular_weight > 40


# =====================================================================
#  Scoring tests
# =====================================================================

class TestScoring:
    """Tests for _score_template."""

    def test_score_returns_float_and_reasons(self):
        mol, fgs, substrate = _parse_and_analyze("CCO")
        fg_names = [fg.name for fg in fgs]
        templates = _gather_candidate_templates(fg_names, "oxidation")
        if templates:
            sc, reasons = _score_template(templates[0], fg_names, substrate, "oxidation")
            assert isinstance(sc, float)
            assert 0 <= sc <= 100
            assert len(reasons) > 0

    def test_name_hint_boosts_score(self):
        mol, fgs, substrate = _parse_and_analyze("CCO")
        fg_names = [fg.name for fg in fgs]
        templates = _gather_candidate_templates(fg_names, "oxidation")
        if templates:
            tmpl = templates[0]
            sc_with, _ = _score_template(tmpl, fg_names, substrate, "oxidation")
            sc_without, _ = _score_template(tmpl, fg_names, substrate, None)
            # Score with hint should be >= score without
            assert sc_with >= sc_without

    def test_incompatible_fg_penalty(self):
        mol, fgs, substrate = _parse_and_analyze("CCO")
        fg_names = [fg.name for fg in fgs]
        templates = _gather_candidate_templates(fg_names, None)
        # Find a template with incompatible FGs if possible
        for tmpl in templates:
            if tmpl.functional_group_incompatible:
                # Add a fake incompatible FG to the substrate
                fake_names = fg_names + [tmpl.functional_group_incompatible[0].lower()]
                sc_bad, reasons = _score_template(tmpl, fake_names, substrate, None)
                sc_good, _ = _score_template(tmpl, fg_names, substrate, None)
                assert sc_bad < sc_good
                break


# =====================================================================
#  Yield adjustment tests
# =====================================================================

class TestYieldAdjustment:
    """Tests for _adjust_yield_range."""

    def test_hindered_reduces_yield(self):
        from molbuilder.reactions.knowledge_base import REACTION_TEMPLATES
        if not REACTION_TEMPLATES:
            pytest.skip("No templates available")
        tmpl = REACTION_TEMPLATES[0]
        # Fake substrate analyses
        unhindered = SubstrateAnalysis(
            detected_functional_groups=[], molecular_weight=100,
            heavy_atom_count=10, steric_class="unhindered",
            electronic_character="neutral", sensitive_groups=[],
        )
        hindered = SubstrateAnalysis(
            detected_functional_groups=[], molecular_weight=100,
            heavy_atom_count=10, steric_class="hindered",
            electronic_character="neutral", sensitive_groups=[],
        )
        lo_u, hi_u = _adjust_yield_range(tmpl, unhindered)
        lo_h, hi_h = _adjust_yield_range(tmpl, hindered)
        assert lo_h <= lo_u
        assert hi_h <= hi_u


# =====================================================================
#  Warnings tests
# =====================================================================

class TestWarnings:
    """Tests for _build_warnings."""

    def test_hindered_warning(self):
        from molbuilder.reactions.knowledge_base import REACTION_TEMPLATES
        if not REACTION_TEMPLATES:
            pytest.skip("No templates available")
        tmpl = REACTION_TEMPLATES[0]
        substrate = SubstrateAnalysis(
            detected_functional_groups=[], molecular_weight=100,
            heavy_atom_count=10, steric_class="hindered",
            electronic_character="neutral", sensitive_groups=[],
        )
        warnings = _build_warnings(tmpl, substrate)
        assert any("hindered" in w.lower() for w in warnings)

    def test_high_mw_warning(self):
        from molbuilder.reactions.knowledge_base import REACTION_TEMPLATES
        if not REACTION_TEMPLATES:
            pytest.skip("No templates available")
        tmpl = REACTION_TEMPLATES[0]
        substrate = SubstrateAnalysis(
            detected_functional_groups=[], molecular_weight=600,
            heavy_atom_count=40, steric_class="unhindered",
            electronic_character="neutral", sensitive_groups=[],
        )
        warnings = _build_warnings(tmpl, substrate)
        assert any("mw" in w.lower() or "high" in w.lower() for w in warnings)

    def test_sensitive_group_warning(self):
        from molbuilder.reactions.knowledge_base import REACTION_TEMPLATES
        if not REACTION_TEMPLATES:
            pytest.skip("No templates available")
        tmpl = REACTION_TEMPLATES[0]
        substrate = SubstrateAnalysis(
            detected_functional_groups=["aldehyde"], molecular_weight=100,
            heavy_atom_count=5, steric_class="unhindered",
            electronic_character="neutral", sensitive_groups=["aldehyde"],
        )
        warnings = _build_warnings(tmpl, substrate)
        assert any("aldehyde" in w for w in warnings)
