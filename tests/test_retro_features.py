"""Tests for retrosynthetic disconnection feature extraction."""

import pytest

from molbuilder.smiles.parser import parse
from molbuilder.reactions.retro_features import (
    extract_retro_features,
    ALL_RETRO_FEATURE_NAMES,
    _TARGET_DESCRIPTOR_NAMES,
    _PRECURSOR_AGGREGATE_NAMES,
    _RELATIONSHIP_NAMES,
    _TARGET_FG_NAMES,
    _TEMPLATE_CATEGORY_NAMES,
    _STRATEGIC_FLAG_NAMES,
)
from molbuilder.reactions.retrosynthesis import Precursor
from molbuilder.reactions.functional_group_detect import detect_functional_groups
from molbuilder.reactions.reaction_types import ReactionTemplate, ReactionCategory


# =====================================================================
#  Helpers
# =====================================================================

def _make_template(**overrides):
    """Create a minimal ReactionTemplate for testing."""
    defaults = dict(
        name="Test Oxidation",
        named_reaction=None,
        category=ReactionCategory.OXIDATION,
        reagents=["KMnO4"],
        solvents=["water"],
        catalysts=[],
        temperature_range=(0, 100),
        typical_yield=(60, 80),
        functional_group_required=["alcohol"],
        functional_group_produced=["ketone"],
        functional_group_incompatible=[],
        mechanism="Oxidation mechanism.",
        reverse_transform="Reduce ketone to alcohol.",
        scale_notes="",
        safety_notes="",
    )
    defaults.update(overrides)
    return ReactionTemplate(**defaults)


def _make_precursor(smiles, name="test", cost=50.0):
    return Precursor(smiles=smiles, molecule=None, name=name,
                     cost_per_kg=cost)


def _ethanol_features(**kwargs):
    """Extract features for ethanol -> oxidation template."""
    mol = parse("CCO")
    tmpl = _make_template()
    precursors = [_make_precursor("CC=O", "acetaldehyde", 30.0)]
    return extract_retro_features(mol, "CCO", tmpl, precursors, **kwargs)


# =====================================================================
#  Feature structure tests
# =====================================================================


class TestFeatureStructure:
    """Tests for feature vector shape and naming."""

    def test_feature_count_is_90(self):
        """ALL_RETRO_FEATURE_NAMES should have exactly 90 entries."""
        assert len(ALL_RETRO_FEATURE_NAMES) == 90

    def test_feature_names_unique(self):
        """All feature names should be unique."""
        assert len(ALL_RETRO_FEATURE_NAMES) == len(set(ALL_RETRO_FEATURE_NAMES))

    def test_keys_match_all_names(self):
        """Extracted features must match ALL_RETRO_FEATURE_NAMES exactly."""
        features = _ethanol_features()
        assert sorted(features.keys()) == sorted(ALL_RETRO_FEATURE_NAMES)

    def test_all_values_numeric(self):
        """All feature values must be int or float."""
        features = _ethanol_features()
        for key, val in features.items():
            assert isinstance(val, (int, float)), f"{key} = {val!r}"


# =====================================================================
#  Block A: Target descriptor tests
# =====================================================================


class TestTargetDescriptors:
    """Tests for target molecule descriptor features."""

    def test_ethanol_mw_in_range(self):
        """Ethanol MW should be approximately 46."""
        features = _ethanol_features()
        assert 44.0 < features["target_mw"] < 48.0

    def test_ethanol_heavy_atoms(self):
        """Ethanol has 3 heavy atoms (C, C, O)."""
        features = _ethanol_features()
        assert features["target_heavy_atoms"] == 3.0

    def test_ethanol_sa_score_in_range(self):
        """SA score should be in [1, 10]."""
        features = _ethanol_features()
        assert 1.0 <= features["target_sa_score"] <= 10.0

    def test_ethanol_fg_alcohol_detected(self):
        """Ethanol should have tfg_alcohol = 1.0."""
        features = _ethanol_features()
        assert features["tfg_alcohol"] == 1.0


# =====================================================================
#  Block B: Precursor aggregate tests
# =====================================================================


class TestPrecursorAggregates:
    """Tests for precursor aggregate features."""

    def test_single_precursor_count(self):
        features = _ethanol_features()
        assert features["precursor_count"] == 1.0

    def test_multi_precursor_count(self):
        """Two precursors should give count = 2."""
        mol = parse("CCO")
        tmpl = _make_template(
            name="Test Coupling",
            category=ReactionCategory.COUPLING,
        )
        precursors = [
            _make_precursor("CC", "ethane", 20.0),
            _make_precursor("O", "water", 1.0),
        ]
        features = extract_retro_features(mol, "CCO", tmpl, precursors)
        assert features["precursor_count"] == 2.0

    def test_purchasable_fraction(self):
        """For ethanol precursors, purchasable fraction should be in [0, 1]."""
        features = _ethanol_features()
        assert 0.0 <= features["precursor_purchasable_frac"] <= 1.0


# =====================================================================
#  Block C: Relationship tests
# =====================================================================


class TestRelationshipFeatures:
    """Tests for relationship / complexity features."""

    def test_depth_reflects_input(self):
        """depth_in_tree should match the depth parameter."""
        f0 = _ethanol_features(depth=0)
        f3 = _ethanol_features(depth=3)
        assert f0["depth_in_tree"] == 0.0
        assert f3["depth_in_tree"] == 3.0


# =====================================================================
#  Block D: Template tests
# =====================================================================


class TestTemplateFeatures:
    """Tests for template features."""

    def test_category_one_hot_at_most_one(self):
        """Category one-hot should have at most one 1.0."""
        features = _ethanol_features()
        active = [n for n in _TEMPLATE_CATEGORY_NAMES if features[n] == 1.0]
        assert len(active) == 1
        assert active[0] == "tcat_oxidation"

    def test_yield_features(self):
        """Template yield features should match template values."""
        features = _ethanol_features()
        assert features["template_yield_lo"] == 60.0
        assert features["template_yield_hi"] == 80.0
        assert features["template_yield_mid"] == 70.0


# =====================================================================
#  Block E: FG one-hot tests
# =====================================================================


class TestFGOneHot:
    """Tests for target FG one-hot features."""

    def test_fg_features_binary(self):
        """All tfg_* features should be 0.0 or 1.0."""
        features = _ethanol_features()
        for name in _TARGET_FG_NAMES:
            assert features[name] in (0.0, 1.0), f"{name} = {features[name]}"


# =====================================================================
#  Block F: Strategic flag tests
# =====================================================================


class TestStrategicFlags:
    """Tests for strategic flag features."""

    def test_cc_coupling_flag_for_coupling(self):
        """Suzuki coupling should set is_cc_coupling = 1.0."""
        mol = parse("c1ccccc1")
        tmpl = _make_template(
            name="Suzuki coupling",
            named_reaction="Suzuki coupling",
            category=ReactionCategory.COUPLING,
        )
        precursors = [_make_precursor("c1ccc(Br)cc1", "bromobenzene")]
        features = extract_retro_features(mol, "c1ccccc1", tmpl, precursors)
        assert features["is_cc_coupling"] == 1.0

    def test_named_reaction_flag(self):
        """Template with named_reaction should set flag to 1.0."""
        mol = parse("CCO")
        tmpl = _make_template(named_reaction="Jones oxidation")
        precursors = [_make_precursor("CC=O")]
        features = extract_retro_features(mol, "CCO", tmpl, precursors)
        assert features["is_named_reaction"] == 1.0

    def test_no_named_reaction_flag(self):
        """Template without named_reaction should set flag to 0.0."""
        features = _ethanol_features()
        assert features["is_named_reaction"] == 0.0


# =====================================================================
#  Determinism test
# =====================================================================


class TestDeterminism:
    """Ensure feature extraction is deterministic."""

    def test_same_input_same_output(self):
        """Identical inputs should produce identical feature dicts."""
        f1 = _ethanol_features()
        f2 = _ethanol_features()
        assert f1 == f2


# =====================================================================
#  Pre-computed FG passthrough test
# =====================================================================


class TestPrecomputedFGs:
    """Test that pre-computed target_fgs are used correctly."""

    def test_precomputed_fgs_match(self):
        """Passing target_fgs should produce same result as auto-detection."""
        mol = parse("CCO")
        fgs = detect_functional_groups(mol)
        tmpl = _make_template()
        precursors = [_make_precursor("CC=O")]

        f_auto = extract_retro_features(mol, "CCO", tmpl, precursors)
        f_pre = extract_retro_features(mol, "CCO", tmpl, precursors,
                                       target_fgs=fgs)
        assert f_auto == f_pre
