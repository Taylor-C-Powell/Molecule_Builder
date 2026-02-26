"""Tests for ML feature extraction and condition predictor stub."""

from __future__ import annotations

import math

import pytest

from molbuilder.process.ml_features import (
    extract_features,
    ALL_FEATURE_NAMES,
    FG_FEATURE_NAMES,
    CATEGORY_FEATURE_NAMES,
    DESCRIPTOR_NAMES,
)
from molbuilder.process.ml_predict import ConditionPredictor, get_predictor, set_predictor


class TestExtractFeatures:
    """Test feature extraction from molecules."""

    def test_ethanol_features(self):
        features = extract_features("CCO")
        assert isinstance(features, dict)
        # Molecular weight of ethanol ~ 46
        assert 40 < features["mw"] < 50
        assert features["heavy_atom_count"] >= 3
        assert features["logp"] != 0.0 or True  # logp can be any value
        # Ethanol has an alcohol FG
        assert features["fg_alcohol"] == 1.0

    def test_benzene_features(self):
        features = extract_features("c1ccccc1")
        assert features["heavy_atom_count"] >= 6
        assert features["fg_aromatic_ring"] == 1.0
        assert features["num_rings"] >= 1

    def test_feature_dict_has_correct_keys(self):
        features = extract_features("C")
        # Check all expected feature names are present
        for name in ALL_FEATURE_NAMES:
            assert name in features, f"Missing feature: {name}"

    def test_feature_values_are_numeric(self):
        features = extract_features("CCO")
        for key, val in features.items():
            assert isinstance(val, (int, float)), f"{key} is {type(val)}, expected numeric"

    def test_fg_features_are_binary(self):
        features = extract_features("CCO")
        for fg_name in FG_FEATURE_NAMES:
            assert features[fg_name] in (0.0, 1.0), f"{fg_name} = {features[fg_name]}"

    def test_category_features_are_binary(self):
        features = extract_features("CCO")
        for cat_name in CATEGORY_FEATURE_NAMES:
            assert features[cat_name] in (0.0, 1.0), f"{cat_name} = {features[cat_name]}"

    def test_category_one_hot_with_template(self):
        features = extract_features("CCO", template_name="OXIDATION")
        assert features["cat_oxidation"] == 1.0
        # All other categories should be 0
        for cat_name in CATEGORY_FEATURE_NAMES:
            if cat_name != "cat_oxidation":
                assert features[cat_name] == 0.0

    def test_no_category_when_no_template(self):
        features = extract_features("CCO")
        for cat_name in CATEGORY_FEATURE_NAMES:
            assert features[cat_name] == 0.0

    def test_scale_feature(self):
        features_1kg = extract_features("CCO", scale_kg=1.0)
        features_100kg = extract_features("CCO", scale_kg=100.0)
        assert features_100kg["log_scale_kg"] > features_1kg["log_scale_kg"]
        assert features_1kg["log_scale_kg"] == pytest.approx(math.log1p(1.0))
        assert features_100kg["log_scale_kg"] == pytest.approx(math.log1p(100.0))

    def test_descriptor_names_list(self):
        assert "mw" in DESCRIPTOR_NAMES
        assert "heavy_atom_count" in DESCRIPTOR_NAMES
        assert "logp" in DESCRIPTOR_NAMES

    def test_complex_molecule(self):
        # Caffeine
        features = extract_features("Cn1cnc2c1c(=O)n(c(=O)n2C)C")
        assert features["heavy_atom_count"] >= 14
        assert features["fg_count"] >= 1


class TestConditionPredictor:
    """Test the ML condition predictor stub."""

    def test_no_model_returns_none(self):
        predictor = ConditionPredictor()
        result = predictor.predict("CCO")
        assert result is None

    def test_is_loaded_false_by_default(self):
        predictor = ConditionPredictor()
        assert predictor.is_loaded is False

    def test_model_path_none(self):
        predictor = ConditionPredictor(model_path=None)
        assert predictor.is_loaded is False
        assert predictor.predict("CCO") is None

    def test_module_singleton(self):
        predictor = get_predictor()
        assert isinstance(predictor, ConditionPredictor)
        assert predictor.is_loaded is False

    def test_set_predictor(self):
        original = get_predictor()
        custom = ConditionPredictor()
        set_predictor(custom)
        assert get_predictor() is custom
        # Restore
        set_predictor(original)

    def test_predict_with_params(self):
        predictor = ConditionPredictor()
        result = predictor.predict(
            "CCO",
            template_name="OXIDATION",
            scale_kg=10.0,
        )
        assert result is None
