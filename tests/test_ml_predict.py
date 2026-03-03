"""Tests for ML-based condition prediction."""

import os
import sys
from unittest.mock import patch, MagicMock

import pytest

from molbuilder.process.ml_features import (
    extract_features,
    ALL_FEATURE_NAMES,
    FG_FEATURE_NAMES,
    CATEGORY_FEATURE_NAMES,
    DESCRIPTOR_NAMES,
)
from molbuilder.process.ml_predict import (
    ConditionPredictor,
    get_predictor,
    set_predictor,
)
from molbuilder.process.conditions import ReactionConditions


# =====================================================================
#  Feature extraction tests
# =====================================================================


class TestFeatureExtraction:
    """Tests for ml_features.extract_features()."""

    def test_feature_count(self):
        """extract_features should produce exactly len(ALL_FEATURE_NAMES) features."""
        features = extract_features("CCO")
        assert len(features) == len(ALL_FEATURE_NAMES)
        # 11 descriptors + 28 FG one-hot + 1 fg_count + 14 category one-hot + 1 log_scale
        assert len(ALL_FEATURE_NAMES) == 55

    def test_feature_keys_match(self):
        """Feature dict keys must match ALL_FEATURE_NAMES exactly."""
        features = extract_features("CCO")
        assert sorted(features.keys()) == sorted(ALL_FEATURE_NAMES)

    def test_descriptor_features_present(self):
        """Molecular descriptor features should be present and numeric."""
        features = extract_features("c1ccccc1")
        for name in DESCRIPTOR_NAMES:
            assert name in features
            assert isinstance(features[name], (int, float))

    def test_fg_features_binary(self):
        """Functional group features should be 0.0 or 1.0."""
        features = extract_features("CCO")
        for name in FG_FEATURE_NAMES:
            assert features[name] in (0.0, 1.0), f"{name} = {features[name]}"

    def test_alcohol_detected(self):
        """Ethanol should have fg_alcohol = 1.0."""
        features = extract_features("CCO")
        assert features["fg_alcohol"] == 1.0

    def test_category_features_one_hot(self):
        """Category one-hot should have at most one 1.0."""
        features = extract_features("CCO", template_name="OXIDATION")
        active = [n for n in CATEGORY_FEATURE_NAMES if features[n] == 1.0]
        assert len(active) == 1
        assert active[0] == "cat_oxidation"

    def test_category_none(self):
        """Without template_name, all category features should be 0.0."""
        features = extract_features("CCO")
        for name in CATEGORY_FEATURE_NAMES:
            assert features[name] == 0.0

    def test_scale_feature(self):
        """log_scale_kg should increase with scale."""
        f1 = extract_features("CCO", scale_kg=1.0)
        f10 = extract_features("CCO", scale_kg=10.0)
        assert f10["log_scale_kg"] > f1["log_scale_kg"]

    def test_tpsa_feature(self):
        """TPSA should be present and positive for ethanol."""
        features = extract_features("CCO")
        assert "tpsa" in features
        assert features["tpsa"] > 0

    def test_chiral_centers_feature(self):
        """Chiral SMILES should report nonzero num_chiral_centers."""
        features = extract_features("C[C@H](O)F")
        assert "num_chiral_centers" in features
        assert features["num_chiral_centers"] >= 1.0

    def test_chiral_centers_achiral(self):
        """Achiral molecule should have 0 chiral centers."""
        features = extract_features("CCO")
        assert features["num_chiral_centers"] == 0.0

    def test_ionizable_groups_feature(self):
        """Carboxylic acid should have ionizable groups."""
        features = extract_features("CC(=O)O")
        assert "num_ionizable_groups" in features
        assert features["num_ionizable_groups"] >= 1.0


# =====================================================================
#  ConditionPredictor tests
# =====================================================================


class TestConditionPredictor:
    """Tests for the ConditionPredictor class."""

    def test_no_model_returns_none(self):
        """Predictor with explicit nonexistent path should return None."""
        pred = ConditionPredictor.__new__(ConditionPredictor)
        pred._model = None
        pred._model_path = None
        assert not pred.is_loaded
        result = pred.predict("CCO")
        assert result is None

    def test_nonexistent_path_stays_unloaded(self):
        """Predictor with bad path should stay unloaded."""
        pred = ConditionPredictor(model_path="/nonexistent/model.pkl")
        assert not pred.is_loaded
        assert pred.predict("CCO") is None

    def test_import_error_graceful(self):
        """If joblib is not importable, predictor should stay unloaded."""
        import builtins
        original_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name == "joblib":
                raise ImportError("mocked")
            return original_import(name, *args, **kwargs)

        pred = ConditionPredictor.__new__(ConditionPredictor)
        pred._model = None
        pred._model_path = None

        with patch("builtins.__import__", side_effect=mock_import):
            pred._load_model("/some/path.pkl")

        assert not pred.is_loaded

    def test_feature_name_mismatch_rejects_model(self):
        """Model with wrong feature_names should be rejected."""
        mock_model = {
            "feature_names": ["wrong_feature_1", "wrong_feature_2"],
            "temperature_model": MagicMock(),
            "solvent_model": MagicMock(),
            "catalyst_model": MagicMock(),
            "yield_model": MagicMock(),
            "solvent_classes": ["THF"],
            "catalyst_classes": ["none"],
            "version": "1.0",
        }
        pred = ConditionPredictor.__new__(ConditionPredictor)
        pred._model = None
        pred._model_path = None

        with patch("joblib.load", return_value=mock_model):
            pred._load_model("/fake/model.pkl")

        assert not pred.is_loaded

    def test_valid_model_loads(self):
        """A model dict with correct feature_names should load."""
        temp_mock = MagicMock()
        temp_mock.predict.return_value = [80.0]
        solvent_mock = MagicMock()
        solvent_mock.predict.return_value = [0]
        catalyst_mock = MagicMock()
        catalyst_mock.predict.return_value = [0]
        yield_mock = MagicMock()
        yield_mock.predict.return_value = [75.0]

        mock_model = {
            "feature_names": ALL_FEATURE_NAMES,
            "temperature_model": temp_mock,
            "solvent_model": solvent_mock,
            "catalyst_model": catalyst_mock,
            "yield_model": yield_mock,
            "solvent_classes": ["THF", "DCM"],
            "catalyst_classes": ["Pd(PPh3)4", "none"],
            "version": "1.0",
        }
        pred = ConditionPredictor.__new__(ConditionPredictor)
        pred._model = None
        pred._model_path = None

        with patch("joblib.load", return_value=mock_model):
            pred._load_model("/fake/model.pkl")

        assert pred.is_loaded

    def test_predict_returns_conditions(self):
        """Loaded predictor should return ReactionConditions."""
        temp_mock = MagicMock()
        temp_mock.predict.return_value = [80.0]
        solvent_mock = MagicMock()
        solvent_mock.predict.return_value = [0]
        catalyst_mock = MagicMock()
        catalyst_mock.predict.return_value = [1]
        yield_mock = MagicMock()
        yield_mock.predict.return_value = [75.0]

        pred = ConditionPredictor.__new__(ConditionPredictor)
        pred._model = {
            "feature_names": ALL_FEATURE_NAMES,
            "temperature_model": temp_mock,
            "solvent_model": solvent_mock,
            "catalyst_model": catalyst_mock,
            "yield_model": yield_mock,
            "solvent_classes": ["THF", "DCM"],
            "catalyst_classes": ["Pd(PPh3)4", "none"],
            "version": "1.0",
        }
        pred._model_path = None

        result = pred.predict("CCO")
        assert isinstance(result, ReactionConditions)
        # Temperature model is gated (disabled); defaults to 25.0 C
        assert result.temperature_C == 25.0
        assert result.solvent == "THF"
        assert result.data_source.startswith("ML model")

    def test_temperature_gated(self):
        """Temperature model is disabled; should use heuristic default."""
        temp_mock = MagicMock()
        temp_mock.predict.return_value = [500.0]
        solvent_mock = MagicMock()
        solvent_mock.predict.return_value = [0]
        catalyst_mock = MagicMock()
        catalyst_mock.predict.return_value = [0]
        yield_mock = MagicMock()
        yield_mock.predict.return_value = [75.0]

        pred = ConditionPredictor.__new__(ConditionPredictor)
        pred._model = {
            "feature_names": ALL_FEATURE_NAMES,
            "temperature_model": temp_mock,
            "solvent_model": solvent_mock,
            "catalyst_model": catalyst_mock,
            "yield_model": yield_mock,
            "solvent_classes": ["THF"],
            "catalyst_classes": ["none"],
            "version": "1.0",
        }
        pred._model_path = None

        result = pred.predict("CCO")
        assert result is not None
        # Temperature model gated; temp_mock.predict should NOT be called
        temp_mock.predict.assert_not_called()
        assert result.temperature_C == 25.0

    def test_predict_exception_returns_none(self):
        """If model.predict raises, should return None gracefully."""
        solvent_mock = MagicMock()
        solvent_mock.predict.side_effect = RuntimeError("boom")

        pred = ConditionPredictor.__new__(ConditionPredictor)
        pred._model = {
            "feature_names": ALL_FEATURE_NAMES,
            "temperature_model": MagicMock(),
            "solvent_model": solvent_mock,
            "catalyst_model": MagicMock(),
            "yield_model": MagicMock(),
            "solvent_classes": ["THF"],
            "catalyst_classes": ["none"],
            "version": "1.0",
        }
        pred._model_path = None

        result = pred.predict("CCO")
        assert result is None


# =====================================================================
#  Singleton tests
# =====================================================================


class TestSingleton:
    """Tests for get_predictor/set_predictor."""

    def teardown_method(self):
        set_predictor(None)

    def test_get_predictor_returns_instance(self):
        """get_predictor should return a ConditionPredictor."""
        pred = get_predictor()
        assert isinstance(pred, ConditionPredictor)

    def test_set_predictor_overrides(self):
        """set_predictor should override the singleton."""
        custom = ConditionPredictor()
        set_predictor(custom)
        assert get_predictor() is custom


# =====================================================================
#  Integration with bundled model (if present)
# =====================================================================


class TestBundledModel:
    """Integration tests that run only if the trained model exists."""

    MODEL_PATH = os.path.join(
        os.path.dirname(__file__), "..", "molbuilder", "data", "condition_model.pkl"
    )

    def teardown_method(self):
        set_predictor(None)

    @pytest.mark.skipif(
        not os.path.isfile(
            os.path.join(
                os.path.dirname(__file__),
                "..", "molbuilder", "data", "condition_model.pkl",
            )
        ),
        reason="Trained model not present",
    )
    def test_bundled_model_loads(self):
        """Bundled model should auto-load via get_predictor."""
        set_predictor(None)
        pred = get_predictor()
        assert pred.is_loaded

    @pytest.mark.skipif(
        not os.path.isfile(
            os.path.join(
                os.path.dirname(__file__),
                "..", "molbuilder", "data", "condition_model.pkl",
            )
        ),
        reason="Trained model not present",
    )
    def test_bundled_predict_ethanol(self):
        """Bundled model should predict reasonable conditions for ethanol."""
        set_predictor(None)
        pred = get_predictor()
        result = pred.predict("CCO", template_name="OXIDATION")
        assert result is not None
        assert isinstance(result, ReactionConditions)
        assert -80 <= result.temperature_C <= 300
        assert result.solvent  # not empty
        assert result.data_source.startswith("ML model")

    @pytest.mark.skipif(
        not os.path.isfile(
            os.path.join(
                os.path.dirname(__file__),
                "..", "molbuilder", "data", "condition_model.pkl",
            )
        ),
        reason="Trained model not present",
    )
    def test_bundled_predict_benzene(self):
        """Bundled model should predict conditions for benzene."""
        set_predictor(None)
        pred = get_predictor()
        result = pred.predict("c1ccccc1", template_name="COUPLING")
        assert result is not None
        assert -80 <= result.temperature_C <= 300
