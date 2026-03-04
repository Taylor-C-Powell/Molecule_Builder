"""Tests for ML-based retrosynthetic disconnection scoring."""

import os
from unittest.mock import patch, MagicMock

import pytest

from molbuilder.smiles.parser import parse
from molbuilder.reactions.retro_features import ALL_RETRO_FEATURE_NAMES
from molbuilder.reactions.ml_scoring import (
    DisconnectionScorer,
    get_scorer,
    set_scorer,
)
from molbuilder.reactions.retrosynthesis import (
    Precursor,
    Disconnection,
    retrosynthesis,
)
from molbuilder.reactions.reaction_types import ReactionTemplate, ReactionCategory


# =====================================================================
#  Helpers
# =====================================================================

def _make_template(**overrides):
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


def _mock_model_dict():
    """Create a valid mock model dict with correct feature names."""
    score_mock = MagicMock()
    score_mock.predict.return_value = [72.5]
    return {
        "score_model": score_mock,
        "feature_names": list(ALL_RETRO_FEATURE_NAMES),
        "version": "1.0",
        "r2_score": 0.85,
        "mae": 3.2,
    }


# =====================================================================
#  DisconnectionScorer class tests
# =====================================================================


class TestDisconnectionScorer:
    """Tests for the DisconnectionScorer class."""

    def test_no_model_returns_none(self):
        """Scorer with no model should return None."""
        scorer = DisconnectionScorer.__new__(DisconnectionScorer)
        scorer._model = None
        scorer._model_path = None
        assert not scorer.is_loaded
        mol = parse("CCO")
        result = scorer.score(mol, "CCO", _make_template(),
                              [_make_precursor("CC=O")])
        assert result is None

    def test_nonexistent_path_stays_unloaded(self):
        """Scorer with bad path should stay unloaded."""
        scorer = DisconnectionScorer(model_path="/nonexistent/model.pkl")
        assert not scorer.is_loaded

    def test_import_error_graceful(self):
        """If joblib is not importable, scorer should stay unloaded."""
        import builtins
        original_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name == "joblib":
                raise ImportError("mocked")
            return original_import(name, *args, **kwargs)

        scorer = DisconnectionScorer.__new__(DisconnectionScorer)
        scorer._model = None
        scorer._model_path = None

        with patch("builtins.__import__", side_effect=mock_import):
            scorer._load_model("/some/path.pkl")

        assert not scorer.is_loaded

    def test_feature_name_mismatch_rejects(self):
        """Model with wrong feature_names should be rejected."""
        mock_model = {
            "feature_names": ["wrong_feature_1", "wrong_feature_2"],
            "score_model": MagicMock(),
            "version": "1.0",
        }
        scorer = DisconnectionScorer.__new__(DisconnectionScorer)
        scorer._model = None
        scorer._model_path = None

        with patch("joblib.load", return_value=mock_model):
            scorer._load_model("/fake/model.pkl")

        assert not scorer.is_loaded

    def test_valid_model_loads(self):
        """A model dict with correct feature_names should load."""
        mock_model = _mock_model_dict()
        scorer = DisconnectionScorer.__new__(DisconnectionScorer)
        scorer._model = None
        scorer._model_path = "/fake/model.pkl"  # custom path -> no sidecar required

        with patch("joblib.load", return_value=mock_model):
            scorer._load_model("/fake/model.pkl")

        assert scorer.is_loaded

    def test_score_returns_clamped_float(self):
        """Loaded scorer should return a clamped float score."""
        mock_model = _mock_model_dict()
        scorer = DisconnectionScorer.__new__(DisconnectionScorer)
        scorer._model = mock_model

        mol = parse("CCO")
        result = scorer.score(mol, "CCO", _make_template(),
                              [_make_precursor("CC=O")])
        assert isinstance(result, float)
        assert 0.0 <= result <= 100.0

    def test_score_clamped_high(self):
        """Scores above 100 should be clamped to 100."""
        mock_model = _mock_model_dict()
        mock_model["score_model"].predict.return_value = [150.0]
        scorer = DisconnectionScorer.__new__(DisconnectionScorer)
        scorer._model = mock_model

        mol = parse("CCO")
        result = scorer.score(mol, "CCO", _make_template(),
                              [_make_precursor("CC=O")])
        assert result == 100.0

    def test_score_clamped_low(self):
        """Scores below 0 should be clamped to 0."""
        mock_model = _mock_model_dict()
        mock_model["score_model"].predict.return_value = [-10.0]
        scorer = DisconnectionScorer.__new__(DisconnectionScorer)
        scorer._model = mock_model

        mol = parse("CCO")
        result = scorer.score(mol, "CCO", _make_template(),
                              [_make_precursor("CC=O")])
        assert result == 0.0

    def test_prediction_exception_returns_none(self):
        """If model.predict raises, should return None gracefully."""
        mock_model = _mock_model_dict()
        mock_model["score_model"].predict.side_effect = RuntimeError("boom")
        scorer = DisconnectionScorer.__new__(DisconnectionScorer)
        scorer._model = mock_model

        mol = parse("CCO")
        result = scorer.score(mol, "CCO", _make_template(),
                              [_make_precursor("CC=O")])
        assert result is None


# =====================================================================
#  Singleton tests
# =====================================================================


class TestSingleton:
    """Tests for get_scorer/set_scorer."""

    def teardown_method(self):
        set_scorer(None)

    def test_get_scorer_returns_instance(self):
        """get_scorer should return a DisconnectionScorer."""
        scorer = get_scorer()
        assert isinstance(scorer, DisconnectionScorer)

    def test_set_scorer_overrides(self):
        """set_scorer should override the singleton."""
        custom = DisconnectionScorer.__new__(DisconnectionScorer)
        custom._model = None
        custom._model_path = None
        set_scorer(custom)
        assert get_scorer() is custom


# =====================================================================
#  Integration: scoring_method field on Disconnection
# =====================================================================


class TestScoringMethodField:
    """Tests that scoring_method is populated correctly."""

    def test_disconnection_default_method(self):
        """Disconnection dataclass defaults to 'heuristic'."""
        d = Disconnection(
            template=_make_template(),
            precursors=[_make_precursor("CC=O")],
            score=50.0,
        )
        assert d.scoring_method == "heuristic"

    def test_disconnection_ml_method(self):
        """Disconnection can be set to 'ml'."""
        d = Disconnection(
            template=_make_template(),
            precursors=[_make_precursor("CC=O")],
            score=72.5,
            scoring_method="ml",
        )
        assert d.scoring_method == "ml"

    def test_retro_heuristic_fallback(self):
        """Without ML model, all disconnections should use 'heuristic'."""
        set_scorer(None)
        # Force fresh singleton without model
        scorer = DisconnectionScorer.__new__(DisconnectionScorer)
        scorer._model = None
        scorer._model_path = None
        set_scorer(scorer)

        mol = parse("CC(=O)O")  # acetic acid
        tree = retrosynthesis(mol, max_depth=1, beam_width=3)
        for disc in tree.target.disconnections:
            assert disc.scoring_method == "heuristic"

        set_scorer(None)


# =====================================================================
#  Integration with bundled model (if present)
# =====================================================================


class TestBundledModel:
    """Integration tests that run only if the trained model exists."""

    MODEL_PATH = os.path.join(
        os.path.dirname(__file__), "..", "molbuilder", "data", "retro_scorer.pkl"
    )

    def teardown_method(self):
        set_scorer(None)

    @pytest.mark.skipif(
        not os.path.isfile(
            os.path.join(
                os.path.dirname(__file__),
                "..", "molbuilder", "data", "retro_scorer.pkl",
            )
        ),
        reason="Trained retro scorer model not present",
    )
    def test_bundled_model_loads(self):
        """Bundled model should auto-load via get_scorer."""
        set_scorer(None)
        scorer = get_scorer()
        assert scorer.is_loaded

    @pytest.mark.skipif(
        not os.path.isfile(
            os.path.join(
                os.path.dirname(__file__),
                "..", "molbuilder", "data", "retro_scorer.pkl",
            )
        ),
        reason="Trained retro scorer model not present",
    )
    def test_bundled_score_ethanol(self):
        """Bundled model should score ethanol disconnection."""
        set_scorer(None)
        scorer = get_scorer()
        mol = parse("CCO")
        result = scorer.score(
            mol, "CCO", _make_template(), [_make_precursor("CC=O")]
        )
        assert result is not None
        assert 0.0 <= result <= 100.0
