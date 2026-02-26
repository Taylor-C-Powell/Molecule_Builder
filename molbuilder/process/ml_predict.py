"""ML-based reaction condition prediction (stub).

This module provides the inference interface for ML-based condition
prediction. Currently it is a stub that returns None (no model loaded),
causing callers to fall back to the existing heuristic approach.

When a trained model is available (e.g. sklearn, torch), it can be
loaded via ConditionPredictor(model_path="path/to/model.pkl") and the
predict() method will return ReactionConditions instead of None.

Integration point: condition_prediction.py checks the ML predictor
first and falls back to heuristics if it returns None.
"""

from __future__ import annotations

import logging

from molbuilder.process.conditions import ReactionConditions
from molbuilder.process.ml_features import extract_features

logger = logging.getLogger("molbuilder.ml_predict")


class ConditionPredictor:
    """ML-based condition prediction.

    When no model is loaded (the default), predict() returns None and
    the caller should fall back to the heuristic approach.

    Parameters
    ----------
    model_path : str | None
        Path to a serialized model file. If None, no model is loaded
        and predict() always returns None.
    """

    def __init__(self, model_path: str | None = None):
        self._model = None
        self._model_path = model_path
        if model_path is not None:
            self._load_model(model_path)

    def _load_model(self, path: str) -> None:
        """Attempt to load a serialized model.

        Currently a no-op stub. When a real model format is chosen
        (e.g. joblib for sklearn, torch.load for PyTorch), this method
        will deserialize the model into self._model.
        """
        logger.info("ML model loading not yet implemented (path: %s)", path)
        # Future implementation:
        # import joblib
        # self._model = joblib.load(path)

    @property
    def is_loaded(self) -> bool:
        """Return True if a model is loaded and ready for inference."""
        return self._model is not None

    def predict(
        self,
        smiles: str,
        template_name: str | None = None,
        scale_kg: float = 1.0,
    ) -> ReactionConditions | None:
        """Predict reaction conditions using the ML model.

        Returns None if no model is loaded, signalling the caller to
        fall back to heuristic prediction.

        Parameters
        ----------
        smiles : str
            Substrate SMILES string.
        template_name : str | None
            Optional reaction category name for context.
        scale_kg : float
            Target production scale in kilograms.

        Returns
        -------
        ReactionConditions | None
            Predicted conditions, or None if no model is available.
        """
        if self._model is None:
            return None

        features = extract_features(smiles, template_name, scale_kg)

        # Future implementation:
        # feature_vector = [features[k] for k in ALL_FEATURE_NAMES]
        # prediction = self._model.predict([feature_vector])[0]
        # return ReactionConditions(
        #     temperature_celsius=prediction["temperature"],
        #     solvent=prediction["solvent"],
        #     ...
        # )

        return None  # pragma: no cover


# Module-level singleton (optional, for shared use)
_predictor: ConditionPredictor | None = None


def get_predictor() -> ConditionPredictor:
    """Get or create the module-level ConditionPredictor singleton."""
    global _predictor
    if _predictor is None:
        _predictor = ConditionPredictor()
    return _predictor


def set_predictor(predictor: ConditionPredictor | None) -> None:
    """Override the module-level predictor (useful for testing)."""
    global _predictor
    _predictor = predictor
