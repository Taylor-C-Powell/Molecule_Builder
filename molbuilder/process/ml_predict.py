"""ML-based reaction condition prediction.

This module provides the inference interface for ML-based condition
prediction backed by a trained scikit-learn gradient boosting ensemble.
When a trained model is available at ``molbuilder/data/condition_model.pkl``,
the predictor loads it automatically and returns ML-predicted conditions.
If the model file is missing or scikit-learn is not installed, the
predictor gracefully falls back (returns None) and callers use heuristics.

Integration point: condition_prediction.py checks the ML predictor
first and falls back to heuristics if it returns None.
"""

from __future__ import annotations

import hashlib
import logging
import os

from molbuilder.process.conditions import ReactionConditions, addition_rate_for_scale
from molbuilder.process.ml_features import extract_features, ALL_FEATURE_NAMES

logger = logging.getLogger("molbuilder.ml_predict")

# Path to the bundled model file (ships in molbuilder/data/)
_MODEL_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "data")
_DEFAULT_MODEL_PATH = os.path.join(_MODEL_DIR, "condition_model.pkl")

# Per-target confidence gates.  Targets with known poor performance
# (R-squared < 0) are disabled; predictions fall back to heuristics.
# This is safer than returning ML predictions that are worse than the mean.
_DISABLED_TARGETS: frozenset[str] = frozenset({"temperature"})


class ConditionPredictor:
    """ML-based condition prediction.

    When no model is loaded (the default), predict() returns None and
    the caller should fall back to the heuristic approach.

    Parameters
    ----------
    model_path : str | None
        Path to a serialized model file. If None, attempts to load
        the bundled model from ``molbuilder/data/condition_model.pkl``.
    """

    def __init__(self, model_path: str | None = None):
        self._model: dict | None = None
        self._model_path = model_path
        path = model_path if model_path is not None else _DEFAULT_MODEL_PATH
        if os.path.isfile(path):
            self._load_model(path)

    @staticmethod
    def _verify_checksum(path: str) -> bool:
        """Verify SHA-256 checksum of model file against sidecar .sha256 file.

        Returns True if the checksum matches or no sidecar file exists
        (allowing custom models without checksums). Returns False if the
        sidecar exists but the hash does not match (possible tampering).
        """
        sha_path = path + ".sha256"
        if not os.path.isfile(sha_path):
            logger.debug("No checksum file at %s; skipping verification", sha_path)
            return True

        try:
            with open(sha_path, "r") as f:
                expected = f.read().strip().lower()
        except OSError:
            logger.warning("Could not read checksum file %s", sha_path)
            return False

        sha = hashlib.sha256()
        try:
            with open(path, "rb") as f:
                for chunk in iter(lambda: f.read(65536), b""):
                    sha.update(chunk)
        except OSError:
            logger.warning("Could not read model file %s for checksum", path)
            return False

        actual = sha.hexdigest().lower()
        if actual != expected:
            logger.error(
                "Model checksum MISMATCH for %s: expected %s, got %s. "
                "The model file may have been tampered with.",
                path, expected, actual,
            )
            return False
        return True

    def _load_model(self, path: str) -> None:
        """Load a serialized model dict from a joblib pickle file.

        The expected dict keys are: temperature_model, solvent_model,
        catalyst_model, yield_model, solvent_classes, catalyst_classes,
        feature_names, version.

        If scikit-learn/joblib is not installed, logs a warning and
        leaves ``self._model`` as None.
        """
        # Verify integrity before deserializing
        if not self._verify_checksum(path):
            logger.error("Refusing to load model with failed checksum: %s", path)
            return

        try:
            import joblib
        except ImportError:
            logger.warning(
                "scikit-learn/joblib not installed; ML prediction unavailable. "
                "Install with: pip install molbuilder[ml]"
            )
            return

        try:
            model_dict = joblib.load(path)
        except Exception:
            logger.warning("Failed to load ML model from %s", path, exc_info=True)
            return

        # Validate feature names match
        stored_names = model_dict.get("feature_names", [])
        if list(stored_names) != list(ALL_FEATURE_NAMES):
            logger.warning(
                "Model feature names mismatch: expected %d features, got %d. "
                "Model may be outdated.",
                len(ALL_FEATURE_NAMES), len(stored_names),
            )
            return

        self._model = model_dict
        logger.info(
            "ML condition model loaded (version=%s) from %s",
            model_dict.get("version", "unknown"), path,
        )

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
        feature_vector = [features[k] for k in ALL_FEATURE_NAMES]

        try:
            # Temperature model is gated: skip if in _DISABLED_TARGETS
            # (current model has negative R-squared for temperature).
            # Use a safe heuristic default of 25 C (room temperature).
            if "temperature" not in _DISABLED_TARGETS:
                temp_pred = float(
                    self._model["temperature_model"].predict([feature_vector])[0]
                )
                temp_pred = max(-80.0, min(300.0, temp_pred))
            else:
                temp_pred = 25.0

            solvent_idx = int(
                self._model["solvent_model"].predict([feature_vector])[0]
            )
            solvent_classes = self._model["solvent_classes"]
            solvent = (
                solvent_classes[solvent_idx]
                if 0 <= solvent_idx < len(solvent_classes)
                else "THF"
            )

            catalyst_idx = int(
                self._model["catalyst_model"].predict([feature_vector])[0]
            )
            catalyst_classes = self._model["catalyst_classes"]
            catalyst = (
                catalyst_classes[catalyst_idx]
                if 0 <= catalyst_idx < len(catalyst_classes)
                else ""
            )

            yield_pred = float(
                self._model["yield_model"].predict([feature_vector])[0]
            )
            yield_pred = max(0.0, min(100.0, yield_pred))
        except Exception:
            logger.warning("ML prediction failed for %s", smiles, exc_info=True)
            return None

        addition_rate = addition_rate_for_scale(scale_kg)

        # Scale-dependent reaction time
        import math as _math
        base_hours = 0.5 + (yield_pred / 100.0) * 2.0
        if temp_pred < -40:
            base_hours *= 1.5
            base_hours += 1.0
        elif temp_pred > 120:
            base_hours *= 0.7
        if scale_kg > 1.0:
            decades = _math.log10(scale_kg)
            base_hours *= (1.0 + 0.2 * decades)
        reaction_time = round(max(0.5, base_hours), 1)

        # Build notes string
        notes_parts = []
        if catalyst and catalyst != "none":
            notes_parts.append(f"Catalyst: {catalyst}")
        notes_parts.append(f"Predicted yield: {yield_pred:.0f}%")
        notes = "  ".join(notes_parts) if notes_parts else "No special notes."

        # Indicate which targets came from ML vs heuristic
        if _DISABLED_TARGETS:
            source = "ML model (temperature=heuristic)"
        else:
            source = "ML model"

        return ReactionConditions(
            temperature_C=round(temp_pred, 1),
            pressure_atm=1.0,
            solvent=solvent,
            concentration_M=0.3,
            addition_rate=addition_rate,
            reaction_time_hours=reaction_time,
            atmosphere="N2",
            workup_procedure="Standard aqueous workup: dilute, extract, wash, dry, concentrate.",
            notes=notes,
            data_source=source,
        )


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
