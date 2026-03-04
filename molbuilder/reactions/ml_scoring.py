"""ML-based retrosynthetic disconnection scoring.

This module provides the inference interface for ML-based disconnection
scoring backed by a trained scikit-learn gradient boosting regressor.
When a trained model is available at ``molbuilder/data/retro_scorer.pkl``,
the scorer loads it automatically and returns ML-predicted scores.
If the model file is missing or scikit-learn is not installed, the
scorer gracefully falls back (returns None) and callers use heuristics.

Integration point: retrosynthesis._build_retro_node() checks the ML
scorer first and falls back to score_disconnection() if it returns None.
"""

from __future__ import annotations

import hashlib
import logging
import os
from typing import TYPE_CHECKING

from molbuilder.reactions.retro_features import (
    extract_retro_features,
    ALL_RETRO_FEATURE_NAMES,
)

if TYPE_CHECKING:
    from molbuilder.molecule.graph import Molecule
    from molbuilder.reactions.reaction_types import ReactionTemplate
    from molbuilder.reactions.retrosynthesis import Precursor
    from molbuilder.reactions.functional_group_detect import FunctionalGroup

logger = logging.getLogger("molbuilder.ml_scoring")

# Path to the bundled model file (ships in molbuilder/data/)
_MODEL_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "data")
_DEFAULT_MODEL_PATH = os.path.join(_MODEL_DIR, "retro_scorer.pkl")


class DisconnectionScorer:
    """ML-based disconnection scoring.

    When no model is loaded (the default), score() returns None and
    the caller should fall back to the heuristic approach.

    Parameters
    ----------
    model_path : str | None
        Path to a serialized model file. If None, attempts to load
        the bundled model from ``molbuilder/data/retro_scorer.pkl``.
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

        The expected dict keys are: score_model, feature_names, version.

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
                "scikit-learn/joblib not installed; ML scoring unavailable. "
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
        if list(stored_names) != list(ALL_RETRO_FEATURE_NAMES):
            logger.warning(
                "Model feature names mismatch: expected %d features, got %d. "
                "Model may be outdated.",
                len(ALL_RETRO_FEATURE_NAMES), len(stored_names),
            )
            return

        self._model = model_dict
        logger.info(
            "ML retro scorer loaded (version=%s) from %s",
            model_dict.get("version", "unknown"), path,
        )

    @property
    def is_loaded(self) -> bool:
        """Return True if a model is loaded and ready for inference."""
        return self._model is not None

    def score(
        self,
        target_mol: "Molecule",
        target_smiles: str,
        template: "ReactionTemplate",
        precursors: list["Precursor"],
        depth: int = 0,
        target_fgs: list["FunctionalGroup"] | None = None,
    ) -> float | None:
        """Score a disconnection using the ML model.

        Returns None if no model is loaded, signalling the caller to
        fall back to heuristic scoring.

        Parameters
        ----------
        target_mol : Molecule
            The parsed target molecule.
        target_smiles : str
            SMILES string for the target.
        template : ReactionTemplate
            The reaction template applied in reverse.
        precursors : list[Precursor]
            Precursor molecules from this disconnection.
        depth : int
            Depth in the retro tree (root = 0).
        target_fgs : list[FunctionalGroup] | None
            Pre-computed functional groups (avoids recomputation).

        Returns
        -------
        float | None
            Predicted score clamped to [0.0, 100.0], or None.
        """
        if self._model is None:
            return None

        try:
            features = extract_retro_features(
                target_mol, target_smiles, template, precursors,
                depth=depth, target_fgs=target_fgs,
            )
            feature_vector = [features[k] for k in ALL_RETRO_FEATURE_NAMES]

            raw_score = float(
                self._model["score_model"].predict([feature_vector])[0]
            )
            return max(0.0, min(100.0, raw_score))
        except Exception:
            logger.warning(
                "ML scoring failed for %s", target_smiles, exc_info=True
            )
            return None


# Module-level singleton
_scorer: DisconnectionScorer | None = None


def get_scorer() -> DisconnectionScorer:
    """Get or create the module-level DisconnectionScorer singleton."""
    global _scorer
    if _scorer is None:
        _scorer = DisconnectionScorer()
    return _scorer


def set_scorer(scorer: DisconnectionScorer | None) -> None:
    """Override the module-level scorer (useful for testing)."""
    global _scorer
    _scorer = scorer
