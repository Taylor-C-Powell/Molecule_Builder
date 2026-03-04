"""Shared utilities for ML model loading and verification.

Provides checksum verification logic used by both the condition
predictor (ml_predict.py) and the disconnection scorer (ml_scoring.py).
"""

from __future__ import annotations

import hashlib
import logging
import os

logger = logging.getLogger("molbuilder.model_utils")


def verify_model_checksum(path: str, require_sidecar: bool = False) -> bool:
    """Verify SHA-256 checksum of model file against sidecar .sha256 file.

    Parameters
    ----------
    path : str
        Path to the model file.
    require_sidecar : bool
        If True, return False when the sidecar is missing (use for
        bundled models where the sidecar should always be present).
        If False, skip verification when sidecar is absent (use for
        user-provided custom model paths).

    Returns
    -------
    bool
        True if checksum matches or verification was skipped.
        False if checksum mismatches or required sidecar is missing.
    """
    sha_path = path + ".sha256"
    if not os.path.isfile(sha_path):
        if require_sidecar:
            logger.error("Required checksum file missing: %s", sha_path)
            return False
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
