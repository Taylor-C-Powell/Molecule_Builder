"""Bundled data files for MolBuilder.

Provides :func:`load_ord_conditions` which lazy-loads the pre-computed
Open Reaction Database condition statistics shipped with the package.
"""

from __future__ import annotations

import json
from importlib import resources
from typing import Any

_ORD_CACHE: dict[str, Any] | None = None


def load_ord_conditions() -> dict[str, Any]:
    """Load and cache the pre-computed ORD condition statistics.

    Returns a dict with keys ``"_meta"`` and ``"reactions"``.
    If the bundled JSON file is missing or empty, returns a dict with
    an empty ``"reactions"`` mapping so callers always get a valid structure.
    """
    global _ORD_CACHE
    if _ORD_CACHE is not None:
        return _ORD_CACHE

    try:
        ref = resources.files("molbuilder.data").joinpath("ord_conditions.json")
        text = ref.read_text(encoding="utf-8")
        data = json.loads(text)
        if not isinstance(data, dict) or "reactions" not in data:
            data = {"_meta": {}, "reactions": {}}
    except Exception:
        data = {"_meta": {}, "reactions": {}}

    _ORD_CACHE = data
    return _ORD_CACHE


def _clear_ord_cache() -> None:
    """Reset the cached ORD data (used by tests)."""
    global _ORD_CACHE
    _ORD_CACHE = None
