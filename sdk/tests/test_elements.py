"""Tests for the elements endpoint."""

from __future__ import annotations

import httpx
import respx

from molbuilder_client import Element, MolBuilder


def test_get_element(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/elements/Fe").mock(
        return_value=httpx.Response(
            200,
            json={
                "atomic_number": 26,
                "symbol": "Fe",
                "name": "Iron",
                "atomic_weight": 55.845,
            },
        )
    )

    result = client.get_element("Fe")

    assert isinstance(result, Element)
    assert result.atomic_number == 26
    assert result.symbol == "Fe"
    assert result.name == "Iron"
    assert result.atomic_weight == 55.845


def test_get_element_forward_compat(client: MolBuilder, mock_api: respx.Router) -> None:
    """Unknown fields from future API versions are silently ignored."""
    mock_api.get("/elements/C").mock(
        return_value=httpx.Response(
            200,
            json={
                "atomic_number": 6,
                "symbol": "C",
                "name": "Carbon",
                "atomic_weight": 12.011,
                "electronegativity": 2.55,  # future field
            },
        )
    )

    result = client.get_element("C")
    assert result.symbol == "C"
    assert not hasattr(result, "electronegativity")
