"""Tests for the async client -- mirrors key sync tests with async/await."""

from __future__ import annotations

import httpx
import pytest
import respx

from molbuilder_client import (
    AsyncMolBuilder,
    AuthenticationError,
    Element,
    MoleculeInfo,
    RateLimitError,
    RetrosynthesisPlan,
    Token,
)


async def test_async_register(async_client: AsyncMolBuilder, mock_api: respx.Router) -> None:
    mock_api.post("/auth/register").mock(
        return_value=httpx.Response(
            200,
            json={"api_key": "mb_a", "email": "a@b.com", "tier": "free", "role": "chemist"},
        )
    )

    result = await async_client.register("a@b.com")
    assert result.api_key == "mb_a"


async def test_async_get_token(async_client: AsyncMolBuilder, mock_api: respx.Router) -> None:
    mock_api.post("/auth/token").mock(
        return_value=httpx.Response(
            200,
            json={"access_token": "jwt_x", "token_type": "bearer", "expires_in": 3600},
        )
    )

    result = await async_client.get_token()
    assert isinstance(result, Token)
    assert result.access_token == "jwt_x"


async def test_async_from_smiles(async_client: AsyncMolBuilder, mock_api: respx.Router) -> None:
    mock_api.post("/molecule/from-smiles").mock(
        return_value=httpx.Response(
            200,
            json={"id": "m1", "name": "", "smiles": "C", "num_atoms": 5, "num_bonds": 4},
        )
    )

    result = await async_client.from_smiles("C")
    assert isinstance(result, MoleculeInfo)
    assert result.id == "m1"


async def test_async_get_element(async_client: AsyncMolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/elements/O").mock(
        return_value=httpx.Response(
            200,
            json={"atomic_number": 8, "symbol": "O", "name": "Oxygen", "atomic_weight": 15.999},
        )
    )

    result = await async_client.get_element("O")
    assert isinstance(result, Element)
    assert result.name == "Oxygen"


async def test_async_retrosynthesis(
    async_client: AsyncMolBuilder, mock_api: respx.Router
) -> None:
    mock_api.post("/retrosynthesis/plan").mock(
        return_value=httpx.Response(
            200,
            json={
                "tree": {
                    "smiles": "CCO",
                    "is_purchasable": True,
                    "functional_groups": [],
                    "best_disconnection": None,
                    "children": [],
                },
                "routes_found": 0,
                "max_depth": 5,
                "beam_width": 5,
                "best_route": None,
            },
        )
    )

    result = await async_client.retrosynthesis("CCO")
    assert isinstance(result, RetrosynthesisPlan)


async def test_async_error_mapping(
    async_client: AsyncMolBuilder, mock_api: respx.Router
) -> None:
    mock_api.get("/elements/X").mock(
        return_value=httpx.Response(401, json={"error": "Unauthorized"})
    )

    with pytest.raises(AuthenticationError):
        await async_client.get_element("X")


async def test_async_rate_limit(async_client: AsyncMolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/elements/C").mock(
        return_value=httpx.Response(
            429,
            json={"error": "Rate limit"},
            headers={"Retry-After": "30"},
        )
    )

    with pytest.raises(RateLimitError) as exc_info:
        await async_client.get_element("C")

    assert exc_info.value.retry_after == 30


async def test_async_context_manager(mock_api: respx.Router) -> None:
    mock_api.get("/elements/H").mock(
        return_value=httpx.Response(
            200,
            json={"atomic_number": 1, "symbol": "H", "name": "Hydrogen", "atomic_weight": 1.008},
        )
    )

    async with AsyncMolBuilder(api_key="mb_test") as client:
        result = await client.get_element("H")
        assert result.symbol == "H"
