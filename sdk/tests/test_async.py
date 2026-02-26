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


async def test_async_library_save(
    async_client: AsyncMolBuilder, mock_api: respx.Router
) -> None:
    mock_api.post("/library/").mock(
        return_value=httpx.Response(200, json={
            "id": 1, "smiles": "CCO", "name": "ethanol", "tags": [],
            "notes": None, "properties": {}, "created_at": "2026-01-01T00:00:00",
            "updated_at": "2026-01-01T00:00:00",
        })
    )
    from molbuilder_client import LibraryMolecule
    result = await async_client.library_save("CCO", name="ethanol")
    assert isinstance(result, LibraryMolecule)
    assert result.id == 1


async def test_async_batch_submit(
    async_client: AsyncMolBuilder, mock_api: respx.Router
) -> None:
    mock_api.post("/batch/submit").mock(
        return_value=httpx.Response(200, json={
            "job_id": "job_x", "status": "pending", "created_at": "2026-01-01T00:00:00",
        })
    )
    from molbuilder_client import BatchSubmit
    result = await async_client.batch_submit(["CCO"], "properties")
    assert isinstance(result, BatchSubmit)
    assert result.job_id == "job_x"


async def test_async_library_list(
    async_client: AsyncMolBuilder, mock_api: respx.Router
) -> None:
    mock_api.get("/library/").mock(
        return_value=httpx.Response(200, json={
            "molecules": [], "total": 0, "page": 1, "per_page": 20,
        })
    )
    from molbuilder_client import LibraryList
    result = await async_client.library_list()
    assert isinstance(result, LibraryList)
    assert result.total == 0
