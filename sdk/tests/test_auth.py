"""Tests for auth endpoints: register and get_token."""

from __future__ import annotations

import httpx
import respx

from molbuilder_client import APIKeyInfo, MolBuilder, Token


def test_register(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.post("/auth/register").mock(
        return_value=httpx.Response(
            200,
            json={
                "api_key": "mb_new_key",
                "email": "user@example.com",
                "tier": "free",
                "role": "chemist",
            },
        )
    )

    result = client.register("user@example.com")

    assert isinstance(result, APIKeyInfo)
    assert result.api_key == "mb_new_key"
    assert result.email == "user@example.com"
    assert result.tier == "free"
    assert result.role == "chemist"


def test_register_sends_correct_body(client: MolBuilder, mock_api: respx.Router) -> None:
    route = mock_api.post("/auth/register").mock(
        return_value=httpx.Response(
            200,
            json={"api_key": "mb_k", "email": "a@b.com", "tier": "free", "role": "chemist"},
        )
    )

    client.register("a@b.com")

    assert route.called
    request = route.calls.last.request
    assert request.headers["X-API-Key"] == "mb_test_key_abc123"


def test_get_token(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.post("/auth/token").mock(
        return_value=httpx.Response(
            200,
            json={
                "access_token": "jwt_abc",
                "token_type": "bearer",
                "expires_in": 3600,
            },
        )
    )

    result = client.get_token()

    assert isinstance(result, Token)
    assert result.access_token == "jwt_abc"
    assert result.token_type == "bearer"
    assert result.expires_in == 3600
