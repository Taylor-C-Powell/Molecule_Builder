"""Tests for HTTP error → exception mapping across all status codes."""

from __future__ import annotations

import httpx
import pytest
import respx

from molbuilder_client import (
    AuthenticationError,
    ForbiddenError,
    MolBuilder,
    MolBuilderError,
    NotFoundError,
    RateLimitError,
    ServerError,
    ServiceUnavailableError,
    ValidationError,
)


def test_401_raises_authentication_error(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/elements/X").mock(
        return_value=httpx.Response(401, json={"error": "Invalid API key"})
    )

    with pytest.raises(AuthenticationError) as exc_info:
        client.get_element("X")

    assert exc_info.value.status_code == 401
    assert "Invalid API key" in exc_info.value.message


def test_403_raises_forbidden_error(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/elements/X").mock(
        return_value=httpx.Response(403, json={"error": "Admin only"})
    )

    with pytest.raises(ForbiddenError) as exc_info:
        client.get_element("X")

    assert exc_info.value.status_code == 403


def test_404_raises_not_found_error(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/elements/Xx").mock(
        return_value=httpx.Response(404, json={"error": "Element not found"})
    )

    with pytest.raises(NotFoundError) as exc_info:
        client.get_element("Xx")

    assert exc_info.value.status_code == 404


def test_422_raises_validation_error(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.post("/molecule/from-smiles").mock(
        return_value=httpx.Response(422, json={"error": "Invalid SMILES"})
    )

    with pytest.raises(ValidationError) as exc_info:
        client.from_smiles("not_a_smiles!!!")

    assert exc_info.value.status_code == 422


def test_429_raises_rate_limit_error_with_retry_after(
    client: MolBuilder, mock_api: respx.Router
) -> None:
    mock_api.get("/elements/C").mock(
        return_value=httpx.Response(
            429,
            json={"error": "Rate limit exceeded"},
            headers={"Retry-After": "60"},
        )
    )

    with pytest.raises(RateLimitError) as exc_info:
        client.get_element("C")

    assert exc_info.value.status_code == 429
    assert exc_info.value.retry_after == 60


def test_429_without_retry_after_header(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/elements/C").mock(
        return_value=httpx.Response(429, json={"error": "Too many requests"})
    )

    with pytest.raises(RateLimitError) as exc_info:
        client.get_element("C")

    assert exc_info.value.retry_after is None


def test_500_raises_server_error(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/elements/C").mock(
        return_value=httpx.Response(500, json={"error": "Internal error"})
    )

    with pytest.raises(ServerError) as exc_info:
        client.get_element("C")

    assert exc_info.value.status_code == 500


def test_501_raises_service_unavailable(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.post("/billing/checkout").mock(
        return_value=httpx.Response(501, json={"error": "Billing not configured"})
    )

    with pytest.raises(ServiceUnavailableError):
        client.create_checkout("pro_monthly")


def test_503_raises_service_unavailable(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/elements/C").mock(
        return_value=httpx.Response(503, text="Service Unavailable")
    )

    with pytest.raises(ServiceUnavailableError):
        client.get_element("C")


def test_unknown_4xx_raises_base_error(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/elements/C").mock(
        return_value=httpx.Response(418, text="I'm a teapot")
    )

    with pytest.raises(MolBuilderError) as exc_info:
        client.get_element("C")

    assert exc_info.value.status_code == 418


def test_non_json_error_body(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/elements/C").mock(
        return_value=httpx.Response(500, text="<html>Server Error</html>")
    )

    with pytest.raises(ServerError) as exc_info:
        client.get_element("C")

    assert "<html>" in exc_info.value.message


def test_exception_hierarchy() -> None:
    """All specific exceptions inherit from MolBuilderError."""
    assert issubclass(AuthenticationError, MolBuilderError)
    assert issubclass(ForbiddenError, MolBuilderError)
    assert issubclass(NotFoundError, MolBuilderError)
    assert issubclass(ValidationError, MolBuilderError)
    assert issubclass(RateLimitError, MolBuilderError)
    assert issubclass(ServerError, MolBuilderError)
    assert issubclass(ServiceUnavailableError, MolBuilderError)


def test_catch_all_with_base_class(client: MolBuilder, mock_api: respx.Router) -> None:
    """Users can catch MolBuilderError to handle any API failure."""
    mock_api.get("/elements/C").mock(
        return_value=httpx.Response(401, json={"error": "bad key"})
    )

    with pytest.raises(MolBuilderError):
        client.get_element("C")
