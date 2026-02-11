"""Shared fixtures for the molbuilder-client test suite."""

from __future__ import annotations

import pytest
import respx

from molbuilder_client import AsyncMolBuilder, MolBuilder
from molbuilder_client._base import DEFAULT_BASE_URL

API_KEY = "mb_test_key_abc123"
BASE = DEFAULT_BASE_URL + "/api/v1"


@pytest.fixture()
def base_url() -> str:
    return BASE


@pytest.fixture()
def mock_api():
    """Yield a started ``respx`` mock router scoped to the API base URL."""
    with respx.mock(base_url=BASE) as router:
        yield router


@pytest.fixture()
def client(mock_api) -> MolBuilder:
    """Sync client wired to the mocked transport."""
    c = MolBuilder(api_key=API_KEY)
    yield c
    c.close()


@pytest.fixture()
def async_client(mock_api) -> AsyncMolBuilder:
    """Async client wired to the mocked transport."""
    return AsyncMolBuilder(api_key=API_KEY)
