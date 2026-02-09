"""Shared test fixtures."""

import os
import uuid
import tempfile
import pytest
from fastapi.testclient import TestClient
from app.main import app
from app.auth.api_keys import api_key_store
from app.auth.rate_limiter import rate_limiter
from app.auth.roles import Role
from app.config import Tier


@pytest.fixture(autouse=True)
def _reset_rate_limiter():
    """Clear rate limiter state between tests."""
    rate_limiter._windows.clear()
    rate_limiter._expensive_windows.clear()


@pytest.fixture(autouse=True)
def _temp_audit_db(tmp_path):
    """Use a temporary audit DB for each test."""
    from app.services.audit_db import AuditDB, set_audit_db
    db_path = str(tmp_path / "test_audit.db")
    db = AuditDB(db_path)
    set_audit_db(db)
    yield db
    db.close()
    set_audit_db(None)


@pytest.fixture(autouse=True)
def _temp_user_db(tmp_path):
    """Use a temporary user DB for each test."""
    from app.services.user_db import UserDB, set_user_db
    db_path = str(tmp_path / "test_users.db")
    db = UserDB(db_path)
    set_user_db(db)
    yield db
    db.close()
    set_user_db(None)


@pytest.fixture(autouse=True)
def _reset_api_key_store():
    """Reset API key store between tests."""
    api_key_store._keys.clear()
    api_key_store._db = None
    yield
    api_key_store._keys.clear()
    api_key_store._db = None


@pytest.fixture()
def client():
    with TestClient(app) as c:
        yield c


@pytest.fixture()
def api_key():
    """Create a free-tier API key with a unique email."""
    email = f"test-{uuid.uuid4().hex[:8]}@example.com"
    return api_key_store.create(email=email, tier=Tier.FREE)


@pytest.fixture()
def auth_headers(api_key):
    return {"X-API-Key": api_key}


@pytest.fixture()
def pro_api_key():
    email = f"pro-{uuid.uuid4().hex[:8]}@example.com"
    return api_key_store.create(email=email, tier=Tier.PRO)


@pytest.fixture()
def pro_headers(pro_api_key):
    return {"X-API-Key": pro_api_key}


@pytest.fixture()
def admin_api_key():
    email = f"admin-{uuid.uuid4().hex[:8]}@example.com"
    return api_key_store.create(email=email, tier=Tier.PRO, role=Role.ADMIN)


@pytest.fixture()
def admin_headers(admin_api_key):
    return {"X-API-Key": admin_api_key}


@pytest.fixture()
def viewer_api_key():
    email = f"viewer-{uuid.uuid4().hex[:8]}@example.com"
    return api_key_store.create(email=email, tier=Tier.FREE, role=Role.VIEWER)


@pytest.fixture()
def viewer_headers(viewer_api_key):
    return {"X-API-Key": viewer_api_key}


@pytest.fixture()
def enterprise_api_key():
    email = f"ent-{uuid.uuid4().hex[:8]}@example.com"
    return api_key_store.create(email=email, tier=Tier.ENTERPRISE, role=Role.CHEMIST)


@pytest.fixture()
def enterprise_headers(enterprise_api_key):
    return {"X-API-Key": enterprise_api_key}
