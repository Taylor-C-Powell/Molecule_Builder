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


def _using_postgres() -> bool:
    """True when tests should run against PostgreSQL."""
    return os.environ.get("DATABASE_BACKEND", "sqlite").lower() == "postgresql"


def _pg_truncate(table: str) -> None:
    """Truncate a table using a direct psycopg connection (not the pool).

    Uses a short-lived connection outside the pool so truncation
    cannot be blocked by or interfere with pooled connections.
    """
    import psycopg
    dsn = os.environ.get("DATABASE_URL", "")
    with psycopg.connect(dsn) as conn:
        conn.execute(f"DELETE FROM {table}")
        conn.commit()


# ------------------------------------------------------------------ #
# Session-scoped PG pool (created once, shared by ALL tests)
# ------------------------------------------------------------------ #

@pytest.fixture(scope="session")
def _pg_backend():
    """Create a single PG backend for the entire test session.

    The pool's close() method is monkey-patched to a no-op so that
    nothing (lifespan shutdown, DB .close() calls, GC __del__) can
    accidentally close the pool while tests are running.
    """
    if not _using_postgres():
        yield None
        return
    from app.services.database import PostgresBackend
    dsn = os.environ.get("DATABASE_URL", "")
    backend = PostgresBackend(dsn=dsn, min_size=2, max_size=5)
    # Prevent accidental closure during the test session
    real_close = backend.close
    backend.close = lambda: None
    yield backend
    real_close()


# ------------------------------------------------------------------ #
# Per-test backend singleton management
# ------------------------------------------------------------------ #

@pytest.fixture(autouse=True)
def _ensure_pg_backend(_pg_backend):
    """Pin the database backend singleton to the session pool.

    Runs before all other autouse fixtures (defined first).  On
    teardown it re-pins the singleton so that the next test starts
    with the pool available even if something cleared it mid-test.
    """
    if not _using_postgres():
        yield
        return
    from app.services.database import set_backend
    set_backend(_pg_backend)
    yield
    # Re-pin after test in case lifespan or other code cleared it
    set_backend(_pg_backend)


# ------------------------------------------------------------------ #
# Per-test state reset
# ------------------------------------------------------------------ #

@pytest.fixture(autouse=True)
def _reset_rate_limiter():
    """Clear rate limiter state between tests."""
    rate_limiter._windows.clear()
    rate_limiter._expensive_windows.clear()
    rate_limiter._ip_windows.clear()


@pytest.fixture(autouse=True)
def _temp_audit_db(tmp_path, _ensure_pg_backend):
    """Use a temporary audit DB for each test."""
    from app.services.audit_db import AuditDB, set_audit_db
    if _using_postgres():
        _pg_truncate("audit_log")
        from app.services.database import get_backend
        db = AuditDB(backend=get_backend())
    else:
        db_path = str(tmp_path / "test_audit.db")
        db = AuditDB(db_path)
    set_audit_db(db)
    yield db
    if not _using_postgres():
        db.close()
    set_audit_db(None)


@pytest.fixture(autouse=True)
def _temp_user_db(tmp_path, _ensure_pg_backend):
    """Use a temporary user DB for each test."""
    from app.services.user_db import UserDB, set_user_db
    if _using_postgres():
        _pg_truncate("api_keys")
        from app.services.database import get_backend
        db = UserDB(backend=get_backend())
    else:
        db_path = str(tmp_path / "test_users.db")
        db = UserDB(db_path)
    set_user_db(db)
    yield db
    if not _using_postgres():
        db.close()
    set_user_db(None)


@pytest.fixture(autouse=True)
def _temp_job_db(tmp_path, _ensure_pg_backend):
    """Use a temporary job DB for each test."""
    from app.services.job_db import JobDB, set_job_db
    if _using_postgres():
        _pg_truncate("jobs")
        from app.services.database import get_backend
        db = JobDB(backend=get_backend())
    else:
        db_path = str(tmp_path / "test_jobs.db")
        db = JobDB(db_path)
    set_job_db(db)
    yield db
    if not _using_postgres():
        db.close()
    set_job_db(None)


@pytest.fixture(autouse=True)
def _temp_library_db(tmp_path, _ensure_pg_backend):
    """Use a temporary library DB for each test."""
    from app.services.library_db import LibraryDB, set_library_db
    if _using_postgres():
        _pg_truncate("library_molecules")
        from app.services.database import get_backend
        db = LibraryDB(backend=get_backend())
    else:
        db_path = str(tmp_path / "test_library.db")
        db = LibraryDB(db_path)
    set_library_db(db)
    yield db
    if not _using_postgres():
        db.close()
    set_library_db(None)


@pytest.fixture(autouse=True)
def _reset_molecule_store():
    """Reset molecule store singleton between tests."""
    from app.services.molecule_store import set_molecule_store
    yield
    set_molecule_store(None)


@pytest.fixture(autouse=True)
def _reset_batch_worker():
    """Reset batch worker between tests to prevent thread leaks."""
    from app.services.batch_worker import get_batch_worker, set_batch_worker
    yield
    try:
        worker = get_batch_worker()
        worker.shutdown()
    except Exception:
        pass
    set_batch_worker(None)


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
