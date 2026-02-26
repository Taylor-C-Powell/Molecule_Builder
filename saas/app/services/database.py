"""Database backend abstraction layer.

Provides a pluggable interface for database operations so that the
application can be migrated from SQLite to PostgreSQL (or other
backends) without changing service-layer code.

Current status: skeleton with SQLiteBackend implementation.
The existing DB classes (user_db, audit_db, job_db, library_db)
continue using direct SQLite connections. This abstraction is ready
for Tier 3 when PostgreSQL migration is actually needed.

Usage::

    from app.services.database import get_backend

    db = get_backend()
    rows = db.execute("SELECT * FROM users WHERE email = ?", (email,))
"""

from __future__ import annotations

import sqlite3
import threading
from abc import ABC, abstractmethod


class DatabaseBackend(ABC):
    """Abstract interface for pluggable database backends.

    All concrete backends must implement execute() and execute_insert().
    Results are returned as lists of dicts for uniform access regardless
    of the underlying driver.
    """

    @abstractmethod
    def execute(self, sql: str, params: tuple = ()) -> list[dict]:
        """Execute a query and return all rows as dicts."""
        ...

    @abstractmethod
    def execute_insert(self, sql: str, params: tuple = ()) -> int:
        """Execute an INSERT and return the last inserted row ID."""
        ...

    @abstractmethod
    def execute_update(self, sql: str, params: tuple = ()) -> int:
        """Execute an UPDATE/DELETE and return the number of affected rows."""
        ...

    @abstractmethod
    def close(self) -> None:
        """Close the database connection."""
        ...


class SQLiteBackend(DatabaseBackend):
    """SQLite implementation of the database backend.

    Uses thread-local connections with WAL mode for concurrent reads.
    This mirrors the pattern used by the existing service DBs.
    """

    def __init__(self, db_path: str):
        self._db_path = db_path
        self._local = threading.local()

    def _get_conn(self) -> sqlite3.Connection:
        if not hasattr(self._local, "conn") or self._local.conn is None:
            conn = sqlite3.connect(self._db_path)
            conn.execute("PRAGMA journal_mode=WAL")
            conn.execute("PRAGMA foreign_keys=ON")
            conn.row_factory = sqlite3.Row
            self._local.conn = conn
        return self._local.conn

    def execute(self, sql: str, params: tuple = ()) -> list[dict]:
        conn = self._get_conn()
        rows = conn.execute(sql, params).fetchall()
        return [dict(r) for r in rows]

    def execute_insert(self, sql: str, params: tuple = ()) -> int:
        conn = self._get_conn()
        cursor = conn.execute(sql, params)
        conn.commit()
        return cursor.lastrowid

    def execute_update(self, sql: str, params: tuple = ()) -> int:
        conn = self._get_conn()
        cursor = conn.execute(sql, params)
        conn.commit()
        return cursor.rowcount

    def close(self) -> None:
        if hasattr(self._local, "conn") and self._local.conn is not None:
            self._local.conn.close()
            self._local.conn = None


# Future: PostgreSQL backend
#
# class PostgresBackend(DatabaseBackend):
#     """PostgreSQL implementation using asyncpg or psycopg2.
#
#     Requirements: pip install psycopg2-binary (or asyncpg for async)
#     """
#
#     def __init__(self, dsn: str):
#         self._dsn = dsn
#         self._conn = None
#
#     def execute(self, sql: str, params: tuple = ()) -> list[dict]:
#         # Convert ? placeholders to %s for PostgreSQL
#         pg_sql = sql.replace("?", "%s")
#         ...
#
#     def execute_insert(self, sql: str, params: tuple = ()) -> int:
#         pg_sql = sql.replace("?", "%s") + " RETURNING id"
#         ...
#
#     def execute_update(self, sql: str, params: tuple = ()) -> int:
#         pg_sql = sql.replace("?", "%s")
#         ...
#
#     def close(self) -> None:
#         if self._conn:
#             self._conn.close()


# Module-level singleton
_backend: DatabaseBackend | None = None


def get_backend() -> DatabaseBackend:
    """Get or create the database backend based on configuration."""
    global _backend
    if _backend is None:
        from app.config import settings
        backend_type = getattr(settings, "database_backend", "sqlite")
        if backend_type == "sqlite":
            _backend = SQLiteBackend(settings.user_db_path)
        # Future: elif backend_type == "postgresql":
        #     _backend = PostgresBackend(settings.database_url)
        else:
            raise ValueError(f"Unknown database backend: {backend_type}")
    return _backend


def set_backend(backend: DatabaseBackend | None) -> None:
    """Override the database backend (useful for testing)."""
    global _backend
    _backend = backend
