"""Database backend abstraction layer.

Provides a pluggable interface for database operations so that the
application can run against SQLite (dev/tests) or PostgreSQL (production)
without changing service-layer SQL.

Usage::

    from app.services.database import get_backend

    db = get_backend()
    rows = db.execute("SELECT * FROM users WHERE email = ?", (email,))

The ``?`` placeholder style used throughout the codebase is automatically
translated to ``%s`` for PostgreSQL at execution time.
"""

from __future__ import annotations

import logging
import sqlite3
import threading
from abc import ABC, abstractmethod

logger = logging.getLogger("molbuilder.database")


class DatabaseIntegrityError(Exception):
    """Raised on unique-constraint or foreign-key violations.

    Service code should catch this instead of driver-specific errors
    (sqlite3.IntegrityError, psycopg.errors.UniqueViolation, etc.).
    """


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

    @property
    def is_postgres(self) -> bool:
        """True when the backend is PostgreSQL."""
        return False

    @property
    def supports_native_json(self) -> bool:
        """True when the backend natively handles dict <-> JSONB."""
        return False


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
        try:
            rows = conn.execute(sql, params).fetchall()
        except sqlite3.IntegrityError as exc:
            raise DatabaseIntegrityError(str(exc)) from exc
        return [dict(r) for r in rows]

    def execute_insert(self, sql: str, params: tuple = ()) -> int:
        conn = self._get_conn()
        try:
            cursor = conn.execute(sql, params)
        except sqlite3.IntegrityError as exc:
            raise DatabaseIntegrityError(str(exc)) from exc
        conn.commit()
        return cursor.lastrowid

    def execute_update(self, sql: str, params: tuple = ()) -> int:
        conn = self._get_conn()
        try:
            cursor = conn.execute(sql, params)
        except sqlite3.IntegrityError as exc:
            raise DatabaseIntegrityError(str(exc)) from exc
        conn.commit()
        return cursor.rowcount

    def close(self) -> None:
        if hasattr(self._local, "conn") and self._local.conn is not None:
            self._local.conn.close()
            self._local.conn = None


def _translate_sql(sql: str) -> str:
    """Convert ``?`` placeholders to ``%s`` for psycopg."""
    # Simple replacement is safe because none of our SQL uses literal '?'
    # inside strings.  If that ever changes, use a proper parser.
    return sql.replace("?", "%s")


class PostgresBackend(DatabaseBackend):
    """PostgreSQL implementation using psycopg v3 + connection pool.

    Requires ``pip install psycopg[binary,pool]>=3.1``.

    All SQL written with ``?`` placeholders is auto-translated to ``%s``.
    Results are returned as ``list[dict]`` via ``dict_row`` row factory.
    """

    def __init__(self, dsn: str, min_size: int = 2, max_size: int = 10):
        try:
            from psycopg_pool import ConnectionPool
            from psycopg.rows import dict_row
        except ImportError as exc:
            raise ImportError(
                "PostgreSQL backend requires psycopg[binary,pool]>=3.1. "
                "Install with: pip install 'psycopg[binary,pool]>=3.1'"
            ) from exc

        self._pool = ConnectionPool(
            conninfo=dsn,
            min_size=min_size,
            max_size=max_size,
            kwargs={"row_factory": dict_row},
        )
        logger.info(
            "PostgreSQL connection pool created (min=%d, max=%d)", min_size, max_size
        )

    @property
    def is_postgres(self) -> bool:
        return True

    @property
    def supports_native_json(self) -> bool:
        return True

    def execute(self, sql: str, params: tuple = ()) -> list[dict]:
        pg_sql = _translate_sql(sql)
        with self._pool.connection() as conn:
            try:
                rows = conn.execute(pg_sql, params).fetchall()
            except Exception as exc:
                if _is_integrity_error(exc):
                    raise DatabaseIntegrityError(str(exc)) from exc
                raise
        return list(rows)

    def execute_insert(self, sql: str, params: tuple = ()) -> int:
        pg_sql = _translate_sql(sql)
        # Append RETURNING id if the table has an integer PK named 'id'
        # and the query doesn't already contain RETURNING.
        if "RETURNING" not in pg_sql.upper():
            pg_sql = pg_sql.rstrip().rstrip(";") + " RETURNING id"
        with self._pool.connection() as conn:
            try:
                row = conn.execute(pg_sql, params).fetchone()
            except Exception as exc:
                if _is_integrity_error(exc):
                    raise DatabaseIntegrityError(str(exc)) from exc
                raise
            conn.commit()
        return row["id"] if row else 0

    def execute_update(self, sql: str, params: tuple = ()) -> int:
        pg_sql = _translate_sql(sql)
        with self._pool.connection() as conn:
            try:
                cursor = conn.execute(pg_sql, params)
            except Exception as exc:
                if _is_integrity_error(exc):
                    raise DatabaseIntegrityError(str(exc)) from exc
                raise
            conn.commit()
        return cursor.rowcount

    def close(self) -> None:
        self._pool.close()
        logger.info("PostgreSQL connection pool closed")


def _is_integrity_error(exc: Exception) -> bool:
    """Check if an exception is a database integrity violation."""
    # psycopg raises subclasses of psycopg.errors.IntegrityError
    exc_type = type(exc).__name__
    exc_module = type(exc).__module__ or ""
    return (
        "IntegrityError" in exc_type
        or "UniqueViolation" in exc_type
        or "psycopg.errors" in exc_module and "IntegrityError" in exc_type
    )


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
        elif backend_type == "postgresql":
            if not settings.database_url:
                raise ValueError(
                    "DATABASE_URL must be set when database_backend='postgresql'"
                )
            _backend = PostgresBackend(
                dsn=settings.database_url,
                min_size=settings.pg_pool_min_size,
                max_size=settings.pg_pool_max_size,
            )
        else:
            raise ValueError(f"Unknown database backend: {backend_type}")
    return _backend


def set_backend(backend: DatabaseBackend | None) -> None:
    """Override the database backend (useful for testing)."""
    global _backend
    _backend = backend
