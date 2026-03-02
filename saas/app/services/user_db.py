"""Persistence for API keys -- SQLite or PostgreSQL."""

from __future__ import annotations

import sqlite3
import threading
import time

from app.services.database import DatabaseBackend, DatabaseIntegrityError


class UserDB:
    """Thread-safe database for API key persistence.

    When *backend* is supplied the class delegates all queries through
    the :class:`DatabaseBackend` abstraction (works for both SQLite and
    PostgreSQL).  When omitted it falls back to a direct SQLite
    connection for backwards compatibility with tests and local dev.
    """

    def __init__(
        self,
        db_path: str = "molbuilder_users.db",
        backend: DatabaseBackend | None = None,
    ):
        if backend is not None:
            self._backend = backend
            self._direct_sqlite = False
        else:
            from app.services.database import SQLiteBackend
            self._backend = SQLiteBackend(db_path)
            self._direct_sqlite = True
            self._init_sqlite(db_path)

    # ------------------------------------------------------------------ #
    # SQLite-only bootstrap (CREATE TABLE + column migration)
    # ------------------------------------------------------------------ #

    def _init_sqlite(self, db_path: str) -> None:
        conn = sqlite3.connect(db_path)
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA foreign_keys=ON")
        conn.execute("""
            CREATE TABLE IF NOT EXISTS api_keys (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                key_hash TEXT NOT NULL UNIQUE,
                email TEXT NOT NULL,
                tier TEXT NOT NULL DEFAULT 'free',
                role TEXT NOT NULL DEFAULT 'chemist',
                created_at REAL NOT NULL,
                active INTEGER NOT NULL DEFAULT 1
            )
        """)
        conn.execute("CREATE INDEX IF NOT EXISTS idx_keys_email ON api_keys(email)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_keys_active ON api_keys(active)")
        existing = {r[1] for r in conn.execute("PRAGMA table_info(api_keys)").fetchall()}
        if "stripe_customer_id" not in existing:
            conn.execute("ALTER TABLE api_keys ADD COLUMN stripe_customer_id TEXT")
        if "stripe_subscription_id" not in existing:
            conn.execute("ALTER TABLE api_keys ADD COLUMN stripe_subscription_id TEXT")
        if "subscription_status" not in existing:
            conn.execute("ALTER TABLE api_keys ADD COLUMN subscription_status TEXT DEFAULT 'none'")
        conn.commit()
        conn.close()

    # ------------------------------------------------------------------ #
    # CRUD operations (backend-agnostic)
    # ------------------------------------------------------------------ #

    def insert_key(self, key_hash: str, email: str, tier: str, role: str) -> int:
        """Insert a new API key record. Returns the record ID."""
        return self._backend.execute_insert(
            "INSERT INTO api_keys (key_hash, email, tier, role, created_at) "
            "VALUES (?, ?, ?, ?, ?)",
            (key_hash, email, tier, role, time.time()),
        )

    def delete_by_email(self, email: str) -> int:
        """Soft-delete all keys for a given email. Returns number affected."""
        return self._backend.execute_update(
            "UPDATE api_keys SET active = 0 WHERE email = ? AND active = 1",
            (email,),
        )

    def delete_by_hash(self, key_hash: str) -> int:
        """Soft-delete a specific key by its hash. Returns number affected."""
        return self._backend.execute_update(
            "UPDATE api_keys SET active = 0 WHERE key_hash = ? AND active = 1",
            (key_hash,),
        )

    def update_tier(self, email: str, tier: str) -> int:
        """Update the tier for all active keys belonging to an email."""
        return self._backend.execute_update(
            "UPDATE api_keys SET tier = ? WHERE email = ? AND active = 1",
            (tier, email),
        )

    def update_stripe_customer(self, email: str, customer_id: str) -> None:
        """Associate a Stripe customer ID with a user's active keys."""
        self._backend.execute_update(
            "UPDATE api_keys SET stripe_customer_id = ? WHERE email = ? AND active = 1",
            (customer_id, email),
        )

    def update_subscription(
        self, email: str, subscription_id: str | None, status: str, tier: str
    ) -> None:
        """Update subscription info and tier for a user."""
        self._backend.execute_update(
            "UPDATE api_keys "
            "SET stripe_subscription_id = ?, subscription_status = ?, tier = ? "
            "WHERE email = ? AND active = 1",
            (subscription_id, status, tier, email),
        )

    def get_by_stripe_customer(self, customer_id: str) -> dict | None:
        """Look up user by Stripe customer ID."""
        rows = self._backend.execute(
            "SELECT email, tier, role, stripe_customer_id, "
            "stripe_subscription_id, subscription_status "
            "FROM api_keys WHERE stripe_customer_id = ? AND active = 1 LIMIT 1",
            (customer_id,),
        )
        return rows[0] if rows else None

    def get_stripe_info(self, email: str) -> dict | None:
        """Get Stripe-related info for a user."""
        rows = self._backend.execute(
            "SELECT email, tier, stripe_customer_id, "
            "stripe_subscription_id, subscription_status "
            "FROM api_keys WHERE email = ? AND active = 1 LIMIT 1",
            (email,),
        )
        return rows[0] if rows else None

    def load_all_active(self) -> list[dict]:
        """Load all active API key records."""
        return self._backend.execute(
            "SELECT key_hash, email, tier, role FROM api_keys WHERE active = 1"
        )

    def close(self) -> None:
        if self._direct_sqlite:
            self._backend.close()


# Module-level singleton (overridable in tests)
_user_db: UserDB | None = None


def get_user_db() -> UserDB:
    global _user_db
    if _user_db is None:
        from app.config import settings
        if settings.database_backend == "postgresql":
            from app.services.database import get_backend
            _user_db = UserDB(backend=get_backend())
        else:
            _user_db = UserDB(settings.user_db_path)
    return _user_db


def set_user_db(db: UserDB | None) -> None:
    global _user_db
    _user_db = db
