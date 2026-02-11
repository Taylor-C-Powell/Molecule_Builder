"""SQLite persistence for API keys."""

import sqlite3
import threading
import time


class UserDB:
    """Thread-safe SQLite database for API key persistence."""

    def __init__(self, db_path: str = "molbuilder_users.db"):
        self._db_path = db_path
        self._local = threading.local()
        self._init_db()

    def _get_conn(self) -> sqlite3.Connection:
        if not hasattr(self._local, "conn") or self._local.conn is None:
            conn = sqlite3.connect(self._db_path)
            conn.execute("PRAGMA journal_mode=WAL")
            conn.execute("PRAGMA foreign_keys=ON")
            conn.row_factory = sqlite3.Row
            self._local.conn = conn
        return self._local.conn

    def _init_db(self) -> None:
        conn = self._get_conn()
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
        # Migrate: add Stripe columns if missing
        existing = {r[1] for r in conn.execute("PRAGMA table_info(api_keys)").fetchall()}
        if "stripe_customer_id" not in existing:
            conn.execute("ALTER TABLE api_keys ADD COLUMN stripe_customer_id TEXT")
        if "stripe_subscription_id" not in existing:
            conn.execute("ALTER TABLE api_keys ADD COLUMN stripe_subscription_id TEXT")
        if "subscription_status" not in existing:
            conn.execute("ALTER TABLE api_keys ADD COLUMN subscription_status TEXT DEFAULT 'none'")
        conn.commit()

    def insert_key(self, key_hash: str, email: str, tier: str, role: str) -> int:
        """Insert a new API key record. Returns the record ID."""
        conn = self._get_conn()
        cursor = conn.execute(
            """INSERT INTO api_keys (key_hash, email, tier, role, created_at)
               VALUES (?, ?, ?, ?, ?)""",
            (key_hash, email, tier, role, time.time()),
        )
        conn.commit()
        return cursor.lastrowid

    def delete_by_email(self, email: str) -> int:
        """Soft-delete all keys for a given email. Returns number affected."""
        conn = self._get_conn()
        cursor = conn.execute(
            "UPDATE api_keys SET active = 0 WHERE email = ? AND active = 1",
            (email,),
        )
        conn.commit()
        return cursor.rowcount

    def update_tier(self, email: str, tier: str) -> int:
        """Update the tier for all active keys belonging to an email."""
        conn = self._get_conn()
        cursor = conn.execute(
            "UPDATE api_keys SET tier = ? WHERE email = ? AND active = 1",
            (tier, email),
        )
        conn.commit()
        return cursor.rowcount

    def update_stripe_customer(self, email: str, customer_id: str) -> None:
        """Associate a Stripe customer ID with a user's active keys."""
        conn = self._get_conn()
        conn.execute(
            "UPDATE api_keys SET stripe_customer_id = ? WHERE email = ? AND active = 1",
            (customer_id, email),
        )
        conn.commit()

    def update_subscription(
        self, email: str, subscription_id: str | None, status: str, tier: str
    ) -> None:
        """Update subscription info and tier for a user."""
        conn = self._get_conn()
        conn.execute(
            """UPDATE api_keys
               SET stripe_subscription_id = ?, subscription_status = ?, tier = ?
               WHERE email = ? AND active = 1""",
            (subscription_id, status, tier, email),
        )
        conn.commit()

    def get_by_stripe_customer(self, customer_id: str) -> dict | None:
        """Look up user by Stripe customer ID."""
        conn = self._get_conn()
        row = conn.execute(
            """SELECT email, tier, role, stripe_customer_id,
                      stripe_subscription_id, subscription_status
               FROM api_keys WHERE stripe_customer_id = ? AND active = 1 LIMIT 1""",
            (customer_id,),
        ).fetchone()
        return dict(row) if row else None

    def get_stripe_info(self, email: str) -> dict | None:
        """Get Stripe-related info for a user."""
        conn = self._get_conn()
        row = conn.execute(
            """SELECT email, tier, stripe_customer_id,
                      stripe_subscription_id, subscription_status
               FROM api_keys WHERE email = ? AND active = 1 LIMIT 1""",
            (email,),
        ).fetchone()
        return dict(row) if row else None

    def load_all_active(self) -> list[dict]:
        """Load all active API key records."""
        conn = self._get_conn()
        rows = conn.execute(
            "SELECT key_hash, email, tier, role FROM api_keys WHERE active = 1"
        ).fetchall()
        return [dict(r) for r in rows]

    def close(self) -> None:
        if hasattr(self._local, "conn") and self._local.conn is not None:
            self._local.conn.close()
            self._local.conn = None


# Module-level singleton (overridable in tests)
_user_db: UserDB | None = None


def get_user_db() -> UserDB:
    global _user_db
    if _user_db is None:
        from app.config import settings
        _user_db = UserDB(settings.user_db_path)
    return _user_db


def set_user_db(db: UserDB | None) -> None:
    global _user_db
    _user_db = db
