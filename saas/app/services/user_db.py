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
