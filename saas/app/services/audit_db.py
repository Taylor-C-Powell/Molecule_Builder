"""Immutable audit trail for 21 CFR Part 11 compliance -- SQLite or PostgreSQL.

HMAC signatures are always computed from the epoch float timestamp to
ensure portability between SQLite (REAL) and PostgreSQL (TIMESTAMPTZ).
"""

from __future__ import annotations

import hashlib
import hmac as hmac_mod
import sqlite3
import threading
import time
from datetime import datetime, timezone

from app.services.database import DatabaseBackend


class AuditDB:
    """Thread-safe, insert-only audit database."""

    def __init__(
        self,
        db_path: str = "molbuilder_audit.db",
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

    def _init_sqlite(self, db_path: str) -> None:
        conn = sqlite3.connect(db_path)
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA foreign_keys=ON")
        conn.execute("""
            CREATE TABLE IF NOT EXISTS audit_log (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                timestamp REAL NOT NULL,
                user_email TEXT NOT NULL,
                user_role TEXT NOT NULL DEFAULT 'chemist',
                user_tier TEXT NOT NULL DEFAULT 'free',
                action TEXT NOT NULL,
                input_summary TEXT NOT NULL DEFAULT '',
                output_hash TEXT NOT NULL DEFAULT '',
                status_code INTEGER NOT NULL DEFAULT 0,
                latency_ms REAL NOT NULL DEFAULT 0,
                ip_address TEXT NOT NULL DEFAULT '',
                signature_hash TEXT NOT NULL
            )
        """)
        conn.execute("CREATE INDEX IF NOT EXISTS idx_audit_timestamp ON audit_log(timestamp)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_audit_user ON audit_log(user_email)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_audit_action ON audit_log(action)")
        conn.commit()
        conn.close()

    # ------------------------------------------------------------------ #
    # Signature helpers (static -- shared by both backends)
    # ------------------------------------------------------------------ #

    @staticmethod
    def compute_signature(user_email: str, timestamp: float, action: str, input_summary: str) -> str:
        """Compute HMAC-SHA256 electronic signature using server secret."""
        from app.config import settings
        data = f"{user_email}|{timestamp}|{action}|{input_summary}"
        secret = settings.audit_hmac_secret
        if secret:
            return hmac_mod.new(secret.encode(), data.encode(), hashlib.sha256).hexdigest()
        return hashlib.sha256(data.encode()).hexdigest()

    @staticmethod
    def compute_output_hash(response_body: str) -> str:
        """Compute SHA-256 hash of response body."""
        return hashlib.sha256(response_body.encode()).hexdigest()

    # ------------------------------------------------------------------ #
    # Write
    # ------------------------------------------------------------------ #

    def record(
        self,
        user_email: str,
        user_role: str,
        user_tier: str,
        action: str,
        input_summary: str = "",
        output_hash: str = "",
        status_code: int = 0,
        latency_ms: float = 0,
        ip_address: str = "",
    ) -> int:
        """Insert an immutable audit record. Returns the record ID."""
        ts = time.time()
        signature = self.compute_signature(user_email, ts, action, input_summary)

        if self._backend.is_postgres:
            ts_dt = datetime.fromtimestamp(ts, tz=timezone.utc).isoformat()
            return self._backend.execute_insert(
                "INSERT INTO audit_log "
                "(timestamp, timestamp_epoch, user_email, user_role, user_tier, "
                "action, input_summary, output_hash, status_code, latency_ms, "
                "ip_address, signature_hash) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (ts_dt, ts, user_email, user_role, user_tier, action,
                 input_summary, output_hash, status_code, latency_ms,
                 ip_address, signature),
            )
        else:
            return self._backend.execute_insert(
                "INSERT INTO audit_log "
                "(timestamp, user_email, user_role, user_tier, action, "
                "input_summary, output_hash, status_code, latency_ms, "
                "ip_address, signature_hash) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (ts, user_email, user_role, user_tier, action,
                 input_summary, output_hash, status_code, latency_ms,
                 ip_address, signature),
            )

    # ------------------------------------------------------------------ #
    # Read
    # ------------------------------------------------------------------ #

    def _normalise_record(self, row: dict) -> dict:
        """Ensure ``timestamp`` is always an epoch float for signature verification."""
        if self._backend.is_postgres:
            # PG stores the epoch in timestamp_epoch; timestamp is TIMESTAMPTZ
            if "timestamp_epoch" in row:
                row["timestamp"] = row["timestamp_epoch"]
        return row

    def get_record(self, record_id: int) -> dict | None:
        rows = self._backend.execute(
            "SELECT * FROM audit_log WHERE id = ?", (record_id,)
        )
        if not rows:
            return None
        return self._normalise_record(rows[0])

    def query(
        self,
        start_time: float | None = None,
        end_time: float | None = None,
        user_email: str | None = None,
        action: str | None = None,
        limit: int = 100,
        offset: int = 0,
    ) -> list[dict]:
        ts_col = "timestamp_epoch" if self._backend.is_postgres else "timestamp"
        conditions: list[str] = []
        params: list = []
        if start_time is not None:
            conditions.append(f"{ts_col} >= ?")
            params.append(start_time)
        if end_time is not None:
            conditions.append(f"{ts_col} <= ?")
            params.append(end_time)
        if user_email is not None:
            conditions.append("user_email = ?")
            params.append(user_email)
        if action is not None:
            conditions.append("action = ?")
            params.append(action)

        where = ""
        if conditions:
            where = "WHERE " + " AND ".join(conditions)

        rows = self._backend.execute(
            f"SELECT * FROM audit_log {where} ORDER BY {ts_col} DESC LIMIT ? OFFSET ?",
            tuple(params + [limit, offset]),
        )
        return [self._normalise_record(r) for r in rows]

    def count(
        self,
        start_time: float | None = None,
        end_time: float | None = None,
        user_email: str | None = None,
        action: str | None = None,
    ) -> int:
        ts_col = "timestamp_epoch" if self._backend.is_postgres else "timestamp"
        conditions: list[str] = []
        params: list = []
        if start_time is not None:
            conditions.append(f"{ts_col} >= ?")
            params.append(start_time)
        if end_time is not None:
            conditions.append(f"{ts_col} <= ?")
            params.append(end_time)
        if user_email is not None:
            conditions.append("user_email = ?")
            params.append(user_email)
        if action is not None:
            conditions.append("action = ?")
            params.append(action)

        where = ""
        if conditions:
            where = "WHERE " + " AND ".join(conditions)

        rows = self._backend.execute(
            f"SELECT COUNT(*) as cnt FROM audit_log {where}", tuple(params)
        )
        return rows[0]["cnt"]

    def verify_integrity(self, record_id: int) -> dict:
        """Verify that a record's signature matches its data."""
        record = self.get_record(record_id)
        if record is None:
            return {"valid": False, "error": "Record not found"}

        expected = self.compute_signature(
            record["user_email"],
            record["timestamp"],
            record["action"],
            record["input_summary"],
        )
        valid = expected == record["signature_hash"]
        result = {
            "valid": valid,
            "record_id": record_id,
            "expected_hash": expected,
            "stored_hash": record["signature_hash"],
        }
        if not valid:
            result["error"] = "Signature mismatch: record may have been tampered with"
        return result

    def export_all(self, limit: int = 10000) -> list[dict]:
        rows = self._backend.execute(
            "SELECT * FROM audit_log ORDER BY id ASC LIMIT ?", (limit,)
        )
        return [self._normalise_record(r) for r in rows]

    def close(self) -> None:
        if self._direct_sqlite:
            self._backend.close()


# Module-level singleton (overridable in tests)
audit_db: AuditDB | None = None


def get_audit_db() -> AuditDB:
    global audit_db
    if audit_db is None:
        from app.config import settings
        if settings.database_backend == "postgresql":
            from app.services.database import get_backend
            audit_db = AuditDB(backend=get_backend())
        else:
            audit_db = AuditDB(settings.audit_db_path)
    return audit_db


def set_audit_db(db: AuditDB | None) -> None:
    global audit_db
    audit_db = db
