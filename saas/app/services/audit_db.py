"""Immutable SQLite audit trail for 21 CFR Part 11 compliance."""

import hashlib
import hmac as hmac_mod
import sqlite3
import threading
import time


class AuditDB:
    """Thread-safe, insert-only SQLite audit database."""

    def __init__(self, db_path: str = "molbuilder_audit.db"):
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

    @staticmethod
    def compute_signature(user_email: str, timestamp: float, action: str, input_summary: str) -> str:
        """Compute HMAC-SHA256 electronic signature using server secret."""
        from app.config import settings
        data = f"{user_email}|{timestamp}|{action}|{input_summary}"
        secret = settings.audit_hmac_secret
        if secret:
            return hmac_mod.new(secret.encode(), data.encode(), hashlib.sha256).hexdigest()
        # Fallback for tests where config isn't fully initialized
        return hashlib.sha256(data.encode()).hexdigest()

    @staticmethod
    def compute_output_hash(response_body: str) -> str:
        """Compute SHA-256 hash of response body."""
        return hashlib.sha256(response_body.encode()).hexdigest()

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
        conn = self._get_conn()
        cursor = conn.execute(
            """INSERT INTO audit_log
               (timestamp, user_email, user_role, user_tier, action, input_summary,
                output_hash, status_code, latency_ms, ip_address, signature_hash)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (ts, user_email, user_role, user_tier, action, input_summary,
             output_hash, status_code, latency_ms, ip_address, signature),
        )
        conn.commit()
        return cursor.lastrowid

    def get_record(self, record_id: int) -> dict | None:
        conn = self._get_conn()
        row = conn.execute("SELECT * FROM audit_log WHERE id = ?", (record_id,)).fetchone()
        if row is None:
            return None
        return dict(row)

    def query(
        self,
        start_time: float | None = None,
        end_time: float | None = None,
        user_email: str | None = None,
        action: str | None = None,
        limit: int = 100,
        offset: int = 0,
    ) -> list[dict]:
        conditions = []
        params: list = []
        if start_time is not None:
            conditions.append("timestamp >= ?")
            params.append(start_time)
        if end_time is not None:
            conditions.append("timestamp <= ?")
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

        conn = self._get_conn()
        rows = conn.execute(
            f"SELECT * FROM audit_log {where} ORDER BY timestamp DESC LIMIT ? OFFSET ?",
            params + [limit, offset],
        ).fetchall()
        return [dict(r) for r in rows]

    def count(
        self,
        start_time: float | None = None,
        end_time: float | None = None,
        user_email: str | None = None,
        action: str | None = None,
    ) -> int:
        conditions = []
        params: list = []
        if start_time is not None:
            conditions.append("timestamp >= ?")
            params.append(start_time)
        if end_time is not None:
            conditions.append("timestamp <= ?")
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

        conn = self._get_conn()
        row = conn.execute(f"SELECT COUNT(*) FROM audit_log {where}", params).fetchone()
        return row[0]

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
        conn = self._get_conn()
        rows = conn.execute(
            "SELECT * FROM audit_log ORDER BY id ASC LIMIT ?", (limit,)
        ).fetchall()
        return [dict(r) for r in rows]

    def close(self) -> None:
        if hasattr(self._local, "conn") and self._local.conn is not None:
            self._local.conn.close()
            self._local.conn = None


# Module-level singleton (overridable in tests)
audit_db: AuditDB | None = None


def get_audit_db() -> AuditDB:
    global audit_db
    if audit_db is None:
        from app.config import settings
        audit_db = AuditDB(settings.audit_db_path)
    return audit_db


def set_audit_db(db: AuditDB) -> None:
    global audit_db
    audit_db = db
