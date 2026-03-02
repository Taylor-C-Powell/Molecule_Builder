"""Persistence for batch jobs -- SQLite or PostgreSQL."""

from __future__ import annotations

import json
import sqlite3
import threading
import time
import uuid

from app.services.database import DatabaseBackend


class JobDB:
    """Thread-safe database for batch job persistence."""

    def __init__(
        self,
        db_path: str = "molbuilder_jobs.db",
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
        conn = sqlite3.connect(db_path, timeout=10)
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA foreign_keys=ON")
        conn.execute("""
            CREATE TABLE IF NOT EXISTS jobs (
                job_id TEXT PRIMARY KEY,
                user_email TEXT NOT NULL,
                status TEXT NOT NULL DEFAULT 'pending',
                job_type TEXT NOT NULL,
                input_data TEXT NOT NULL,
                result_data TEXT,
                error TEXT,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL,
                progress_pct REAL NOT NULL DEFAULT 0.0
            )
        """)
        conn.execute("CREATE INDEX IF NOT EXISTS idx_jobs_user ON jobs(user_email)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_jobs_status ON jobs(status)")
        conn.commit()
        conn.close()

    # ------------------------------------------------------------------ #
    # JSON helpers
    # ------------------------------------------------------------------ #

    def _encode_json(self, data: dict) -> str | dict:
        """Encode a dict for storage -- raw dict for PG JSONB, JSON string for SQLite."""
        if self._backend.supports_native_json:
            return data
        return json.dumps(data)

    def _decode_json(self, value) -> dict | None:
        """Decode stored JSON -- PG returns dict already, SQLite returns string."""
        if value is None:
            return None
        if isinstance(value, dict):
            return value
        return json.loads(value)

    # ------------------------------------------------------------------ #
    # CRUD operations
    # ------------------------------------------------------------------ #

    def create_job(self, user_email: str, job_type: str, input_data: dict) -> str:
        job_id = uuid.uuid4().hex
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        encoded = self._encode_json(input_data)
        # JobDB uses TEXT PK (job_id), not an auto-increment id.
        # execute_insert appends RETURNING id on PG which won't match,
        # so use execute_update which just runs the statement.
        self._backend.execute_update(
            "INSERT INTO jobs (job_id, user_email, status, job_type, "
            "input_data, created_at, updated_at) "
            "VALUES (?, ?, 'pending', ?, ?, ?, ?)",
            (job_id, user_email, job_type, encoded, now, now),
        )
        return job_id

    def update_status(self, job_id: str, status: str, progress_pct: float | None = None) -> None:
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        if progress_pct is not None:
            self._backend.execute_update(
                "UPDATE jobs SET status = ?, progress_pct = ?, updated_at = ? WHERE job_id = ?",
                (status, progress_pct, now, job_id),
            )
        else:
            self._backend.execute_update(
                "UPDATE jobs SET status = ?, updated_at = ? WHERE job_id = ?",
                (status, now, job_id),
            )

    def set_result(self, job_id: str, result: dict) -> None:
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        encoded = self._encode_json(result)
        self._backend.execute_update(
            "UPDATE jobs SET status = 'completed', result_data = ?, "
            "progress_pct = 100.0, updated_at = ? WHERE job_id = ?",
            (encoded, now, job_id),
        )

    def set_error(self, job_id: str, error: str) -> None:
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        self._backend.execute_update(
            "UPDATE jobs SET status = 'failed', error = ?, updated_at = ? WHERE job_id = ?",
            (error, now, job_id),
        )

    def get_job(self, job_id: str) -> dict | None:
        rows = self._backend.execute(
            "SELECT * FROM jobs WHERE job_id = ?", (job_id,)
        )
        if not rows:
            return None
        d = rows[0]
        d["result_data"] = self._decode_json(d.get("result_data"))
        d["input_data"] = self._decode_json(d.get("input_data"))
        return d

    def list_jobs(self, user_email: str, limit: int = 20, offset: int = 0) -> tuple[list[dict], int]:
        total_rows = self._backend.execute(
            "SELECT COUNT(*) as cnt FROM jobs WHERE user_email = ?", (user_email,)
        )
        total = total_rows[0]["cnt"]
        rows = self._backend.execute(
            "SELECT job_id, status, job_type, progress_pct, created_at, updated_at "
            "FROM jobs WHERE user_email = ? ORDER BY created_at DESC LIMIT ? OFFSET ?",
            (user_email, limit, offset),
        )
        return rows, total

    def cancel_job(self, job_id: str, user_email: str) -> bool:
        """Cancel a pending or running job. Returns True if cancelled."""
        affected = self._backend.execute_update(
            "UPDATE jobs SET status = 'cancelled' "
            "WHERE job_id = ? AND user_email = ? AND status IN ('pending', 'running')",
            (job_id, user_email),
        )
        return affected > 0

    def close(self) -> None:
        if self._direct_sqlite:
            self._backend.close()


_job_db: JobDB | None = None


def get_job_db() -> JobDB:
    global _job_db
    if _job_db is None:
        from app.config import settings
        if settings.database_backend == "postgresql":
            from app.services.database import get_backend
            _job_db = JobDB(backend=get_backend())
        else:
            _job_db = JobDB(settings.job_db_path)
    return _job_db


def set_job_db(db: JobDB | None) -> None:
    global _job_db
    _job_db = db
