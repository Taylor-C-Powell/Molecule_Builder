"""SQLite persistence for batch jobs."""

import json
import sqlite3
import threading
import time
import uuid


class JobDB:
    """Thread-safe SQLite database for batch job persistence."""

    def __init__(self, db_path: str = "molbuilder_jobs.db"):
        self._db_path = db_path
        self._local = threading.local()
        self._init_db()

    def _get_conn(self) -> sqlite3.Connection:
        if not hasattr(self._local, "conn") or self._local.conn is None:
            conn = sqlite3.connect(self._db_path, timeout=10)
            conn.execute("PRAGMA journal_mode=WAL")
            conn.execute("PRAGMA foreign_keys=ON")
            conn.row_factory = sqlite3.Row
            self._local.conn = conn
        return self._local.conn

    def _init_db(self) -> None:
        conn = self._get_conn()
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

    def create_job(self, user_email: str, job_type: str, input_data: dict) -> str:
        job_id = uuid.uuid4().hex
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        conn = self._get_conn()
        conn.execute(
            """INSERT INTO jobs (job_id, user_email, status, job_type, input_data, created_at, updated_at)
               VALUES (?, ?, 'pending', ?, ?, ?, ?)""",
            (job_id, user_email, job_type, json.dumps(input_data), now, now),
        )
        conn.commit()
        return job_id

    def update_status(self, job_id: str, status: str, progress_pct: float | None = None) -> None:
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        conn = self._get_conn()
        if progress_pct is not None:
            conn.execute(
                "UPDATE jobs SET status = ?, progress_pct = ?, updated_at = ? WHERE job_id = ?",
                (status, progress_pct, now, job_id),
            )
        else:
            conn.execute(
                "UPDATE jobs SET status = ?, updated_at = ? WHERE job_id = ?",
                (status, now, job_id),
            )
        conn.commit()

    def set_result(self, job_id: str, result: dict) -> None:
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        conn = self._get_conn()
        conn.execute(
            "UPDATE jobs SET status = 'completed', result_data = ?, progress_pct = 100.0, updated_at = ? WHERE job_id = ?",
            (json.dumps(result), now, job_id),
        )
        conn.commit()

    def set_error(self, job_id: str, error: str) -> None:
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        conn = self._get_conn()
        conn.execute(
            "UPDATE jobs SET status = 'failed', error = ?, updated_at = ? WHERE job_id = ?",
            (error, now, job_id),
        )
        conn.commit()

    def get_job(self, job_id: str) -> dict | None:
        conn = self._get_conn()
        row = conn.execute("SELECT * FROM jobs WHERE job_id = ?", (job_id,)).fetchone()
        if row is None:
            return None
        d = dict(row)
        if d.get("result_data"):
            d["result_data"] = json.loads(d["result_data"])
        if d.get("input_data"):
            d["input_data"] = json.loads(d["input_data"])
        return d

    def list_jobs(self, user_email: str, limit: int = 20, offset: int = 0) -> tuple[list[dict], int]:
        conn = self._get_conn()
        total = conn.execute(
            "SELECT COUNT(*) FROM jobs WHERE user_email = ?", (user_email,)
        ).fetchone()[0]
        rows = conn.execute(
            "SELECT job_id, status, job_type, progress_pct, created_at, updated_at FROM jobs WHERE user_email = ? ORDER BY created_at DESC LIMIT ? OFFSET ?",
            (user_email, limit, offset),
        ).fetchall()
        return [dict(r) for r in rows], total

    def cancel_job(self, job_id: str, user_email: str) -> bool:
        """Cancel a pending or running job. Returns True if cancelled."""
        conn = self._get_conn()
        cursor = conn.execute(
            "UPDATE jobs SET status = 'cancelled' WHERE job_id = ? AND user_email = ? AND status IN ('pending', 'running')",
            (job_id, user_email),
        )
        conn.commit()
        return cursor.rowcount > 0

    def close(self) -> None:
        if hasattr(self._local, "conn") and self._local.conn is not None:
            self._local.conn.close()
            self._local.conn = None


_job_db: JobDB | None = None


def get_job_db() -> JobDB:
    global _job_db
    if _job_db is None:
        _job_db = JobDB()
    return _job_db


def set_job_db(db: JobDB | None) -> None:
    global _job_db
    _job_db = db
