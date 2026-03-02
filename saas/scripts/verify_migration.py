#!/usr/bin/env python3
"""Post-migration verification: compares row counts between SQLite and
PostgreSQL, and spot-checks audit signature integrity.

Usage::

    cd saas
    python scripts/verify_migration.py
"""

from __future__ import annotations

import json
import os
import sqlite3
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def _sqlite_count(path: str, table: str) -> int | None:
    if not os.path.exists(path):
        return None
    conn = sqlite3.connect(path)
    try:
        count = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
    except sqlite3.OperationalError:
        count = None
    conn.close()
    return count


def _pg_conn():
    import psycopg
    from psycopg.rows import dict_row

    dsn = os.environ.get("DATABASE_URL")
    if not dsn:
        print("ERROR: DATABASE_URL is required")
        sys.exit(1)
    return psycopg.connect(dsn, row_factory=dict_row)


def main():
    pg = _pg_conn()
    ok = True

    tables = [
        ("api_keys", os.environ.get("USER_DB_PATH", "molbuilder_users.db")),
        ("audit_log", os.environ.get("AUDIT_DB_PATH", "molbuilder_audit.db")),
        ("jobs", os.environ.get("JOB_DB_PATH", "molbuilder_jobs.db")),
        ("library_molecules", os.environ.get("LIBRARY_DB_PATH", "molbuilder_library.db")),
        ("molecules", os.environ.get("MOLECULE_DB_PATH", "molbuilder_molecules.db")),
    ]

    print("=== Row Count Comparison ===\n")
    for table, sqlite_path in tables:
        sqlite_n = _sqlite_count(sqlite_path, table)
        pg_n = pg.execute(f"SELECT COUNT(*) as cnt FROM {table}").fetchone()["cnt"]
        match = "OK" if sqlite_n == pg_n else "MISMATCH"
        if sqlite_n is None:
            match = "SKIP (no SQLite file)"
        if match == "MISMATCH":
            ok = False
        print(f"  {table:25s}  SQLite={sqlite_n}  PG={pg_n}  [{match}]")

    # Spot-check audit signatures
    print("\n=== Audit Signature Spot-Check ===\n")
    rows = pg.execute(
        "SELECT id, timestamp_epoch, user_email, action, input_summary, signature_hash "
        "FROM audit_log ORDER BY id LIMIT 5"
    ).fetchall()

    if not rows:
        print("  No audit records to verify.")
    else:
        from app.services.audit_db import AuditDB
        for r in rows:
            expected = AuditDB.compute_signature(
                r["user_email"], r["timestamp_epoch"], r["action"], r["input_summary"]
            )
            valid = expected == r["signature_hash"]
            status = "VALID" if valid else "INVALID"
            if not valid:
                ok = False
            print(f"  Record {r['id']:5d}: {status}")

    pg.close()

    if ok:
        print("\nAll checks passed.")
    else:
        print("\nSome checks FAILED -- review above output.")
        sys.exit(1)


if __name__ == "__main__":
    main()
