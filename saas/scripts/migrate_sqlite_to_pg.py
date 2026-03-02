#!/usr/bin/env python3
"""One-time migration script: reads all 5 SQLite databases and bulk-inserts
their data into the PostgreSQL database.

Prerequisites:
  - PostgreSQL schema already created via ``alembic upgrade head``
  - Environment variables set: DATABASE_URL, plus the 5 *_DB_PATH vars

Usage::

    cd saas
    python scripts/migrate_sqlite_to_pg.py

The script is idempotent -- it will skip rows that already exist (by PK).
"""

from __future__ import annotations

import json
import os
import sqlite3
import sys
import time
from datetime import datetime, timezone

# Ensure app is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def _connect_sqlite(path: str) -> sqlite3.Connection:
    if not os.path.exists(path):
        print(f"  SKIP {path} (file not found)")
        return None
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    return conn


def _pg_conn():
    """Return a psycopg connection using DATABASE_URL."""
    import psycopg
    from psycopg.rows import dict_row

    dsn = os.environ.get("DATABASE_URL")
    if not dsn:
        print("ERROR: DATABASE_URL environment variable is required")
        sys.exit(1)
    return psycopg.connect(dsn, row_factory=dict_row)


def migrate_users(pg) -> int:
    path = os.environ.get("USER_DB_PATH", "molbuilder_users.db")
    conn = _connect_sqlite(path)
    if conn is None:
        return 0
    rows = conn.execute("SELECT * FROM api_keys").fetchall()
    count = 0
    for r in rows:
        d = dict(r)
        ts = datetime.fromtimestamp(d["created_at"], tz=timezone.utc)
        active = bool(d["active"])
        try:
            pg.execute(
                "INSERT INTO api_keys "
                "(id, key_hash, email, tier, role, created_at, active, "
                "stripe_customer_id, stripe_subscription_id, subscription_status) "
                "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s) "
                "ON CONFLICT (id) DO NOTHING",
                (d["id"], d["key_hash"], d["email"], d["tier"], d["role"],
                 ts, active, d.get("stripe_customer_id"),
                 d.get("stripe_subscription_id"),
                 d.get("subscription_status", "none")),
            )
            count += 1
        except Exception as exc:
            print(f"  WARN user row {d['id']}: {exc}")
    # Reset sequence
    pg.execute(
        "SELECT setval('api_keys_id_seq', COALESCE((SELECT MAX(id) FROM api_keys), 0))"
    )
    conn.close()
    return count


def migrate_audit(pg) -> int:
    path = os.environ.get("AUDIT_DB_PATH", "molbuilder_audit.db")
    conn = _connect_sqlite(path)
    if conn is None:
        return 0
    rows = conn.execute("SELECT * FROM audit_log").fetchall()
    count = 0
    for r in rows:
        d = dict(r)
        ts_float = d["timestamp"]
        ts_dt = datetime.fromtimestamp(ts_float, tz=timezone.utc)
        try:
            pg.execute(
                "INSERT INTO audit_log "
                "(id, timestamp, timestamp_epoch, user_email, user_role, user_tier, "
                "action, input_summary, output_hash, status_code, latency_ms, "
                "ip_address, signature_hash) "
                "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) "
                "ON CONFLICT (id) DO NOTHING",
                (d["id"], ts_dt, ts_float, d["user_email"], d["user_role"],
                 d["user_tier"], d["action"], d["input_summary"],
                 d["output_hash"], d["status_code"], d["latency_ms"],
                 d["ip_address"], d["signature_hash"]),
            )
            count += 1
        except Exception as exc:
            print(f"  WARN audit row {d['id']}: {exc}")
    pg.execute(
        "SELECT setval('audit_log_id_seq', COALESCE((SELECT MAX(id) FROM audit_log), 0))"
    )
    conn.close()
    return count


def migrate_jobs(pg) -> int:
    path = os.environ.get("JOB_DB_PATH", "molbuilder_jobs.db")
    conn = _connect_sqlite(path)
    if conn is None:
        return 0
    rows = conn.execute("SELECT * FROM jobs").fetchall()
    count = 0
    for r in rows:
        d = dict(r)
        input_data = json.loads(d["input_data"]) if d["input_data"] else {}
        result_data = json.loads(d["result_data"]) if d["result_data"] else None
        try:
            pg.execute(
                "INSERT INTO jobs "
                "(job_id, user_email, status, job_type, input_data, result_data, "
                "error, created_at, updated_at, progress_pct) "
                "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s) "
                "ON CONFLICT (job_id) DO NOTHING",
                (d["job_id"], d["user_email"], d["status"], d["job_type"],
                 json.dumps(input_data), json.dumps(result_data) if result_data else None,
                 d.get("error"), d["created_at"], d["updated_at"], d["progress_pct"]),
            )
            count += 1
        except Exception as exc:
            print(f"  WARN job {d['job_id']}: {exc}")
    conn.close()
    return count


def migrate_library(pg) -> int:
    path = os.environ.get("LIBRARY_DB_PATH", "molbuilder_library.db")
    conn = _connect_sqlite(path)
    if conn is None:
        return 0
    rows = conn.execute("SELECT * FROM library_molecules").fetchall()
    count = 0
    for r in rows:
        d = dict(r)
        tags = json.loads(d["tags"]) if d["tags"] else []
        props = json.loads(d["properties_json"]) if d.get("properties_json") else {}
        try:
            pg.execute(
                "INSERT INTO library_molecules "
                "(id, user_email, smiles, name, tags, notes, properties_json, "
                "created_at, updated_at) "
                "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s) "
                "ON CONFLICT (id) DO NOTHING",
                (d["id"], d["user_email"], d["smiles"], d.get("name"),
                 json.dumps(tags), d.get("notes"), json.dumps(props),
                 d["created_at"], d["updated_at"]),
            )
            count += 1
        except Exception as exc:
            print(f"  WARN library row {d['id']}: {exc}")
    pg.execute(
        "SELECT setval('library_molecules_id_seq', "
        "COALESCE((SELECT MAX(id) FROM library_molecules), 0))"
    )
    conn.close()
    return count


def migrate_molecules(pg) -> int:
    path = os.environ.get("MOLECULE_DB_PATH", "molbuilder_molecules.db")
    conn = _connect_sqlite(path)
    if conn is None:
        return 0
    rows = conn.execute("SELECT * FROM molecules").fetchall()
    count = 0
    for r in rows:
        d = dict(r)
        created = datetime.fromtimestamp(d["created_at"], tz=timezone.utc)
        accessed = datetime.fromtimestamp(d["accessed_at"], tz=timezone.utc)
        try:
            pg.execute(
                "INSERT INTO molecules "
                "(mol_id, smiles, created_at, accessed_at) "
                "VALUES (%s, %s, %s, %s) "
                "ON CONFLICT (mol_id) DO NOTHING",
                (d["mol_id"], d["smiles"], created, accessed),
            )
            count += 1
        except Exception as exc:
            print(f"  WARN molecule {d['mol_id']}: {exc}")
    conn.close()
    return count


def main():
    print("=== SQLite -> PostgreSQL Migration ===")
    print()

    pg = _pg_conn()
    start = time.time()

    print("[1/5] Migrating api_keys (user_db)...")
    n = migrate_users(pg)
    print(f"  -> {n} rows")

    print("[2/5] Migrating audit_log (audit_db)...")
    n = migrate_audit(pg)
    print(f"  -> {n} rows")

    print("[3/5] Migrating jobs (job_db)...")
    n = migrate_jobs(pg)
    print(f"  -> {n} rows")

    print("[4/5] Migrating library_molecules (library_db)...")
    n = migrate_library(pg)
    print(f"  -> {n} rows")

    print("[5/5] Migrating molecules (molecule_store)...")
    n = migrate_molecules(pg)
    print(f"  -> {n} rows")

    pg.commit()
    pg.close()

    elapsed = time.time() - start
    print(f"\nMigration completed in {elapsed:.1f}s")
    print("Run scripts/verify_migration.py to confirm row counts and audit signatures.")


if __name__ == "__main__":
    main()
