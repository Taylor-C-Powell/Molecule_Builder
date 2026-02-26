"""SQLite persistence for molecule libraries."""

from __future__ import annotations

import json
import sqlite3
import threading
import time


class LibraryDB:
    """Thread-safe SQLite database for molecule library persistence.

    Uses a single shared connection with check_same_thread=False and a
    threading lock to serialise writes.  This avoids stale thread-local
    connections that cannot be cleaned up across test runs on Windows.
    """

    def __init__(self, db_path: str = "molbuilder_library.db"):
        self._db_path = db_path
        self._lock = threading.Lock()
        self._conn = sqlite3.connect(
            db_path, timeout=30, check_same_thread=False,
        )
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("PRAGMA foreign_keys=ON")
        self._conn.execute("PRAGMA busy_timeout=10000")
        self._conn.row_factory = sqlite3.Row
        self._init_db()

    def _init_db(self) -> None:
        self._conn.execute("""
            CREATE TABLE IF NOT EXISTS library_molecules (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                user_email TEXT NOT NULL,
                smiles TEXT NOT NULL,
                name TEXT,
                tags TEXT NOT NULL DEFAULT '[]',
                notes TEXT,
                properties_json TEXT,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL
            )
        """)
        self._conn.execute(
            "CREATE UNIQUE INDEX IF NOT EXISTS idx_lib_user_smiles "
            "ON library_molecules(user_email, smiles)"
        )
        self._conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_lib_user ON library_molecules(user_email)"
        )
        self._conn.commit()

    def save_molecule(self, user_email: str, smiles: str, name: str | None = None,
                      tags: list[str] | None = None, notes: str | None = None,
                      properties: dict | None = None) -> int:
        """Save a molecule. Returns the row ID. Raises IntegrityError on duplicate."""
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        with self._lock:
            cursor = self._conn.execute(
                """INSERT INTO library_molecules
                   (user_email, smiles, name, tags, notes, properties_json, created_at, updated_at)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
                (user_email, smiles, name, json.dumps(tags or []),
                 notes, json.dumps(properties or {}), now, now),
            )
            self._conn.commit()
            return cursor.lastrowid

    def get_molecule(self, mol_id: int, user_email: str) -> dict | None:
        with self._lock:
            row = self._conn.execute(
                "SELECT * FROM library_molecules WHERE id = ? AND user_email = ?",
                (mol_id, user_email),
            ).fetchone()
        if row is None:
            return None
        return self._row_to_dict(row)

    def list_molecules(self, user_email: str, tag: str | None = None,
                       search: str | None = None,
                       limit: int = 20, offset: int = 0) -> tuple[list[dict], int]:
        conditions = ["user_email = ?"]
        params: list = [user_email]

        if tag:
            # JSON array contains check
            conditions.append("tags LIKE ?")
            params.append(f'%"{tag}"%')

        if search:
            conditions.append("(name LIKE ? OR smiles LIKE ?)")
            params.append(f"%{search}%")
            params.append(f"%{search}%")

        where = " AND ".join(conditions)

        with self._lock:
            total = self._conn.execute(
                f"SELECT COUNT(*) FROM library_molecules WHERE {where}", params
            ).fetchone()[0]

            rows = self._conn.execute(
                f"SELECT * FROM library_molecules WHERE {where} "
                f"ORDER BY updated_at DESC LIMIT ? OFFSET ?",
                params + [limit, offset],
            ).fetchall()

        return [self._row_to_dict(r) for r in rows], total

    def update_molecule(self, mol_id: int, user_email: str,
                        name: str | None = None, tags: list[str] | None = None,
                        notes: str | None = None) -> bool:
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        updates = ["updated_at = ?"]
        params: list = [now]
        if name is not None:
            updates.append("name = ?")
            params.append(name)
        if tags is not None:
            updates.append("tags = ?")
            params.append(json.dumps(tags))
        if notes is not None:
            updates.append("notes = ?")
            params.append(notes)

        params.extend([mol_id, user_email])
        with self._lock:
            cursor = self._conn.execute(
                f"UPDATE library_molecules SET {', '.join(updates)} "
                f"WHERE id = ? AND user_email = ?",
                params,
            )
            self._conn.commit()
        return cursor.rowcount > 0

    def delete_molecule(self, mol_id: int, user_email: str) -> bool:
        with self._lock:
            cursor = self._conn.execute(
                "DELETE FROM library_molecules WHERE id = ? AND user_email = ?",
                (mol_id, user_email),
            )
            self._conn.commit()
        return cursor.rowcount > 0

    def count_molecules(self, user_email: str) -> int:
        with self._lock:
            return self._conn.execute(
                "SELECT COUNT(*) FROM library_molecules WHERE user_email = ?",
                (user_email,),
            ).fetchone()[0]

    def _row_to_dict(self, row: sqlite3.Row) -> dict:
        d = dict(row)
        if d.get("tags"):
            d["tags"] = json.loads(d["tags"])
        else:
            d["tags"] = []
        if d.get("properties_json"):
            d["properties"] = json.loads(d["properties_json"])
        else:
            d["properties"] = {}
        d.pop("properties_json", None)
        return d

    def close(self) -> None:
        try:
            self._conn.close()
        except Exception:
            pass


_library_db: LibraryDB | None = None


def get_library_db() -> LibraryDB:
    global _library_db
    if _library_db is None:
        _library_db = LibraryDB()
    return _library_db


def set_library_db(db: LibraryDB | None) -> None:
    global _library_db
    _library_db = db
