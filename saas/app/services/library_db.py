"""Persistence for molecule libraries -- SQLite or PostgreSQL."""

from __future__ import annotations

import json
import sqlite3
import threading
import time

from app.services.database import DatabaseBackend, DatabaseIntegrityError


class LibraryDB:
    """Thread-safe database for molecule library persistence.

    SQLite mode uses a single shared connection with a threading lock.
    PostgreSQL mode delegates concurrency to the connection pool.
    """

    def __init__(
        self,
        db_path: str = "molbuilder_library.db",
        backend: DatabaseBackend | None = None,
    ):
        if backend is not None:
            self._backend = backend
            self._direct_sqlite = False
            self._lock: threading.Lock | None = None
        else:
            # Direct SQLite path -- shared connection with lock
            self._direct_sqlite = True
            self._lock = threading.Lock()
            self._conn = sqlite3.connect(
                db_path, timeout=30, check_same_thread=False,
            )
            self._conn.execute("PRAGMA journal_mode=WAL")
            self._conn.execute("PRAGMA foreign_keys=ON")
            self._conn.execute("PRAGMA busy_timeout=10000")
            self._conn.row_factory = sqlite3.Row
            self._init_sqlite()
            # Wrap in a backend for the query methods
            from app.services.database import SQLiteBackend
            self._backend = SQLiteBackend(db_path)

    def _init_sqlite(self) -> None:
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

    # ------------------------------------------------------------------ #
    # JSON helpers
    # ------------------------------------------------------------------ #

    def _encode_json(self, data) -> str | list | dict:
        if self._backend.supports_native_json:
            return data
        return json.dumps(data)

    def _decode_json(self, value, default=None):
        if value is None:
            return default
        if isinstance(value, (dict, list)):
            return value
        return json.loads(value)

    # ------------------------------------------------------------------ #
    # CRUD operations
    # ------------------------------------------------------------------ #

    def save_molecule(
        self,
        user_email: str,
        smiles: str,
        name: str | None = None,
        tags: list[str] | None = None,
        notes: str | None = None,
        properties: dict | None = None,
    ) -> int:
        """Save a molecule. Returns the row ID. Raises DatabaseIntegrityError on duplicate."""
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        tags_val = self._encode_json(tags or [])
        props_val = self._encode_json(properties or {})

        if self._direct_sqlite:
            with self._lock:
                try:
                    cursor = self._conn.execute(
                        "INSERT INTO library_molecules "
                        "(user_email, smiles, name, tags, notes, properties_json, "
                        "created_at, updated_at) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                        (user_email, smiles, name, tags_val, notes, props_val, now, now),
                    )
                except sqlite3.IntegrityError as exc:
                    raise DatabaseIntegrityError(str(exc)) from exc
                self._conn.commit()
                return cursor.lastrowid
        else:
            return self._backend.execute_insert(
                "INSERT INTO library_molecules "
                "(user_email, smiles, name, tags, notes, properties_json, "
                "created_at, updated_at) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                (user_email, smiles, name, tags_val, notes, props_val, now, now),
            )

    def get_molecule(self, mol_id: int, user_email: str) -> dict | None:
        if self._direct_sqlite:
            with self._lock:
                row = self._conn.execute(
                    "SELECT * FROM library_molecules WHERE id = ? AND user_email = ?",
                    (mol_id, user_email),
                ).fetchone()
            if row is None:
                return None
            return self._row_to_dict(dict(row))
        else:
            rows = self._backend.execute(
                "SELECT * FROM library_molecules WHERE id = ? AND user_email = ?",
                (mol_id, user_email),
            )
            if not rows:
                return None
            return self._row_to_dict(rows[0])

    def list_molecules(
        self,
        user_email: str,
        tag: str | None = None,
        search: str | None = None,
        limit: int = 20,
        offset: int = 0,
    ) -> tuple[list[dict], int]:
        conditions = ["user_email = ?"]
        params: list = [user_email]

        if tag:
            if self._backend.is_postgres:
                conditions.append("tags @> ?::jsonb")
                params.append(json.dumps([tag]))
            else:
                conditions.append("tags LIKE ?")
                params.append(f'%"{tag}"%')

        if search:
            conditions.append("(name LIKE ? OR smiles LIKE ?)")
            params.append(f"%{search}%")
            params.append(f"%{search}%")

        where = " AND ".join(conditions)

        if self._direct_sqlite:
            with self._lock:
                total = self._conn.execute(
                    f"SELECT COUNT(*) FROM library_molecules WHERE {where}", params
                ).fetchone()[0]
                rows = self._conn.execute(
                    f"SELECT * FROM library_molecules WHERE {where} "
                    f"ORDER BY updated_at DESC LIMIT ? OFFSET ?",
                    params + [limit, offset],
                ).fetchall()
            return [self._row_to_dict(dict(r)) for r in rows], total
        else:
            total_rows = self._backend.execute(
                f"SELECT COUNT(*) as cnt FROM library_molecules WHERE {where}",
                tuple(params),
            )
            total = total_rows[0]["cnt"]
            rows = self._backend.execute(
                f"SELECT * FROM library_molecules WHERE {where} "
                f"ORDER BY updated_at DESC LIMIT ? OFFSET ?",
                tuple(params + [limit, offset]),
            )
            return [self._row_to_dict(r) for r in rows], total

    def update_molecule(
        self,
        mol_id: int,
        user_email: str,
        name: str | None = None,
        tags: list[str] | None = None,
        notes: str | None = None,
    ) -> bool:
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        updates = ["updated_at = ?"]
        params: list = [now]
        if name is not None:
            updates.append("name = ?")
            params.append(name)
        if tags is not None:
            updates.append("tags = ?")
            params.append(self._encode_json(tags))
        if notes is not None:
            updates.append("notes = ?")
            params.append(notes)

        params.extend([mol_id, user_email])

        if self._direct_sqlite:
            with self._lock:
                cursor = self._conn.execute(
                    f"UPDATE library_molecules SET {', '.join(updates)} "
                    f"WHERE id = ? AND user_email = ?",
                    params,
                )
                self._conn.commit()
            return cursor.rowcount > 0
        else:
            affected = self._backend.execute_update(
                f"UPDATE library_molecules SET {', '.join(updates)} "
                f"WHERE id = ? AND user_email = ?",
                tuple(params),
            )
            return affected > 0

    def delete_molecule(self, mol_id: int, user_email: str) -> bool:
        if self._direct_sqlite:
            with self._lock:
                cursor = self._conn.execute(
                    "DELETE FROM library_molecules WHERE id = ? AND user_email = ?",
                    (mol_id, user_email),
                )
                self._conn.commit()
            return cursor.rowcount > 0
        else:
            affected = self._backend.execute_update(
                "DELETE FROM library_molecules WHERE id = ? AND user_email = ?",
                (mol_id, user_email),
            )
            return affected > 0

    def count_molecules(self, user_email: str) -> int:
        if self._direct_sqlite:
            with self._lock:
                return self._conn.execute(
                    "SELECT COUNT(*) FROM library_molecules WHERE user_email = ?",
                    (user_email,),
                ).fetchone()[0]
        else:
            rows = self._backend.execute(
                "SELECT COUNT(*) as cnt FROM library_molecules WHERE user_email = ?",
                (user_email,),
            )
            return rows[0]["cnt"]

    def _row_to_dict(self, d: dict) -> dict:
        d["tags"] = self._decode_json(d.get("tags"), default=[])
        props = d.pop("properties_json", None)
        d["properties"] = self._decode_json(props, default={})
        return d

    def close(self) -> None:
        if self._direct_sqlite:
            try:
                self._conn.close()
            except Exception:
                pass
        else:
            self._backend.close()


_library_db: LibraryDB | None = None


def get_library_db() -> LibraryDB:
    global _library_db
    if _library_db is None:
        from app.config import settings
        if settings.database_backend == "postgresql":
            from app.services.database import get_backend
            _library_db = LibraryDB(backend=get_backend())
        else:
            _library_db = LibraryDB(settings.library_db_path)
    return _library_db


def set_library_db(db: LibraryDB | None) -> None:
    global _library_db
    _library_db = db
