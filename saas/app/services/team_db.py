"""Persistence for teams, members, and team molecule libraries."""

from __future__ import annotations

import json
import sqlite3
import threading
import time

from app.services.database import DatabaseBackend, DatabaseIntegrityError


def _escape_like(s: str) -> str:
    """Escape LIKE wildcard characters in user input."""
    return s.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")


class TeamDB:
    """Thread-safe database for team management.

    SQLite mode uses a single shared connection with a threading lock.
    PostgreSQL mode delegates concurrency to the connection pool.
    """

    def __init__(
        self,
        db_path: str = "molbuilder_teams.db",
        backend: DatabaseBackend | None = None,
    ):
        if backend is not None:
            self._backend = backend
            self._direct_sqlite = False
            self._lock: threading.Lock | None = None
        else:
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
            from app.services.database import SQLiteBackend
            self._backend = SQLiteBackend(db_path)

    def _init_sqlite(self) -> None:
        self._conn.executescript("""
            CREATE TABLE IF NOT EXISTS teams (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT NOT NULL,
                slug TEXT NOT NULL UNIQUE,
                owner_email TEXT NOT NULL,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL
            );
            CREATE INDEX IF NOT EXISTS idx_teams_owner
                ON teams(owner_email);

            CREATE TABLE IF NOT EXISTS team_members (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                team_id INTEGER NOT NULL REFERENCES teams(id) ON DELETE CASCADE,
                user_email TEXT NOT NULL,
                team_role TEXT NOT NULL DEFAULT 'member',
                joined_at TEXT NOT NULL,
                UNIQUE(team_id, user_email)
            );

            CREATE TABLE IF NOT EXISTS team_library_molecules (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                team_id INTEGER NOT NULL REFERENCES teams(id) ON DELETE CASCADE,
                smiles TEXT NOT NULL,
                name TEXT,
                tags TEXT NOT NULL DEFAULT '[]',
                notes TEXT,
                properties_json TEXT,
                added_by TEXT NOT NULL,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL,
                UNIQUE(team_id, smiles)
            );
        """)
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
    # Timestamp helper
    # ------------------------------------------------------------------ #

    @staticmethod
    def _now() -> str:
        return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    # ================================================================== #
    # Team CRUD
    # ================================================================== #

    def create_team(self, name: str, slug: str, owner_email: str) -> int:
        """Create a team and add the owner as a member. Returns team ID."""
        now = self._now()
        if self._direct_sqlite:
            with self._lock:
                try:
                    cursor = self._conn.execute(
                        "INSERT INTO teams (name, slug, owner_email, created_at, updated_at) "
                        "VALUES (?, ?, ?, ?, ?)",
                        (name, slug, owner_email, now, now),
                    )
                except sqlite3.IntegrityError as exc:
                    raise DatabaseIntegrityError(str(exc)) from exc
                team_id = cursor.lastrowid
                self._conn.execute(
                    "INSERT INTO team_members (team_id, user_email, team_role, joined_at) "
                    "VALUES (?, ?, ?, ?)",
                    (team_id, owner_email, "owner", now),
                )
                self._conn.commit()
            return team_id
        else:
            team_id = self._backend.execute_insert(
                "INSERT INTO teams (name, slug, owner_email, created_at, updated_at) "
                "VALUES (?, ?, ?, ?, ?)",
                (name, slug, owner_email, now, now),
            )
            self._backend.execute_insert(
                "INSERT INTO team_members (team_id, user_email, team_role, joined_at) "
                "VALUES (?, ?, ?, ?)",
                (team_id, owner_email, "owner", now),
            )
            return team_id

    def get_team(self, team_id: int) -> dict | None:
        if self._direct_sqlite:
            with self._lock:
                row = self._conn.execute(
                    "SELECT * FROM teams WHERE id = ?", (team_id,),
                ).fetchone()
            return dict(row) if row else None
        else:
            rows = self._backend.execute(
                "SELECT * FROM teams WHERE id = ?", (team_id,),
            )
            return rows[0] if rows else None

    def get_team_by_slug(self, slug: str) -> dict | None:
        if self._direct_sqlite:
            with self._lock:
                row = self._conn.execute(
                    "SELECT * FROM teams WHERE slug = ?", (slug,),
                ).fetchone()
            return dict(row) if row else None
        else:
            rows = self._backend.execute(
                "SELECT * FROM teams WHERE slug = ?", (slug,),
            )
            return rows[0] if rows else None

    def list_teams_for_user(self, user_email: str) -> list[dict]:
        sql = (
            "SELECT t.* FROM teams t "
            "JOIN team_members m ON t.id = m.team_id "
            "WHERE m.user_email = ? "
            "ORDER BY t.name"
        )
        if self._direct_sqlite:
            with self._lock:
                rows = self._conn.execute(sql, (user_email,)).fetchall()
            return [dict(r) for r in rows]
        else:
            return self._backend.execute(sql, (user_email,))

    def update_team(self, team_id: int, name: str) -> bool:
        now = self._now()
        if self._direct_sqlite:
            with self._lock:
                cursor = self._conn.execute(
                    "UPDATE teams SET name = ?, updated_at = ? WHERE id = ?",
                    (name, now, team_id),
                )
                self._conn.commit()
            return cursor.rowcount > 0
        else:
            affected = self._backend.execute_update(
                "UPDATE teams SET name = ?, updated_at = ? WHERE id = ?",
                (name, now, team_id),
            )
            return affected > 0

    def delete_team(self, team_id: int) -> bool:
        if self._direct_sqlite:
            with self._lock:
                cursor = self._conn.execute(
                    "DELETE FROM teams WHERE id = ?", (team_id,),
                )
                self._conn.commit()
            return cursor.rowcount > 0
        else:
            affected = self._backend.execute_update(
                "DELETE FROM teams WHERE id = ?", (team_id,),
            )
            return affected > 0

    # ================================================================== #
    # Member management
    # ================================================================== #

    def add_member(
        self, team_id: int, user_email: str, team_role: str = "member",
    ) -> int:
        """Add a member to a team. Returns the membership row ID."""
        now = self._now()
        if self._direct_sqlite:
            with self._lock:
                try:
                    cursor = self._conn.execute(
                        "INSERT INTO team_members (team_id, user_email, team_role, joined_at) "
                        "VALUES (?, ?, ?, ?)",
                        (team_id, user_email, team_role, now),
                    )
                except sqlite3.IntegrityError as exc:
                    raise DatabaseIntegrityError(str(exc)) from exc
                self._conn.commit()
                return cursor.lastrowid
        else:
            return self._backend.execute_insert(
                "INSERT INTO team_members (team_id, user_email, team_role, joined_at) "
                "VALUES (?, ?, ?, ?)",
                (team_id, user_email, team_role, now),
            )

    def get_member(self, team_id: int, user_email: str) -> dict | None:
        if self._direct_sqlite:
            with self._lock:
                row = self._conn.execute(
                    "SELECT * FROM team_members WHERE team_id = ? AND user_email = ?",
                    (team_id, user_email),
                ).fetchone()
            return dict(row) if row else None
        else:
            rows = self._backend.execute(
                "SELECT * FROM team_members WHERE team_id = ? AND user_email = ?",
                (team_id, user_email),
            )
            return rows[0] if rows else None

    def list_members(self, team_id: int) -> list[dict]:
        if self._direct_sqlite:
            with self._lock:
                rows = self._conn.execute(
                    "SELECT * FROM team_members WHERE team_id = ? ORDER BY joined_at",
                    (team_id,),
                ).fetchall()
            return [dict(r) for r in rows]
        else:
            return self._backend.execute(
                "SELECT * FROM team_members WHERE team_id = ? ORDER BY joined_at",
                (team_id,),
            )

    def update_member_role(
        self, team_id: int, user_email: str, new_role: str,
    ) -> bool:
        if self._direct_sqlite:
            with self._lock:
                cursor = self._conn.execute(
                    "UPDATE team_members SET team_role = ? "
                    "WHERE team_id = ? AND user_email = ?",
                    (new_role, team_id, user_email),
                )
                self._conn.commit()
            return cursor.rowcount > 0
        else:
            affected = self._backend.execute_update(
                "UPDATE team_members SET team_role = ? "
                "WHERE team_id = ? AND user_email = ?",
                (new_role, team_id, user_email),
            )
            return affected > 0

    def remove_member(self, team_id: int, user_email: str) -> bool:
        if self._direct_sqlite:
            with self._lock:
                cursor = self._conn.execute(
                    "DELETE FROM team_members WHERE team_id = ? AND user_email = ?",
                    (team_id, user_email),
                )
                self._conn.commit()
            return cursor.rowcount > 0
        else:
            affected = self._backend.execute_update(
                "DELETE FROM team_members WHERE team_id = ? AND user_email = ?",
                (team_id, user_email),
            )
            return affected > 0

    # ================================================================== #
    # Team library
    # ================================================================== #

    def save_molecule(
        self,
        team_id: int,
        smiles: str,
        added_by: str,
        name: str | None = None,
        tags: list[str] | None = None,
        notes: str | None = None,
        properties: dict | None = None,
    ) -> int:
        """Save a molecule to the team library. Returns row ID."""
        now = self._now()
        tags_val = self._encode_json(tags or [])
        props_val = self._encode_json(properties or {})

        if self._direct_sqlite:
            with self._lock:
                try:
                    cursor = self._conn.execute(
                        "INSERT INTO team_library_molecules "
                        "(team_id, smiles, name, tags, notes, properties_json, "
                        "added_by, created_at, updated_at) "
                        "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                        (team_id, smiles, name, tags_val, notes, props_val,
                         added_by, now, now),
                    )
                except sqlite3.IntegrityError as exc:
                    raise DatabaseIntegrityError(str(exc)) from exc
                self._conn.commit()
                return cursor.lastrowid
        else:
            return self._backend.execute_insert(
                "INSERT INTO team_library_molecules "
                "(team_id, smiles, name, tags, notes, properties_json, "
                "added_by, created_at, updated_at) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (team_id, smiles, name, tags_val, notes, props_val,
                 added_by, now, now),
            )

    def get_molecule(self, mol_id: int, team_id: int) -> dict | None:
        if self._direct_sqlite:
            with self._lock:
                row = self._conn.execute(
                    "SELECT * FROM team_library_molecules "
                    "WHERE id = ? AND team_id = ?",
                    (mol_id, team_id),
                ).fetchone()
            if row is None:
                return None
            return self._mol_to_dict(dict(row))
        else:
            rows = self._backend.execute(
                "SELECT * FROM team_library_molecules "
                "WHERE id = ? AND team_id = ?",
                (mol_id, team_id),
            )
            if not rows:
                return None
            return self._mol_to_dict(rows[0])

    def list_molecules(
        self,
        team_id: int,
        tag: str | None = None,
        search: str | None = None,
        limit: int = 20,
        offset: int = 0,
    ) -> tuple[list[dict], int]:
        conditions = ["team_id = ?"]
        params: list = [team_id]

        if tag:
            if self._backend.is_postgres:
                conditions.append("tags @> ?::jsonb")
                params.append(json.dumps([tag]))
            else:
                conditions.append("tags LIKE ? ESCAPE '\\'")
                params.append(f'%"{_escape_like(tag)}"%')

        if search:
            escaped = _escape_like(search)
            conditions.append("(name LIKE ? ESCAPE '\\' OR smiles LIKE ? ESCAPE '\\')")
            params.append(f"%{escaped}%")
            params.append(f"%{escaped}%")

        where = " AND ".join(conditions)

        if self._direct_sqlite:
            with self._lock:
                total = self._conn.execute(
                    f"SELECT COUNT(*) FROM team_library_molecules WHERE {where}",
                    params,
                ).fetchone()[0]
                rows = self._conn.execute(
                    f"SELECT * FROM team_library_molecules WHERE {where} "
                    f"ORDER BY updated_at DESC LIMIT ? OFFSET ?",
                    params + [limit, offset],
                ).fetchall()
            return [self._mol_to_dict(dict(r)) for r in rows], total
        else:
            total_rows = self._backend.execute(
                f"SELECT COUNT(*) as cnt FROM team_library_molecules WHERE {where}",
                tuple(params),
            )
            total = total_rows[0]["cnt"]
            rows = self._backend.execute(
                f"SELECT * FROM team_library_molecules WHERE {where} "
                f"ORDER BY updated_at DESC LIMIT ? OFFSET ?",
                tuple(params + [limit, offset]),
            )
            return [self._mol_to_dict(r) for r in rows], total

    def update_molecule(
        self,
        mol_id: int,
        team_id: int,
        name: str | None = None,
        tags: list[str] | None = None,
        notes: str | None = None,
    ) -> bool:
        now = self._now()
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

        params.extend([mol_id, team_id])

        if self._direct_sqlite:
            with self._lock:
                cursor = self._conn.execute(
                    f"UPDATE team_library_molecules SET {', '.join(updates)} "
                    f"WHERE id = ? AND team_id = ?",
                    params,
                )
                self._conn.commit()
            return cursor.rowcount > 0
        else:
            affected = self._backend.execute_update(
                f"UPDATE team_library_molecules SET {', '.join(updates)} "
                f"WHERE id = ? AND team_id = ?",
                tuple(params),
            )
            return affected > 0

    def delete_molecule(self, mol_id: int, team_id: int) -> bool:
        if self._direct_sqlite:
            with self._lock:
                cursor = self._conn.execute(
                    "DELETE FROM team_library_molecules "
                    "WHERE id = ? AND team_id = ?",
                    (mol_id, team_id),
                )
                self._conn.commit()
            return cursor.rowcount > 0
        else:
            affected = self._backend.execute_update(
                "DELETE FROM team_library_molecules "
                "WHERE id = ? AND team_id = ?",
                (mol_id, team_id),
            )
            return affected > 0

    def _mol_to_dict(self, d: dict) -> dict:
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
        # When using a shared PG backend, don't close the pool -- it's
        # managed by the application lifecycle, not individual DB classes.


_team_db: TeamDB | None = None


def get_team_db() -> TeamDB:
    global _team_db
    if _team_db is None:
        from app.config import settings
        if settings.database_backend == "postgresql":
            from app.services.database import get_backend
            _team_db = TeamDB(backend=get_backend())
        else:
            _team_db = TeamDB(settings.team_db_path)
    return _team_db


def set_team_db(db: TeamDB | None) -> None:
    global _team_db
    _team_db = db
