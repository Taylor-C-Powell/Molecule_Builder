"""Molecule store with LRU eviction -- SQLite or PostgreSQL.

Stores SMILES strings and reconstructs Molecule objects on get().
Falls back to in-memory SQLite when molecule_db_path is empty or ':memory:'.
"""

from __future__ import annotations

import sqlite3
import threading
import time
import uuid

from molbuilder.molecule.graph import Molecule
from molbuilder.smiles.parser import parse

from app.services.database import DatabaseBackend


class MoleculeStore:
    """Molecule storage with LRU eviction at capacity."""

    def __init__(
        self,
        db_path: str = ":memory:",
        max_size: int | None = None,
        backend: DatabaseBackend | None = None,
    ) -> None:
        from app.config import settings

        self._max = max_size or settings.molecule_store_max

        if backend is not None:
            self._backend = backend
            self._direct_sqlite = False
        else:
            self._db_path = db_path if db_path else ":memory:"
            self._direct_sqlite = True
            self._local = threading.local()
            self._init_sqlite()
            from app.services.database import SQLiteBackend
            self._backend = SQLiteBackend(self._db_path)

    # ------------------------------------------------------------------ #
    # SQLite-only bootstrap
    # ------------------------------------------------------------------ #

    def _get_conn(self) -> sqlite3.Connection:
        if not hasattr(self._local, "conn") or self._local.conn is None:
            conn = sqlite3.connect(self._db_path)
            conn.execute("PRAGMA journal_mode=WAL")
            conn.row_factory = sqlite3.Row
            self._local.conn = conn
        return self._local.conn

    def _init_sqlite(self) -> None:
        conn = self._get_conn()
        conn.execute("""
            CREATE TABLE IF NOT EXISTS molecules (
                mol_id TEXT PRIMARY KEY,
                smiles TEXT NOT NULL,
                created_at REAL NOT NULL,
                accessed_at REAL NOT NULL
            )
        """)
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_molecules_accessed "
            "ON molecules(accessed_at)"
        )
        conn.commit()

    # ------------------------------------------------------------------ #
    # Public API
    # ------------------------------------------------------------------ #

    def put(self, mol: Molecule, smiles: str) -> str:
        """Store a molecule's SMILES. Returns 12-char hex ID."""
        mol_id = uuid.uuid4().hex[:12]
        now = time.time()

        if self._direct_sqlite:
            conn = self._get_conn()
            conn.execute(
                "INSERT OR REPLACE INTO molecules (mol_id, smiles, created_at, accessed_at) "
                "VALUES (?, ?, ?, ?)",
                (mol_id, smiles, now, now),
            )
            count = conn.execute("SELECT COUNT(*) FROM molecules").fetchone()[0]
            if count > self._max:
                excess = count - self._max
                conn.execute(
                    "DELETE FROM molecules WHERE mol_id IN ("
                    "  SELECT mol_id FROM molecules ORDER BY accessed_at ASC LIMIT ?"
                    ")",
                    (excess,),
                )
            conn.commit()
        else:
            # PostgreSQL upsert
            self._backend.execute_update(
                "INSERT INTO molecules (mol_id, smiles, created_at, accessed_at) "
                "VALUES (?, ?, ?, ?) "
                "ON CONFLICT (mol_id) DO UPDATE "
                "SET smiles = EXCLUDED.smiles, accessed_at = EXCLUDED.accessed_at",
                (mol_id, smiles, now, now),
            )
            # LRU eviction
            count_rows = self._backend.execute(
                "SELECT COUNT(*) as cnt FROM molecules"
            )
            count = count_rows[0]["cnt"]
            if count > self._max:
                excess = count - self._max
                self._backend.execute_update(
                    "DELETE FROM molecules WHERE mol_id IN ("
                    "  SELECT mol_id FROM molecules ORDER BY accessed_at ASC LIMIT ?"
                    ")",
                    (excess,),
                )

        return mol_id

    def get(self, mol_id: str) -> tuple[Molecule, str] | None:
        """Retrieve a molecule by ID. Updates LRU timestamp. Returns None if missing."""
        if self._direct_sqlite:
            conn = self._get_conn()
            row = conn.execute(
                "SELECT smiles FROM molecules WHERE mol_id = ?", (mol_id,)
            ).fetchone()
            if row is None:
                return None
            conn.execute(
                "UPDATE molecules SET accessed_at = ? WHERE mol_id = ?",
                (time.time(), mol_id),
            )
            conn.commit()
            smiles = row["smiles"]
        else:
            rows = self._backend.execute(
                "SELECT smiles FROM molecules WHERE mol_id = ?", (mol_id,)
            )
            if not rows:
                return None
            self._backend.execute_update(
                "UPDATE molecules SET accessed_at = ? WHERE mol_id = ?",
                (time.time(), mol_id),
            )
            smiles = rows[0]["smiles"]

        mol = parse(smiles)
        return mol, smiles

    def close(self) -> None:
        if self._direct_sqlite:
            if hasattr(self, "_local") and hasattr(self._local, "conn") and self._local.conn is not None:
                self._local.conn.close()
                self._local.conn = None
        else:
            self._backend.close()


# Module-level singleton (overridable in tests)
_molecule_store: MoleculeStore | None = None


def get_molecule_store() -> MoleculeStore:
    global _molecule_store
    if _molecule_store is None:
        from app.config import settings

        if settings.database_backend == "postgresql":
            from app.services.database import get_backend
            _molecule_store = MoleculeStore(
                max_size=settings.molecule_store_max,
                backend=get_backend(),
            )
        else:
            _molecule_store = MoleculeStore(
                db_path=settings.molecule_db_path,
                max_size=settings.molecule_store_max,
            )
    return _molecule_store


def set_molecule_store(store: MoleculeStore | None) -> None:
    global _molecule_store
    _molecule_store = store


class _StoreProxy:
    """Lazy proxy so ``molecule_store.put(...)`` works as a drop-in replacement."""

    def __getattr__(self, name: str):
        return getattr(get_molecule_store(), name)


molecule_store: MoleculeStore = _StoreProxy()  # type: ignore[assignment]
