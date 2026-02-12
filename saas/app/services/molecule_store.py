"""SQLite-backed molecule store with LRU eviction.

Stores SMILES strings in SQLite and reconstructs Molecule objects on get().
Falls back to in-memory SQLite when molecule_db_path is empty or ':memory:'.
"""

import sqlite3
import threading
import time
import uuid

from molbuilder.molecule.graph import Molecule
from molbuilder.smiles.parser import parse


class MoleculeStore:
    """SQLite molecule storage with LRU eviction at capacity."""

    def __init__(self, db_path: str = ":memory:", max_size: int | None = None) -> None:
        from app.config import settings

        self._db_path = db_path if db_path else ":memory:"
        self._max = max_size or settings.molecule_store_max
        self._local = threading.local()
        self._init_db()

    def _get_conn(self) -> sqlite3.Connection:
        if not hasattr(self._local, "conn") or self._local.conn is None:
            conn = sqlite3.connect(self._db_path)
            conn.execute("PRAGMA journal_mode=WAL")
            conn.row_factory = sqlite3.Row
            self._local.conn = conn
        return self._local.conn

    def _init_db(self) -> None:
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

    def put(self, mol: Molecule, smiles: str) -> str:
        """Store a molecule's SMILES. Returns 12-char hex ID."""
        mol_id = uuid.uuid4().hex[:12]
        now = time.time()
        conn = self._get_conn()
        conn.execute(
            "INSERT OR REPLACE INTO molecules (mol_id, smiles, created_at, accessed_at) "
            "VALUES (?, ?, ?, ?)",
            (mol_id, smiles, now, now),
        )
        # LRU eviction: delete oldest entries beyond max_size
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
        return mol_id

    def get(self, mol_id: str) -> tuple[Molecule, str] | None:
        """Retrieve a molecule by ID. Updates LRU timestamp. Returns None if missing."""
        conn = self._get_conn()
        row = conn.execute(
            "SELECT smiles FROM molecules WHERE mol_id = ?", (mol_id,)
        ).fetchone()
        if row is None:
            return None
        # LRU touch
        conn.execute(
            "UPDATE molecules SET accessed_at = ? WHERE mol_id = ?",
            (time.time(), mol_id),
        )
        conn.commit()
        smiles = row["smiles"]
        mol = parse(smiles)
        return mol, smiles

    def close(self) -> None:
        if hasattr(self._local, "conn") and self._local.conn is not None:
            self._local.conn.close()
            self._local.conn = None


# Module-level singleton (overridable in tests)
_molecule_store: MoleculeStore | None = None


def get_molecule_store() -> MoleculeStore:
    global _molecule_store
    if _molecule_store is None:
        from app.config import settings

        _molecule_store = MoleculeStore(
            db_path=settings.molecule_db_path,
            max_size=settings.molecule_store_max,
        )
    return _molecule_store


def set_molecule_store(store: MoleculeStore | None) -> None:
    global _molecule_store
    _molecule_store = store


class _StoreProxy:
    """Lazy proxy so ``molecule_store.put(...)`` works as a drop-in replacement.

    The router imports ``molecule_store`` at module level.  We need the
    singleton to be created lazily (after settings are loaded), so this
    proxy forwards attribute access to :func:`get_molecule_store`.
    """

    def __getattr__(self, name: str):
        return getattr(get_molecule_store(), name)


# Drop-in replacement: ``from app.services.molecule_store import molecule_store``
molecule_store: MoleculeStore = _StoreProxy()  # type: ignore[assignment]
