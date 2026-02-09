"""Thread-safe in-memory molecule store with LRU eviction."""

import threading
import uuid
from collections import OrderedDict
from molbuilder.molecule.graph import Molecule
from app.config import settings


class MoleculeStore:
    """In-memory molecule storage with LRU eviction at capacity."""

    def __init__(self, max_size: int | None = None) -> None:
        self._max = max_size or settings.molecule_store_max
        self._lock = threading.Lock()
        self._store: OrderedDict[str, tuple[Molecule, str]] = OrderedDict()

    def put(self, mol: Molecule, smiles: str) -> str:
        mol_id = uuid.uuid4().hex[:12]
        with self._lock:
            if mol_id in self._store:
                self._store.move_to_end(mol_id)
            self._store[mol_id] = (mol, smiles)
            while len(self._store) > self._max:
                self._store.popitem(last=False)
        return mol_id

    def get(self, mol_id: str) -> tuple[Molecule, str] | None:
        with self._lock:
            item = self._store.get(mol_id)
            if item is not None:
                self._store.move_to_end(mol_id)
            return item


molecule_store = MoleculeStore()
