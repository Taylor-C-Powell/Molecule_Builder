"""JSON serialisation / deserialisation for Molecule objects.

Schema::

    {
      "name": "ethanol",
      "atoms": [
        {"index": 0, "symbol": "C", "position": [0.0, 0.0, 0.0],
         "hybridization": "SP3"},
        ...
      ],
      "bonds": [
        {"atom_i": 0, "atom_j": 1, "order": 1, "rotatable": true},
        ...
      ],
      "properties": {
        "formula": "C2H6O",
        "atom_count": 9,
        "bond_count": 8
      }
    }
"""

from __future__ import annotations

import json
from collections import Counter

from molbuilder.molecule.graph import Molecule, Hybridization


# ── Helpers ───────────────────────────────────────────────────────────

# Map for serialising Hybridization to a string and back.
_HYB_TO_STR = {
    Hybridization.SP3: "SP3",
    Hybridization.SP2: "SP2",
    Hybridization.SP:  "SP",
}
_STR_TO_HYB = {v: k for k, v in _HYB_TO_STR.items()}


def _molecular_formula(mol: Molecule) -> str:
    """Return a Hill-system molecular formula (C first, H second, then
    alphabetical)."""
    counts: Counter[str] = Counter()
    for atom in mol.atoms:
        counts[atom.symbol] += 1

    parts: list[str] = []
    for sym in ("C", "H"):
        if sym in counts:
            n = counts.pop(sym)
            parts.append(sym if n == 1 else f"{sym}{n}")
    for sym in sorted(counts):
        n = counts[sym]
        parts.append(sym if n == 1 else f"{sym}{n}")
    return "".join(parts)


# ── Molecule -> dict -> JSON ──────────────────────────────────────────

def _mol_to_dict(mol: Molecule) -> dict:
    """Convert a Molecule to a JSON-serialisable dictionary."""
    atoms = []
    for atom in mol.atoms:
        hyb_str = _HYB_TO_STR.get(atom.hybridization) if atom.hybridization else None
        atoms.append({
            "index": atom.index,
            "symbol": atom.symbol,
            "position": atom.position.tolist(),
            "hybridization": hyb_str,
        })

    bonds = []
    for bond in mol.bonds:
        bonds.append({
            "atom_i": bond.atom_i,
            "atom_j": bond.atom_j,
            "order": bond.order,
            "rotatable": bond.rotatable,
        })

    return {
        "name": mol.name,
        "atoms": atoms,
        "bonds": bonds,
        "properties": {
            "formula": _molecular_formula(mol),
            "atom_count": len(mol.atoms),
            "bond_count": len(mol.bonds),
        },
    }


def _dict_to_mol(data: dict) -> Molecule:
    """Reconstruct a Molecule from a JSON-derived dictionary."""
    mol = Molecule(name=data.get("name", ""))

    for atom_data in data["atoms"]:
        hyb = None
        hyb_str = atom_data.get("hybridization")
        if hyb_str is not None:
            hyb = _STR_TO_HYB.get(hyb_str)
        mol.add_atom(
            symbol=atom_data["symbol"],
            position=atom_data["position"],
            hybridization=hyb,
        )

    for bond_data in data["bonds"]:
        mol.add_bond(
            i=bond_data["atom_i"],
            j=bond_data["atom_j"],
            order=bond_data.get("order", 1),
            rotatable=bond_data.get("rotatable", True),
        )

    return mol


# ── String serialisation ─────────────────────────────────────────────

def to_json_string(mol: Molecule, indent: int = 2) -> str:
    """Serialise a Molecule to a JSON string."""
    return json.dumps(_mol_to_dict(mol), indent=indent)


def from_json_string(content: str) -> Molecule:
    """Deserialise a Molecule from a JSON string."""
    data = json.loads(content)
    return _dict_to_mol(data)


# ── File I/O ──────────────────────────────────────────────────────────

def write_json(mol: Molecule, filepath: str) -> None:
    """Write a Molecule to a JSON file."""
    with open(filepath, "w") as f:
        f.write(to_json_string(mol))

def read_json(filepath: str) -> Molecule:
    """Read a Molecule from a JSON file."""
    with open(filepath, "r") as f:
        content = f.read()
    return from_json_string(content)
