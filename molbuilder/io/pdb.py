"""PDB file format reader/writer.

Supports ATOM/HETATM records for coordinates and CONECT records for
bond connectivity.  Intended for small-molecule use (not full protein
PDB support).

Record formats used::

    ATOM      1  C1  MOL A   1       0.000   0.000   0.000  1.00  0.00           C
    CONECT    1    2    3    4    5
    END
"""

from __future__ import annotations

from collections import Counter, defaultdict

import numpy as np

from molbuilder.molecule.graph import Molecule


# ── Helpers ───────────────────────────────────────────────────────────

def _atom_names(mol: Molecule) -> list[str]:
    """Generate PDB atom names (e.g. C1, C2, H1, H2, ...) for each atom
    based on its element and order of occurrence within that element."""
    counts: Counter[str] = Counter()
    names: list[str] = []
    for atom in mol.atoms:
        counts[atom.symbol] += 1
        names.append(f"{atom.symbol}{counts[atom.symbol]}")
    return names


def _residue_name(mol: Molecule) -> str:
    """Return a 3-character residue name derived from the molecule name."""
    name = mol.name.strip() if mol.name else "MOL"
    if len(name) == 0:
        return "MOL"
    # Uppercase, first 3 characters
    return name[:3].upper().ljust(3)


# ── String serialisation ─────────────────────────────────────────────

def to_pdb_string(mol: Molecule) -> str:
    """Serialise a Molecule to PDB-format string with ATOM and CONECT
    records."""
    lines: list[str] = []
    atom_names = _atom_names(mol)
    res_name = _residue_name(mol)
    chain = "A"
    res_seq = 1

    # ATOM records
    for atom in mol.atoms:
        serial = atom.index + 1  # 1-based
        aname = atom_names[atom.index]
        x, y, z = atom.position
        element = atom.symbol.rjust(2)

        # PDB ATOM record (fixed-width columns)
        #  1- 6  Record type
        #  7-11  Serial
        # 13-16  Atom name
        # 17     Alternate location
        # 18-20  Residue name
        # 22     Chain ID
        # 23-26  Residue sequence number
        # 31-38  x (8.3f)
        # 39-46  y (8.3f)
        # 47-54  z (8.3f)
        # 55-60  Occupancy (6.2f)
        # 61-66  Temp factor (6.2f)
        # 77-78  Element symbol
        line = (
            f"HETATM{serial:5d} {aname:<4s} {res_name:3s} {chain:1s}"
            f"{res_seq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
            f"{1.0:6.2f}{0.0:6.2f}"
            f"          {element:>2s}"
        )
        lines.append(line)

    # CONECT records
    # Build adjacency: for each atom, list of bonded atom serials (1-based)
    adj: dict[int, list[int]] = defaultdict(list)
    for bond in mol.bonds:
        adj[bond.atom_i + 1].append(bond.atom_j + 1)
        adj[bond.atom_j + 1].append(bond.atom_i + 1)

    for serial in sorted(adj):
        neighbours = sorted(adj[serial])
        # PDB CONECT records can hold up to 4 bonded atoms per line
        for chunk_start in range(0, len(neighbours), 4):
            chunk = neighbours[chunk_start:chunk_start + 4]
            parts = "".join(f"{n:5d}" for n in chunk)
            lines.append(f"CONECT{serial:5d}{parts}")

    lines.append("END")
    return "\n".join(lines) + "\n"


def from_pdb_string(content: str) -> Molecule:
    """Parse a Molecule from PDB-format string.

    Reads HETATM / ATOM records for coordinates and CONECT records for
    bond connectivity.
    """
    mol = Molecule()

    serial_to_index: dict[int, int] = {}
    conect_records: list[tuple[int, list[int]]] = []

    for line in content.splitlines():
        record = line[:6].strip()

        if record in ("ATOM", "HETATM"):
            serial = int(line[6:11])
            # Element symbol: columns 77-78 (preferred), fallback to atom name
            element = line[76:78].strip() if len(line) >= 78 else ""
            if not element:
                # Fallback: strip digits from atom name (cols 12-16)
                raw_name = line[12:16].strip()
                element = "".join(c for c in raw_name if c.isalpha())

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            idx = mol.add_atom(element, [x, y, z])
            serial_to_index[serial] = idx

            # Extract molecule name from residue name on first atom
            if idx == 0:
                res = line[17:20].strip()
                mol.name = res

        elif record == "CONECT":
            serial = int(line[6:11])
            neighbours: list[int] = []
            col = 11
            while col + 5 <= len(line):
                token = line[col:col + 5].strip()
                if token:
                    try:
                        neighbours.append(int(token))
                    except ValueError:
                        pass
                col += 5
            conect_records.append((serial, neighbours))

    # Build bonds from CONECT records (avoid duplicates)
    added_bonds: set[tuple[int, int]] = set()
    for serial, neighbours in conect_records:
        if serial not in serial_to_index:
            continue
        i = serial_to_index[serial]
        for nb_serial in neighbours:
            if nb_serial not in serial_to_index:
                continue
            j = serial_to_index[nb_serial]
            bond_key = (min(i, j), max(i, j))
            if bond_key not in added_bonds:
                added_bonds.add(bond_key)
                mol.add_bond(i, j, order=1, rotatable=True)

    return mol


# ── File I/O ──────────────────────────────────────────────────────────

def write_pdb(mol: Molecule, filepath: str) -> None:
    """Write a Molecule to a PDB file."""
    with open(filepath, "w") as f:
        f.write(to_pdb_string(mol))


def read_pdb(filepath: str) -> Molecule:
    """Read a Molecule from a PDB file."""
    with open(filepath, "r") as f:
        content = f.read()
    return from_pdb_string(content)
