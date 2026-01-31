"""V2000 MOL/SDF file format reader/writer.

MOL file layout (V2000)::

    <molecule name>
      molbuilder          3D

     <atom_count> <bond_count>  0  0  0  0  0  0  0  0999 V2000
    <x10.4><y10.4><z10.4> <sym3> 0  0  0  0  0  0  0  0  0  0  0  0
    ...
    <i3><j3><type3>  0  0  0  0
    ...
    M  END

SDF files contain one or more MOL blocks separated by ``$$$$``.

Note: MOL files use **1-based** atom indices.
Bond type: 1 = single, 2 = double, 3 = triple.
"""

from __future__ import annotations

import numpy as np

from molbuilder.molecule.graph import Molecule, Hybridization


# ── MOL string serialisation ─────────────────────────────────────────

def to_mol_string(mol: Molecule) -> str:
    """Serialise a Molecule to a V2000 MOL block string."""
    lines: list[str] = []

    # Header block (3 lines)
    lines.append(mol.name if mol.name else "")
    lines.append("  molbuilder          3D")
    lines.append("")

    # Counts line
    n_atoms = len(mol.atoms)
    n_bonds = len(mol.bonds)
    lines.append(
        f"{n_atoms:3d}{n_bonds:3d}"
        f"  0  0  0  0  0  0  0  0999 V2000"
    )

    # Atom block
    for atom in mol.atoms:
        x, y, z = atom.position
        symbol = atom.symbol
        lines.append(
            f"{x:10.4f}{y:10.4f}{z:10.4f} {symbol:<3s} 0  0  0  0  0  0  0  0  0  0  0  0"
        )

    # Bond block
    for bond in mol.bonds:
        i = bond.atom_i + 1  # 1-based
        j = bond.atom_j + 1
        order = bond.order
        lines.append(f"{i:3d}{j:3d}{order:3d}  0  0  0  0")

    lines.append("M  END")
    return "\n".join(lines) + "\n"


def from_mol_string(content: str) -> Molecule:
    """Parse a Molecule from a V2000 MOL block string."""
    lines = content.splitlines()

    # Header
    name = lines[0].strip() if lines[0].strip() else ""
    # lines[1] = program/timestamp, lines[2] = comment (both ignored)

    # Counts line
    counts_line = lines[3]
    n_atoms = int(counts_line[0:3])
    n_bonds = int(counts_line[3:6])

    mol = Molecule(name=name)

    # Atom block: starts at line 4
    for i in range(n_atoms):
        line = lines[4 + i]
        x = float(line[0:10])
        y = float(line[10:20])
        z = float(line[20:30])
        symbol = line[31:34].strip()
        mol.add_atom(symbol, [x, y, z])

    # Bond block: starts after atom block
    bond_start = 4 + n_atoms
    for i in range(n_bonds):
        line = lines[bond_start + i]
        ai = int(line[0:3]) - 1   # convert to 0-based
        aj = int(line[3:6]) - 1
        order = int(line[6:9])
        rotatable = (order == 1)
        mol.add_bond(ai, aj, order=order, rotatable=rotatable)

    return mol


# ── MOL file I/O ─────────────────────────────────────────────────────

def write_mol(mol: Molecule, filepath: str) -> None:
    """Write a Molecule to a V2000 MOL file."""
    with open(filepath, "w") as f:
        f.write(to_mol_string(mol))


def read_mol(filepath: str) -> Molecule:
    """Read a Molecule from a V2000 MOL file."""
    with open(filepath, "r") as f:
        content = f.read()
    return from_mol_string(content)


# ── SDF multi-molecule I/O ───────────────────────────────────────────

def write_sdf(molecules: list[Molecule], filepath: str) -> None:
    """Write multiple Molecules to an SDF file.

    Each MOL block is followed by ``$$$$`` as the record separator.
    """
    with open(filepath, "w") as f:
        for mol in molecules:
            f.write(to_mol_string(mol))
            f.write("$$$$\n")


def read_sdf(filepath: str) -> list[Molecule]:
    """Read all Molecules from an SDF file.

    SDF files contain one or more MOL blocks separated by ``$$$$``.
    Data items between ``M  END`` and ``$$$$`` are silently ignored.
    """
    with open(filepath, "r") as f:
        content = f.read()

    molecules: list[Molecule] = []
    blocks = content.split("$$$$")

    for block in blocks:
        block = block.strip()
        if not block:
            continue

        # Ensure the block contains a valid MOL section
        if "V2000" not in block:
            continue

        # Trim anything after "M  END" (SDF data items)
        end_idx = block.find("M  END")
        if end_idx != -1:
            mol_block = block[:end_idx + len("M  END")]
        else:
            mol_block = block

        try:
            mol = from_mol_string(mol_block)
            molecules.append(mol)
        except (ValueError, IndexError):
            # Skip malformed blocks rather than crashing
            continue

    return molecules
