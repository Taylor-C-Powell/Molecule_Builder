"""XYZ file format reader/writer.

The XYZ format is one of the simplest molecular geometry formats::

    <atom_count>
    <comment line>
    <symbol> <x> <y> <z>
    <symbol> <x> <y> <z>
    ...

XYZ files do not store bond connectivity, so bonds are inferred from
interatomic distances and covalent radii when reading.
"""

from __future__ import annotations


import numpy as np

from molbuilder.molecule.graph import Molecule
from molbuilder.core.element_properties import covalent_radius_pm


# ── Bond inference ────────────────────────────────────────────────────

_BOND_TOLERANCE = 1.3  # multiplier on sum-of-covalent-radii


def _infer_bonds(mol: Molecule) -> None:
    """Add bonds between atoms whose distance is within tolerance of
    the sum of their covalent radii (converted from pm to Angstroms).
    """
    n = len(mol.atoms)
    for i in range(n):
        for j in range(i + 1, n):
            ri = covalent_radius_pm(mol.atoms[i].symbol) / 100.0
            rj = covalent_radius_pm(mol.atoms[j].symbol) / 100.0
            max_dist = _BOND_TOLERANCE * (ri + rj)
            dist = float(np.linalg.norm(
                mol.atoms[i].position - mol.atoms[j].position))
            if dist < max_dist:
                mol.add_bond(i, j, order=1, rotatable=True)


# ── String serialisation ─────────────────────────────────────────────

def to_xyz_string(mol: Molecule) -> str:
    """Return the molecule as an XYZ-format string."""
    lines: list[str] = []
    lines.append(str(len(mol.atoms)))
    lines.append(mol.name if mol.name else "")
    for atom in mol.atoms:
        x, y, z = atom.position
        lines.append(f"{atom.symbol:<4s} {x:15.8f} {y:15.8f} {z:15.8f}")
    return "\n".join(lines) + "\n"


def from_xyz_string(content: str) -> Molecule:
    """Parse a Molecule from an XYZ-format string.

    Bonds are inferred from interatomic distances using covalent radii.
    """
    lines = content.strip().splitlines()
    if len(lines) < 2:
        raise ValueError("XYZ content must have at least two lines "
                         "(atom count and comment).")

    atom_count = int(lines[0].strip())
    comment = lines[1].strip()

    if len(lines) < atom_count + 2:
        raise ValueError(f"XYZ file declares {atom_count} atoms but has only "
                         f"{len(lines) - 2} atom lines")

    mol = Molecule(name=comment)

    for i in range(atom_count):
        parts = lines[2 + i].split()
        symbol = parts[0]
        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        mol.add_atom(symbol, [x, y, z])

    _infer_bonds(mol)
    return mol


# ── File I/O ──────────────────────────────────────────────────────────

def write_xyz(mol: Molecule, filepath: str) -> None:
    """Write a Molecule to an XYZ file."""
    with open(filepath, "w") as f:
        f.write(to_xyz_string(mol))


def read_xyz(filepath: str) -> Molecule:
    """Read a Molecule from an XYZ file.

    Bonds are inferred from interatomic distances using covalent radii.
    """
    with open(filepath, "r") as f:
        content = f.read()
    return from_xyz_string(content)
