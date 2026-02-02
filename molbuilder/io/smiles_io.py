"""SMILES file I/O wrapper.

Reads and writes SMILES files where each line contains a SMILES string
optionally followed by a molecule name::

    CCO ethanol
    c1ccccc1 benzene
"""

import warnings

from molbuilder.smiles.parser import parse
from molbuilder.smiles.writer import to_smiles
from molbuilder.molecule.graph import Molecule


def write_smiles(mol: Molecule, filepath: str) -> None:
    """Write a SMILES string to a file."""
    smi = to_smiles(mol)
    with open(filepath, "w") as f:
        f.write(f"{smi} {mol.name}\n")


def read_smiles(filepath: str) -> list[Molecule]:
    """Read molecules from a SMILES file (one per line).

    Blank lines and lines starting with ``#`` are skipped.  Each line
    may contain a SMILES string followed by an optional name separated
    by whitespace.
    """
    molecules = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            smi = parts[0]
            name = parts[1] if len(parts) > 1 else smi
            try:
                mol = parse(smi)
            except (ValueError, IndexError) as e:
                warnings.warn(f"Skipping invalid SMILES '{smi}': {e}")
                continue
            mol.name = name
            molecules.append(mol)
    return molecules
