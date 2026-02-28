"""
molbuilder -- Professional-grade molecular engineering toolkit.

Build, analyze, and plan synthesis of molecules from atoms up through
industrial-scale manufacturing processes.

Quick start::

    from molbuilder import Molecule, parse, to_smiles

    mol = parse("CCO")  # ethanol
    print(f"{len(mol.atoms)} atoms, {len(mol.bonds)} bonds")
    print(to_smiles(mol))
"""

__version__ = "1.2.0"

# Core data types
from molbuilder.molecule.graph import Molecule, Atom, Bond, Hybridization

# SMILES I/O
from molbuilder.smiles.parser import parse
from molbuilder.smiles.writer import to_smiles

__all__ = [
    "__version__",
    "Molecule",
    "Atom",
    "Bond",
    "Hybridization",
    "parse",
    "to_smiles",
]
