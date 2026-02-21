"""3D coordinate generation for molecules.

Public API
----------
generate_3d(mol, backend="auto", seed=42)
    Generate 3D coordinates for a Molecule, modifying positions in-place.

Backends
--------
- ``"auto"``: try RDKit, fall back to builtin
- ``"rdkit"``: require RDKit (ImportError if missing)
- ``"builtin"``: pure Python distance geometry + force field pipeline
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from molbuilder.molecule.graph import Molecule


def generate_3d(mol: Molecule, backend: str = "auto", seed: int = 42) -> None:
    """Generate 3D coordinates for mol, modifying positions in-place.

    Parameters
    ----------
    mol : Molecule
        Molecule with atoms and bonds defined. Positions will be overwritten.
    backend : str
        ``"auto"`` tries RDKit first, then builtin.
        ``"rdkit"`` requires RDKit (raises ImportError if missing).
        ``"builtin"`` uses pure-Python DG + FF pipeline.
    seed : int
        Random seed for reproducibility.

    Raises
    ------
    ImportError
        If backend="rdkit" and rdkit is not installed.
    ValueError
        If backend is not recognized.
    """
    if backend == "rdkit":
        from molbuilder.coords.rdkit_backend import generate_rdkit
        generate_rdkit(mol, seed=seed)
        return

    if backend == "builtin":
        from molbuilder.coords.builtin import generate_builtin
        generate_builtin(mol, seed=seed)
        return

    if backend == "auto":
        try:
            from molbuilder.coords.rdkit_backend import generate_rdkit
            generate_rdkit(mol, seed=seed)
        except ImportError:
            from molbuilder.coords.builtin import generate_builtin
            generate_builtin(mol, seed=seed)
        return

    raise ValueError(
        f"Unknown backend {backend!r}. Use 'auto', 'builtin', or 'rdkit'."
    )
