"""Chirality and E/Z stereo enforcement in 3D coordinates.

After FF optimization, checks that stereochemistry parsed from SMILES
(@/@@ for R/S, /\\ for E/Z) is correctly reflected in atom positions.
If wrong, swaps substituent positions to fix.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from molbuilder.molecule.graph import Molecule


def enforce_stereo(mol: Molecule) -> None:
    """Enforce chirality and E/Z stereochemistry in 3D coordinates.

    Modifies atom positions in-place.
    """
    _enforce_tetrahedral(mol)
    _enforce_ez(mol)


def _enforce_tetrahedral(mol: Molecule) -> None:
    """Enforce @/@@ chirality at tetrahedral centers.

    SMILES convention:
        @ = anticlockwise looking from first neighbor
        @@ = clockwise looking from first neighbor

    We compute the signed tetrahedral volume and swap two substituents
    if the sign doesn't match the expected chirality.
    """
    for atom in mol.atoms:
        if atom.chirality not in ("@", "@@"):
            continue

        nbrs = mol.neighbors(atom.index)
        if len(nbrs) != 4:
            continue

        # Compute signed volume: det([v1-v0, v2-v0, v3-v0])
        # where v0..v3 are positions of the 4 neighbors in connection order
        p = [mol.atoms[n].position for n in nbrs]
        v1 = p[1] - p[0]
        v2 = p[2] - p[0]
        v3 = p[3] - p[0]
        signed_vol = float(np.dot(v1, np.cross(v2, v3)))

        # @ -> negative volume (anticlockwise), @@ -> positive (clockwise)
        # in the standard SMILES convention for the first 3 neighbors
        want_positive = (atom.chirality == "@@")

        if (signed_vol > 0) != want_positive:
            # Swap positions of the last two neighbors to invert chirality
            i, j = nbrs[2], nbrs[3]
            mol.atoms[i].position, mol.atoms[j].position = (
                mol.atoms[j].position.copy(),
                mol.atoms[i].position.copy(),
            )


def _enforce_ez(mol: Molecule) -> None:
    """Enforce E/Z at double bonds with directional bond markers.

    Directional bonds (/ and \\) in SMILES specify E/Z geometry.
    After optimization, verify the dihedral sign is correct and
    flip if needed.

    Note: This is a basic implementation. Full E/Z enforcement would
    require tracking which bonds were marked / vs \\ during parsing,
    which is not currently stored on the Bond objects. For now, we
    rely on the optimizer preserving the initial E/Z orientation
    from the DG embedding.
    """
    # Future: when directional bond markers are stored, check and fix here
    pass
