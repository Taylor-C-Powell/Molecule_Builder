"""Builtin pure-Python 3D coordinate generation pipeline.

Orchestrates: DG embed -> ring templates -> FF optimize -> stereo enforce.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from molbuilder.molecule.graph import Molecule


def generate_builtin(mol: Molecule, seed: int = 42) -> None:
    """Generate 3D coordinates using the builtin DG + FF pipeline.

    Modifies atom positions in-place.

    Pipeline
    --------
    1. Build distance bounds matrix from topology
    2. DG embed via eigendecomposition
    3. Apply ring templates (Kabsch-align ideal ring geometries)
    4. FF optimize with ideal equilibrium values (L-BFGS-B)
    5. Enforce chirality (swap atoms if signed volume is wrong)
    """
    from molbuilder.coords._bounds import build_bounds
    from molbuilder.coords._embed import embed
    from molbuilder.coords._ring_templates import apply_ring_templates
    from molbuilder.coords._optimize import optimize_geometry
    from molbuilder.coords._stereo import enforce_stereo

    n = len(mol.atoms)
    if n == 0:
        return

    if n == 1:
        import numpy as np
        mol.atoms[0].position = np.array([0.0, 0.0, 0.0])
        return

    # Step 1-2: Distance geometry embedding
    lower, upper = build_bounds(mol)
    coords = embed(lower, upper, n, seed=seed)

    # Set initial positions
    for i, atom in enumerate(mol.atoms):
        atom.position = coords[i]

    # Step 3: Ring template correction
    apply_ring_templates(mol)

    # Step 4: Force field optimization with ideal angles
    optimize_geometry(mol, max_iter=200)

    # Step 5: Chirality enforcement
    enforce_stereo(mol)
