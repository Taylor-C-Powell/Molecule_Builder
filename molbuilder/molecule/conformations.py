"""
Conformational analysis utilities.

Provides classification of dihedral angles into named conformations
and torsional energy scanning.
"""

from __future__ import annotations


from molbuilder.molecule.graph import Molecule, ConformationType


def classify_conformation(dihedral_deg: float) -> ConformationType:
    """Classify a dihedral angle into a named conformation.

    Ranges (normalised to -180..180):
        |d| < 10          -> ECLIPSED  (syn-periplanar)
        |d| ~ 60          -> GAUCHE
        |d| ~ 120         -> ECLIPSED  (anti-clinal)
        |d| > 170         -> ANTI      (antiperiplanar)
        otherwise         -> CUSTOM
    """
    d = dihedral_deg % 360
    if d > 180:
        d -= 360

    if abs(d) < 10:
        return ConformationType.ECLIPSED
    if abs(abs(d) - 60) < 15:
        return ConformationType.GAUCHE
    if abs(abs(d) - 180) < 15:
        return ConformationType.ANTI
    if abs(abs(d) - 120) < 15:
        return ConformationType.ECLIPSED
    return ConformationType.CUSTOM


def scan_torsion(mol: Molecule, j: int, k: int,
                 ref_i: int, ref_l: int,
                 steps: int = 36) -> list[tuple[float, float]]:
    """Scan torsional energy as a function of dihedral angle.

    Rotates the k-side of bond j-k through 360 degrees, computing
    strain at each step.  Restores original geometry between each
    step to avoid cumulative numerical drift.

    Returns list of (angle_deg, energy_kJ_per_mol).
    """
    originals = [a.position.copy() for a in mol.atoms]

    results: list[tuple[float, float]] = []
    step_size = 360.0 / steps

    for s in range(steps):
        # Reset to original before each step
        for i, pos in enumerate(originals):
            mol.atoms[i].position = pos.copy()

        target = -180.0 + s * step_size
        mol.set_dihedral(ref_i, j, k, ref_l, target)
        energy = mol.torsional_energy(j, k)
        results.append((target, energy.total_kj_per_mol))

    # Restore original
    for i, pos in enumerate(originals):
        mol.atoms[i].position = pos

    return results
