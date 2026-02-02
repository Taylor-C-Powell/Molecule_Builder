"""
3D geometry utilities for molecular coordinate construction.

Public API for functions that were previously private (_normalize,
_rotation_matrix, etc.) in molecular_conformations.py.  Also includes
coordinate conversions from quantum_wavefunctions.py.
"""

import math
import numpy as np

# Tolerance constants for near-zero comparisons
_ZERO_VECTOR_TOL = 1e-12     # For zero-length vector detection
_COLLINEAR_TOL = 1e-10       # For collinear atom detection


# ===================================================================
# Vector utilities
# ===================================================================

def normalize(v: np.ndarray) -> np.ndarray:
    """Return unit vector, or zero vector if input is near-zero."""
    n = np.linalg.norm(v)
    if n < _ZERO_VECTOR_TOL:
        return np.zeros(3)
    return v / n


def rotation_matrix(axis: np.ndarray, theta: float) -> np.ndarray:
    """Rodrigues rotation matrix: rotate by *theta* radians about *axis*."""
    u = normalize(axis)
    K = np.array([
        [0.0,   -u[2],  u[1]],
        [u[2],   0.0,  -u[0]],
        [-u[1],  u[0],  0.0],
    ])
    return np.eye(3) + math.sin(theta) * K + (1.0 - math.cos(theta)) * (K @ K)


# ===================================================================
# Z-matrix atom placement
# ===================================================================

def place_atom_zmatrix(pos_ref: np.ndarray,
                       pos_angle_ref: np.ndarray,
                       pos_dihedral_ref: np.ndarray,
                       bond_length: float,
                       bond_angle_deg: float,
                       dihedral_deg: float) -> np.ndarray:
    """Place an atom using internal (z-matrix) coordinates.

    Given three reference positions (j, i, k), place a new atom m so that
    distance(j,m)=bond_length, angle(i-j-m)=bond_angle_deg, and
    dihedral(k-i-j-m)=dihedral_deg.

    Parameters
    ----------
    pos_ref : position of bonded atom j.
    pos_angle_ref : position of atom i (defines bond angle).
    pos_dihedral_ref : position of atom k (defines dihedral).
    bond_length, bond_angle_deg, dihedral_deg : internal coordinates.
    """
    theta = math.radians(bond_angle_deg)
    phi = math.radians(dihedral_deg)

    v_ij = normalize(pos_ref - pos_angle_ref)
    v_ki = normalize(pos_angle_ref - pos_dihedral_ref)

    # Plane normal
    n = np.cross(v_ki, v_ij)
    if np.linalg.norm(n) < _COLLINEAR_TOL:
        # Collinear fallback
        perp = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(v_ij, perp)) > 0.9:
            perp = np.array([0.0, 1.0, 0.0])
        n = np.cross(v_ij, perp)
    n = normalize(n)

    # In-plane perpendicular
    d = normalize(np.cross(n, v_ij))

    new_dir = (
        -v_ij * math.cos(theta)
        + d * math.sin(theta) * math.cos(phi)
        + n * math.sin(theta) * math.sin(phi)
    )
    return pos_ref + bond_length * normalize(new_dir)


# ===================================================================
# Tetrahedral direction computation
# ===================================================================

def available_tetrahedral_dirs(
        existing_dirs: list[np.ndarray],
        count: int,
) -> list[np.ndarray]:
    """Compute *count* tetrahedral directions that avoid *existing_dirs*.

    Given the unit-direction vectors of bonds already on a tetrahedral
    centre, returns *count* new unit-direction vectors that complete the
    tetrahedral arrangement.
    """
    n = len(existing_dirs)
    tet_angle = math.acos(-1.0 / 3.0)  # ~109.47 deg

    if n == 0:
        dirs = [
            np.array([1, 1, 1]) / math.sqrt(3),
            np.array([1, -1, -1]) / math.sqrt(3),
            np.array([-1, 1, -1]) / math.sqrt(3),
            np.array([-1, -1, 1]) / math.sqrt(3),
        ]
        return dirs[:count]

    if n == 1:
        v0 = normalize(existing_dirs[0])
        perp = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(v0, perp)) > 0.9:
            perp = np.array([0.0, 1.0, 0.0])
        p = normalize(np.cross(v0, perp))
        q = normalize(np.cross(v0, p))

        ct = math.cos(tet_angle)
        st = math.sin(tet_angle)

        dirs = []
        for i in range(count):
            phi = math.radians(120.0 * i)
            d = v0 * ct + (p * math.cos(phi) + q * math.sin(phi)) * st
            dirs.append(normalize(d))
        return dirs

    if n == 2:
        v0 = normalize(existing_dirs[0])
        v1 = normalize(existing_dirs[1])
        bisector = normalize(v0 + v1)
        out_of_plane = normalize(np.cross(v0, v1))

        bv0 = float(np.dot(bisector, v0))
        if abs(bv0) > _COLLINEAR_TOL:
            cb = (1.0 / 3.0) / bv0
        else:
            cb = 0.0
        cb = np.clip(cb, -1.0, 1.0)
        sb = math.sqrt(max(0.0, 1.0 - cb * cb))

        dirs = []
        if count >= 1:
            dirs.append(normalize(-bisector * cb + out_of_plane * sb))
        if count >= 2:
            dirs.append(normalize(-bisector * cb - out_of_plane * sb))
        return dirs

    if n >= 3:
        centroid = normalize(sum(existing_dirs))
        return [normalize(-centroid)][:count]

    return []


# ===================================================================
# Hydrogen placement helper
# ===================================================================

def add_sp3_hydrogens(mol, carbon_idx: int, count: int):
    """Add *count* hydrogens to fill remaining tetrahedral positions.

    Uses direct 3D geometry to find available tetrahedral directions
    relative to the existing bonds on the carbon.

    Parameters
    ----------
    mol : molecule object with .atoms, .neighbors(), .add_atom(), .add_bond()
    carbon_idx : index of the atom to add hydrogens to
    count : number of hydrogens to add
    """
    from molbuilder.core.bond_data import bond_length as _bond_length

    c_pos = mol.atoms[carbon_idx].position
    existing = mol.neighbors(carbon_idx)
    if not existing:
        return

    existing_dirs = [
        normalize(mol.atoms[n].position - c_pos) for n in existing
    ]

    CH = _bond_length("C", "H", 1)
    new_dirs = available_tetrahedral_dirs(existing_dirs, count)

    for d in new_dirs:
        h_pos = c_pos + CH * d
        h_idx = mol.add_atom("H", h_pos)
        mol.add_bond(carbon_idx, h_idx, order=1, rotatable=False)


# ===================================================================
# Coordinate conversions (from quantum_wavefunctions.py)
# ===================================================================

def cartesian_to_spherical(x, y, z):
    """Convert (x, y, z) to (r, theta, phi).

    theta: polar angle [0, pi],  phi: azimuthal [0, 2*pi].
    """
    x, y, z = np.asarray(x, float), np.asarray(y, float), np.asarray(z, float)
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.where(r > 0, np.arccos(np.clip(z / np.where(r > 0, r, 1.0), -1, 1)), 0.0)
    phi = np.arctan2(y, x) % (2 * np.pi)
    return r, theta, phi


def spherical_to_cartesian(r, theta, phi):
    """Convert (r, theta, phi) to (x, y, z)."""
    r = np.asarray(r, float)
    theta = np.asarray(theta, float)
    phi = np.asarray(phi, float)
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z
