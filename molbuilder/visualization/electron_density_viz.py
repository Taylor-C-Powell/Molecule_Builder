"""Electron cloud rendering for bond events.

Uses an LCAO (Linear Combination of Atomic Orbitals) approximation to
visualize electron density during bond formation and breaking.  The
bonding orbital is modelled as:

    psi_bond = c_a * phi_a(r - R_a) + c_b * phi_b(r - R_b)

where phi are hydrogen-like orbitals evaluated with Slater effective
nuclear charges, and c_a, c_b are mixing coefficients determined by
the bond order.

Reuses radial_wavefunction, real_spherical_harmonic, and
wavefunction_real from molbuilder.atomic.wavefunctions.
"""

from __future__ import annotations

import math

import numpy as np

from molbuilder.core.constants import BOHR_RADIUS_M
from molbuilder.visualization.theme import ELECTRON_CLOUD_COLOR


# Bohr radius in Angstroms
_A0_ANGSTROM = BOHR_RADIUS_M * 1e10  # ~0.529 A


def _hex_to_rgba(hex_color: str, alpha: float = 0.3) -> tuple:
    """Convert hex color string to RGBA tuple (0-1 scale)."""
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4))
    return (r, g, b, alpha)


def _slater_zeff_simple(Z: int) -> float:
    """Simplified Slater Z_eff for outermost valence s/p orbital.

    Uses a rough approximation: Z_eff ~ Z - S where S is the total
    screening.  For quick rendering this is sufficient.
    """
    # Very simplified: use known values for common elements
    _ZEFF = {
        1: 1.0, 6: 3.25, 7: 3.9, 8: 4.55, 9: 5.2,
        15: 4.8, 16: 5.45, 17: 6.1, 35: 9.0, 53: 12.0,
    }
    return _ZEFF.get(Z, max(1.0, Z * 0.3))


def _valence_nl(Z: int) -> tuple[int, int]:
    """Return (n, l) for the outermost valence orbital."""
    if Z <= 2:
        return (1, 0)
    if Z <= 10:
        if Z <= 4:
            return (2, 0)
        return (2, 1)
    if Z <= 18:
        if Z <= 12:
            return (3, 0)
        return (3, 1)
    if Z <= 36:
        return (4, 0)
    return (5, 0)


class ElectronDensityRenderer:
    """Renders electron density clouds using Monte Carlo point sampling.

    Generates scatter points weighted by the LCAO bonding orbital
    probability density, producing a visual representation of the
    electron cloud during bond events.

    Parameters
    ----------
    n_points : int
        Number of scatter points to generate per bond.  More points
        give smoother clouds but slower rendering.
    """

    def __init__(self, n_points: int = 5000):
        self.n_points = n_points
        self._rng = np.random.default_rng(42)

    def compute_bond_density(self, pos_a: np.ndarray, pos_b: np.ndarray,
                              z_a: int, z_b: int,
                              bond_order: float = 1.0,
                              ) -> tuple[np.ndarray, np.ndarray]:
        """Compute electron cloud scatter points for a bond.

        Uses rejection sampling: generates random points in the bonding
        region and accepts them with probability proportional to
        |psi_bond|^2.

        Parameters
        ----------
        pos_a, pos_b : ndarray of shape (3,)
            Atomic positions in Angstroms.
        z_a, z_b : int
            Atomic numbers.
        bond_order : float
            Fractional bond order (0 to 3).  Controls cloud density.

        Returns
        -------
        points : ndarray of shape (n_accepted, 3)
            Accepted scatter point positions.
        colors : ndarray of shape (n_accepted, 4)
            RGBA colors for each point.
        """
        if bond_order < 0.05:
            return np.empty((0, 3)), np.empty((0, 4))

        midpoint = 0.5 * (pos_a + pos_b)
        bond_vec = pos_b - pos_a
        bond_len = np.linalg.norm(bond_vec)
        if bond_len < 1e-6:
            return np.empty((0, 3)), np.empty((0, 4))

        # Effective nuclear charges
        zeff_a = _slater_zeff_simple(z_a)
        zeff_b = _slater_zeff_simple(z_b)

        # Valence orbital quantum numbers
        n_a, l_a = _valence_nl(z_a)
        n_b, l_b = _valence_nl(z_b)

        # Sampling region: ellipsoid around the bond
        spread = max(bond_len * 0.8, 1.5)

        # Generate candidate points
        pts = self._rng.normal(0, spread * 0.4, (self.n_points, 3))
        pts += midpoint

        # Evaluate approximate orbital values at each point
        # Use simplified 1s-like orbitals for speed
        r_a = np.linalg.norm(pts - pos_a, axis=1)
        r_b = np.linalg.norm(pts - pos_b, axis=1)

        # Slater-type orbital: phi ~ r^(n-1) * exp(-zeff * r / (n * a0))
        decay_a = zeff_a / (n_a * _A0_ANGSTROM)
        decay_b = zeff_b / (n_b * _A0_ANGSTROM)

        phi_a = np.exp(-decay_a * r_a)
        phi_b = np.exp(-decay_b * r_b)

        # LCAO bonding combination
        c_a = 1.0 / math.sqrt(2.0)
        c_b = 1.0 / math.sqrt(2.0)
        psi_bond = c_a * phi_a + c_b * phi_b
        density = psi_bond ** 2

        # Normalize and apply bond order scaling
        max_dens = np.max(density)
        if max_dens < 1e-30:
            return np.empty((0, 3)), np.empty((0, 4))
        prob = density / max_dens * min(bond_order, 1.5)

        # Rejection sampling
        rng_vals = self._rng.uniform(0, 1, self.n_points)
        accept = rng_vals < prob

        accepted_pts = pts[accept]
        n_accepted = len(accepted_pts)

        if n_accepted == 0:
            return np.empty((0, 3)), np.empty((0, 4))

        # Colors: base color with alpha proportional to density
        base_rgba = _hex_to_rgba(ELECTRON_CLOUD_COLOR, 0.4)
        colors = np.tile(base_rgba, (n_accepted, 1))
        # Modulate alpha by density
        accepted_density = density[accept]
        alpha_scale = accepted_density / max_dens
        colors[:, 3] = 0.1 + 0.4 * alpha_scale * min(bond_order, 1.5)

        return accepted_pts, colors

    def render_on_axis(self, ax, pos_a: np.ndarray, pos_b: np.ndarray,
                       z_a: int, z_b: int,
                       bond_order: float = 1.0,
                       point_size: float = 2.0):
        """Draw electron density cloud on a matplotlib 3D axis.

        Parameters
        ----------
        ax : Axes3D
            Matplotlib 3D axis to draw on.
        pos_a, pos_b : ndarray of shape (3,)
            Atomic positions in Angstroms.
        z_a, z_b : int
            Atomic numbers.
        bond_order : float
            Fractional bond order.
        point_size : float
            Size of scatter points.
        """
        points, colors = self.compute_bond_density(
            pos_a, pos_b, z_a, z_b, bond_order)
        if len(points) == 0:
            return

        ax.scatter(
            points[:, 0], points[:, 1], points[:, 2],
            c=colors,
            s=point_size,
            alpha=0.3,
            depthshade=True,
            edgecolors="none",
        )

    def render_on_axis_2d(self, ax, pos_a: np.ndarray, pos_b: np.ndarray,
                           z_a: int, z_b: int,
                           bond_order: float = 1.0,
                           point_size: float = 2.0):
        """Draw electron density cloud on a 2D matplotlib axis.

        Uses the first two coordinates (x, y) for 2D projection.

        Parameters
        ----------
        ax : Axes
            Matplotlib 2D axis.
        pos_a, pos_b : ndarray
            Atomic positions.
        z_a, z_b : int
            Atomic numbers.
        bond_order : float
            Fractional bond order.
        point_size : float
            Scatter point size.
        """
        points, colors = self.compute_bond_density(
            pos_a, pos_b, z_a, z_b, bond_order)
        if len(points) == 0:
            return

        ax.scatter(
            points[:, 0], points[:, 1],
            c=colors,
            s=point_size,
            edgecolors="none",
        )
