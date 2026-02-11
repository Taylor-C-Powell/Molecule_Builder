"""Trajectory storage and sub-femtosecond interpolation.

Stores MD frames and provides CubicSpline interpolation for
requesting atomic positions at arbitrary time resolution, enabling
extreme slow-motion visualization of atomic interactions.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline


@dataclass
class TrajectoryFrame:
    """A single snapshot from an MD simulation.

    Attributes
    ----------
    time_fs : float
        Simulation time in femtoseconds.
    positions : ndarray of shape (n_atoms, 3)
        Atomic positions in Angstroms.
    velocities : ndarray of shape (n_atoms, 3) or None
        Atomic velocities in A/fs (optional).
    energy_kinetic : float
        Kinetic energy in kJ/mol.
    energy_potential : float
        Potential energy in kJ/mol.
    bond_orders : dict[tuple[int, int], float] or None
        Fractional bond orders (for mechanism animations).
    """
    time_fs: float
    positions: np.ndarray
    velocities: np.ndarray | None = None
    energy_kinetic: float = 0.0
    energy_potential: float = 0.0
    bond_orders: dict[tuple[int, int], float] | None = None

    @property
    def energy_total(self) -> float:
        return self.energy_kinetic + self.energy_potential


class Trajectory:
    """Ordered collection of MD frames with CubicSpline interpolation.

    The CubicSpline interpolation is the key to sub-femtosecond slow
    motion: MD runs at discrete timesteps (e.g. 0.5 fs), but the
    animation pipeline can request positions at any intermediate time.

    Parameters
    ----------
    n_atoms : int
        Number of atoms (must be consistent across all frames).
    """

    def __init__(self, n_atoms: int):
        self.n_atoms = n_atoms
        self.frames: list[TrajectoryFrame] = []
        self._spline: CubicSpline | None = None
        self._spline_dirty = True

    def add_frame(self, frame: TrajectoryFrame):
        """Append a frame to the trajectory."""
        self.frames.append(frame)
        self._spline_dirty = True

    @property
    def n_frames(self) -> int:
        return len(self.frames)

    @property
    def times(self) -> np.ndarray:
        """Array of frame times in fs."""
        return np.array([f.time_fs for f in self.frames])

    @property
    def t_start(self) -> float:
        """Start time of the trajectory in fs."""
        if not self.frames:
            return 0.0
        return self.frames[0].time_fs

    @property
    def t_end(self) -> float:
        """End time of the trajectory in fs."""
        if not self.frames:
            return 0.0
        return self.frames[-1].time_fs

    @property
    def duration(self) -> float:
        """Total duration in fs."""
        return self.t_end - self.t_start

    def _build_spline(self):
        """Build or rebuild the CubicSpline interpolator."""
        if len(self.frames) < 2:
            self._spline = None
            self._spline_dirty = False
            return

        times = self.times
        # Flatten positions to (n_frames, n_atoms*3) for spline
        pos_flat = np.array([f.positions.ravel() for f in self.frames])
        self._spline = CubicSpline(times, pos_flat)
        self._spline_dirty = False

    def at_time(self, t_fs: float) -> np.ndarray:
        """Interpolate atomic positions at an arbitrary time.

        Parameters
        ----------
        t_fs : float
            Time in femtoseconds.  Clamped to trajectory bounds.

        Returns
        -------
        ndarray of shape (n_atoms, 3)
            Interpolated positions in Angstroms.
        """
        if not self.frames:
            raise ValueError("Trajectory is empty")

        if len(self.frames) == 1:
            return self.frames[0].positions.copy()

        if self._spline_dirty:
            self._build_spline()

        t_fs = np.clip(t_fs, self.t_start, self.t_end)
        flat = self._spline(t_fs)
        return flat.reshape(self.n_atoms, 3)

    def bond_order_at_time(self, t_fs: float, i: int, j: int) -> float:
        """Interpolate fractional bond order between atoms i and j.

        Parameters
        ----------
        t_fs : float
            Time in femtoseconds.
        i, j : int
            Atom indices.

        Returns
        -------
        float
            Interpolated bond order (0.0 = no bond, 1.0 = single, etc.).
        """
        if not self.frames:
            return 0.0

        key = (min(i, j), max(i, j))
        times = []
        orders = []

        for frame in self.frames:
            if frame.bond_orders is not None and key in frame.bond_orders:
                times.append(frame.time_fs)
                orders.append(frame.bond_orders[key])

        if len(times) < 2:
            if orders:
                return orders[0]
            return 0.0

        t_fs = np.clip(t_fs, times[0], times[-1])
        spline = CubicSpline(np.array(times), np.array(orders))
        return float(np.clip(spline(t_fs), 0.0, 4.0))

    def energies(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return energy arrays.

        Returns
        -------
        times : ndarray
            Frame times in fs.
        ke : ndarray
            Kinetic energies in kJ/mol.
        pe : ndarray
            Potential energies in kJ/mol.
        """
        times = self.times
        ke = np.array([f.energy_kinetic for f in self.frames])
        pe = np.array([f.energy_potential for f in self.frames])
        return times, ke, pe

    def positions_at_times(self, time_array: np.ndarray) -> np.ndarray:
        """Interpolate positions at multiple times.

        Parameters
        ----------
        time_array : ndarray of shape (n_times,)
            Times in fs.

        Returns
        -------
        ndarray of shape (n_times, n_atoms, 3)
            Interpolated positions.
        """
        if self._spline_dirty:
            self._build_spline()

        if self._spline is None:
            if self.frames:
                pos = self.frames[0].positions
                return np.tile(pos, (len(time_array), 1, 1))
            raise ValueError("Trajectory is empty")

        t_clamped = np.clip(time_array, self.t_start, self.t_end)
        flat = self._spline(t_clamped)  # (n_times, n_atoms*3)
        return flat.reshape(len(time_array), self.n_atoms, 3)
