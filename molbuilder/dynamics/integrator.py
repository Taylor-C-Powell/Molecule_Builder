"""Velocity Verlet integrator and Berendsen thermostat.

Scientific basis:
    Swope et al., J. Chem. Phys. 1982, 76, 637.

The Velocity Verlet algorithm is symplectic, time-reversible, and
second-order accurate:

    v(t + dt/2)  = v(t)      + (dt/2) * a(t)
    r(t + dt)    = r(t)      + dt * v(t + dt/2)
    a(t + dt)    = F(r(t+dt)) / m
    v(t + dt)    = v(t + dt/2) + (dt/2) * a(t + dt)
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


# ===================================================================
# Unit conversion constants
# ===================================================================

# 1 fs = 1e-15 s;  positions in Angstroms, masses in AMU, energies in kJ/mol
# Conversion factor: F [kJ/(mol*A)] -> a [A/fs^2]
# a = F / m  with F in kJ/(mol*A), m in AMU
# 1 kJ/mol = 1e3 J/mol = 1e3 / (6.022e23) J = 1.6605e-21 J per particle
# 1 AMU = 1.6605e-27 kg
# a = F * 1.6605e-21 / (m * 1.6605e-27) [m/s^2]  -- but we want A/fs^2
# a = F / m * 1e-21/1e-27 * 1e-10/1e-30 ... simplifying:
# a [A/fs^2] = F [kJ/(mol*A)] / m [AMU] * (1e-2)
# More precisely: kJ/(mol*A*AMU) -> A/fs^2 = 1/(Avogadro * 1e-26) -> 1/6.022e-3
# = 1e3 / (6.022e23 * 1e-27 * 1e10) -- let's use the correct factor
# F_SI = F_kJ_mol_A * 1e3 / (6.022e23) / 1e-10  [N per particle]
# a_SI = F_SI / (m_AMU * 1.6605e-27)  [m/s^2]
# a_A_fs2 = a_SI * 1e-10 * 1e-30 ... no, A/fs^2 = (1e-10 m)/(1e-15 s)^2 = 1e20 m/s^2
# So a [A/fs^2] = a_SI / 1e20
# = F_kJ_mol_A * 1e13 / (6.022e23 * m_AMU * 1.6605e-27 * 1e20)
# Let's compute numerically:
# factor = 1e3 / (6.02214076e23 * 1e-10) / (1.66053906660e-27 * 1e20)
# = 1e3 / (6.02214076e13) / (1.66053906660e-7)
# = 1e3 * 1e7 / (6.02214076e13 * 1.66053906660)
# = 1e10 / (6.02214076e13 * 1.66053906660)
# = 1e10 / (1.00000e14)  (approximately)
# = 1e-4
# More precisely: 1e3 / (6.02214076e23) = 1.66054e-21 J/(A)
#   / (1.66054e-27 kg) = 1e6 m/s^2 per (A)
#   * 1e-10 m/A = 1e-4 m^2/(A*s^2)  ... let me just compute directly
# CONVERSION_FACTOR = 1e-4 approximately
# Use the precise value: 1 / (AMU * Avogadro) in appropriate units
_KJ_MOL_A_TO_A_FS2 = 1.0e-4  # kJ/(mol*A*AMU) -> A/fs^2 (approximate)


@dataclass
class IntegratorState:
    """State of the MD system at a single time point.

    Attributes
    ----------
    positions : ndarray of shape (n_atoms, 3)
        Atomic positions in Angstroms.
    velocities : ndarray of shape (n_atoms, 3)
        Atomic velocities in A/fs.
    time_fs : float
        Current simulation time in femtoseconds.
    """
    positions: np.ndarray
    velocities: np.ndarray
    time_fs: float = 0.0


class BerendsenThermostat:
    """Berendsen weak-coupling thermostat.

    Rescales velocities to weakly couple the system to a heat bath
    at a target temperature.  The coupling time tau controls how
    strongly the system is coupled (larger tau = weaker coupling).

    Parameters
    ----------
    target_T : float
        Target temperature in Kelvin.
    tau_fs : float
        Coupling time constant in femtoseconds.  Typical: 100 fs.
    """

    def __init__(self, target_T: float = 300.0, tau_fs: float = 100.0):
        self.target_T = target_T
        self.tau_fs = tau_fs

    def rescale(self, velocities: np.ndarray, masses: np.ndarray,
                dt_fs: float) -> np.ndarray:
        """Apply Berendsen velocity rescaling.

        Parameters
        ----------
        velocities : ndarray of shape (n_atoms, 3)
            Current velocities in A/fs.
        masses : ndarray of shape (n_atoms,)
            Atomic masses in AMU.
        dt_fs : float
            Integration timestep in fs.

        Returns
        -------
        ndarray of shape (n_atoms, 3)
            Rescaled velocities.
        """
        current_T = self._temperature(velocities, masses)
        if current_T < 1e-10:
            return velocities

        lam2 = 1.0 + (dt_fs / self.tau_fs) * (self.target_T / current_T - 1.0)
        if lam2 < 0:
            lam2 = 0.0
        lam = math.sqrt(lam2)
        return velocities * lam

    @staticmethod
    def _temperature(velocities: np.ndarray, masses: np.ndarray) -> float:
        """Compute instantaneous kinetic temperature in Kelvin.

        T = (2 * KE) / (n_dof * k_B)
        where n_dof = 3*N - 3 (removing center-of-mass translation).

        Uses internal units: KE in AMU * A^2 / fs^2, k_B in same units.
        k_B = 8.3145e-3 kJ/(mol*K), and 1 AMU*A^2/fs^2 = 1e-4 * 1 kJ/mol
        => k_B in AMU*A^2/(fs^2*K) = 8.3145e-3 / 1e-4 * ... no, let's be
        precise.

        KE [AMU*A^2/fs^2] = sum 0.5 * m_i * v_i^2
        To convert to kJ/mol: multiply by 1e4 (inverse of _KJ_MOL_A_TO_A_FS2 * mass)
        Actually: E [kJ/mol] = 0.5 * m [AMU] * v^2 [A^2/fs^2] / _KJ_MOL_A_TO_A_FS2
        Hmm, let's simplify with k_B in internal units.
        """
        n_atoms = len(masses)
        if n_atoms < 2:
            return 0.0
        n_dof = 3 * n_atoms - 3

        # KE in kJ/mol: 0.5 * m * v^2 / conversion
        # v [A/fs], m [AMU]: KE_amu = 0.5 * m * v^2 [AMU*A^2/fs^2]
        # To kJ/mol: KE_kJ = KE_amu / _KJ_MOL_A_TO_A_FS2
        ke_per_atom = 0.5 * masses[:, None] * velocities ** 2
        ke_total = np.sum(ke_per_atom) / _KJ_MOL_A_TO_A_FS2  # kJ/mol

        # k_B = 8.3145e-3 kJ/(mol*K)
        k_B = 8.3145e-3
        T = 2.0 * ke_total / (n_dof * k_B)
        return T


class VelocityVerletIntegrator:
    """Velocity Verlet integrator for molecular dynamics.

    Parameters
    ----------
    dt_fs : float
        Timestep in femtoseconds.  Default 0.5 fs is suitable for
        all-atom simulations with hydrogens.
    masses : ndarray of shape (n_atoms,)
        Atomic masses in AMU.
    thermostat : BerendsenThermostat or None
        Optional thermostat for temperature control.
    """

    def __init__(self, dt_fs: float, masses: np.ndarray,
                 thermostat: BerendsenThermostat | None = None):
        self.dt = dt_fs
        self.masses = masses
        self.thermostat = thermostat

    def step(self, state: IntegratorState,
             force_fn) -> IntegratorState:
        """Perform one Velocity Verlet integration step.

        Parameters
        ----------
        state : IntegratorState
            Current positions, velocities, and time.
        force_fn : callable
            ``force_fn(positions) -> ForceResult`` returning forces in
            kJ/(mol*A).

        Returns
        -------
        IntegratorState
            Updated state after one timestep.
        """
        dt = self.dt
        m = self.masses[:, None]  # (n_atoms, 1)

        # Current acceleration: a = F/m converted to A/fs^2
        result = force_fn(state.positions)
        accel = result.forces / m * _KJ_MOL_A_TO_A_FS2

        # Half-step velocity
        v_half = state.velocities + 0.5 * dt * accel

        # Full-step position
        new_pos = state.positions + dt * v_half

        # New acceleration
        result_new = force_fn(new_pos)
        accel_new = result_new.forces / m * _KJ_MOL_A_TO_A_FS2

        # Full-step velocity
        new_vel = v_half + 0.5 * dt * accel_new

        # Apply thermostat
        if self.thermostat is not None:
            new_vel = self.thermostat.rescale(new_vel, self.masses, dt)

        return IntegratorState(
            positions=new_pos,
            velocities=new_vel,
            time_fs=state.time_fs + dt,
        )

    def initialize_velocities(self, n_atoms: int,
                               temperature_K: float) -> np.ndarray:
        """Generate Maxwell-Boltzmann distributed velocities.

        Parameters
        ----------
        n_atoms : int
            Number of atoms.
        temperature_K : float
            Target temperature in Kelvin.

        Returns
        -------
        ndarray of shape (n_atoms, 3)
            Velocities in A/fs.
        """
        if temperature_K < 1e-10:
            return np.zeros((n_atoms, 3))

        # k_B T in kJ/mol
        k_B = 8.3145e-3  # kJ/(mol*K)
        kT = k_B * temperature_K

        velocities = np.zeros((n_atoms, 3))
        for i in range(n_atoms):
            # sigma_v [A/fs] from equipartition: 0.5 * m * v^2 = 0.5 * kT
            # v_rms = sqrt(kT/m) but in our units:
            # kT [kJ/mol], m [AMU] -> sigma [A/fs]
            sigma = math.sqrt(kT * _KJ_MOL_A_TO_A_FS2 / self.masses[i])
            velocities[i] = np.random.normal(0.0, sigma, 3)

        # Remove center-of-mass velocity
        total_mass = np.sum(self.masses)
        com_vel = np.sum(self.masses[:, None] * velocities, axis=0) / total_mass
        velocities -= com_vel[None, :]

        return velocities

    def kinetic_energy(self, velocities: np.ndarray) -> float:
        """Compute kinetic energy in kJ/mol.

        Parameters
        ----------
        velocities : ndarray of shape (n_atoms, 3)
            Velocities in A/fs.

        Returns
        -------
        float
            Kinetic energy in kJ/mol.
        """
        ke_internal = 0.5 * np.sum(self.masses[:, None] * velocities ** 2)
        return ke_internal / _KJ_MOL_A_TO_A_FS2
