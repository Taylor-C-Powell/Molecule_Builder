"""Tests for the molecular dynamics engine.

Covers ForceField parameterization, force computation, Velocity Verlet
integration, trajectory interpolation, and MDSimulation orchestration.
"""

import math
import numpy as np
import pytest

from molbuilder.molecule.graph import Molecule, Hybridization
from molbuilder.molecule.builders import build_ethane
from molbuilder.core.bond_data import bond_length, SP3_ANGLE


# ===================================================================
# ForceField tests
# ===================================================================

class TestForceFieldParams:
    """Test ForceField.from_molecule parameterization."""

    def test_from_ethane(self):
        """ForceField should parameterize from ethane without error."""
        from molbuilder.dynamics.forcefield import ForceField
        mol = build_ethane(60.0)
        ff = ForceField.from_molecule(mol)
        assert ff.params.n_atoms == len(mol.atoms)
        assert len(ff.params.masses) == len(mol.atoms)
        assert ff.params.bond_indices.shape[0] == len(mol.bonds)

    def test_masses_positive(self):
        """All masses should be positive."""
        from molbuilder.dynamics.forcefield import ForceField
        mol = build_ethane(60.0)
        ff = ForceField.from_molecule(mol)
        assert np.all(ff.params.masses > 0)

    def test_lj_params_positive(self):
        """LJ sigma and epsilon should be non-negative."""
        from molbuilder.dynamics.forcefield import ForceField
        mol = build_ethane(60.0)
        ff = ForceField.from_molecule(mol)
        assert np.all(ff.params.sigma > 0)
        assert np.all(ff.params.epsilon >= 0)

    def test_bond_r0_reasonable(self):
        """Equilibrium bond lengths should be reasonable (0.5-3 A)."""
        from molbuilder.dynamics.forcefield import ForceField
        mol = build_ethane(60.0)
        ff = ForceField.from_molecule(mol)
        assert np.all(ff.params.bond_r0 > 0.5)
        assert np.all(ff.params.bond_r0 < 3.0)

    def test_diatomic_h2(self):
        """H2 molecule should have exactly 1 bond."""
        from molbuilder.dynamics.forcefield import ForceField
        mol = Molecule("H2")
        mol.add_atom("H", [0.0, 0.0, 0.0])
        mol.add_atom("H", [0.74, 0.0, 0.0])
        mol.add_bond(0, 1, order=1)
        ff = ForceField.from_molecule(mol)
        assert ff.params.n_atoms == 2
        assert ff.params.bond_indices.shape[0] == 1
        assert abs(ff.params.bond_r0[0] - 0.74) < 0.01

    def test_charges_sum_near_zero(self):
        """Partial charges should sum to approximately zero."""
        from molbuilder.dynamics.forcefield import ForceField
        mol = build_ethane(60.0)
        ff = ForceField.from_molecule(mol)
        assert abs(np.sum(ff.params.charges)) < 0.1


class TestForceComputation:
    """Test force computation correctness."""

    def test_forces_shape(self):
        """Forces should have shape (n_atoms, 3)."""
        from molbuilder.dynamics.forcefield import ForceField
        mol = build_ethane(60.0)
        ff = ForceField.from_molecule(mol)
        pos = np.array([a.position for a in mol.atoms])
        result = ff.compute(pos)
        assert result.forces.shape == (len(mol.atoms), 3)

    def test_energy_finite(self):
        """All energy terms should be finite."""
        from molbuilder.dynamics.forcefield import ForceField
        mol = build_ethane(60.0)
        ff = ForceField.from_molecule(mol)
        pos = np.array([a.position for a in mol.atoms])
        result = ff.compute(pos)
        assert np.isfinite(result.energy_bond)
        assert np.isfinite(result.energy_angle)
        assert np.isfinite(result.energy_torsion)
        assert np.isfinite(result.energy_lj)
        assert np.isfinite(result.energy_coulomb)
        assert np.isfinite(result.energy_total)

    def test_equilibrium_forces_small(self):
        """At equilibrium geometry, net forces should be small."""
        from molbuilder.dynamics.forcefield import ForceField
        mol = build_ethane(60.0)
        ff = ForceField.from_molecule(mol)
        pos = np.array([a.position for a in mol.atoms])
        result = ff.compute(pos)
        # Forces shouldn't be huge at equilibrium
        max_force = np.max(np.abs(result.forces))
        assert max_force < 5000.0  # kJ/(mol*A)

    def test_stretched_bond_restoring_force(self):
        """Stretching a bond should produce a restoring force."""
        from molbuilder.dynamics.forcefield import ForceField
        mol = Molecule("H2 stretched")
        mol.add_atom("H", [0.0, 0.0, 0.0])
        mol.add_atom("H", [1.5, 0.0, 0.0])  # stretched beyond 0.74 A
        mol.add_bond(0, 1, order=1)
        ff = ForceField.from_molecule(mol)
        pos = np.array([a.position for a in mol.atoms])
        result = ff.compute(pos)
        # Atom 0 should be pushed toward atom 1 (positive x force)
        assert result.forces[0, 0] > 0
        # Atom 1 should be pushed toward atom 0 (negative x force)
        assert result.forces[1, 0] < 0

    def test_compressed_bond_repulsive_force(self):
        """Compressing a bond should produce a repulsive force."""
        from molbuilder.dynamics.forcefield import ForceField
        mol = Molecule("H2 compressed")
        mol.add_atom("H", [0.0, 0.0, 0.0])
        mol.add_atom("H", [0.5, 0.0, 0.0])  # compressed below 0.74 A
        mol.add_bond(0, 1, order=1)
        ff = ForceField.from_molecule(mol)
        pos = np.array([a.position for a in mol.atoms])
        result = ff.compute(pos)
        # Atom 0 pushed away from atom 1 (negative x force)
        assert result.forces[0, 0] < 0
        # Atom 1 pushed away from atom 0 (positive x force)
        assert result.forces[1, 0] > 0

    def test_bond_energy_increases_with_stretch(self):
        """Bond energy should increase when bond is stretched."""
        from molbuilder.dynamics.forcefield import ForceField
        mol_eq = Molecule("H2 eq")
        mol_eq.add_atom("H", [0.0, 0.0, 0.0])
        mol_eq.add_atom("H", [0.74, 0.0, 0.0])
        mol_eq.add_bond(0, 1, order=1)
        ff = ForceField.from_molecule(mol_eq)

        pos_eq = np.array([[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]])
        pos_str = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
        e_eq = ff.compute(pos_eq).energy_bond
        e_str = ff.compute(pos_str).energy_bond
        assert e_str > e_eq


# ===================================================================
# Integrator tests
# ===================================================================

class TestIntegrator:
    """Test the Velocity Verlet integrator."""

    def test_create_integrator(self):
        """Integrator should initialize without error."""
        from molbuilder.dynamics.integrator import VelocityVerletIntegrator
        masses = np.array([1.008, 1.008])
        integ = VelocityVerletIntegrator(dt_fs=0.5, masses=masses)
        assert integ.dt == 0.5

    def test_initialize_velocities_shape(self):
        """Initialized velocities should have correct shape."""
        from molbuilder.dynamics.integrator import VelocityVerletIntegrator
        masses = np.array([12.0, 12.0, 1.008, 1.008])
        integ = VelocityVerletIntegrator(dt_fs=0.5, masses=masses)
        vel = integ.initialize_velocities(4, 300.0)
        assert vel.shape == (4, 3)

    def test_zero_com_velocity(self):
        """Center-of-mass velocity should be zero after initialization."""
        from molbuilder.dynamics.integrator import VelocityVerletIntegrator
        masses = np.array([12.0, 12.0, 1.008, 1.008, 1.008, 1.008])
        integ = VelocityVerletIntegrator(dt_fs=0.5, masses=masses)
        vel = integ.initialize_velocities(6, 300.0)
        com_vel = np.sum(masses[:, None] * vel, axis=0) / np.sum(masses)
        assert np.allclose(com_vel, 0.0, atol=1e-12)

    def test_zero_temperature_gives_zero_velocity(self):
        """T=0 should produce zero velocities."""
        from molbuilder.dynamics.integrator import VelocityVerletIntegrator
        masses = np.array([12.0, 1.008])
        integ = VelocityVerletIntegrator(dt_fs=0.5, masses=masses)
        vel = integ.initialize_velocities(2, 0.0)
        assert np.allclose(vel, 0.0)

    def test_single_step_positions_change(self):
        """After one step, positions should change (unless at exact equilibrium)."""
        from molbuilder.dynamics.integrator import (
            VelocityVerletIntegrator, IntegratorState,
        )
        from molbuilder.dynamics.forcefield import ForceField, ForceResult

        masses = np.array([1.008, 1.008])
        integ = VelocityVerletIntegrator(dt_fs=0.5, masses=masses)

        pos = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
        vel = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        state = IntegratorState(positions=pos, velocities=vel)

        # Simple force function: harmonic spring
        def force_fn(p):
            f = np.zeros_like(p)
            r = p[1] - p[0]
            d = np.linalg.norm(r)
            k = 500.0
            r0 = 0.74
            f_mag = -k * (d - r0)
            f_vec = f_mag * r / d
            f[0] -= f_vec
            f[1] += f_vec
            return ForceResult(forces=f, energy_bond=0.5 * k * (d - r0) ** 2)

        new_state = integ.step(state, force_fn)
        assert not np.allclose(new_state.positions, pos)
        assert new_state.time_fs == 0.5

    def test_kinetic_energy_nonnegative(self):
        """Kinetic energy should always be non-negative."""
        from molbuilder.dynamics.integrator import VelocityVerletIntegrator
        masses = np.array([12.0, 1.008, 1.008])
        integ = VelocityVerletIntegrator(dt_fs=0.5, masses=masses)
        vel = integ.initialize_velocities(3, 300.0)
        ke = integ.kinetic_energy(vel)
        assert ke >= 0.0


class TestThermostat:
    """Test the Berendsen thermostat."""

    def test_thermostat_creation(self):
        """Thermostat should initialize without error."""
        from molbuilder.dynamics.integrator import BerendsenThermostat
        thermo = BerendsenThermostat(300.0, 100.0)
        assert thermo.target_T == 300.0
        assert thermo.tau_fs == 100.0

    def test_rescale_preserves_shape(self):
        """Rescaling should preserve velocity shape."""
        from molbuilder.dynamics.integrator import BerendsenThermostat
        thermo = BerendsenThermostat(300.0, 100.0)
        vel = np.random.randn(5, 3) * 0.01
        masses = np.array([12.0, 1.0, 1.0, 1.0, 1.0])
        new_vel = thermo.rescale(vel, masses, 0.5)
        assert new_vel.shape == vel.shape


# ===================================================================
# Energy conservation test
# ===================================================================

class TestEnergyConservation:
    """Test energy conservation in NVE (no thermostat) dynamics."""

    def test_harmonic_diatomic_energy_conservation(self):
        """A harmonic diatomic should conserve total energy to ~5%
        over 1000 steps with 0.5 fs timestep (no thermostat)."""
        from molbuilder.dynamics.forcefield import ForceField
        from molbuilder.dynamics.integrator import (
            VelocityVerletIntegrator, IntegratorState,
        )

        mol = Molecule("H2")
        mol.add_atom("H", [0.0, 0.0, 0.0])
        mol.add_atom("H", [0.74, 0.0, 0.0])
        mol.add_bond(0, 1, order=1)

        ff = ForceField.from_molecule(mol)
        masses = ff.params.masses
        integ = VelocityVerletIntegrator(dt_fs=0.5, masses=masses)

        pos = np.array([a.position for a in mol.atoms], dtype=float)
        # Give a small initial velocity for vibration
        vel = np.array([[0.001, 0.0, 0.0], [-0.001, 0.0, 0.0]])
        state = IntegratorState(positions=pos, velocities=vel)

        energies = []
        for _ in range(1000):
            result = ff.compute(state.positions)
            ke = integ.kinetic_energy(state.velocities)
            energies.append(ke + result.energy_total)
            state = integ.step(state, ff.compute)

        energies = np.array(energies)
        # Energy should be conserved to within 5%
        mean_e = np.mean(energies)
        if abs(mean_e) > 1e-6:
            max_deviation = np.max(np.abs(energies - mean_e)) / abs(mean_e)
            assert max_deviation < 0.05, (
                f"Energy not conserved: max deviation {max_deviation:.2%}")


# ===================================================================
# Trajectory tests
# ===================================================================

class TestTrajectory:
    """Test trajectory storage and interpolation."""

    def test_add_frames(self):
        """Should be able to add frames to a trajectory."""
        from molbuilder.dynamics.trajectory import Trajectory, TrajectoryFrame
        traj = Trajectory(n_atoms=2)
        for t in [0.0, 0.5, 1.0]:
            traj.add_frame(TrajectoryFrame(
                time_fs=t,
                positions=np.array([[0.0, 0.0, 0.0],
                                     [t * 0.1, 0.0, 0.0]]),
            ))
        assert traj.n_frames == 3
        assert traj.t_start == 0.0
        assert traj.t_end == 1.0

    def test_interpolation_at_recorded_time(self):
        """Interpolation at a recorded time should return the exact frame."""
        from molbuilder.dynamics.trajectory import Trajectory, TrajectoryFrame
        traj = Trajectory(n_atoms=2)
        for t in [0.0, 1.0, 2.0]:
            traj.add_frame(TrajectoryFrame(
                time_fs=t,
                positions=np.array([[0.0, 0.0, 0.0],
                                     [t, 0.0, 0.0]]),
            ))
        pos_at_1 = traj.at_time(1.0)
        assert abs(pos_at_1[1, 0] - 1.0) < 1e-6

    def test_interpolation_between_frames(self):
        """Interpolation between frames should give intermediate values."""
        from molbuilder.dynamics.trajectory import Trajectory, TrajectoryFrame
        traj = Trajectory(n_atoms=1)
        traj.add_frame(TrajectoryFrame(
            time_fs=0.0,
            positions=np.array([[0.0, 0.0, 0.0]])))
        traj.add_frame(TrajectoryFrame(
            time_fs=1.0,
            positions=np.array([[2.0, 0.0, 0.0]])))

        pos_mid = traj.at_time(0.5)
        # CubicSpline through 2 points = linear => midpoint = 1.0
        assert abs(pos_mid[0, 0] - 1.0) < 0.1

    def test_clamping_beyond_bounds(self):
        """Times outside trajectory bounds should be clamped."""
        from molbuilder.dynamics.trajectory import Trajectory, TrajectoryFrame
        traj = Trajectory(n_atoms=1)
        traj.add_frame(TrajectoryFrame(
            time_fs=0.0,
            positions=np.array([[0.0, 0.0, 0.0]])))
        traj.add_frame(TrajectoryFrame(
            time_fs=1.0,
            positions=np.array([[1.0, 0.0, 0.0]])))

        pos_before = traj.at_time(-1.0)
        pos_after = traj.at_time(5.0)
        assert abs(pos_before[0, 0] - 0.0) < 0.1
        assert abs(pos_after[0, 0] - 1.0) < 0.1

    def test_bond_order_interpolation(self):
        """Bond order interpolation should work with recorded values."""
        from molbuilder.dynamics.trajectory import Trajectory, TrajectoryFrame
        traj = Trajectory(n_atoms=2)
        traj.add_frame(TrajectoryFrame(
            time_fs=0.0,
            positions=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
            bond_orders={(0, 1): 0.0}))
        traj.add_frame(TrajectoryFrame(
            time_fs=1.0,
            positions=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
            bond_orders={(0, 1): 1.0}))

        bo_mid = traj.bond_order_at_time(0.5, 0, 1)
        assert 0.3 < bo_mid < 0.7

    def test_energies_method(self):
        """energies() should return correct arrays."""
        from molbuilder.dynamics.trajectory import Trajectory, TrajectoryFrame
        traj = Trajectory(n_atoms=1)
        traj.add_frame(TrajectoryFrame(
            time_fs=0.0,
            positions=np.array([[0.0, 0.0, 0.0]]),
            energy_kinetic=10.0,
            energy_potential=20.0))
        traj.add_frame(TrajectoryFrame(
            time_fs=1.0,
            positions=np.array([[1.0, 0.0, 0.0]]),
            energy_kinetic=15.0,
            energy_potential=18.0))

        times, ke, pe = traj.energies()
        assert len(times) == 2
        assert ke[0] == 10.0
        assert pe[1] == 18.0

    def test_positions_at_times(self):
        """positions_at_times should return correct shape."""
        from molbuilder.dynamics.trajectory import Trajectory, TrajectoryFrame
        traj = Trajectory(n_atoms=2)
        for t in [0.0, 1.0, 2.0]:
            traj.add_frame(TrajectoryFrame(
                time_fs=t,
                positions=np.random.randn(2, 3)))

        t_array = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
        positions = traj.positions_at_times(t_array)
        assert positions.shape == (5, 2, 3)


# ===================================================================
# Simulation tests
# ===================================================================

class TestMDSimulation:
    """Test the MDSimulation orchestrator."""

    def test_run_short_simulation(self):
        """Should run a short simulation and return a trajectory."""
        from molbuilder.dynamics.simulation import MDSimulation
        mol = Molecule("H2")
        mol.add_atom("H", [0.0, 0.0, 0.0])
        mol.add_atom("H", [0.74, 0.0, 0.0])
        mol.add_bond(0, 1, order=1)

        sim = MDSimulation(mol, dt_fs=0.5, temperature_K=100.0,
                           thermostat=True, record_interval=5)
        traj = sim.run(50)
        assert traj.n_frames >= 10  # 50/5 = 10 frames + initial
        assert traj.duration > 0

    def test_trajectory_has_positions(self):
        """Trajectory frames should have valid positions."""
        from molbuilder.dynamics.simulation import MDSimulation
        mol = Molecule("H2")
        mol.add_atom("H", [0.0, 0.0, 0.0])
        mol.add_atom("H", [0.74, 0.0, 0.0])
        mol.add_bond(0, 1, order=1)

        sim = MDSimulation(mol, dt_fs=0.5, temperature_K=100.0)
        traj = sim.run(20)
        for frame in traj.frames:
            assert frame.positions.shape == (2, 3)
            assert np.all(np.isfinite(frame.positions))

    def test_ethane_simulation(self):
        """Ethane simulation should run without error."""
        from molbuilder.dynamics.simulation import MDSimulation
        mol = build_ethane(60.0)
        sim = MDSimulation(mol, dt_fs=0.5, temperature_K=300.0,
                           record_interval=10)
        traj = sim.run(100)
        assert traj.n_frames >= 10

    def test_interpolation_after_run(self):
        """Should be able to interpolate positions after simulation."""
        from molbuilder.dynamics.simulation import MDSimulation
        mol = Molecule("H2")
        mol.add_atom("H", [0.0, 0.0, 0.0])
        mol.add_atom("H", [0.74, 0.0, 0.0])
        mol.add_bond(0, 1, order=1)

        sim = MDSimulation(mol, dt_fs=0.5, temperature_K=100.0)
        traj = sim.run(20)

        # Interpolate at midpoint
        t_mid = (traj.t_start + traj.t_end) / 2
        pos = traj.at_time(t_mid)
        assert pos.shape == (2, 3)
        assert np.all(np.isfinite(pos))


class TestPhysicalValidation:
    """Test physical validity of simulation outputs."""

    def test_cc_vibration_period(self):
        """C-C stretch vibration period should be approximately 30-50 fs.

        A C-C single bond has a stretching frequency of ~1000 cm^-1,
        corresponding to a period of ~33 fs.  We test that the
        dominant vibration period is in a reasonable range.
        """
        from molbuilder.dynamics.simulation import MDSimulation

        mol = Molecule("CC diatomic")
        mol.add_atom("C", [0.0, 0.0, 0.0], Hybridization.SP3)
        mol.add_atom("C", [1.54, 0.0, 0.0], Hybridization.SP3)
        mol.add_bond(0, 1, order=1)

        sim = MDSimulation(mol, dt_fs=0.1, temperature_K=100.0,
                           thermostat=False, record_interval=1)
        traj = sim.run(1000)

        # Extract C-C distances over time
        distances = []
        for frame in traj.frames:
            d = np.linalg.norm(frame.positions[1] - frame.positions[0])
            distances.append(d)

        distances = np.array(distances)
        times = traj.times

        # Find period from zero crossings of (d - mean)
        mean_d = np.mean(distances)
        deviations = distances - mean_d
        crossings = []
        for i in range(1, len(deviations)):
            if deviations[i-1] * deviations[i] < 0:
                # Linear interpolation for crossing time
                frac = abs(deviations[i-1]) / (
                    abs(deviations[i-1]) + abs(deviations[i]))
                t_cross = times[i-1] + frac * (times[i] - times[i-1])
                crossings.append(t_cross)

        if len(crossings) >= 4:
            # Period = 2 * (average half-period)
            half_periods = np.diff(crossings)
            period = 2 * np.mean(half_periods)
            # C-C stretch should be roughly 20-80 fs
            assert 15.0 < period < 100.0, (
                f"C-C vibration period {period:.1f} fs out of expected range")
