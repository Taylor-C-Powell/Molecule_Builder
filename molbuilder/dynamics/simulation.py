"""MDSimulation orchestrator.

High-level class that wires together ForceField, VelocityVerletIntegrator,
and Trajectory to run plain or steered molecular dynamics from a Molecule.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from molbuilder.dynamics.forcefield import ForceField
from molbuilder.dynamics.integrator import (
    VelocityVerletIntegrator,
    IntegratorState,
    BerendsenThermostat,
)
from molbuilder.dynamics.trajectory import Trajectory, TrajectoryFrame

if TYPE_CHECKING:
    from molbuilder.molecule.graph import Molecule
    from molbuilder.dynamics.mechanisms import ReactionMechanism


class MDSimulation:
    """Molecular dynamics simulation orchestrator.

    Builds a force field from a Molecule, initializes velocities at a
    target temperature, and runs either plain or mechanism-steered MD,
    recording every frame for visualization.

    Parameters
    ----------
    molecule : Molecule
        Starting molecular structure.
    dt_fs : float
        Integration timestep in femtoseconds.
    temperature_K : float
        Target temperature in Kelvin.
    thermostat : bool
        Whether to enable Berendsen thermostat.
    thermostat_tau_fs : float
        Coupling time for the thermostat.
    record_interval : int
        Record a trajectory frame every N steps.
    """

    def __init__(self, molecule: Molecule,
                 dt_fs: float = 0.5,
                 temperature_K: float = 300.0,
                 thermostat: bool = True,
                 thermostat_tau_fs: float = 100.0,
                 record_interval: int = 1):
        self.molecule = molecule
        self.dt_fs = dt_fs
        self.temperature_K = temperature_K
        self.record_interval = record_interval

        # Build force field
        self.ff = ForceField.from_molecule(molecule)

        # Set up thermostat
        thermo = None
        if thermostat:
            thermo = BerendsenThermostat(temperature_K, thermostat_tau_fs)

        # Set up integrator
        self.integrator = VelocityVerletIntegrator(
            dt_fs=dt_fs,
            masses=self.ff.params.masses,
            thermostat=thermo,
        )

        # Initialize state from molecule positions
        positions = np.array([a.position for a in molecule.atoms], dtype=float)
        velocities = self.integrator.initialize_velocities(
            len(molecule.atoms), temperature_K)

        self.state = IntegratorState(
            positions=positions,
            velocities=velocities,
            time_fs=0.0,
        )

    def run(self, n_steps: int) -> Trajectory:
        """Run plain molecular dynamics for n_steps.

        Parameters
        ----------
        n_steps : int
            Number of integration steps.

        Returns
        -------
        Trajectory
            Recorded trajectory with frames at record_interval spacing.
        """
        n_atoms = len(self.molecule.atoms)
        traj = Trajectory(n_atoms)

        def force_fn(pos):
            return self.ff.compute(pos)

        # Record initial frame
        result = force_fn(self.state.positions)
        ke = self.integrator.kinetic_energy(self.state.velocities)
        traj.add_frame(TrajectoryFrame(
            time_fs=self.state.time_fs,
            positions=self.state.positions.copy(),
            velocities=self.state.velocities.copy(),
            energy_kinetic=ke,
            energy_potential=result.energy_total,
        ))

        for step in range(1, n_steps + 1):
            self.state = self.integrator.step(self.state, force_fn)

            if step % self.record_interval == 0:
                result = force_fn(self.state.positions)
                ke = self.integrator.kinetic_energy(self.state.velocities)
                traj.add_frame(TrajectoryFrame(
                    time_fs=self.state.time_fs,
                    positions=self.state.positions.copy(),
                    velocities=self.state.velocities.copy(),
                    energy_kinetic=ke,
                    energy_potential=result.energy_total,
                ))

        return traj

    def run_mechanism(self, mechanism: ReactionMechanism,
                       n_steps_per_stage: int = 200) -> Trajectory:
        """Run steered MD driven by a reaction mechanism.

        The MechanismChoreographer applies time-varying restraint forces
        to guide atoms through the mechanism stages, producing smooth
        bond formation/breaking events.

        Parameters
        ----------
        mechanism : ReactionMechanism
            Reaction mechanism template.
        n_steps_per_stage : int
            Number of MD steps per mechanism stage.

        Returns
        -------
        Trajectory
            Recorded trajectory with fractional bond orders.
        """
        from molbuilder.dynamics.mechanism_choreography import MechanismChoreographer

        n_atoms = len(self.molecule.atoms)
        traj = Trajectory(n_atoms)
        choreographer = MechanismChoreographer(
            mechanism, self.ff, n_steps_per_stage)

        n_stages = len(mechanism.stages)
        total_steps = n_stages * n_steps_per_stage

        for stage_idx in range(n_stages):
            for step_in_stage in range(n_steps_per_stage):
                progress = step_in_stage / max(1, n_steps_per_stage - 1)
                global_step = stage_idx * n_steps_per_stage + step_in_stage

                def force_fn(pos, _si=stage_idx, _pr=progress):
                    ff_result = self.ff.compute(pos)
                    restraint = choreographer.restraint_forces(
                        pos, _si, _pr)
                    ff_result.forces = ff_result.forces + restraint
                    return ff_result

                self.state = self.integrator.step(self.state, force_fn)

                if global_step % self.record_interval == 0:
                    result = self.ff.compute(self.state.positions)
                    ke = self.integrator.kinetic_energy(self.state.velocities)
                    bond_orders = choreographer.bond_orders_at(
                        stage_idx, progress)
                    traj.add_frame(TrajectoryFrame(
                        time_fs=self.state.time_fs,
                        positions=self.state.positions.copy(),
                        velocities=self.state.velocities.copy(),
                        energy_kinetic=ke,
                        energy_potential=result.energy_total,
                        bond_orders=bond_orders,
                    ))

        return traj

    def get_positions_at_time(self, t_fs: float,
                               trajectory: Trajectory) -> np.ndarray:
        """Convenience interpolation of positions from a trajectory.

        Parameters
        ----------
        t_fs : float
            Time in femtoseconds.
        trajectory : Trajectory
            A previously recorded trajectory.

        Returns
        -------
        ndarray of shape (n_atoms, 3)
            Interpolated positions.
        """
        return trajectory.at_time(t_fs)
