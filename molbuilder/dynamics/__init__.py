"""Molecular dynamics engine for atomic interaction visualization.

Provides a force field, Velocity Verlet integrator, trajectory storage
with sub-femtosecond interpolation, reaction mechanism templates, and
steered-MD choreography for animating bond formation/breaking events.

Public API
----------
ForceField, ForceFieldParams, ForceResult
    Vectorized force computation (LJ + Coulomb + harmonic bond/angle +
    OPLS-AA torsion).

VelocityVerletIntegrator, IntegratorState, BerendsenThermostat
    Symplectic integration and temperature control.

Trajectory, TrajectoryFrame
    Frame storage with CubicSpline interpolation for sub-fs resolution.

MDSimulation
    High-level orchestrator: build from a ``Molecule``, run plain or
    steered MD, and record trajectories.

MechanismType, ReactionMechanism, MechanismStage, ElectronFlow
    Reaction mechanism data model and predefined templates.

MechanismChoreographer
    Converts mechanism stages into time-varying restraint forces.
"""

from molbuilder.dynamics.forcefield import ForceField, ForceFieldParams, ForceResult
from molbuilder.dynamics.integrator import (
    VelocityVerletIntegrator,
    IntegratorState,
    BerendsenThermostat,
)
from molbuilder.dynamics.trajectory import Trajectory, TrajectoryFrame
from molbuilder.dynamics.simulation import MDSimulation
from molbuilder.dynamics.mechanisms import (
    MechanismType,
    ReactionMechanism,
    MechanismStage,
    ElectronFlow,
    FlowType,
    sn2_mechanism,
    e2_mechanism,
    radical_substitution_mechanism,
    nucleophilic_addition_mechanism,
)
from molbuilder.dynamics.mechanism_choreography import MechanismChoreographer
