"""Steered-MD choreography for reaction mechanisms.

Converts MechanismStage targets (distance, bond order, angle) into
time-varying harmonic restraint forces that guide the MD simulation
smoothly through each stage of the mechanism.  A sigmoid ramp function
ensures gradual transitions.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from molbuilder.dynamics.mechanisms import (
    ReactionMechanism,
    MechanismStage,
    ElectronFlow,
)

if TYPE_CHECKING:
    from molbuilder.dynamics.forcefield import ForceField


def _sigmoid(x: float) -> float:
    """Smooth sigmoid ramp: 0 at x=0, ~1 at x=1."""
    t = np.clip(x, 0.0, 1.0)
    return float(1.0 / (1.0 + np.exp(-10.0 * (t - 0.5))))


class MechanismChoreographer:
    """Converts mechanism stages into restraint forces for steered MD.

    For each stage, harmonic distance restraints with a sigmoid ramp
    are applied to guide atoms toward their target distances, producing
    smooth bond formation and breaking events.

    Parameters
    ----------
    mechanism : ReactionMechanism
        The reaction mechanism template.
    forcefield : ForceField
        The underlying force field (used for reference).
    n_steps_per_stage : int
        Number of MD steps per mechanism stage.
    restraint_k : float
        Restraint force constant in kJ/(mol*A^2).  Higher values
        make the steering stronger.
    """

    def __init__(self, mechanism: ReactionMechanism,
                 forcefield: ForceField,
                 n_steps_per_stage: int = 200,
                 restraint_k: float = 500.0):
        self.mechanism = mechanism
        self.ff = forcefield
        self.n_steps_per_stage = n_steps_per_stage
        self.restraint_k = restraint_k

        # Precompute initial distances for each target
        self._initial_distances: dict[int, dict[tuple[int, int], float]] = {}

    def _get_initial_distance(self, positions: np.ndarray,
                               i: int, j: int) -> float:
        """Compute distance between atoms i and j."""
        return float(np.linalg.norm(positions[j] - positions[i]))

    def restraint_forces(self, positions: np.ndarray,
                          stage_idx: int,
                          progress: float) -> np.ndarray:
        """Compute restraint forces for the current stage and progress.

        Parameters
        ----------
        positions : ndarray of shape (n_atoms, 3)
            Current atomic positions.
        stage_idx : int
            Index of the current mechanism stage.
        progress : float
            Progress through the current stage (0.0 to 1.0).

        Returns
        -------
        ndarray of shape (n_atoms, 3)
            Restraint forces in kJ/(mol*A).
        """
        forces = np.zeros_like(positions)
        stage = self.mechanism.stages[stage_idx]
        ramp = _sigmoid(progress)

        # Distance restraints
        for (i, j), target_dist in stage.distance_targets.items():
            rij = positions[j] - positions[i]
            r = np.linalg.norm(rij)
            if r < 1e-12:
                continue

            # Interpolate target distance with sigmoid ramp
            current_target = target_dist  # could interpolate from previous
            dr = r - current_target
            f_mag = -self.restraint_k * ramp * dr
            f_vec = f_mag * (rij / r)

            forces[i] -= f_vec
            forces[j] += f_vec

        # Angle restraints
        for (i, j, k), target_angle_deg in stage.angle_targets.items():
            rji = positions[i] - positions[j]
            rjk = positions[k] - positions[j]
            nji = np.linalg.norm(rji)
            njk = np.linalg.norm(rjk)
            if nji < 1e-12 or njk < 1e-12:
                continue

            cos_theta = np.clip(
                np.dot(rji, rjk) / (nji * njk), -1.0, 1.0)
            theta = math.acos(cos_theta)
            target_rad = math.radians(target_angle_deg)
            d_theta = theta - target_rad

            sin_theta = math.sin(theta)
            if abs(sin_theta) < 1e-12:
                continue

            # Angle restraint force constant (kJ/(mol*rad^2))
            angle_k = self.restraint_k * 0.5 * ramp
            dE_dtheta = angle_k * d_theta

            rji_hat = rji / nji
            rjk_hat = rjk / njk
            fi = (dE_dtheta / (nji * sin_theta)) * (
                cos_theta * rji_hat - rjk_hat)
            fk = (dE_dtheta / (njk * sin_theta)) * (
                cos_theta * rjk_hat - rji_hat)

            forces[i] -= fi
            forces[k] -= fk
            forces[j] += fi + fk

        return forces

    def bond_orders_at(self, stage_idx: int,
                        progress: float) -> dict[tuple[int, int], float]:
        """Compute interpolated bond orders at the current stage/progress.

        Parameters
        ----------
        stage_idx : int
            Index of the current mechanism stage.
        progress : float
            Progress through the current stage (0.0 to 1.0).

        Returns
        -------
        dict[tuple[int, int], float]
            Bond order for each atom pair that has a target.
        """
        stage = self.mechanism.stages[stage_idx]
        ramp = _sigmoid(progress)

        # Get previous stage bond orders as starting point
        prev_orders: dict[tuple[int, int], float] = {}
        if stage_idx > 0:
            prev_stage = self.mechanism.stages[stage_idx - 1]
            prev_orders = dict(prev_stage.bond_order_changes)

        result: dict[tuple[int, int], float] = {}
        for key, target_order in stage.bond_order_changes.items():
            normalized_key = (min(key), max(key))
            prev = prev_orders.get(key, 0.0)
            result[normalized_key] = prev + ramp * (target_order - prev)

        return result

    def electron_flows_at(self, stage_idx: int,
                           progress: float) -> list[ElectronFlow]:
        """Return the electron flows active at the current stage/progress.

        Flows are returned when progress is between 0.1 and 0.9 to
        avoid showing arrows at the very start/end of a stage.

        Parameters
        ----------
        stage_idx : int
            Index of the current mechanism stage.
        progress : float
            Progress through the current stage.

        Returns
        -------
        list[ElectronFlow]
            Active electron flow arrows.
        """
        if progress < 0.1 or progress > 0.9:
            return []
        stage = self.mechanism.stages[stage_idx]
        return list(stage.electron_flows)

    def stage_annotation(self, stage_idx: int) -> str:
        """Return the text annotation for the current stage.

        Parameters
        ----------
        stage_idx : int
            Index of the current mechanism stage.

        Returns
        -------
        str
            Annotation text.
        """
        if 0 <= stage_idx < len(self.mechanism.stages):
            return self.mechanism.stages[stage_idx].annotation
        return ""
