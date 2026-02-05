"""Main animation renderer for extreme slow-motion atomic interactions.

Provides ``InteractionVisualizer`` which builds matplotlib FuncAnimation
pipelines that render MD trajectories at sub-femtosecond resolution,
with optional overlays for electron density, curly arrows, energy bars,
and time labels.

Convenience functions:
    - ``visualize_md_trajectory(mol, n_steps, ...)``
    - ``visualize_reaction(mol, mechanism, ...)``
    - ``visualize_bond_formation(mol, atom_i, atom_j, ...)``
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

try:
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    from matplotlib.patches import FancyArrowPatch
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from molbuilder.core.element_properties import cpk_color, covalent_radius_pm
from molbuilder.core.elements import SYMBOL_TO_Z
from molbuilder.visualization.theme import (
    BG_COLOR, TEXT_COLOR, BOND_COLOR,
    FORMING_BOND_COLOR, BREAKING_BOND_COLOR,
    ELECTRON_CLOUD_COLOR, ARROW_COLOR, TRANSITION_STATE_COLOR,
)

if TYPE_CHECKING:
    from molbuilder.molecule.graph import Molecule
    from molbuilder.dynamics.trajectory import Trajectory
    from molbuilder.dynamics.mechanisms import ReactionMechanism, ElectronFlow
    from molbuilder.dynamics.mechanism_choreography import MechanismChoreographer


# ===================================================================
# Configuration
# ===================================================================

@dataclass
class PlaybackConfig:
    """Configuration for the interaction visualizer.

    Attributes
    ----------
    slowmo_factor : float
        Slowdown factor.  Default 1e15 means 1 fs of simulation time
        maps to 1 second of animation (extreme slow motion).
    fps : int
        Animation frame rate.
    show_electron_density : bool
        Render electron density clouds during bond events.
    show_electron_flows : bool
        Render curly/fishhook arrows for electron flow.
    show_energy_bar : bool
        Show an energy bar overlay.
    show_time_label : bool
        Show a time label with SI prefix.
    show_bond_orders : bool
        Show fractional bond orders as annotations.
    export_path : str or None
        If set, export animation to this file path (.mp4 or .gif).
    camera_follow : bool
        Auto-track the reaction center.
    dpi : int
        Resolution for export.
    figsize : tuple[float, float]
        Figure size in inches.
    """
    slowmo_factor: float = 1e15
    fps: int = 30
    show_electron_density: bool = True
    show_electron_flows: bool = True
    show_energy_bar: bool = True
    show_time_label: bool = True
    show_bond_orders: bool = False
    export_path: str | None = None
    camera_follow: bool = False
    dpi: int = 150
    figsize: tuple[float, float] = (10, 8)


def _time_label(t_fs: float) -> str:
    """Format time with appropriate SI prefix."""
    if abs(t_fs) < 1e-3:
        return f"{t_fs * 1e3:.2f} as"
    if abs(t_fs) < 1.0:
        return f"{t_fs:.3f} fs"
    if abs(t_fs) < 1000.0:
        return f"{t_fs:.1f} fs"
    return f"{t_fs / 1000.0:.2f} ps"


def _bond_color_from_order(order: float) -> str:
    """Choose bond color based on fractional bond order."""
    if order < 0.3:
        return BREAKING_BOND_COLOR
    if order > 0.7 and order < 1.3:
        return BOND_COLOR
    if order >= 0.3 and order < 0.7:
        return TRANSITION_STATE_COLOR
    return FORMING_BOND_COLOR


# ===================================================================
# InteractionVisualizer
# ===================================================================

class InteractionVisualizer:
    """Renders MD trajectories as slow-motion matplotlib animations.

    Parameters
    ----------
    trajectory : Trajectory
        The MD trajectory to visualize.
    molecule : Molecule
        Original molecule (for bond connectivity and atom symbols).
    mechanism : ReactionMechanism or None
        Optional mechanism for electron flow overlays.
    choreographer : MechanismChoreographer or None
        Optional choreographer for mechanism annotations.
    config : PlaybackConfig or None
        Playback configuration.
    """

    def __init__(self, trajectory, molecule,
                 mechanism=None, choreographer=None,
                 config: PlaybackConfig | None = None):
        if not HAS_MATPLOTLIB:
            raise ImportError("matplotlib is required for visualization")

        self.trajectory = trajectory
        self.molecule = molecule
        self.mechanism = mechanism
        self.choreographer = choreographer
        self.config = config or PlaybackConfig()

        # Precompute animation frames
        duration_fs = trajectory.duration
        duration_anim_s = duration_fs * self.config.slowmo_factor * 1e-15
        # Minimum 2 seconds of animation
        duration_anim_s = max(duration_anim_s, 2.0)
        self.n_frames = max(int(duration_anim_s * self.config.fps), 2)

        # Time array for each animation frame
        self.frame_times = np.linspace(
            trajectory.t_start, trajectory.t_end, self.n_frames)

        # Electron density renderer
        self._edensity = None
        if self.config.show_electron_density:
            from molbuilder.visualization.electron_density_viz import (
                ElectronDensityRenderer,
            )
            self._edensity = ElectronDensityRenderer(n_points=2000)

        self._fig = None
        self._ax = None
        self._anim = None
        self._paused = False

    def build_animation(self) -> FuncAnimation:
        """Build and return the matplotlib FuncAnimation.

        Returns
        -------
        FuncAnimation
            The animation object.  Call ``plt.show()`` to display, or
            use ``export()`` to save.
        """
        self._fig, self._ax = plt.subplots(
            figsize=self.config.figsize,
            facecolor=BG_COLOR,
        )
        self._ax.set_facecolor(BG_COLOR)
        self._ax.set_aspect("equal")
        self._ax.tick_params(colors=TEXT_COLOR)
        for spine in self._ax.spines.values():
            spine.set_color(TEXT_COLOR)

        self._anim = FuncAnimation(
            self._fig,
            self._render_frame,
            frames=self.n_frames,
            interval=1000 // self.config.fps,
            blit=False,
            repeat=True,
        )

        return self._anim

    def _render_frame(self, frame_idx: int):
        """Render a single animation frame."""
        ax = self._ax
        ax.clear()
        ax.set_facecolor(BG_COLOR)

        t_fs = self.frame_times[frame_idx]
        positions = self.trajectory.at_time(t_fs)

        mol = self.molecule
        symbols = [a.symbol for a in mol.atoms]

        # Compute bounds for axis
        margin = 2.0
        x_min = positions[:, 0].min() - margin
        x_max = positions[:, 0].max() + margin
        y_min = positions[:, 1].min() - margin
        y_max = positions[:, 1].max() + margin
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_aspect("equal")

        # Draw bonds
        for bond in mol.bonds:
            i, j = bond.atom_i, bond.atom_j
            key = (min(i, j), max(i, j))
            order = float(bond.order)

            # Check for fractional bond order from trajectory
            frac_order = self.trajectory.bond_order_at_time(t_fs, i, j)
            if frac_order > 0.01:
                order = frac_order

            # Bond styling
            color = BOND_COLOR
            linestyle = "-"
            linewidth = 1.5 * max(order, 0.2)

            if frac_order > 0.01:
                color = _bond_color_from_order(frac_order)
                if frac_order < 0.5:
                    linestyle = ":"
                elif frac_order < 0.8:
                    linestyle = "--"

            pi = positions[i]
            pj = positions[j]
            ax.plot(
                [pi[0], pj[0]], [pi[1], pj[1]],
                color=color, linewidth=linewidth,
                linestyle=linestyle, zorder=1,
            )

            # Bond order annotation
            if self.config.show_bond_orders and frac_order > 0.05:
                mid = 0.5 * (pi + pj)
                ax.text(mid[0], mid[1] + 0.15, f"{frac_order:.2f}",
                        color=TEXT_COLOR, fontsize=7, ha="center",
                        va="bottom", zorder=5)

        # Draw electron density
        if self._edensity is not None and self.config.show_electron_density:
            for bond in mol.bonds:
                i, j = bond.atom_i, bond.atom_j
                frac_order = self.trajectory.bond_order_at_time(t_fs, i, j)
                if 0.1 < frac_order < 0.9:
                    z_a = SYMBOL_TO_Z.get(symbols[i], 1)
                    z_b = SYMBOL_TO_Z.get(symbols[j], 1)
                    self._edensity.render_on_axis_2d(
                        ax, positions[i], positions[j],
                        z_a, z_b, frac_order, point_size=1.5)

        # Draw electron flow arrows
        if (self.config.show_electron_flows
                and self.mechanism is not None
                and self.choreographer is not None):
            self._draw_electron_flows(ax, positions, t_fs, frame_idx)

        # Draw atoms (on top of bonds)
        for idx, atom in enumerate(mol.atoms):
            x, y = positions[idx, 0], positions[idx, 1]
            color = cpk_color(atom.symbol)
            radius = covalent_radius_pm(atom.symbol) / 100.0 * 0.3
            radius = max(radius, 0.15)

            circle = plt.Circle(
                (x, y), radius, color=color,
                zorder=3, ec="#222222", linewidth=0.5)
            ax.add_patch(circle)

            # Atom label for non-hydrogen atoms
            if atom.symbol != "H":
                ax.text(x, y, atom.symbol,
                        color="#000000", fontsize=7, fontweight="bold",
                        ha="center", va="center", zorder=4)

        # Overlays
        if self.config.show_time_label:
            ax.text(0.02, 0.98, _time_label(t_fs),
                    transform=ax.transAxes,
                    color=TEXT_COLOR, fontsize=12, fontweight="bold",
                    va="top", ha="left",
                    bbox=dict(boxstyle="round,pad=0.3",
                              facecolor=BG_COLOR, alpha=0.8))

        if self.config.show_energy_bar:
            self._draw_energy_bar(ax, t_fs)

        # Mechanism annotation
        if self.choreographer is not None:
            stage_idx = self._stage_for_time(t_fs)
            if stage_idx is not None:
                annotation = self.choreographer.stage_annotation(stage_idx)
                if annotation:
                    ax.text(0.5, 0.02, annotation,
                            transform=ax.transAxes,
                            color=TRANSITION_STATE_COLOR,
                            fontsize=10, ha="center", va="bottom",
                            bbox=dict(boxstyle="round,pad=0.3",
                                      facecolor=BG_COLOR, alpha=0.8))

        ax.set_xlabel("x (A)", color=TEXT_COLOR, fontsize=8)
        ax.set_ylabel("y (A)", color=TEXT_COLOR, fontsize=8)

    def _draw_electron_flows(self, ax, positions, t_fs, frame_idx):
        """Draw curly/fishhook arrows for electron flows."""
        from molbuilder.dynamics.mechanisms import FlowType

        stage_idx = self._stage_for_time(t_fs)
        if stage_idx is None:
            return

        total_steps = len(self.mechanism.stages)
        progress = (frame_idx / max(1, self.n_frames - 1))
        stage_progress = (progress * total_steps) - stage_idx

        flows = self.choreographer.electron_flows_at(
            stage_idx, stage_progress)

        for flow in flows:
            start = positions[flow.from_atom][:2]
            mid_bond = 0.5 * (
                positions[flow.to_bond[0]][:2]
                + positions[flow.to_bond[1]][:2])

            color = ARROW_COLOR
            style = "Simple,tail_width=2,head_width=8,head_length=6"
            if flow.flow_type == FlowType.FISHHOOK_ARROW:
                style = "Simple,tail_width=1,head_width=5,head_length=4"

            arrow = FancyArrowPatch(
                posA=tuple(start), posB=tuple(mid_bond),
                arrowstyle=style,
                color=color, alpha=0.8,
                connectionstyle="arc3,rad=0.3",
                zorder=6,
            )
            ax.add_patch(arrow)

            if flow.label:
                label_pos = 0.5 * (start + mid_bond)
                ax.text(label_pos[0], label_pos[1] + 0.3,
                        flow.label, color=ARROW_COLOR,
                        fontsize=6, ha="center", va="bottom",
                        zorder=7)

    def _draw_energy_bar(self, ax, t_fs):
        """Draw a small energy bar overlay."""
        times, ke, pe = self.trajectory.energies()
        if len(times) < 2:
            return

        # Find closest frame
        idx = np.searchsorted(times, t_fs)
        idx = min(idx, len(times) - 1)
        ke_val = ke[idx]
        pe_val = pe[idx]
        total = ke_val + pe_val

        bar_text = f"KE: {ke_val:.1f}  PE: {pe_val:.1f}  Total: {total:.1f} kJ/mol"
        ax.text(0.98, 0.98, bar_text,
                transform=ax.transAxes,
                color=TEXT_COLOR, fontsize=8,
                va="top", ha="right",
                bbox=dict(boxstyle="round,pad=0.3",
                          facecolor=BG_COLOR, alpha=0.8))

    def _stage_for_time(self, t_fs):
        """Determine which mechanism stage corresponds to a time."""
        if self.mechanism is None:
            return None
        n_stages = len(self.mechanism.stages)
        if n_stages == 0:
            return None
        frac = (t_fs - self.trajectory.t_start) / max(
            self.trajectory.duration, 1e-10)
        stage = int(frac * n_stages)
        return min(stage, n_stages - 1)

    def show(self):
        """Build and display the animation interactively."""
        self.build_animation()
        plt.show()

    def export(self, path: str | None = None):
        """Export the animation to a file.

        Parameters
        ----------
        path : str or None
            Output path.  If None, uses config.export_path.
            Supports .mp4 (requires ffmpeg) and .gif (Pillow).
        """
        path = path or self.config.export_path
        if path is None:
            raise ValueError("No export path specified")

        if self._anim is None:
            self.build_animation()

        if path.endswith(".mp4"):
            try:
                from matplotlib.animation import FFMpegWriter
                writer = FFMpegWriter(fps=self.config.fps)
                self._anim.save(path, writer=writer, dpi=self.config.dpi)
            except Exception:
                # Fallback to gif
                path = path.replace(".mp4", ".gif")
                from matplotlib.animation import PillowWriter
                writer = PillowWriter(fps=self.config.fps)
                self._anim.save(path, writer=writer, dpi=self.config.dpi)
        elif path.endswith(".gif"):
            from matplotlib.animation import PillowWriter
            writer = PillowWriter(fps=self.config.fps)
            self._anim.save(path, writer=writer, dpi=self.config.dpi)
        else:
            self._anim.save(path, dpi=self.config.dpi)

        plt.close(self._fig)
        return path

    @property
    def fig(self):
        return self._fig


# ===================================================================
# Convenience functions
# ===================================================================

def visualize_md_trajectory(molecule,
                            n_steps: int = 500,
                            dt_fs: float = 0.5,
                            temperature_K: float = 300.0,
                            config: PlaybackConfig | None = None,
                            show: bool = True):
    """Quick MD simulation + visualization.

    Parameters
    ----------
    molecule : Molecule
        Molecule to simulate.
    n_steps : int
        Number of MD steps.
    dt_fs : float
        Timestep in fs.
    temperature_K : float
        Target temperature.
    config : PlaybackConfig or None
        Visualization config.
    show : bool
        If True, display interactively.

    Returns
    -------
    InteractionVisualizer
        The visualizer (can be used for export).
    """
    from molbuilder.dynamics.simulation import MDSimulation

    sim = MDSimulation(molecule, dt_fs=dt_fs,
                       temperature_K=temperature_K)
    traj = sim.run(n_steps)

    viz = InteractionVisualizer(traj, molecule, config=config)
    if show:
        viz.show()
    return viz


def visualize_reaction(molecule,
                       mechanism,
                       n_steps_per_stage: int = 200,
                       dt_fs: float = 0.5,
                       temperature_K: float = 50.0,
                       config: PlaybackConfig | None = None,
                       show: bool = True):
    """Mechanism-steered MD + visualization.

    Parameters
    ----------
    molecule : Molecule
        Starting molecule.
    mechanism : ReactionMechanism
        Reaction mechanism template.
    n_steps_per_stage : int
        MD steps per mechanism stage.
    dt_fs : float
        Timestep.
    temperature_K : float
        Temperature (low for cleaner mechanism).
    config : PlaybackConfig or None
        Visualization config.
    show : bool
        Display interactively.

    Returns
    -------
    InteractionVisualizer
    """
    from molbuilder.dynamics.simulation import MDSimulation
    from molbuilder.dynamics.mechanism_choreography import MechanismChoreographer
    from molbuilder.dynamics.forcefield import ForceField

    sim = MDSimulation(molecule, dt_fs=dt_fs,
                       temperature_K=temperature_K,
                       thermostat=True,
                       thermostat_tau_fs=50.0)
    traj = sim.run_mechanism(mechanism, n_steps_per_stage)

    choreographer = MechanismChoreographer(
        mechanism, sim.ff, n_steps_per_stage)
    viz = InteractionVisualizer(
        traj, molecule,
        mechanism=mechanism,
        choreographer=choreographer,
        config=config)
    if show:
        viz.show()
    return viz


def visualize_bond_formation(molecule,
                             atom_i: int, atom_j: int,
                             n_steps: int = 400,
                             config: PlaybackConfig | None = None,
                             show: bool = True):
    """Visualize a single bond formation event.

    Creates a simple mechanism that brings two atoms together and
    forms a bond.

    Parameters
    ----------
    molecule : Molecule
        Molecule containing both atoms.
    atom_i, atom_j : int
        Atom indices for the bond to form.
    n_steps : int
        Total MD steps.
    config : PlaybackConfig or None
        Visualization config.
    show : bool
        Display interactively.

    Returns
    -------
    InteractionVisualizer
    """
    from molbuilder.dynamics.mechanisms import (
        ReactionMechanism, MechanismStage, MechanismType,
        ElectronFlow, FlowType,
    )
    from molbuilder.core.bond_data import bond_length

    sym_i = molecule.atoms[atom_i].symbol
    sym_j = molecule.atoms[atom_j].symbol
    target_bl = bond_length(sym_i, sym_j, 1)

    key = (min(atom_i, atom_j), max(atom_i, atom_j))
    mechanism = ReactionMechanism(
        name=f"Bond formation {sym_i}-{sym_j}",
        mechanism_type=MechanismType.SN2,
        stages=[
            MechanismStage(
                name="Approach",
                distance_targets={(atom_i, atom_j): target_bl * 1.5},
                bond_order_changes={key: 0.3},
                electron_flows=[
                    ElectronFlow(atom_i, (atom_i, atom_j),
                                 FlowType.CURLY_ARROW,
                                 "Bond forming"),
                ],
                duration_weight=1.0,
                annotation=f"Atoms approaching: {sym_i}...{sym_j}",
            ),
            MechanismStage(
                name="Bond formation",
                distance_targets={(atom_i, atom_j): target_bl},
                bond_order_changes={key: 1.0},
                electron_flows=[],
                duration_weight=1.5,
                annotation=f"Bond formed: {sym_i}-{sym_j}",
            ),
        ],
        atom_roles={"atom_i": atom_i, "atom_j": atom_j},
    )

    return visualize_reaction(
        molecule, mechanism,
        n_steps_per_stage=n_steps // 2,
        temperature_K=50.0,
        config=config,
        show=show)
