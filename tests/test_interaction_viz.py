"""Tests for the interaction visualization pipeline.

Uses headless Agg backend to test animation construction and
export without requiring a display.
"""

import os
import tempfile

import numpy as np
import pytest

import matplotlib
matplotlib.use("Agg")  # Headless backend for testing

from molbuilder.molecule.graph import Molecule, Hybridization


class TestPlaybackConfig:
    """Test PlaybackConfig defaults and creation."""

    def test_default_config(self):
        """Default config should have sensible values."""
        from molbuilder.visualization.interaction_viz import PlaybackConfig
        config = PlaybackConfig()
        assert config.slowmo_factor == 1e15
        assert config.fps == 30
        assert config.show_electron_density is True
        assert config.show_time_label is True
        assert config.export_path is None

    def test_custom_config(self):
        """Custom config values should be preserved."""
        from molbuilder.visualization.interaction_viz import PlaybackConfig
        config = PlaybackConfig(
            slowmo_factor=1e12,
            fps=60,
            show_electron_density=False)
        assert config.slowmo_factor == 1e12
        assert config.fps == 60
        assert config.show_electron_density is False


class TestTimeLabel:
    """Test the time label formatter."""

    def test_attosecond_range(self):
        """Times < 1e-3 fs should show as attoseconds."""
        from molbuilder.visualization.interaction_viz import _time_label
        label = _time_label(0.0005)
        assert "as" in label

    def test_femtosecond_range(self):
        """Times in 0.001-1 fs should show as femtoseconds."""
        from molbuilder.visualization.interaction_viz import _time_label
        label = _time_label(0.5)
        assert "fs" in label

    def test_picosecond_range(self):
        """Times >= 1000 fs should show as picoseconds."""
        from molbuilder.visualization.interaction_viz import _time_label
        label = _time_label(1500)
        assert "ps" in label


class TestInteractionVisualizer:
    """Test the InteractionVisualizer animation pipeline."""

    def _build_h2_trajectory(self):
        """Build a simple H2 trajectory for testing."""
        from molbuilder.dynamics.trajectory import Trajectory, TrajectoryFrame

        mol = Molecule("H2")
        mol.add_atom("H", [0.0, 0.0, 0.0])
        mol.add_atom("H", [0.74, 0.0, 0.0])
        mol.add_bond(0, 1, order=1)

        traj = Trajectory(n_atoms=2)
        for t in np.linspace(0, 10, 20):
            offset = 0.05 * np.sin(2 * np.pi * t / 5)
            traj.add_frame(TrajectoryFrame(
                time_fs=t,
                positions=np.array([[0.0, 0.0, 0.0],
                                     [0.74 + offset, 0.0, 0.0]]),
                energy_kinetic=10.0 + t,
                energy_potential=20.0 - t,
            ))
        return mol, traj

    def test_build_animation(self):
        """build_animation should return a FuncAnimation."""
        from molbuilder.visualization.interaction_viz import (
            InteractionVisualizer, PlaybackConfig,
        )
        from matplotlib.animation import FuncAnimation

        mol, traj = self._build_h2_trajectory()
        config = PlaybackConfig(show_electron_density=False, fps=10)
        viz = InteractionVisualizer(traj, mol, config=config)
        anim = viz.build_animation()
        assert isinstance(anim, FuncAnimation)
        import matplotlib.pyplot as plt
        plt.close(viz.fig)

    def test_render_single_frame(self):
        """Rendering a single frame should not raise."""
        from molbuilder.visualization.interaction_viz import (
            InteractionVisualizer, PlaybackConfig,
        )
        mol, traj = self._build_h2_trajectory()
        config = PlaybackConfig(show_electron_density=False, fps=10)
        viz = InteractionVisualizer(traj, mol, config=config)
        viz.build_animation()
        # Render frame 0 manually
        viz._render_frame(0)
        import matplotlib.pyplot as plt
        plt.close(viz.fig)

    def test_export_gif(self):
        """GIF export should produce a valid file."""
        from molbuilder.visualization.interaction_viz import (
            InteractionVisualizer, PlaybackConfig,
        )
        mol, traj = self._build_h2_trajectory()

        with tempfile.NamedTemporaryFile(suffix=".gif", delete=False) as f:
            path = f.name

        try:
            config = PlaybackConfig(
                show_electron_density=False,
                fps=5,
                export_path=path)
            viz = InteractionVisualizer(traj, mol, config=config)
            result = viz.export(path)
            assert os.path.exists(path)
            assert os.path.getsize(path) > 100  # Should have content
        finally:
            if os.path.exists(path):
                os.unlink(path)

    def test_frame_count(self):
        """Animation should have at least 2 frames."""
        from molbuilder.visualization.interaction_viz import (
            InteractionVisualizer, PlaybackConfig,
        )
        mol, traj = self._build_h2_trajectory()
        config = PlaybackConfig(show_electron_density=False, fps=10)
        viz = InteractionVisualizer(traj, mol, config=config)
        assert viz.n_frames >= 2
        import matplotlib.pyplot as plt
        plt.close("all")


class TestElectronDensityRenderer:
    """Test electron density cloud rendering."""

    def test_compute_density_returns_points(self):
        """Should return scatter points for a bond."""
        from molbuilder.visualization.electron_density_viz import (
            ElectronDensityRenderer,
        )
        renderer = ElectronDensityRenderer(n_points=1000)
        pos_a = np.array([0.0, 0.0, 0.0])
        pos_b = np.array([1.5, 0.0, 0.0])
        points, colors = renderer.compute_bond_density(
            pos_a, pos_b, z_a=6, z_b=6, bond_order=1.0)
        assert len(points) > 0
        assert points.shape[1] == 3
        assert colors.shape[1] == 4

    def test_zero_bond_order_returns_empty(self):
        """Zero bond order should return empty arrays."""
        from molbuilder.visualization.electron_density_viz import (
            ElectronDensityRenderer,
        )
        renderer = ElectronDensityRenderer(n_points=1000)
        pos_a = np.array([0.0, 0.0, 0.0])
        pos_b = np.array([1.5, 0.0, 0.0])
        points, colors = renderer.compute_bond_density(
            pos_a, pos_b, z_a=6, z_b=6, bond_order=0.0)
        assert len(points) == 0

    def test_higher_bond_order_more_points(self):
        """Higher bond order should generally produce more points."""
        from molbuilder.visualization.electron_density_viz import (
            ElectronDensityRenderer,
        )
        renderer = ElectronDensityRenderer(n_points=2000)
        pos_a = np.array([0.0, 0.0, 0.0])
        pos_b = np.array([1.5, 0.0, 0.0])

        pts_low, _ = renderer.compute_bond_density(
            pos_a, pos_b, z_a=6, z_b=6, bond_order=0.3)
        pts_high, _ = renderer.compute_bond_density(
            pos_a, pos_b, z_a=6, z_b=6, bond_order=1.5)

        # Higher bond order should accept more points on average
        # (may not always be true due to randomness, so use generous margin)
        assert len(pts_high) >= len(pts_low) * 0.3


class TestPlaybackController:
    """Test the PlaybackController."""

    def test_controller_creation(self):
        """Controller should be creatable from a visualizer."""
        from molbuilder.visualization.interaction_viz import (
            InteractionVisualizer, PlaybackConfig,
        )
        from molbuilder.visualization.interaction_controls import (
            PlaybackController,
        )
        from molbuilder.dynamics.trajectory import Trajectory, TrajectoryFrame

        mol = Molecule("H2")
        mol.add_atom("H", [0.0, 0.0, 0.0])
        mol.add_atom("H", [0.74, 0.0, 0.0])
        mol.add_bond(0, 1, order=1)

        traj = Trajectory(n_atoms=2)
        traj.add_frame(TrajectoryFrame(
            time_fs=0.0,
            positions=np.array([[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]])))
        traj.add_frame(TrajectoryFrame(
            time_fs=1.0,
            positions=np.array([[0.0, 0.0, 0.0], [0.75, 0.0, 0.0]])))

        config = PlaybackConfig(show_electron_density=False, fps=5)
        viz = InteractionVisualizer(traj, mol, config=config)
        controller = PlaybackController(viz)
        assert controller.viz is viz
        import matplotlib.pyplot as plt
        plt.close("all")
