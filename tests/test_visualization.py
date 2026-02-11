"""Tests for the visualization subpackage.

Covers pure-computation helpers, PlaybackController state methods,
and smoke tests for all plotting functions using headless Agg backend.
"""

import math
import warnings

import numpy as np
import pytest

import matplotlib
matplotlib.use("Agg")  # Headless backend -- must be set before pyplot import
import matplotlib.pyplot as plt

from molbuilder.molecule.graph import Molecule


# ===================================================================
# 1. wavelength_to_rgb  (bohr_viz.py)
# ===================================================================

class TestWavelengthToRGB:
    """Pure-computation tests for bohr_viz.wavelength_to_rgb."""

    def setup_method(self):
        from molbuilder.visualization.bohr_viz import wavelength_to_rgb
        self.fn = wavelength_to_rgb

    def test_below_visible_returns_white(self):
        assert self.fn(300) == (1.0, 1.0, 1.0)

    def test_above_visible_returns_white(self):
        assert self.fn(800) == (1.0, 1.0, 1.0)

    def test_blue_region(self):
        r, g, b = self.fn(440)
        assert b == 1.0
        assert r == pytest.approx(0.0, abs=0.01)

    def test_green_region(self):
        r, g, b = self.fn(510)
        assert g == 1.0

    def test_red_region(self):
        r, g, b = self.fn(650)
        assert r == 1.0
        assert b == 0.0

    def test_edge_falloff_low(self):
        """At 390 nm, factor should be < 1.0, dimming the result."""
        r, g, b = self.fn(390)
        # At 390 nm: r = -(390-440)/(440-380) = 50/60 ≈ 0.833, b = 1.0
        # factor = 0.3 + 0.7*(390-380)/(420-380) = 0.3 + 0.175 = 0.475
        assert b < 1.0  # factor applied

    def test_edge_falloff_high(self):
        """At 750 nm, factor should be < 1.0, dimming the result."""
        r, g, b = self.fn(750)
        assert r < 1.0  # factor applied


# ===================================================================
# 2. _perpendicular_offset  (molecule_viz.py)
# ===================================================================

class TestPerpendicularOffset:
    """Pure-computation tests for molecule_viz._perpendicular_offset."""

    def setup_method(self):
        from molbuilder.visualization.molecule_viz import _perpendicular_offset
        self.fn = _perpendicular_offset

    def test_perpendicular_to_bond(self):
        pos_a = np.array([0.0, 0.0, 0.0])
        pos_b = np.array([0.0, 0.0, 2.0])
        offset = self.fn(pos_a, pos_b, 0.1)
        bond_vec = pos_b - pos_a
        dot = np.dot(offset, bond_vec)
        assert abs(dot) < 1e-10

    def test_correct_magnitude(self):
        pos_a = np.array([0.0, 0.0, 0.0])
        pos_b = np.array([1.0, 1.0, 0.0])
        mag = 0.5
        offset = self.fn(pos_a, pos_b, mag)
        assert np.linalg.norm(offset) == pytest.approx(mag, abs=1e-10)

    def test_zero_length_bond(self):
        pos_a = np.array([1.0, 2.0, 3.0])
        pos_b = np.array([1.0, 2.0, 3.0])
        offset = self.fn(pos_a, pos_b, 0.1)
        np.testing.assert_allclose(offset, [0.1, 0.0, 0.0])

    def test_bond_along_x_uses_y_reference(self):
        """When bond is parallel to x, the function should pick y as reference."""
        pos_a = np.array([0.0, 0.0, 0.0])
        pos_b = np.array([5.0, 0.0, 0.0])
        offset = self.fn(pos_a, pos_b, 0.2)
        # Should be perpendicular to x-axis
        assert abs(offset[0]) < 1e-10
        assert np.linalg.norm(offset) == pytest.approx(0.2, abs=1e-10)


# ===================================================================
# 3. _slerp  (molecule_viz.py)
# ===================================================================

class TestSlerp:
    """Pure-computation tests for molecule_viz._slerp."""

    def setup_method(self):
        from molbuilder.visualization.molecule_viz import _slerp
        self.fn = _slerp

    def test_t0_returns_va(self):
        va = np.array([1.0, 0.0, 0.0])
        vb = np.array([0.0, 1.0, 0.0])
        result = self.fn(va, vb, 0.0)
        np.testing.assert_allclose(result, va, atol=1e-10)

    def test_t1_returns_vb(self):
        va = np.array([1.0, 0.0, 0.0])
        vb = np.array([0.0, 1.0, 0.0])
        result = self.fn(va, vb, 1.0)
        np.testing.assert_allclose(result, vb, atol=1e-10)

    def test_midpoint_is_unit_vector(self):
        va = np.array([1.0, 0.0, 0.0])
        vb = np.array([0.0, 1.0, 0.0])
        result = self.fn(va, vb, 0.5)
        # Midpoint on unit sphere should have unit length
        assert np.linalg.norm(result) == pytest.approx(1.0, abs=1e-10)
        # Should be at 45 degrees from both
        expected = np.array([1.0, 1.0, 0.0]) / math.sqrt(2)
        np.testing.assert_allclose(result, expected, atol=1e-10)

    def test_nearly_parallel_fallback(self):
        va = np.array([1.0, 0.0, 0.0])
        vb = np.array([1.0, 1e-8, 0.0])
        vb = vb / np.linalg.norm(vb)
        # Should not crash -- falls back to linear interpolation
        result = self.fn(va, vb, 0.5)
        assert np.linalg.norm(result) > 0


# ===================================================================
# 4. Electron density helpers  (electron_density_viz.py)
# ===================================================================

class TestElectronDensityHelpers:
    """Pure-computation tests for electron_density_viz helpers."""

    def test_hex_to_rgba_red(self):
        from molbuilder.visualization.electron_density_viz import _hex_to_rgba
        result = _hex_to_rgba("#ff0000", 0.5)
        assert result == pytest.approx((1.0, 0.0, 0.0, 0.5), abs=1e-3)

    def test_hex_to_rgba_green(self):
        from molbuilder.visualization.electron_density_viz import _hex_to_rgba
        result = _hex_to_rgba("#00ff00", 0.3)
        assert result == pytest.approx((0.0, 1.0, 0.0, 0.3), abs=1e-3)

    def test_hex_to_rgba_no_hash(self):
        from molbuilder.visualization.electron_density_viz import _hex_to_rgba
        result = _hex_to_rgba("0000ff", 1.0)
        assert result == pytest.approx((0.0, 0.0, 1.0, 1.0), abs=1e-3)

    def test_slater_zeff_hydrogen(self):
        from molbuilder.visualization.electron_density_viz import _slater_zeff_simple
        assert _slater_zeff_simple(1) == 1.0

    def test_slater_zeff_carbon(self):
        from molbuilder.visualization.electron_density_viz import _slater_zeff_simple
        assert _slater_zeff_simple(6) == 3.25

    def test_slater_zeff_fallback(self):
        from molbuilder.visualization.electron_density_viz import _slater_zeff_simple
        result = _slater_zeff_simple(999)
        assert result == max(1.0, 999 * 0.3)

    def test_valence_nl_hydrogen(self):
        from molbuilder.visualization.electron_density_viz import _valence_nl
        assert _valence_nl(1) == (1, 0)

    def test_valence_nl_carbon(self):
        from molbuilder.visualization.electron_density_viz import _valence_nl
        assert _valence_nl(6) == (2, 1)

    def test_valence_nl_calcium(self):
        from molbuilder.visualization.electron_density_viz import _valence_nl
        assert _valence_nl(20) == (4, 0)

    def test_valence_nl_helium(self):
        from molbuilder.visualization.electron_density_viz import _valence_nl
        assert _valence_nl(2) == (1, 0)

    def test_valence_nl_sodium(self):
        from molbuilder.visualization.electron_density_viz import _valence_nl
        assert _valence_nl(11) == (3, 0)


# ===================================================================
# 5. _bond_color_from_order  (interaction_viz.py)
# ===================================================================

class TestBondColorFromOrder:
    """Pure-computation tests for interaction_viz._bond_color_from_order."""

    def setup_method(self):
        from molbuilder.visualization.interaction_viz import _bond_color_from_order
        from molbuilder.visualization.theme import (
            BOND_COLOR, FORMING_BOND_COLOR, BREAKING_BOND_COLOR,
            TRANSITION_STATE_COLOR,
        )
        self.fn = _bond_color_from_order
        self.BOND = BOND_COLOR
        self.FORMING = FORMING_BOND_COLOR
        self.BREAKING = BREAKING_BOND_COLOR
        self.TRANSITION = TRANSITION_STATE_COLOR

    def test_very_low_order_is_breaking(self):
        assert self.fn(0.1) == self.BREAKING

    def test_transition_state(self):
        assert self.fn(0.5) == self.TRANSITION

    def test_normal_bond(self):
        assert self.fn(1.0) == self.BOND

    def test_high_order_is_forming(self):
        assert self.fn(1.5) == self.FORMING


# ===================================================================
# 6. PlaybackController  (interaction_controls.py)
# ===================================================================

class TestPlaybackControllerMethods:
    """State-method tests for PlaybackController with a real visualizer."""

    def _make_viz(self):
        """Build a minimal InteractionVisualizer with H2 trajectory."""
        from molbuilder.visualization.interaction_viz import (
            InteractionVisualizer, PlaybackConfig,
        )
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

        config = PlaybackConfig(show_electron_density=False, fps=10)
        viz = InteractionVisualizer(traj, mol, config=config)
        viz.build_animation()
        return viz

    def setup_method(self):
        from molbuilder.visualization.interaction_controls import PlaybackController
        self.viz = self._make_viz()
        self.ctrl = PlaybackController(self.viz)

    def teardown_method(self):
        plt.close("all")

    def test_toggle_pause(self):
        assert not self.viz._paused
        self.ctrl._toggle_pause()
        assert self.viz._paused

    def test_toggle_pause_twice(self):
        self.ctrl._toggle_pause()
        self.ctrl._toggle_pause()
        assert not self.viz._paused

    def test_step_forward(self):
        self.ctrl._step_forward()
        assert self.ctrl._current_frame == 1

    def test_step_backward_clamped(self):
        self.ctrl._step_backward()
        assert self.ctrl._current_frame == 0

    def test_step_forward_then_backward(self):
        self.ctrl._step_forward()
        self.ctrl._step_forward()
        assert self.ctrl._current_frame == 2
        self.ctrl._step_backward()
        assert self.ctrl._current_frame == 1

    def test_speed_up(self):
        self.ctrl._speed_up()
        assert self.ctrl._speed_multiplier == pytest.approx(2.0)

    def test_slow_down(self):
        self.ctrl._slow_down()
        assert self.ctrl._speed_multiplier == pytest.approx(0.5)

    def test_toggle_electron_density(self):
        original = self.viz.config.show_electron_density
        self.ctrl._toggle_electron_density()
        assert self.viz.config.show_electron_density is (not original)

    def test_toggle_labels(self):
        original_time = self.viz.config.show_time_label
        original_energy = self.viz.config.show_energy_bar
        self.ctrl._toggle_labels()
        assert self.viz.config.show_time_label is (not original_time)
        assert self.viz.config.show_energy_bar is (not original_energy)

    def test_reset(self):
        self.ctrl._step_forward()
        self.ctrl._step_forward()
        self.ctrl._reset()
        assert self.ctrl._current_frame == 0

    def test_connect(self):
        self.ctrl.connect()
        assert self.ctrl._connected is True


# ===================================================================
# 7. quantum_viz smoke tests
# ===================================================================

class TestQuantumVizSmoke:
    """Smoke tests for quantum_viz plotting functions (Agg backend)."""

    def teardown_method(self):
        plt.close("all")

    def test_plot_orbital_3d(self):
        from molbuilder.visualization.quantum_viz import plot_orbital_3d
        plot_orbital_3d(1, 0, 0, num_points=500)
        assert len(plt.get_fignums()) >= 1

    def test_plot_radial_wavefunction(self):
        from molbuilder.visualization.quantum_viz import plot_radial_wavefunction
        plot_radial_wavefunction([(1, 0), (2, 0)])
        assert len(plt.get_fignums()) >= 1

    def test_plot_radial_probability(self):
        from molbuilder.visualization.quantum_viz import plot_radial_probability
        plot_radial_probability([(1, 0), (2, 1)])
        assert len(plt.get_fignums()) >= 1

    def test_plot_angular_distribution(self):
        from molbuilder.visualization.quantum_viz import plot_angular_distribution
        plot_angular_distribution(1, 0)
        assert len(plt.get_fignums()) >= 1

    def test_plot_electron_configuration(self):
        from molbuilder.visualization.quantum_viz import plot_electron_configuration
        from molbuilder.atomic.quantum_atom import QuantumAtom
        atom = QuantumAtom(6)  # Carbon
        plot_electron_configuration(atom)
        assert len(plt.get_fignums()) >= 1

    def test_plot_energy_levels(self):
        from molbuilder.visualization.quantum_viz import plot_energy_levels
        from molbuilder.atomic.quantum_atom import QuantumAtom
        atom = QuantumAtom(6)  # Carbon
        plot_energy_levels(atom)
        assert len(plt.get_fignums()) >= 1


# ===================================================================
# 8. bohr_viz and molecule_viz smoke tests
# ===================================================================

class TestBohrAndMoleculeVizSmoke:
    """Smoke tests for bohr_viz.visualize and molecule_viz functions."""

    def teardown_method(self):
        plt.close("all")

    def test_bohr_visualize_static(self):
        from molbuilder.visualization.bohr_viz import visualize
        from molbuilder.atomic.bohr import BohrAtom
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            atom = BohrAtom(1)  # Hydrogen
        visualize(atom, animate=False)
        assert len(plt.get_fignums()) >= 1

    def test_draw_angle_arc(self):
        from molbuilder.visualization.molecule_viz import _draw_angle_arc
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        center = np.array([0.0, 0.0, 0.0])
        vec_a = np.array([1.0, 0.0, 0.0])
        vec_b = np.array([0.0, 1.0, 0.0])
        _draw_angle_arc(ax, center, vec_a, vec_b, 90.0)
        # No crash means success

    def test_visualize_molecule(self):
        from molbuilder.visualization.molecule_viz import visualize_molecule
        from molbuilder.bonding.vsepr import VSEPRMolecule
        mol = VSEPRMolecule("H2O")
        visualize_molecule(mol)
        assert len(plt.get_fignums()) >= 1

    def test_visualize_gallery(self):
        from molbuilder.visualization.molecule_viz import visualize_gallery
        from molbuilder.bonding.vsepr import VSEPRMolecule
        mols = [VSEPRMolecule("H2O"), VSEPRMolecule("CH4")]
        visualize_gallery(mols)
        assert len(plt.get_fignums()) >= 1
