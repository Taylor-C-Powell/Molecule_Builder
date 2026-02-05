"""Playback controls for the interaction visualizer.

Provides keyboard bindings for interactive animation control and an
optional tkinter panel for GUI embedding.

Keyboard bindings:
    Space     : Play / Pause
    Right     : Step forward
    Left      : Step backward
    Up        : Increase speed (2x)
    Down      : Decrease speed (0.5x)
    E         : Toggle electron density
    L         : Toggle labels
    R         : Reset to beginning
    Q / Esc   : Close
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from molbuilder.visualization.interaction_viz import InteractionVisualizer


class PlaybackController:
    """Keyboard-driven playback controller for InteractionVisualizer.

    Connects to matplotlib key_press_event to provide interactive
    playback controls.

    Parameters
    ----------
    visualizer : InteractionVisualizer
        The visualizer to control.
    """

    def __init__(self, visualizer: InteractionVisualizer):
        self.viz = visualizer
        self._speed_multiplier = 1.0
        self._current_frame = 0
        self._connected = False

    def connect(self):
        """Connect keyboard event handlers to the visualizer's figure."""
        if self.viz.fig is None:
            return
        self.viz.fig.canvas.mpl_connect("key_press_event", self._on_key)
        self._connected = True

    def _on_key(self, event):
        """Handle a key press event."""
        if event.key == " ":
            self._toggle_pause()
        elif event.key == "right":
            self._step_forward()
        elif event.key == "left":
            self._step_backward()
        elif event.key == "up":
            self._speed_up()
        elif event.key == "down":
            self._slow_down()
        elif event.key == "e":
            self._toggle_electron_density()
        elif event.key == "l":
            self._toggle_labels()
        elif event.key == "r":
            self._reset()
        elif event.key in ("q", "escape"):
            self._close()

    def _toggle_pause(self):
        """Toggle play/pause."""
        anim = self.viz._anim
        if anim is None:
            return
        if self.viz._paused:
            anim.resume()
            self.viz._paused = False
        else:
            anim.pause()
            self.viz._paused = True

    def _step_forward(self):
        """Step one frame forward (pauses first)."""
        anim = self.viz._anim
        if anim is None:
            return
        if not self.viz._paused:
            anim.pause()
            self.viz._paused = True
        self._current_frame = min(
            self._current_frame + 1, self.viz.n_frames - 1)
        self.viz._render_frame(self._current_frame)
        self.viz.fig.canvas.draw_idle()

    def _step_backward(self):
        """Step one frame backward."""
        anim = self.viz._anim
        if anim is None:
            return
        if not self.viz._paused:
            anim.pause()
            self.viz._paused = True
        self._current_frame = max(self._current_frame - 1, 0)
        self.viz._render_frame(self._current_frame)
        self.viz.fig.canvas.draw_idle()

    def _speed_up(self):
        """Double the playback speed."""
        self._speed_multiplier *= 2.0
        anim = self.viz._anim
        if anim is not None:
            new_interval = max(
                1, int(1000 / (self.viz.config.fps * self._speed_multiplier)))
            anim.event_source.interval = new_interval

    def _slow_down(self):
        """Halve the playback speed."""
        self._speed_multiplier *= 0.5
        anim = self.viz._anim
        if anim is not None:
            new_interval = int(
                1000 / (self.viz.config.fps * self._speed_multiplier))
            anim.event_source.interval = new_interval

    def _toggle_electron_density(self):
        """Toggle electron density cloud rendering."""
        self.viz.config.show_electron_density = (
            not self.viz.config.show_electron_density)

    def _toggle_labels(self):
        """Toggle time and energy labels."""
        self.viz.config.show_time_label = not self.viz.config.show_time_label
        self.viz.config.show_energy_bar = not self.viz.config.show_energy_bar

    def _reset(self):
        """Reset to the first frame."""
        self._current_frame = 0
        if self.viz._paused:
            self.viz._render_frame(0)
            self.viz.fig.canvas.draw_idle()

    def _close(self):
        """Close the figure."""
        import matplotlib.pyplot as plt
        plt.close(self.viz.fig)


def create_tkinter_panel(parent, visualizer: InteractionVisualizer):
    """Create an optional tkinter control panel for GUI embedding.

    Parameters
    ----------
    parent : tk.Widget
        Parent tkinter widget.
    visualizer : InteractionVisualizer
        The visualizer to control.

    Returns
    -------
    tk.Frame
        The control panel frame.
    """
    import tkinter as tk
    from tkinter import ttk

    controller = PlaybackController(visualizer)

    frame = ttk.Frame(parent)

    play_btn = ttk.Button(
        frame, text="Play/Pause",
        command=controller._toggle_pause)
    play_btn.pack(side="left", padx=2)

    step_back_btn = ttk.Button(
        frame, text="<<",
        command=controller._step_backward)
    step_back_btn.pack(side="left", padx=2)

    step_fwd_btn = ttk.Button(
        frame, text=">>",
        command=controller._step_forward)
    step_fwd_btn.pack(side="left", padx=2)

    slower_btn = ttk.Button(
        frame, text="Slower",
        command=controller._slow_down)
    slower_btn.pack(side="left", padx=2)

    faster_btn = ttk.Button(
        frame, text="Faster",
        command=controller._speed_up)
    faster_btn.pack(side="left", padx=2)

    edensity_var = tk.BooleanVar(value=visualizer.config.show_electron_density)
    edensity_cb = ttk.Checkbutton(
        frame, text="e- Density",
        variable=edensity_var,
        command=controller._toggle_electron_density)
    edensity_cb.pack(side="left", padx=2)

    labels_var = tk.BooleanVar(value=visualizer.config.show_time_label)
    labels_cb = ttk.Checkbutton(
        frame, text="Labels",
        variable=labels_var,
        command=controller._toggle_labels)
    labels_cb.pack(side="left", padx=2)

    return frame
