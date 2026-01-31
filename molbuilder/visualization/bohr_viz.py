"""Bohr model visualization: animated orbital diagram.

Migrated from legacy/bohr_model.py -- wavelength_to_rgb and visualize functions.
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
from molbuilder.visualization.theme import BG_COLOR, TEXT_COLOR, GRID_COLOR


def wavelength_to_rgb(wavelength_nm: float) -> tuple:
    """Convert a visible-light wavelength (380-780 nm) to an RGB tuple.
    Returns white for wavelengths outside visible range."""
    if wavelength_nm < 380 or wavelength_nm > 780:
        return (1.0, 1.0, 1.0)

    if wavelength_nm < 440:
        r = -(wavelength_nm - 440) / (440 - 380)
        g = 0.0
        b = 1.0
    elif wavelength_nm < 490:
        r = 0.0
        g = (wavelength_nm - 440) / (490 - 440)
        b = 1.0
    elif wavelength_nm < 510:
        r = 0.0
        g = 1.0
        b = -(wavelength_nm - 510) / (510 - 490)
    elif wavelength_nm < 580:
        r = (wavelength_nm - 510) / (580 - 510)
        g = 1.0
        b = 0.0
    elif wavelength_nm < 645:
        r = 1.0
        g = -(wavelength_nm - 645) / (645 - 580)
        b = 0.0
    else:
        r = 1.0
        g = 0.0
        b = 0.0

    # intensity fall-off at edges of visible spectrum
    if wavelength_nm < 420:
        factor = 0.3 + 0.7 * (wavelength_nm - 380) / (420 - 380)
    elif wavelength_nm > 700:
        factor = 0.3 + 0.7 * (780 - wavelength_nm) / (780 - 700)
    else:
        factor = 1.0

    return (r * factor, g * factor, b * factor)


def visualize(atom, animate: bool = True, interval_ms: int = 30):
    """Render an animated Bohr model diagram of the atom.

    Parameters
    ----------
    atom : BohrAtom
        The atom to visualize.
    animate : bool
        If True, electrons orbit the nucleus. If False, show a static frame.
    interval_ms : int
        Milliseconds between animation frames.
    """
    n_shells = atom.num_shells
    if n_shells == 0:
        print("No electrons to visualize.")
        return

    fig, ax = plt.subplots(1, 1, figsize=(8, 8), facecolor="black")
    ax.set_facecolor("black")
    ax.set_aspect("equal")
    ax.set_xlim(-n_shells - 1, n_shells + 1)
    ax.set_ylim(-n_shells - 1, n_shells + 1)
    ax.axis("off")

    # Title
    charge_label = ""
    if atom.charge > 0:
        charge_label = f"$^{{+{atom.charge}}}$"
    elif atom.charge < 0:
        charge_label = f"$^{{{atom.charge}}}$"
    ax.set_title(
        f"Bohr Model - {atom.name} ({atom.symbol}{charge_label})  "
        f"Z={atom.atomic_number}  e={atom.num_electrons}",
        color="white", fontsize=14, pad=12,
    )

    # ----- Nucleus -----
    nucleus_radius = 0.25 + 0.03 * atom.atomic_number**0.33
    nucleus = Circle((0, 0), nucleus_radius, color="#ff6633", zorder=10)
    ax.add_patch(nucleus)
    ax.text(0, 0, f"{atom.protons}p\n{atom.neutrons}n",
            ha="center", va="center", fontsize=7, color="white",
            fontweight="bold", zorder=11)

    # ----- Orbital rings -----
    shell_radii = []
    for n in range(1, n_shells + 1):
        r = n  # use integer spacing for visual clarity
        shell_radii.append(r)
        orbit = Circle((0, 0), r, fill=False, edgecolor="#334466",
                        linewidth=0.8, linestyle="--", zorder=1)
        ax.add_patch(orbit)
        ax.text(r + 0.15, 0.15, f"n={n}", fontsize=7, color="#5588aa", zorder=2)

    # ----- Electron dots (initial positions) -----
    electron_artists = []
    electron_positions = []  # (shell_index, angle_offset, shell_radius)

    for shell_idx, count in enumerate(atom.shell_config):
        r = shell_radii[shell_idx]
        for e in range(count):
            angle = 2 * math.pi * e / count
            x = r * math.cos(angle)
            y = r * math.sin(angle)
            dot = ax.plot(x, y, 'o', color="#44ccff", markersize=5, zorder=5)[0]
            electron_artists.append(dot)
            electron_positions.append((shell_idx, angle, r, count))

    # ----- Shell electron count labels -----
    for shell_idx, count in enumerate(atom.shell_config):
        r = shell_radii[shell_idx]
        ax.text(-r - 0.15, -0.25, str(count), fontsize=8, color="#88bbdd",
                ha="right", zorder=2)

    # ----- Energy level sidebar -----
    sidebar_x = n_shells + 0.6
    e_min = atom.energy_level(1)
    e_max = atom.energy_level(n_shells) if n_shells > 1 else e_min * 0.1
    e_range = abs(e_max - e_min) if abs(e_max - e_min) > 0 else 1.0
    bar_bottom = -n_shells
    bar_height = 2 * n_shells

    for n in range(1, n_shells + 1):
        e = atom.energy_level(n)
        if n_shells > 1:
            y_pos = bar_bottom + bar_height * (e - e_min) / e_range
        else:
            y_pos = 0
        ax.plot([sidebar_x, sidebar_x + 0.5], [y_pos, y_pos],
                color="#ffaa33", linewidth=1.5, zorder=3)
        ax.text(sidebar_x + 0.6, y_pos, f"{e:.2f} eV",
                fontsize=6, color="#ffcc66", va="center", zorder=3)

    # ----- Animation -----
    def update(frame):
        for i, (shell_idx, angle0, r, count) in enumerate(electron_positions):
            n = shell_idx + 1
            angular_speed = 0.06 / n  # outer shells orbit slower
            angle = angle0 + angular_speed * frame
            x = r * math.cos(angle)
            y = r * math.sin(angle)
            electron_artists[i].set_data([x], [y])
        return electron_artists

    if animate:
        anim = animation.FuncAnimation(
            fig, update, frames=None, interval=interval_ms, blit=True,
        )

    plt.tight_layout()
    plt.show()
