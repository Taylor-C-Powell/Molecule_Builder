"""Quantum Mechanical Atom Visualization

Provides plotting functions for:
    - 3D orbital probability clouds (electron density scatter)
    - Radial wave functions  R_nl(r)
    - Radial probability distributions  r^2 |R_nl|^2
    - Angular probability distributions  |Y_l^m|^2
    - Orbital energy level diagrams
    - Electron configuration box diagrams (with arrows)

Migrated from legacy/quantum_visualization.py.
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D

from molbuilder.core.constants import BOHR_RADIUS_M as BOHR_RADIUS
from molbuilder.atomic.quantum_numbers import SUBSHELL_LETTER, ORBITAL_NAMES
from molbuilder.atomic.wavefunctions import (
    radial_wavefunction,
    real_spherical_harmonic,
    wavefunction_real,
    radial_probability_density,
    expectation_r,
    orbital_label,
)
from molbuilder.core.geometry import cartesian_to_spherical, spherical_to_cartesian
from molbuilder.visualization.theme import (
    BG_COLOR, TEXT_COLOR, GRID_COLOR,
    POSITIVE_COLOR, NEGATIVE_COLOR, NEUTRAL_COLOR, ENERGY_COLOR,
)


# ===================================================================
# 3D Orbital probability cloud
# ===================================================================

def plot_orbital_3d(n: int, l: int, m: int, Z: int = 1,
                    num_points: int = 30000, threshold: float = 0.3,
                    figsize: tuple = (8, 8)):
    """Render a 3D probability cloud for orbital (n, l, m).

    Uses Monte Carlo rejection sampling to generate points distributed
    according to |psi|^2, then plots as a 3D scatter with positive/negative
    lobes coloured differently.

    Parameters
    ----------
    n, l, m   : quantum numbers
    Z         : nuclear charge
    num_points: number of candidate random points
    threshold : fraction of max density below which points are culled
    figsize   : figure size
    """
    # Determine radial extent: ~5x the expectation value of r
    r_extent = 5.0 * expectation_r(n, l, Z)
    r_extent_a0 = r_extent / BOHR_RADIUS  # in units of a_0 for display

    # Generate random points in a cube, then convert to spherical
    side = r_extent
    rng = np.random.default_rng(42)
    x = rng.uniform(-side, side, num_points)
    y = rng.uniform(-side, side, num_points)
    z_coord = rng.uniform(-side, side, num_points)

    r, theta, phi = cartesian_to_spherical(x, y, z_coord)

    # Compute real wave function and density
    psi = wavefunction_real(n, l, m, r, theta, phi, Z)
    density = psi**2

    # Rejection sampling: keep points with probability proportional to density
    max_density = np.max(density)
    if max_density == 0:
        print("Wave function is zero everywhere in sampled region.")
        return

    accept_prob = density / max_density
    randoms = rng.uniform(0, 1, num_points)
    mask = randoms < accept_prob

    # Also cull very low-density points for visual clarity
    mask &= density > threshold * max_density

    x_keep = x[mask] / BOHR_RADIUS  # convert to a_0 units for display
    y_keep = y[mask] / BOHR_RADIUS
    z_keep = z_coord[mask] / BOHR_RADIUS
    psi_keep = psi[mask]

    # Color by sign of wave function
    colors = np.where(psi_keep >= 0, POSITIVE_COLOR, NEGATIVE_COLOR)

    # Plot
    fig = plt.figure(figsize=figsize, facecolor=BG_COLOR)
    ax = fig.add_subplot(111, projection="3d", facecolor=BG_COLOR)

    ax.scatter(x_keep, y_keep, z_keep, c=colors, s=1.0, alpha=0.5,
               depthshade=True)

    lim = r_extent_a0 * 0.8
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)

    label = orbital_label(n, l, m)
    ax.set_title(f"Orbital {label}  (Z={Z})", color=TEXT_COLOR, fontsize=14)
    ax.set_xlabel("x / a0", color=TEXT_COLOR, fontsize=9)
    ax.set_ylabel("y / a0", color=TEXT_COLOR, fontsize=9)
    ax.set_zlabel("z / a0", color=TEXT_COLOR, fontsize=9)
    ax.tick_params(colors=TEXT_COLOR, labelsize=7)

    # Style pane colors
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.set_facecolor(BG_COLOR)
        pane.set_edgecolor(GRID_COLOR)

    plt.tight_layout()
    plt.show()


# ===================================================================
# Radial wave function plot
# ===================================================================

def plot_radial_wavefunction(n_l_pairs: list[tuple[int, int]],
                             Z: int = 1, r_max_a0: float = None,
                             figsize: tuple = (9, 5)):
    """Plot R_nl(r) for one or more (n, l) pairs.

    Parameters
    ----------
    n_l_pairs : list of (n, l) tuples
    Z         : nuclear charge
    r_max_a0  : maximum r in units of a_0 (auto-scaled if None)
    """
    if r_max_a0 is None:
        r_max_a0 = max(
            5 * expectation_r(n, l, Z) / BOHR_RADIUS for n, l in n_l_pairs
        )

    r_a0 = np.linspace(1e-6, r_max_a0, 2000)
    r_m = r_a0 * BOHR_RADIUS

    fig, ax = plt.subplots(figsize=figsize, facecolor=BG_COLOR)
    ax.set_facecolor(BG_COLOR)

    for n, l in n_l_pairs:
        R = radial_wavefunction(n, l, r_m, Z)
        # Scale to a_0 units for display: R has units of m^{-3/2}
        R_scaled = R * BOHR_RADIUS**1.5
        label_str = f"R({n},{SUBSHELL_LETTER.get(l,'?')})"
        ax.plot(r_a0, R_scaled, linewidth=1.5, label=label_str)

    ax.axhline(0, color=GRID_COLOR, linewidth=0.5)
    ax.set_xlabel("r / a0", color=TEXT_COLOR, fontsize=11)
    ax.set_ylabel("R(r) * a0^(3/2)", color=TEXT_COLOR, fontsize=11)
    ax.set_title(f"Radial Wave Functions (Z={Z})", color=TEXT_COLOR, fontsize=13)
    ax.legend(facecolor="#111122", edgecolor=GRID_COLOR, labelcolor=TEXT_COLOR)
    ax.tick_params(colors=TEXT_COLOR)
    for spine in ax.spines.values():
        spine.set_color(GRID_COLOR)
    ax.grid(True, color=GRID_COLOR, alpha=0.3)

    plt.tight_layout()
    plt.show()


# ===================================================================
# Radial probability distribution plot
# ===================================================================

def plot_radial_probability(n_l_pairs: list[tuple[int, int]],
                            Z: int = 1, r_max_a0: float = None,
                            figsize: tuple = (9, 5)):
    """Plot r^2 |R_nl(r)|^2 for one or more (n, l) pairs."""
    if r_max_a0 is None:
        r_max_a0 = max(
            5 * expectation_r(n, l, Z) / BOHR_RADIUS for n, l in n_l_pairs
        )

    r_a0 = np.linspace(1e-6, r_max_a0, 2000)
    r_m = r_a0 * BOHR_RADIUS

    fig, ax = plt.subplots(figsize=figsize, facecolor=BG_COLOR)
    ax.set_facecolor(BG_COLOR)

    for n, l in n_l_pairs:
        P = radial_probability_density(n, l, r_m, Z)
        # Scale: P has units of m^{-1}, multiply by a_0 to get dimensionless
        P_scaled = P * BOHR_RADIUS
        label_str = f"P({n},{SUBSHELL_LETTER.get(l,'?')})"
        ax.plot(r_a0, P_scaled, linewidth=1.5, label=label_str)
        ax.fill_between(r_a0, P_scaled, alpha=0.15)

    ax.set_xlabel("r / a0", color=TEXT_COLOR, fontsize=11)
    ax.set_ylabel("P(r) * a0", color=TEXT_COLOR, fontsize=11)
    ax.set_title(
        f"Radial Probability Distribution (Z={Z})",
        color=TEXT_COLOR, fontsize=13,
    )
    ax.legend(facecolor="#111122", edgecolor=GRID_COLOR, labelcolor=TEXT_COLOR)
    ax.tick_params(colors=TEXT_COLOR)
    for spine in ax.spines.values():
        spine.set_color(GRID_COLOR)
    ax.grid(True, color=GRID_COLOR, alpha=0.3)

    plt.tight_layout()
    plt.show()


# ===================================================================
# Angular distribution cross-section
# ===================================================================

def plot_angular_distribution(l: int, m: int, figsize: tuple = (6, 6)):
    """Polar plot of |Y_l^m(theta, phi=0)|^2 in the xz-plane."""
    theta = np.linspace(0, 2 * np.pi, 500)
    phi_fixed = np.zeros_like(theta)

    Y = real_spherical_harmonic(l, m, theta, phi_fixed)
    mag = np.abs(Y)

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"},
                           figsize=figsize, facecolor=BG_COLOR)
    ax.set_facecolor(BG_COLOR)

    # Color by sign
    pos_mask = Y >= 0
    neg_mask = ~pos_mask

    ax.plot(theta[pos_mask], mag[pos_mask], ".", color=POSITIVE_COLOR,
            markersize=1.5)
    ax.plot(theta[neg_mask], mag[neg_mask], ".", color=NEGATIVE_COLOR,
            markersize=1.5)
    ax.fill_between(theta, mag, alpha=0.15, color=NEUTRAL_COLOR)

    name = ORBITAL_NAMES.get((l, m), f"l={l},m={m}")
    ax.set_title(f"|Y({name})|  (xz plane)", color=TEXT_COLOR,
                 fontsize=12, pad=15)
    ax.tick_params(colors=TEXT_COLOR, labelsize=7)
    ax.grid(True, color=GRID_COLOR, alpha=0.3)

    plt.tight_layout()
    plt.show()


# ===================================================================
# Electron configuration box diagram
# ===================================================================

def plot_electron_configuration(atom, figsize: tuple = None):
    """Draw an orbital box (arrow) diagram for an atom's electron config.

    Each orbital is a box. Spin-up electrons are shown as up-arrows,
    spin-down as down-arrows.
    """
    from molbuilder.atomic.quantum_atom import QuantumAtom

    subshells = atom.subshells
    if not subshells:
        print("No electrons to display.")
        return

    # Layout: one row per subshell, boxes for each m_l
    total_orbitals = sum(2 * ss.l + 1 for ss in subshells)
    if figsize is None:
        max_boxes = max(2 * ss.l + 1 for ss in subshells)
        figsize = (max(6, max_boxes * 1.2 + 3), len(subshells) * 0.9 + 1)

    fig, ax = plt.subplots(figsize=figsize, facecolor=BG_COLOR)
    ax.set_facecolor(BG_COLOR)
    ax.set_xlim(-1, max(2 * ss.l + 1 for ss in subshells) + 2)
    ax.set_ylim(-0.5, len(subshells) * 1.0 + 0.5)
    ax.axis("off")

    charge_label = ""
    if atom.charge > 0:
        charge_label = f" (+{atom.charge})"
    elif atom.charge < 0:
        charge_label = f" ({atom.charge})"
    ax.set_title(
        f"Electron Configuration -- {atom.name} ({atom.symbol}{charge_label})",
        color=TEXT_COLOR, fontsize=13, pad=10,
    )

    box_w, box_h = 0.8, 0.6

    for row_idx, ss in enumerate(reversed(subshells)):
        y = row_idx * 1.0 + 0.3
        num_orbitals = 2 * ss.l + 1
        ml_values = list(range(-ss.l, ss.l + 1))

        # Determine which orbitals have spin-up, spin-down
        states = ss.quantum_states()
        occupied = {}
        for st in states:
            occupied.setdefault(st.ml, []).append(st.ms)

        # Label
        ax.text(-0.8, y, ss.label, fontsize=11, color=TEXT_COLOR,
                ha="right", va="center", fontweight="bold")

        for j, ml in enumerate(ml_values):
            bx = j * 1.0 + 0.2

            # Draw box
            rect = plt.Rectangle((bx, y - box_h / 2), box_w, box_h,
                                  fill=False, edgecolor="#4466aa",
                                  linewidth=1.2)
            ax.add_patch(rect)

            # Draw arrows for electrons
            spins = occupied.get(ml, [])
            for k, ms in enumerate(sorted(spins, reverse=True)):
                arrow_x = bx + box_w * (0.3 + 0.4 * k)
                if ms > 0:
                    ax.annotate("", xy=(arrow_x, y + 0.2),
                                xytext=(arrow_x, y - 0.15),
                                arrowprops=dict(arrowstyle="->",
                                                color=POSITIVE_COLOR, lw=1.8))
                else:
                    ax.annotate("", xy=(arrow_x, y - 0.15),
                                xytext=(arrow_x, y + 0.2),
                                arrowprops=dict(arrowstyle="->",
                                                color=NEGATIVE_COLOR, lw=1.8))

    plt.tight_layout()
    plt.show()


# ===================================================================
# Energy level diagram
# ===================================================================

def plot_energy_levels(atom, figsize: tuple = (7, 8)):
    """Plot an energy level diagram for the atom's subshells."""
    from molbuilder.atomic.quantum_atom import QuantumAtom

    subshells = atom.subshells
    if not subshells:
        print("No subshells to display.")
        return

    fig, ax = plt.subplots(figsize=figsize, facecolor=BG_COLOR)
    ax.set_facecolor(BG_COLOR)

    # Compute energies
    energies = []
    labels = []
    counts = []
    for ss in subshells:
        e = atom.orbital_energy_eV(ss.n, ss.l)
        energies.append(e)
        labels.append(ss.label)
        counts.append(ss.electron_count)

    # Group by n for horizontal positioning
    n_values = sorted(set(ss.n for ss in subshells))
    n_to_x = {n: i * 2.0 for i, n in enumerate(n_values)}

    for i, ss in enumerate(subshells):
        e = energies[i]
        x_center = n_to_x[ss.n] + ss.l * 0.5
        x_left = x_center - 0.3
        x_right = x_center + 0.3

        ax.plot([x_left, x_right], [e, e], color=ENERGY_COLOR,
                linewidth=2.5, solid_capstyle="round")
        ax.text(x_center, e + 0.3, f"{labels[i]}",
                ha="center", va="bottom", fontsize=9, color=TEXT_COLOR)
        ax.text(x_center, e - 0.5, f"({counts[i]}e-)",
                ha="center", va="top", fontsize=7, color="#7799bb")

    ax.set_ylabel("Energy (eV)", color=TEXT_COLOR, fontsize=12)
    ax.set_title(
        f"Orbital Energy Levels -- {atom.name} ({atom.symbol})",
        color=TEXT_COLOR, fontsize=13, pad=12,
    )
    ax.tick_params(axis="y", colors=TEXT_COLOR)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    for spine in ax.spines.values():
        spine.set_color(GRID_COLOR)
    ax.grid(True, axis="y", color=GRID_COLOR, alpha=0.3)

    plt.tight_layout()
    plt.show()


# ===================================================================
# Combined atom overview
# ===================================================================

def visualize_atom(atom, orbital_nlm: tuple = None):
    """Show a multi-panel overview of an atom.

    Panels:
        1. Radial probability distributions for all occupied subshells
        2. Energy level diagram
        3. Electron configuration box diagram
        4. 3D orbital for the specified (n, l, m) or outermost subshell

    Parameters
    ----------
    atom        : QuantumAtom instance
    orbital_nlm : (n, l, m) to visualize in 3D; defaults to outermost subshell m=0
    """
    from molbuilder.atomic.quantum_atom import QuantumAtom

    print(atom.summary())

    # Determine orbital to show in 3D
    if orbital_nlm is None:
        ss = atom.subshells[-1]
        orbital_nlm = (ss.n, ss.l, 0)

    n, l, m = orbital_nlm

    # 1) Radial probability
    n_l_pairs = list(set((ss.n, ss.l) for ss in atom.subshells))
    n_l_pairs.sort()
    plot_radial_probability(n_l_pairs, Z=atom.atomic_number)

    # 2) Energy levels
    plot_energy_levels(atom)

    # 3) Configuration diagram
    plot_electron_configuration(atom)

    # 4) 3D orbital
    plot_orbital_3d(n, l, m, Z=atom.atomic_number)
