"""Molecular Geometry Visualization

Renders 3D ball-and-stick models of VSEPR-predicted molecular
geometries using matplotlib.  Matches the dark theme from the
existing quantum visualization module.

Features:
    - Atom spheres with CPK colours
    - Single / double / triple bond lines
    - Lone pair lobes
    - Bond angle arcs with labels
    - Informational overlay text
    - Multi-molecule gallery view

Migrated from legacy/molecule_visualization.py.
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from molbuilder.core.element_properties import cpk_color, covalent_radius_pm
from molbuilder.bonding.vsepr import VSEPRMolecule
from molbuilder.visualization.theme import (
    BG_COLOR, TEXT_COLOR, GRID_COLOR,
    LONE_PAIR_COLOR, BOND_COLOR, ANGLE_ARC_COLOR,
)


# ===================================================================
# Helpers
# ===================================================================

def _perpendicular_offset(pos_a, pos_b, magnitude):
    """Compute a perpendicular offset vector for drawing multi-bonds."""
    bond_vec = pos_b - pos_a
    norm = np.linalg.norm(bond_vec)
    if norm < 1e-10:
        return np.array([magnitude, 0.0, 0.0])
    bond_dir = bond_vec / norm
    # Choose a reference not parallel to bond_dir
    ref = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(bond_dir, ref)) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])
    perp = np.cross(bond_dir, ref)
    perp = perp / np.linalg.norm(perp) * magnitude
    return perp


def _slerp(va, vb, t):
    """Spherical linear interpolation between unit vectors."""
    dot = np.clip(np.dot(va, vb), -1.0, 1.0)
    omega = math.acos(dot)
    if omega < 1e-6:
        return va * (1.0 - t) + vb * t
    return (math.sin((1 - t) * omega) * va + math.sin(t * omega) * vb) / math.sin(omega)


def _draw_angle_arc(ax, center, vec_a, vec_b, angle_deg,
                    radius=0.3, n_points=30, label=True):
    """Draw a circular arc between two bond vectors to show angle."""
    na = np.linalg.norm(vec_a)
    nb = np.linalg.norm(vec_b)
    if na < 1e-10 or nb < 1e-10:
        return
    va = vec_a / na
    vb = vec_b / nb

    arc_points = []
    for i in range(n_points + 1):
        t = i / n_points
        v = _slerp(va, vb, t)
        arc_points.append(center + radius * v)
    arc_points = np.array(arc_points)

    ax.plot(arc_points[:, 0], arc_points[:, 1], arc_points[:, 2],
            color=ANGLE_ARC_COLOR, linewidth=0.8, alpha=0.6)

    if label:
        mid = arc_points[n_points // 2]
        ax.text(mid[0], mid[1], mid[2], f" {angle_deg:.1f}",
                color=ANGLE_ARC_COLOR, fontsize=7, alpha=0.8)


# ===================================================================
# Single molecule visualisation
# ===================================================================

def visualize_molecule(molecule: VSEPRMolecule,
                       show_lone_pairs: bool = True,
                       show_labels: bool = True,
                       show_angles: bool = True,
                       figsize: tuple = (9, 8)):
    """Render a 3D ball-and-stick model of a molecule.

    Parameters
    ----------
    molecule       : VSEPRMolecule
    show_lone_pairs: draw lone pair lobes
    show_labels    : label each atom
    show_angles    : draw bond angle arcs
    figsize        : figure size
    """
    coords = molecule.coordinates
    atom_positions = coords['atom_positions']
    bonds = coords['bonds']
    lp_positions = coords['lone_pair_positions']
    central_idx = coords['central_index']

    fig = plt.figure(figsize=figsize, facecolor=BG_COLOR)
    ax = fig.add_subplot(111, projection='3d', facecolor=BG_COLOR)

    # ---- Bonds ----
    for idx_a, idx_b, order in bonds:
        sym_a, pos_a = atom_positions[idx_a]
        sym_b, pos_b = atom_positions[idx_b]
        if pos_a is None or pos_b is None:
            continue

        if order == 1:
            ax.plot([pos_a[0], pos_b[0]],
                    [pos_a[1], pos_b[1]],
                    [pos_a[2], pos_b[2]],
                    color=BOND_COLOR, linewidth=2.5, zorder=3)
        elif order == 2:
            offset = _perpendicular_offset(pos_a, pos_b, 0.06)
            for sign in [1, -1]:
                ax.plot([pos_a[0] + sign*offset[0], pos_b[0] + sign*offset[0]],
                        [pos_a[1] + sign*offset[1], pos_b[1] + sign*offset[1]],
                        [pos_a[2] + sign*offset[2], pos_b[2] + sign*offset[2]],
                        color=BOND_COLOR, linewidth=2.0, zorder=3)
        elif order == 3:
            ax.plot([pos_a[0], pos_b[0]],
                    [pos_a[1], pos_b[1]],
                    [pos_a[2], pos_b[2]],
                    color=BOND_COLOR, linewidth=2.5, zorder=3)
            offset = _perpendicular_offset(pos_a, pos_b, 0.07)
            for sign in [1, -1]:
                ax.plot([pos_a[0] + sign*offset[0], pos_b[0] + sign*offset[0]],
                        [pos_a[1] + sign*offset[1], pos_b[1] + sign*offset[1]],
                        [pos_a[2] + sign*offset[2], pos_b[2] + sign*offset[2]],
                        color=BOND_COLOR, linewidth=1.5, zorder=3)

    # ---- Atoms ----
    for i, (sym, pos) in enumerate(atom_positions):
        if sym is None or pos is None:
            continue
        color = cpk_color(sym)
        radius = covalent_radius_pm(sym)
        size = 180 + radius * 0.8
        if i == central_idx:
            size *= 1.15
        ax.scatter(pos[0], pos[1], pos[2],
                   c=color, s=size, edgecolors='white', linewidths=0.5,
                   alpha=0.92, zorder=5, depthshade=True)

    # ---- Atom labels ----
    if show_labels:
        for i, (sym, pos) in enumerate(atom_positions):
            if sym is None or pos is None:
                continue
            ax.text(pos[0] + 0.08, pos[1] + 0.08, pos[2] + 0.08,
                    sym, color=TEXT_COLOR, fontsize=11, fontweight='bold',
                    zorder=10)

    # ---- Lone pairs ----
    if show_lone_pairs and lp_positions:
        rng = np.random.default_rng(42)
        for atom_idx, direction in lp_positions:
            sym, atom_pos = atom_positions[atom_idx]
            if atom_pos is None:
                continue
            lobe_length = 0.45
            t = np.linspace(0.15, lobe_length, 25)
            lobe_pts = atom_pos[np.newaxis, :] + direction[np.newaxis, :] * t[:, np.newaxis]
            spread = 0.03
            lobe_pts += rng.normal(0, spread, lobe_pts.shape)
            ax.scatter(lobe_pts[:, 0], lobe_pts[:, 1], lobe_pts[:, 2],
                       c=LONE_PAIR_COLOR, s=12, alpha=0.4, zorder=4,
                       depthshade=True)

    # ---- Bond angles ----
    if show_angles:
        central_sym, central_pos = atom_positions[central_idx]
        if central_pos is not None:
            terminal_data = []
            for idx_a, idx_b, order in bonds:
                ti = idx_b if idx_a == central_idx else idx_a
                _, tpos = atom_positions[ti]
                if tpos is not None:
                    terminal_data.append(tpos)

            # Draw arcs for adjacent bond pairs (limit to avoid clutter)
            drawn = set()
            for i in range(len(terminal_data)):
                for j in range(i + 1, len(terminal_data)):
                    va = terminal_data[i] - central_pos
                    vb = terminal_data[j] - central_pos
                    na = np.linalg.norm(va)
                    nb = np.linalg.norm(vb)
                    if na < 1e-10 or nb < 1e-10:
                        continue
                    cos_a = np.clip(np.dot(va, vb) / (na * nb), -1, 1)
                    angle = math.degrees(math.acos(cos_a))
                    # Only draw ~90 or ~120 degree angles (skip 180)
                    angle_key = round(angle)
                    if angle_key > 170:
                        continue
                    if angle_key in drawn and len(terminal_data) > 3:
                        continue
                    drawn.add(angle_key)
                    _draw_angle_arc(ax, central_pos, va, vb, angle,
                                    radius=0.3, label=True)

    # ---- Axis limits ----
    all_pos = [pos for _, pos in atom_positions if pos is not None]
    if all_pos:
        all_pos = np.array(all_pos)
        max_range = np.max(np.abs(all_pos)) * 1.5
        max_range = max(max_range, 1.0)
    else:
        max_range = 2.0
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.set_zlim(-max_range, max_range)

    # ---- Style ----
    axe = molecule.axe
    ax.set_title(
        f"{molecule.formula} -- {axe.molecular_geometry} ({axe.axe_notation})",
        color=TEXT_COLOR, fontsize=14, pad=10,
    )
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.set_facecolor(BG_COLOR)
        pane.set_edgecolor(GRID_COLOR)
    ax.tick_params(colors=TEXT_COLOR, labelsize=7)
    ax.set_xlabel("x (A)", color=TEXT_COLOR, fontsize=8)
    ax.set_ylabel("y (A)", color=TEXT_COLOR, fontsize=8)
    ax.set_zlabel("z (A)", color=TEXT_COLOR, fontsize=8)

    # Info text box
    angles_str = ", ".join(f"{a:.0f}" for a in axe.ideal_bond_angles) if axe.ideal_bond_angles else "N/A"
    info = (
        f"Geometry: {axe.molecular_geometry}\n"
        f"Hybridization: {axe.hybridization}\n"
        f"Bond angle(s): {angles_str} deg\n"
        f"Lone pairs: {axe.lone_pairs}"
    )
    ax.text2D(0.02, 0.95, info, transform=ax.transAxes,
              color=TEXT_COLOR, fontsize=9, verticalalignment='top',
              fontfamily='monospace',
              bbox=dict(boxstyle='round', facecolor='#111122',
                        edgecolor=GRID_COLOR, alpha=0.85))

    plt.tight_layout()
    plt.show()


# ===================================================================
# Gallery view
# ===================================================================

def _render_on_axis(ax, molecule: VSEPRMolecule):
    """Render a molecule onto a given 3D axis (no plt.show)."""
    coords = molecule.coordinates
    atom_positions = coords['atom_positions']
    bonds = coords['bonds']
    lp_positions = coords['lone_pair_positions']
    central_idx = coords['central_index']

    # Bonds
    for idx_a, idx_b, order in bonds:
        _, pos_a = atom_positions[idx_a]
        _, pos_b = atom_positions[idx_b]
        if pos_a is None or pos_b is None:
            continue
        if order >= 2:
            offset = _perpendicular_offset(pos_a, pos_b, 0.05)
            for sign in ([0] if order == 1 else [1, -1]):
                ax.plot([pos_a[0]+sign*offset[0], pos_b[0]+sign*offset[0]],
                        [pos_a[1]+sign*offset[1], pos_b[1]+sign*offset[1]],
                        [pos_a[2]+sign*offset[2], pos_b[2]+sign*offset[2]],
                        color=BOND_COLOR, linewidth=1.5, zorder=3)
            if order == 3:
                ax.plot([pos_a[0], pos_b[0]],
                        [pos_a[1], pos_b[1]],
                        [pos_a[2], pos_b[2]],
                        color=BOND_COLOR, linewidth=1.8, zorder=3)
        else:
            ax.plot([pos_a[0], pos_b[0]],
                    [pos_a[1], pos_b[1]],
                    [pos_a[2], pos_b[2]],
                    color=BOND_COLOR, linewidth=1.8, zorder=3)

    # Atoms
    for i, (sym, pos) in enumerate(atom_positions):
        if sym is None or pos is None:
            continue
        color = cpk_color(sym)
        size = 100 + covalent_radius_pm(sym) * 0.4
        ax.scatter(pos[0], pos[1], pos[2],
                   c=color, s=size, edgecolors='white', linewidths=0.3,
                   alpha=0.9, zorder=5, depthshade=True)
        ax.text(pos[0]+0.05, pos[1]+0.05, pos[2]+0.05,
                sym, color=TEXT_COLOR, fontsize=7, fontweight='bold', zorder=10)

    # Lone pairs
    if lp_positions:
        rng = np.random.default_rng(42)
        for atom_idx, direction in lp_positions:
            _, atom_pos = atom_positions[atom_idx]
            if atom_pos is None:
                continue
            t = np.linspace(0.12, 0.35, 15)
            lobe_pts = atom_pos[np.newaxis, :] + direction[np.newaxis, :] * t[:, np.newaxis]
            lobe_pts += rng.normal(0, 0.025, lobe_pts.shape)
            ax.scatter(lobe_pts[:, 0], lobe_pts[:, 1], lobe_pts[:, 2],
                       c=LONE_PAIR_COLOR, s=6, alpha=0.35, zorder=4,
                       depthshade=True)

    # Limits
    all_pos = [pos for _, pos in atom_positions if pos is not None]
    if all_pos:
        max_r = np.max(np.abs(np.array(all_pos))) * 1.5
        max_r = max(max_r, 0.8)
    else:
        max_r = 1.5
    ax.set_xlim(-max_r, max_r)
    ax.set_ylim(-max_r, max_r)
    ax.set_zlim(-max_r, max_r)

    axe = molecule.axe
    ax.set_title(f"{molecule.formula}\n{axe.molecular_geometry}",
                 color=TEXT_COLOR, fontsize=9, pad=2)
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.set_facecolor(BG_COLOR)
        pane.set_edgecolor(GRID_COLOR)
    ax.tick_params(colors=TEXT_COLOR, labelsize=5)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_zlabel("")


def visualize_gallery(molecules: list[VSEPRMolecule], cols: int = 3,
                      figsize: tuple = (16, 14)):
    """Show multiple molecules in a grid layout for comparison.

    Parameters
    ----------
    molecules : list of VSEPRMolecule
    cols      : columns in the grid
    figsize   : figure size
    """
    n = len(molecules)
    rows = math.ceil(n / cols)

    fig = plt.figure(figsize=figsize, facecolor=BG_COLOR)
    fig.suptitle("VSEPR Molecular Geometry Gallery",
                 color=TEXT_COLOR, fontsize=16, y=0.98)

    for i, mol in enumerate(molecules):
        ax = fig.add_subplot(rows, cols, i + 1, projection='3d',
                             facecolor=BG_COLOR)
        _render_on_axis(ax, mol)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()
