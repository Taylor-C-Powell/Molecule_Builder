"""Embedded matplotlib 3D canvas for tkinter."""

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from molbuilder.core.element_properties import cpk_color, covalent_radius_pm
from molbuilder.visualization.theme import BG_COLOR, TEXT_COLOR, GRID_COLOR, BOND_COLOR


class MolCanvas3D:
    """3D molecule viewport embedded in a tkinter frame."""

    def __init__(self, parent_frame):
        self.parent = parent_frame
        self.fig = Figure(figsize=(7, 6), facecolor=BG_COLOR)
        self.ax = self.fig.add_subplot(111, projection="3d", facecolor=BG_COLOR)

        self.canvas = FigureCanvasTkAgg(self.fig, master=parent_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        self.toolbar = NavigationToolbar2Tk(self.canvas, parent_frame)
        self.toolbar.update()

        self._style_axes()
        self.mol = None
        self._pick_callback = None

    def _style_axes(self):
        """Apply dark theme to 3D axes."""
        for pane in [self.ax.xaxis.pane, self.ax.yaxis.pane, self.ax.zaxis.pane]:
            pane.set_facecolor(BG_COLOR)
            pane.set_edgecolor(GRID_COLOR)
        self.ax.tick_params(colors=TEXT_COLOR, labelsize=7)
        self.ax.set_xlabel("x (A)", color=TEXT_COLOR, fontsize=8)
        self.ax.set_ylabel("y (A)", color=TEXT_COLOR, fontsize=8)
        self.ax.set_zlabel("z (A)", color=TEXT_COLOR, fontsize=8)

    def set_pick_callback(self, callback):
        """Set callback for atom picking: callback(atom_index)."""
        self._pick_callback = callback

    def render(self, mol):
        """Render a Molecule object."""
        self.mol = mol
        self.ax.clear()
        self._style_axes()

        if mol is None or len(mol.atoms) == 0:
            self.canvas.draw()
            return

        # Draw bonds
        for bond in mol.bonds:
            pa = mol.atoms[bond.atom_i].position
            pb = mol.atoms[bond.atom_j].position
            if bond.order == 1:
                self.ax.plot([pa[0], pb[0]], [pa[1], pb[1]], [pa[2], pb[2]],
                             color=BOND_COLOR, linewidth=2.0, zorder=3)
            elif bond.order == 2:
                mid = (pa + pb) / 2
                perp = self._perp_offset(pa, pb, 0.06)
                for s in [1, -1]:
                    self.ax.plot([pa[0] + s * perp[0], pb[0] + s * perp[0]],
                                [pa[1] + s * perp[1], pb[1] + s * perp[1]],
                                [pa[2] + s * perp[2], pb[2] + s * perp[2]],
                                color=BOND_COLOR, linewidth=1.8, zorder=3)
            elif bond.order == 3:
                self.ax.plot([pa[0], pb[0]], [pa[1], pb[1]], [pa[2], pb[2]],
                             color=BOND_COLOR, linewidth=2.0, zorder=3)
                perp = self._perp_offset(pa, pb, 0.07)
                for s in [1, -1]:
                    self.ax.plot([pa[0] + s * perp[0], pb[0] + s * perp[0]],
                                [pa[1] + s * perp[1], pb[1] + s * perp[1]],
                                [pa[2] + s * perp[2], pb[2] + s * perp[2]],
                                color=BOND_COLOR, linewidth=1.2, zorder=3)

        # Draw atoms
        for atom in mol.atoms:
            color = cpk_color(atom.symbol)
            size = 120 + covalent_radius_pm(atom.symbol) * 0.6
            self.ax.scatter(*atom.position, c=color, s=size,
                            edgecolors="white", linewidths=0.3,
                            alpha=0.92, zorder=5, depthshade=True)
            self.ax.text(atom.position[0] + 0.06, atom.position[1] + 0.06,
                         atom.position[2] + 0.06, atom.symbol,
                         color=TEXT_COLOR, fontsize=8, fontweight="bold", zorder=10)

        # Auto-scale
        all_pos = np.array([a.position for a in mol.atoms])
        if len(all_pos) > 0:
            mx = max(np.max(np.abs(all_pos)) * 1.4, 1.0)
            self.ax.set_xlim(-mx, mx)
            self.ax.set_ylim(-mx, mx)
            self.ax.set_zlim(-mx, mx)

        self.ax.set_title(mol.name or "Molecule", color=TEXT_COLOR, fontsize=11)
        self.canvas.draw()

    @staticmethod
    def _perp_offset(pa, pb, mag):
        """Compute a perpendicular offset vector for drawing multi-bonds."""
        d = pb - pa
        n = np.linalg.norm(d)
        if n < 1e-10:
            return np.array([mag, 0, 0])
        d = d / n
        ref = np.array([1.0, 0, 0])
        if abs(np.dot(d, ref)) > 0.9:
            ref = np.array([0, 1.0, 0])
        p = np.cross(d, ref)
        return p / np.linalg.norm(p) * mag
