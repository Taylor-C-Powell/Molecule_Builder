"""
VSEPR (Valence Shell Electron Pair Repulsion) Theory

Predicts molecular geometry from electron pair repulsion around a
central atom.  Given a molecular formula, this module:

    1. Builds a Lewis structure (bonding pairs + lone pairs)
    2. Classifies as AXnEm  (A=central, X=bonds, E=lone pairs)
    3. Predicts electron geometry, molecular geometry, hybridisation
    4. Generates 3D atomic coordinates for visualisation
"""

import math
import numpy as np
from dataclasses import dataclass

from molbuilder.core.element_properties import estimated_bond_length_pm
from molbuilder.bonding.lewis import LewisStructure


# ===================================================================
# Geometry reference tables
# ===================================================================

ELECTRON_GEOMETRY = {
    1: "terminal",
    2: "linear",
    3: "trigonal planar",
    4: "tetrahedral",
    5: "trigonal bipyramidal",
    6: "octahedral",
    7: "pentagonal bipyramidal",
}

# (steric_number, lone_pairs) -> (molecular_geometry, ideal_bond_angles)
MOLECULAR_GEOMETRY = {
    (1, 0): ("terminal",               [0.0]),
    (2, 0): ("linear",                 [180.0]),
    (2, 1): ("linear",                 [180.0]),
    (3, 0): ("trigonal planar",        [120.0]),
    (3, 1): ("bent",                   [117.0]),
    (4, 0): ("tetrahedral",            [109.5]),
    (4, 1): ("trigonal pyramidal",     [107.0]),
    (4, 2): ("bent",                   [104.5]),
    (5, 0): ("trigonal bipyramidal",   [90.0, 120.0]),
    (5, 1): ("seesaw",                 [90.0, 117.0]),
    (5, 2): ("T-shaped",              [90.0]),
    (5, 3): ("linear",                 [180.0]),
    (6, 0): ("octahedral",            [90.0]),
    (6, 1): ("square pyramidal",       [90.0]),
    (6, 2): ("square planar",         [90.0]),
    (6, 3): ("T-shaped",              [90.0]),
    (6, 4): ("linear",                 [180.0]),
    (7, 0): ("pentagonal bipyramidal", [72.0, 90.0]),
    (7, 1): ("pentagonal pyramidal",   [72.0, 90.0]),
    (7, 2): ("pentagonal planar",      [72.0]),
}

HYBRIDIZATION = {
    1: "s",
    2: "sp",
    3: "sp2",
    4: "sp3",
    5: "sp3d",
    6: "sp3d2",
    7: "sp3d3",
}


# ===================================================================
# AXE classification
# ===================================================================

@dataclass
class AXEClassification:
    """AXnEm classification of a molecule's central atom."""
    central_symbol: str
    bonding_groups: int       # X
    lone_pairs: int           # E
    steric_number: int        # X + E
    axe_notation: str         # e.g. "AX2E2"
    electron_geometry: str
    molecular_geometry: str
    ideal_bond_angles: list
    hybridization: str


# ===================================================================
# Ideal direction vectors for each electron geometry
# ===================================================================

def _linear_directions():
    return [
        np.array([0.0, 0.0,  1.0]),
        np.array([0.0, 0.0, -1.0]),
    ]


def _trigonal_planar_directions():
    return [
        np.array([1.0, 0.0, 0.0]),
        np.array([-0.5,  math.sqrt(3)/2, 0.0]),
        np.array([-0.5, -math.sqrt(3)/2, 0.0]),
    ]


def _tetrahedral_directions():
    return [
        np.array([ 1.0,  1.0,  1.0]) / math.sqrt(3),
        np.array([ 1.0, -1.0, -1.0]) / math.sqrt(3),
        np.array([-1.0,  1.0, -1.0]) / math.sqrt(3),
        np.array([-1.0, -1.0,  1.0]) / math.sqrt(3),
    ]


def _trigonal_bipyramidal_directions():
    """3 equatorial (xy-plane, 120 deg) + 2 axial (+z, -z)."""
    eq = _trigonal_planar_directions()
    ax = [np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, -1.0])]
    # Return equatorial first, then axial (lone pairs prefer equatorial)
    return eq + ax


def _octahedral_directions():
    return [
        np.array([ 1.0,  0.0,  0.0]),
        np.array([-1.0,  0.0,  0.0]),
        np.array([ 0.0,  1.0,  0.0]),
        np.array([ 0.0, -1.0,  0.0]),
        np.array([ 0.0,  0.0,  1.0]),
        np.array([ 0.0,  0.0, -1.0]),
    ]


def _pentagonal_bipyramidal_directions():
    """5 equatorial (xy-plane, 72 deg) + 2 axial."""
    eq = []
    for i in range(5):
        angle = 2 * math.pi * i / 5
        eq.append(np.array([math.cos(angle), math.sin(angle), 0.0]))
    ax = [np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, -1.0])]
    return eq + ax


_DIRECTION_GENERATORS = {
    2: _linear_directions,
    3: _trigonal_planar_directions,
    4: _tetrahedral_directions,
    5: _trigonal_bipyramidal_directions,
    6: _octahedral_directions,
    7: _pentagonal_bipyramidal_directions,
}


def _assign_positions(steric_number: int, num_lone_pairs: int):
    """Assign direction vectors to bonding groups and lone pairs.

    Lone pairs are placed at positions that minimise 90-degree repulsions:
        SN=3: lone pair at any equatorial position
        SN=4: lone pair at any tetrahedral vertex
        SN=5: lone pairs at equatorial positions first
        SN=6: 1 LP at any position; 2 LPs at trans (opposite) positions

    Returns
    -------
    bonding_dirs : list of np.array  -- unit vectors for bonded atoms
    lone_pair_dirs : list of np.array -- unit vectors for lone pairs
    """
    gen = _DIRECTION_GENERATORS.get(steric_number)
    if gen is None:
        # Single atom or unknown: just put along +z
        return [np.array([0.0, 0.0, 1.0])] * (steric_number - num_lone_pairs), \
               [np.array([0.0, 0.0, -1.0])] * num_lone_pairs

    all_dirs = gen()

    if steric_number == 5:
        # Equatorial positions are indices 0,1,2; axial are 3,4
        # Lone pairs prefer equatorial
        lp_indices = list(range(min(num_lone_pairs, 3)))
        if num_lone_pairs > 3:
            lp_indices.extend([3, 4][:num_lone_pairs - 3])
    elif steric_number == 6:
        # For 1 LP: remove index 0 (arbitrary, all equivalent)
        # For 2 LPs: remove trans pair (indices 0,1 are +x/-x)
        if num_lone_pairs == 1:
            lp_indices = [0]
        elif num_lone_pairs == 2:
            lp_indices = [0, 1]  # trans pair
        elif num_lone_pairs == 3:
            lp_indices = [0, 1, 2]  # fac arrangement
        elif num_lone_pairs == 4:
            lp_indices = [0, 1, 2, 3]
        else:
            lp_indices = list(range(num_lone_pairs))
    else:
        # For SN=2,3,4: lone pairs take the last positions
        lp_indices = list(range(steric_number - num_lone_pairs, steric_number))

    bonding_dirs = [all_dirs[i] for i in range(len(all_dirs)) if i not in lp_indices]
    lone_pair_dirs = [all_dirs[i] for i in lp_indices]

    # Apply lone-pair compression: when lone pairs are present, bonding
    # groups are pushed closer together.  Use the ideal bond angles from
    # MOLECULAR_GEOMETRY to compute the compression.
    if num_lone_pairs > 0 and len(bonding_dirs) >= 2:
        key = (steric_number, num_lone_pairs)
        if key in MOLECULAR_GEOMETRY:
            target_angle_deg = MOLECULAR_GEOMETRY[key][1][0]
            # Current angle between first two bonding directions
            d0 = bonding_dirs[0] / np.linalg.norm(bonding_dirs[0])
            d1 = bonding_dirs[1] / np.linalg.norm(bonding_dirs[1])
            cos_current = float(np.clip(np.dot(d0, d1), -1.0, 1.0))
            current_angle = math.degrees(math.acos(cos_current))
            if abs(current_angle - target_angle_deg) > 0.5:
                # Centroid of bonding directions
                centroid = sum(bonding_dirs) / len(bonding_dirs)
                c_norm = np.linalg.norm(centroid)
                if c_norm > 1e-8:
                    centroid = centroid / c_norm
                    # Iteratively scale compression toward centroid
                    half_target = math.radians(target_angle_deg / 2.0)
                    half_current = math.radians(current_angle / 2.0)
                    # Blend factor: how much to move toward centroid
                    # sin(half_target)/sin(half_current) gives the
                    # ratio needed, but a simple slerp-like blend works
                    t = 1.0 - math.sin(half_target) / max(math.sin(half_current), 1e-8)
                    t = max(0.0, min(t, 0.5))
                    compressed = []
                    for d in bonding_dirs:
                        blended = d * (1.0 - t) + centroid * t
                        blended = blended / np.linalg.norm(blended)
                        compressed.append(blended)
                    bonding_dirs = compressed

    return bonding_dirs, lone_pair_dirs


# ===================================================================
# 3D coordinate generation
# ===================================================================

def generate_3d_coordinates(lewis: LewisStructure) -> dict:
    """Generate 3D atomic coordinates from a Lewis structure.

    Central atom is placed at the origin.  Terminal atoms are placed
    along ideal geometry direction vectors, scaled by estimated bond
    lengths.

    Returns
    -------
    dict with keys:
        'atom_positions'     : list of (symbol, np.array([x,y,z])) in Angstroms
        'bonds'              : list of (idx_a, idx_b, bond_order)
        'lone_pair_positions': list of (atom_idx, np.array direction)
        'central_index'      : int
    """
    sn = lewis.steric_number()
    lp_count = lewis.lone_pairs_on_central()
    bp_count = lewis.bonding_pairs_on_central()

    # Handle single-atom or diatomic edge case
    if sn < 2 and len(lewis.atoms) == 1:
        return {
            'atom_positions': [(lewis.atoms[0], np.array([0.0, 0.0, 0.0]))],
            'bonds': [],
            'lone_pair_positions': [],
            'central_index': 0,
        }

    # Get direction vectors
    if sn >= 2:
        bonding_dirs, lp_dirs = _assign_positions(sn, lp_count)
    else:
        # Diatomic with SN=1
        bonding_dirs = [np.array([0.0, 0.0, 1.0])]
        lp_dirs = []

    # Build atom positions
    atom_positions = [(None, None)] * len(lewis.atoms)
    atom_positions[lewis.central_index] = (
        lewis.central_symbol,
        np.array([0.0, 0.0, 0.0]),
    )

    # Assign each bonded terminal to a direction
    bonds_out = []
    bond_dir_idx = 0
    for bond in lewis.bonds:
        ti = bond.atom_b if bond.atom_a == lewis.central_index else bond.atom_a
        sym = lewis.atoms[ti]

        bl = estimated_bond_length_pm(lewis.central_symbol, sym, bond.order) / 100.0

        if bond_dir_idx < len(bonding_dirs):
            direction = bonding_dirs[bond_dir_idx]
        else:
            direction = np.array([0.0, 0.0, 1.0])
        bond_dir_idx += 1

        pos = direction * bl
        atom_positions[ti] = (sym, pos)
        bonds_out.append((lewis.central_index, ti, bond.order))

    # Lone pair direction info (for visualisation)
    lp_positions = []
    for i, lp_dir in enumerate(lp_dirs):
        lp_positions.append((lewis.central_index, lp_dir))

    return {
        'atom_positions': atom_positions,
        'bonds': bonds_out,
        'lone_pair_positions': lp_positions,
        'central_index': lewis.central_index,
    }


# ===================================================================
# VSEPRMolecule
# ===================================================================

class VSEPRMolecule:
    """Complete VSEPR analysis of a molecule.

    Parameters
    ----------
    formula : str
        Molecular formula (e.g., 'H2O', 'CH4').
    charge : int
        Net charge (default 0).
    """

    def __init__(self, formula: str, charge: int = 0):
        self.formula = formula
        self.charge = charge
        self.lewis = LewisStructure(formula, charge)
        self.axe = self._classify()
        self.coordinates = generate_3d_coordinates(self.lewis)

    def _classify(self) -> AXEClassification:
        """Perform AXnEm classification."""
        X = self.lewis.bonding_pairs_on_central()
        E = self.lewis.lone_pairs_on_central()
        sn = X + E

        notation = f"AX{X}" + (f"E{E}" if E > 0 else "")
        e_geom = ELECTRON_GEOMETRY.get(sn, "unknown")

        entry = MOLECULAR_GEOMETRY.get((sn, E))
        if entry:
            m_geom, angles = entry
        else:
            m_geom, angles = "unknown", []

        hybrid = HYBRIDIZATION.get(sn, "unknown")

        return AXEClassification(
            central_symbol=self.lewis.central_symbol,
            bonding_groups=X,
            lone_pairs=E,
            steric_number=sn,
            axe_notation=notation,
            electron_geometry=e_geom,
            molecular_geometry=m_geom,
            ideal_bond_angles=angles,
            hybridization=hybrid,
        )

    def computed_bond_angles(self) -> list[float]:
        """Compute actual bond angles from 3D coordinates (degrees)."""
        coords = self.coordinates
        central_pos = coords['atom_positions'][coords['central_index']][1]
        terminal_positions = []
        for idx_a, idx_b, order in coords['bonds']:
            ti = idx_b if idx_a == coords['central_index'] else idx_a
            terminal_positions.append(coords['atom_positions'][ti][1])

        angles = []
        for i in range(len(terminal_positions)):
            for j in range(i + 1, len(terminal_positions)):
                va = terminal_positions[i] - central_pos
                vb = terminal_positions[j] - central_pos
                na, nb = np.linalg.norm(va), np.linalg.norm(vb)
                if na < 1e-10 or nb < 1e-10:
                    continue
                cos_angle = np.clip(np.dot(va, vb) / (na * nb), -1.0, 1.0)
                angles.append(math.degrees(math.acos(cos_angle)))
        return sorted(angles)

    # ------ display ------

    def __repr__(self):
        return (f"VSEPRMolecule({self.formula}, "
                f"{self.axe.axe_notation}, "
                f"{self.axe.molecular_geometry})")

    def summary(self) -> str:
        """Return a detailed ASCII summary."""
        axe = self.axe
        coords = self.coordinates

        angles_str = ", ".join(f"{a:.1f}" for a in axe.ideal_bond_angles) if axe.ideal_bond_angles else "N/A"
        computed = self.computed_bond_angles()
        computed_str = ", ".join(f"{a:.1f}" for a in computed) if computed else "N/A"

        lines = [
            f"{'='*60}",
            f"  VSEPR Analysis: {self.formula}",
            f"  Charge: {self.charge:+d}",
            f"{'='*60}",
            f"  Lewis Structure:",
            f"    Total valence electrons: {self.lewis.total_valence_electrons}",
            f"    Central atom: {self.lewis.central_symbol}",
        ]

        # Bonds summary
        bond_strs = []
        for bond in self.lewis.bonds:
            sym_a = self.lewis.atoms[bond.atom_a]
            sym_b = self.lewis.atoms[bond.atom_b]
            order_sym = {1: "-", 2: "=", 3: "#"}.get(bond.order, "?")
            bond_strs.append(f"{sym_a}{order_sym}{sym_b}")
        lines.append(f"    Bonds: {', '.join(bond_strs)}")

        lp_central = self.lewis.lone_pairs_on_central()
        lines.append(f"    Lone pairs on {self.lewis.central_symbol}: {lp_central}")

        lines.extend([
            f"{'='*60}",
            f"  VSEPR Classification:",
            f"    AXE notation:        {axe.axe_notation}",
            f"    Steric number:       {axe.steric_number}",
            f"    Electron geometry:   {axe.electron_geometry}",
            f"    Molecular geometry:  {axe.molecular_geometry}",
            f"    Hybridization:       {axe.hybridization}",
            f"    Ideal bond angles:   {angles_str} deg",
            f"    Computed angles:     {computed_str} deg",
        ])

        lines.append(f"{'='*60}")
        lines.append("  3D Coordinates (Angstroms):")
        for i, (sym, pos) in enumerate(coords['atom_positions']):
            if sym is not None and pos is not None:
                tag = " <-- central" if i == coords['central_index'] else ""
                lines.append(
                    f"    {sym:<4} {pos[0]:>8.4f}  {pos[1]:>8.4f}  {pos[2]:>8.4f}{tag}"
                )
        lines.append(f"{'='*60}")
        return "\n".join(lines)


# ===================================================================
# Convenience constructors
# ===================================================================

def from_formula(formula: str, charge: int = 0) -> VSEPRMolecule:
    """Create a VSEPRMolecule from a molecular formula string."""
    return VSEPRMolecule(formula, charge)


def from_atoms(symbols: list[str], charge: int = 0) -> VSEPRMolecule:
    """Create a VSEPRMolecule from a list of element symbols."""
    from collections import Counter
    counts = Counter(symbols)
    formula = "".join(f"{sym}{cnt if cnt > 1 else ''}" for sym, cnt in counts.items())
    return VSEPRMolecule(formula, charge)
