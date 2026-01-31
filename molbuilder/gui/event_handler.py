"""Event handler: atom placement, selection, coordinate computation."""

from __future__ import annotations
import numpy as np
from molbuilder.molecule.graph import Molecule, Hybridization
from molbuilder.core.bond_data import bond_length, SP3_ANGLE
from molbuilder.core.geometry import normalize, available_tetrahedral_dirs


class MolEventHandler:
    """Manages molecule editing state and operations."""

    def __init__(self):
        self.mol = Molecule("untitled")
        self.selected_atoms: list[int] = []  # selected atom indices
        self.current_element = "C"
        self.current_bond_order = 1

    def add_atom_free(self, element: str = None) -> int:
        """Add an atom at a computed position.

        If no atoms exist, place at origin.
        If atoms are selected, place bonded to the first selected atom.
        Otherwise place bonded to the last atom.
        """
        sym = element or self.current_element

        if len(self.mol.atoms) == 0:
            idx = self.mol.add_atom(sym, [0.0, 0.0, 0.0], Hybridization.SP3)
            return idx

        # Determine attachment point
        if self.selected_atoms:
            parent = self.selected_atoms[0]
        else:
            parent = len(self.mol.atoms) - 1

        idx = self.mol.add_atom_bonded(
            sym, parent,
            bond_order=self.current_bond_order,
            hybridization=Hybridization.SP3,
        )
        return idx

    def add_bond_between_selected(self) -> bool:
        """Add a bond between two selected atoms."""
        if len(self.selected_atoms) < 2:
            return False
        i, j = self.selected_atoms[0], self.selected_atoms[1]
        existing = self.mol.get_bond(i, j)
        if existing:
            return False
        self.mol.add_bond(i, j, order=self.current_bond_order)
        return True

    def delete_selected(self):
        """Delete selected atoms (rebuild molecule without them)."""
        if not self.selected_atoms:
            return
        keep = [i for i in range(len(self.mol.atoms)) if i not in self.selected_atoms]
        if not keep:
            self.mol = Molecule(self.mol.name)
            self.selected_atoms.clear()
            return

        remap = {}
        new_mol = Molecule(self.mol.name)
        for new_i, old_i in enumerate(keep):
            a = self.mol.atoms[old_i]
            new_mol.add_atom(a.symbol, a.position.copy(), a.hybridization)
            remap[old_i] = new_i

        for bond in self.mol.bonds:
            if bond.atom_i in remap and bond.atom_j in remap:
                new_mol.add_bond(remap[bond.atom_i], remap[bond.atom_j],
                                 bond.order, bond.rotatable)

        self.mol = new_mol
        self.selected_atoms.clear()

    def add_hydrogens(self):
        """Add missing hydrogens to all heavy atoms (simple valence fill)."""
        from molbuilder.core.geometry import add_sp3_hydrogens

        target_valence = {"C": 4, "N": 3, "O": 2, "S": 2, "P": 3, "B": 3}
        for i, atom in enumerate(list(self.mol.atoms)):
            if atom.symbol == "H":
                continue
            target = target_valence.get(atom.symbol)
            if target is None:
                continue
            current = sum(b.order for b in self.mol.bonds
                          if b.atom_i == i or b.atom_j == i)
            needed = target - current
            if needed > 0 and atom.symbol == "C":
                add_sp3_hydrogens(self.mol, i, needed)
            elif needed > 0:
                # Generic hydrogen addition
                from molbuilder.core.bond_data import bond_length as bl
                pos = atom.position + np.array([0.0, 0.0, bl(atom.symbol, "H", 1)])
                for _ in range(needed):
                    h_idx = self.mol.add_atom_bonded("H", i, bond_order=1, rotatable=False)

    def clear(self):
        """Remove all atoms."""
        self.mol = Molecule(self.mol.name)
        self.selected_atoms.clear()

    def toggle_selection(self, atom_idx: int):
        """Toggle selection of an atom."""
        if atom_idx in self.selected_atoms:
            self.selected_atoms.remove(atom_idx)
        else:
            self.selected_atoms.append(atom_idx)

    def clear_selection(self):
        """Clear atom selection."""
        self.selected_atoms.clear()
