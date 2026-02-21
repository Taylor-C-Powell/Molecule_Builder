"""Optional RDKit backend for high-quality 3D coordinate generation.

Uses RDKit's ETKDG + MMFF94 for superior geometry. Only available when
rdkit-pypi is installed. Converts MolBuilder topology to RDKit mol
via explicit atom/bond construction (not SMILES round-trip) to preserve
atom ordering.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from molbuilder.molecule.graph import Molecule


def generate_rdkit(mol: Molecule, seed: int = 42) -> None:
    """Generate 3D coordinates using RDKit ETKDG + MMFF94.

    Modifies atom positions in-place.

    Parameters
    ----------
    mol : Molecule
        Molecule with atoms and bonds defined.
    seed : int
        Random seed for ETKDG embedding.

    Raises
    ------
    ImportError
        If rdkit is not installed.
    RuntimeError
        If RDKit embedding fails.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        raise ImportError(
            "RDKit is required for the 'rdkit' backend. "
            "Install with: pip install rdkit-pypi"
        )

    rd_mol = Chem.RWMol()

    # Add atoms - build explicit index mapping
    idx_map: dict[int, int] = {}
    for atom in mol.atoms:
        rd_atom = Chem.Atom(atom.symbol)
        if atom.formal_charge:
            rd_atom.SetFormalCharge(atom.formal_charge)
        if atom.chirality == "@":
            rd_atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif atom.chirality == "@@":
            rd_atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
        rd_idx = rd_mol.AddAtom(rd_atom)
        idx_map[atom.index] = rd_idx

    # Add bonds
    bond_type_map = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
    }
    for bond in mol.bonds:
        rd_i = idx_map[bond.atom_i]
        rd_j = idx_map[bond.atom_j]
        bt = bond_type_map.get(bond.order, Chem.BondType.SINGLE)
        rd_mol.AddBond(rd_i, rd_j, bt)

    rd_mol = rd_mol.GetMol()

    # Embed with ETKDG v3
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    status = AllChem.EmbedMolecule(rd_mol, params)
    if status != 0:
        # Fallback: try with less strict parameters
        params.useRandomCoords = True
        status = AllChem.EmbedMolecule(rd_mol, params)
        if status != 0:
            raise RuntimeError(
                f"RDKit embedding failed for molecule '{mol.name}'"
            )

    # MMFF94 optimization
    try:
        AllChem.MMFFOptimizeMolecule(rd_mol, maxIters=200)
    except Exception:
        pass  # Use unoptimized embedding if MMFF fails

    # Extract coordinates back into MolBuilder atoms
    conf = rd_mol.GetConformer()
    for mb_idx, rd_idx in idx_map.items():
        pos = conf.GetAtomPosition(rd_idx)
        mol.atoms[mb_idx].position = np.array([pos.x, pos.y, pos.z])
