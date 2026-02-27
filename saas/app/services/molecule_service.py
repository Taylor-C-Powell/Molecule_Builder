"""Molecule parsing and property extraction."""

from molbuilder.smiles.parser import parse
from molbuilder.smiles.writer import to_smiles
from molbuilder.molecule.graph import Molecule
from molbuilder.io.json_io import _molecular_formula
from molbuilder.reactions.functional_group_detect import detect_functional_groups
from molbuilder.core.elements import atomic_weight
from molbuilder.molecule.properties import lipinski_properties
from molbuilder.molecule.sa_score import sa_score
from app.models.molecule import (
    AtomResponse,
    BondResponse,
    MoleculePropertiesResponse,
    Molecule3DResponse,
)


def parse_smiles(smiles: str, name: str = "") -> tuple[Molecule, str]:
    """Parse SMILES and return (Molecule, canonical_smiles). Raises on invalid input."""
    mol = parse(smiles)
    mol.name = name or smiles
    canonical = to_smiles(mol)
    return mol, canonical


def get_formula(mol: Molecule) -> str:
    return _molecular_formula(mol)


def get_molecular_weight(mol: Molecule) -> float:
    total = 0.0
    for atom in mol.atoms:
        total += atomic_weight(atom.symbol)
    return round(total, 4)


def get_functional_groups(mol: Molecule) -> list[str]:
    fgs = detect_functional_groups(mol)
    return [fg.name for fg in fgs]


def build_properties(mol_id: str, mol: Molecule, smiles: str) -> MoleculePropertiesResponse:
    props = lipinski_properties(mol)
    sa = sa_score(mol)
    return MoleculePropertiesResponse(
        id=mol_id,
        smiles=smiles,
        formula=get_formula(mol),
        molecular_weight=get_molecular_weight(mol),
        num_atoms=len(mol.atoms),
        num_bonds=len(mol.bonds),
        functional_groups=get_functional_groups(mol),
        logp=props.logp,
        hbd=props.hbd,
        hba=props.hba,
        rotatable_bonds=props.rotatable_bonds,
        tpsa=props.tpsa,
        heavy_atom_count=props.heavy_atom_count,
        lipinski_violations=props.lipinski_violations,
        lipinski_pass=props.lipinski_pass,
        sa_score=sa.sa_score,
    )


def build_3d(mol_id: str, mol: Molecule) -> Molecule3DResponse:
    atoms = [
        AtomResponse(
            index=a.index,
            symbol=a.symbol,
            position=a.position.tolist(),
            hybridization=a.hybridization.name if a.hybridization else None,
            formal_charge=a.formal_charge,
        )
        for a in mol.atoms
    ]
    bonds = [
        BondResponse(
            atom_i=b.atom_i,
            atom_j=b.atom_j,
            order=b.order,
            rotatable=b.rotatable,
        )
        for b in mol.bonds
    ]
    return Molecule3DResponse(id=mol_id, atoms=atoms, bonds=bonds)
