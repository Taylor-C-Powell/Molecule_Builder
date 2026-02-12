"""Molecule endpoints: parse SMILES, get properties, get 3D structure."""

from fastapi import APIRouter, Depends
from app.dependencies import UserContext, check_rate_limit
from app.exceptions import InvalidSMILES, MoleculeNotFound
from app.models.molecule import (
    ParseSmilesRequest,
    MoleculeResponse,
    MoleculePropertiesResponse,
    Molecule3DResponse,
)
from app.services.molecule_store import molecule_store
from app.services.molecule_service import parse_smiles, build_properties, build_3d

router = APIRouter(prefix="/api/v1/molecule", tags=["molecule"])


@router.post("/from-smiles", response_model=MoleculeResponse)
def from_smiles(
    body: ParseSmilesRequest,
    user: UserContext = Depends(check_rate_limit),
):
    try:
        mol, canonical = parse_smiles(body.smiles, body.name)
    except (ValueError, KeyError) as e:
        raise InvalidSMILES(body.smiles, str(e))
    mol_id = molecule_store.put(mol, canonical)
    return MoleculeResponse(
        id=mol_id,
        name=mol.name,
        smiles=canonical,
        num_atoms=len(mol.atoms),
        num_bonds=len(mol.bonds),
    )


@router.get("/{mol_id}/properties", response_model=MoleculePropertiesResponse)
def get_properties(
    mol_id: str,
    user: UserContext = Depends(check_rate_limit),
):
    item = molecule_store.get(mol_id)
    if item is None:
        raise MoleculeNotFound(mol_id)
    mol, smiles = item
    return build_properties(mol_id, mol, smiles)


@router.get("/{mol_id}/3d", response_model=Molecule3DResponse)
def get_3d(
    mol_id: str,
    user: UserContext = Depends(check_rate_limit),
):
    item = molecule_store.get(mol_id)
    if item is None:
        raise MoleculeNotFound(mol_id)
    mol, _ = item
    return build_3d(mol_id, mol)
