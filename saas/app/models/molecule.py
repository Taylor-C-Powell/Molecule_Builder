"""Molecule request/response models."""

from pydantic import BaseModel, Field


class ParseSmilesRequest(BaseModel):
    smiles: str = Field(..., description="SMILES string to parse")
    name: str = Field("", description="Optional molecule name")


class AtomResponse(BaseModel):
    index: int
    symbol: str
    position: list[float] = Field(..., description="[x, y, z] coordinates")
    hybridization: str | None = None
    formal_charge: int = 0


class BondResponse(BaseModel):
    atom_i: int
    atom_j: int
    order: int
    rotatable: bool


class MoleculeResponse(BaseModel):
    id: str
    name: str
    smiles: str
    num_atoms: int
    num_bonds: int


class MoleculePropertiesResponse(BaseModel):
    id: str
    smiles: str
    formula: str
    molecular_weight: float
    num_atoms: int
    num_bonds: int
    functional_groups: list[str]


class Molecule3DResponse(BaseModel):
    id: str
    atoms: list[AtomResponse]
    bonds: list[BondResponse]
