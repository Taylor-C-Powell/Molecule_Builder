"""Molecule request/response models."""

from pydantic import BaseModel, Field


class ParseSmilesRequest(BaseModel):
    smiles: str = Field(..., description="SMILES string to parse", min_length=1, max_length=2000)
    name: str = Field("", description="Optional molecule name", max_length=200)


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
    logp: float | None = None
    hbd: int | None = None
    hba: int | None = None
    rotatable_bonds: int | None = None
    tpsa: float | None = None
    heavy_atom_count: int | None = None
    lipinski_violations: int | None = None
    lipinski_pass: bool | None = None
    sa_score: float | None = None


class ADMETResponse(BaseModel):
    id: str
    smiles: str
    oral_bioavailability: str
    intestinal_absorption: str
    caco2_permeability: str
    pgp_substrate: bool
    bbb_penetrant: bool
    plasma_protein_binding: str
    vd_class: str
    cyp_inhibition: dict[str, bool]
    metabolic_stability: str
    renal_clearance: str
    half_life_class: str
    herg_risk: str
    ames_mutagenicity: bool
    hepatotoxicity_risk: str
    structural_alerts: list[str]
    overall_score: float
    warnings: list[str]
    flags: list[str]


class SolubilityResponse(BaseModel):
    id: str
    smiles: str
    log_s_esol: float
    log_s_gse: float
    solubility_mg_ml: float
    solubility_class: str
    estimated_melting_point_c: float
    crystallization_risk: str
    polymorph_risk: str


class Molecule3DResponse(BaseModel):
    id: str
    atoms: list[AtomResponse]
    bonds: list[BondResponse]
