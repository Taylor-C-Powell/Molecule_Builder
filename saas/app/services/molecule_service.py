"""Molecule parsing and property extraction."""

from functools import lru_cache

from molbuilder.smiles.parser import parse
from molbuilder.smiles.writer import to_smiles
from molbuilder.molecule.graph import Molecule
from molbuilder.io.json_io import _molecular_formula
from molbuilder.reactions.functional_group_detect import detect_functional_groups
from molbuilder.reactions.fg_smarts_validation import cross_validate_fg
from molbuilder.core.elements import atomic_weight
from molbuilder.molecule.properties import lipinski_properties
from molbuilder.molecule.sa_score import sa_score
from molbuilder.molecule.solubility import predict_solubility
from molbuilder.molecule.admet import predict_admet
from app.models.molecule import (
    ADMETResponse,
    AtomResponse,
    BondResponse,
    MoleculePropertiesResponse,
    Molecule3DResponse,
    SolubilityResponse,
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


def get_fg_confidence(mol: Molecule) -> float:
    """Return heuristic/SMARTS agreement confidence (0.0-1.0)."""
    result = cross_validate_fg(mol)
    return round(result.confidence, 3)


@lru_cache(maxsize=512)
def _compute_properties(smiles: str) -> dict:
    """Compute and cache molecule properties keyed by canonical SMILES."""
    mol = parse(smiles)
    props = lipinski_properties(mol)
    sa = sa_score(mol)
    return {
        "formula": get_formula(mol),
        "molecular_weight": get_molecular_weight(mol),
        "num_atoms": len(mol.atoms),
        "num_bonds": len(mol.bonds),
        "functional_groups": get_functional_groups(mol),
        "logp": props.logp,
        "hbd": props.hbd,
        "hba": props.hba,
        "rotatable_bonds": props.rotatable_bonds,
        "tpsa": props.tpsa,
        "heavy_atom_count": props.heavy_atom_count,
        "lipinski_violations": props.lipinski_violations,
        "lipinski_pass": props.lipinski_pass,
        "sa_score": sa.sa_score,
    }


def build_properties(mol_id: str, mol: Molecule, smiles: str) -> MoleculePropertiesResponse:
    cached = _compute_properties(smiles)
    return MoleculePropertiesResponse(id=mol_id, smiles=smiles, **cached)


@lru_cache(maxsize=512)
def _compute_solubility(smiles: str) -> dict:
    """Compute and cache solubility properties keyed by canonical SMILES."""
    mol = parse(smiles)
    result = predict_solubility(mol)
    return {
        "log_s_esol": result.log_s_esol,
        "log_s_gse": result.log_s_gse,
        "solubility_mg_ml": result.solubility_mg_ml,
        "solubility_class": result.solubility_class,
        "estimated_melting_point_c": result.estimated_melting_point_c,
        "crystallization_risk": result.crystallization_risk,
        "polymorph_risk": result.polymorph_risk,
    }


@lru_cache(maxsize=512)
def _compute_admet(smiles: str) -> dict:
    """Compute and cache ADMET profile keyed by canonical SMILES."""
    mol = parse(smiles)
    result = predict_admet(mol)
    return {
        "oral_bioavailability": result.oral_bioavailability,
        "intestinal_absorption": result.intestinal_absorption,
        "caco2_permeability": result.caco2_permeability,
        "pgp_substrate": result.pgp_substrate,
        "bbb_penetrant": result.bbb_penetrant,
        "plasma_protein_binding": result.plasma_protein_binding,
        "vd_class": result.vd_class,
        "cyp_inhibition": result.cyp_inhibition,
        "metabolic_stability": result.metabolic_stability,
        "renal_clearance": result.renal_clearance,
        "half_life_class": result.half_life_class,
        "herg_risk": result.herg_risk,
        "ames_mutagenicity": result.ames_mutagenicity,
        "hepatotoxicity_risk": result.hepatotoxicity_risk,
        "structural_alerts": result.structural_alerts,
        "overall_score": result.overall_score,
        "warnings": result.warnings,
        "flags": result.flags,
    }


def build_admet(mol_id: str, smiles: str) -> ADMETResponse:
    cached = _compute_admet(smiles)
    return ADMETResponse(id=mol_id, smiles=smiles, **cached)


def build_solubility(mol_id: str, smiles: str) -> SolubilityResponse:
    cached = _compute_solubility(smiles)
    return SolubilityResponse(id=mol_id, smiles=smiles, **cached)


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
