"""Molecule library endpoints."""

from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException

from app.services.database import DatabaseIntegrityError

from app.config import Tier
from app.dependencies import UserContext, check_rate_limit
from app.models.library import (
    LibrarySaveRequest,
    LibraryUpdateRequest,
    LibraryMoleculeResponse,
    LibraryListResponse,
    LibraryImportRequest,
    LibraryImportResponse,
)
from app.services.library_db import get_library_db

router = APIRouter(prefix="/api/v1/library", tags=["library"])

_TIER_LIBRARY_LIMITS = {
    Tier.FREE: 50,
    Tier.PRO: 1000,
    Tier.TEAM: 5000,
    Tier.ACADEMIC: 500,
    Tier.ENTERPRISE: 100000,
}


def _compute_properties(smiles: str) -> dict:
    """Compute and return molecule properties for caching."""
    try:
        from molbuilder.smiles.parser import parse
        from molbuilder.core.elements import atomic_weight
        from molbuilder.io.json_io import _molecular_formula
        from molbuilder.molecule.properties import (
            crippen_logp, heavy_atom_count,
            hydrogen_bond_donors, hydrogen_bond_acceptors,
        )
        from molbuilder.reactions.functional_group_detect import detect_functional_groups

        mol = parse(smiles)
        fgs = detect_functional_groups(mol)
        mw = round(sum(atomic_weight(a.symbol) for a in mol.atoms), 4)
        hac = heavy_atom_count(mol)
        logp = crippen_logp(mol)
        hbd = hydrogen_bond_donors(mol)
        hba = hydrogen_bond_acceptors(mol)
        lipinski_pass = (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10)

        return {
            "molecular_weight": mw,
            "formula": _molecular_formula(mol),
            "logp": round(logp, 4),
            "heavy_atom_count": hac,
            "functional_groups": [fg.name for fg in fgs],
            "hbd": hbd,
            "hba": hba,
            "lipinski_pass": lipinski_pass,
        }
    except Exception:
        return {}


def _check_library_limit(user: UserContext) -> None:
    """Raise 403 if user has reached their tier's library limit."""
    limit = _TIER_LIBRARY_LIMITS.get(user.tier, 50)
    db = get_library_db()
    count = db.count_molecules(user.email)
    if count >= limit:
        raise HTTPException(
            status_code=403,
            detail=f"Library limit of {limit} molecules reached for {user.tier.value} tier",
        )


@router.post("/", response_model=LibraryMoleculeResponse)
def save_molecule(
    body: LibrarySaveRequest,
    user: UserContext = Depends(check_rate_limit),
):
    _check_library_limit(user)

    # Validate SMILES first
    try:
        from molbuilder.smiles.parser import parse
        from molbuilder.smiles.writer import to_smiles
        mol = parse(body.smiles)
        canonical = to_smiles(mol)
    except (ValueError, KeyError) as e:
        raise HTTPException(status_code=400, detail=f"Invalid SMILES: {e}")

    props = _compute_properties(canonical)
    db = get_library_db()

    try:
        mol_id = db.save_molecule(
            user_email=user.email,
            smiles=canonical,
            name=body.name,
            tags=body.tags,
            notes=body.notes,
            properties=props,
        )
    except DatabaseIntegrityError:
        raise HTTPException(status_code=409, detail="Molecule already in library")

    saved = db.get_molecule(mol_id, user.email)
    return LibraryMoleculeResponse(**saved)


@router.get("/", response_model=LibraryListResponse)
def list_molecules(
    page: int = 1,
    per_page: int = 20,
    tag: str | None = None,
    search: str | None = None,
    user: UserContext = Depends(check_rate_limit),
):
    db = get_library_db()
    offset = (page - 1) * per_page
    molecules, total = db.list_molecules(
        user.email, tag=tag, search=search, limit=per_page, offset=offset,
    )
    return LibraryListResponse(
        molecules=[LibraryMoleculeResponse(**m) for m in molecules],
        total=total,
        page=page,
        per_page=per_page,
    )


@router.get("/{mol_id}", response_model=LibraryMoleculeResponse)
def get_molecule(
    mol_id: int,
    user: UserContext = Depends(check_rate_limit),
):
    db = get_library_db()
    mol = db.get_molecule(mol_id, user.email)
    if mol is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return LibraryMoleculeResponse(**mol)


@router.put("/{mol_id}", response_model=LibraryMoleculeResponse)
def update_molecule(
    mol_id: int,
    body: LibraryUpdateRequest,
    user: UserContext = Depends(check_rate_limit),
):
    db = get_library_db()
    updated = db.update_molecule(
        mol_id, user.email,
        name=body.name, tags=body.tags, notes=body.notes,
    )
    if not updated:
        raise HTTPException(status_code=404, detail="Molecule not found")
    mol = db.get_molecule(mol_id, user.email)
    return LibraryMoleculeResponse(**mol)


@router.delete("/{mol_id}")
def delete_molecule(
    mol_id: int,
    user: UserContext = Depends(check_rate_limit),
):
    db = get_library_db()
    deleted = db.delete_molecule(mol_id, user.email)
    if not deleted:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return {"status": "deleted", "id": mol_id}


@router.post("/import", response_model=LibraryImportResponse)
def import_molecules(
    body: LibraryImportRequest,
    user: UserContext = Depends(check_rate_limit),
):
    db = get_library_db()
    saved = 0
    duplicates = 0
    errors = []

    for smiles in body.smiles_list:
        try:
            from molbuilder.smiles.parser import parse
            from molbuilder.smiles.writer import to_smiles
            mol = parse(smiles)
            canonical = to_smiles(mol)
        except (ValueError, KeyError) as e:
            errors.append(f"{smiles}: {e}")
            continue

        # Check library limit
        limit = _TIER_LIBRARY_LIMITS.get(user.tier, 50)
        count = db.count_molecules(user.email)
        if count >= limit:
            errors.append(f"{smiles}: library limit reached")
            continue

        props = _compute_properties(canonical)
        tags = [body.tag] if body.tag else []

        try:
            db.save_molecule(
                user_email=user.email,
                smiles=canonical,
                tags=tags,
                properties=props,
            )
            saved += 1
        except DatabaseIntegrityError:
            duplicates += 1

    return LibraryImportResponse(saved=saved, duplicates=duplicates, errors=errors)
