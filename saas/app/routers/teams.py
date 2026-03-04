"""Team management endpoints -- CRUD, members, and shared library."""

from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException, Query

from app.auth.roles import Role
from app.auth.team_roles import TeamRole, can_manage_members, can_delete_team
from app.config import Tier
from app.dependencies import UserContext, check_rate_limit
from app.models.team import (
    TeamCreateRequest,
    TeamUpdateRequest,
    TeamResponse,
    TeamListResponse,
    MemberAddRequest,
    MemberRoleUpdateRequest,
    MemberResponse,
    MemberListResponse,
    TeamLibrarySaveRequest,
    TeamLibraryUpdateRequest,
    TeamLibraryMoleculeResponse,
    TeamLibraryListResponse,
    TeamLibraryImportRequest,
    TeamLibraryImportResponse,
)
from app.services.database import DatabaseIntegrityError
from app.services.team_db import get_team_db

router = APIRouter(prefix="/api/v1/teams", tags=["teams"])

_TEAM_TIERS = {Tier.TEAM, Tier.ENTERPRISE}


# ------------------------------------------------------------------ #
# Authorization helpers
# ------------------------------------------------------------------ #

def _require_membership(
    team_id: int, user: UserContext,
) -> dict:
    """Return the membership record or raise 403/404.

    System ADMINs get virtual owner access without a real membership row.
    """
    db = get_team_db()
    team = db.get_team(team_id)
    if team is None:
        raise HTTPException(status_code=404, detail="Team not found")

    # System admin gets virtual owner access
    if user.role == Role.ADMIN:
        return {
            "team_id": team_id,
            "user_email": user.email,
            "team_role": TeamRole.OWNER.value,
            "joined_at": "",
        }

    member = db.get_member(team_id, user.email)
    if member is None:
        raise HTTPException(status_code=403, detail="Not a member of this team")
    return member


def _require_manager(
    team_id: int, user: UserContext,
) -> dict:
    """Return membership record if user is owner or admin, else 403."""
    member = _require_membership(team_id, user)
    role = TeamRole(member["team_role"])
    if not can_manage_members(role):
        raise HTTPException(
            status_code=403, detail="Team admin or owner access required",
        )
    return member


# ------------------------------------------------------------------ #
# Team CRUD
# ------------------------------------------------------------------ #

@router.post("/", response_model=TeamResponse, status_code=201)
def create_team(
    body: TeamCreateRequest,
    user: UserContext = Depends(check_rate_limit),
):
    if user.tier not in _TEAM_TIERS and user.role != Role.ADMIN:
        raise HTTPException(
            status_code=403,
            detail="Team creation requires team or enterprise tier",
        )
    db = get_team_db()
    try:
        team_id = db.create_team(body.name, body.slug, user.email)
    except DatabaseIntegrityError:
        raise HTTPException(status_code=409, detail="Team slug already exists")
    team = db.get_team(team_id)
    return TeamResponse(**team)


@router.get("/", response_model=TeamListResponse)
def list_teams(user: UserContext = Depends(check_rate_limit)):
    db = get_team_db()
    teams = db.list_teams_for_user(user.email)
    return TeamListResponse(teams=[TeamResponse(**t) for t in teams])


@router.get("/{team_id}", response_model=TeamResponse)
def get_team(team_id: int, user: UserContext = Depends(check_rate_limit)):
    _require_membership(team_id, user)
    db = get_team_db()
    team = db.get_team(team_id)
    return TeamResponse(**team)


@router.patch("/{team_id}", response_model=TeamResponse)
def update_team(
    team_id: int,
    body: TeamUpdateRequest,
    user: UserContext = Depends(check_rate_limit),
):
    _require_manager(team_id, user)
    db = get_team_db()
    db.update_team(team_id, body.name)
    team = db.get_team(team_id)
    return TeamResponse(**team)


@router.delete("/{team_id}")
def delete_team(team_id: int, user: UserContext = Depends(check_rate_limit)):
    member = _require_membership(team_id, user)
    role = TeamRole(member["team_role"])
    if not can_delete_team(role):
        raise HTTPException(status_code=403, detail="Only the team owner can delete the team")
    db = get_team_db()
    db.delete_team(team_id)
    return {"status": "deleted", "id": team_id}


# ------------------------------------------------------------------ #
# Member management
# ------------------------------------------------------------------ #

@router.get("/{team_id}/members", response_model=MemberListResponse)
def list_members(team_id: int, user: UserContext = Depends(check_rate_limit)):
    _require_membership(team_id, user)
    db = get_team_db()
    members = db.list_members(team_id)
    return MemberListResponse(members=[MemberResponse(**m) for m in members])


@router.post("/{team_id}/members", response_model=MemberResponse, status_code=201)
def add_member(
    team_id: int,
    body: MemberAddRequest,
    user: UserContext = Depends(check_rate_limit),
):
    _require_manager(team_id, user)

    # Validate requested role
    try:
        requested_role = TeamRole(body.role)
    except ValueError:
        raise HTTPException(status_code=400, detail=f"Invalid role: {body.role}")

    if requested_role == TeamRole.OWNER:
        raise HTTPException(status_code=400, detail="Cannot add a member as owner")

    db = get_team_db()
    try:
        member_id = db.add_member(team_id, body.email, requested_role.value)
    except DatabaseIntegrityError:
        raise HTTPException(status_code=409, detail="User is already a member")

    member = db.get_member(team_id, body.email)
    return MemberResponse(**member)


@router.patch("/{team_id}/members/{email}", response_model=MemberResponse)
def update_member_role(
    team_id: int,
    email: str,
    body: MemberRoleUpdateRequest,
    user: UserContext = Depends(check_rate_limit),
):
    _require_manager(team_id, user)

    # Validate the new role
    try:
        new_role = TeamRole(body.role)
    except ValueError:
        raise HTTPException(status_code=400, detail=f"Invalid role: {body.role}")

    if new_role == TeamRole.OWNER:
        raise HTTPException(status_code=400, detail="Cannot promote to owner")

    db = get_team_db()
    target = db.get_member(team_id, email)
    if target is None:
        raise HTTPException(status_code=404, detail="Member not found")
    if target["team_role"] == TeamRole.OWNER.value:
        raise HTTPException(status_code=400, detail="Cannot change the owner's role")

    db.update_member_role(team_id, email, new_role.value)
    member = db.get_member(team_id, email)
    return MemberResponse(**member)


@router.delete("/{team_id}/members/{email}")
def remove_member(
    team_id: int,
    email: str,
    user: UserContext = Depends(check_rate_limit),
):
    db = get_team_db()
    target = db.get_member(team_id, email)
    if target is None:
        raise HTTPException(status_code=404, detail="Member not found")
    if target["team_role"] == TeamRole.OWNER.value:
        raise HTTPException(status_code=400, detail="Owner cannot be removed")

    # Self-removal: any member can leave
    is_self = (email == user.email)
    if not is_self:
        _require_manager(team_id, user)

    db.remove_member(team_id, email)
    return {"status": "removed", "email": email}


# ------------------------------------------------------------------ #
# Team library
# ------------------------------------------------------------------ #

def _compute_properties(smiles: str) -> dict:
    """Compute molecule properties (reuses logic from library router)."""
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


@router.post("/{team_id}/library/", response_model=TeamLibraryMoleculeResponse)
def save_team_molecule(
    team_id: int,
    body: TeamLibrarySaveRequest,
    user: UserContext = Depends(check_rate_limit),
):
    _require_membership(team_id, user)

    try:
        from molbuilder.smiles.parser import parse
        from molbuilder.smiles.writer import to_smiles
        mol = parse(body.smiles)
        canonical = to_smiles(mol)
    except (ValueError, KeyError) as e:
        raise HTTPException(status_code=400, detail=f"Invalid SMILES: {e}")

    props = _compute_properties(canonical)
    db = get_team_db()
    try:
        mol_id = db.save_molecule(
            team_id=team_id,
            smiles=canonical,
            added_by=user.email,
            name=body.name,
            tags=body.tags,
            notes=body.notes,
            properties=props,
        )
    except DatabaseIntegrityError:
        raise HTTPException(status_code=409, detail="Molecule already in team library")

    saved = db.get_molecule(mol_id, team_id)
    return TeamLibraryMoleculeResponse(**saved)


@router.get("/{team_id}/library/", response_model=TeamLibraryListResponse)
def list_team_molecules(
    team_id: int,
    page: int = Query(default=1, ge=1),
    per_page: int = Query(default=20, ge=1, le=100),
    tag: str | None = None,
    search: str | None = None,
    user: UserContext = Depends(check_rate_limit),
):
    _require_membership(team_id, user)
    db = get_team_db()
    offset = (page - 1) * per_page
    molecules, total = db.list_molecules(
        team_id, tag=tag, search=search, limit=per_page, offset=offset,
    )
    return TeamLibraryListResponse(
        molecules=[TeamLibraryMoleculeResponse(**m) for m in molecules],
        total=total,
        page=page,
        per_page=per_page,
    )


@router.get("/{team_id}/library/{mol_id}", response_model=TeamLibraryMoleculeResponse)
def get_team_molecule(
    team_id: int,
    mol_id: int,
    user: UserContext = Depends(check_rate_limit),
):
    _require_membership(team_id, user)
    db = get_team_db()
    mol = db.get_molecule(mol_id, team_id)
    if mol is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return TeamLibraryMoleculeResponse(**mol)


@router.patch("/{team_id}/library/{mol_id}", response_model=TeamLibraryMoleculeResponse)
def update_team_molecule(
    team_id: int,
    mol_id: int,
    body: TeamLibraryUpdateRequest,
    user: UserContext = Depends(check_rate_limit),
):
    _require_membership(team_id, user)
    db = get_team_db()
    updated = db.update_molecule(
        mol_id, team_id,
        name=body.name, tags=body.tags, notes=body.notes,
    )
    if not updated:
        raise HTTPException(status_code=404, detail="Molecule not found")
    mol = db.get_molecule(mol_id, team_id)
    return TeamLibraryMoleculeResponse(**mol)


@router.delete("/{team_id}/library/{mol_id}")
def delete_team_molecule(
    team_id: int,
    mol_id: int,
    user: UserContext = Depends(check_rate_limit),
):
    _require_membership(team_id, user)
    db = get_team_db()
    deleted = db.delete_molecule(mol_id, team_id)
    if not deleted:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return {"status": "deleted", "id": mol_id}


@router.post("/{team_id}/library/import", response_model=TeamLibraryImportResponse)
def import_team_molecules(
    team_id: int,
    body: TeamLibraryImportRequest,
    user: UserContext = Depends(check_rate_limit),
):
    _require_membership(team_id, user)
    db = get_team_db()
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

        props = _compute_properties(canonical)
        tags = [body.tag] if body.tag else []

        try:
            db.save_molecule(
                team_id=team_id,
                smiles=canonical,
                added_by=user.email,
                tags=tags,
                properties=props,
            )
            saved += 1
        except DatabaseIntegrityError:
            duplicates += 1

    return TeamLibraryImportResponse(saved=saved, duplicates=duplicates, errors=errors)
