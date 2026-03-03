"""Team management request/response models."""

from __future__ import annotations

from pydantic import BaseModel, Field


# ------------------------------------------------------------------ #
# Team CRUD
# ------------------------------------------------------------------ #

class TeamCreateRequest(BaseModel):
    name: str = Field(..., min_length=1, max_length=100, description="Team display name")
    slug: str = Field(
        ..., min_length=1, max_length=50,
        pattern=r"^[a-z0-9][a-z0-9\-]*[a-z0-9]$|^[a-z0-9]$",
        description="URL-safe team identifier (lowercase, hyphens)",
    )


class TeamUpdateRequest(BaseModel):
    name: str = Field(..., min_length=1, max_length=100)


class TeamResponse(BaseModel):
    id: int
    name: str
    slug: str
    owner_email: str
    created_at: str
    updated_at: str


class TeamListResponse(BaseModel):
    teams: list[TeamResponse]


# ------------------------------------------------------------------ #
# Member management
# ------------------------------------------------------------------ #

class MemberAddRequest(BaseModel):
    email: str = Field(..., min_length=1, description="Email of user to add")
    role: str = Field(default="member", description="Team role: owner, admin, or member")


class MemberRoleUpdateRequest(BaseModel):
    role: str = Field(..., description="New team role: admin or member")


class MemberResponse(BaseModel):
    id: int
    team_id: int
    user_email: str
    team_role: str
    joined_at: str


class MemberListResponse(BaseModel):
    members: list[MemberResponse]


# ------------------------------------------------------------------ #
# Team library
# ------------------------------------------------------------------ #

class TeamLibrarySaveRequest(BaseModel):
    smiles: str = Field(..., min_length=1, description="SMILES string")
    name: str | None = Field(default=None, description="Optional molecule name")
    tags: list[str] = Field(default_factory=list, description="Tags for categorization")
    notes: str | None = Field(default=None, description="Optional notes")


class TeamLibraryUpdateRequest(BaseModel):
    name: str | None = None
    tags: list[str] | None = None
    notes: str | None = None


class TeamLibraryMoleculeResponse(BaseModel):
    id: int
    team_id: int
    smiles: str
    name: str | None
    tags: list[str]
    notes: str | None
    properties: dict
    added_by: str
    created_at: str
    updated_at: str


class TeamLibraryListResponse(BaseModel):
    molecules: list[TeamLibraryMoleculeResponse]
    total: int
    page: int
    per_page: int


class TeamLibraryImportRequest(BaseModel):
    smiles_list: list[str] = Field(..., min_length=1, description="SMILES strings to import")
    tag: str | None = Field(default=None, description="Optional tag for all imported molecules")


class TeamLibraryImportResponse(BaseModel):
    saved: int
    duplicates: int
    errors: list[str]
