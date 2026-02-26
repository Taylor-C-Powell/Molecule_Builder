"""Molecule library request/response models."""

from __future__ import annotations

from pydantic import BaseModel, Field


class LibrarySaveRequest(BaseModel):
    smiles: str = Field(..., min_length=1, description="SMILES string")
    name: str | None = Field(default=None, description="Optional molecule name")
    tags: list[str] = Field(default_factory=list, description="Tags for categorization")
    notes: str | None = Field(default=None, description="Optional notes")


class LibraryUpdateRequest(BaseModel):
    name: str | None = None
    tags: list[str] | None = None
    notes: str | None = None


class LibraryMoleculeResponse(BaseModel):
    id: int
    smiles: str
    name: str | None
    tags: list[str]
    notes: str | None
    properties: dict
    created_at: str
    updated_at: str


class LibraryListResponse(BaseModel):
    molecules: list[LibraryMoleculeResponse]
    total: int
    page: int
    per_page: int


class LibraryImportRequest(BaseModel):
    smiles_list: list[str] = Field(..., min_length=1, description="SMILES strings to import")
    tag: str | None = Field(default=None, description="Optional tag for all imported molecules")


class LibraryImportResponse(BaseModel):
    saved: int
    duplicates: int
    errors: list[str]
