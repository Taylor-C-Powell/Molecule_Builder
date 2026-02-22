"""Pydantic models for feasibility scoring API."""

from __future__ import annotations

from pydantic import BaseModel, Field


class FeasibilityRequest(BaseModel):
    """Request for feasibility scoring."""

    smiles: str = Field(..., description="SMILES string of target molecule")
    max_depth: int = Field(5, ge=1, le=10, description="Max retrosynthesis depth")
    beam_width: int = Field(5, ge=1, le=10, description="Beam search width")


class FeasibilityCompareRequest(BaseModel):
    """Request for comparing feasibility of multiple molecules."""

    smiles_list: list[str] = Field(
        ..., min_length=1, max_length=10,
        description="List of SMILES strings to compare",
    )
    max_depth: int = Field(5, ge=1, le=10)
    beam_width: int = Field(5, ge=1, le=10)


class DimensionScores(BaseModel):
    """Per-dimension feasibility scores (0-100)."""

    availability: float
    cost: float
    green_chemistry: float
    safety: float
    complexity: float
    regulatory: float


class RouteSummary(BaseModel):
    """Summary of the best synthesis route."""

    total_steps: int
    overall_yield: float
    starting_materials: list[str]
    routes_found: int


class FeasibilityResponse(BaseModel):
    """Response with feasibility assessment."""

    target_smiles: str
    composite_score: float = Field(..., description="Overall feasibility score (0-100)")
    grade: str = Field(..., description="Letter grade: A, B, C, D, F")
    is_feasible: bool
    dimensions: DimensionScores
    route_summary: RouteSummary


class FeasibilityCompareResponse(BaseModel):
    """Response comparing multiple molecules."""

    results: list[FeasibilityResponse]
    best_target: str
    best_score: float
