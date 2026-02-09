"""Retrosynthesis request/response models."""

from pydantic import BaseModel, Field


class RetroRequest(BaseModel):
    smiles: str = Field(..., description="Target SMILES to retro-analyze")
    max_depth: int = Field(5, ge=1, le=10, description="Max retrosynthetic steps")
    beam_width: int = Field(5, ge=1, le=10, description="Beam search width")


class PrecursorResponse(BaseModel):
    smiles: str
    name: str
    cost_per_kg: float


class DisconnectionResponse(BaseModel):
    reaction_name: str
    named_reaction: str | None
    category: str
    score: float
    precursors: list[PrecursorResponse]


class RetroNodeResponse(BaseModel):
    smiles: str
    is_purchasable: bool
    functional_groups: list[str]
    best_disconnection: DisconnectionResponse | None = None
    children: list["RetroNodeResponse"] = []


class SynthesisStepResponse(BaseModel):
    step_number: int
    reaction_name: str
    named_reaction: str | None
    category: str
    precursor_smiles: list[str]
    product_smiles: str
    product_name: str
    conditions: str
    expected_yield: float
    notes: str


class SynthesisRouteResponse(BaseModel):
    target_smiles: str
    target_name: str
    total_steps: int
    overall_yield: float
    longest_linear_sequence: int
    starting_materials: list[PrecursorResponse]
    steps: list[SynthesisStepResponse]


class RetroResponse(BaseModel):
    tree: RetroNodeResponse
    routes_found: int
    max_depth: int
    beam_width: int
    best_route: SynthesisRouteResponse | None = None
