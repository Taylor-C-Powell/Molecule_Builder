"""Process engineering request/response models."""

from pydantic import BaseModel, Field


class ProcessEvaluateRequest(BaseModel):
    smiles: str = Field(..., description="Target SMILES")
    scale_kg: float = Field(1.0, gt=0, description="Production scale in kg")
    max_depth: int = Field(5, ge=1, le=10)
    beam_width: int = Field(5, ge=1, le=10)


class ReactorResponse(BaseModel):
    reactor_type: str
    volume_L: float
    temperature_C: float
    pressure_atm: float
    residence_time_min: float
    mixing_type: str
    heat_transfer: str
    material: str
    estimated_cost_usd: float
    notes: str


class ConditionsResponse(BaseModel):
    temperature_C: float
    pressure_atm: float
    solvent: str
    concentration_M: float
    addition_rate: str
    reaction_time_hours: float
    atmosphere: str
    workup_procedure: str
    notes: str
    data_source: str = "heuristic"


class PurificationStepResponse(BaseModel):
    method: str
    description: str
    estimated_recovery: float
    estimated_purity: float
    scale_appropriate: bool
    notes: str


class HazardInfoResponse(BaseModel):
    reagent_name: str
    ghs_hazards: list[str]
    ghs_pictograms: list[str]
    hazard_descriptions: list[str]
    pictogram_descriptions: list[str]


class SafetyAssessmentResponse(BaseModel):
    step_number: int
    step_name: str
    risk_level: str
    ppe_required: list[str]
    engineering_controls: list[str]
    emergency_procedures: list[str]
    incompatible_materials: list[str]
    waste_classification: str
    hazards: list[HazardInfoResponse]


class CostBreakdownResponse(BaseModel):
    raw_materials_usd: float
    labor_usd: float
    equipment_usd: float
    energy_usd: float
    waste_disposal_usd: float
    overhead_usd: float


class CostEstimateResponse(BaseModel):
    total_usd: float
    per_kg_usd: float
    scale_kg: float
    breakdown: CostBreakdownResponse
    notes: list[str]


class ScaleUpResponse(BaseModel):
    target_annual_kg: float
    recommended_mode: str
    batch_size_kg: float | None
    batches_per_year: int | None
    cycle_time_hours: float
    annual_capacity_kg: float
    capital_cost_usd: float
    operating_cost_annual_usd: float
    scale_up_risks: list[str]
    recommendations: list[str]


class StepProcessDetail(BaseModel):
    step_number: int
    reaction_name: str
    reactor: ReactorResponse
    conditions: ConditionsResponse
    purification: list[PurificationStepResponse]


class ProcessEvaluateResponse(BaseModel):
    smiles: str
    scale_kg: float
    route_found: bool
    total_steps: int = 0
    overall_yield: float = 0.0
    step_details: list[StepProcessDetail] = []
    safety: list[SafetyAssessmentResponse] = []
    cost: CostEstimateResponse | None = None
    scale_up: ScaleUpResponse | None = None


# -----------------------------------------------------------------
#  Condition prediction models
# -----------------------------------------------------------------

class PredictConditionsRequest(BaseModel):
    smiles: str = Field(..., description="Substrate SMILES")
    reaction_name: str | None = Field(None, description="Optional reaction type hint")
    scale_kg: float = Field(1.0, gt=0, description="Production scale in kg")
    max_candidates: int = Field(5, ge=1, le=20, description="Max ranked results")


class SolventScoreResponse(BaseModel):
    solvent_name: str
    composite_score: float
    green_score: int
    reasons: list[str]


class SubstrateAnalysisResponse(BaseModel):
    detected_functional_groups: list[str]
    molecular_weight: float
    heavy_atom_count: int
    steric_class: str
    electronic_character: str
    sensitive_groups: list[str]


class TemplateMatchResponse(BaseModel):
    template_name: str
    named_reaction: str | None
    category: str
    match_score: float
    match_reasons: list[str]
    conditions: ConditionsResponse
    recommended_solvents: list[SolventScoreResponse]
    adjusted_yield_range: list[float]
    warnings: list[str]
    data_source: str = "heuristic"


class PredictConditionsResponse(BaseModel):
    smiles: str
    substrate_analysis: SubstrateAnalysisResponse
    candidates: list[TemplateMatchResponse]
    best_match: TemplateMatchResponse | None
    overall_confidence: str
