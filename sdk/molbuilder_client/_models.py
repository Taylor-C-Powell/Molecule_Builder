"""Frozen dataclass models for all API response types.

All models are immutable (frozen=True) and use slots for memory efficiency.
No external dependencies -- stdlib only.
"""

from __future__ import annotations

from dataclasses import dataclass, field


# -- Auth ---------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class APIKeyInfo:
    api_key: str
    email: str
    tier: str
    role: str


@dataclass(frozen=True, slots=True)
class Token:
    access_token: str
    token_type: str
    expires_in: int


# -- Molecule -----------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class MoleculeInfo:
    id: str
    name: str
    smiles: str
    num_atoms: int
    num_bonds: int


@dataclass(frozen=True, slots=True)
class MoleculeProperties:
    id: str
    smiles: str
    formula: str
    molecular_weight: float
    num_atoms: int
    num_bonds: int
    functional_groups: list[str] = field(default_factory=list)


@dataclass(frozen=True, slots=True)
class Atom3D:
    index: int
    symbol: str
    position: tuple[float, float, float]
    hybridization: str | None = None
    formal_charge: int = 0


@dataclass(frozen=True, slots=True)
class Bond3D:
    atom_i: int
    atom_j: int
    order: float
    rotatable: bool


@dataclass(frozen=True, slots=True)
class Molecule3D:
    id: str
    atoms: list[Atom3D] = field(default_factory=list)
    bonds: list[Bond3D] = field(default_factory=list)


# -- Elements -----------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class Element:
    atomic_number: int
    symbol: str
    name: str
    atomic_weight: float


# -- Retrosynthesis -----------------------------------------------------------

@dataclass(frozen=True, slots=True)
class Precursor:
    smiles: str
    name: str
    cost_per_kg: float


@dataclass(frozen=True, slots=True)
class Disconnection:
    reaction_name: str
    named_reaction: str | None
    category: str
    score: float
    precursors: list[Precursor] = field(default_factory=list)


@dataclass(frozen=True, slots=True)
class RetroNode:
    smiles: str
    is_purchasable: bool
    functional_groups: list[str] = field(default_factory=list)
    best_disconnection: Disconnection | None = None
    children: list[RetroNode] = field(default_factory=list)


@dataclass(frozen=True, slots=True)
class RouteStep:
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


@dataclass(frozen=True, slots=True)
class BestRoute:
    target_smiles: str
    target_name: str
    total_steps: int
    overall_yield: float
    longest_linear_sequence: int
    starting_materials: list[Precursor] = field(default_factory=list)
    steps: list[RouteStep] = field(default_factory=list)


@dataclass(frozen=True, slots=True)
class RetrosynthesisPlan:
    tree: RetroNode
    routes_found: int
    max_depth: int
    beam_width: int
    best_route: BestRoute | None = None


# -- Process Evaluation -------------------------------------------------------

@dataclass(frozen=True, slots=True)
class Reactor:
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


@dataclass(frozen=True, slots=True)
class Conditions:
    temperature_C: float
    pressure_atm: float
    solvent: str
    concentration_M: float
    addition_rate: str
    reaction_time_hours: float
    atmosphere: str
    workup_procedure: str
    notes: str


@dataclass(frozen=True, slots=True)
class PurificationMethod:
    method: str
    description: str
    estimated_recovery: float
    estimated_purity: float
    scale_appropriate: bool
    notes: str


@dataclass(frozen=True, slots=True)
class StepDetail:
    step_number: int
    reaction_name: str
    reactor: Reactor
    conditions: Conditions
    purification: list[PurificationMethod] = field(default_factory=list)


@dataclass(frozen=True, slots=True)
class Hazard:
    reagent_name: str
    ghs_hazards: list[str] = field(default_factory=list)
    ghs_pictograms: list[str] = field(default_factory=list)
    hazard_descriptions: list[str] = field(default_factory=list)
    pictogram_descriptions: list[str] = field(default_factory=list)


@dataclass(frozen=True, slots=True)
class SafetyAssessment:
    step_number: int
    step_name: str
    risk_level: str
    ppe_required: list[str] = field(default_factory=list)
    engineering_controls: list[str] = field(default_factory=list)
    emergency_procedures: list[str] = field(default_factory=list)
    incompatible_materials: list[str] = field(default_factory=list)
    waste_classification: str = ""
    hazards: list[Hazard] = field(default_factory=list)


@dataclass(frozen=True, slots=True)
class CostBreakdown:
    raw_materials_usd: float
    labor_usd: float
    equipment_usd: float
    energy_usd: float
    waste_disposal_usd: float
    overhead_usd: float


@dataclass(frozen=True, slots=True)
class CostEstimate:
    total_usd: float
    per_kg_usd: float
    scale_kg: float
    breakdown: CostBreakdown
    notes: list[str] = field(default_factory=list)


@dataclass(frozen=True, slots=True)
class ScaleUp:
    target_annual_kg: float
    recommended_mode: str
    batch_size_kg: float | None
    batches_per_year: float | None
    cycle_time_hours: float
    annual_capacity_kg: float
    capital_cost_usd: float
    operating_cost_annual_usd: float
    scale_up_risks: list[str] = field(default_factory=list)
    recommendations: list[str] = field(default_factory=list)


@dataclass(frozen=True, slots=True)
class ProcessEvaluation:
    smiles: str
    scale_kg: float
    route_found: bool
    total_steps: int = 0
    overall_yield: float = 0.0
    step_details: list[StepDetail] = field(default_factory=list)
    safety: list[SafetyAssessment] = field(default_factory=list)
    cost: CostEstimate | None = None
    scale_up: ScaleUp | None = None


# -- Billing ------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class CheckoutSession:
    checkout_url: str


@dataclass(frozen=True, slots=True)
class BillingStatus:
    email: str
    tier: str
    subscription_status: str
    stripe_customer_id: str | None = None
    stripe_subscription_id: str | None = None
