"""Retrosynthesis Feasibility Scoring Engine.

Extends the retrosynthesis engine with comprehensive feasibility assessment
combining synthetic accessibility, precursor availability, cost estimation,
environmental impact, and regulatory classification into a unified score.

Patent-relevant novelty: Multi-dimensional feasibility scoring that integrates
retrosynthetic route analysis with real-world constraints (cost, availability,
regulation, green chemistry) into a single actionable composite score with
per-dimension breakdowns.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

import numpy as np

from molbuilder.reactions.retrosynthesis import (
    RetrosynthesisTree,
    RetroNode,
    Disconnection,
    Precursor,
    retrosynthesis,
    PURCHASABLE_MATERIALS,
    is_purchasable,
)
from molbuilder.reactions.synthesis_route import (
    SynthesisRoute,
    SynthesisStep,
    extract_best_route,
)
from molbuilder.reactions.reaction_types import ReactionCategory
from molbuilder.reactions.reagent_data import REAGENT_DB, SOLVENT_DB


class RegulatoryClass(Enum):
    """Regulatory classification for synthetic targets."""

    UNCLASSIFIED = "unclassified"
    GRAS = "gras"  # Generally Recognized as Safe
    OTC = "otc"  # Over-the-counter drug
    PRESCRIPTION = "prescription"
    CONTROLLED = "controlled"
    RESTRICTED = "restricted"


class HazardLevel(Enum):
    """Safety hazard level for synthetic routes."""

    LOW = "low"
    MODERATE = "moderate"
    HIGH = "high"
    EXTREME = "extreme"


# Green chemistry solvent scores (higher = greener)
GREEN_SOLVENT_SCORES = {
    "water": 1.0,
    "ethanol": 0.9,
    "isopropanol": 0.85,
    "ethyl acetate": 0.8,
    "acetone": 0.75,
    "methanol": 0.7,
    "toluene": 0.4,
    "dichloromethane": 0.2,
    "chloroform": 0.15,
    "dmf": 0.3,
    "dmso": 0.5,
    "thf": 0.35,
    "diethyl ether": 0.3,
    "hexane": 0.25,
    "benzene": 0.1,
    "carbon tetrachloride": 0.05,
}

# Known controlled/restricted substances (simplified - SMILES patterns)
CONTROLLED_PATTERNS = {
    "ephedrine": "C[C@@H](O)[C@H](NC)c1ccccc1",
    "pseudoephedrine": "C[C@@H](O)[C@@H](NC)c1ccccc1",
    "ergotamine": None,  # Complex SMILES
    "lysergic_acid": None,
}


@dataclass
class AvailabilityScore:
    """Precursor availability assessment."""

    score: float  # 0.0-1.0
    n_purchasable: int
    n_total_precursors: int
    fraction_purchasable: float
    estimated_lead_time_days: int
    supply_chain_risk: str  # low, moderate, high
    details: list[dict[str, Any]] = field(default_factory=list)


@dataclass
class CostEstimate:
    """Cost estimate for a synthetic route."""

    score: float  # 0.0-1.0 (inverse of cost, normalized)
    total_material_cost_usd: float
    cost_per_gram_target: float
    reagent_costs: dict[str, float] = field(default_factory=dict)
    precursor_costs: dict[str, float] = field(default_factory=dict)
    labor_cost_estimate: float = 0.0
    scale: str = "lab"  # lab, pilot, production


@dataclass
class GreenChemistryScore:
    """Environmental impact assessment."""

    score: float  # 0.0-1.0 (higher = greener)
    atom_economy: float
    solvent_greenness: float
    step_economy: float  # fewer steps = greener
    hazardous_reagent_fraction: float
    e_factor_estimate: float  # kg waste / kg product
    details: str = ""


@dataclass
class SafetyAssessment:
    """Safety assessment for the synthetic route."""

    score: float  # 0.0-1.0 (higher = safer)
    hazard_level: HazardLevel
    hazardous_steps: list[int]  # Step indices with safety concerns
    ghs_hazards_encountered: list[str]
    requires_special_equipment: bool
    details: str = ""


@dataclass
class SyntheticComplexityScore:
    """Synthetic complexity assessment."""

    score: float  # 0.0-1.0 (higher = more feasible / less complex)
    total_steps: int
    longest_linear_sequence: int
    overall_yield: float
    convergent: bool  # Is the route convergent?
    max_step_difficulty: float
    details: str = ""


@dataclass
class RegulatoryAssessment:
    """Regulatory classification assessment."""

    score: float  # 0.0-1.0 (higher = fewer regulatory concerns)
    classification: RegulatoryClass
    controlled_precursors: list[str]
    requires_dea_license: bool
    export_controlled: bool
    details: str = ""


@dataclass
class FeasibilityResult:
    """Complete feasibility assessment for a retrosynthetic route."""

    target_smiles: str
    composite_score: float  # 0.0-100.0
    grade: str  # A, B, C, D, F
    availability: AvailabilityScore
    cost: CostEstimate
    green_chemistry: GreenChemistryScore
    safety: SafetyAssessment
    complexity: SyntheticComplexityScore
    regulatory: RegulatoryAssessment
    route_summary: dict[str, Any] = field(default_factory=dict)

    @property
    def is_feasible(self) -> bool:
        return self.composite_score >= 40.0

    def to_dict(self) -> dict[str, Any]:
        return {
            "target_smiles": self.target_smiles,
            "composite_score": round(self.composite_score, 2),
            "grade": self.grade,
            "is_feasible": self.is_feasible,
            "dimensions": {
                "availability": round(self.availability.score * 100, 1),
                "cost": round(self.cost.score * 100, 1),
                "green_chemistry": round(self.green_chemistry.score * 100, 1),
                "safety": round(self.safety.score * 100, 1),
                "complexity": round(self.complexity.score * 100, 1),
                "regulatory": round(self.regulatory.score * 100, 1),
            },
            "route_summary": self.route_summary,
        }


class FeasibilityEngine:
    """Multi-dimensional retrosynthesis feasibility scoring engine.

    Combines retrosynthetic route analysis with real-world feasibility
    constraints to produce a composite score.

    Scoring dimensions and weights:
    - Precursor availability (25%): Are starting materials obtainable?
    - Cost (20%): Total route cost per gram of target
    - Synthetic complexity (20%): Steps, yield, convergence
    - Green chemistry (15%): Atom economy, solvent impact, waste
    - Safety (10%): Hazardous reagents and conditions
    - Regulatory (10%): Controlled substances, export restrictions
    """

    DIMENSION_WEIGHTS = {
        "availability": 0.25,
        "cost": 0.20,
        "complexity": 0.20,
        "green_chemistry": 0.15,
        "safety": 0.10,
        "regulatory": 0.10,
    }

    def __init__(
        self,
        max_depth: int = 5,
        beam_width: int = 5,
        cost_threshold_per_gram: float = 10000.0,
        custom_weights: dict[str, float] | None = None,
    ) -> None:
        self.max_depth = max_depth
        self.beam_width = beam_width
        self.cost_threshold = cost_threshold_per_gram
        if custom_weights:
            self.DIMENSION_WEIGHTS = {**self.DIMENSION_WEIGHTS, **custom_weights}

    def score_availability(self, route: SynthesisRoute) -> AvailabilityScore:
        """Assess precursor availability."""
        if not route.starting_materials:
            return AvailabilityScore(
                score=0.0, n_purchasable=0, n_total_precursors=0,
                fraction_purchasable=0.0, estimated_lead_time_days=365,
                supply_chain_risk="high",
            )

        n_total = len(route.starting_materials)
        purchasable_count = 0
        details = []

        for sm in route.starting_materials:
            smiles = sm.smiles if hasattr(sm, "smiles") else str(sm)
            available = is_purchasable(smiles)
            purchasable_count += int(available)
            details.append({
                "smiles": smiles,
                "purchasable": available,
                "name": sm.name if hasattr(sm, "name") else "unknown",
            })

        fraction = purchasable_count / max(n_total, 1)

        # Lead time estimate
        if fraction >= 1.0:
            lead_time = 3
            risk = "low"
        elif fraction >= 0.7:
            lead_time = 14
            risk = "low"
        elif fraction >= 0.4:
            lead_time = 30
            risk = "moderate"
        else:
            lead_time = 90
            risk = "high"

        # Score: heavily penalize if key precursors unavailable
        score = fraction ** 0.7  # Concave: partial availability still gets decent score

        return AvailabilityScore(
            score=score,
            n_purchasable=purchasable_count,
            n_total_precursors=n_total,
            fraction_purchasable=fraction,
            estimated_lead_time_days=lead_time,
            supply_chain_risk=risk,
            details=details,
        )

    def score_cost(self, route: SynthesisRoute) -> CostEstimate:
        """Estimate route cost."""
        precursor_costs: dict[str, float] = {}
        total_material = 0.0

        for sm in route.starting_materials:
            smiles = sm.smiles if hasattr(sm, "smiles") else str(sm)
            cost = sm.cost_per_kg if hasattr(sm, "cost_per_kg") else 100.0
            precursor_costs[smiles] = cost
            total_material += cost * 0.001  # Assume 1g scale -> 0.001 kg

        # Reagent costs from steps
        reagent_costs: dict[str, float] = {}
        for step in route.steps:
            conditions = step.conditions if hasattr(step, "conditions") else ""
            if isinstance(conditions, str) and conditions:
                # Rough estimate per reagent
                reagent_costs[conditions] = reagent_costs.get(conditions, 0.0) + 10.0
                total_material += 10.0

        # Labor cost estimate (simple: $50/step at lab scale)
        labor = len(route.steps) * 50.0

        # Cost per gram of target (accounting for yield)
        overall_yield = max(route.overall_yield, 0.01)
        cost_per_gram = (total_material + labor) / overall_yield

        # Score: normalized inverse of cost
        score = max(0.0, 1.0 - cost_per_gram / self.cost_threshold)
        score = score ** 0.5  # Soften the penalty curve

        return CostEstimate(
            score=score,
            total_material_cost_usd=total_material,
            cost_per_gram_target=cost_per_gram,
            reagent_costs=reagent_costs,
            precursor_costs=precursor_costs,
            labor_cost_estimate=labor,
        )

    def score_green_chemistry(self, route: SynthesisRoute) -> GreenChemistryScore:
        """Assess environmental impact using green chemistry principles."""
        # Step economy: fewer steps = greener
        n_steps = max(len(route.steps), 1)
        step_economy = max(0.0, 1.0 - (n_steps - 1) * 0.1)  # Decreases with steps

        # Atom economy: overall yield as proxy
        atom_economy = min(1.0, route.overall_yield * 1.5)

        # Solvent greenness
        solvent_scores = []
        for step in route.steps:
            conditions = step.conditions if hasattr(step, "conditions") else ""
            if isinstance(conditions, str):
                for solvent, green_score in GREEN_SOLVENT_SCORES.items():
                    if solvent in conditions.lower():
                        solvent_scores.append(green_score)
                        break
                else:
                    solvent_scores.append(0.5)  # Unknown solvent

        solvent_greenness = float(np.mean(solvent_scores)) if solvent_scores else 0.5

        # Hazardous reagent fraction (simplified)
        hazardous_count = 0
        for step in route.steps:
            conditions = step.conditions if hasattr(step, "conditions") else ""
            if isinstance(conditions, str):
                hazardous_keywords = ["lialh4", "nah", "buli", "hf", "hcn", "osmium", "cr"]
                if any(kw in conditions.lower() for kw in hazardous_keywords):
                    hazardous_count += 1

        hazardous_fraction = hazardous_count / max(n_steps, 1)

        # E-factor estimate (kg waste / kg product)
        e_factor = (1.0 / max(route.overall_yield, 0.01) - 1.0) * n_steps
        e_factor = min(100.0, e_factor)

        # Composite green score
        green_score = (
            step_economy * 0.25
            + atom_economy * 0.30
            + solvent_greenness * 0.25
            + (1.0 - hazardous_fraction) * 0.20
        )

        return GreenChemistryScore(
            score=green_score,
            atom_economy=atom_economy,
            solvent_greenness=solvent_greenness,
            step_economy=step_economy,
            hazardous_reagent_fraction=hazardous_fraction,
            e_factor_estimate=e_factor,
        )

    def score_safety(self, route: SynthesisRoute) -> SafetyAssessment:
        """Assess safety of the synthetic route."""
        hazardous_steps: list[int] = []
        ghs_hazards: list[str] = []
        requires_special = False

        for i, step in enumerate(route.steps):
            conditions = step.conditions if hasattr(step, "conditions") else ""
            if isinstance(conditions, str):
                # Check for high-hazard conditions
                if any(kw in conditions.lower() for kw in [
                    "lialh4", "buli", "grignard", "sodium hydride",
                ]):
                    hazardous_steps.append(i)
                    ghs_hazards.append("GHS02-Flammable")

                if any(kw in conditions.lower() for kw in [
                    "hcn", "phosgene", "osmium", "azide",
                ]):
                    hazardous_steps.append(i)
                    ghs_hazards.append("GHS06-Toxic")
                    requires_special = True

                if any(kw in conditions.lower() for kw in [
                    "pressure", "autoclave", "hydrogenation",
                ]):
                    hazardous_steps.append(i)
                    ghs_hazards.append("GHS04-CompressedGas")
                    requires_special = True

        n_steps = max(len(route.steps), 1)
        hazard_fraction = len(set(hazardous_steps)) / n_steps

        if hazard_fraction == 0:
            hazard_level = HazardLevel.LOW
            score = 1.0
        elif hazard_fraction < 0.3:
            hazard_level = HazardLevel.MODERATE
            score = 0.7
        elif hazard_fraction < 0.6:
            hazard_level = HazardLevel.HIGH
            score = 0.4
        else:
            hazard_level = HazardLevel.EXTREME
            score = 0.15

        if requires_special:
            score *= 0.8

        return SafetyAssessment(
            score=score,
            hazard_level=hazard_level,
            hazardous_steps=list(set(hazardous_steps)),
            ghs_hazards_encountered=list(set(ghs_hazards)),
            requires_special_equipment=requires_special,
        )

    def score_complexity(self, route: SynthesisRoute) -> SyntheticComplexityScore:
        """Assess synthetic complexity."""
        n_steps = route.total_steps
        lls = route.longest_linear_sequence
        overall_yield = route.overall_yield
        convergent = lls < n_steps  # If LLS < total steps, route has convergent branches

        # Step count score (fewer = better)
        step_score = max(0.0, 1.0 - (n_steps - 1) * 0.08)

        # Yield score
        yield_score = min(1.0, overall_yield * 2.0)

        # Convergence bonus
        convergence_bonus = 0.1 if convergent else 0.0

        # Max step difficulty (simplified)
        max_difficulty = 0.0
        for step in route.steps:
            yield_val = step.expected_yield if hasattr(step, "expected_yield") else 0.7
            step_difficulty = 1.0 - yield_val
            max_difficulty = max(max_difficulty, step_difficulty)

        difficulty_penalty = max_difficulty * 0.3

        complexity_score = (
            step_score * 0.35
            + yield_score * 0.35
            + (1.0 - difficulty_penalty) * 0.20
            + convergence_bonus
        )
        complexity_score = max(0.0, min(1.0, complexity_score))

        return SyntheticComplexityScore(
            score=complexity_score,
            total_steps=n_steps,
            longest_linear_sequence=lls,
            overall_yield=overall_yield,
            convergent=convergent,
            max_step_difficulty=max_difficulty,
        )

    def score_regulatory(
        self, route: SynthesisRoute, target_smiles: str
    ) -> RegulatoryAssessment:
        """Assess regulatory constraints."""
        controlled_precursors: list[str] = []
        requires_dea = False

        # Check starting materials against controlled list
        for sm in route.starting_materials:
            smiles = sm.smiles if hasattr(sm, "smiles") else str(sm)
            for name, pattern in CONTROLLED_PATTERNS.items():
                if pattern and smiles == pattern:
                    controlled_precursors.append(name)
                    requires_dea = True

        # Classification
        if requires_dea:
            classification = RegulatoryClass.CONTROLLED
            score = 0.2
        elif controlled_precursors:
            classification = RegulatoryClass.RESTRICTED
            score = 0.4
        else:
            classification = RegulatoryClass.UNCLASSIFIED
            score = 1.0

        return RegulatoryAssessment(
            score=score,
            classification=classification,
            controlled_precursors=controlled_precursors,
            requires_dea_license=requires_dea,
            export_controlled=requires_dea,
        )

    def _compute_composite(self, scores: dict[str, float]) -> tuple[float, str]:
        """Compute weighted composite score and grade."""
        total = sum(
            scores[dim] * weight
            for dim, weight in self.DIMENSION_WEIGHTS.items()
            if dim in scores
        )
        composite = total * 100  # Scale to 0-100

        if composite >= 80:
            grade = "A"
        elif composite >= 65:
            grade = "B"
        elif composite >= 50:
            grade = "C"
        elif composite >= 35:
            grade = "D"
        else:
            grade = "F"

        return composite, grade

    def score(
        self,
        target_smiles: str,
        tree: RetrosynthesisTree | None = None,
        route: SynthesisRoute | None = None,
    ) -> FeasibilityResult:
        """Run complete feasibility assessment.

        Can accept pre-computed retrosynthesis tree/route, or will
        compute them from the target SMILES.
        """
        # Run retrosynthesis if not provided
        if tree is None:
            from molbuilder.smiles.parser import parse
            mol = parse(target_smiles)
            tree = retrosynthesis(mol, self.max_depth, self.beam_width)

        if route is None:
            route = extract_best_route(tree)

        # Score each dimension
        availability = self.score_availability(route)
        cost = self.score_cost(route)
        green = self.score_green_chemistry(route)
        safety = self.score_safety(route)
        complexity = self.score_complexity(route)
        regulatory = self.score_regulatory(route, target_smiles)

        # Composite
        dim_scores = {
            "availability": availability.score,
            "cost": cost.score,
            "green_chemistry": green.score,
            "safety": safety.score,
            "complexity": complexity.score,
            "regulatory": regulatory.score,
        }
        composite, grade = self._compute_composite(dim_scores)

        route_summary = {
            "total_steps": route.total_steps,
            "overall_yield": round(route.overall_yield, 4),
            "starting_materials": [
                sm.name if hasattr(sm, "name") else str(sm)
                for sm in route.starting_materials
            ],
            "routes_found": tree.routes_found,
        }

        return FeasibilityResult(
            target_smiles=target_smiles,
            composite_score=composite,
            grade=grade,
            availability=availability,
            cost=cost,
            green_chemistry=green,
            safety=safety,
            complexity=complexity,
            regulatory=regulatory,
            route_summary=route_summary,
        )

    def compare(
        self,
        smiles_list: list[str],
    ) -> list[FeasibilityResult]:
        """Compare feasibility of multiple target molecules."""
        results = []
        for smiles in smiles_list:
            try:
                result = self.score(smiles)
                results.append(result)
            except Exception:
                # Return a zero-score result for failed analyses
                results.append(FeasibilityResult(
                    target_smiles=smiles,
                    composite_score=0.0,
                    grade="F",
                    availability=AvailabilityScore(0, 0, 0, 0, 365, "high"),
                    cost=CostEstimate(0, 0, 0),
                    green_chemistry=GreenChemistryScore(0, 0, 0, 0, 0, 100),
                    safety=SafetyAssessment(0, HazardLevel.EXTREME, [], [], True),
                    complexity=SyntheticComplexityScore(0, 0, 0, 0, False, 1.0),
                    regulatory=RegulatoryAssessment(0, RegulatoryClass.UNCLASSIFIED, [], False, False),
                ))
        return sorted(results, key=lambda r: r.composite_score, reverse=True)
