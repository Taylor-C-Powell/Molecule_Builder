"""Tests for retrosynthesis feasibility scoring engine.

Tests the FeasibilityEngine and its component scoring functions
using synthetic route data and known molecules.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pytest

from molbuilder.reactions.feasibility import (
    FeasibilityEngine,
    FeasibilityResult,
    AvailabilityScore,
    CostEstimate,
    GreenChemistryScore,
    SafetyAssessment,
    SyntheticComplexityScore,
    RegulatoryAssessment,
    RegulatoryClass,
    HazardLevel,
    GREEN_SOLVENT_SCORES,
)


# Mock objects for testing without full retrosynthesis
@dataclass
class MockPrecursor:
    smiles: str
    name: str = "unknown"
    cost_per_kg: float = 50.0


@dataclass
class MockStep:
    step_number: int = 1
    conditions: str = ""
    expected_yield: float = 0.7
    template: Any = None
    precursors: list = field(default_factory=list)
    product_smiles: str = ""
    product_name: str = ""
    notes: str = ""


@dataclass
class MockRoute:
    target_smiles: str = "CCO"
    target_name: str = "ethanol"
    steps: list = field(default_factory=list)
    overall_yield: float = 0.5
    starting_materials: list = field(default_factory=list)
    total_steps: int = 3
    longest_linear_sequence: int = 3


def _make_simple_route(n_steps=3, yield_per_step=0.8, n_purchasable=3, n_total=4):
    """Create a mock route for testing."""
    steps = []
    for i in range(n_steps):
        steps.append(MockStep(
            step_number=i + 1,
            conditions="ethanol, room temperature",
            expected_yield=yield_per_step,
        ))

    materials = []
    # Use actual purchasable SMILES from the database
    purchasable_smiles = ["O", "CCO", "CC=O", "CC(=O)C"]  # water, ethanol, acetaldehyde, acetone
    for i in range(n_total):
        if i < n_purchasable and i < len(purchasable_smiles):
            materials.append(MockPrecursor(purchasable_smiles[i], f"reagent_{i}", 30.0))
        else:
            materials.append(MockPrecursor(f"CCCCCCCCCCCCCC(=O)N{i}", f"custom_{i}", 500.0))

    return MockRoute(
        steps=steps,
        overall_yield=yield_per_step ** n_steps,
        starting_materials=materials,
        total_steps=n_steps,
        longest_linear_sequence=n_steps,
    )


def _make_hazardous_route():
    """Create a route with hazardous steps."""
    return MockRoute(
        steps=[
            MockStep(1, "LiAlH4 in THF, 0C"),
            MockStep(2, "BuLi in THF, -78C"),
            MockStep(3, "HCN, catalytic KCN"),
        ],
        overall_yield=0.3,
        starting_materials=[MockPrecursor("O", "water", 1.0)],
        total_steps=3,
        longest_linear_sequence=3,
    )


class TestFeasibilityEngine:
    def test_init_default_weights(self):
        engine = FeasibilityEngine()
        assert engine.DIMENSION_WEIGHTS["availability"] == 0.25
        assert sum(engine.DIMENSION_WEIGHTS.values()) == pytest.approx(1.0)

    def test_init_custom_weights(self):
        engine = FeasibilityEngine(custom_weights={"availability": 0.5, "cost": 0.5})
        assert engine.DIMENSION_WEIGHTS["availability"] == 0.5


class TestAvailabilityScoring:
    def test_all_purchasable(self):
        engine = FeasibilityEngine()
        route = _make_simple_route(n_purchasable=3, n_total=3)
        score = engine.score_availability(route)
        assert score.fraction_purchasable == 1.0
        assert score.score >= 0.9
        assert score.supply_chain_risk == "low"

    def test_none_purchasable(self):
        engine = FeasibilityEngine()
        route = _make_simple_route(n_purchasable=0, n_total=3)
        score = engine.score_availability(route)
        assert score.fraction_purchasable == 0.0
        assert score.score == 0.0
        assert score.supply_chain_risk == "high"

    def test_partial_purchasable(self):
        engine = FeasibilityEngine()
        route = _make_simple_route(n_purchasable=2, n_total=4)
        score = engine.score_availability(route)
        assert 0 < score.score < 1
        assert score.n_purchasable == 2
        assert score.n_total_precursors == 4

    def test_empty_materials(self):
        engine = FeasibilityEngine()
        route = MockRoute(starting_materials=[])
        score = engine.score_availability(route)
        assert score.score == 0.0

    def test_lead_time_scales(self):
        engine = FeasibilityEngine()
        full = engine.score_availability(_make_simple_route(n_purchasable=4, n_total=4))
        partial = engine.score_availability(_make_simple_route(n_purchasable=1, n_total=4))
        assert full.estimated_lead_time_days < partial.estimated_lead_time_days


class TestCostScoring:
    def test_cheap_route(self):
        engine = FeasibilityEngine()
        route = _make_simple_route(n_steps=2)
        score = engine.score_cost(route)
        assert score.score > 0
        assert score.total_material_cost_usd >= 0

    def test_expensive_route(self):
        engine = FeasibilityEngine(cost_threshold_per_gram=1.0)
        route = _make_simple_route(n_steps=10)
        route.overall_yield = 0.001
        score = engine.score_cost(route)
        assert score.score < 0.5

    def test_labor_scales_with_steps(self):
        engine = FeasibilityEngine()
        short = engine.score_cost(_make_simple_route(n_steps=2))
        long = engine.score_cost(_make_simple_route(n_steps=8))
        assert short.labor_cost_estimate < long.labor_cost_estimate


class TestGreenChemistryScoring:
    def test_green_route(self):
        engine = FeasibilityEngine()
        route = MockRoute(
            steps=[MockStep(1, "water, room temperature", 0.9)],
            overall_yield=0.9,
            starting_materials=[MockPrecursor("O")],
            total_steps=1,
            longest_linear_sequence=1,
        )
        score = engine.score_green_chemistry(route)
        assert score.score > 0.6

    def test_dirty_route(self):
        engine = FeasibilityEngine()
        route = MockRoute(
            steps=[
                MockStep(1, "dichloromethane, reflux", 0.5),
                MockStep(2, "chloroform, LiAlH4", 0.4),
                MockStep(3, "benzene, HCN", 0.3),
            ],
            overall_yield=0.06,
            starting_materials=[MockPrecursor("C")],
            total_steps=3,
            longest_linear_sequence=3,
        )
        score = engine.score_green_chemistry(route)
        assert score.score < 0.6
        assert score.hazardous_reagent_fraction > 0

    def test_e_factor(self):
        engine = FeasibilityEngine()
        high_yield = MockRoute(steps=[MockStep()], overall_yield=0.9,
                               starting_materials=[], total_steps=1, longest_linear_sequence=1)
        low_yield = MockRoute(steps=[MockStep()], overall_yield=0.1,
                              starting_materials=[], total_steps=1, longest_linear_sequence=1)

        green_high = engine.score_green_chemistry(high_yield)
        green_low = engine.score_green_chemistry(low_yield)
        assert green_high.e_factor_estimate < green_low.e_factor_estimate


class TestSafetyScoring:
    def test_safe_route(self):
        engine = FeasibilityEngine()
        route = _make_simple_route()
        score = engine.score_safety(route)
        assert score.hazard_level == HazardLevel.LOW
        assert score.score >= 0.9

    def test_hazardous_route(self):
        engine = FeasibilityEngine()
        route = _make_hazardous_route()
        score = engine.score_safety(route)
        assert score.hazard_level in (HazardLevel.HIGH, HazardLevel.EXTREME)
        assert score.score < 0.5
        assert len(score.ghs_hazards_encountered) > 0
        assert score.requires_special_equipment is True

    def test_partial_hazard(self):
        engine = FeasibilityEngine()
        route = MockRoute(
            steps=[
                MockStep(1, "ethanol, room temperature"),
                MockStep(2, "LiAlH4 in THF"),
                MockStep(3, "water, room temperature"),
                MockStep(4, "ethanol, reflux"),
            ],
            overall_yield=0.5,
            starting_materials=[],
            total_steps=4,
            longest_linear_sequence=4,
        )
        score = engine.score_safety(route)
        # 1/4 = 25% hazardous, which is < 0.3 threshold for MODERATE
        assert score.hazard_level == HazardLevel.MODERATE


class TestComplexityScoring:
    def test_simple_route(self):
        engine = FeasibilityEngine()
        route = MockRoute(total_steps=1, longest_linear_sequence=1,
                          overall_yield=0.9, steps=[MockStep(expected_yield=0.9)],
                          starting_materials=[])
        score = engine.score_complexity(route)
        assert score.score > 0.7

    def test_complex_route(self):
        engine = FeasibilityEngine()
        route = MockRoute(total_steps=12, longest_linear_sequence=12,
                          overall_yield=0.01,
                          steps=[MockStep(expected_yield=0.5) for _ in range(12)],
                          starting_materials=[])
        score = engine.score_complexity(route)
        assert score.score < 0.3

    def test_convergent_bonus(self):
        engine = FeasibilityEngine()
        linear = MockRoute(total_steps=5, longest_linear_sequence=5,
                           overall_yield=0.3, steps=[MockStep() for _ in range(5)],
                           starting_materials=[])
        convergent = MockRoute(total_steps=5, longest_linear_sequence=3,
                               overall_yield=0.3, steps=[MockStep() for _ in range(5)],
                               starting_materials=[])

        s_linear = engine.score_complexity(linear)
        s_conv = engine.score_complexity(convergent)
        assert s_conv.convergent is True
        assert s_linear.convergent is False
        assert s_conv.score >= s_linear.score


class TestRegulatoryScoring:
    def test_clean_route(self):
        engine = FeasibilityEngine()
        route = _make_simple_route()
        score = engine.score_regulatory(route, "CCO")
        assert score.score >= 0.9
        assert score.classification == RegulatoryClass.UNCLASSIFIED

    def test_controlled_precursor(self):
        engine = FeasibilityEngine()
        route = MockRoute(
            starting_materials=[
                MockPrecursor("C[C@@H](O)[C@H](NC)c1ccccc1", "ephedrine", 100.0),
            ],
            steps=[], total_steps=0, longest_linear_sequence=0, overall_yield=1.0,
        )
        score = engine.score_regulatory(route, "target")
        assert score.score < 0.5
        assert score.requires_dea_license is True
        assert "ephedrine" in score.controlled_precursors


class TestCompositeScoring:
    def test_compute_composite(self):
        engine = FeasibilityEngine()
        scores = {
            "availability": 0.8,
            "cost": 0.7,
            "green_chemistry": 0.6,
            "safety": 0.9,
            "complexity": 0.7,
            "regulatory": 1.0,
        }
        composite, grade = engine._compute_composite(scores)
        assert 0 <= composite <= 100
        assert grade in ("A", "B", "C", "D", "F")

    def test_all_perfect(self):
        engine = FeasibilityEngine()
        scores = {d: 1.0 for d in engine.DIMENSION_WEIGHTS}
        composite, grade = engine._compute_composite(scores)
        assert composite == pytest.approx(100.0)
        assert grade == "A"

    def test_all_zero(self):
        engine = FeasibilityEngine()
        scores = {d: 0.0 for d in engine.DIMENSION_WEIGHTS}
        composite, grade = engine._compute_composite(scores)
        assert composite == pytest.approx(0.0)
        assert grade == "F"

    def test_grade_boundaries(self):
        engine = FeasibilityEngine()
        # Test each grade boundary
        for target_score, expected_grade in [(85, "A"), (70, "B"), (55, "C"), (40, "D"), (20, "F")]:
            val = target_score / 100.0
            scores = {d: val for d in engine.DIMENSION_WEIGHTS}
            _, grade = engine._compute_composite(scores)
            assert grade == expected_grade, f"Expected {expected_grade} for score {target_score}"


class TestFeasibilityResult:
    def test_is_feasible(self):
        result = FeasibilityResult(
            target_smiles="CCO",
            composite_score=60.0,
            grade="B",
            availability=AvailabilityScore(0.8, 3, 3, 1.0, 3, "low"),
            cost=CostEstimate(0.7, 10, 50),
            green_chemistry=GreenChemistryScore(0.6, 0.8, 0.7, 0.8, 0.1, 5),
            safety=SafetyAssessment(0.9, HazardLevel.LOW, [], [], False),
            complexity=SyntheticComplexityScore(0.7, 3, 3, 0.5, False, 0.3),
            regulatory=RegulatoryAssessment(1.0, RegulatoryClass.UNCLASSIFIED, [], False, False),
        )
        assert result.is_feasible is True

    def test_not_feasible(self):
        result = FeasibilityResult(
            target_smiles="CCO",
            composite_score=20.0,
            grade="F",
            availability=AvailabilityScore(0.1, 0, 3, 0.0, 90, "high"),
            cost=CostEstimate(0.1, 1000, 5000),
            green_chemistry=GreenChemistryScore(0.2, 0.1, 0.2, 0.3, 0.8, 80),
            safety=SafetyAssessment(0.2, HazardLevel.EXTREME, [0, 1], ["GHS06"], True),
            complexity=SyntheticComplexityScore(0.1, 12, 12, 0.01, False, 0.9),
            regulatory=RegulatoryAssessment(0.2, RegulatoryClass.CONTROLLED, ["X"], True, True),
        )
        assert result.is_feasible is False

    def test_to_dict(self):
        result = FeasibilityResult(
            target_smiles="CCO",
            composite_score=65.0,
            grade="B",
            availability=AvailabilityScore(0.8, 3, 3, 1.0, 3, "low"),
            cost=CostEstimate(0.7, 10, 50),
            green_chemistry=GreenChemistryScore(0.6, 0.8, 0.7, 0.8, 0.1, 5),
            safety=SafetyAssessment(0.9, HazardLevel.LOW, [], [], False),
            complexity=SyntheticComplexityScore(0.7, 3, 3, 0.5, False, 0.3),
            regulatory=RegulatoryAssessment(1.0, RegulatoryClass.UNCLASSIFIED, [], False, False),
        )
        d = result.to_dict()
        assert d["composite_score"] == 65.0
        assert d["grade"] == "B"
        assert "dimensions" in d
        assert "availability" in d["dimensions"]


class TestGreenSolventScores:
    def test_water_is_greenest(self):
        assert GREEN_SOLVENT_SCORES["water"] == 1.0

    def test_benzene_is_worst(self):
        assert GREEN_SOLVENT_SCORES["benzene"] <= 0.15

    def test_all_in_range(self):
        for solvent, score in GREEN_SOLVENT_SCORES.items():
            assert 0.0 <= score <= 1.0, f"{solvent} has invalid score {score}"


class TestFullScoreWithMockRoute:
    """Test full score() with mock route (no retrosynthesis)."""

    def test_score_with_route(self):
        engine = FeasibilityEngine()
        route = _make_simple_route()

        # We need a mock tree too
        @dataclass
        class MockTree:
            routes_found: int = 1

        result = engine.score("CCO", tree=MockTree(), route=route)
        assert isinstance(result, FeasibilityResult)
        assert 0 <= result.composite_score <= 100
        assert result.grade in ("A", "B", "C", "D", "F")

    def test_score_dimensions_present(self):
        engine = FeasibilityEngine()
        route = _make_simple_route()

        @dataclass
        class MockTree:
            routes_found: int = 1

        result = engine.score("CCO", tree=MockTree(), route=route)
        assert result.availability.score >= 0
        assert result.cost.score >= 0
        assert result.green_chemistry.score >= 0
        assert result.safety.score >= 0
        assert result.complexity.score >= 0
        assert result.regulatory.score >= 0
