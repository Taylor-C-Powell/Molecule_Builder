"""Tests for the process evaluation endpoint."""

from __future__ import annotations

import httpx
import respx

from molbuilder_client import (
    Conditions,
    CostBreakdown,
    CostEstimate,
    MolBuilder,
    ProcessEvaluation,
    PurificationMethod,
    Reactor,
    SafetyAssessment,
    ScaleUp,
    StepDetail,
)


PROCESS_RESPONSE = {
    "smiles": "CCO",
    "scale_kg": 10.0,
    "route_found": True,
    "total_steps": 1,
    "overall_yield": 0.85,
    "step_details": [
        {
            "step_number": 1,
            "reaction_name": "hydration",
            "reactor": {
                "reactor_type": "CSTR",
                "volume_L": 100.0,
                "temperature_C": 80.0,
                "pressure_atm": 1.0,
                "residence_time_min": 60.0,
                "mixing_type": "impeller",
                "heat_transfer": "jacket",
                "material": "stainless_steel",
                "estimated_cost_usd": 5000.0,
                "notes": "",
            },
            "conditions": {
                "temperature_C": 80.0,
                "pressure_atm": 1.0,
                "solvent": "water",
                "concentration_M": 1.0,
                "addition_rate": "dropwise",
                "reaction_time_hours": 2.0,
                "atmosphere": "air",
                "workup_procedure": "extraction",
                "notes": "",
            },
            "purification": [
                {
                    "method": "distillation",
                    "description": "Simple distillation at 78C",
                    "estimated_recovery": 0.95,
                    "estimated_purity": 0.99,
                    "scale_appropriate": True,
                    "notes": "",
                },
            ],
        },
    ],
    "safety": [
        {
            "step_number": 1,
            "step_name": "hydration",
            "risk_level": "low",
            "ppe_required": ["goggles", "gloves"],
            "engineering_controls": ["fume hood"],
            "emergency_procedures": ["eye wash"],
            "incompatible_materials": [],
            "waste_classification": "non-hazardous",
            "hazards": [],
        },
    ],
    "cost": {
        "total_usd": 500.0,
        "per_kg_usd": 50.0,
        "scale_kg": 10.0,
        "breakdown": {
            "raw_materials_usd": 100.0,
            "labor_usd": 150.0,
            "equipment_usd": 50.0,
            "energy_usd": 30.0,
            "waste_disposal_usd": 20.0,
            "overhead_usd": 150.0,
        },
        "notes": ["Estimate based on lab scale"],
    },
    "scale_up": {
        "target_annual_kg": 1000.0,
        "recommended_mode": "batch",
        "batch_size_kg": 50.0,
        "batches_per_year": 20.0,
        "cycle_time_hours": 8.0,
        "annual_capacity_kg": 1000.0,
        "capital_cost_usd": 100000.0,
        "operating_cost_annual_usd": 50000.0,
        "scale_up_risks": ["exotherm control"],
        "recommendations": ["install temp monitoring"],
    },
}


def test_process_evaluate(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.post("/process/evaluate").mock(
        return_value=httpx.Response(200, json=PROCESS_RESPONSE)
    )

    result = client.process_evaluate("CCO", scale_kg=10.0)

    assert isinstance(result, ProcessEvaluation)
    assert result.smiles == "CCO"
    assert result.scale_kg == 10.0
    assert result.route_found is True
    assert result.total_steps == 1

    # Step details
    assert len(result.step_details) == 1
    step = result.step_details[0]
    assert isinstance(step, StepDetail)
    assert isinstance(step.reactor, Reactor)
    assert step.reactor.reactor_type == "CSTR"
    assert isinstance(step.conditions, Conditions)
    assert step.conditions.solvent == "water"
    assert len(step.purification) == 1
    assert isinstance(step.purification[0], PurificationMethod)
    assert step.purification[0].estimated_purity == 0.99

    # Safety
    assert len(result.safety) == 1
    assert isinstance(result.safety[0], SafetyAssessment)
    assert result.safety[0].risk_level == "low"
    assert "goggles" in result.safety[0].ppe_required

    # Cost
    assert isinstance(result.cost, CostEstimate)
    assert result.cost.total_usd == 500.0
    assert isinstance(result.cost.breakdown, CostBreakdown)
    assert result.cost.breakdown.labor_usd == 150.0

    # Scale-up
    assert isinstance(result.scale_up, ScaleUp)
    assert result.scale_up.recommended_mode == "batch"
    assert result.scale_up.batch_size_kg == 50.0


def test_process_evaluate_no_route(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.post("/process/evaluate").mock(
        return_value=httpx.Response(
            200,
            json={
                "smiles": "X",
                "scale_kg": 1.0,
                "route_found": False,
                "total_steps": 0,
                "overall_yield": 0.0,
                "step_details": [],
                "safety": [],
                "cost": None,
                "scale_up": None,
            },
        )
    )

    result = client.process_evaluate("X")
    assert result.route_found is False
    assert result.cost is None
    assert result.scale_up is None
