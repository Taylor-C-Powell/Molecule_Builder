"""Process engineering service - full pipeline: retro + reactor/conditions/safety/cost/scaleup."""

from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeoutError

from molbuilder.smiles.parser import parse
from molbuilder.reactions.retrosynthesis import retrosynthesis
from molbuilder.reactions.synthesis_route import extract_best_route
from molbuilder.process.reactor import select_reactor
from molbuilder.process.conditions import optimize_conditions
from molbuilder.process.purification import recommend_purification
from molbuilder.process.safety import assess_safety
from molbuilder.process.costing import estimate_cost
from molbuilder.process.scale_up import analyze_scale_up
from app.models.process import (
    ProcessEvaluateResponse,
    ReactorResponse,
    ConditionsResponse,
    PurificationStepResponse,
    SafetyAssessmentResponse,
    HazardInfoResponse,
    CostEstimateResponse,
    CostBreakdownResponse,
    ScaleUpResponse,
    StepProcessDetail,
)


def _serialize_reactor(spec) -> ReactorResponse:
    return ReactorResponse(
        reactor_type=spec.reactor_type.name,
        volume_L=spec.volume_L,
        temperature_C=spec.temperature_C,
        pressure_atm=spec.pressure_atm,
        residence_time_min=spec.residence_time_min,
        mixing_type=spec.mixing_type,
        heat_transfer=spec.heat_transfer,
        material=spec.material,
        estimated_cost_usd=spec.estimated_cost_usd,
        notes=spec.notes,
    )


def _serialize_conditions(cond) -> ConditionsResponse:
    return ConditionsResponse(
        temperature_C=cond.temperature_C,
        pressure_atm=cond.pressure_atm,
        solvent=cond.solvent,
        concentration_M=cond.concentration_M,
        addition_rate=cond.addition_rate,
        reaction_time_hours=cond.reaction_time_hours,
        atmosphere=cond.atmosphere,
        workup_procedure=cond.workup_procedure,
        notes=cond.notes,
    )


def _serialize_purification(steps) -> list[PurificationStepResponse]:
    return [
        PurificationStepResponse(
            method=p.method.name,
            description=p.description,
            estimated_recovery=p.estimated_recovery,
            estimated_purity=p.estimated_purity,
            scale_appropriate=p.scale_appropriate,
            notes=p.notes,
        )
        for p in steps
    ]


def _serialize_safety(assessments) -> list[SafetyAssessmentResponse]:
    return [
        SafetyAssessmentResponse(
            step_number=a.step_number,
            step_name=a.step_name,
            risk_level=a.risk_level,
            ppe_required=a.ppe_required,
            engineering_controls=a.engineering_controls,
            emergency_procedures=a.emergency_procedures,
            incompatible_materials=a.incompatible_materials,
            waste_classification=a.waste_classification,
            hazards=[
                HazardInfoResponse(
                    reagent_name=h.reagent_name,
                    ghs_hazards=h.ghs_hazards,
                    ghs_pictograms=h.ghs_pictograms,
                    hazard_descriptions=h.hazard_descriptions,
                    pictogram_descriptions=h.pictogram_descriptions,
                )
                for h in a.hazards
            ],
        )
        for a in assessments
    ]


def _serialize_cost(cost) -> CostEstimateResponse:
    return CostEstimateResponse(
        total_usd=cost.total_usd,
        per_kg_usd=cost.per_kg_usd,
        scale_kg=cost.scale_kg,
        breakdown=CostBreakdownResponse(
            raw_materials_usd=cost.breakdown.raw_materials_usd,
            labor_usd=cost.breakdown.labor_usd,
            equipment_usd=cost.breakdown.equipment_usd,
            energy_usd=cost.breakdown.energy_usd,
            waste_disposal_usd=cost.breakdown.waste_disposal_usd,
            overhead_usd=cost.breakdown.overhead_usd,
        ),
        notes=cost.notes,
    )


def _serialize_scaleup(su) -> ScaleUpResponse:
    return ScaleUpResponse(
        target_annual_kg=su.target_annual_kg,
        recommended_mode=su.recommended_mode,
        batch_size_kg=su.batch_size_kg,
        batches_per_year=su.batches_per_year,
        cycle_time_hours=su.cycle_time_hours,
        annual_capacity_kg=su.annual_capacity_kg,
        capital_cost_usd=su.capital_cost_usd,
        operating_cost_annual_usd=su.operating_cost_annual_usd,
        scale_up_risks=su.scale_up_risks,
        recommendations=su.recommendations,
    )


PROCESS_TIMEOUT_SECONDS = 60

_executor = ThreadPoolExecutor(max_workers=2)


def _evaluate_process_sync(
    smiles: str,
    scale_kg: float,
    max_depth: int,
    beam_width: int,
) -> ProcessEvaluateResponse:
    """Core process evaluation logic (runs in thread pool)."""
    mol = parse(smiles)
    tree = retrosynthesis(mol, max_depth=max_depth, beam_width=beam_width)

    if not tree.target.best_disconnection:
        return ProcessEvaluateResponse(
            smiles=smiles, scale_kg=scale_kg, route_found=False
        )

    route = extract_best_route(tree)

    # Per-step process details
    step_details = []
    for step in route.steps:
        t = step.template
        reactor = select_reactor(t, scale_kg)
        conditions = optimize_conditions(t, scale_kg)
        purif = recommend_purification(t, scale_kg)
        step_details.append(
            StepProcessDetail(
                step_number=step.step_number,
                reaction_name=t.name,
                reactor=_serialize_reactor(reactor),
                conditions=_serialize_conditions(conditions),
                purification=_serialize_purification(purif),
            )
        )

    # Aggregate assessments
    safety = _serialize_safety(assess_safety(route.steps))
    cost = _serialize_cost(estimate_cost(route.steps, scale_kg))
    annual_kg = scale_kg * 100
    scaleup = _serialize_scaleup(analyze_scale_up(route.steps, annual_kg))

    return ProcessEvaluateResponse(
        smiles=smiles,
        scale_kg=scale_kg,
        route_found=True,
        total_steps=route.total_steps,
        overall_yield=route.overall_yield,
        step_details=step_details,
        safety=safety,
        cost=cost,
        scale_up=scaleup,
    )


def evaluate_process(
    smiles: str,
    scale_kg: float = 1.0,
    max_depth: int = 5,
    beam_width: int = 5,
) -> ProcessEvaluateResponse:
    """Full pipeline with timeout protection: parse -> retro -> process engineering."""
    future = _executor.submit(
        _evaluate_process_sync, smiles, scale_kg, max_depth, beam_width
    )
    try:
        return future.result(timeout=PROCESS_TIMEOUT_SECONDS)
    except FuturesTimeoutError:
        future.cancel()
        raise ValueError(
            f"Process evaluation timed out after {PROCESS_TIMEOUT_SECONDS}s "
            f"for SMILES: {smiles}"
        )
