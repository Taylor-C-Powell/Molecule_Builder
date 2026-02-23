"""Condition prediction service -- thin wrapper around core library."""

from molbuilder.process.condition_prediction import (
    predict_conditions,
    ConditionPrediction,
    TemplateMatch,
)
from app.models.process import (
    PredictConditionsResponse,
    SubstrateAnalysisResponse,
    TemplateMatchResponse,
    ConditionsResponse,
    SolventScoreResponse,
)


def _serialize_match(m: TemplateMatch) -> TemplateMatchResponse:
    return TemplateMatchResponse(
        template_name=m.template_name,
        named_reaction=m.named_reaction,
        category=m.category,
        match_score=m.match_score,
        match_reasons=m.match_reasons,
        conditions=ConditionsResponse(
            temperature_C=m.conditions.temperature_C,
            pressure_atm=m.conditions.pressure_atm,
            solvent=m.conditions.solvent,
            concentration_M=m.conditions.concentration_M,
            addition_rate=m.conditions.addition_rate,
            reaction_time_hours=m.conditions.reaction_time_hours,
            atmosphere=m.conditions.atmosphere,
            workup_procedure=m.conditions.workup_procedure,
            notes=m.conditions.notes,
        ),
        recommended_solvents=[
            SolventScoreResponse(
                solvent_name=s.solvent_name,
                composite_score=s.composite_score,
                green_score=s.green_score,
                reasons=s.reasons,
            )
            for s in m.recommended_solvents
        ],
        adjusted_yield_range=list(m.adjusted_yield_range),
        warnings=m.warnings,
    )


def run_prediction(
    smiles: str,
    reaction_name: str | None = None,
    scale_kg: float = 1.0,
    max_candidates: int = 5,
) -> PredictConditionsResponse:
    """Run condition prediction and serialize to Pydantic response."""
    result = predict_conditions(
        smiles,
        reaction_name=reaction_name,
        scale_kg=scale_kg,
        max_candidates=max_candidates,
    )
    sa = result.substrate_analysis
    return PredictConditionsResponse(
        smiles=result.smiles,
        substrate_analysis=SubstrateAnalysisResponse(
            detected_functional_groups=sa.detected_functional_groups,
            molecular_weight=sa.molecular_weight,
            heavy_atom_count=sa.heavy_atom_count,
            steric_class=sa.steric_class,
            electronic_character=sa.electronic_character,
            sensitive_groups=sa.sensitive_groups,
        ),
        candidates=[_serialize_match(m) for m in result.candidates],
        best_match=_serialize_match(result.best_match) if result.best_match else None,
        overall_confidence=result.overall_confidence,
    )
