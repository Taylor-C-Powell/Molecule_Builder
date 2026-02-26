"""Process engineering endpoints."""

from fastapi import APIRouter, Depends, Query
from app.dependencies import UserContext, check_expensive_rate_limit
from app.exceptions import InvalidSMILES
from app.models.process import (
    ProcessEvaluateRequest,
    ProcessEvaluateResponse,
    PredictConditionsRequest,
    PredictConditionsResponse,
)
from app.services.process_service import evaluate_process
from app.services.prediction_service import run_prediction

router = APIRouter(prefix="/api/v1/process", tags=["process"])


@router.post("/evaluate", response_model=ProcessEvaluateResponse)
def evaluate(
    body: ProcessEvaluateRequest,
    cost_region: str = Query("us", pattern="^(us|eu|india|china)$"),
    user: UserContext = Depends(check_expensive_rate_limit),
):
    try:
        return evaluate_process(
            body.smiles,
            scale_kg=body.scale_kg,
            max_depth=body.max_depth,
            beam_width=body.beam_width,
            cost_region=cost_region,
        )
    except (ValueError, KeyError) as e:
        raise InvalidSMILES(body.smiles, str(e))


@router.post("/predict-conditions", response_model=PredictConditionsResponse)
def predict_conditions_endpoint(
    body: PredictConditionsRequest,
    user: UserContext = Depends(check_expensive_rate_limit),
):
    try:
        return run_prediction(
            body.smiles,
            reaction_name=body.reaction_name,
            scale_kg=body.scale_kg,
            max_candidates=body.max_candidates,
        )
    except (ValueError, KeyError) as e:
        raise InvalidSMILES(body.smiles, str(e))
