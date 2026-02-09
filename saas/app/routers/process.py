"""Process engineering endpoint."""

from fastapi import APIRouter, Depends
from app.dependencies import UserContext, check_expensive_rate_limit
from app.exceptions import InvalidSMILES
from app.models.process import ProcessEvaluateRequest, ProcessEvaluateResponse
from app.services.process_service import evaluate_process

router = APIRouter(prefix="/api/v1/process", tags=["process"])


@router.post("/evaluate", response_model=ProcessEvaluateResponse)
def evaluate(
    body: ProcessEvaluateRequest,
    user: UserContext = Depends(check_expensive_rate_limit),
):
    try:
        return evaluate_process(
            body.smiles,
            scale_kg=body.scale_kg,
            max_depth=body.max_depth,
            beam_width=body.beam_width,
        )
    except Exception as e:
        raise InvalidSMILES(body.smiles, str(e))
