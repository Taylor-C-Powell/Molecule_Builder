"""Retrosynthesis endpoint."""

from fastapi import APIRouter, Depends
from app.dependencies import UserContext, check_expensive_rate_limit
from app.exceptions import InvalidSMILES
from app.models.retrosynthesis import RetroRequest, RetroResponse
from app.services.retro_service import run_retrosynthesis

router = APIRouter(prefix="/api/v1/retrosynthesis", tags=["retrosynthesis"])


@router.post("/plan", response_model=RetroResponse)
def plan(
    body: RetroRequest,
    user: UserContext = Depends(check_expensive_rate_limit),
):
    try:
        return run_retrosynthesis(
            body.smiles,
            max_depth=body.max_depth,
            beam_width=body.beam_width,
        )
    except (ValueError, KeyError) as e:
        raise InvalidSMILES(body.smiles, str(e))
