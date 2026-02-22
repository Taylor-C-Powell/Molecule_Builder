"""Feasibility scoring API router.

POST /api/v1/feasibility/score - Score a single target
POST /api/v1/feasibility/compare - Compare multiple targets
"""

from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException

from app.dependencies import UserContext, check_expensive_rate_limit
from app.models.feasibility import (
    FeasibilityRequest,
    FeasibilityCompareRequest,
    FeasibilityResponse,
    FeasibilityCompareResponse,
)
from app.services.feasibility_service import (
    run_feasibility,
    run_feasibility_compare,
)

router = APIRouter(prefix="/api/v1/feasibility", tags=["feasibility"])


@router.post("/score", response_model=FeasibilityResponse)
def score_feasibility(
    body: FeasibilityRequest,
    user: UserContext = Depends(check_expensive_rate_limit),
) -> FeasibilityResponse:
    """Score the synthetic feasibility of a target molecule.

    Runs retrosynthetic analysis and evaluates:
    - Precursor availability
    - Estimated cost
    - Green chemistry metrics
    - Safety assessment
    - Synthetic complexity
    - Regulatory classification
    """
    try:
        return run_feasibility(
            smiles=body.smiles,
            max_depth=body.max_depth,
            beam_width=body.beam_width,
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except KeyError as exc:
        raise HTTPException(
            status_code=400, detail=f"Invalid SMILES: {exc}"
        ) from exc


@router.post("/compare", response_model=FeasibilityCompareResponse)
def compare_feasibility(
    body: FeasibilityCompareRequest,
    user: UserContext = Depends(check_expensive_rate_limit),
) -> FeasibilityCompareResponse:
    """Compare synthetic feasibility of multiple target molecules.

    Returns results sorted by composite feasibility score (best first).
    """
    try:
        return run_feasibility_compare(
            smiles_list=body.smiles_list,
            max_depth=body.max_depth,
            beam_width=body.beam_width,
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
