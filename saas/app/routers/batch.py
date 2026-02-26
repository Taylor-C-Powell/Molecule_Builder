"""Batch processing endpoints."""

from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException

from app.config import Tier
from app.dependencies import UserContext, check_rate_limit, check_expensive_rate_limit
from app.models.batch import (
    BatchSubmitRequest,
    BatchSubmitResponse,
    BatchStatusResponse,
    BatchListResponse,
    BatchJobSummary,
)
from app.services.job_db import get_job_db
from app.services.batch_worker import get_batch_worker

router = APIRouter(prefix="/api/v1/batch", tags=["batch"])

_TIER_BATCH_LIMITS = {
    Tier.FREE: 10,
    Tier.PRO: 100,
    Tier.TEAM: 500,
    Tier.ACADEMIC: 50,
    Tier.ENTERPRISE: 10000,
}


@router.post("/submit", response_model=BatchSubmitResponse)
def submit_batch(
    body: BatchSubmitRequest,
    user: UserContext = Depends(check_expensive_rate_limit),
):
    limit = _TIER_BATCH_LIMITS.get(user.tier, 10)
    if len(body.smiles_list) > limit:
        raise HTTPException(
            status_code=403,
            detail=f"Batch size {len(body.smiles_list)} exceeds {user.tier.value} tier limit of {limit}",
        )

    db = get_job_db()
    job_id = db.create_job(
        user_email=user.email,
        job_type=body.job_type.value,
        input_data={"smiles_list": body.smiles_list, "params": body.params or {}},
    )

    worker = get_batch_worker()
    worker.submit_batch(job_id, body.smiles_list, body.job_type.value, body.params or {})

    job = db.get_job(job_id)
    return BatchSubmitResponse(
        job_id=job_id,
        status=job["status"],
        created_at=job["created_at"],
    )


@router.get("/{job_id}", response_model=BatchStatusResponse)
def get_job_status(
    job_id: str,
    user: UserContext = Depends(check_rate_limit),
):
    db = get_job_db()
    job = db.get_job(job_id)
    if job is None or job["user_email"] != user.email:
        raise HTTPException(status_code=404, detail="Job not found")

    return BatchStatusResponse(
        job_id=job["job_id"],
        status=job["status"],
        job_type=job["job_type"],
        progress_pct=job["progress_pct"],
        result=job.get("result_data"),
        error=job.get("error"),
        created_at=job["created_at"],
        updated_at=job["updated_at"],
    )


@router.get("/", response_model=BatchListResponse)
def list_jobs(
    page: int = 1,
    per_page: int = 20,
    user: UserContext = Depends(check_rate_limit),
):
    db = get_job_db()
    offset = (page - 1) * per_page
    jobs, total = db.list_jobs(user.email, limit=per_page, offset=offset)
    return BatchListResponse(
        jobs=[BatchJobSummary(**j) for j in jobs],
        total=total,
        page=page,
        per_page=per_page,
    )


@router.delete("/{job_id}")
def cancel_job(
    job_id: str,
    user: UserContext = Depends(check_rate_limit),
):
    db = get_job_db()
    cancelled = db.cancel_job(job_id, user.email)
    if not cancelled:
        raise HTTPException(status_code=404, detail="Job not found or not cancellable")
    worker = get_batch_worker()
    worker.cancel(job_id)
    return {"status": "cancelled", "job_id": job_id}
