"""Batch API request/response models."""

from __future__ import annotations

from enum import Enum
from pydantic import BaseModel, Field


class BatchJobType(str, Enum):
    PROPERTIES = "properties"
    RETROSYNTHESIS = "retrosynthesis"
    CONDITIONS = "conditions"
    EVALUATE = "evaluate"


class BatchSubmitRequest(BaseModel):
    smiles_list: list[str] = Field(..., min_length=1, description="List of SMILES strings to process")
    job_type: BatchJobType = Field(..., description="Type of analysis to run")
    params: dict | None = Field(default=None, description="Optional parameters for the job type")


class BatchSubmitResponse(BaseModel):
    job_id: str
    status: str
    created_at: str


class BatchStatusResponse(BaseModel):
    job_id: str
    status: str
    job_type: str
    progress_pct: float
    result: dict | None = None
    error: str | None = None
    created_at: str
    updated_at: str


class BatchJobSummary(BaseModel):
    job_id: str
    status: str
    job_type: str
    progress_pct: float
    created_at: str
    updated_at: str


class BatchListResponse(BaseModel):
    jobs: list[BatchJobSummary]
    total: int
    page: int
    per_page: int
