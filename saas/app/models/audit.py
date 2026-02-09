"""Audit trail request/response models."""

from pydantic import BaseModel, Field


class AuditRecordResponse(BaseModel):
    id: int
    timestamp: float
    user_email: str
    user_role: str
    user_tier: str
    action: str
    input_summary: str
    output_hash: str
    status_code: int
    latency_ms: float
    ip_address: str
    signature_hash: str


class AuditQueryResponse(BaseModel):
    records: list[AuditRecordResponse]
    total: int
    limit: int
    offset: int


class AuditVerifyResponse(BaseModel):
    valid: bool
    record_id: int
    expected_hash: str = ""
    stored_hash: str = ""
    error: str = ""
