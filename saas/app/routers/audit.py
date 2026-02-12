"""Audit trail endpoints (admin only)."""

from fastapi import APIRouter, Depends, Query
from app.dependencies import UserContext, require_admin
from app.models.audit import AuditRecordResponse, AuditQueryResponse, AuditVerifyResponse
from app.services.audit_db import get_audit_db

router = APIRouter(prefix="/api/v1/audit", tags=["audit"])


@router.get("/log", response_model=AuditQueryResponse)
def query_audit_log(
    start_time: float | None = Query(None, description="Unix timestamp lower bound"),
    end_time: float | None = Query(None, description="Unix timestamp upper bound"),
    user_email: str | None = Query(None),
    action: str | None = Query(None),
    limit: int = Query(100, ge=1, le=1000),
    offset: int = Query(0, ge=0),
    admin: UserContext = Depends(require_admin),
):
    db = get_audit_db()
    records = db.query(
        start_time=start_time,
        end_time=end_time,
        user_email=user_email,
        action=action,
        limit=limit,
        offset=offset,
    )
    total = db.count(
        start_time=start_time,
        end_time=end_time,
        user_email=user_email,
        action=action,
    )
    return AuditQueryResponse(
        records=[AuditRecordResponse(**r) for r in records],
        total=total,
        limit=limit,
        offset=offset,
    )


@router.get("/log/{record_id}", response_model=AuditRecordResponse)
def get_audit_record(
    record_id: int,
    admin: UserContext = Depends(require_admin),
):
    db = get_audit_db()
    record = db.get_record(record_id)
    if record is None:
        from app.exceptions import MolBuilderAPIError
        raise MolBuilderAPIError(404, f"Audit record {record_id} not found")
    return AuditRecordResponse(**record)


@router.get("/log/{record_id}/verify", response_model=AuditVerifyResponse)
def verify_audit_record(
    record_id: int,
    admin: UserContext = Depends(require_admin),
):
    db = get_audit_db()
    result = db.verify_integrity(record_id)
    return AuditVerifyResponse(**result)


@router.get("/export")
def export_audit_log(
    limit: int = Query(10000, ge=1, le=10000),
    admin: UserContext = Depends(require_admin),
):
    db = get_audit_db()
    return db.export_all(limit=limit)
