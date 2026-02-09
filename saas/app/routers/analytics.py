"""Analytics endpoints for usage tracking (admin only)."""

from fastapi import APIRouter, Depends, Query
from app.dependencies import UserContext, require_admin
from app.services.usage_tracker import usage_tracker

router = APIRouter(prefix="/api/v1/analytics", tags=["analytics"])


@router.get("/summary")
def get_summary(user: UserContext = Depends(require_admin)):
    return usage_tracker.get_summary()


@router.get("/recent")
def get_recent(
    limit: int = Query(50, ge=1, le=500),
    user: UserContext = Depends(require_admin),
):
    return usage_tracker.get_recent(limit=limit)


@router.get("/user/{email}")
def get_user_usage(
    email: str,
    user: UserContext = Depends(require_admin),
):
    return usage_tracker.get_user_usage(email)
