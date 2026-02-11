"""Auth endpoints: register API key, exchange for JWT, user management."""

from fastapi import APIRouter, Depends
from app.auth.api_keys import api_key_store
from app.auth.jwt_handler import create_token
from app.auth.roles import Role
from app.config import Tier, settings
from app.dependencies import UserContext, require_admin
from app.exceptions import AuthenticationError
from app.models.auth import (
    APIKeyCreate, AdminKeyCreate, APIKeyResponse, TokenRequest, TokenResponse, UserInfo,
)

router = APIRouter(prefix="/api/v1/auth", tags=["auth"])


@router.post("/register", response_model=APIKeyResponse)
def register(body: APIKeyCreate):
    """Public registration: always creates a free-tier chemist key."""
    raw_key = api_key_store.create(email=body.email, tier=Tier.FREE, role=Role.CHEMIST)
    return APIKeyResponse(
        api_key=raw_key, email=body.email, tier=Tier.FREE, role=Role.CHEMIST,
    )


@router.post("/provision", response_model=APIKeyResponse)
def provision(body: AdminKeyCreate, admin: UserContext = Depends(require_admin)):
    """Admin-only: create an API key with any tier and role."""
    raw_key = api_key_store.create(email=body.email, tier=body.tier, role=body.role)
    return APIKeyResponse(
        api_key=raw_key, email=body.email, tier=body.tier, role=body.role,
    )


@router.post("/token", response_model=TokenResponse)
def get_token(body: TokenRequest):
    record = api_key_store.validate(body.api_key)
    if record is None:
        raise AuthenticationError("Invalid API key")
    token = create_token(email=record.email, tier=record.tier, role=record.role)
    return TokenResponse(
        access_token=token,
        expires_in=settings.jwt_expiry_minutes * 60,
    )


@router.get("/users", response_model=list[UserInfo])
def list_users(admin: UserContext = Depends(require_admin)):
    return [UserInfo(**u) for u in api_key_store.list_users()]


@router.delete("/users/{email}")
def delete_user(email: str, admin: UserContext = Depends(require_admin)):
    count = api_key_store.delete_by_email(email)
    return {"deleted": count, "email": email}
