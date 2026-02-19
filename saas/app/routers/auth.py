"""Auth endpoints: register API key, exchange for JWT, user management."""

import re

from fastapi import APIRouter, Depends, Request
from app.auth.api_keys import api_key_store
from app.auth.jwt_handler import create_token
from app.auth.rate_limiter import rate_limiter
from app.auth.roles import Role
from app.config import Tier, settings
from app.dependencies import UserContext, require_admin
from app.exceptions import AuthenticationError, RateLimitExceeded, MolBuilderAPIError
from app.models.auth import (
    APIKeyCreate, AdminKeyCreate, APIKeyResponse, TokenRequest, TokenResponse, UserInfo,
)

router = APIRouter(prefix="/api/v1/auth", tags=["auth"])

# RFC 5322 simplified email regex
_EMAIL_RE = re.compile(r"^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$")
_MAX_EMAIL_LEN = 254


def _validate_email(email: str) -> str:
    """Validate and normalize email address."""
    email = email.strip().lower()
    if not email or len(email) > _MAX_EMAIL_LEN:
        raise MolBuilderAPIError(422, "Invalid email address")
    if not _EMAIL_RE.match(email):
        raise MolBuilderAPIError(422, "Invalid email format")
    return email


def _get_client_ip(request: Request) -> str:
    """Extract client IP, respecting X-Forwarded-For behind proxy."""
    forwarded = request.headers.get("x-forwarded-for")
    if forwarded:
        return forwarded.split(",")[0].strip()
    return request.client.host if request.client else "unknown"


@router.post("/register", response_model=APIKeyResponse)
def register(body: APIKeyCreate, request: Request):
    """Public registration: always creates a free-tier chemist key."""
    # IP-based rate limit for unauthenticated registration
    ip = _get_client_ip(request)
    if not rate_limiter.check_ip(ip, "register", settings.register_rpm):
        raise RateLimitExceeded()

    email = _validate_email(body.email)
    raw_key = api_key_store.create(email=email, tier=Tier.FREE, role=Role.CHEMIST)
    return APIKeyResponse(
        api_key=raw_key, email=email, tier=Tier.FREE, role=Role.CHEMIST,
    )


@router.post("/provision", response_model=APIKeyResponse)
def provision(body: AdminKeyCreate, admin: UserContext = Depends(require_admin)):
    """Admin-only: create an API key with any tier and role."""
    email = _validate_email(body.email)
    raw_key = api_key_store.create(email=email, tier=body.tier, role=body.role)
    return APIKeyResponse(
        api_key=raw_key, email=email, tier=body.tier, role=body.role,
    )


@router.post("/token", response_model=TokenResponse)
def get_token(body: TokenRequest, request: Request):
    # IP-based rate limit for token exchange
    ip = _get_client_ip(request)
    if not rate_limiter.check_ip(ip, "token", settings.token_rpm):
        raise RateLimitExceeded()

    record = api_key_store.validate(body.api_key)
    if record is None:
        raise AuthenticationError("Invalid API key")
    token = create_token(email=record.email, tier=record.tier, role=record.role)
    return TokenResponse(
        access_token=token,
        expires_in=settings.jwt_expiry_minutes * 60,
        email=record.email,
        tier=record.tier.value,
    )


@router.post("/rotate", response_model=APIKeyResponse)
def rotate_key(body: TokenRequest):
    """Exchange a valid API key for a new one. The old key is revoked."""
    record = api_key_store.validate(body.api_key)
    if record is None:
        raise AuthenticationError("Invalid API key")

    # Revoke old key and create new one with same tier/role
    api_key_store.revoke_key(body.api_key)
    new_key = api_key_store.create(
        email=record.email, tier=record.tier, role=record.role,
    )
    return APIKeyResponse(
        api_key=new_key, email=record.email, tier=record.tier, role=record.role,
    )


@router.get("/users", response_model=list[UserInfo])
def list_users(admin: UserContext = Depends(require_admin)):
    return [UserInfo(**u) for u in api_key_store.list_users()]


@router.patch("/users/{email}/tier")
def update_user_tier(email: str, body: AdminKeyCreate, admin: UserContext = Depends(require_admin)):
    """Admin-only: update tier for all keys belonging to an email."""
    count = api_key_store.update_tier(email, body.tier)
    if count == 0:
        from fastapi import HTTPException
        raise HTTPException(status_code=404, detail=f"No keys found for {email}")
    return {"updated": count, "email": email, "tier": body.tier.value}


@router.delete("/users/{email}")
def delete_user(email: str, admin: UserContext = Depends(require_admin)):
    count = api_key_store.delete_by_email(email)
    return {"deleted": count, "email": email}
