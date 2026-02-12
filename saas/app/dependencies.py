"""FastAPI dependencies for auth and rate limiting."""

import logging

import jwt as pyjwt
from fastapi import Depends, Header
from app.auth.api_keys import api_key_store
from app.auth.jwt_handler import decode_token
from app.auth.rate_limiter import rate_limiter
from app.auth.roles import Role
from app.config import Tier
from app.exceptions import AuthenticationError, AuthorizationError, RateLimitExceeded

logger = logging.getLogger("molbuilder.dependencies")


class UserContext:
    __slots__ = ("email", "tier", "role")

    def __init__(self, email: str, tier: Tier, role: Role = Role.CHEMIST):
        self.email = email
        self.tier = tier
        self.role = role


def get_current_user(
    x_api_key: str | None = Header(None),
    authorization: str | None = Header(None),
) -> UserContext:
    """Extract user from API key header or Bearer JWT."""
    # Try JWT first
    if authorization and authorization.startswith("Bearer "):
        token = authorization[7:]
        try:
            payload = decode_token(token)
            role_str = payload.get("role", Role.CHEMIST.value)
            try:
                role = Role(role_str)
            except ValueError:
                logger.warning("Invalid role %r in JWT, defaulting to CHEMIST", role_str)
                role = Role.CHEMIST
            return UserContext(
                email=payload["sub"],
                tier=Tier(payload["tier"]),
                role=role,
            )
        except (pyjwt.PyJWTError, KeyError, ValueError):
            raise AuthenticationError("Invalid or expired JWT token")

    # Fall back to API key
    if x_api_key:
        record = api_key_store.validate(x_api_key)
        if record is not None:
            return UserContext(email=record.email, tier=record.tier, role=record.role)
        raise AuthenticationError("Invalid API key")

    raise AuthenticationError("Missing authentication: provide X-API-Key or Authorization header")


def require_admin(user: UserContext = Depends(get_current_user)) -> UserContext:
    """Dependency that requires admin role."""
    if user.role != Role.ADMIN:
        raise AuthorizationError("Admin access required")
    return user


def require_role(*allowed_roles: Role):
    """Factory returning a dependency that requires one of the allowed roles."""
    def _check(user: UserContext = Depends(get_current_user)) -> UserContext:
        if user.role not in allowed_roles:
            raise AuthorizationError(
                f"Requires one of: {', '.join(r.value for r in allowed_roles)}"
            )
        return user
    return _check


def check_rate_limit(user: UserContext = Depends(get_current_user)) -> UserContext:
    """Enforce per-minute rate limit."""
    if not rate_limiter.check(user.email, user.tier):
        raise RateLimitExceeded()
    return user


def check_expensive_rate_limit(user: UserContext = Depends(get_current_user)) -> UserContext:
    """Enforce both per-minute and hourly expensive-endpoint rate limits."""
    if not rate_limiter.check(user.email, user.tier):
        raise RateLimitExceeded()
    if not rate_limiter.check_expensive(user.email, user.tier):
        raise RateLimitExceeded(retry_after=3600)
    return user
