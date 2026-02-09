"""JWT creation and verification."""

from datetime import datetime, timedelta, timezone
import jwt
from app.config import settings, Tier
from app.auth.roles import Role


def create_token(email: str, tier: Tier, role: Role = Role.CHEMIST) -> str:
    now = datetime.now(timezone.utc)
    payload = {
        "sub": email,
        "tier": tier.value,
        "role": role.value,
        "iat": now,
        "exp": now + timedelta(minutes=settings.jwt_expiry_minutes),
    }
    return jwt.encode(payload, settings.jwt_secret_key, algorithm=settings.jwt_algorithm)


def decode_token(token: str) -> dict:
    """Decode and validate a JWT. Raises jwt.PyJWTError on failure."""
    return jwt.decode(
        token,
        settings.jwt_secret_key,
        algorithms=[settings.jwt_algorithm],
    )
