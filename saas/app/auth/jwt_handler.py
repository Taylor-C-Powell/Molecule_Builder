"""JWT creation, verification, and revocation."""

import threading
import time
from datetime import datetime, timedelta, timezone

import jwt
from app.config import settings, Tier
from app.auth.roles import Role


class TokenBlocklist:
    """In-memory JWT blocklist with automatic expiry cleanup."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        # jti -> expiry timestamp
        self._blocked: dict[str, float] = {}
        self._last_cleanup = time.monotonic()

    def block(self, jti: str, exp: float) -> None:
        """Block a token by its JTI until its expiry time."""
        with self._lock:
            self._blocked[jti] = exp
            self._cleanup_if_needed()

    def is_blocked(self, jti: str) -> bool:
        with self._lock:
            self._cleanup_if_needed()
            return jti in self._blocked

    def _cleanup_if_needed(self) -> None:
        """Remove expired entries every 5 minutes."""
        now = time.monotonic()
        if now - self._last_cleanup < 300:
            return
        self._last_cleanup = now
        current_time = time.time()
        self._blocked = {
            jti: exp for jti, exp in self._blocked.items()
            if exp > current_time
        }


token_blocklist = TokenBlocklist()


def create_token(email: str, tier: Tier, role: Role = Role.CHEMIST) -> str:
    import secrets
    now = datetime.now(timezone.utc)
    payload = {
        "sub": email,
        "tier": tier.value,
        "role": role.value,
        "iat": now,
        "exp": now + timedelta(minutes=settings.jwt_expiry_minutes),
        "jti": secrets.token_hex(16),
    }
    return jwt.encode(payload, settings.jwt_secret_key, algorithm=settings.jwt_algorithm)


def decode_token(token: str) -> dict:
    """Decode and validate a JWT. Raises jwt.PyJWTError on failure."""
    payload = jwt.decode(
        token,
        settings.jwt_secret_key,
        algorithms=[settings.jwt_algorithm],
    )
    # Check blocklist
    jti = payload.get("jti")
    if jti and token_blocklist.is_blocked(jti):
        raise jwt.InvalidTokenError("Token has been revoked")
    return payload


def revoke_token(token: str) -> bool:
    """Revoke a JWT by adding its JTI to the blocklist."""
    try:
        payload = jwt.decode(
            token,
            settings.jwt_secret_key,
            algorithms=[settings.jwt_algorithm],
        )
        jti = payload.get("jti")
        exp = payload.get("exp", 0)
        if jti:
            token_blocklist.block(jti, exp)
            return True
    except jwt.PyJWTError:
        pass
    return False
