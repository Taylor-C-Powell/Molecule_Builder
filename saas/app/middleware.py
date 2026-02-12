"""Usage tracking and audit middleware."""

import hashlib
import json
import logging
import time
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.requests import Request
from starlette.responses import Response
import jwt as pyjwt
from app.auth.api_keys import api_key_store
from app.auth.jwt_handler import decode_token
from app.auth.roles import Role, check_permission
from app.services.usage_tracker import usage_tracker, RequestRecord

logger = logging.getLogger("molbuilder.middleware")

# Paths exempt from RBAC checks (public, health, auth)
_RBAC_EXEMPT_PREFIXES = (
    "/",
    "/health",
    "/docs",
    "/openapi.json",
    "/api/v1/auth/register",
    "/api/v1/auth/token",
    "/api/v1/billing/webhook",
)


def _hash_ip(ip: str) -> str:
    """One-way hash of IP address for privacy-preserving audit trail."""
    return hashlib.sha256(ip.encode()).hexdigest()[:16]


def _extract_user(request: Request) -> tuple[str, str, str]:
    """Best-effort extraction of user email, tier, and role from request headers."""
    # Try JWT
    auth = request.headers.get("authorization", "")
    if auth.startswith("Bearer "):
        try:
            payload = decode_token(auth[7:])
            return (
                payload.get("sub", ""),
                payload.get("tier", ""),
                payload.get("role", "chemist"),
            )
        except pyjwt.PyJWTError:
            pass

    # Try API key
    api_key = request.headers.get("x-api-key", "")
    if api_key:
        record = api_key_store.validate(api_key)
        if record is not None:
            return record.email, record.tier.value, record.role.value

    return "", "", ""


def _extract_body_summary(body: bytes, path: str) -> dict:
    """Extract key fields from request body for tracking (no PII)."""
    if not body:
        return {}
    try:
        data = json.loads(body)
    except (json.JSONDecodeError, UnicodeDecodeError):
        return {}

    summary = {}
    if "smiles" in data:
        summary["smiles"] = data["smiles"]
    if "scale_kg" in data:
        summary["scale_kg"] = data["scale_kg"]
    # Deliberately omit email and other PII from tracking summaries
    return summary


def _is_rbac_exempt(path: str) -> bool:
    """Check if a path is exempt from RBAC enforcement."""
    for prefix in _RBAC_EXEMPT_PREFIXES:
        if path == prefix or (prefix != "/" and path.startswith(prefix)):
            return True
    return False


class UsageTrackingMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request: Request, call_next) -> Response:
        start = time.monotonic()

        # Read body for POST requests (need to cache it for downstream)
        body = b""
        if request.method == "POST":
            body = await request.body()

        email, tier, role = _extract_user(request)
        body_summary = _extract_body_summary(body, request.url.path)

        # RBAC enforcement for authenticated requests
        if email and role and not _is_rbac_exempt(request.url.path):
            try:
                role_enum = Role(role)
            except ValueError:
                role_enum = Role.CHEMIST
            if not check_permission(role_enum, request.method, request.url.path):
                from starlette.responses import JSONResponse
                return JSONResponse(
                    status_code=403,
                    content={"error": "Insufficient permissions for this resource"},
                )

        response = await call_next(request)

        latency_ms = (time.monotonic() - start) * 1000

        # Normalize path: replace molecule IDs with {id}
        path = request.url.path
        parts = path.split("/")
        if len(parts) >= 5 and parts[3] == "molecule" and len(parts[4]) == 12:
            parts[4] = "{id}"
            path = "/".join(parts)

        usage_tracker.record(RequestRecord(
            timestamp=time.time(),
            method=request.method,
            path=path,
            status_code=response.status_code,
            latency_ms=latency_ms,
            user_email=email,
            user_tier=tier,
            body_summary=body_summary,
        ))

        # Audit trail recording
        if email:
            try:
                from app.services.audit_db import get_audit_db
                db = get_audit_db()
                action = f"{request.method} {path}"
                input_str = json.dumps(body_summary) if body_summary else ""
                # Hash IP for privacy
                raw_ip = request.client.host if request.client else ""
                ip_hash = _hash_ip(raw_ip) if raw_ip else ""
                db.record(
                    user_email=email,
                    user_role=role,
                    user_tier=tier,
                    action=action,
                    input_summary=input_str,
                    output_hash="",
                    status_code=response.status_code,
                    latency_ms=latency_ms,
                    ip_address=ip_hash,
                )
            except Exception:
                logger.warning("Audit trail write failed", exc_info=True)

        return response
