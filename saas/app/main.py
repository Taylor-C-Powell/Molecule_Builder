"""MolBuilder SaaS API - FastAPI application."""

import logging
import secrets
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.responses import Response

from app.exceptions import register_exception_handlers
from app.middleware import UsageTrackingMiddleware
from app.middleware_versioning import VersioningMiddleware
from app.models.common import HealthResponse
from app.routers import auth, billing, legal, molecule, retrosynthesis, process, elements, analytics, audit, version, feasibility

logger = logging.getLogger("molbuilder.startup")

_LANDING_HTML = (Path(__file__).parent / "landing.html").read_text(encoding="utf-8")


class SecurityHeadersMiddleware(BaseHTTPMiddleware):
    """Add standard security headers to all responses."""

    async def dispatch(self, request: Request, call_next) -> Response:
        response = await call_next(request)
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["X-Frame-Options"] = "DENY"
        response.headers["X-XSS-Protection"] = "1; mode=block"
        response.headers["Referrer-Policy"] = "strict-origin-when-cross-origin"
        response.headers["Permissions-Policy"] = "geolocation=(), camera=(), microphone=()"
        if request.url.scheme == "https":
            response.headers["Strict-Transport-Security"] = "max-age=63072000; includeSubDomains"
        return response


@asynccontextmanager
async def lifespan(app: FastAPI):
    from app.config import settings as _cfg

    # --- Sentry initialization ---
    if _cfg.sentry_dsn:
        try:
            import sentry_sdk
            sentry_sdk.init(
                dsn=_cfg.sentry_dsn,
                traces_sample_rate=_cfg.sentry_traces_sample_rate,
                environment=_cfg.sentry_environment,
                release=f"molbuilder-api@{app.version}",
                send_default_pii=False,
            )
            logger.info("Sentry error monitoring initialized")
        except ImportError:
            logger.warning("sentry-sdk not installed, skipping error monitoring")

    # --- JWT secret validation ---
    import os
    is_production = os.environ.get("RAILWAY_ENVIRONMENT") or os.environ.get("PRODUCTION")

    if not _cfg.jwt_secret_key:
        if is_production:
            logger.critical(
                "JWT_SECRET_KEY is not set! Set a secure secret "
                "in production via environment variable JWT_SECRET_KEY"
            )
            raise SystemExit(1)
        else:
            # Generate ephemeral random secret for local dev
            _cfg.jwt_secret_key = f"dev-ephemeral-{secrets.token_hex(32)}"
            logger.warning(
                "JWT_SECRET_KEY not set - using ephemeral random secret (local dev only)"
            )

    # Generate ephemeral HMAC secrets if not configured (prod should set these)
    if not _cfg.api_key_hmac_secret:
        if is_production:
            logger.warning("API_KEY_HMAC_SECRET not set in production - generating ephemeral secret")
        _cfg.api_key_hmac_secret = secrets.token_hex(32)

    if not _cfg.audit_hmac_secret:
        if is_production:
            logger.warning("AUDIT_HMAC_SECRET not set in production - generating ephemeral secret")
        _cfg.audit_hmac_secret = secrets.token_hex(32)

    # Warm up molbuilder data structures
    from molbuilder.smiles.parser import parse
    parse("C")

    # Load persisted API keys from SQLite
    from app.services.user_db import get_user_db
    from app.auth.api_keys import api_key_store
    api_key_store.load_from_db(get_user_db())

    # Bootstrap admin account if configured and not already present
    from app.config import settings, Tier
    from app.auth.roles import Role
    if settings.admin_bootstrap_email:
        existing = [u for u in api_key_store.list_users()
                    if u["email"] == settings.admin_bootstrap_email]
        if not existing:
            raw_key = api_key_store.create(
                email=settings.admin_bootstrap_email,
                tier=Tier.ENTERPRISE,
                role=Role.ADMIN,
            )
            # SECURITY: Only log a masked version of the key
            masked = raw_key[:12] + "..." + raw_key[-4:]
            logger.warning(
                "ADMIN BOOTSTRAP: Created admin key for %s: %s "
                "(store this securely, it will not be shown again)",
                settings.admin_bootstrap_email, masked,
            )
            # Write full key to a file in /data if available, otherwise log masked only
            import os
            data_dir = os.path.dirname(settings.user_db_path)
            key_file = os.path.join(data_dir, ".admin_bootstrap_key")
            try:
                old_umask = os.umask(0o077)
                try:
                    with open(key_file, "w") as f:
                        f.write(raw_key)
                finally:
                    os.umask(old_umask)
                logger.info("Admin key written to %s (retrieve and delete this file)", key_file)
            except OSError:
                logger.warning(
                    "Could not write admin key to file. "
                    "Key prefix: %s - contact support if you need recovery.",
                    raw_key[:12],
                )
        else:
            logger.info(
                "Admin key for %s already exists, skipping bootstrap",
                settings.admin_bootstrap_email,
            )

    yield


app = FastAPI(
    title="MolBuilder API",
    description="REST API for molecular engineering: SMILES parsing, retrosynthesis, and process engineering",
    version="0.1.0",
    lifespan=lifespan,
)

# Security headers (outermost middleware, runs last)
app.add_middleware(SecurityHeadersMiddleware)
app.add_middleware(VersioningMiddleware)
app.add_middleware(UsageTrackingMiddleware)

# CORS configuration - never allow credentials with wildcard origins
from app.config import settings as _settings
_cors_origins = [o.strip() for o in _settings.cors_origins.split(",") if o.strip()]
_allow_credentials = bool(_cors_origins) and "*" not in _cors_origins

if "*" in [o.strip() for o in _settings.cors_origins.split(",") if o.strip()]:
    logger.warning(
        "CORS_ORIGINS contains '*' - credentials will be disabled. "
        "Set specific origins to enable credentials."
    )
    _cors_origins = ["*"]
    _allow_credentials = False

app.add_middleware(
    CORSMiddleware,
    allow_origins=_cors_origins if _cors_origins else [],
    allow_credentials=_allow_credentials,
    allow_methods=["GET", "POST", "PUT", "PATCH", "DELETE", "OPTIONS"],
    allow_headers=["Authorization", "X-API-Key", "X-API-Version", "Content-Type"],
)

register_exception_handlers(app)

app.include_router(auth.router)
app.include_router(billing.router)
app.include_router(legal.router)
app.include_router(molecule.router)
app.include_router(retrosynthesis.router)
app.include_router(process.router)
app.include_router(elements.router)
app.include_router(analytics.router)
app.include_router(audit.router)
app.include_router(version.router)
app.include_router(feasibility.router)


@app.get("/", response_class=HTMLResponse, include_in_schema=False)
async def landing_page():
    return _LANDING_HTML


@app.get("/health", response_model=HealthResponse, tags=["health"])
async def health():
    import molbuilder
    return HealthResponse(version=molbuilder.__version__)
