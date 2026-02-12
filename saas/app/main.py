"""MolBuilder SaaS API - FastAPI application."""

from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse

from app.exceptions import register_exception_handlers
from app.middleware import UsageTrackingMiddleware
from app.middleware_versioning import VersioningMiddleware
from app.models.common import HealthResponse
from app.routers import auth, billing, molecule, retrosynthesis, process, elements, analytics, audit, version

_LANDING_HTML = (Path(__file__).parent / "landing.html").read_text(encoding="utf-8")


@asynccontextmanager
async def lifespan(app: FastAPI):
    import logging
    logger = logging.getLogger("molbuilder.bootstrap")

    # Validate JWT secret
    from app.config import settings as _cfg
    if _cfg.jwt_secret_key == "dev-secret-change-in-production":
        import os
        if os.environ.get("RAILWAY_ENVIRONMENT") or os.environ.get("PRODUCTION"):
            logger.critical(
                "JWT_SECRET_KEY is the default value! Set a secure secret "
                "in production via environment variable JWT_SECRET_KEY"
            )
            raise SystemExit(1)
        else:
            logger.warning(
                "JWT_SECRET_KEY is the default value - acceptable for local dev only"
            )

    # Warm up molbuilder data structures by parsing a trivial molecule
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
            logger.warning(
                "ADMIN BOOTSTRAP: Created admin key for %s: %s "
                "(store this securely, it will not be shown again)",
                settings.admin_bootstrap_email, raw_key,
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

app.add_middleware(VersioningMiddleware)
app.add_middleware(UsageTrackingMiddleware)
from app.config import settings as _settings
_cors_origins = [o.strip() for o in _settings.cors_origins.split(",")]
app.add_middleware(
    CORSMiddleware,
    allow_origins=_cors_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

register_exception_handlers(app)

app.include_router(auth.router)
app.include_router(billing.router)
app.include_router(molecule.router)
app.include_router(retrosynthesis.router)
app.include_router(process.router)
app.include_router(elements.router)
app.include_router(analytics.router)
app.include_router(audit.router)
app.include_router(version.router)


@app.get("/", response_class=HTMLResponse, include_in_schema=False)
async def landing_page():
    return _LANDING_HTML


@app.get("/health", response_model=HealthResponse, tags=["health"])
async def health():
    import molbuilder
    return HealthResponse(version=molbuilder.__version__)
