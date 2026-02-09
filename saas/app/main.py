"""MolBuilder SaaS API - FastAPI application."""

from contextlib import asynccontextmanager
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.exceptions import register_exception_handlers
from app.middleware import UsageTrackingMiddleware
from app.middleware_versioning import VersioningMiddleware
from app.models.common import HealthResponse
from app.routers import auth, molecule, retrosynthesis, process, elements, analytics, audit, version


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Warm up molbuilder data structures by parsing a trivial molecule
    from molbuilder.smiles.parser import parse
    parse("C")
    yield


app = FastAPI(
    title="MolBuilder API",
    description="REST API for molecular engineering: SMILES parsing, retrosynthesis, and process engineering",
    version="0.1.0",
    lifespan=lifespan,
)

app.add_middleware(VersioningMiddleware)
app.add_middleware(UsageTrackingMiddleware)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

register_exception_handlers(app)

app.include_router(auth.router)
app.include_router(molecule.router)
app.include_router(retrosynthesis.router)
app.include_router(process.router)
app.include_router(elements.router)
app.include_router(analytics.router)
app.include_router(audit.router)
app.include_router(version.router)


@app.get("/health", response_model=HealthResponse, tags=["health"])
async def health():
    import molbuilder
    return HealthResponse(version=molbuilder.__version__)
