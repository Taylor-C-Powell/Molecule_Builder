"""Version info endpoint."""

from fastapi import APIRouter
from pydantic import BaseModel

from app.middleware_versioning import API_VERSION, MIN_SUPPORTED_VERSION

router = APIRouter(prefix="/api/v1", tags=["version"])


class VersionResponse(BaseModel):
    api_version: str
    min_supported_version: str
    deprecation_policy: str


@router.get("/version", response_model=VersionResponse)
async def get_version():
    return VersionResponse(
        api_version=API_VERSION,
        min_supported_version=MIN_SUPPORTED_VERSION,
        deprecation_policy=(
            "API versions are supported for a minimum of 12 months after "
            "a new version is released. Deprecated versions will return an "
            "X-Deprecation-Warning header. Clients should migrate to the "
            "latest version within the deprecation window."
        ),
    )
