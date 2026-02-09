"""API versioning middleware — adds version headers to all responses."""

from starlette.middleware.base import BaseHTTPMiddleware
from starlette.requests import Request
from starlette.responses import Response

API_VERSION = "1.0.0"
MIN_SUPPORTED_VERSION = "1.0.0"


class VersioningMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request: Request, call_next) -> Response:
        response = await call_next(request)

        response.headers["X-API-Version"] = API_VERSION
        response.headers["X-Min-Supported-Version"] = MIN_SUPPORTED_VERSION

        # If client sends an outdated version header, warn about deprecation
        client_version = request.headers.get("X-API-Version", "")
        if client_version and client_version < MIN_SUPPORTED_VERSION:
            response.headers["X-Deprecation-Warning"] = (
                f"API version {client_version} is deprecated. "
                f"Minimum supported version is {MIN_SUPPORTED_VERSION}. "
                "Please upgrade within 12 months."
            )

        return response
