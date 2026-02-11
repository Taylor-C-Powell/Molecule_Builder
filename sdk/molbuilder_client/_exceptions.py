"""Exception hierarchy for MolBuilder API errors.

Each exception maps 1:1 to an HTTP status code returned by the server.
All exceptions carry the raw status_code and server message for debugging.
"""

from __future__ import annotations


class MolBuilderError(Exception):
    """Base exception for all MolBuilder SDK errors."""

    def __init__(self, message: str, status_code: int = 0) -> None:
        self.message = message
        self.status_code = status_code
        super().__init__(message)


# -- 4xx Client Errors -------------------------------------------------------

class AuthenticationError(MolBuilderError):
    """401 - Invalid or missing API key / token."""

    def __init__(self, message: str = "Invalid or missing credentials") -> None:
        super().__init__(message, status_code=401)


class ForbiddenError(MolBuilderError):
    """403 - Insufficient permissions (e.g. viewer calling POST)."""

    def __init__(self, message: str = "Insufficient permissions") -> None:
        super().__init__(message, status_code=403)


class NotFoundError(MolBuilderError):
    """404 - Resource not found (molecule ID, element symbol, etc.)."""

    def __init__(self, message: str = "Resource not found") -> None:
        super().__init__(message, status_code=404)


class ValidationError(MolBuilderError):
    """422 - Request validation failed (invalid SMILES, bad parameters)."""

    def __init__(self, message: str = "Validation error") -> None:
        super().__init__(message, status_code=422)


class RateLimitError(MolBuilderError):
    """429 - Rate limit exceeded.  Check ``retry_after`` for seconds to wait."""

    def __init__(
        self,
        message: str = "Rate limit exceeded",
        *,
        retry_after: int | None = None,
    ) -> None:
        self.retry_after = retry_after
        super().__init__(message, status_code=429)


# -- 5xx Server Errors -------------------------------------------------------

class ServerError(MolBuilderError):
    """500 - Unexpected server-side failure."""

    def __init__(self, message: str = "Internal server error") -> None:
        super().__init__(message, status_code=500)


class ServiceUnavailableError(MolBuilderError):
    """501/503 - Feature not configured or service temporarily unavailable."""

    def __init__(self, message: str = "Service unavailable") -> None:
        super().__init__(message, status_code=503)
