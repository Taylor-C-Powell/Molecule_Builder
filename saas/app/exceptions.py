"""Custom exceptions and global error handlers."""

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse


class MolBuilderAPIError(Exception):
    def __init__(self, status_code: int, detail: str):
        self.status_code = status_code
        self.detail = detail


class MoleculeNotFound(MolBuilderAPIError):
    def __init__(self, molecule_id: str):
        super().__init__(404, f"Molecule '{molecule_id}' not found")


class InvalidSMILES(MolBuilderAPIError):
    def __init__(self, smiles: str, reason: str = ""):
        msg = f"Invalid SMILES: '{smiles}'"
        if reason:
            msg += f" - {reason}"
        super().__init__(422, msg)


class RateLimitExceeded(MolBuilderAPIError):
    def __init__(self, retry_after: int = 60):
        self.retry_after = retry_after
        super().__init__(429, "Rate limit exceeded")


class AuthenticationError(MolBuilderAPIError):
    def __init__(self, detail: str = "Invalid or missing credentials"):
        super().__init__(401, detail)


class AuthorizationError(MolBuilderAPIError):
    def __init__(self, detail: str = "Insufficient permissions"):
        super().__init__(403, detail)


class InvalidFileFormat(MolBuilderAPIError):
    def __init__(self, filename_or_format: str):
        super().__init__(422, f"Unsupported or invalid file format: '{filename_or_format}'")


def register_exception_handlers(app: FastAPI) -> None:
    @app.exception_handler(MolBuilderAPIError)
    async def molbuilder_error_handler(
        request: Request, exc: MolBuilderAPIError
    ) -> JSONResponse:
        headers = {}
        if isinstance(exc, RateLimitExceeded):
            headers["Retry-After"] = str(exc.retry_after)
        return JSONResponse(
            status_code=exc.status_code,
            content={"error": exc.detail},
            headers=headers,
        )

    @app.exception_handler(Exception)
    async def general_error_handler(
        request: Request, exc: Exception
    ) -> JSONResponse:
        return JSONResponse(
            status_code=500,
            content={"error": "Internal server error"},
        )
