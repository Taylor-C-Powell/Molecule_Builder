"""Common response models."""

from pydantic import BaseModel


class ErrorResponse(BaseModel):
    error: str


class HealthResponse(BaseModel):
    status: str = "ok"
    version: str
