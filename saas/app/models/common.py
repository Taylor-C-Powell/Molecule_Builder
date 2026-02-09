"""Common response models."""

from pydantic import BaseModel
from app.config import Tier


class ErrorResponse(BaseModel):
    error: str


class HealthResponse(BaseModel):
    status: str = "ok"
    version: str
