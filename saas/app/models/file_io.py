"""File I/O request/response models."""

from pydantic import BaseModel
from app.models.molecule import MoleculeResponse


class FileImportResponse(BaseModel):
    molecules: list[MoleculeResponse]
    format: str
    count: int
