"""Element lookup response models."""

from pydantic import BaseModel


class ElementResponse(BaseModel):
    atomic_number: int
    symbol: str
    name: str
    atomic_weight: float
