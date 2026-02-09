"""Element lookup endpoint."""

from fastapi import APIRouter, Depends
from app.dependencies import UserContext, check_rate_limit
from app.exceptions import MolBuilderAPIError
from app.models.elements import ElementResponse
from molbuilder.core.elements import from_symbol

router = APIRouter(prefix="/api/v1/elements", tags=["elements"])


@router.get("/{symbol}", response_model=ElementResponse)
def get_element(
    symbol: str,
    user: UserContext = Depends(check_rate_limit),
):
    try:
        z, sym, name, weight = from_symbol(symbol)
    except ValueError:
        raise MolBuilderAPIError(404, f"Element '{symbol}' not found")
    return ElementResponse(
        atomic_number=z, symbol=sym, name=name, atomic_weight=weight
    )
