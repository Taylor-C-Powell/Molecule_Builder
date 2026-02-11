"""Asynchronous MolBuilder API client backed by ``httpx.AsyncClient``."""

from __future__ import annotations

from typing import Any

import httpx

from molbuilder_client._base import (
    DEFAULT_BASE_URL,
    build_headers,
    build_url,
    from_dict,
    raise_for_status,
)
from molbuilder_client._models import (
    APIKeyInfo,
    BillingStatus,
    CheckoutSession,
    Element,
    Molecule3D,
    MoleculeInfo,
    MoleculeProperties,
    ProcessEvaluation,
    RetrosynthesisPlan,
    Token,
)


class AsyncMolBuilder:
    """Asynchronous client for the MolBuilder REST API.

    Parameters
    ----------
    api_key:
        Your MolBuilder API key (starts with ``mb_``).
    base_url:
        Override the default production URL (useful for local dev).
    timeout:
        Request timeout in seconds.  Defaults to 30.
    """

    def __init__(
        self,
        api_key: str,
        *,
        base_url: str = DEFAULT_BASE_URL,
        timeout: float = 30.0,
    ) -> None:
        self._api_key = api_key
        self._base_url = base_url
        self._client = httpx.AsyncClient(
            headers=build_headers(api_key),
            timeout=timeout,
        )

    # -- Async context manager ------------------------------------------------

    async def __aenter__(self) -> AsyncMolBuilder:
        return self

    async def __aexit__(self, *exc: Any) -> None:
        await self.close()

    async def close(self) -> None:
        """Release the underlying connection pool."""
        await self._client.aclose()

    # -- Private helpers ------------------------------------------------------

    def _url(self, path: str) -> str:
        return build_url(self._base_url, path)

    async def _get(self, path: str, *, params: dict[str, Any] | None = None) -> dict[str, Any]:
        resp = await self._client.get(self._url(path), params=params)
        raise_for_status(resp)
        return resp.json()

    async def _post(self, path: str, *, json: dict[str, Any] | None = None) -> dict[str, Any]:
        resp = await self._client.post(self._url(path), json=json)
        raise_for_status(resp)
        return resp.json()

    # -- Auth -----------------------------------------------------------------

    async def register(self, email: str) -> APIKeyInfo:
        """Register a new free-tier API key for *email*."""
        data = await self._post("auth/register", json={"email": email})
        return from_dict(APIKeyInfo, data)

    async def get_token(self) -> Token:
        """Exchange the current API key for a short-lived JWT."""
        data = await self._post("auth/token", json={"api_key": self._api_key})
        return from_dict(Token, data)

    # -- Molecule -------------------------------------------------------------

    async def from_smiles(self, smiles: str, *, name: str = "") -> MoleculeInfo:
        """Parse a SMILES string and return basic molecule info."""
        payload: dict[str, Any] = {"smiles": smiles}
        if name:
            payload["name"] = name
        data = await self._post("molecule/from-smiles", json=payload)
        return from_dict(MoleculeInfo, data)

    async def get_properties(self, mol_id: str) -> MoleculeProperties:
        """Retrieve computed molecular properties for *mol_id*."""
        data = await self._get(f"molecule/{mol_id}/properties")
        return from_dict(MoleculeProperties, data)

    async def get_3d(self, mol_id: str) -> Molecule3D:
        """Retrieve 3D atomic coordinates and bonds for *mol_id*."""
        data = await self._get(f"molecule/{mol_id}/3d")
        return from_dict(Molecule3D, data)

    # -- Elements -------------------------------------------------------------

    async def get_element(self, symbol: str) -> Element:
        """Look up an element by its chemical symbol (e.g. ``"Fe"``)."""
        data = await self._get(f"elements/{symbol}")
        return from_dict(Element, data)

    # -- Retrosynthesis -------------------------------------------------------

    async def retrosynthesis(
        self,
        smiles: str,
        *,
        max_depth: int = 5,
        beam_width: int = 5,
    ) -> RetrosynthesisPlan:
        """Plan retrosynthetic routes for a target SMILES."""
        data = await self._post(
            "retrosynthesis/plan",
            json={
                "smiles": smiles,
                "max_depth": max_depth,
                "beam_width": beam_width,
            },
        )
        return from_dict(RetrosynthesisPlan, data)

    # -- Process Evaluation ---------------------------------------------------

    async def process_evaluate(
        self,
        smiles: str,
        *,
        scale_kg: float = 1.0,
        max_depth: int = 5,
        beam_width: int = 5,
    ) -> ProcessEvaluation:
        """Evaluate process engineering for a synthesis route."""
        data = await self._post(
            "process/evaluate",
            json={
                "smiles": smiles,
                "scale_kg": scale_kg,
                "max_depth": max_depth,
                "beam_width": beam_width,
            },
        )
        return from_dict(ProcessEvaluation, data)

    # -- Billing --------------------------------------------------------------

    async def create_checkout(self, plan: str) -> CheckoutSession:
        """Create a Stripe Checkout session for *plan* (e.g. ``"pro_monthly"``)."""
        data = await self._post("billing/checkout", json={"plan": plan})
        return from_dict(CheckoutSession, data)

    async def billing_status(self) -> BillingStatus:
        """Return the current billing / subscription status."""
        data = await self._get("billing/status")
        return from_dict(BillingStatus, data)
