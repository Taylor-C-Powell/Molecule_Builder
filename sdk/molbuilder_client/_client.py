"""Synchronous MolBuilder API client backed by ``httpx.Client``."""

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


class MolBuilder:
    """Synchronous client for the MolBuilder REST API.

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
        self._client = httpx.Client(
            headers=build_headers(api_key),
            timeout=timeout,
        )

    # -- Context manager ------------------------------------------------------

    def __enter__(self) -> MolBuilder:
        return self

    def __exit__(self, *exc: Any) -> None:
        self.close()

    def close(self) -> None:
        """Release the underlying connection pool."""
        self._client.close()

    # -- Private helpers ------------------------------------------------------

    def _url(self, path: str) -> str:
        return build_url(self._base_url, path)

    def _get(self, path: str, *, params: dict[str, Any] | None = None) -> dict[str, Any]:
        resp = self._client.get(self._url(path), params=params)
        raise_for_status(resp)
        return resp.json()

    def _post(self, path: str, *, json: dict[str, Any] | None = None) -> dict[str, Any]:
        resp = self._client.post(self._url(path), json=json)
        raise_for_status(resp)
        return resp.json()

    # -- Auth -----------------------------------------------------------------

    def register(self, email: str) -> APIKeyInfo:
        """Register a new free-tier API key for *email*."""
        data = self._post("auth/register", json={"email": email})
        return from_dict(APIKeyInfo, data)

    def get_token(self) -> Token:
        """Exchange the current API key for a short-lived JWT."""
        data = self._post("auth/token", json={"api_key": self._api_key})
        return from_dict(Token, data)

    # -- Molecule -------------------------------------------------------------

    def from_smiles(self, smiles: str, *, name: str = "") -> MoleculeInfo:
        """Parse a SMILES string and return basic molecule info."""
        payload: dict[str, Any] = {"smiles": smiles}
        if name:
            payload["name"] = name
        data = self._post("molecule/from-smiles", json=payload)
        return from_dict(MoleculeInfo, data)

    def get_properties(self, mol_id: str) -> MoleculeProperties:
        """Retrieve computed molecular properties for *mol_id*."""
        data = self._get(f"molecule/{mol_id}/properties")
        return from_dict(MoleculeProperties, data)

    def get_3d(self, mol_id: str) -> Molecule3D:
        """Retrieve 3D atomic coordinates and bonds for *mol_id*."""
        data = self._get(f"molecule/{mol_id}/3d")
        return from_dict(Molecule3D, data)

    # -- Elements -------------------------------------------------------------

    def get_element(self, symbol: str) -> Element:
        """Look up an element by its chemical symbol (e.g. ``"Fe"``)."""
        data = self._get(f"elements/{symbol}")
        return from_dict(Element, data)

    # -- Retrosynthesis -------------------------------------------------------

    def retrosynthesis(
        self,
        smiles: str,
        *,
        max_depth: int = 5,
        beam_width: int = 5,
    ) -> RetrosynthesisPlan:
        """Plan retrosynthetic routes for a target SMILES."""
        data = self._post(
            "retrosynthesis/plan",
            json={
                "smiles": smiles,
                "max_depth": max_depth,
                "beam_width": beam_width,
            },
        )
        return from_dict(RetrosynthesisPlan, data)

    # -- Process Evaluation ---------------------------------------------------

    def process_evaluate(
        self,
        smiles: str,
        *,
        scale_kg: float = 1.0,
        max_depth: int = 5,
        beam_width: int = 5,
    ) -> ProcessEvaluation:
        """Evaluate process engineering for a synthesis route."""
        data = self._post(
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

    def create_checkout(self, plan: str) -> CheckoutSession:
        """Create a Stripe Checkout session for *plan* (e.g. ``"pro_monthly"``)."""
        data = self._post("billing/checkout", json={"plan": plan})
        return from_dict(CheckoutSession, data)

    def billing_status(self) -> BillingStatus:
        """Return the current billing / subscription status."""
        data = self._get("billing/status")
        return from_dict(BillingStatus, data)
