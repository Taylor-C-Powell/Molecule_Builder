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
    ADMETProfile,
    APIKeyInfo,
    BatchList,
    BatchStatus,
    BatchSubmit,
    BillingStatus,
    CheckoutSession,
    Element,
    FileImportResult,
    LibraryImport,
    LibraryList,
    LibraryMolecule,
    Molecule3D,
    MoleculeInfo,
    MoleculeProperties,
    ProcessEvaluation,
    RetrosynthesisPlan,
    SolubilityResult,
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

    async def _put(self, path: str, *, json: dict[str, Any] | None = None) -> dict[str, Any]:
        resp = await self._client.put(self._url(path), json=json)
        raise_for_status(resp)
        return resp.json()

    async def _delete(self, path: str) -> dict[str, Any]:
        resp = await self._client.delete(self._url(path))
        raise_for_status(resp)
        return resp.json()

    async def _post_file(self, path: str, filename: str, content: bytes) -> dict[str, Any]:
        resp = await self._client.post(
            self._url(path),
            files={"file": (filename, content, "application/octet-stream")},
        )
        raise_for_status(resp)
        return resp.json()

    async def _get_raw(self, path: str, *, params: dict[str, Any] | None = None) -> httpx.Response:
        resp = await self._client.get(self._url(path), params=params)
        raise_for_status(resp)
        return resp

    async def _post_raw(self, path: str, *, params: dict[str, Any] | None = None) -> httpx.Response:
        resp = await self._client.post(self._url(path), params=params)
        raise_for_status(resp)
        return resp

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

    async def get_admet(self, mol_id: str) -> ADMETProfile:
        """Retrieve ADMET predictions for *mol_id*."""
        data = await self._get(f"molecule/{mol_id}/admet")
        return from_dict(ADMETProfile, data)

    async def get_solubility(self, mol_id: str) -> SolubilityResult:
        """Retrieve solubility and crystallization predictions for *mol_id*."""
        data = await self._get(f"molecule/{mol_id}/solubility")
        return from_dict(SolubilityResult, data)

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

    # -- Library ---------------------------------------------------------------

    async def library_save(
        self,
        smiles: str,
        *,
        name: str | None = None,
        tags: list[str] | None = None,
        notes: str | None = None,
    ) -> LibraryMolecule:
        """Save a molecule to the personal library."""
        payload: dict[str, Any] = {"smiles": smiles}
        if name is not None:
            payload["name"] = name
        if tags is not None:
            payload["tags"] = tags
        if notes is not None:
            payload["notes"] = notes
        data = await self._post("library/", json=payload)
        return from_dict(LibraryMolecule, data)

    async def library_get(self, mol_id: int) -> LibraryMolecule:
        """Retrieve a single library molecule by ID."""
        data = await self._get(f"library/{mol_id}")
        return from_dict(LibraryMolecule, data)

    async def library_list(
        self,
        *,
        tag: str | None = None,
        search: str | None = None,
        page: int = 1,
        per_page: int = 20,
    ) -> LibraryList:
        """List molecules in the personal library."""
        params: dict[str, Any] = {"page": page, "per_page": per_page}
        if tag is not None:
            params["tag"] = tag
        if search is not None:
            params["search"] = search
        data = await self._get("library/", params=params)
        return from_dict(LibraryList, data)

    async def library_update(
        self,
        mol_id: int,
        *,
        name: str | None = None,
        tags: list[str] | None = None,
        notes: str | None = None,
    ) -> LibraryMolecule:
        """Update metadata on a library molecule."""
        payload: dict[str, Any] = {}
        if name is not None:
            payload["name"] = name
        if tags is not None:
            payload["tags"] = tags
        if notes is not None:
            payload["notes"] = notes
        data = await self._put(f"library/{mol_id}", json=payload)
        return from_dict(LibraryMolecule, data)

    async def library_delete(self, mol_id: int) -> None:
        """Delete a molecule from the library."""
        await self._delete(f"library/{mol_id}")

    async def library_import(
        self,
        smiles_list: list[str],
        *,
        tag: str | None = None,
    ) -> LibraryImport:
        """Bulk-import SMILES into the library."""
        payload: dict[str, Any] = {"smiles_list": smiles_list}
        if tag is not None:
            payload["tag"] = tag
        data = await self._post("library/import", json=payload)
        return from_dict(LibraryImport, data)

    # -- Batch -----------------------------------------------------------------

    async def batch_submit(
        self,
        smiles_list: list[str],
        job_type: str,
        *,
        params: dict[str, Any] | None = None,
    ) -> BatchSubmit:
        """Submit a batch processing job."""
        payload: dict[str, Any] = {
            "smiles_list": smiles_list,
            "job_type": job_type,
        }
        if params is not None:
            payload["params"] = params
        data = await self._post("batch/submit", json=payload)
        return from_dict(BatchSubmit, data)

    async def batch_status(self, job_id: str) -> BatchStatus:
        """Check the status of a batch job."""
        data = await self._get(f"batch/{job_id}")
        return from_dict(BatchStatus, data)

    async def batch_list(self, *, page: int = 1, per_page: int = 20) -> BatchList:
        """List batch jobs."""
        data = await self._get("batch/", params={"page": page, "per_page": per_page})
        return from_dict(BatchList, data)

    async def batch_cancel(self, job_id: str) -> bool:
        """Cancel a running batch job. Returns True on success."""
        await self._delete(f"batch/{job_id}")
        return True

    # -- File I/O -------------------------------------------------------------

    async def import_file(self, file_path: str, *, format: str | None = None) -> FileImportResult:
        """Import a molecule file (XYZ, MOL, SDF, PDB, JSON)."""
        import pathlib
        p = pathlib.Path(file_path)
        content = p.read_bytes()
        filename = p.name
        if format:
            filename = f"{p.stem}.{format}"
        data = await self._post_file("molecule/import-file", filename, content)
        return from_dict(FileImportResult, data)

    async def export_file(self, mol_id: str, format: str, *, save_to: str | None = None) -> str:
        """Export a molecule in the specified format."""
        resp = await self._get_raw(f"molecule/{mol_id}/export/{format}")
        content = resp.text
        if save_to:
            import pathlib
            pathlib.Path(save_to).write_text(content, encoding="utf-8")
            return save_to
        return content

    # -- Reports --------------------------------------------------------------

    async def download_report(
        self,
        smiles: str,
        *,
        scale_kg: float = 1.0,
        save_to: str | None = None,
    ) -> bytes:
        """Download a process engineering PDF report."""
        resp = await self._post_raw(
            "reports/process-pdf",
            params={"smiles": smiles, "scale_kg": scale_kg},
        )
        pdf_bytes = resp.content
        if save_to:
            import pathlib
            pathlib.Path(save_to).write_bytes(pdf_bytes)
        return pdf_bytes
