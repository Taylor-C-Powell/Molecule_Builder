"""Tests for SDK report download methods."""

from __future__ import annotations

import httpx
import respx

from molbuilder_client import MolBuilder


def test_download_report(client: MolBuilder, mock_api: respx.Router) -> None:
    pdf_bytes = b"%PDF-1.4 fake pdf content here"
    mock_api.post("/reports/process-pdf").mock(
        return_value=httpx.Response(
            200,
            content=pdf_bytes,
            headers={"Content-Type": "application/pdf"},
        )
    )

    result = client.download_report("CCO", scale_kg=1.0)

    assert isinstance(result, bytes)
    assert result == pdf_bytes


def test_download_report_save_to(client: MolBuilder, mock_api: respx.Router, tmp_path) -> None:
    pdf_bytes = b"%PDF-1.4 fake pdf content"
    mock_api.post("/reports/process-pdf").mock(
        return_value=httpx.Response(
            200,
            content=pdf_bytes,
            headers={"Content-Type": "application/pdf"},
        )
    )

    out_path = tmp_path / "report.pdf"
    result = client.download_report("CCO", save_to=str(out_path))

    assert result == pdf_bytes
    assert out_path.read_bytes() == pdf_bytes


def test_download_report_params(client: MolBuilder, mock_api: respx.Router) -> None:
    route = mock_api.post("/reports/process-pdf").mock(
        return_value=httpx.Response(
            200,
            content=b"%PDF-1.4",
            headers={"Content-Type": "application/pdf"},
        )
    )

    client.download_report("c1ccccc1", scale_kg=10.0)

    # Verify query params were sent
    request = route.calls.last.request
    assert "smiles=c1ccccc1" in str(request.url)
    assert "scale_kg=10" in str(request.url)
