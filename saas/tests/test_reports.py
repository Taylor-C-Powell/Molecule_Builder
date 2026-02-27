"""Tests for PDF report endpoints."""


def test_process_pdf_returns_pdf(client, auth_headers):
    """POST /reports/process-pdf returns a PDF for a valid SMILES."""
    resp = client.post(
        "/api/v1/reports/process-pdf",
        params={"smiles": "CCO", "scale_kg": 1.0},
        headers=auth_headers,
    )
    assert resp.status_code == 200
    assert resp.headers["content-type"] == "application/pdf"
    assert resp.content[:5] == b"%PDF-"


def test_process_pdf_requires_auth(client):
    """PDF endpoint requires authentication."""
    resp = client.post(
        "/api/v1/reports/process-pdf",
        params={"smiles": "CCO"},
    )
    assert resp.status_code == 401


def test_process_pdf_invalid_smiles(client, auth_headers):
    """Invalid SMILES returns 422."""
    resp = client.post(
        "/api/v1/reports/process-pdf",
        params={"smiles": "INVALID!!!"},
        headers=auth_headers,
    )
    assert resp.status_code == 422
