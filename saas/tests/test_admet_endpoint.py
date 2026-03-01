"""Tests for the /molecule/{id}/admet endpoint."""

import pytest


def _create_molecule(client, auth_headers, smiles="CCO"):
    resp = client.post(
        "/api/v1/molecule/from-smiles",
        json={"smiles": smiles},
        headers=auth_headers,
    )
    assert resp.status_code == 200
    return resp.json()["id"]


class TestADMETEndpoint:
    def test_requires_auth(self, client):
        resp = client.get("/api/v1/molecule/fake-id/admet")
        assert resp.status_code == 401

    def test_not_found(self, client, auth_headers):
        resp = client.get(
            "/api/v1/molecule/nonexistent/admet",
            headers=auth_headers,
        )
        assert resp.status_code == 404

    def test_returns_profile(self, client, auth_headers):
        mol_id = _create_molecule(client, auth_headers)
        resp = client.get(
            f"/api/v1/molecule/{mol_id}/admet",
            headers=auth_headers,
        )
        assert resp.status_code == 200
        data = resp.json()
        assert data["id"] == mol_id
        assert data["oral_bioavailability"] in ("high", "moderate", "low")
        assert isinstance(data["cyp_inhibition"], dict)
        assert len(data["cyp_inhibition"]) == 5
        assert isinstance(data["overall_score"], float)
        assert 0.0 <= data["overall_score"] <= 10.0
        assert isinstance(data["warnings"], list)
        assert isinstance(data["flags"], list)
        assert isinstance(data["structural_alerts"], list)
