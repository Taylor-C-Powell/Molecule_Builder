"""Tests for the /molecule/{id}/solubility endpoint."""

import pytest


def _create_molecule(client, auth_headers, smiles="CCO"):
    resp = client.post(
        "/api/v1/molecule/from-smiles",
        json={"smiles": smiles},
        headers=auth_headers,
    )
    assert resp.status_code == 200
    return resp.json()["id"]


class TestSolubilityEndpoint:
    def test_requires_auth(self, client):
        resp = client.get("/api/v1/molecule/fake-id/solubility")
        assert resp.status_code == 401

    def test_not_found(self, client, auth_headers):
        resp = client.get(
            "/api/v1/molecule/nonexistent/solubility",
            headers=auth_headers,
        )
        assert resp.status_code == 404

    def test_returns_result(self, client, auth_headers):
        mol_id = _create_molecule(client, auth_headers)
        resp = client.get(
            f"/api/v1/molecule/{mol_id}/solubility",
            headers=auth_headers,
        )
        assert resp.status_code == 200
        data = resp.json()
        assert data["id"] == mol_id
        assert "log_s_esol" in data
        assert "log_s_gse" in data
        assert "solubility_mg_ml" in data
        assert data["solubility_class"] in {
            "highly soluble", "soluble", "moderately soluble",
            "poorly soluble", "insoluble",
        }
        assert data["crystallization_risk"] in {"low", "moderate", "high"}
        assert data["polymorph_risk"] in {"low", "moderate", "high"}
