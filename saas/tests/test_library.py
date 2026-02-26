"""Tests for molecule library endpoints."""

import uuid
import pytest
from app.auth.api_keys import api_key_store
from app.config import Tier


@pytest.fixture()
def free_headers():
    email = f"lib-free-{uuid.uuid4().hex[:8]}@example.com"
    key = api_key_store.create(email=email, tier=Tier.FREE)
    return {"X-API-Key": key}


@pytest.fixture()
def pro_headers():
    email = f"lib-pro-{uuid.uuid4().hex[:8]}@example.com"
    key = api_key_store.create(email=email, tier=Tier.PRO)
    return {"X-API-Key": key}


class TestLibrarySave:
    def test_save_molecule(self, client, free_headers):
        resp = client.post("/api/v1/library/", json={
            "smiles": "CCO",
            "name": "ethanol",
            "tags": ["alcohol", "solvent"],
            "notes": "Common solvent",
        }, headers=free_headers)
        assert resp.status_code == 200
        data = resp.json()
        assert data["name"] == "ethanol"
        assert "alcohol" in data["tags"]
        assert data["properties"]["molecular_weight"] > 0

    def test_duplicate_detection(self, client, free_headers):
        client.post("/api/v1/library/", json={
            "smiles": "CCO",
        }, headers=free_headers)
        resp = client.post("/api/v1/library/", json={
            "smiles": "CCO",
        }, headers=free_headers)
        assert resp.status_code == 409

    def test_invalid_smiles(self, client, free_headers):
        resp = client.post("/api/v1/library/", json={
            "smiles": "INVALID_SMILES_XYZ",
        }, headers=free_headers)
        assert resp.status_code == 400


class TestLibraryRetrieve:
    def test_get_molecule(self, client, free_headers):
        resp = client.post("/api/v1/library/", json={
            "smiles": "C",
            "name": "methane",
        }, headers=free_headers)
        assert resp.status_code == 200, resp.text
        mol_id = resp.json()["id"]

        resp = client.get(f"/api/v1/library/{mol_id}", headers=free_headers)
        assert resp.status_code == 200
        assert resp.json()["name"] == "methane"

    def test_get_nonexistent(self, client, free_headers):
        resp = client.get("/api/v1/library/99999", headers=free_headers)
        assert resp.status_code == 404


class TestLibraryList:
    def test_list_molecules(self, client, free_headers):
        client.post("/api/v1/library/", json={"smiles": "CCO", "name": "ethanol"}, headers=free_headers)
        client.post("/api/v1/library/", json={"smiles": "CC", "name": "ethane"}, headers=free_headers)

        resp = client.get("/api/v1/library/", headers=free_headers)
        assert resp.status_code == 200
        data = resp.json()
        assert data["total"] == 2

    def test_tag_filter(self, client, free_headers):
        client.post("/api/v1/library/", json={"smiles": "CCO", "tags": ["solvent"]}, headers=free_headers)
        client.post("/api/v1/library/", json={"smiles": "CC", "tags": ["hydrocarbon"]}, headers=free_headers)

        resp = client.get("/api/v1/library/?tag=solvent", headers=free_headers)
        assert resp.status_code == 200
        data = resp.json()
        assert data["total"] == 1

    def test_search_by_name(self, client, free_headers):
        client.post("/api/v1/library/", json={"smiles": "CCO", "name": "ethanol"}, headers=free_headers)
        client.post("/api/v1/library/", json={"smiles": "CC", "name": "ethane"}, headers=free_headers)

        resp = client.get("/api/v1/library/?search=ethanol", headers=free_headers)
        assert resp.status_code == 200
        data = resp.json()
        assert data["total"] == 1


class TestLibraryUpdate:
    def test_update_name_and_tags(self, client, free_headers):
        resp = client.post("/api/v1/library/", json={"smiles": "CCO"}, headers=free_headers)
        mol_id = resp.json()["id"]

        resp = client.put(f"/api/v1/library/{mol_id}", json={
            "name": "ethyl alcohol",
            "tags": ["updated"],
        }, headers=free_headers)
        assert resp.status_code == 200
        assert resp.json()["name"] == "ethyl alcohol"
        assert "updated" in resp.json()["tags"]


class TestLibraryDelete:
    def test_delete_molecule(self, client, free_headers):
        resp = client.post("/api/v1/library/", json={"smiles": "CCO"}, headers=free_headers)
        mol_id = resp.json()["id"]

        resp = client.delete(f"/api/v1/library/{mol_id}", headers=free_headers)
        assert resp.status_code == 200

        resp = client.get(f"/api/v1/library/{mol_id}", headers=free_headers)
        assert resp.status_code == 404


class TestLibraryImport:
    def test_bulk_import(self, client, free_headers):
        resp = client.post("/api/v1/library/import", json={
            "smiles_list": ["CCO", "CC", "C"],
            "tag": "imported",
        }, headers=free_headers)
        assert resp.status_code == 200
        data = resp.json()
        assert data["saved"] == 3
        assert data["duplicates"] == 0

    def test_import_with_invalid(self, client, free_headers):
        resp = client.post("/api/v1/library/import", json={
            "smiles_list": ["CCO", "INVALID_XYZ", "CC"],
        }, headers=free_headers)
        assert resp.status_code == 200
        data = resp.json()
        assert data["saved"] == 2
        assert len(data["errors"]) == 1

    def test_import_with_duplicates(self, client, free_headers):
        client.post("/api/v1/library/", json={"smiles": "CCO"}, headers=free_headers)
        resp = client.post("/api/v1/library/import", json={
            "smiles_list": ["CCO", "CC"],
        }, headers=free_headers)
        assert resp.status_code == 200
        data = resp.json()
        assert data["saved"] == 1
        assert data["duplicates"] == 1


class TestLibraryTierLimit:
    def test_free_tier_save_succeeds(self, client, free_headers):
        """Free tier can save molecules (up to 50 limit)."""
        resp = client.post("/api/v1/library/", json={"smiles": "C"}, headers=free_headers)
        assert resp.status_code == 200
