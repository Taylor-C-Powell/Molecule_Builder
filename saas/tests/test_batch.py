"""Tests for batch processing endpoints."""

import time
import uuid
import pytest
from app.auth.api_keys import api_key_store
from app.config import Tier


@pytest.fixture()
def free_headers():
    email = f"batch-free-{uuid.uuid4().hex[:8]}@example.com"
    key = api_key_store.create(email=email, tier=Tier.FREE)
    return {"X-API-Key": key}


@pytest.fixture()
def pro_headers():
    email = f"batch-pro-{uuid.uuid4().hex[:8]}@example.com"
    key = api_key_store.create(email=email, tier=Tier.PRO)
    return {"X-API-Key": key}


class TestBatchSubmit:
    def test_submit_and_poll(self, client, free_headers):
        resp = client.post("/api/v1/batch/submit", json={
            "smiles_list": ["CCO", "CC", "C"],
            "job_type": "properties",
        }, headers=free_headers)
        assert resp.status_code == 200
        data = resp.json()
        assert "job_id" in data
        assert data["status"] in ("pending", "running")

        # Poll until complete (max 10s)
        job_id = data["job_id"]
        for _ in range(20):
            status = client.get(f"/api/v1/batch/{job_id}", headers=free_headers)
            assert status.status_code == 200
            sdata = status.json()
            if sdata["status"] in ("completed", "failed"):
                break
            time.sleep(0.5)

        assert sdata["status"] == "completed"
        assert sdata["result"] is not None
        assert sdata["result"]["total"] == 3
        assert sdata["result"]["succeeded"] == 3

    def test_submit_with_invalid_smiles(self, client, free_headers):
        resp = client.post("/api/v1/batch/submit", json={
            "smiles_list": ["CCO", "INVALID_SMILES", "CC"],
            "job_type": "evaluate",
        }, headers=free_headers)
        assert resp.status_code == 200
        job_id = resp.json()["job_id"]

        for _ in range(20):
            status = client.get(f"/api/v1/batch/{job_id}", headers=free_headers)
            sdata = status.json()
            if sdata["status"] in ("completed", "failed"):
                break
            time.sleep(0.5)

        assert sdata["status"] == "completed"
        result = sdata["result"]
        assert result["failed"] >= 1

    def test_tier_limit_enforcement(self, client, free_headers):
        """Free tier limited to 10 molecules per batch."""
        resp = client.post("/api/v1/batch/submit", json={
            "smiles_list": ["C"] * 11,
            "job_type": "evaluate",
        }, headers=free_headers)
        assert resp.status_code == 403

    def test_pro_tier_higher_limit(self, client, pro_headers):
        """Pro tier allows up to 100."""
        resp = client.post("/api/v1/batch/submit", json={
            "smiles_list": ["C"] * 20,
            "job_type": "evaluate",
        }, headers=pro_headers)
        assert resp.status_code == 200


class TestBatchList:
    def test_list_jobs(self, client, free_headers):
        # Submit a job first
        client.post("/api/v1/batch/submit", json={
            "smiles_list": ["C"],
            "job_type": "evaluate",
        }, headers=free_headers)

        resp = client.get("/api/v1/batch/", headers=free_headers)
        assert resp.status_code == 200
        data = resp.json()
        assert data["total"] >= 1
        assert len(data["jobs"]) >= 1


class TestBatchCancel:
    def test_cancel_job(self, client, free_headers):
        resp = client.post("/api/v1/batch/submit", json={
            "smiles_list": ["C"] * 5,
            "job_type": "properties",
        }, headers=free_headers)
        job_id = resp.json()["job_id"]

        # Try to cancel (may or may not succeed depending on timing/rate limits)
        cancel = client.delete(f"/api/v1/batch/{job_id}", headers=free_headers)
        # 200 cancelled, 404 already done, 429 rate limited
        assert cancel.status_code in (200, 404, 429)
