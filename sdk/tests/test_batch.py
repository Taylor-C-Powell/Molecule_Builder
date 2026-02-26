"""Tests for batch SDK methods."""

from __future__ import annotations

import httpx
import pytest
import respx

from molbuilder_client import MolBuilder, BatchSubmit, BatchStatus, BatchList


_JOB_JSON = {
    "job_id": "job_abc123",
    "status": "pending",
    "created_at": "2026-01-01T00:00:00",
}

_STATUS_JSON = {
    "job_id": "job_abc123",
    "status": "running",
    "job_type": "properties",
    "progress_pct": 50.0,
    "result": None,
    "error": None,
    "created_at": "2026-01-01T00:00:00",
    "updated_at": "2026-01-01T00:00:01",
}


def test_batch_submit(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.post("/batch/submit").mock(
        return_value=httpx.Response(200, json=_JOB_JSON)
    )
    result = client.batch_submit(["CCO", "C"], "properties")
    assert isinstance(result, BatchSubmit)
    assert result.job_id == "job_abc123"


def test_batch_status(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/batch/job_abc123").mock(
        return_value=httpx.Response(200, json=_STATUS_JSON)
    )
    result = client.batch_status("job_abc123")
    assert isinstance(result, BatchStatus)
    assert result.progress_pct == 50.0


def test_batch_list(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/batch/").mock(
        return_value=httpx.Response(200, json={
            "jobs": [{
                "job_id": "job_abc123",
                "status": "completed",
                "job_type": "properties",
                "progress_pct": 100.0,
                "created_at": "2026-01-01T00:00:00",
                "updated_at": "2026-01-01T00:00:05",
            }],
            "total": 1,
            "page": 1,
            "per_page": 20,
        })
    )
    result = client.batch_list()
    assert isinstance(result, BatchList)
    assert result.total == 1


def test_batch_cancel(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.delete("/batch/job_abc123").mock(
        return_value=httpx.Response(200, json={"status": "cancelled", "job_id": "job_abc123"})
    )
    result = client.batch_cancel("job_abc123")
    assert result is True


def test_batch_tier_limit_403(client: MolBuilder, mock_api: respx.Router) -> None:
    from molbuilder_client import ForbiddenError
    mock_api.post("/batch/submit").mock(
        return_value=httpx.Response(403, json={"detail": "Batch size exceeds tier limit"})
    )
    with pytest.raises(ForbiddenError):
        client.batch_submit(["CCO"] * 100, "properties")


def test_batch_submit_with_params(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.post("/batch/submit").mock(
        return_value=httpx.Response(200, json=_JOB_JSON)
    )
    result = client.batch_submit(
        ["CCO"], "retrosynthesis", params={"max_depth": 3}
    )
    assert result.status == "pending"
