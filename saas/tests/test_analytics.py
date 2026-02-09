"""Tests for usage tracking and analytics endpoints."""

from app.services.usage_tracker import usage_tracker


def test_analytics_summary_after_requests(client, admin_headers):
    # Make some tracked requests
    client.get("/api/v1/elements/C", headers=admin_headers)
    client.post(
        "/api/v1/molecule/from-smiles",
        json={"smiles": "CCO"},
        headers=admin_headers,
    )
    client.post(
        "/api/v1/molecule/from-smiles",
        json={"smiles": "CCO"},
        headers=admin_headers,
    )

    resp = client.get("/api/v1/analytics/summary", headers=admin_headers)
    assert resp.status_code == 200
    data = resp.json()
    assert data["total_requests"] >= 3
    assert "by_endpoint" in data
    assert "by_tier" in data
    assert "top_smiles" in data
    # CCO should appear in top SMILES
    assert "CCO" in data["top_smiles"]


def test_analytics_recent(client, admin_headers):
    client.get("/api/v1/elements/O", headers=admin_headers)

    resp = client.get("/api/v1/analytics/recent?limit=5", headers=admin_headers)
    assert resp.status_code == 200
    records = resp.json()
    assert isinstance(records, list)
    assert len(records) >= 1
    rec = records[0]
    assert "method" in rec
    assert "path" in rec
    assert "status_code" in rec
    assert "latency_ms" in rec


def test_analytics_user_usage(client, admin_headers):
    # Make a request so there's data
    client.get("/api/v1/elements/H", headers=admin_headers)

    # Get the email from the summary
    summary = client.get("/api/v1/analytics/summary", headers=admin_headers).json()
    if summary["by_user"]:
        email = list(summary["by_user"].keys())[0]
        resp = client.get(f"/api/v1/analytics/user/{email}", headers=admin_headers)
        assert resp.status_code == 200
        data = resp.json()
        assert data["email"] == email
        assert data["total_requests"] >= 1
        assert "endpoints" in data


def test_analytics_tracks_smiles(client, admin_headers):
    client.post(
        "/api/v1/retrosynthesis/plan",
        json={"smiles": "c1ccccc1"},
        headers=admin_headers,
    )

    summary = client.get("/api/v1/analytics/summary", headers=admin_headers).json()
    assert "c1ccccc1" in summary["top_smiles"]


def test_analytics_tracks_errors(client, admin_headers):
    # Trigger a 404
    client.get("/api/v1/molecule/nonexistent/properties", headers=admin_headers)

    summary = client.get("/api/v1/analytics/summary", headers=admin_headers).json()
    assert summary["total_errors"] >= 1


def test_analytics_requires_auth(client):
    resp = client.get("/api/v1/analytics/summary")
    assert resp.status_code == 401


def test_analytics_forbidden_for_non_admin(client, auth_headers):
    resp = client.get("/api/v1/analytics/summary", headers=auth_headers)
    assert resp.status_code == 403


def test_analytics_latency_tracked(client, admin_headers):
    client.get("/api/v1/elements/N", headers=admin_headers)

    summary = client.get("/api/v1/analytics/summary", headers=admin_headers).json()
    # At least one endpoint should have latency data
    for endpoint_data in summary["by_endpoint"].values():
        if endpoint_data["calls"] > 0:
            assert endpoint_data["avg_latency_ms"] >= 0
            assert endpoint_data["max_latency_ms"] >= 0
            break
