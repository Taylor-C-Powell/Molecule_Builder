"""Tests for API versioning middleware and endpoint."""

from fastapi.testclient import TestClient

from app.main import app
from app.middleware_versioning import API_VERSION, MIN_SUPPORTED_VERSION


client = TestClient(app)


class TestVersionHeaders:
    """Verify version headers are present on all responses."""

    def test_health_has_version_headers(self):
        resp = client.get("/health")
        assert resp.headers["X-API-Version"] == API_VERSION
        assert resp.headers["X-Min-Supported-Version"] == MIN_SUPPORTED_VERSION

    def test_version_endpoint(self):
        resp = client.get("/api/v1/version")
        assert resp.status_code == 200
        data = resp.json()
        assert data["api_version"] == API_VERSION
        assert data["min_supported_version"] == MIN_SUPPORTED_VERSION
        assert "12 months" in data["deprecation_policy"]

    def test_deprecation_warning_for_old_version(self):
        resp = client.get(
            "/health",
            headers={"X-API-Version": "0.0.1"},
        )
        assert "X-Deprecation-Warning" in resp.headers

    def test_no_deprecation_warning_for_current_version(self):
        resp = client.get(
            "/health",
            headers={"X-API-Version": API_VERSION},
        )
        assert "X-Deprecation-Warning" not in resp.headers
