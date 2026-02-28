"""Tests for middleware (Request ID, Security Headers)."""

import uuid


class TestRequestIDMiddleware:
    def test_response_has_request_id_header(self, client):
        resp = client.get("/health")
        assert "X-Request-ID" in resp.headers

    def test_request_id_is_valid_uuid(self, client):
        resp = client.get("/health")
        request_id = resp.headers["X-Request-ID"]
        # Should parse as a valid UUID
        parsed = uuid.UUID(request_id)
        assert str(parsed) == request_id

    def test_request_id_unique_per_request(self, client):
        r1 = client.get("/health")
        r2 = client.get("/health")
        assert r1.headers["X-Request-ID"] != r2.headers["X-Request-ID"]


class TestSecurityHeaders:
    def test_nosniff_header(self, client):
        resp = client.get("/health")
        assert resp.headers.get("X-Content-Type-Options") == "nosniff"

    def test_frame_deny_header(self, client):
        resp = client.get("/health")
        assert resp.headers.get("X-Frame-Options") == "DENY"

    def test_xss_protection_header(self, client):
        resp = client.get("/health")
        assert resp.headers.get("X-XSS-Protection") == "1; mode=block"

    def test_referrer_policy_header(self, client):
        resp = client.get("/health")
        assert resp.headers.get("Referrer-Policy") == "strict-origin-when-cross-origin"
