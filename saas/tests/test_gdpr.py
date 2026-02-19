"""Tests for GDPR self-service endpoints (data export and account deletion)."""


class TestDataExport:
    def test_export_requires_auth(self, client):
        resp = client.get("/api/v1/auth/me/export")
        assert resp.status_code == 401

    def test_export_returns_user_data(self, client, auth_headers):
        resp = client.get("/api/v1/auth/me/export", headers=auth_headers)
        assert resp.status_code == 200
        data = resp.json()
        assert "email" in data
        assert "tier" in data
        assert "role" in data
        assert "billing" in data

    def test_export_contains_billing_info(self, client, auth_headers):
        resp = client.get("/api/v1/auth/me/export", headers=auth_headers)
        data = resp.json()
        assert "subscription_status" in data["billing"]
        assert "has_stripe_account" in data["billing"]


class TestAccountDeletion:
    def test_delete_requires_auth(self, client):
        resp = client.delete("/api/v1/auth/me")
        assert resp.status_code == 401

    def test_delete_revokes_keys(self, client):
        # Register a user
        reg = client.post("/api/v1/auth/register", json={"email": "gdpr@test.com"})
        key = reg.json()["api_key"]

        # Get a token
        tok = client.post("/api/v1/auth/token", json={"api_key": key})
        token = tok.json()["access_token"]
        headers = {"Authorization": f"Bearer {token}"}

        # Delete account
        resp = client.delete("/api/v1/auth/me", headers=headers)
        assert resp.status_code == 200
        data = resp.json()
        assert data["deleted"] is True
        assert data["keys_revoked"] >= 1

        # Verify the old key no longer works
        verify = client.post("/api/v1/auth/token", json={"api_key": key})
        assert verify.status_code == 401
