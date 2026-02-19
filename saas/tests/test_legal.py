"""Tests for legal endpoints (Terms of Service and Privacy Policy)."""


class TestTermsOfService:
    def test_terms_returns_200(self, client):
        resp = client.get("/api/v1/legal/terms")
        assert resp.status_code == 200

    def test_terms_has_title(self, client):
        data = client.get("/api/v1/legal/terms").json()
        assert "title" in data
        assert "Terms of Service" in data["title"]

    def test_terms_has_sections(self, client):
        data = client.get("/api/v1/legal/terms").json()
        assert "sections" in data
        assert len(data["sections"]) >= 10

    def test_terms_no_auth_required(self, client):
        """Legal endpoints must be publicly accessible."""
        resp = client.get("/api/v1/legal/terms")
        assert resp.status_code == 200


class TestPrivacyPolicy:
    def test_privacy_returns_200(self, client):
        resp = client.get("/api/v1/legal/privacy")
        assert resp.status_code == 200

    def test_privacy_has_title(self, client):
        data = client.get("/api/v1/legal/privacy").json()
        assert "Privacy Policy" in data["title"]

    def test_privacy_has_gdpr_section(self, client):
        data = client.get("/api/v1/legal/privacy").json()
        headings = [s["heading"] for s in data["sections"]]
        assert any("GDPR" in h for h in headings)

    def test_privacy_has_molecular_data_section(self, client):
        data = client.get("/api/v1/legal/privacy").json()
        headings = [s["heading"] for s in data["sections"]]
        assert any("Molecular" in h for h in headings)
