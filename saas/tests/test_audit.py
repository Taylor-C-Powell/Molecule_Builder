"""Tests for audit trail (21 CFR Part 11 compliance)."""

import time
from app.services.audit_db import AuditDB


class TestAuditDB:
    def test_record_creates_entry(self, _temp_audit_db):
        db = _temp_audit_db
        rid = db.record(
            user_email="test@example.com",
            user_role="chemist",
            user_tier="free",
            action="GET /api/v1/elements/C",
        )
        assert rid >= 1
        rec = db.get_record(rid)
        assert rec is not None
        assert rec["user_email"] == "test@example.com"
        assert rec["action"] == "GET /api/v1/elements/C"

    def test_signature_verification_valid(self, _temp_audit_db):
        db = _temp_audit_db
        rid = db.record(
            user_email="test@example.com",
            user_role="chemist",
            user_tier="free",
            action="POST /api/v1/molecule/from-smiles",
            input_summary='{"smiles": "CCO"}',
        )
        result = db.verify_integrity(rid)
        assert result["valid"] is True

    def test_signature_verification_tampered(self, _temp_audit_db):
        db = _temp_audit_db
        rid = db.record(
            user_email="test@example.com",
            user_role="chemist",
            user_tier="free",
            action="POST /api/v1/molecule/from-smiles",
        )
        # Tamper with the record via the backend
        db._backend.execute_update(
            "UPDATE audit_log SET user_email = 'hacker@evil.com' WHERE id = ?",
            (rid,),
        )
        result = db.verify_integrity(rid)
        assert result["valid"] is False
        assert "tampered" in result.get("error", "").lower()

    def test_query_by_user(self, _temp_audit_db):
        db = _temp_audit_db
        db.record(user_email="alice@co.com", user_role="chemist", user_tier="pro", action="GET /health")
        db.record(user_email="bob@co.com", user_role="viewer", user_tier="free", action="GET /health")
        db.record(user_email="alice@co.com", user_role="chemist", user_tier="pro", action="POST /api/v1/molecule/from-smiles")

        results = db.query(user_email="alice@co.com")
        assert len(results) == 2
        assert all(r["user_email"] == "alice@co.com" for r in results)

    def test_query_by_action(self, _temp_audit_db):
        db = _temp_audit_db
        db.record(user_email="u@co.com", user_role="chemist", user_tier="free", action="GET /health")
        db.record(user_email="u@co.com", user_role="chemist", user_tier="free", action="POST /api/v1/molecule/from-smiles")

        results = db.query(action="GET /health")
        assert len(results) == 1

    def test_query_by_date_range(self, _temp_audit_db):
        db = _temp_audit_db
        before = time.time()
        db.record(user_email="u@co.com", user_role="chemist", user_tier="free", action="GET /test")
        after = time.time()

        results = db.query(start_time=before, end_time=after)
        assert len(results) >= 1

    def test_count(self, _temp_audit_db):
        db = _temp_audit_db
        db.record(user_email="u@co.com", user_role="chemist", user_tier="free", action="A")
        db.record(user_email="u@co.com", user_role="chemist", user_tier="free", action="B")
        db.record(user_email="u@co.com", user_role="chemist", user_tier="free", action="C")
        assert db.count() == 3

    def test_export(self, _temp_audit_db):
        db = _temp_audit_db
        db.record(user_email="u@co.com", user_role="chemist", user_tier="free", action="A")
        db.record(user_email="u@co.com", user_role="chemist", user_tier="free", action="B")
        exported = db.export_all()
        assert len(exported) == 2
        assert exported[0]["id"] < exported[1]["id"]  # ordered by id ASC


class TestAuditEndpoints:
    def test_audit_log_admin_access(self, client, admin_headers):
        # Make a request to generate audit data
        client.get("/api/v1/elements/C", headers=admin_headers)

        resp = client.get("/api/v1/audit/log", headers=admin_headers)
        assert resp.status_code == 200
        data = resp.json()
        assert "records" in data
        assert "total" in data
        assert data["total"] >= 1

    def test_audit_log_non_admin_denied(self, client, auth_headers):
        resp = client.get("/api/v1/audit/log", headers=auth_headers)
        assert resp.status_code == 403

    def test_audit_single_record(self, client, admin_headers):
        # Generate some requests
        client.get("/api/v1/elements/C", headers=admin_headers)

        log = client.get("/api/v1/audit/log", headers=admin_headers).json()
        if log["records"]:
            rid = log["records"][0]["id"]
            resp = client.get(f"/api/v1/audit/log/{rid}", headers=admin_headers)
            assert resp.status_code == 200
            assert resp.json()["id"] == rid

    def test_audit_verify_integrity(self, client, admin_headers):
        client.get("/api/v1/elements/C", headers=admin_headers)

        log = client.get("/api/v1/audit/log", headers=admin_headers).json()
        if log["records"]:
            rid = log["records"][0]["id"]
            resp = client.get(f"/api/v1/audit/log/{rid}/verify", headers=admin_headers)
            assert resp.status_code == 200
            assert resp.json()["valid"] is True

    def test_audit_export(self, client, admin_headers):
        client.get("/api/v1/elements/C", headers=admin_headers)

        resp = client.get("/api/v1/audit/export", headers=admin_headers)
        assert resp.status_code == 200
        data = resp.json()
        assert isinstance(data, list)
        assert len(data) >= 1

    def test_audit_records_created_per_request(self, client, admin_headers, _temp_audit_db):
        db = _temp_audit_db
        initial = db.count()
        client.get("/api/v1/elements/C", headers=admin_headers)
        assert db.count() > initial
