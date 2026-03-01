"""Tests for auth endpoints and middleware."""

import pytest


def test_register_returns_key(client):
    resp = client.post("/api/v1/auth/register", json={"email": "u@test.com"})
    assert resp.status_code == 200
    data = resp.json()
    assert data["api_key"].startswith("mb_free_")
    assert data["email"] == "u@test.com"
    assert data["tier"] == "free"


def test_register_ignores_tier(client):
    """Self-registration always produces free-tier keys regardless of input."""
    resp = client.post(
        "/api/v1/auth/register", json={"email": "p@test.com", "tier": "pro"}
    )
    assert resp.status_code == 200
    assert resp.json()["api_key"].startswith("mb_free_")
    assert resp.json()["tier"] == "free"


def test_provision_requires_admin(client):
    """The /provision endpoint rejects non-admin users."""
    # Register a free user first
    reg = client.post("/api/v1/auth/register", json={"email": "user@test.com"})
    key = reg.json()["api_key"]
    resp = client.post(
        "/api/v1/auth/provision",
        json={"email": "cust@test.com", "tier": "pro"},
        headers={"X-API-Key": key},
    )
    assert resp.status_code == 403


def test_provision_by_admin(client, admin_headers):
    """Admins can provision keys with custom tiers and roles."""
    resp = client.post(
        "/api/v1/auth/provision",
        json={"email": "cust@test.com", "tier": "pro", "role": "chemist"},
        headers=admin_headers,
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["api_key"].startswith("mb_pro_")
    assert data["tier"] == "pro"
    assert data["role"] == "chemist"


def test_token_exchange(client):
    # Register first
    reg = client.post("/api/v1/auth/register", json={"email": "tok@test.com"})
    key = reg.json()["api_key"]
    # Exchange for JWT
    resp = client.post("/api/v1/auth/token", json={"api_key": key})
    assert resp.status_code == 200
    data = resp.json()
    assert data["token_type"] == "bearer"
    assert "access_token" in data
    assert data["expires_in"] == 3600
    assert data["email"] == "tok@test.com"
    assert data["tier"] == "free"


def test_token_invalid_key(client):
    resp = client.post("/api/v1/auth/token", json={"api_key": "mb_free_bogus"})
    assert resp.status_code == 401


def test_jwt_auth_works(client):
    # Register + get token
    reg = client.post("/api/v1/auth/register", json={"email": "jwt@test.com"})
    key = reg.json()["api_key"]
    tok = client.post("/api/v1/auth/token", json={"api_key": key})
    token = tok.json()["access_token"]
    # Use JWT to hit authenticated endpoint
    resp = client.get(
        "/api/v1/elements/C",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.status_code == 200


def test_no_auth_returns_401(client):
    resp = client.get("/api/v1/elements/C")
    assert resp.status_code == 401


def test_bad_api_key_returns_401(client):
    resp = client.get("/api/v1/elements/C", headers={"X-API-Key": "invalid"})
    assert resp.status_code == 401


def test_bad_jwt_returns_401(client):
    resp = client.get(
        "/api/v1/elements/C",
        headers={"Authorization": "Bearer invalid.jwt.token"},
    )
    assert resp.status_code == 401


def test_register_academic_edu(client):
    """Academic .edu emails get academic tier automatically."""
    resp = client.post("/api/v1/auth/register", json={"email": "prof@mit.edu"})
    assert resp.status_code == 200
    data = resp.json()
    assert data["tier"] == "academic"
    assert data["api_key"].startswith("mb_acad")


def test_register_academic_ac_uk(client):
    """Academic .ac.uk emails get academic tier."""
    resp = client.post("/api/v1/auth/register", json={"email": "researcher@ox.ac.uk"})
    assert resp.status_code == 200
    assert resp.json()["tier"] == "academic"


def test_register_non_academic(client):
    """Non-academic emails remain free tier."""
    resp = client.post("/api/v1/auth/register", json={"email": "user@gmail.com"})
    assert resp.status_code == 200
    assert resp.json()["tier"] == "free"
    assert resp.json()["api_key"].startswith("mb_free_")


def test_register_academic_edu_au(client):
    """Australian .edu.au emails get academic tier."""
    resp = client.post("/api/v1/auth/register", json={"email": "student@unsw.edu.au"})
    assert resp.status_code == 200
    assert resp.json()["tier"] == "academic"


def test_academic_detection_case_insensitive(client):
    """Academic detection is case-insensitive."""
    resp = client.post("/api/v1/auth/register", json={"email": "Prof@MIT.EDU"})
    assert resp.status_code == 200
    assert resp.json()["tier"] == "academic"


def test_rate_limit_returns_429(client, auth_headers):
    # Free tier: 10 rpm. Exhaust the limit then check 429.
    for _ in range(10):
        client.get("/api/v1/elements/C", headers=auth_headers)
    resp = client.get("/api/v1/elements/C", headers=auth_headers)
    assert resp.status_code == 429
    assert "Retry-After" in resp.headers
