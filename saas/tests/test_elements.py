"""Tests for elements endpoint."""


def test_get_carbon(client, auth_headers):
    resp = client.get("/api/v1/elements/C", headers=auth_headers)
    assert resp.status_code == 200
    data = resp.json()
    assert data["symbol"] == "C"
    assert data["name"] == "Carbon"
    assert data["atomic_number"] == 6
    assert data["atomic_weight"] > 12.0


def test_get_oxygen(client, auth_headers):
    resp = client.get("/api/v1/elements/O", headers=auth_headers)
    assert resp.status_code == 200
    assert resp.json()["name"] == "Oxygen"


def test_element_not_found(client, auth_headers):
    resp = client.get("/api/v1/elements/Xx", headers=auth_headers)
    assert resp.status_code == 404


def test_element_requires_auth(client):
    resp = client.get("/api/v1/elements/C")
    assert resp.status_code == 401
