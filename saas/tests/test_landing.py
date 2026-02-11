"""Tests for landing page."""


def test_landing_returns_200(client):
    resp = client.get("/")
    assert resp.status_code == 200


def test_landing_content_type_is_html(client):
    resp = client.get("/")
    assert "text/html" in resp.headers["content-type"]


def test_landing_contains_product_name(client):
    resp = client.get("/")
    assert "MolBuilder" in resp.text


def test_landing_contains_api_endpoints(client):
    resp = client.get("/")
    assert "/api/v1/molecule/from-smiles" in resp.text
    assert "/api/v1/retrosynthesis/plan" in resp.text


def test_landing_contains_pricing(client):
    resp = client.get("/")
    assert "$0" in resp.text
    assert "$49" in resp.text


def test_landing_not_in_openapi_schema(client):
    resp = client.get("/openapi.json")
    paths = resp.json()["paths"]
    assert "/" not in paths
