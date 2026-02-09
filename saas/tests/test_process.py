"""Tests for process engineering endpoint."""


def test_process_evaluate(client, auth_headers):
    resp = client.post(
        "/api/v1/process/evaluate",
        json={"smiles": "CC(=O)OCC", "scale_kg": 1.0},
        headers=auth_headers,
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["smiles"] == "CC(=O)OCC"
    assert data["scale_kg"] == 1.0
    assert isinstance(data["route_found"], bool)


def test_process_with_route(client, auth_headers):
    resp = client.post(
        "/api/v1/process/evaluate",
        json={"smiles": "CC(=O)OCC", "scale_kg": 1.0},
        headers=auth_headers,
    )
    data = resp.json()
    if data["route_found"]:
        assert data["total_steps"] > 0
        assert len(data["step_details"]) > 0
        assert data["cost"] is not None
        assert data["cost"]["total_usd"] > 0
        assert data["scale_up"] is not None
        assert len(data["safety"]) > 0
        # Check step detail structure
        step = data["step_details"][0]
        assert "reactor" in step
        assert "conditions" in step
        assert "purification" in step


def test_process_simple_molecule(client, auth_headers):
    resp = client.post(
        "/api/v1/process/evaluate",
        json={"smiles": "CCO", "scale_kg": 0.5},
        headers=auth_headers,
    )
    assert resp.status_code == 200


def test_process_invalid_smiles(client, auth_headers):
    resp = client.post(
        "/api/v1/process/evaluate",
        json={"smiles": "NOPE!!!", "scale_kg": 1.0},
        headers=auth_headers,
    )
    assert resp.status_code == 422


def test_process_requires_auth(client):
    resp = client.post(
        "/api/v1/process/evaluate",
        json={"smiles": "CCO", "scale_kg": 1.0},
    )
    assert resp.status_code == 401
