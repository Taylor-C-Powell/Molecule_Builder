"""Tests for retrosynthesis endpoint."""


def test_retro_plan(client, auth_headers):
    resp = client.post(
        "/api/v1/retrosynthesis/plan",
        json={"smiles": "CC(=O)OCC"},
        headers=auth_headers,
    )
    assert resp.status_code == 200
    data = resp.json()
    assert "tree" in data
    assert "routes_found" in data
    assert data["tree"]["smiles"]


def test_retro_with_params(client, auth_headers):
    resp = client.post(
        "/api/v1/retrosynthesis/plan",
        json={"smiles": "CC(=O)OCC", "max_depth": 3, "beam_width": 3},
        headers=auth_headers,
    )
    assert resp.status_code == 200
    assert resp.json()["max_depth"] == 3
    assert resp.json()["beam_width"] == 3


def test_retro_simple_molecule(client, auth_headers):
    resp = client.post(
        "/api/v1/retrosynthesis/plan",
        json={"smiles": "CCO"},
        headers=auth_headers,
    )
    assert resp.status_code == 200


def test_retro_invalid_smiles(client, auth_headers):
    resp = client.post(
        "/api/v1/retrosynthesis/plan",
        json={"smiles": "NOT_SMILES"},
        headers=auth_headers,
    )
    assert resp.status_code == 422


def test_retro_requires_auth(client):
    resp = client.post(
        "/api/v1/retrosynthesis/plan",
        json={"smiles": "CCO"},
    )
    assert resp.status_code == 401


def test_retro_disconnections_field_present(client, auth_headers):
    """The tree node should include a disconnections list."""
    resp = client.post(
        "/api/v1/retrosynthesis/plan",
        json={"smiles": "CC(=O)OCC", "max_depth": 2, "beam_width": 3},
        headers=auth_headers,
    )
    assert resp.status_code == 200
    tree = resp.json()["tree"]
    assert "disconnections" in tree
    assert isinstance(tree["disconnections"], list)


def test_retro_disconnections_contain_best(client, auth_headers):
    """When best_disconnection exists, it should appear in disconnections too."""
    resp = client.post(
        "/api/v1/retrosynthesis/plan",
        json={"smiles": "CC(=O)OCC", "max_depth": 3, "beam_width": 5},
        headers=auth_headers,
    )
    assert resp.status_code == 200
    tree = resp.json()["tree"]
    if tree["best_disconnection"] is not None:
        best_name = tree["best_disconnection"]["reaction_name"]
        disc_names = [d["reaction_name"] for d in tree["disconnections"]]
        assert best_name in disc_names
        # Each disconnection should have required fields
        for d in tree["disconnections"]:
            assert "reaction_name" in d
            assert "score" in d
            assert "precursors" in d
            assert isinstance(d["score"], (int, float))
