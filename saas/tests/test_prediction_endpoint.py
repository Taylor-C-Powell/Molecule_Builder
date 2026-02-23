"""Tests for POST /api/v1/process/predict-conditions endpoint."""


def test_predict_conditions_basic(client, auth_headers):
    resp = client.post(
        "/api/v1/process/predict-conditions",
        json={"smiles": "CCO"},
        headers=auth_headers,
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["smiles"] == "CCO"
    assert "substrate_analysis" in data
    assert "candidates" in data
    assert data["overall_confidence"] in ("high", "medium", "low")


def test_predict_conditions_with_hint(client, auth_headers):
    resp = client.post(
        "/api/v1/process/predict-conditions",
        json={"smiles": "CCO", "reaction_name": "oxidation", "scale_kg": 1.0},
        headers=auth_headers,
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["best_match"] is not None
    best = data["best_match"]
    assert "conditions" in best
    assert "match_score" in best
    assert best["match_score"] >= 0


def test_predict_conditions_substrate_analysis(client, auth_headers):
    resp = client.post(
        "/api/v1/process/predict-conditions",
        json={"smiles": "CCO"},
        headers=auth_headers,
    )
    data = resp.json()
    sa = data["substrate_analysis"]
    assert "detected_functional_groups" in sa
    assert sa["heavy_atom_count"] == 3
    assert sa["steric_class"] in ("unhindered", "moderately_hindered", "hindered")
    assert sa["electronic_character"] in ("electron_rich", "neutral", "electron_poor")


def test_predict_conditions_max_candidates(client, auth_headers):
    resp = client.post(
        "/api/v1/process/predict-conditions",
        json={"smiles": "CCO", "max_candidates": 2},
        headers=auth_headers,
    )
    data = resp.json()
    assert len(data["candidates"]) <= 2


def test_predict_conditions_candidate_structure(client, auth_headers):
    resp = client.post(
        "/api/v1/process/predict-conditions",
        json={"smiles": "CCO", "reaction_name": "oxidation"},
        headers=auth_headers,
    )
    data = resp.json()
    if data["candidates"]:
        cand = data["candidates"][0]
        assert "template_name" in cand
        assert "category" in cand
        assert "conditions" in cand
        assert "recommended_solvents" in cand
        assert "adjusted_yield_range" in cand
        assert len(cand["adjusted_yield_range"]) == 2
        assert "warnings" in cand


def test_predict_conditions_invalid_smiles(client, auth_headers):
    resp = client.post(
        "/api/v1/process/predict-conditions",
        json={"smiles": "NOPE!!!"},
        headers=auth_headers,
    )
    assert resp.status_code == 422


def test_predict_conditions_requires_auth(client):
    resp = client.post(
        "/api/v1/process/predict-conditions",
        json={"smiles": "CCO"},
    )
    assert resp.status_code == 401


def test_predict_conditions_scale_param(client, auth_headers):
    resp = client.post(
        "/api/v1/process/predict-conditions",
        json={"smiles": "CCO", "scale_kg": 100.0},
        headers=auth_headers,
    )
    assert resp.status_code == 200
    data = resp.json()
    if data["best_match"]:
        cond = data["best_match"]["conditions"]
        assert cond["temperature_C"] is not None
        assert cond["solvent"] is not None


def test_predict_conditions_alkyl_halide(client, auth_headers):
    resp = client.post(
        "/api/v1/process/predict-conditions",
        json={"smiles": "CCCl", "reaction_name": "substitution"},
        headers=auth_headers,
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["best_match"] is not None
    assert data["best_match"]["category"] == "SUBSTITUTION"
