"""Tests for molecule endpoints."""


def test_parse_ethanol(client, auth_headers):
    resp = client.post(
        "/api/v1/molecule/from-smiles",
        json={"smiles": "CCO"},
        headers=auth_headers,
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["id"]
    assert data["num_atoms"] > 0
    assert data["num_bonds"] > 0


def test_parse_with_name(client, auth_headers):
    resp = client.post(
        "/api/v1/molecule/from-smiles",
        json={"smiles": "CCO", "name": "ethanol"},
        headers=auth_headers,
    )
    assert resp.status_code == 200
    assert resp.json()["name"] == "ethanol"


def test_parse_invalid_smiles(client, auth_headers):
    resp = client.post(
        "/api/v1/molecule/from-smiles",
        json={"smiles": "INVALID!!!"},
        headers=auth_headers,
    )
    assert resp.status_code == 422


def test_get_properties(client, auth_headers):
    # Create molecule first
    create = client.post(
        "/api/v1/molecule/from-smiles",
        json={"smiles": "CCO"},
        headers=auth_headers,
    )
    mol_id = create.json()["id"]
    # Get properties
    resp = client.get(f"/api/v1/molecule/{mol_id}/properties", headers=auth_headers)
    assert resp.status_code == 200
    data = resp.json()
    assert data["formula"]
    assert data["molecular_weight"] > 0
    assert isinstance(data["functional_groups"], list)


def test_get_3d(client, auth_headers):
    create = client.post(
        "/api/v1/molecule/from-smiles",
        json={"smiles": "CCO"},
        headers=auth_headers,
    )
    mol_id = create.json()["id"]
    resp = client.get(f"/api/v1/molecule/{mol_id}/3d", headers=auth_headers)
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["atoms"]) > 0
    assert len(data["bonds"]) > 0
    # Each atom has 3D position
    for atom in data["atoms"]:
        assert len(atom["position"]) == 3


def test_molecule_not_found(client, auth_headers):
    resp = client.get("/api/v1/molecule/nonexistent/properties", headers=auth_headers)
    assert resp.status_code == 404


def test_molecule_requires_auth(client):
    resp = client.post(
        "/api/v1/molecule/from-smiles",
        json={"smiles": "CCO"},
    )
    assert resp.status_code == 401


def test_properties_lipinski_fields(client, auth_headers):
    """Lipinski fields appear in properties response for ethanol (CCO)."""
    create = client.post(
        "/api/v1/molecule/from-smiles",
        json={"smiles": "CCO"},
        headers=auth_headers,
    )
    mol_id = create.json()["id"]
    resp = client.get(f"/api/v1/molecule/{mol_id}/properties", headers=auth_headers)
    assert resp.status_code == 200
    data = resp.json()
    # Check types
    assert isinstance(data["logp"], float)
    assert isinstance(data["hbd"], int)
    assert isinstance(data["hba"], int)
    assert isinstance(data["rotatable_bonds"], int)
    assert isinstance(data["tpsa"], float)
    assert isinstance(data["heavy_atom_count"], int)
    assert isinstance(data["lipinski_violations"], int)
    assert isinstance(data["lipinski_pass"], bool)
    # Check expected values for ethanol
    assert data["hba"] == 1
    assert data["heavy_atom_count"] == 3
    assert data["lipinski_pass"] is True
