"""Tests for file I/O endpoints (import/export)."""

import io

# -- Sample file content for tests --

SAMPLE_XYZ = """\
3
water
O  0.000  0.000  0.000
H  0.757  0.586  0.000
H -0.757  0.586  0.000
"""

SAMPLE_MOL = """\
water
     MolBuilder

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7570    0.5860    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7570    0.5860    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END
"""

SAMPLE_SDF_MULTI = SAMPLE_MOL + "\n$$$$\n" + SAMPLE_MOL + "\n$$$$\n"

SAMPLE_PDB = """\
HETATM    1  O1  UNK A   1       0.000   0.000   0.000  1.00  0.00           O
HETATM    2  H1  UNK A   1       0.757   0.586   0.000  1.00  0.00           H
HETATM    3  H2  UNK A   1      -0.757   0.586   0.000  1.00  0.00           H
CONECT    1    2    3
CONECT    2    1
CONECT    3    1
END
"""


def _upload(client, headers, content, filename):
    """Helper to POST a file upload."""
    return client.post(
        "/api/v1/molecule/import-file",
        files={"file": (filename, io.BytesIO(content.encode()), "application/octet-stream")},
        headers=headers,
    )


# ---- Import Tests ----


def test_import_xyz(client, auth_headers):
    resp = _upload(client, auth_headers, SAMPLE_XYZ, "water.xyz")
    assert resp.status_code == 200
    data = resp.json()
    assert data["format"] == "xyz"
    assert data["count"] == 1
    assert len(data["molecules"]) == 1
    assert data["molecules"][0]["num_atoms"] == 3


def test_import_mol(client, auth_headers):
    resp = _upload(client, auth_headers, SAMPLE_MOL, "water.mol")
    assert resp.status_code == 200
    data = resp.json()
    assert data["format"] == "mol"
    assert data["count"] == 1


def test_import_sdf_multi(client, auth_headers):
    resp = _upload(client, auth_headers, SAMPLE_SDF_MULTI, "multi.sdf")
    assert resp.status_code == 200
    data = resp.json()
    assert data["format"] == "sdf"
    assert data["count"] == 2


def test_import_pdb(client, auth_headers):
    resp = _upload(client, auth_headers, SAMPLE_PDB, "water.pdb")
    assert resp.status_code == 200
    data = resp.json()
    assert data["format"] == "pdb"
    assert data["count"] == 1


def test_import_unsupported_format(client, auth_headers):
    resp = _upload(client, auth_headers, "random content", "data.csv")
    assert resp.status_code == 422


def test_import_requires_auth(client):
    resp = client.post(
        "/api/v1/molecule/import-file",
        files={"file": ("test.xyz", io.BytesIO(SAMPLE_XYZ.encode()), "application/octet-stream")},
    )
    assert resp.status_code == 401


# ---- Export Tests ----


def _create_mol(client, headers):
    """Helper: parse a SMILES and return mol_id."""
    resp = client.post(
        "/api/v1/molecule/from-smiles",
        json={"smiles": "CCO"},
        headers=headers,
    )
    return resp.json()["id"]


def test_export_xyz(client, auth_headers):
    mol_id = _create_mol(client, auth_headers)
    resp = client.get(f"/api/v1/molecule/{mol_id}/export/xyz", headers=auth_headers)
    assert resp.status_code == 200
    assert "chemical/x-xyz" in resp.headers["content-type"]
    # XYZ first line is atom count
    lines = resp.text.strip().splitlines()
    assert lines[0].strip().isdigit()


def test_export_mol(client, auth_headers):
    mol_id = _create_mol(client, auth_headers)
    resp = client.get(f"/api/v1/molecule/{mol_id}/export/mol", headers=auth_headers)
    assert resp.status_code == 200
    assert "V2000" in resp.text


def test_export_pdb(client, auth_headers):
    mol_id = _create_mol(client, auth_headers)
    resp = client.get(f"/api/v1/molecule/{mol_id}/export/pdb", headers=auth_headers)
    assert resp.status_code == 200
    assert "HETATM" in resp.text or "ATOM" in resp.text


def test_export_json(client, auth_headers):
    mol_id = _create_mol(client, auth_headers)
    resp = client.get(f"/api/v1/molecule/{mol_id}/export/json", headers=auth_headers)
    assert resp.status_code == 200
    data = resp.json()
    assert "atoms" in data


def test_export_not_found(client, auth_headers):
    resp = client.get("/api/v1/molecule/nonexistent/export/xyz", headers=auth_headers)
    assert resp.status_code == 404


def test_export_requires_auth(client):
    resp = client.get("/api/v1/molecule/someid/export/xyz")
    assert resp.status_code == 401


def test_export_unsupported_format(client, auth_headers):
    mol_id = _create_mol(client, auth_headers)
    resp = client.get(f"/api/v1/molecule/{mol_id}/export/csv", headers=auth_headers)
    assert resp.status_code == 422
