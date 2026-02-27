"""Tests for SDK file I/O methods: import_file, export_file."""

from __future__ import annotations

import httpx
import respx

from molbuilder_client import FileImportResult, MolBuilder


def test_import_file(client: MolBuilder, mock_api: respx.Router, tmp_path) -> None:
    mock_api.post("/molecule/import-file").mock(
        return_value=httpx.Response(
            200,
            json={
                "molecules": [
                    {
                        "id": "mol_001",
                        "name": "water",
                        "smiles": "O",
                        "num_atoms": 3,
                        "num_bonds": 2,
                    }
                ],
                "format": "xyz",
                "count": 1,
            },
        )
    )

    # Write a dummy XYZ file
    xyz_file = tmp_path / "water.xyz"
    xyz_file.write_text("3\nwater\nO 0 0 0\nH 1 0 0\nH 0 1 0\n")

    result = client.import_file(str(xyz_file))

    assert isinstance(result, FileImportResult)
    assert result.count == 1
    assert result.format == "xyz"
    assert len(result.molecules) == 1
    assert result.molecules[0].id == "mol_001"


def test_export_file(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/molecule/mol_001/export/xyz").mock(
        return_value=httpx.Response(
            200,
            text="3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n",
            headers={"Content-Type": "chemical/x-xyz"},
        )
    )

    content = client.export_file("mol_001", "xyz")

    assert "O 0.0 0.0 0.0" in content


def test_export_file_save_to(client: MolBuilder, mock_api: respx.Router, tmp_path) -> None:
    xyz_text = "3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH 0.0 1.0 0.0\n"
    mock_api.get("/molecule/mol_001/export/xyz").mock(
        return_value=httpx.Response(
            200,
            text=xyz_text,
            headers={"Content-Type": "chemical/x-xyz"},
        )
    )

    out_path = tmp_path / "output.xyz"
    result = client.export_file("mol_001", "xyz", save_to=str(out_path))

    assert result == str(out_path)
    assert out_path.read_text() == xyz_text


def test_import_file_with_format_override(client: MolBuilder, mock_api: respx.Router, tmp_path) -> None:
    route = mock_api.post("/molecule/import-file").mock(
        return_value=httpx.Response(
            200,
            json={
                "molecules": [],
                "format": "mol",
                "count": 0,
            },
        )
    )

    mol_file = tmp_path / "data.txt"
    mol_file.write_text("some mol content")

    client.import_file(str(mol_file), format="mol")

    # The filename sent should have .mol extension
    request = route.calls.last.request
    assert b"data.mol" in request.content
