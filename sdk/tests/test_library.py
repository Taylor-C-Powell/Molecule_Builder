"""Tests for library SDK methods."""

from __future__ import annotations

import httpx
import pytest
import respx

from molbuilder_client import MolBuilder, LibraryMolecule, LibraryList, LibraryImport


_MOL_JSON = {
    "id": 1,
    "smiles": "CCO",
    "name": "ethanol",
    "tags": ["alcohol"],
    "notes": "test note",
    "properties": {"molecular_weight": 46.07},
    "created_at": "2026-01-01T00:00:00",
    "updated_at": "2026-01-01T00:00:00",
}


def test_library_save(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.post("/library/").mock(
        return_value=httpx.Response(200, json=_MOL_JSON)
    )
    result = client.library_save("CCO", name="ethanol", tags=["alcohol"], notes="test note")
    assert isinstance(result, LibraryMolecule)
    assert result.id == 1
    assert result.smiles == "CCO"


def test_library_get(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/library/1").mock(
        return_value=httpx.Response(200, json=_MOL_JSON)
    )
    result = client.library_get(1)
    assert isinstance(result, LibraryMolecule)
    assert result.name == "ethanol"


def test_library_list(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/library/").mock(
        return_value=httpx.Response(200, json={
            "molecules": [_MOL_JSON],
            "total": 1,
            "page": 1,
            "per_page": 20,
        })
    )
    result = client.library_list()
    assert isinstance(result, LibraryList)
    assert result.total == 1
    assert len(result.molecules) == 1


def test_library_list_with_filters(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/library/").mock(
        return_value=httpx.Response(200, json={
            "molecules": [_MOL_JSON],
            "total": 1,
            "page": 1,
            "per_page": 10,
        })
    )
    result = client.library_list(tag="alcohol", search="ethanol", page=1, per_page=10)
    assert result.total == 1


def test_library_update(client: MolBuilder, mock_api: respx.Router) -> None:
    updated = {**_MOL_JSON, "name": "EtOH"}
    mock_api.put("/library/1").mock(
        return_value=httpx.Response(200, json=updated)
    )
    result = client.library_update(1, name="EtOH")
    assert result.name == "EtOH"


def test_library_delete(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.delete("/library/1").mock(
        return_value=httpx.Response(200, json={"status": "deleted", "id": 1})
    )
    client.library_delete(1)  # should not raise


def test_library_import(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.post("/library/import").mock(
        return_value=httpx.Response(200, json={"saved": 2, "duplicates": 1, "errors": []})
    )
    result = client.library_import(["CCO", "C", "CCO"], tag="batch1")
    assert isinstance(result, LibraryImport)
    assert result.saved == 2
    assert result.duplicates == 1


def test_library_duplicate_409(client: MolBuilder, mock_api: respx.Router) -> None:
    from molbuilder_client import MolBuilderError
    mock_api.post("/library/").mock(
        return_value=httpx.Response(409, json={"detail": "Molecule already in library"})
    )
    with pytest.raises(MolBuilderError):
        client.library_save("CCO")
