"""Tests for molecule endpoints: from_smiles, get_properties, get_3d."""

from __future__ import annotations

import httpx
import respx

from molbuilder_client import (
    Atom3D,
    Bond3D,
    MolBuilder,
    Molecule3D,
    MoleculeInfo,
    MoleculeProperties,
)


def test_from_smiles(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.post("/molecule/from-smiles").mock(
        return_value=httpx.Response(
            200,
            json={
                "id": "mol_001",
                "name": "ethanol",
                "smiles": "CCO",
                "num_atoms": 9,
                "num_bonds": 8,
            },
        )
    )

    result = client.from_smiles("CCO", name="ethanol")

    assert isinstance(result, MoleculeInfo)
    assert result.id == "mol_001"
    assert result.name == "ethanol"
    assert result.smiles == "CCO"
    assert result.num_atoms == 9
    assert result.num_bonds == 8


def test_from_smiles_default_name(client: MolBuilder, mock_api: respx.Router) -> None:
    route = mock_api.post("/molecule/from-smiles").mock(
        return_value=httpx.Response(
            200,
            json={"id": "m1", "name": "", "smiles": "C", "num_atoms": 5, "num_bonds": 4},
        )
    )

    client.from_smiles("C")

    # When name is empty string, it should not be included in the payload
    import json
    body = json.loads(route.calls.last.request.content)
    assert "name" not in body


def test_get_properties(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/molecule/mol_001/properties").mock(
        return_value=httpx.Response(
            200,
            json={
                "id": "mol_001",
                "smiles": "CCO",
                "formula": "C2H6O",
                "molecular_weight": 46.07,
                "num_atoms": 9,
                "num_bonds": 8,
                "functional_groups": ["alcohol"],
            },
        )
    )

    result = client.get_properties("mol_001")

    assert isinstance(result, MoleculeProperties)
    assert result.formula == "C2H6O"
    assert result.molecular_weight == 46.07
    assert result.functional_groups == ["alcohol"]


def test_get_3d(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/molecule/mol_001/3d").mock(
        return_value=httpx.Response(
            200,
            json={
                "id": "mol_001",
                "atoms": [
                    {
                        "index": 0,
                        "symbol": "C",
                        "position": [0.0, 0.0, 0.0],
                        "hybridization": "sp3",
                        "formal_charge": 0,
                    },
                    {
                        "index": 1,
                        "symbol": "O",
                        "position": [1.43, 0.0, 0.0],
                        "hybridization": "sp3",
                        "formal_charge": 0,
                    },
                ],
                "bonds": [
                    {"atom_i": 0, "atom_j": 1, "order": 1.0, "rotatable": True},
                ],
            },
        )
    )

    result = client.get_3d("mol_001")

    assert isinstance(result, Molecule3D)
    assert len(result.atoms) == 2
    assert isinstance(result.atoms[0], Atom3D)
    assert result.atoms[0].symbol == "C"
    assert result.atoms[0].position == (0.0, 0.0, 0.0)
    assert result.atoms[1].hybridization == "sp3"
    assert len(result.bonds) == 1
    assert isinstance(result.bonds[0], Bond3D)
    assert result.bonds[0].rotatable is True


def test_get_3d_null_hybridization(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/molecule/mol_002/3d").mock(
        return_value=httpx.Response(
            200,
            json={
                "id": "mol_002",
                "atoms": [
                    {
                        "index": 0,
                        "symbol": "H",
                        "position": [0.0, 0.0, 0.0],
                        "hybridization": None,
                        "formal_charge": 0,
                    },
                ],
                "bonds": [],
            },
        )
    )

    result = client.get_3d("mol_002")
    assert result.atoms[0].hybridization is None
