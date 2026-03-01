"""Tests for ADMET SDK client method."""

from __future__ import annotations

import httpx
import respx

from molbuilder_client import ADMETProfile, MolBuilder


_MOCK_ADMET = {
    "id": "mol_001",
    "smiles": "CCO",
    "oral_bioavailability": "high",
    "intestinal_absorption": "high",
    "caco2_permeability": "moderate",
    "pgp_substrate": False,
    "bbb_penetrant": True,
    "plasma_protein_binding": "low",
    "vd_class": "low",
    "cyp_inhibition": {
        "CYP1A2": False,
        "CYP2C9": False,
        "CYP2C19": False,
        "CYP2D6": False,
        "CYP3A4": False,
    },
    "metabolic_stability": "high",
    "renal_clearance": "high",
    "half_life_class": "moderate",
    "herg_risk": "low",
    "ames_mutagenicity": False,
    "hepatotoxicity_risk": "low",
    "structural_alerts": [],
    "overall_score": 9.5,
    "warnings": [],
    "flags": ["Predicted BBB penetrant -- consider CNS side effects"],
}


def test_get_admet(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/molecule/mol_001/admet").mock(
        return_value=httpx.Response(200, json=_MOCK_ADMET)
    )

    result = client.get_admet("mol_001")

    assert isinstance(result, ADMETProfile)
    assert result.id == "mol_001"
    assert result.oral_bioavailability == "high"
    assert result.bbb_penetrant is True
    assert len(result.cyp_inhibition) == 5
    assert result.overall_score == 9.5
    assert isinstance(result.warnings, list)
    assert isinstance(result.flags, list)
    assert isinstance(result.structural_alerts, list)
