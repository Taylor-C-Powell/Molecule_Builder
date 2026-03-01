"""Tests for ADMET prediction module."""

import pytest

from molbuilder.smiles.parser import parse
from molbuilder.molecule.admet import predict_admet, ADMETProfile


def _admet(smiles: str) -> ADMETProfile:
    mol = parse(smiles)
    return predict_admet(mol)


class TestADMETBasic:
    def test_result_type(self):
        result = _admet("CCO")
        assert isinstance(result, ADMETProfile)

    def test_ethanol_good_absorption(self):
        """Ethanol: small, low logP, few HBD -> high oral bioavailability."""
        result = _admet("CCO")
        assert result.oral_bioavailability == "high"

    def test_ethanol_not_bbb_penetrant(self):
        """Ethanol: logP < 1, too hydrophilic for BBB penetration rule."""
        result = _admet("CCO")
        assert result.bbb_penetrant is False

    def test_toluene_bbb_penetrant(self):
        """Toluene: small MW, logP ~2, low TPSA -> BBB penetrant."""
        result = _admet("Cc1ccccc1")
        assert result.bbb_penetrant is True

    def test_overall_score_range(self):
        for smi in ["CCO", "c1ccccc1", "CC(=O)O", "CCCCCCCCCCCCCCCCCCCC"]:
            result = _admet(smi)
            assert 0.0 <= result.overall_score <= 10.0

    def test_cyp_inhibition_keys(self):
        """CYP inhibition dict has exactly 5 isoforms."""
        result = _admet("CCO")
        expected = {"CYP1A2", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A4"}
        assert set(result.cyp_inhibition.keys()) == expected
        assert all(isinstance(v, bool) for v in result.cyp_inhibition.values())


class TestADMETAbsorption:
    def test_large_lipophilic_poor_absorption(self):
        """Long alkane chain -- high logP, may violate Lipinski."""
        result = _admet("CCCCCCCCCCCCCCCCCCCC")
        assert result.oral_bioavailability in ("moderate", "low")


class TestADMETToxicity:
    def test_nitro_compound_ames_positive(self):
        """Nitrobenzene has nitro group -> Ames positive."""
        result = _admet("[O-][N+](=O)c1ccccc1")
        assert result.ames_mutagenicity is True

    def test_herg_risk_lipophilic_basic(self):
        """Haloperidol-like: lipophilic, basic N, aromatic -> hERG risk."""
        # Simplified haloperidol-like scaffold
        result = _admet("OC(CCCN1CCC(=O)CC1)(c1ccc(F)cc1)c1ccc(Cl)cc1")
        assert result.herg_risk in ("moderate", "high")

    def test_small_molecule_no_structural_alerts(self):
        """Ethanol should have no structural alerts."""
        result = _admet("CCO")
        assert result.structural_alerts == []


class TestADMETMisc:
    def test_warnings_and_flags_are_lists(self):
        result = _admet("CCO")
        assert isinstance(result.warnings, list)
        assert isinstance(result.flags, list)
        assert all(isinstance(w, str) for w in result.warnings)
        assert all(isinstance(f, str) for f in result.flags)
