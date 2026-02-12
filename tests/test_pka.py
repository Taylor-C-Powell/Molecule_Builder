"""Tests for pKa prediction functionality."""

import pytest

from molbuilder.smiles.parser import parse
from molbuilder.molecule.properties import predict_pka, pKaPrediction


# =====================================================================
#  Helper
# =====================================================================

def mol(smiles: str):
    """Parse a SMILES string into a Molecule."""
    return parse(smiles)


# =====================================================================
#  pKaPrediction dataclass
# =====================================================================

class TestpKaPredictionDataclass:
    def test_fields(self):
        pred = pKaPrediction(
            group_name="carboxylic acid",
            atom_index=3,
            pka_value=4.8,
            acidic=True,
        )
        assert pred.group_name == "carboxylic acid"
        assert pred.atom_index == 3
        assert pred.pka_value == 4.8
        assert pred.acidic is True

    def test_amine_not_acidic(self):
        pred = pKaPrediction(
            group_name="primary amine",
            atom_index=0,
            pka_value=10.6,
            acidic=False,
        )
        assert pred.acidic is False


# =====================================================================
#  Acetic acid (CH3COOH) ~ 4.76
# =====================================================================

class TestAceticAcid:
    def test_finds_carboxylic_acid(self):
        results = predict_pka(mol("CC(=O)O"))
        groups = [r.group_name for r in results]
        assert "carboxylic acid" in groups

    def test_pka_near_4_76(self):
        results = predict_pka(mol("CC(=O)O"))
        acid = [r for r in results if r.group_name == "carboxylic acid"][0]
        # Base pKa is 4.76; the C=O on the same carbon may cause a small
        # correction, but the O-H oxygen's own carbonyl is the defining
        # feature, not an extra EWG.  Accept 3.0 - 6.0.
        assert 3.0 <= acid.pka_value <= 6.0, f"acetic acid pKa={acid.pka_value}"

    def test_acidic_flag(self):
        results = predict_pka(mol("CC(=O)O"))
        acid = [r for r in results if r.group_name == "carboxylic acid"][0]
        assert acid.acidic is True


# =====================================================================
#  Benzoic acid ~ 4.2 (aromatic correction)
# =====================================================================

class TestBenzoicAcid:
    def test_finds_carboxylic_acid(self):
        results = predict_pka(mol("c1ccc(cc1)C(=O)O"))
        groups = [r.group_name for r in results]
        assert "carboxylic acid" in groups

    def test_pka_lower_than_acetic(self):
        acetic = predict_pka(mol("CC(=O)O"))
        benzoic = predict_pka(mol("c1ccc(cc1)C(=O)O"))
        acetic_pka = [r for r in acetic if r.group_name == "carboxylic acid"][0].pka_value
        benzoic_pka = [r for r in benzoic if r.group_name == "carboxylic acid"][0].pka_value
        assert benzoic_pka < acetic_pka, (
            f"benzoic ({benzoic_pka}) should be more acidic than acetic ({acetic_pka})"
        )

    def test_pka_in_range(self):
        results = predict_pka(mol("c1ccc(cc1)C(=O)O"))
        acid = [r for r in results if r.group_name == "carboxylic acid"][0]
        # Expect roughly 3.0 - 5.0 with aromatic correction
        assert 2.0 <= acid.pka_value <= 5.5, f"benzoic acid pKa={acid.pka_value}"


# =====================================================================
#  Phenol ~ 10.0
# =====================================================================

class TestPhenol:
    def test_finds_phenol(self):
        results = predict_pka(mol("c1ccc(cc1)O"))
        groups = [r.group_name for r in results]
        assert "phenol" in groups

    def test_pka_near_10(self):
        results = predict_pka(mol("c1ccc(cc1)O"))
        phenol = [r for r in results if r.group_name == "phenol"][0]
        assert 8.0 <= phenol.pka_value <= 11.0, f"phenol pKa={phenol.pka_value}"

    def test_acidic_flag(self):
        results = predict_pka(mol("c1ccc(cc1)O"))
        phenol = [r for r in results if r.group_name == "phenol"][0]
        assert phenol.acidic is True


# =====================================================================
#  Ethylamine (CCN) ~ 10.6
# =====================================================================

class TestEthylamine:
    def test_finds_primary_amine(self):
        results = predict_pka(mol("CCN"))
        groups = [r.group_name for r in results]
        assert "primary amine" in groups

    def test_pka_near_10_6(self):
        results = predict_pka(mol("CCN"))
        amine = [r for r in results if r.group_name == "primary amine"][0]
        assert 9.0 <= amine.pka_value <= 12.0, f"ethylamine pKa={amine.pka_value}"

    def test_not_acidic(self):
        results = predict_pka(mol("CCN"))
        amine = [r for r in results if r.group_name == "primary amine"][0]
        assert amine.acidic is False


# =====================================================================
#  Ethanol (CCO) ~ 16.0
# =====================================================================

class TestEthanol:
    def test_finds_alcohol(self):
        results = predict_pka(mol("CCO"))
        groups = [r.group_name for r in results]
        assert "alcohol" in groups

    def test_pka_near_16(self):
        results = predict_pka(mol("CCO"))
        alcohol = [r for r in results if r.group_name == "alcohol"][0]
        assert 14.0 <= alcohol.pka_value <= 18.0, f"ethanol pKa={alcohol.pka_value}"

    def test_acidic_flag(self):
        results = predict_pka(mol("CCO"))
        alcohol = [r for r in results if r.group_name == "alcohol"][0]
        assert alcohol.acidic is True


# =====================================================================
#  Glycine (NCC(=O)O) - both amine and carboxylic acid
# =====================================================================

class TestGlycine:
    def test_finds_both_groups(self):
        results = predict_pka(mol("NCC(=O)O"))
        groups = {r.group_name for r in results}
        assert "primary amine" in groups
        assert "carboxylic acid" in groups

    def test_two_predictions(self):
        results = predict_pka(mol("NCC(=O)O"))
        assert len(results) >= 2

    def test_acid_more_acidic_than_amine(self):
        results = predict_pka(mol("NCC(=O)O"))
        acid = [r for r in results if r.group_name == "carboxylic acid"][0]
        amine = [r for r in results if r.group_name == "primary amine"][0]
        assert acid.pka_value < amine.pka_value

    def test_sorted_by_pka(self):
        results = predict_pka(mol("NCC(=O)O"))
        pkas = [r.pka_value for r in results]
        assert pkas == sorted(pkas)


# =====================================================================
#  Methane - no ionizable groups
# =====================================================================

class TestMethane:
    def test_empty_list(self):
        results = predict_pka(mol("C"))
        assert results == []

    def test_returns_list_type(self):
        results = predict_pka(mol("C"))
        assert isinstance(results, list)


# =====================================================================
#  Additional edge cases
# =====================================================================

class TestEdgeCases:
    def test_thiol(self):
        # Ethanethiol (CCS)
        results = predict_pka(mol("CCS"))
        groups = [r.group_name for r in results]
        assert "thiol" in groups
        thiol = [r for r in results if r.group_name == "thiol"][0]
        assert 9.0 <= thiol.pka_value <= 12.0, f"thiol pKa={thiol.pka_value}"

    def test_secondary_amine(self):
        # Diethylamine: CCNCC
        results = predict_pka(mol("CCNCC"))
        groups = [r.group_name for r in results]
        assert "secondary amine" in groups

    def test_results_sorted(self):
        # Glycine has multiple groups; verify global sort
        results = predict_pka(mol("NCC(=O)O"))
        for i in range(len(results) - 1):
            assert results[i].pka_value <= results[i + 1].pka_value

    def test_pka_values_rounded(self):
        # All predicted pKa values should be rounded to 1 decimal place
        results = predict_pka(mol("NCC(=O)O"))
        for r in results:
            assert r.pka_value == round(r.pka_value, 1)

    def test_benzene_no_ionizable(self):
        results = predict_pka(mol("c1ccccc1"))
        assert results == []

    def test_water(self):
        results = predict_pka(mol("O"))
        # Water is essentially an alcohol (O-H)
        assert len(results) >= 1
        groups = [r.group_name for r in results]
        assert "alcohol" in groups
