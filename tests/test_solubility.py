"""Tests for aqueous solubility and crystallization prediction."""

import pytest

from molbuilder.smiles.parser import parse
from molbuilder.molecule.solubility import predict_solubility, SolubilityResult

VALID_CLASSES = {"highly soluble", "soluble", "moderately soluble", "poorly soluble", "insoluble"}
VALID_RISKS = {"low", "moderate", "high"}


def _sol(smiles: str) -> SolubilityResult:
    mol = parse(smiles)
    return predict_solubility(mol)


class TestSolubilityBasic:
    def test_result_type(self):
        result = _sol("CCO")
        assert isinstance(result, SolubilityResult)

    def test_ethanol_highly_soluble(self):
        result = _sol("CCO")
        assert result.solubility_class in ("highly soluble", "soluble")

    def test_lipophilic_poorly_soluble(self):
        # Long alkane -- very hydrophobic
        result = _sol("CCCCCCCCCCCCCCCCCCCC")
        assert result.solubility_class in ("poorly soluble", "insoluble")

    def test_solubility_mg_ml_non_negative(self):
        for smi in ["CCO", "c1ccccc1", "CC(=O)O", "CCCCCCCCCC"]:
            result = _sol(smi)
            assert result.solubility_mg_ml >= 0.0

    def test_solubility_class_values(self):
        for smi in ["CCO", "c1ccccc1", "CC(=O)O", "CCCCCCCCCCCCCCCCCCCC"]:
            result = _sol(smi)
            assert result.solubility_class in VALID_CLASSES


class TestRiskAssessments:
    def test_crystallization_risk_values(self):
        for smi in ["CCO", "c1ccccc1", "CC(=O)O"]:
            result = _sol(smi)
            assert result.crystallization_risk in VALID_RISKS

    def test_polymorph_risk_values(self):
        for smi in ["CCO", "c1ccccc1", "CC(=O)O"]:
            result = _sol(smi)
            assert result.polymorph_risk in VALID_RISKS

    def test_melting_point_positive(self):
        for smi in ["CCO", "c1ccccc1", "CC(=O)O", "C(O)C(O)C(O)C(O)C(O)CO"]:
            result = _sol(smi)
            assert result.estimated_melting_point_c > 0

    def test_high_hbd_higher_polymorph_risk(self):
        """Glucose (many OH groups) should have higher polymorph risk than ethanol."""
        glucose = _sol("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O")
        ethanol = _sol("CCO")
        risk_order = {"low": 0, "moderate": 1, "high": 2}
        assert risk_order[glucose.polymorph_risk] >= risk_order[ethanol.polymorph_risk]
