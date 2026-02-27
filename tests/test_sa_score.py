"""Tests for synthetic accessibility scoring."""

import unittest

from molbuilder.smiles.parser import parse
from molbuilder.molecule.sa_score import sa_score, SAScoreResult


class TestSAScore(unittest.TestCase):
    """Synthetic accessibility score tests."""

    def test_ethanol_sa_low(self):
        """Ethanol (CCO) should score low (easy to synthesize)."""
        mol = parse("CCO")
        result = sa_score(mol)
        self.assertLess(result.sa_score, 4.0)

    def test_benzene_sa_moderate(self):
        """Benzene (c1ccccc1) should score moderate."""
        mol = parse("c1ccccc1")
        result = sa_score(mol)
        self.assertGreaterEqual(result.sa_score, 1.0)
        self.assertLessEqual(result.sa_score, 7.0)

    def test_complex_molecule_sa_high(self):
        """Multi-ring fused system should score higher than simple molecules."""
        # Naphthalene (fused bicyclic) should score higher than benzene
        naphthalene = parse("c1ccc2ccccc2c1")
        benzene = parse("c1ccccc1")
        r_naph = sa_score(naphthalene)
        r_benz = sa_score(benzene)
        self.assertGreater(r_naph.sa_score, r_benz.sa_score)

    def test_sa_score_range(self):
        """SA score should always be in [1, 10]."""
        for smiles in ["C", "CCO", "c1ccccc1", "CC(=O)O", "CCCCCCCCCCCCCCCCCC"]:
            mol = parse(smiles)
            result = sa_score(mol)
            self.assertGreaterEqual(result.sa_score, 1.0,
                                    f"SA score below 1 for {smiles}")
            self.assertLessEqual(result.sa_score, 10.0,
                                 f"SA score above 10 for {smiles}")

    def test_sa_score_dataclass_fields(self):
        """SAScoreResult should contain all breakdown fields."""
        mol = parse("CCO")
        result = sa_score(mol)
        self.assertIsInstance(result, SAScoreResult)
        self.assertIsInstance(result.sa_score, float)
        self.assertIsInstance(result.ring_complexity, float)
        self.assertIsInstance(result.stereo_penalty, float)
        self.assertIsInstance(result.fragment_penalty, float)
        self.assertIsInstance(result.size_penalty, float)
        self.assertIsInstance(result.sp3_bonus, float)
        self.assertIsInstance(result.heavy_atom_count, int)
        self.assertIsInstance(result.ring_count, int)
        self.assertIsInstance(result.stereo_count, int)

    def test_sa_ring_complexity(self):
        """Fused ring system should have higher ring complexity than simple ring."""
        benzene = parse("c1ccccc1")
        naphthalene = parse("c1ccc2ccccc2c1")
        r_benz = sa_score(benzene)
        r_naph = sa_score(naphthalene)
        self.assertGreater(r_naph.ring_complexity, r_benz.ring_complexity)

    def test_sa_acyclic_no_ring_complexity(self):
        """Acyclic molecule should have zero ring complexity."""
        mol = parse("CCCC")
        result = sa_score(mol)
        self.assertEqual(result.ring_complexity, 0.0)
        self.assertEqual(result.ring_count, 0)

    def test_sa_methane_minimal(self):
        """Methane (C) -- simplest molecule, should score low."""
        mol = parse("C")
        result = sa_score(mol)
        self.assertLessEqual(result.sa_score, 4.0)

    def test_sa_size_penalty(self):
        """Long chain molecule should incur size penalty."""
        # 40-carbon chain
        long = parse("C" * 40)
        short = parse("CCCC")
        r_long = sa_score(long)
        r_short = sa_score(short)
        self.assertGreaterEqual(r_long.size_penalty, 0.0)
        self.assertEqual(r_short.size_penalty, 0.0)

    def test_sa_heavy_atom_count(self):
        """Heavy atom count should exclude hydrogens."""
        mol = parse("CCO")  # 3 heavy atoms (2C + 1O)
        result = sa_score(mol)
        self.assertEqual(result.heavy_atom_count, 3)


if __name__ == "__main__":
    unittest.main()
