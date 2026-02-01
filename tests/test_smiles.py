"""Unit tests for the molbuilder SMILES tokenizer, parser, and writer."""

import unittest

from molbuilder.smiles.tokenizer import tokenize, TokenType
from molbuilder.smiles.parser import parse
from molbuilder.smiles.writer import to_smiles


class TestTokenizer(unittest.TestCase):
    """Tests for molbuilder.smiles.tokenizer.tokenize."""

    def test_tokenize_ethanol(self):
        """tokenize('CCO') should produce three ATOM tokens: C, C, O."""
        tokens = tokenize("CCO")
        self.assertEqual(len(tokens), 3)
        self.assertTrue(all(t.type == TokenType.ATOM for t in tokens))
        self.assertEqual(tokens[0].value, "C")
        self.assertEqual(tokens[1].value, "C")
        self.assertEqual(tokens[2].value, "O")

    def test_tokenize_bracket_atom(self):
        """tokenize('[NH4+]') should yield one ATOM token with correct fields."""
        tokens = tokenize("[NH4+]")
        self.assertEqual(len(tokens), 1)
        tok = tokens[0]
        self.assertEqual(tok.type, TokenType.ATOM)
        self.assertEqual(tok.value, "N")
        self.assertEqual(tok.hcount, 4)
        self.assertEqual(tok.charge, 1)

    def test_tokenize_ring_digits(self):
        """tokenize('c1ccccc1') should contain RING_DIGIT tokens for '1'."""
        tokens = tokenize("c1ccccc1")
        ring_tokens = [t for t in tokens if t.type == TokenType.RING_DIGIT]
        self.assertEqual(len(ring_tokens), 2)
        self.assertTrue(all(t.value == "1" for t in ring_tokens))

    def test_tokenize_branch_tokens(self):
        """tokenize('CC(=O)O') should contain BRANCH_OPEN and BRANCH_CLOSE."""
        tokens = tokenize("CC(=O)O")
        types = [t.type for t in tokens]
        self.assertIn(TokenType.BRANCH_OPEN, types)
        self.assertIn(TokenType.BRANCH_CLOSE, types)

    def test_tokenize_bond_symbol(self):
        """tokenize('C=C') should contain a BOND token for '='."""
        tokens = tokenize("C=C")
        bond_tokens = [t for t in tokens if t.type == TokenType.BOND]
        self.assertEqual(len(bond_tokens), 1)
        self.assertEqual(bond_tokens[0].value, "=")

    def test_tokenize_two_letter_atom(self):
        """tokenize('ClC') should produce Cl and C ATOM tokens."""
        tokens = tokenize("ClC")
        self.assertEqual(len(tokens), 2)
        self.assertEqual(tokens[0].value, "Cl")
        self.assertEqual(tokens[1].value, "C")


class TestParser(unittest.TestCase):
    """Tests for molbuilder.smiles.parser.parse."""

    def test_parse_methane(self):
        """parse('C') should yield 5 atoms: 1 C + 4 H."""
        mol = parse("C")
        symbols = [a.symbol for a in mol.atoms]
        self.assertEqual(symbols.count("C"), 1)
        self.assertEqual(symbols.count("H"), 4)
        self.assertEqual(len(mol.atoms), 5)

    def test_parse_ethanol(self):
        """parse('CCO') should yield 9 atoms: 2 C + 1 O + 6 H."""
        mol = parse("CCO")
        symbols = [a.symbol for a in mol.atoms]
        self.assertEqual(symbols.count("C"), 2)
        self.assertEqual(symbols.count("O"), 1)
        self.assertEqual(symbols.count("H"), 6)
        self.assertEqual(len(mol.atoms), 9)

    def test_parse_ethene(self):
        """parse('C=C') should yield 6 atoms: 2 C + 4 H."""
        mol = parse("C=C")
        symbols = [a.symbol for a in mol.atoms]
        self.assertEqual(symbols.count("C"), 2)
        self.assertEqual(symbols.count("H"), 4)
        self.assertEqual(len(mol.atoms), 6)

    def test_parse_benzene(self):
        """parse('c1ccccc1') should yield 12 atoms: 6 C + 6 H."""
        mol = parse("c1ccccc1")
        symbols = [a.symbol for a in mol.atoms]
        self.assertEqual(symbols.count("C"), 6)
        self.assertEqual(symbols.count("H"), 6)
        self.assertEqual(len(mol.atoms), 12)

    def test_parse_acetic_acid(self):
        """parse('CC(=O)O') should yield 8 atoms: 2C + 2O + 4H."""
        mol = parse("CC(=O)O")
        symbols = [a.symbol for a in mol.atoms]
        self.assertEqual(symbols.count("C"), 2)
        self.assertEqual(symbols.count("O"), 2)
        self.assertEqual(symbols.count("H"), 4)
        self.assertEqual(len(mol.atoms), 8)

    def test_all_atoms_have_3d_positions(self):
        """Every atom produced by parse should have a 3D position array."""
        for smi in ("C", "CCO", "C=C", "c1ccccc1", "CC(=O)O"):
            mol = parse(smi)
            for atom in mol.atoms:
                pos = atom.position
                self.assertEqual(len(pos), 3,
                                 f"Atom {atom.symbol} in {smi} missing 3D pos")


class TestWriter(unittest.TestCase):
    """Tests for molbuilder.smiles.writer.to_smiles."""

    def test_to_smiles_non_empty(self):
        """to_smiles should return a non-empty string for simple molecules."""
        for smi in ("C", "CCO", "C=C", "CC(=O)O"):
            mol = parse(smi)
            result = to_smiles(mol)
            self.assertIsInstance(result, str)
            self.assertTrue(len(result) > 0,
                            f"to_smiles returned empty for {smi}")

    def test_round_trip_methane(self):
        """parse then to_smiles then re-parse should preserve atom count."""
        mol1 = parse("C")
        smi = to_smiles(mol1)
        mol2 = parse(smi)
        self.assertEqual(len(mol1.atoms), len(mol2.atoms))

    def test_round_trip_ethanol(self):
        """Round-trip CCO: atom count should be preserved."""
        mol1 = parse("CCO")
        smi = to_smiles(mol1)
        mol2 = parse(smi)
        self.assertEqual(len(mol1.atoms), len(mol2.atoms))

    def test_round_trip_ethene(self):
        """Round-trip C=C: atom count should be preserved."""
        mol1 = parse("C=C")
        smi = to_smiles(mol1)
        mol2 = parse(smi)
        self.assertEqual(len(mol1.atoms), len(mol2.atoms))


class TestErrorHandling(unittest.TestCase):
    """Tests that invalid SMILES raise ValueError."""

    def test_invalid_character(self):
        """An unexpected character in SMILES should raise ValueError."""
        with self.assertRaises(ValueError):
            parse("C&C")

    def test_unclosed_ring(self):
        """An unclosed ring digit should raise ValueError."""
        with self.assertRaises(ValueError):
            parse("C1CC")

    def test_unclosed_bracket(self):
        """An unclosed bracket atom should raise an error."""
        with self.assertRaises((ValueError, Exception)):
            parse("[NH4")


if __name__ == "__main__":
    unittest.main()
