"""Unit tests for the molbuilder SMARTS tokenizer, parser, and matcher."""

from __future__ import annotations

import unittest

from molbuilder.smarts.tokenizer import tokenize_smarts, SmartsTokenType
from molbuilder.smarts.parser import parse_smarts, SmartsPattern
from molbuilder.smarts.matcher import find_matches, has_match, detect_by_smarts
from molbuilder.smiles.parser import parse


# ===================================================================
# Tokenizer tests
# ===================================================================

class TestSmartsTokenizer(unittest.TestCase):
    """Tests for molbuilder.smarts.tokenizer.tokenize_smarts."""

    def test_simple_organic(self):
        """tokenize_smarts('CCO') produces 3 ATOM tokens."""
        tokens = tokenize_smarts("CCO")
        self.assertEqual(len(tokens), 3)
        self.assertTrue(all(t.type == SmartsTokenType.ATOM for t in tokens))
        self.assertEqual(tokens[0].symbol, "C")
        self.assertEqual(tokens[1].symbol, "C")
        self.assertEqual(tokens[2].symbol, "O")

    def test_bracket_atomic_number(self):
        """tokenize_smarts('[#6]') produces one ATOM token with atomic_number=6."""
        tokens = tokenize_smarts("[#6]")
        self.assertEqual(len(tokens), 1)
        tok = tokens[0]
        self.assertEqual(tok.type, SmartsTokenType.ATOM)
        self.assertEqual(tok.atomic_number, 6)

    def test_bracket_degree(self):
        """tokenize_smarts('[CD3]') produces atom with symbol=C, degree=3."""
        tokens = tokenize_smarts("[CD3]")
        self.assertEqual(len(tokens), 1)
        tok = tokens[0]
        self.assertEqual(tok.symbol, "C")
        self.assertEqual(tok.degree, 3)

    def test_bracket_wildcard(self):
        """tokenize_smarts('[*]') produces atom with symbol=None (wildcard)."""
        tokens = tokenize_smarts("[*]")
        self.assertEqual(len(tokens), 1)
        tok = tokens[0]
        self.assertIsNone(tok.symbol)

    def test_bare_wildcard(self):
        """tokenize_smarts('*') produces a WILDCARD token."""
        tokens = tokenize_smarts("*")
        self.assertEqual(len(tokens), 1)
        self.assertEqual(tokens[0].type, SmartsTokenType.WILDCARD)
        self.assertIsNone(tokens[0].symbol)

    def test_not_operator(self):
        """tokenize_smarts('[!#6]') produces atom with is_not=True, atomic_number=6."""
        tokens = tokenize_smarts("[!#6]")
        self.assertEqual(len(tokens), 1)
        tok = tokens[0]
        self.assertTrue(tok.is_not)
        self.assertEqual(tok.atomic_number, 6)

    def test_or_in_brackets(self):
        """tokenize_smarts('[#6,#7]') produces token with or_specs."""
        tokens = tokenize_smarts("[#6,#7]")
        self.assertEqual(len(tokens), 1)
        tok = tokens[0]
        self.assertIsNotNone(tok.or_specs)
        self.assertEqual(len(tok.or_specs), 2)
        self.assertEqual(tok.or_specs[0].atomic_number, 6)
        self.assertEqual(tok.or_specs[1].atomic_number, 7)

    def test_bond_any(self):
        """tokenize_smarts('C~C') produces ATOM, BOND(~), ATOM."""
        tokens = tokenize_smarts("C~C")
        self.assertEqual(len(tokens), 3)
        self.assertEqual(tokens[1].type, SmartsTokenType.BOND)
        self.assertEqual(tokens[1].value, "~")

    def test_bond_double(self):
        """tokenize_smarts('C=O') produces a BOND token for '='."""
        tokens = tokenize_smarts("C=O")
        bond_tokens = [t for t in tokens if t.type == SmartsTokenType.BOND]
        self.assertEqual(len(bond_tokens), 1)
        self.assertEqual(bond_tokens[0].value, "=")

    def test_bond_triple(self):
        """tokenize_smarts('C#N') produces a BOND token for '#'."""
        tokens = tokenize_smarts("C#N")
        bond_tokens = [t for t in tokens if t.type == SmartsTokenType.BOND]
        self.assertEqual(len(bond_tokens), 1)
        self.assertEqual(bond_tokens[0].value, "#")

    def test_aromatic_atoms(self):
        """tokenize_smarts('c') produces aromatic atom token."""
        tokens = tokenize_smarts("c")
        self.assertEqual(len(tokens), 1)
        self.assertTrue(tokens[0].aromatic)
        self.assertEqual(tokens[0].symbol, "C")

    def test_two_letter_organic(self):
        """tokenize_smarts('ClC') produces Cl and C ATOM tokens."""
        tokens = tokenize_smarts("ClC")
        self.assertEqual(len(tokens), 2)
        self.assertEqual(tokens[0].symbol, "Cl")
        self.assertEqual(tokens[1].symbol, "C")

    def test_ring_digits(self):
        """tokenize_smarts('c1ccccc1') contains RING_DIGIT tokens."""
        tokens = tokenize_smarts("c1ccccc1")
        ring_tokens = [t for t in tokens if t.type == SmartsTokenType.RING_DIGIT]
        self.assertEqual(len(ring_tokens), 2)

    def test_branch_tokens(self):
        """tokenize_smarts('[CD3](=O)[OH]') contains branch open/close."""
        tokens = tokenize_smarts("[CD3](=O)[OH]")
        types = [t.type for t in tokens]
        self.assertIn(SmartsTokenType.BRANCH_OPEN, types)
        self.assertIn(SmartsTokenType.BRANCH_CLOSE, types)

    def test_bracket_hcount(self):
        """tokenize_smarts('[OH]') produces atom with hcount=1."""
        tokens = tokenize_smarts("[OH]")
        self.assertEqual(len(tokens), 1)
        tok = tokens[0]
        self.assertEqual(tok.symbol, "O")
        self.assertEqual(tok.hcount, 1)

    def test_bracket_hcount_explicit(self):
        """tokenize_smarts('[NH2]') produces atom with hcount=2."""
        tokens = tokenize_smarts("[NH2]")
        self.assertEqual(len(tokens), 1)
        self.assertEqual(tokens[0].hcount, 2)

    def test_bracket_charge(self):
        """tokenize_smarts('[O-]') produces atom with charge=-1."""
        tokens = tokenize_smarts("[O-]")
        self.assertEqual(len(tokens), 1)
        self.assertEqual(tokens[0].charge, -1)

    def test_bracket_ring_membership(self):
        """tokenize_smarts('[CR1]') produces atom with ring_membership=1."""
        tokens = tokenize_smarts("[CR1]")
        self.assertEqual(len(tokens), 1)
        self.assertEqual(tokens[0].ring_membership, 1)

    def test_bracket_ring_size(self):
        """tokenize_smarts('[Cr6]') produces atom with ring_size=6."""
        tokens = tokenize_smarts("[Cr6]")
        self.assertEqual(len(tokens), 1)
        self.assertEqual(tokens[0].ring_size, 6)

    def test_bracket_valence(self):
        """tokenize_smarts('[Cv4]') produces atom with valence=4."""
        tokens = tokenize_smarts("[Cv4]")
        self.assertEqual(len(tokens), 1)
        self.assertEqual(tokens[0].valence, 4)


# ===================================================================
# Parser tests
# ===================================================================

class TestSmartsParser(unittest.TestCase):
    """Tests for molbuilder.smarts.parser.parse_smarts."""

    def test_single_atom(self):
        """parse_smarts('[OH]') produces 1 atom, 0 bonds."""
        pat = parse_smarts("[OH]")
        self.assertEqual(len(pat.atoms), 1)
        self.assertEqual(len(pat.bonds), 0)
        self.assertEqual(pat.atoms[0].symbol, "O")
        self.assertEqual(pat.atoms[0].hcount, 1)

    def test_double_bond(self):
        """parse_smarts('C=O') produces 2 atoms, 1 bond with order 2.0."""
        pat = parse_smarts("C=O")
        self.assertEqual(len(pat.atoms), 2)
        self.assertEqual(len(pat.bonds), 1)
        self.assertEqual(pat.bonds[0].order, 2.0)

    def test_ring_closure(self):
        """parse_smarts('c1ccccc1') produces 6 atoms, 6 bonds."""
        pat = parse_smarts("c1ccccc1")
        self.assertEqual(len(pat.atoms), 6)
        self.assertEqual(len(pat.bonds), 6)

    def test_branched_pattern(self):
        """parse_smarts('[CD3](=O)[OH]') -- carboxylic acid pattern."""
        pat = parse_smarts("[CD3](=O)[OH]")
        self.assertEqual(len(pat.atoms), 3)
        self.assertEqual(len(pat.bonds), 2)
        # First atom is C with degree 3
        self.assertEqual(pat.atoms[0].symbol, "C")
        self.assertEqual(pat.atoms[0].degree, 3)
        # One bond should be double (=O)
        bond_orders = [b.order for b in pat.bonds if b.order is not None]
        self.assertIn(2.0, bond_orders)

    def test_default_bond_is_any(self):
        """Default bond in SMARTS (no explicit bond symbol) matches any."""
        pat = parse_smarts("CO")
        self.assertEqual(len(pat.bonds), 1)
        bond = pat.bonds[0]
        # Default: order is None, is_any is False -> matches any bond
        self.assertIsNone(bond.order)
        self.assertFalse(bond.is_any)

    def test_any_bond_tilde(self):
        """parse_smarts('C~O') produces a bond with is_any=True."""
        pat = parse_smarts("C~O")
        self.assertEqual(len(pat.bonds), 1)
        self.assertTrue(pat.bonds[0].is_any)

    def test_wildcard_atom(self):
        """parse_smarts('[*]') produces atom with symbol=None."""
        pat = parse_smarts("[*]")
        self.assertEqual(len(pat.atoms), 1)
        self.assertIsNone(pat.atoms[0].symbol)

    def test_bare_wildcard(self):
        """parse_smarts('*') produces atom with symbol=None."""
        pat = parse_smarts("*")
        self.assertEqual(len(pat.atoms), 1)
        self.assertIsNone(pat.atoms[0].symbol)

    def test_or_atom(self):
        """parse_smarts('[#6,#7]') produces atom with or_atoms list."""
        pat = parse_smarts("[#6,#7]")
        self.assertEqual(len(pat.atoms), 1)
        self.assertIsNotNone(pat.atoms[0].or_atoms)
        self.assertEqual(len(pat.atoms[0].or_atoms), 2)

    def test_not_atom(self):
        """parse_smarts('[!#6]') produces atom with is_not=True."""
        pat = parse_smarts("[!#6]")
        self.assertEqual(len(pat.atoms), 1)
        self.assertTrue(pat.atoms[0].is_not)
        self.assertEqual(pat.atoms[0].atomic_number, 6)

    def test_pattern_neighbors(self):
        """SmartsPattern.neighbors() returns correct indices."""
        pat = parse_smarts("C=O")
        self.assertEqual(pat.neighbors(0), [1])
        self.assertEqual(pat.neighbors(1), [0])

    def test_pattern_neighbor_bonds(self):
        """SmartsPattern.neighbor_bonds() returns bonds for an atom."""
        pat = parse_smarts("[CD3](=O)[OH]")
        bonds_0 = pat.neighbor_bonds(0)
        self.assertEqual(len(bonds_0), 2)


# ===================================================================
# Matcher tests
# ===================================================================

class TestSmartsMatcher(unittest.TestCase):
    """Tests for molbuilder.smarts.matcher."""

    def test_ethanol_oxygen(self):
        """Ethanol (CCO) should match [#8] (any oxygen)."""
        mol = parse("CCO")
        pat = parse_smarts("[#8]")
        matches = find_matches(mol, pat)
        self.assertGreaterEqual(len(matches), 1)
        # Every match maps pattern atom 0 to an oxygen
        for m in matches:
            self.assertEqual(mol.atoms[m[0]].symbol, "O")

    def test_ethanol_oh(self):
        """Ethanol (CCO) should match [OH] (oxygen with 1 H)."""
        mol = parse("CCO")
        pat = parse_smarts("[OH]")
        matches = find_matches(mol, pat)
        self.assertGreaterEqual(len(matches), 1)

    def test_benzene_aromatic(self):
        """Benzene (c1ccccc1) should match 'c' (aromatic carbon) 6 times."""
        mol = parse("c1ccccc1")
        pat = parse_smarts("c")
        matches = find_matches(mol, pat)
        # Benzene has 6 aromatic carbons
        aromatic_c_matches = [
            m for m in matches
            if mol.atoms[m[0]].symbol == "C"
        ]
        self.assertEqual(len(aromatic_c_matches), 6)

    def test_aspirin_ester(self):
        """Aspirin should match ester-like pattern [CD3](=O)[OD2]."""
        # Aspirin: CC(=O)Oc1ccccc1C(=O)O
        mol = parse("CC(=O)Oc1ccccc1C(=O)O")
        pat = parse_smarts("[#6](=[#8])[#8]")
        matches = find_matches(mol, pat)
        # Should find at least 1 match (ester and/or acid carbonyl)
        self.assertGreaterEqual(len(matches), 1)

    def test_negative_hexane_no_oh(self):
        """Hexane (CCCCCC) should NOT match [OH]."""
        mol = parse("CCCCCC")
        pat = parse_smarts("[OH]")
        matches = find_matches(mol, pat)
        self.assertEqual(len(matches), 0)

    def test_wildcard_matches_any(self):
        """Wildcard [*] should match heavy atoms in ethanol."""
        mol = parse("CCO")
        pat = parse_smarts("[*]")
        matches = find_matches(mol, pat)
        # Should match at least 3 atoms (C, C, O) -- may also match H
        self.assertGreaterEqual(len(matches), 3)

    def test_has_match_positive(self):
        """has_match returns True for ethanol + [OH]."""
        mol = parse("CCO")
        pat = parse_smarts("[OH]")
        self.assertTrue(has_match(mol, pat))

    def test_has_match_negative(self):
        """has_match returns False for hexane + [OH]."""
        mol = parse("CCCCCC")
        pat = parse_smarts("[OH]")
        self.assertFalse(has_match(mol, pat))

    def test_double_bond_pattern(self):
        """Ethene (C=C) should match C=C pattern."""
        mol = parse("C=C")
        pat = parse_smarts("C=C")
        matches = find_matches(mol, pat)
        self.assertGreaterEqual(len(matches), 1)

    def test_triple_bond_pattern(self):
        """Acetylene (C#C) should match C#C pattern."""
        mol = parse("C#C")
        pat = parse_smarts("C#C")
        matches = find_matches(mol, pat)
        self.assertGreaterEqual(len(matches), 1)

    def test_not_carbon(self):
        """[!#6] should match non-carbon atoms in ethanol."""
        mol = parse("CCO")
        pat = parse_smarts("[!#6]")
        matches = find_matches(mol, pat)
        # Should match O and H atoms (not C)
        for m in matches:
            self.assertNotEqual(mol.atoms[m[0]].symbol, "C")
        self.assertGreaterEqual(len(matches), 1)

    def test_or_atom_match(self):
        """[#6,#8] should match both C and O in ethanol."""
        mol = parse("CCO")
        pat = parse_smarts("[#6,#8]")
        matches = find_matches(mol, pat)
        matched_symbols = {mol.atoms[m[0]].symbol for m in matches}
        self.assertIn("C", matched_symbols)
        self.assertIn("O", matched_symbols)

    def test_any_bond_pattern(self):
        """C~O with any-bond should match C-O in ethanol."""
        mol = parse("CCO")
        pat = parse_smarts("C~O")
        matches = find_matches(mol, pat)
        self.assertGreaterEqual(len(matches), 1)

    def test_atomic_number_carbon(self):
        """[#6] should match all carbon atoms in propane."""
        mol = parse("CCC")
        pat = parse_smarts("[#6]")
        matches = find_matches(mol, pat)
        carbon_count = sum(
            1 for a in mol.atoms if a.symbol == "C"
        )
        carbon_matches = [
            m for m in matches if mol.atoms[m[0]].symbol == "C"
        ]
        self.assertEqual(len(carbon_matches), carbon_count)

    def test_two_atom_pattern(self):
        """[#6][#8] should match C-O connectivity in ethanol."""
        mol = parse("CCO")
        pat = parse_smarts("[#6][#8]")
        matches = find_matches(mol, pat)
        self.assertGreaterEqual(len(matches), 1)
        for m in matches:
            self.assertEqual(mol.atoms[m[0]].symbol, "C")
            self.assertEqual(mol.atoms[m[1]].symbol, "O")

    def test_carbonyl_pattern(self):
        """[#6]=[#8] should match C=O in acetone."""
        mol = parse("CC(=O)C")
        pat = parse_smarts("[#6]=[#8]")
        matches = find_matches(mol, pat)
        self.assertGreaterEqual(len(matches), 1)

    def test_empty_pattern(self):
        """Empty SMARTS pattern should return no matches."""
        mol = parse("CCO")
        pat = SmartsPattern()
        matches = find_matches(mol, pat)
        self.assertEqual(len(matches), 0)


# ===================================================================
# detect_by_smarts tests
# ===================================================================

class TestDetectBySmarts(unittest.TestCase):
    """Tests for the detect_by_smarts convenience function."""

    def test_detect_alcohol(self):
        """detect_by_smarts finds alcohol OH in ethanol."""
        mol = parse("CCO")
        groups = detect_by_smarts(mol, "[OH]", "alcohol")
        self.assertGreaterEqual(len(groups), 1)
        self.assertEqual(groups[0].name, "alcohol")

    def test_detect_returns_functional_group(self):
        """detect_by_smarts returns FunctionalGroup instances."""
        from molbuilder.reactions.functional_group_detect import FunctionalGroup
        mol = parse("CCO")
        groups = detect_by_smarts(mol, "[#8]", "oxygen")
        for g in groups:
            self.assertIsInstance(g, FunctionalGroup)

    def test_detect_no_match(self):
        """detect_by_smarts returns empty list when no match."""
        mol = parse("CCCCCC")
        groups = detect_by_smarts(mol, "[OH]", "alcohol")
        self.assertEqual(len(groups), 0)

    def test_detect_custom_name(self):
        """detect_by_smarts uses the provided group_name."""
        mol = parse("CCO")
        groups = detect_by_smarts(mol, "[#8]", "my_group")
        for g in groups:
            self.assertEqual(g.name, "my_group")

    def test_detect_smarts_like_field(self):
        """detect_by_smarts sets smarts_like to the SMARTS string."""
        mol = parse("CCO")
        groups = detect_by_smarts(mol, "[#8]", "oxygen")
        for g in groups:
            self.assertEqual(g.smarts_like, "[#8]")


# ===================================================================
# Performance test
# ===================================================================

class TestSmartsPerformance(unittest.TestCase):
    """Basic performance sanity check."""

    def test_caffeine_match(self):
        """Matching a 3-atom pattern against caffeine should complete quickly."""
        # Caffeine SMILES
        mol = parse("Cn1cnc2c1c(=O)n(c(=O)n2C)C")
        pat = parse_smarts("[#6]=[#8]")
        matches = find_matches(mol, pat)
        # Caffeine has 2 C=O groups
        self.assertGreaterEqual(len(matches), 1)

    def test_larger_pattern_match(self):
        """3-atom branched pattern on caffeine."""
        mol = parse("Cn1cnc2c1c(=O)n(c(=O)n2C)C")
        pat = parse_smarts("[#6](=[#8])[#7]")
        matches = find_matches(mol, pat)
        # Caffeine has C(=O)N motifs
        self.assertGreaterEqual(len(matches), 1)


if __name__ == "__main__":
    unittest.main()
