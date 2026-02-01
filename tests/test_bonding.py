"""Tests for the molbuilder.bonding package.

Covers:
    - lewis: formula parsing, LewisStructure steric_number() and lone_pairs_on_central()
    - vsepr: VSEPRMolecule geometry via .axe.molecular_geometry and .axe.hybridization
    - covalent: bond constructors, polarity classification, orbital composition
"""

import unittest

from molbuilder.bonding.lewis import parse_formula, LewisStructure
from molbuilder.bonding.vsepr import VSEPRMolecule
from molbuilder.bonding.covalent import (
    CovalentBond,
    BondPolarity,
    OrbitalType,
    single_bond,
    double_bond,
    triple_bond,
    sigma_pi_composition,
    MolecularBondAnalysis,
)


# ===================================================================
# Lewis structure tests
# ===================================================================

class TestParseFormula(unittest.TestCase):
    """Tests for the parse_formula helper."""

    def test_h2o(self):
        atoms = parse_formula("H2O")
        self.assertEqual(atoms.count("H"), 2)
        self.assertEqual(atoms.count("O"), 1)
        self.assertEqual(len(atoms), 3)

    def test_ch4(self):
        atoms = parse_formula("CH4")
        self.assertEqual(atoms.count("C"), 1)
        self.assertEqual(atoms.count("H"), 4)
        self.assertEqual(len(atoms), 5)

    def test_sf6(self):
        atoms = parse_formula("SF6")
        self.assertEqual(atoms.count("S"), 1)
        self.assertEqual(atoms.count("F"), 6)
        self.assertEqual(len(atoms), 7)

    def test_co2(self):
        atoms = parse_formula("CO2")
        self.assertEqual(atoms.count("C"), 1)
        self.assertEqual(atoms.count("O"), 2)

    def test_single_element_no_count(self):
        atoms = parse_formula("O")
        self.assertEqual(atoms, ["O"])

    def test_unknown_element_raises(self):
        with self.assertRaises(ValueError):
            parse_formula("Xx2")


class TestLewisStructureH2O(unittest.TestCase):
    """Water: central O, 2 bonding pairs, 2 lone pairs, steric 4."""

    def setUp(self):
        self.lewis = LewisStructure("H2O")

    def test_central_atom(self):
        self.assertEqual(self.lewis.central_symbol, "O")

    def test_steric_number(self):
        self.assertEqual(self.lewis.steric_number(), 4)

    def test_lone_pairs_on_central(self):
        self.assertEqual(self.lewis.lone_pairs_on_central(), 2)

    def test_bonding_pairs_on_central(self):
        self.assertEqual(self.lewis.bonding_pairs_on_central(), 2)


class TestLewisStructureCO2(unittest.TestCase):
    """Carbon dioxide: central C, 2 double bonds, 0 lone pairs, steric 2."""

    def setUp(self):
        self.lewis = LewisStructure("CO2")

    def test_central_atom(self):
        self.assertEqual(self.lewis.central_symbol, "C")

    def test_steric_number(self):
        self.assertEqual(self.lewis.steric_number(), 2)

    def test_lone_pairs_on_central(self):
        self.assertEqual(self.lewis.lone_pairs_on_central(), 0)

    def test_double_bonds(self):
        for bond in self.lewis.bonds:
            self.assertEqual(bond.order, 2)


class TestLewisStructureNH3(unittest.TestCase):
    """Ammonia: central N, 3 bonding pairs, 1 lone pair, steric 4."""

    def setUp(self):
        self.lewis = LewisStructure("NH3")

    def test_central_atom(self):
        self.assertEqual(self.lewis.central_symbol, "N")

    def test_steric_number(self):
        self.assertEqual(self.lewis.steric_number(), 4)

    def test_lone_pairs_on_central(self):
        self.assertEqual(self.lewis.lone_pairs_on_central(), 1)

    def test_bonding_pairs_on_central(self):
        self.assertEqual(self.lewis.bonding_pairs_on_central(), 3)


class TestLewisStructureSF6(unittest.TestCase):
    """SF6: expanded octet, central S, 6 bonding pairs, 0 lone pairs, steric 6."""

    def setUp(self):
        self.lewis = LewisStructure("SF6")

    def test_central_atom(self):
        self.assertEqual(self.lewis.central_symbol, "S")

    def test_steric_number(self):
        self.assertEqual(self.lewis.steric_number(), 6)

    def test_lone_pairs_on_central(self):
        self.assertEqual(self.lewis.lone_pairs_on_central(), 0)

    def test_bonding_pairs_on_central(self):
        self.assertEqual(self.lewis.bonding_pairs_on_central(), 6)

    def test_all_single_bonds(self):
        for bond in self.lewis.bonds:
            self.assertEqual(bond.order, 1)


class TestLewisStructureCH4(unittest.TestCase):
    """Methane: central C, steric 4, 0 lone pairs."""

    def setUp(self):
        self.lewis = LewisStructure("CH4")

    def test_steric_number(self):
        self.assertEqual(self.lewis.steric_number(), 4)

    def test_lone_pairs_on_central(self):
        self.assertEqual(self.lewis.lone_pairs_on_central(), 0)


# ===================================================================
# VSEPR tests
# ===================================================================

class TestVSEPRMoleculeGeometry(unittest.TestCase):
    """Test molecular geometry classification via .axe.molecular_geometry."""

    def test_h2o_bent(self):
        mol = VSEPRMolecule("H2O")
        self.assertEqual(mol.axe.molecular_geometry, "bent")

    def test_ch4_tetrahedral(self):
        mol = VSEPRMolecule("CH4")
        self.assertEqual(mol.axe.molecular_geometry, "tetrahedral")

    def test_co2_linear(self):
        mol = VSEPRMolecule("CO2")
        self.assertEqual(mol.axe.molecular_geometry, "linear")

    def test_nh3_trigonal_pyramidal(self):
        mol = VSEPRMolecule("NH3")
        self.assertEqual(mol.axe.molecular_geometry, "trigonal pyramidal")

    def test_bf3_trigonal_planar(self):
        mol = VSEPRMolecule("BF3")
        self.assertEqual(mol.axe.molecular_geometry, "trigonal planar")

    def test_sf6_octahedral(self):
        mol = VSEPRMolecule("SF6")
        self.assertEqual(mol.axe.molecular_geometry, "octahedral")


class TestVSEPRMoleculeHybridization(unittest.TestCase):
    """Test hybridization via .axe.hybridization."""

    def test_ch4_sp3(self):
        mol = VSEPRMolecule("CH4")
        self.assertEqual(mol.axe.hybridization, "sp3")

    def test_co2_sp(self):
        mol = VSEPRMolecule("CO2")
        self.assertEqual(mol.axe.hybridization, "sp")

    def test_bf3_sp2(self):
        mol = VSEPRMolecule("BF3")
        self.assertEqual(mol.axe.hybridization, "sp2")

    def test_h2o_sp3(self):
        mol = VSEPRMolecule("H2O")
        self.assertEqual(mol.axe.hybridization, "sp3")

    def test_sf6_sp3d2(self):
        mol = VSEPRMolecule("SF6")
        self.assertEqual(mol.axe.hybridization, "sp3d2")

    def test_nh3_sp3(self):
        mol = VSEPRMolecule("NH3")
        self.assertEqual(mol.axe.hybridization, "sp3")


class TestVSEPRMoleculeAXE(unittest.TestCase):
    """Test AXE notation and steric number on the .axe classification."""

    def test_h2o_axe_notation(self):
        mol = VSEPRMolecule("H2O")
        self.assertEqual(mol.axe.axe_notation, "AX2E2")

    def test_ch4_axe_notation(self):
        mol = VSEPRMolecule("CH4")
        self.assertEqual(mol.axe.axe_notation, "AX4")

    def test_nh3_axe_notation(self):
        mol = VSEPRMolecule("NH3")
        self.assertEqual(mol.axe.axe_notation, "AX3E1")

    def test_sf6_steric_number(self):
        mol = VSEPRMolecule("SF6")
        self.assertEqual(mol.axe.steric_number, 6)

    def test_co2_electron_geometry(self):
        mol = VSEPRMolecule("CO2")
        self.assertEqual(mol.axe.electron_geometry, "linear")


class TestVSEPRCoordinates(unittest.TestCase):
    """Test that 3D coordinates are generated."""

    def test_coordinates_present(self):
        mol = VSEPRMolecule("H2O")
        coords = mol.coordinates
        self.assertIn("atom_positions", coords)
        self.assertIn("bonds", coords)
        self.assertEqual(len(coords["atom_positions"]), 3)

    def test_computed_bond_angles(self):
        mol = VSEPRMolecule("CH4")
        angles = mol.computed_bond_angles()
        self.assertTrue(len(angles) > 0)
        # All angles should be near 109.5 for tetrahedral
        for a in angles:
            self.assertAlmostEqual(a, 109.5, delta=1.0)


# ===================================================================
# Covalent bond tests
# ===================================================================

class TestCovalentBondConstructors(unittest.TestCase):
    """Test the convenience constructors single_bond, double_bond, triple_bond."""

    def test_single_bond(self):
        bond = single_bond("C", "H")
        self.assertEqual(bond.bond_order, 1)
        self.assertEqual(bond.symbol_a, "C")
        self.assertEqual(bond.symbol_b, "H")

    def test_double_bond(self):
        bond = double_bond("C", "O")
        self.assertEqual(bond.bond_order, 2)
        self.assertEqual(bond.symbol_a, "C")
        self.assertEqual(bond.symbol_b, "O")

    def test_triple_bond(self):
        bond = triple_bond("N", "N")
        self.assertEqual(bond.bond_order, 3)
        self.assertEqual(bond.symbol_a, "N")
        self.assertEqual(bond.symbol_b, "N")

    def test_invalid_element_raises(self):
        with self.assertRaises(ValueError):
            CovalentBond("X", "Y", bond_order=1)

    def test_invalid_bond_order_raises(self):
        with self.assertRaises(ValueError):
            CovalentBond("C", "H", bond_order=4)


class TestBondPolarity(unittest.TestCase):
    """Test polarity classification based on electronegativity difference."""

    def test_nonpolar_hh(self):
        bond = single_bond("H", "H")
        self.assertEqual(bond.polarity, BondPolarity.NONPOLAR_COVALENT)

    def test_nonpolar_cc(self):
        bond = single_bond("C", "C")
        self.assertEqual(bond.polarity, BondPolarity.NONPOLAR_COVALENT)

    def test_polar_co(self):
        bond = single_bond("C", "O")
        self.assertEqual(bond.polarity, BondPolarity.POLAR_COVALENT)

    def test_polar_oh(self):
        bond = single_bond("O", "H")
        self.assertEqual(bond.polarity, BondPolarity.POLAR_COVALENT)

    def test_ionic_naf(self):
        bond = single_bond("Na", "F")
        self.assertEqual(bond.polarity, BondPolarity.IONIC)

    def test_delta_en_symmetric(self):
        bond = single_bond("C", "C")
        self.assertAlmostEqual(bond.delta_en, 0.0)


class TestOrbitalComposition(unittest.TestCase):
    """Test sigma and pi orbital counts for different bond orders."""

    def test_single_bond_sigma_pi(self):
        bond = single_bond("C", "H")
        self.assertEqual(bond.sigma_bonds, 1)
        self.assertEqual(bond.pi_bonds, 0)

    def test_double_bond_sigma_pi(self):
        bond = double_bond("C", "O")
        self.assertEqual(bond.sigma_bonds, 1)
        self.assertEqual(bond.pi_bonds, 1)

    def test_triple_bond_sigma_pi(self):
        bond = triple_bond("N", "N")
        self.assertEqual(bond.sigma_bonds, 1)
        self.assertEqual(bond.pi_bonds, 2)

    def test_sigma_pi_composition_single(self):
        orbitals = sigma_pi_composition(1)
        self.assertEqual(len(orbitals), 1)
        self.assertEqual(orbitals[0].orbital_type, OrbitalType.SIGMA)

    def test_sigma_pi_composition_double(self):
        orbitals = sigma_pi_composition(2)
        self.assertEqual(len(orbitals), 2)
        sigma_count = sum(1 for o in orbitals if o.orbital_type == OrbitalType.SIGMA)
        pi_count = sum(1 for o in orbitals if o.orbital_type == OrbitalType.PI)
        self.assertEqual(sigma_count, 1)
        self.assertEqual(pi_count, 1)

    def test_sigma_pi_composition_triple(self):
        orbitals = sigma_pi_composition(3)
        self.assertEqual(len(orbitals), 3)
        sigma_count = sum(1 for o in orbitals if o.orbital_type == OrbitalType.SIGMA)
        pi_count = sum(1 for o in orbitals if o.orbital_type == OrbitalType.PI)
        self.assertEqual(sigma_count, 1)
        self.assertEqual(pi_count, 2)

    def test_invalid_order_raises(self):
        with self.assertRaises(ValueError):
            sigma_pi_composition(0)
        with self.assertRaises(ValueError):
            sigma_pi_composition(4)


class TestCovalentBondProperties(unittest.TestCase):
    """Test computed properties on CovalentBond."""

    def test_bond_length_positive(self):
        bond = single_bond("C", "H")
        self.assertGreater(bond.bond_length_pm, 0)
        self.assertGreater(bond.bond_length_angstrom, 0)

    def test_bond_length_angstrom_conversion(self):
        bond = single_bond("C", "C")
        self.assertAlmostEqual(bond.bond_length_angstrom,
                               bond.bond_length_pm / 100.0)

    def test_partial_charges_polar(self):
        bond = single_bond("C", "O")
        # O is more electronegative, so C is partial positive
        self.assertEqual(bond.partial_positive, "C")
        self.assertEqual(bond.partial_negative, "O")

    def test_percent_ionic_zero_for_same_element(self):
        bond = single_bond("C", "C")
        self.assertAlmostEqual(bond.percent_ionic_character, 0.0)

    def test_dipole_moment_zero_for_nonpolar(self):
        bond = single_bond("C", "C")
        self.assertAlmostEqual(bond.dipole_moment_debye, 0.0)

    def test_orbital_contributions_length(self):
        bond = double_bond("C", "C")
        self.assertEqual(len(bond.orbital_contributions), 2)


class TestMolecularBondAnalysis(unittest.TestCase):
    """Test the MolecularBondAnalysis class."""

    def test_co2_sigma_pi_counts(self):
        analysis = MolecularBondAnalysis("CO2")
        self.assertEqual(analysis.total_sigma_bonds, 2)
        self.assertEqual(analysis.total_pi_bonds, 2)

    def test_ch4_all_sigma(self):
        analysis = MolecularBondAnalysis("CH4")
        self.assertEqual(analysis.total_sigma_bonds, 4)
        self.assertEqual(analysis.total_pi_bonds, 0)

    def test_h2o_has_lone_pairs(self):
        analysis = MolecularBondAnalysis("H2O")
        self.assertTrue(analysis.has_lone_pairs_on_central)

    def test_ch4_no_lone_pairs(self):
        analysis = MolecularBondAnalysis("CH4")
        self.assertFalse(analysis.has_lone_pairs_on_central)


if __name__ == "__main__":
    unittest.main()
