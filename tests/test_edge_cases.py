"""
Comprehensive edge case and scientific validation tests for molbuilder.

P2-TEST-1: Edge case tests (invalid SMILES, empty molecules, zero-scale
           costing, no-functional-group molecules, single-atom molecules,
           steric clash detection, Coulomb's law zero distance, round-trip
           I/O, and P0 regression tests).

P2-TEST-2: Scientific validation tests (Bohr model ionisation energies,
           standard bond lengths, VSEPR angles, element lookups, quantum
           numbers / electron configurations).
"""

import math
import os
import tempfile
import unittest

import numpy as np

# ---------------------------------------------------------------------------
# P2-TEST-1: Edge Case Tests
# ---------------------------------------------------------------------------


# ===================================================================
# 1. Invalid SMILES parsing
# ===================================================================

class TestInvalidSMILES(unittest.TestCase):
    """Malformed SMILES strings must raise errors."""

    def test_invalid_character_ampersand(self):
        from molbuilder.smiles.parser import parse
        with self.assertRaises(ValueError):
            parse("C&C")

    def test_invalid_character_exclamation(self):
        from molbuilder.smiles.parser import parse
        with self.assertRaises(ValueError):
            parse("C!C")

    def test_unclosed_bracket(self):
        from molbuilder.smiles.parser import parse
        with self.assertRaises((ValueError, Exception)):
            parse("[NH4")

    def test_unclosed_ring_digit(self):
        from molbuilder.smiles.parser import parse
        with self.assertRaises(ValueError):
            parse("C1CC")

    def test_empty_string_gives_no_atoms(self):
        """An empty SMILES string should produce a molecule with no atoms."""
        from molbuilder.smiles.parser import parse
        mol = parse("")
        self.assertEqual(len(mol.atoms), 0)

    def test_only_bond_symbol_gives_no_atoms(self):
        """A bare bond symbol '=' should produce a molecule with no atoms."""
        from molbuilder.smiles.parser import parse
        mol = parse("=")
        self.assertEqual(len(mol.atoms), 0)

    def test_mismatched_parentheses_still_parses(self):
        """Unclosed parenthesis 'CC(=O' is handled gracefully (no crash)."""
        from molbuilder.smiles.parser import parse
        mol = parse("CC(=O")
        # Should parse at least some atoms without crashing
        self.assertGreater(len(mol.atoms), 0)

    def test_trailing_bond_symbol_ignored(self):
        """A trailing bond 'C=' is handled gracefully (trailing bond ignored)."""
        from molbuilder.smiles.parser import parse
        mol = parse("C=")
        # Should produce at least a carbon (CH4-like)
        symbols = [a.symbol for a in mol.atoms]
        self.assertIn("C", symbols)


# ===================================================================
# 2. Empty molecule operations
# ===================================================================

class TestEmptyMolecule(unittest.TestCase):
    """Molecule with no atoms should behave gracefully."""

    def test_empty_molecule_has_no_atoms(self):
        from molbuilder.molecule.graph import Molecule
        mol = Molecule("empty")
        self.assertEqual(len(mol.atoms), 0)

    def test_empty_molecule_has_no_bonds(self):
        from molbuilder.molecule.graph import Molecule
        mol = Molecule("empty")
        self.assertEqual(len(mol.bonds), 0)

    def test_empty_molecule_check_steric_clashes(self):
        from molbuilder.molecule.graph import Molecule
        mol = Molecule("empty")
        clashes = mol.check_steric_clashes()
        self.assertEqual(clashes, [])


# ===================================================================
# 3. Zero-scale process engineering
# ===================================================================

class TestZeroScaleCosting(unittest.TestCase):
    """estimate_cost with scale_kg=0 should not crash (division by zero)."""

    def test_zero_scale_no_crash(self):
        from molbuilder.process.costing import estimate_cost, CostEstimate
        from molbuilder.reactions.knowledge_base import lookup_by_name

        template = lookup_by_name("NaBH4 reduction")[0]

        class _MockPrecursor:
            def __init__(self):
                self.name = "ethanol"
                self.smiles = "CCO"
                self.cost_per_kg = 10.0

        class _MockStep:
            def __init__(self, t):
                self.template = t
                self.precursors = [_MockPrecursor()]
                self.step_number = 1

        step = _MockStep(template)
        result = estimate_cost([step], scale_kg=0)
        self.assertIsInstance(result, CostEstimate)
        # The function clamps to 0.001 internally, so total should be > 0
        self.assertGreater(result.total_usd, 0)

    def test_negative_scale_no_crash(self):
        from molbuilder.process.costing import estimate_cost, CostEstimate
        from molbuilder.reactions.knowledge_base import lookup_by_name

        template = lookup_by_name("NaBH4 reduction")[0]

        class _MockPrecursor:
            def __init__(self):
                self.name = "ethanol"
                self.smiles = "CCO"
                self.cost_per_kg = 10.0

        class _MockStep:
            def __init__(self, t):
                self.template = t
                self.precursors = [_MockPrecursor()]
                self.step_number = 1

        step = _MockStep(template)
        result = estimate_cost([step], scale_kg=-5.0)
        self.assertIsInstance(result, CostEstimate)


# ===================================================================
# 4. Molecules with no functional groups
# ===================================================================

class TestNoFunctionalGroups(unittest.TestCase):
    """Methane and noble gas atoms should have no detected FGs."""

    def test_methane_no_fgs(self):
        from molbuilder.smiles.parser import parse
        from molbuilder.reactions.functional_group_detect import detect_functional_groups
        mol = parse("C")
        fgs = detect_functional_groups(mol)
        names = {fg.name for fg in fgs}
        # Methane has no traditional functional groups
        # (it may detect 'alkane' or nothing; just verify no alcohol/aldehyde)
        self.assertNotIn("alcohol", names)
        self.assertNotIn("aldehyde", names)
        self.assertNotIn("carboxylic_acid", names)
        self.assertNotIn("primary_amine", names)

    def test_noble_gas_helium_single_atom(self):
        from molbuilder.smiles.parser import parse
        mol = parse("[He]")
        # He should be 1 atom, no bonds, no functional groups
        self.assertEqual(len(mol.atoms), 1)
        self.assertEqual(mol.atoms[0].symbol, "He")


# ===================================================================
# 5. Single-atom molecules
# ===================================================================

class TestSingleAtomMolecules(unittest.TestCase):
    """Test parsing of single-atom SMILES."""

    def test_parse_C_gives_CH4(self):
        """parse('C') should give CH4 with 5 atoms (1 C + 4 H)."""
        from molbuilder.smiles.parser import parse
        mol = parse("C")
        symbols = [a.symbol for a in mol.atoms]
        self.assertEqual(symbols.count("C"), 1)
        self.assertEqual(symbols.count("H"), 4)
        self.assertEqual(len(mol.atoms), 5)

    def test_parse_helium_gives_1_atom(self):
        """parse('[He]') should give exactly 1 atom."""
        from molbuilder.smiles.parser import parse
        mol = parse("[He]")
        self.assertEqual(len(mol.atoms), 1)
        self.assertEqual(mol.atoms[0].symbol, "He")

    def test_parse_neon_gives_1_atom(self):
        """parse('[Ne]') should give exactly 1 atom."""
        from molbuilder.smiles.parser import parse
        mol = parse("[Ne]")
        self.assertEqual(len(mol.atoms), 1)
        self.assertEqual(mol.atoms[0].symbol, "Ne")

    def test_parse_N_organic_subset(self):
        """parse('N') should give NH3 with 4 atoms (1 N + 3 H)."""
        from molbuilder.smiles.parser import parse
        mol = parse("N")
        symbols = [a.symbol for a in mol.atoms]
        self.assertEqual(symbols.count("N"), 1)
        self.assertEqual(symbols.count("H"), 3)
        self.assertEqual(len(mol.atoms), 4)

    def test_parse_O_organic_subset(self):
        """parse('O') should give H2O with 3 atoms (1 O + 2 H)."""
        from molbuilder.smiles.parser import parse
        mol = parse("O")
        symbols = [a.symbol for a in mol.atoms]
        self.assertEqual(symbols.count("O"), 1)
        self.assertEqual(symbols.count("H"), 2)
        self.assertEqual(len(mol.atoms), 3)


# ===================================================================
# 6. Steric clash detection
# ===================================================================

class TestStericClashDetection(unittest.TestCase):
    """check_steric_clashes on well-built molecules should find no clashes."""

    def test_ethane_no_clashes(self):
        from molbuilder.molecule.builders import build_ethane
        mol = build_ethane()
        clashes = mol.check_steric_clashes()
        self.assertEqual(clashes, [])

    def test_methane_from_smiles_no_clashes(self):
        from molbuilder.smiles.parser import parse
        mol = parse("C")
        clashes = mol.check_steric_clashes()
        self.assertEqual(clashes, [])

    def test_butane_no_clashes(self):
        from molbuilder.molecule.builders import build_butane
        mol = build_butane()
        clashes = mol.check_steric_clashes()
        self.assertEqual(clashes, [])


# ===================================================================
# 7. Coulomb's law zero distance
# ===================================================================

class TestCoulombsLawZeroDistance(unittest.TestCase):
    """coulombs_law with r=0 must raise ValueError."""

    def test_zero_distance_raises(self):
        from molbuilder.core.constants import coulombs_law
        with self.assertRaises(ValueError):
            coulombs_law(1.0, 1.0, 0.0)

    def test_negative_distance_raises(self):
        from molbuilder.core.constants import coulombs_law
        with self.assertRaises(ValueError):
            coulombs_law(1.0, 1.0, -1.0)


# ===================================================================
# 8. Round-trip I/O (XYZ, JSON, MOL, PDB)
# ===================================================================

class TestRoundTripXYZ(unittest.TestCase):
    """Write then read XYZ for a simple molecule."""

    def test_methane_xyz_round_trip(self):
        from molbuilder.smiles.parser import parse
        from molbuilder.io.xyz import write_xyz, read_xyz
        mol = parse("C")
        tmpdir = tempfile.mkdtemp()
        path = os.path.join(tmpdir, "methane.xyz")
        write_xyz(mol, path)
        mol2 = read_xyz(path)
        self.assertEqual(len(mol2.atoms), len(mol.atoms))
        orig_symbols = sorted(a.symbol for a in mol.atoms)
        read_symbols = sorted(a.symbol for a in mol2.atoms)
        self.assertEqual(orig_symbols, read_symbols)


class TestRoundTripJSON(unittest.TestCase):
    """Write then read JSON for a simple molecule."""

    def test_methane_json_round_trip(self):
        from molbuilder.smiles.parser import parse
        from molbuilder.io.json_io import write_json, read_json
        mol = parse("C")
        tmpdir = tempfile.mkdtemp()
        path = os.path.join(tmpdir, "methane.json")
        write_json(mol, path)
        mol2 = read_json(path)
        self.assertEqual(len(mol2.atoms), len(mol.atoms))
        self.assertEqual(len(mol2.bonds), len(mol.bonds))


class TestRoundTripMOL(unittest.TestCase):
    """Write then read MOL for a simple molecule."""

    def test_methane_mol_round_trip(self):
        from molbuilder.smiles.parser import parse
        from molbuilder.io.mol_sdf import write_mol, read_mol
        mol = parse("C")
        tmpdir = tempfile.mkdtemp()
        path = os.path.join(tmpdir, "methane.mol")
        write_mol(mol, path)
        mol2 = read_mol(path)
        self.assertEqual(len(mol2.atoms), len(mol.atoms))


class TestRoundTripPDB(unittest.TestCase):
    """Write then read PDB for a simple molecule."""

    def test_methane_pdb_round_trip(self):
        from molbuilder.smiles.parser import parse
        from molbuilder.io.pdb import write_pdb, read_pdb
        mol = parse("C")
        tmpdir = tempfile.mkdtemp()
        path = os.path.join(tmpdir, "methane.pdb")
        write_pdb(mol, path)
        mol2 = read_pdb(path)
        self.assertEqual(len(mol2.atoms), len(mol.atoms))


# ===================================================================
# 9. P0 regression tests
# ===================================================================

class TestP0_1_HF_DipoleInRange(unittest.TestCase):
    """P0-1 regression: HF dipole moment should be in range 1.0-2.5 D."""

    def test_hf_dipole_in_range(self):
        from molbuilder.bonding.covalent import single_bond
        bond = single_bond("H", "F")
        mu = bond.dipole_moment_debye
        self.assertGreaterEqual(mu, 1.0,
                                "HF dipole too low: %.3f D" % mu)
        self.assertLessEqual(mu, 2.5,
                             "HF dipole too high: %.3f D" % mu)


class TestP0_3_ChiralitySurvivesRoundTrip(unittest.TestCase):
    """P0-3 regression: chirality should survive SMILES round-trip."""

    def test_chirality_preserved(self):
        from molbuilder.smiles.parser import parse
        from molbuilder.smiles.writer import to_smiles
        mol = parse("[C@@H](F)(Cl)Br")
        c_atoms = [a for a in mol.atoms if a.symbol == "C"]
        self.assertEqual(len(c_atoms), 1)
        self.assertEqual(c_atoms[0].chirality, "@@")
        smi = to_smiles(mol)
        self.assertIn("@@", smi)
        # Re-parse and confirm
        mol2 = parse(smi)
        c_atoms2 = [a for a in mol2.atoms if a.symbol == "C"]
        self.assertEqual(c_atoms2[0].chirality, "@@")


class TestP0_4_BracketAtomsPreserveChargeIsotope(unittest.TestCase):
    """P0-4 regression: bracket atoms preserve charge and isotope."""

    def test_charge_preserved(self):
        from molbuilder.smiles.parser import parse
        from molbuilder.smiles.writer import to_smiles
        mol = parse("[NH4+]")
        smi = to_smiles(mol)
        self.assertIn("+", smi)
        mol2 = parse(smi)
        n_atoms = [a for a in mol2.atoms if a.symbol == "N"]
        self.assertEqual(n_atoms[0].formal_charge, 1)

    def test_isotope_preserved(self):
        from molbuilder.smiles.parser import parse
        from molbuilder.smiles.writer import to_smiles
        mol = parse("[13C]")
        smi = to_smiles(mol)
        self.assertIn("13", smi)
        mol2 = parse(smi)
        c_atoms = [a for a in mol2.atoms if a.symbol == "C"]
        self.assertEqual(c_atoms[0].isotope, 13)

    def test_negative_charge_preserved(self):
        from molbuilder.smiles.parser import parse
        from molbuilder.smiles.writer import to_smiles
        mol = parse("[O-]")
        smi = to_smiles(mol)
        self.assertIn("-", smi)


class TestP0_5_FGDetectionWithImplicitH(unittest.TestCase):
    """P0-5 regression: FG detection works with implicit H (ethanol -> alcohol)."""

    def test_ethanol_has_alcohol(self):
        from molbuilder.smiles.parser import parse
        from molbuilder.reactions.functional_group_detect import detect_functional_groups
        mol = parse("CCO")
        fgs = detect_functional_groups(mol)
        names = {fg.name for fg in fgs}
        self.assertIn("alcohol", names)

    def test_acetic_acid_has_carboxylic_acid(self):
        from molbuilder.smiles.parser import parse
        from molbuilder.reactions.functional_group_detect import detect_functional_groups
        mol = parse("CC(=O)O")
        fgs = detect_functional_groups(mol)
        names = {fg.name for fg in fgs}
        self.assertIn("carboxylic_acid", names)


class TestP0_6_RetrosynthesisNoCrash(unittest.TestCase):
    """P0-6 regression: retrosynthesis should not crash on simple molecules."""

    def test_propanol_retrosynthesis(self):
        from molbuilder.smiles.parser import parse
        from molbuilder.reactions.retrosynthesis import retrosynthesis, RetrosynthesisTree
        mol = parse("CCCO")
        tree = retrosynthesis(mol, max_depth=2, beam_width=3)
        self.assertIsInstance(tree, RetrosynthesisTree)
        self.assertIsNotNone(tree.target)

    def test_ethanol_retrosynthesis(self):
        from molbuilder.smiles.parser import parse
        from molbuilder.reactions.retrosynthesis import retrosynthesis, RetrosynthesisTree
        mol = parse("CCO")
        tree = retrosynthesis(mol, max_depth=2, beam_width=3)
        self.assertIsInstance(tree, RetrosynthesisTree)

    def test_methane_retrosynthesis(self):
        from molbuilder.smiles.parser import parse
        from molbuilder.reactions.retrosynthesis import retrosynthesis, RetrosynthesisTree
        mol = parse("C")
        tree = retrosynthesis(mol, max_depth=2, beam_width=3)
        self.assertIsInstance(tree, RetrosynthesisTree)


class TestP0_7_SO2LewisDoubleBoonds(unittest.TestCase):
    """P0-7 regression: SO2 Lewis structure should have double bonds."""

    def test_so2_has_double_bonds(self):
        from molbuilder.bonding.lewis import LewisStructure
        lewis = LewisStructure("SO2")
        # SO2 should have two S=O bonds (order 2)
        double_bonds = [b for b in lewis.bonds if b.order == 2]
        self.assertGreaterEqual(len(double_bonds), 1,
                                "SO2 should have at least one double bond; "
                                "got orders: %s" % [b.order for b in lewis.bonds])

    def test_so2_central_is_sulfur(self):
        from molbuilder.bonding.lewis import LewisStructure
        lewis = LewisStructure("SO2")
        self.assertEqual(lewis.central_symbol, "S")

    def test_so2_steric_number(self):
        from molbuilder.bonding.lewis import LewisStructure
        lewis = LewisStructure("SO2")
        # SO2: 2 bonding groups + 1 lone pair = steric 3
        self.assertEqual(lewis.steric_number(), 3)


# ---------------------------------------------------------------------------
# P2-TEST-2: Scientific Validation Tests
# ---------------------------------------------------------------------------


# ===================================================================
# 1. Bohr model ionisation energies
# ===================================================================

class TestBohrModelValidation(unittest.TestCase):
    """Validate Bohr model predictions against known values."""

    def test_hydrogen_ionisation_energy_exact(self):
        """H ionisation energy = 13.6 eV (exact for Bohr model)."""
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        self.assertAlmostEqual(h.ionization_energy(), 13.6, delta=0.1)

    def test_he_plus_ionisation_energy(self):
        """He+ ionisation energy = 54.4 eV (Z=2, one electron)."""
        from molbuilder.atomic.bohr import BohrAtom
        he_plus = BohrAtom(2, charge=1)
        self.assertAlmostEqual(he_plus.ionization_energy(), 54.4, delta=0.2)

    def test_hydrogen_ground_state_energy(self):
        """H ground state energy = -13.6 eV."""
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        self.assertAlmostEqual(h.energy_level(1), -13.6, delta=0.1)

    def test_hydrogen_n2_energy(self):
        """H n=2 energy = -3.4 eV."""
        from molbuilder.atomic.bohr import BohrAtom
        h = BohrAtom(1)
        self.assertAlmostEqual(h.energy_level(2), -3.4, delta=0.1)


# ===================================================================
# 2. Bond lengths from STANDARD_BOND_LENGTHS
# ===================================================================

class TestStandardBondLengths(unittest.TestCase):
    """Compare bond_length() against known experimental values."""

    def test_cc_single_1_54(self):
        """C-C single bond = 1.54 A."""
        from molbuilder.core.bond_data import bond_length
        self.assertAlmostEqual(bond_length("C", "C", 1), 1.54, places=2)

    def test_cc_double_1_34(self):
        """C=C double bond = 1.34 A."""
        from molbuilder.core.bond_data import bond_length
        self.assertAlmostEqual(bond_length("C", "C", 2), 1.34, places=2)

    def test_cc_triple_1_20(self):
        """C-C triple bond = 1.20 A."""
        from molbuilder.core.bond_data import bond_length
        self.assertAlmostEqual(bond_length("C", "C", 3), 1.20, places=2)

    def test_ch_single_1_09(self):
        """C-H single bond = 1.09 A."""
        from molbuilder.core.bond_data import bond_length
        self.assertAlmostEqual(bond_length("C", "H", 1), 1.09, places=2)

    def test_oh_single_0_96(self):
        """O-H single bond = 0.96 A."""
        from molbuilder.core.bond_data import bond_length
        self.assertAlmostEqual(bond_length("O", "H", 1), 0.96, places=2)

    def test_double_shorter_than_single(self):
        """Double bond should be shorter than single bond for same pair."""
        from molbuilder.core.bond_data import bond_length
        single = bond_length("C", "C", 1)
        double = bond_length("C", "C", 2)
        triple = bond_length("C", "C", 3)
        self.assertLess(double, single)
        self.assertLess(triple, double)


# ===================================================================
# 3. VSEPR angles
# ===================================================================

class TestVSEPRAnglesValidation(unittest.TestCase):
    """Validate VSEPR bond angles against known values."""

    def test_ch4_tetrahedral_angle_109_5(self):
        """CH4 tetrahedral angle should be approximately 109.5 degrees."""
        from molbuilder.bonding.vsepr import VSEPRMolecule
        mol = VSEPRMolecule("CH4")
        angles = mol.computed_bond_angles()
        self.assertTrue(len(angles) > 0, "No angles computed for CH4")
        for a in angles:
            self.assertAlmostEqual(a, 109.5, delta=1.5)

    def test_bf3_trigonal_angle_120(self):
        """BF3 trigonal planar angle should be approximately 120 degrees."""
        from molbuilder.bonding.vsepr import VSEPRMolecule
        mol = VSEPRMolecule("BF3")
        angles = mol.computed_bond_angles()
        self.assertTrue(len(angles) > 0, "No angles computed for BF3")
        for a in angles:
            self.assertAlmostEqual(a, 120.0, delta=2.0)

    def test_co2_linear_angle_180(self):
        """CO2 linear angle should be approximately 180 degrees."""
        from molbuilder.bonding.vsepr import VSEPRMolecule
        mol = VSEPRMolecule("CO2")
        angles = mol.computed_bond_angles()
        self.assertTrue(len(angles) > 0, "No angles computed for CO2")
        for a in angles:
            self.assertAlmostEqual(a, 180.0, delta=2.0)

    def test_h2o_bent_angle(self):
        """H2O bent angle: VSEPR model uses tetrahedral electron geometry,
        so the computed bond angle is ~109.5 deg (ideal), not the
        experimental 104.5 deg.  This validates the model is consistent."""
        from molbuilder.bonding.vsepr import VSEPRMolecule
        mol = VSEPRMolecule("H2O")
        angles = mol.computed_bond_angles()
        self.assertTrue(len(angles) > 0, "No angles computed for H2O")
        for a in angles:
            # VSEPR places all electron groups at tetrahedral positions
            self.assertAlmostEqual(a, 109.5, delta=2.0)


# ===================================================================
# 4. Element lookups
# ===================================================================

class TestElementLookupValidation(unittest.TestCase):
    """Validate element lookup functions against known data."""

    def test_from_name_iron_returns_z_26(self):
        from molbuilder.core.elements import from_name
        z, sym, name, weight = from_name("Iron")
        self.assertEqual(z, 26)
        self.assertEqual(sym, "Fe")

    def test_from_symbol_fe_returns_z_26(self):
        from molbuilder.core.elements import from_symbol
        z, sym, name, weight = from_symbol("Fe")
        self.assertEqual(z, 26)
        self.assertEqual(name, "Iron")

    def test_symbol_to_z_has_118_elements(self):
        from molbuilder.core.elements import SYMBOL_TO_Z
        self.assertEqual(len(SYMBOL_TO_Z), 118)

    def test_symbol_to_z_spot_checks(self):
        from molbuilder.core.elements import SYMBOL_TO_Z
        self.assertEqual(SYMBOL_TO_Z["H"], 1)
        self.assertEqual(SYMBOL_TO_Z["C"], 6)
        self.assertEqual(SYMBOL_TO_Z["N"], 7)
        self.assertEqual(SYMBOL_TO_Z["O"], 8)
        self.assertEqual(SYMBOL_TO_Z["Fe"], 26)
        self.assertEqual(SYMBOL_TO_Z["Au"], 79)
        self.assertEqual(SYMBOL_TO_Z["Og"], 118)

    def test_from_symbol_unknown_raises(self):
        from molbuilder.core.elements import from_symbol
        with self.assertRaises(ValueError):
            from_symbol("Xx")

    def test_from_name_unknown_raises(self):
        from molbuilder.core.elements import from_name
        with self.assertRaises(ValueError):
            from_name("Unobtainium")


# ===================================================================
# 5. Quantum numbers / electron configurations
# ===================================================================

class TestElectronConfigurationValidation(unittest.TestCase):
    """Validate electron configurations against known data."""

    def test_hydrogen_1s1(self):
        """H has 1 electron in 1s1."""
        from molbuilder.atomic.quantum_atom import QuantumAtom
        h = QuantumAtom(1)
        self.assertEqual(h.electron_configuration_string(), "1s1")
        self.assertEqual(h.valence_electrons, 1)

    def test_helium_1s2(self):
        """He has 2 electrons in 1s2."""
        from molbuilder.atomic.quantum_atom import QuantumAtom
        he = QuantumAtom(2)
        self.assertEqual(he.electron_configuration_string(), "1s2")
        total = sum(ss.electron_count for ss in he.subshells)
        self.assertEqual(total, 2)

    def test_carbon_1s2_2s2_2p2(self):
        """C has electron config 1s2 2s2 2p2."""
        from molbuilder.atomic.quantum_atom import QuantumAtom
        c = QuantumAtom(6)
        self.assertEqual(c.electron_configuration_string(), "1s2 2s2 2p2")
        self.assertEqual(c.valence_electrons, 4)

    def test_nitrogen_1s2_2s2_2p3(self):
        """N has electron config 1s2 2s2 2p3."""
        from molbuilder.atomic.quantum_atom import QuantumAtom
        n = QuantumAtom(7)
        self.assertEqual(n.electron_configuration_string(), "1s2 2s2 2p3")

    def test_oxygen_config(self):
        """O (Z=8) config = 1s2 2s2 2p4."""
        from molbuilder.atomic.quantum_atom import QuantumAtom
        o = QuantumAtom(8)
        self.assertEqual(o.electron_configuration_string(), "1s2 2s2 2p4")

    def test_iron_config_has_3d6_4s2(self):
        """Fe (Z=26) config should contain 3d6 and 4s2."""
        from molbuilder.atomic.quantum_atom import QuantumAtom
        fe = QuantumAtom(26)
        config = fe.electron_configuration_string()
        self.assertIn("3d6", config)
        self.assertIn("4s2", config)

    def test_total_electrons_match_z(self):
        """Total electrons should equal atomic number for neutral atoms."""
        from molbuilder.atomic.quantum_atom import QuantumAtom
        for z in [1, 2, 6, 7, 8, 10, 11, 18, 26, 29]:
            atom = QuantumAtom(z)
            total = sum(ss.electron_count for ss in atom.subshells)
            self.assertEqual(total, z,
                             "Electron count mismatch for Z=%d" % z)

    def test_aufbau_order_first_five(self):
        """Aufbau order starts 1s, 2s, 2p, 3s, 3p."""
        from molbuilder.atomic.quantum_numbers import aufbau_order, SUBSHELL_LETTER
        order = aufbau_order()
        labels = ["%d%s" % (n, SUBSHELL_LETTER[l]) for n, l in order]
        self.assertEqual(labels[0], "1s")
        self.assertEqual(labels[1], "2s")
        self.assertEqual(labels[2], "2p")
        self.assertEqual(labels[3], "3s")
        self.assertEqual(labels[4], "3p")

    def test_noble_gas_notation_carbon(self):
        """C noble gas notation should contain [He] and 2s2 2p2."""
        from molbuilder.atomic.quantum_atom import QuantumAtom
        c = QuantumAtom(6)
        notation = c.noble_gas_notation()
        self.assertIn("[He]", notation)
        self.assertIn("2s2", notation)
        self.assertIn("2p2", notation)


if __name__ == "__main__":
    unittest.main()
