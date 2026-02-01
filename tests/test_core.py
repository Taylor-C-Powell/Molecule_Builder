"""
Unit tests for molbuilder.core modules:
    constants, elements, element_properties, geometry, bond_data
"""

import math
import unittest

import numpy as np


# ---------------------------------------------------------------------------
# constants
# ---------------------------------------------------------------------------
class TestConstants(unittest.TestCase):
    """Verify fundamental physical constants and derived values."""

    def test_planck_constant(self):
        from molbuilder.core.constants import PLANCK_CONSTANT
        self.assertAlmostEqual(PLANCK_CONSTANT, 6.62607015e-34, places=42)

    def test_speed_of_light(self):
        from molbuilder.core.constants import SPEED_OF_LIGHT
        self.assertAlmostEqual(SPEED_OF_LIGHT, 2.99792458e8, places=0)

    def test_bohr_radius_metres(self):
        from molbuilder.core.constants import BOHR_RADIUS_M
        self.assertAlmostEqual(BOHR_RADIUS_M, 5.29177210903e-11, delta=1e-20)

    def test_bohr_radius_pm(self):
        from molbuilder.core.constants import BOHR_RADIUS_PM
        self.assertAlmostEqual(BOHR_RADIUS_PM, 52.9177210903, places=5)

    def test_electron_charge(self):
        from molbuilder.core.constants import ELECTRON_CHARGE
        self.assertAlmostEqual(ELECTRON_CHARGE, 1.602176634e-19, delta=1e-28)

    def test_electron_mass(self):
        from molbuilder.core.constants import ELECTRON_MASS
        self.assertAlmostEqual(ELECTRON_MASS, 9.1093837015e-31, delta=1e-40)

    def test_hbar_derived_from_planck(self):
        from molbuilder.core.constants import PLANCK_CONSTANT, HBAR
        expected = PLANCK_CONSTANT / (2 * math.pi)
        self.assertAlmostEqual(HBAR, expected, places=46)

    def test_rydberg_energy_eV(self):
        from molbuilder.core.constants import RYDBERG_ENERGY_EV
        self.assertAlmostEqual(RYDBERG_ENERGY_EV, 13.605693122994, places=6)

    def test_coulomb_constant(self):
        from molbuilder.core.constants import COULOMB_CONSTANT
        self.assertAlmostEqual(COULOMB_CONSTANT, 8.9875517873681764e9, places=0)

    def test_vacuum_permittivity(self):
        from molbuilder.core.constants import VACUUM_PERMITTIVITY
        self.assertAlmostEqual(VACUUM_PERMITTIVITY, 8.8541878128e-12, delta=1e-21)

    def test_ev_to_joules(self):
        from molbuilder.core.constants import EV_TO_JOULES
        self.assertAlmostEqual(EV_TO_JOULES, 1.602176634e-19, delta=1e-28)

    def test_fine_structure_alpha(self):
        from molbuilder.core.constants import FINE_STRUCTURE_ALPHA
        # alpha ~ 1/137.036
        self.assertAlmostEqual(FINE_STRUCTURE_ALPHA, 1.0 / 137.036, places=4)

    def test_max_electrons_per_shell_count(self):
        from molbuilder.core.constants import MAX_ELECTRONS_PER_SHELL
        self.assertEqual(len(MAX_ELECTRONS_PER_SHELL), 7)

    def test_max_electrons_per_shell_values(self):
        from molbuilder.core.constants import MAX_ELECTRONS_PER_SHELL
        expected = [2, 8, 18, 32, 50, 72, 98]
        self.assertEqual(MAX_ELECTRONS_PER_SHELL, expected)

    def test_coulombs_law_positive(self):
        from molbuilder.core.constants import coulombs_law, COULOMB_CONSTANT
        # Two 1-C charges at 1 m
        force = coulombs_law(1.0, 1.0, 1.0)
        self.assertAlmostEqual(force, COULOMB_CONSTANT, places=0)

    def test_coulombs_law_inverse_square(self):
        from molbuilder.core.constants import coulombs_law
        f1 = coulombs_law(1.0, 1.0, 1.0)
        f2 = coulombs_law(1.0, 1.0, 2.0)
        self.assertAlmostEqual(f1 / f2, 4.0, places=10)

    def test_debye_conversion(self):
        from molbuilder.core.constants import DEBYE_PER_E_ANGSTROM
        self.assertAlmostEqual(DEBYE_PER_E_ANGSTROM, 4.8032, places=4)


# ---------------------------------------------------------------------------
# elements
# ---------------------------------------------------------------------------
class TestElements(unittest.TestCase):
    """Verify periodic table data and lookup functions."""

    def test_elements_has_118_entries(self):
        from molbuilder.core.elements import ELEMENTS
        self.assertEqual(len(ELEMENTS), 118)

    def test_elements_keys_1_to_118(self):
        from molbuilder.core.elements import ELEMENTS
        self.assertEqual(set(ELEMENTS.keys()), set(range(1, 119)))

    def test_symbol_to_z_count(self):
        from molbuilder.core.elements import SYMBOL_TO_Z
        self.assertEqual(len(SYMBOL_TO_Z), 118)

    def test_symbol_to_z_hydrogen(self):
        from molbuilder.core.elements import SYMBOL_TO_Z
        self.assertEqual(SYMBOL_TO_Z["H"], 1)

    def test_symbol_to_z_carbon(self):
        from molbuilder.core.elements import SYMBOL_TO_Z
        self.assertEqual(SYMBOL_TO_Z["C"], 6)

    def test_symbol_to_z_oganesson(self):
        from molbuilder.core.elements import SYMBOL_TO_Z
        self.assertEqual(SYMBOL_TO_Z["Og"], 118)

    def test_from_symbol_returns_tuple(self):
        from molbuilder.core.elements import from_symbol
        result = from_symbol("C")
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 4)

    def test_from_symbol_carbon(self):
        from molbuilder.core.elements import from_symbol
        z, sym, name, weight = from_symbol("C")
        self.assertEqual(z, 6)
        self.assertEqual(sym, "C")
        self.assertEqual(name, "Carbon")
        self.assertAlmostEqual(weight, 12.011, places=3)

    def test_from_symbol_case_insensitive(self):
        from molbuilder.core.elements import from_symbol
        z1, _, _, _ = from_symbol("fe")
        z2, _, _, _ = from_symbol("Fe")
        z3, _, _, _ = from_symbol("FE")
        self.assertEqual(z1, 26)
        self.assertEqual(z2, 26)
        self.assertEqual(z3, 26)

    def test_from_symbol_unknown_raises(self):
        from molbuilder.core.elements import from_symbol
        with self.assertRaises(ValueError):
            from_symbol("Xx")

    def test_from_name_oxygen(self):
        from molbuilder.core.elements import from_name
        z, sym, name, weight = from_name("Oxygen")
        self.assertEqual(z, 8)
        self.assertEqual(sym, "O")
        self.assertEqual(name, "Oxygen")
        self.assertAlmostEqual(weight, 15.999, places=3)

    def test_from_name_case_insensitive(self):
        from molbuilder.core.elements import from_name
        z, _, _, _ = from_name("iron")
        self.assertEqual(z, 26)

    def test_from_name_unknown_raises(self):
        from molbuilder.core.elements import from_name
        with self.assertRaises(ValueError):
            from_name("Unobtainium")

    def test_atomic_weight_hydrogen(self):
        from molbuilder.core.elements import atomic_weight
        self.assertAlmostEqual(atomic_weight("H"), 1.008, places=3)

    def test_atomic_weight_gold(self):
        from molbuilder.core.elements import atomic_weight
        self.assertAlmostEqual(atomic_weight("Au"), 196.967, places=3)

    def test_atomic_weight_unknown_raises(self):
        from molbuilder.core.elements import atomic_weight
        with self.assertRaises(ValueError):
            atomic_weight("Zz")

    def test_noble_gases_dict(self):
        from molbuilder.core.elements import NOBLE_GASES
        self.assertIn(2, NOBLE_GASES)
        self.assertEqual(NOBLE_GASES[2], "He")
        self.assertIn(10, NOBLE_GASES)
        self.assertEqual(NOBLE_GASES[10], "Ne")
        self.assertIn(18, NOBLE_GASES)
        self.assertEqual(NOBLE_GASES[18], "Ar")


# ---------------------------------------------------------------------------
# element_properties
# ---------------------------------------------------------------------------
class TestElementProperties(unittest.TestCase):
    """Verify element chemical property data and helper functions."""

    # --- electronegativity ---
    def test_electronegativity_carbon(self):
        from molbuilder.core.element_properties import electronegativity
        self.assertAlmostEqual(electronegativity("C"), 2.55, places=2)

    def test_electronegativity_nitrogen(self):
        from molbuilder.core.element_properties import electronegativity
        self.assertAlmostEqual(electronegativity("N"), 3.04, places=2)

    def test_electronegativity_oxygen(self):
        from molbuilder.core.element_properties import electronegativity
        self.assertAlmostEqual(electronegativity("O"), 3.44, places=2)

    def test_electronegativity_fluorine(self):
        from molbuilder.core.element_properties import electronegativity
        self.assertAlmostEqual(electronegativity("F"), 3.98, places=2)

    def test_electronegativity_ordering(self):
        from molbuilder.core.element_properties import electronegativity
        # F > O > N > C
        self.assertGreater(electronegativity("F"), electronegativity("O"))
        self.assertGreater(electronegativity("O"), electronegativity("N"))
        self.assertGreater(electronegativity("N"), electronegativity("C"))

    def test_electronegativity_noble_gas_returns_zero(self):
        from molbuilder.core.element_properties import electronegativity
        self.assertAlmostEqual(electronegativity("He"), 0.0, places=5)

    def test_electronegativity_dict_raw_value(self):
        from molbuilder.core.element_properties import PAULING_ELECTRONEGATIVITY
        self.assertEqual(PAULING_ELECTRONEGATIVITY["He"], None)
        self.assertAlmostEqual(PAULING_ELECTRONEGATIVITY["H"], 2.20, places=2)

    # --- covalent radii ---
    def test_covalent_radius_hydrogen(self):
        from molbuilder.core.element_properties import covalent_radius_pm
        self.assertAlmostEqual(covalent_radius_pm("H"), 31.0, places=1)

    def test_covalent_radius_carbon(self):
        from molbuilder.core.element_properties import covalent_radius_pm
        self.assertAlmostEqual(covalent_radius_pm("C"), 76.0, places=1)

    def test_covalent_radius_oxygen(self):
        from molbuilder.core.element_properties import covalent_radius_pm
        self.assertAlmostEqual(covalent_radius_pm("O"), 66.0, places=1)

    def test_covalent_radius_fallback(self):
        from molbuilder.core.element_properties import covalent_radius_pm
        # Unknown symbol falls back to 150.0
        self.assertAlmostEqual(covalent_radius_pm("Zz"), 150.0, places=1)

    # --- CPK colours ---
    def test_cpk_color_carbon(self):
        from molbuilder.core.element_properties import cpk_color
        self.assertEqual(cpk_color("C"), "#909090")

    def test_cpk_color_oxygen(self):
        from molbuilder.core.element_properties import cpk_color
        self.assertEqual(cpk_color("O"), "#FF0D0D")

    def test_cpk_color_hydrogen(self):
        from molbuilder.core.element_properties import cpk_color
        self.assertEqual(cpk_color("H"), "#FFFFFF")

    def test_cpk_color_default_fallback(self):
        from molbuilder.core.element_properties import cpk_color
        self.assertEqual(cpk_color("Zz"), "#FF69B4")

    # --- target electrons ---
    def test_target_electrons_hydrogen_duet(self):
        from molbuilder.core.element_properties import target_electrons
        self.assertEqual(target_electrons("H"), 2)

    def test_target_electrons_helium_duet(self):
        from molbuilder.core.element_properties import target_electrons
        self.assertEqual(target_electrons("He"), 2)

    def test_target_electrons_boron(self):
        from molbuilder.core.element_properties import target_electrons
        self.assertEqual(target_electrons("B"), 6)

    def test_target_electrons_carbon_octet(self):
        from molbuilder.core.element_properties import target_electrons
        # Carbon is not in TARGET_ELECTRONS dict, so default 8
        self.assertEqual(target_electrons("C"), 8)

    def test_target_electrons_oxygen_octet(self):
        from molbuilder.core.element_properties import target_electrons
        self.assertEqual(target_electrons("O"), 8)

    # --- period ---
    def test_period_hydrogen(self):
        from molbuilder.core.element_properties import period
        self.assertEqual(period("H"), 1)

    def test_period_carbon(self):
        from molbuilder.core.element_properties import period
        self.assertEqual(period("C"), 2)

    def test_period_sodium(self):
        from molbuilder.core.element_properties import period
        self.assertEqual(period("Na"), 3)

    def test_period_iron(self):
        from molbuilder.core.element_properties import period
        self.assertEqual(period("Fe"), 4)

    def test_period_iodine(self):
        from molbuilder.core.element_properties import period
        self.assertEqual(period("I"), 5)

    def test_period_gold(self):
        from molbuilder.core.element_properties import period
        self.assertEqual(period("Au"), 6)

    def test_period_francium(self):
        from molbuilder.core.element_properties import period
        self.assertEqual(period("Fr"), 7)

    # --- can_expand_octet ---
    def test_can_expand_octet_carbon_false(self):
        from molbuilder.core.element_properties import can_expand_octet
        self.assertFalse(can_expand_octet("C"))

    def test_can_expand_octet_nitrogen_false(self):
        from molbuilder.core.element_properties import can_expand_octet
        self.assertFalse(can_expand_octet("N"))

    def test_can_expand_octet_sulfur_true(self):
        from molbuilder.core.element_properties import can_expand_octet
        self.assertTrue(can_expand_octet("S"))

    def test_can_expand_octet_phosphorus_true(self):
        from molbuilder.core.element_properties import can_expand_octet
        self.assertTrue(can_expand_octet("P"))

    # --- estimated bond lengths ---
    def test_estimated_bond_length_pm_CC_single(self):
        from molbuilder.core.element_properties import estimated_bond_length_pm
        result = estimated_bond_length_pm("C", "C", 1)
        # C covalent radius 76 pm, so single = 152 pm
        self.assertAlmostEqual(result, 152.0, places=1)

    def test_estimated_bond_length_pm_double_shorter(self):
        from molbuilder.core.element_properties import estimated_bond_length_pm
        single = estimated_bond_length_pm("C", "C", 1)
        double = estimated_bond_length_pm("C", "C", 2)
        triple = estimated_bond_length_pm("C", "C", 3)
        self.assertLess(double, single)
        self.assertLess(triple, double)

    def test_estimated_bond_length_angstrom(self):
        from molbuilder.core.element_properties import (
            estimated_bond_length_pm, estimated_bond_length_angstrom
        )
        pm = estimated_bond_length_pm("C", "H", 1)
        ang = estimated_bond_length_angstrom("C", "H", 1)
        self.assertAlmostEqual(ang, pm / 100.0, places=5)


# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------
class TestGeometry(unittest.TestCase):
    """Verify 3D geometry utilities."""

    # --- normalize ---
    def test_normalize_produces_unit_vector(self):
        from molbuilder.core.geometry import normalize
        v = np.array([3.0, 4.0, 0.0])
        result = normalize(v)
        self.assertAlmostEqual(np.linalg.norm(result), 1.0, places=12)

    def test_normalize_direction_preserved(self):
        from molbuilder.core.geometry import normalize
        v = np.array([0.0, 0.0, 5.0])
        result = normalize(v)
        np.testing.assert_allclose(result, [0.0, 0.0, 1.0], atol=1e-12)

    def test_normalize_zero_vector_returns_zero(self):
        from molbuilder.core.geometry import normalize
        v = np.array([0.0, 0.0, 0.0])
        result = normalize(v)
        np.testing.assert_allclose(result, [0.0, 0.0, 0.0], atol=1e-12)

    def test_normalize_arbitrary_vector(self):
        from molbuilder.core.geometry import normalize
        v = np.array([1.0, 2.0, 3.0])
        result = normalize(v)
        self.assertAlmostEqual(np.linalg.norm(result), 1.0, places=12)
        # Direction should be proportional to original
        ratio = v / result
        self.assertAlmostEqual(ratio[0], ratio[1], places=10)
        self.assertAlmostEqual(ratio[1], ratio[2], places=10)

    # --- rotation_matrix ---
    def test_rotation_matrix_preserves_length(self):
        from molbuilder.core.geometry import rotation_matrix
        axis = np.array([0.0, 0.0, 1.0])
        angle = math.radians(45.0)
        R = rotation_matrix(axis, angle)
        v = np.array([1.0, 0.0, 0.0])
        v_rotated = R @ v
        self.assertAlmostEqual(np.linalg.norm(v_rotated), 1.0, places=12)

    def test_rotation_matrix_90_deg_z_axis(self):
        from molbuilder.core.geometry import rotation_matrix
        axis = np.array([0.0, 0.0, 1.0])
        R = rotation_matrix(axis, math.pi / 2)
        v = np.array([1.0, 0.0, 0.0])
        v_rotated = R @ v
        np.testing.assert_allclose(v_rotated, [0.0, 1.0, 0.0], atol=1e-10)

    def test_rotation_matrix_360_deg_identity(self):
        from molbuilder.core.geometry import rotation_matrix
        axis = np.array([1.0, 1.0, 1.0])
        R = rotation_matrix(axis, 2 * math.pi)
        np.testing.assert_allclose(R, np.eye(3), atol=1e-10)

    def test_rotation_matrix_preserves_axis(self):
        from molbuilder.core.geometry import rotation_matrix
        axis = np.array([0.0, 0.0, 1.0])
        R = rotation_matrix(axis, math.radians(73.0))
        v_rotated = R @ axis
        np.testing.assert_allclose(v_rotated, axis, atol=1e-10)

    # --- place_atom_zmatrix ---
    def test_place_atom_zmatrix_correct_distance(self):
        from molbuilder.core.geometry import place_atom_zmatrix
        ref = np.array([0.0, 0.0, 0.0])
        angle_ref = np.array([-1.0, 0.0, 0.0])
        dihedral_ref = np.array([-1.0, -1.0, 0.0])
        bond_length = 1.54
        result = place_atom_zmatrix(ref, angle_ref, dihedral_ref,
                                    bond_length, 109.47, 120.0)
        dist = np.linalg.norm(result - ref)
        self.assertAlmostEqual(dist, bond_length, places=5)

    def test_place_atom_zmatrix_different_bond_length(self):
        from molbuilder.core.geometry import place_atom_zmatrix
        ref = np.array([1.0, 0.0, 0.0])
        angle_ref = np.array([0.0, 0.0, 0.0])
        dihedral_ref = np.array([0.0, 1.0, 0.0])
        bond_length = 2.0
        result = place_atom_zmatrix(ref, angle_ref, dihedral_ref,
                                    bond_length, 120.0, 0.0)
        dist = np.linalg.norm(result - ref)
        self.assertAlmostEqual(dist, bond_length, places=5)

    # --- coordinate conversions ---
    def test_cartesian_to_spherical_on_z_axis(self):
        from molbuilder.core.geometry import cartesian_to_spherical
        r, theta, phi = cartesian_to_spherical(0.0, 0.0, 5.0)
        self.assertAlmostEqual(float(r), 5.0, places=10)
        self.assertAlmostEqual(float(theta), 0.0, places=10)

    def test_spherical_to_cartesian_roundtrip(self):
        from molbuilder.core.geometry import (
            cartesian_to_spherical, spherical_to_cartesian
        )
        x0, y0, z0 = 1.0, 2.0, 3.0
        r, theta, phi = cartesian_to_spherical(x0, y0, z0)
        x1, y1, z1 = spherical_to_cartesian(r, theta, phi)
        self.assertAlmostEqual(float(x1), x0, places=10)
        self.assertAlmostEqual(float(y1), y0, places=10)
        self.assertAlmostEqual(float(z1), z0, places=10)

    # --- available_tetrahedral_dirs ---
    def test_tetrahedral_dirs_from_zero(self):
        from molbuilder.core.geometry import available_tetrahedral_dirs
        dirs = available_tetrahedral_dirs([], 4)
        self.assertEqual(len(dirs), 4)
        # Each should be a unit vector
        for d in dirs:
            self.assertAlmostEqual(np.linalg.norm(d), 1.0, places=10)

    def test_tetrahedral_dirs_from_one(self):
        from molbuilder.core.geometry import available_tetrahedral_dirs
        existing = [np.array([1.0, 0.0, 0.0])]
        dirs = available_tetrahedral_dirs(existing, 3)
        self.assertEqual(len(dirs), 3)
        for d in dirs:
            self.assertAlmostEqual(np.linalg.norm(d), 1.0, places=10)


# ---------------------------------------------------------------------------
# bond_data
# ---------------------------------------------------------------------------
class TestBondData(unittest.TestCase):
    """Verify bond length lookups, angle constants, BDE data, torsion barriers."""

    # --- bond_length lookups ---
    def test_bond_length_CC_single(self):
        from molbuilder.core.bond_data import bond_length
        self.assertAlmostEqual(bond_length("C", "C", 1), 1.54, places=2)

    def test_bond_length_CC_double(self):
        from molbuilder.core.bond_data import bond_length
        self.assertAlmostEqual(bond_length("C", "C", 2), 1.34, places=2)

    def test_bond_length_CC_triple(self):
        from molbuilder.core.bond_data import bond_length
        self.assertAlmostEqual(bond_length("C", "C", 3), 1.20, places=2)

    def test_bond_length_CH(self):
        from molbuilder.core.bond_data import bond_length
        self.assertAlmostEqual(bond_length("C", "H", 1), 1.09, places=2)

    def test_bond_length_reversed_symbols(self):
        from molbuilder.core.bond_data import bond_length
        # bond_length should sort symbols internally
        self.assertAlmostEqual(bond_length("H", "C", 1), 1.09, places=2)

    def test_bond_length_OH(self):
        from molbuilder.core.bond_data import bond_length
        self.assertAlmostEqual(bond_length("O", "H", 1), 0.96, places=2)

    def test_bond_length_fallback_to_covalent_radii(self):
        from molbuilder.core.bond_data import bond_length
        from molbuilder.core.element_properties import estimated_bond_length_angstrom
        # Use a pair not in STANDARD_BOND_LENGTHS
        result = bond_length("Li", "F", 1)
        expected = estimated_bond_length_angstrom("Li", "F", 1)
        self.assertAlmostEqual(result, expected, places=5)

    # --- standard angles ---
    def test_sp3_angle(self):
        from molbuilder.core.bond_data import SP3_ANGLE
        self.assertAlmostEqual(SP3_ANGLE, 109.47, places=2)

    def test_sp2_angle(self):
        from molbuilder.core.bond_data import SP2_ANGLE
        self.assertAlmostEqual(SP2_ANGLE, 120.0, places=1)

    def test_sp_angle(self):
        from molbuilder.core.bond_data import SP_ANGLE
        self.assertAlmostEqual(SP_ANGLE, 180.0, places=1)

    # --- BDE table ---
    def test_bde_HH(self):
        from molbuilder.core.bond_data import BDE_TABLE
        self.assertAlmostEqual(BDE_TABLE[("H", "H", 1)], 436, places=0)

    def test_bde_CH(self):
        from molbuilder.core.bond_data import BDE_TABLE
        self.assertAlmostEqual(BDE_TABLE[("C", "H", 1)], 413, places=0)

    def test_bde_CC_double(self):
        from molbuilder.core.bond_data import BDE_TABLE
        self.assertAlmostEqual(BDE_TABLE[("C", "C", 2)], 614, places=0)

    def test_bde_CC_triple(self):
        from molbuilder.core.bond_data import BDE_TABLE
        self.assertAlmostEqual(BDE_TABLE[("C", "C", 3)], 839, places=0)

    def test_bde_NN_triple(self):
        from molbuilder.core.bond_data import BDE_TABLE
        self.assertAlmostEqual(BDE_TABLE[("N", "N", 3)], 941, places=0)

    def test_bde_lookup_function(self):
        from molbuilder.core.bond_data import bde_lookup
        self.assertAlmostEqual(bde_lookup("C", "H", 1), 413, places=0)

    def test_bde_lookup_sorted_keys(self):
        from molbuilder.core.bond_data import bde_lookup
        # Should sort internally so H,C -> C,H
        self.assertAlmostEqual(bde_lookup("H", "C", 1), 413, places=0)

    def test_bde_lookup_none_for_missing(self):
        from molbuilder.core.bond_data import bde_lookup
        self.assertIsNone(bde_lookup("Li", "F", 1))

    # --- torsion barriers ---
    def test_torsion_barriers_H_sp3_sp3_H(self):
        from molbuilder.core.bond_data import TORSION_BARRIERS
        entry = TORSION_BARRIERS["H_sp3_sp3_H"]
        self.assertAlmostEqual(entry["V1"], 0.0, places=2)
        self.assertAlmostEqual(entry["V2"], 0.0, places=2)
        self.assertAlmostEqual(entry["V3"], 1.39, places=2)

    def test_torsion_barriers_C_sp3_sp3_C(self):
        from molbuilder.core.bond_data import TORSION_BARRIERS
        entry = TORSION_BARRIERS["C_sp3_sp3_C"]
        self.assertAlmostEqual(entry["V1"], 2.73, places=2)
        self.assertAlmostEqual(entry["V2"], -0.53, places=2)
        self.assertAlmostEqual(entry["V3"], 0.84, places=2)

    def test_torsion_barriers_default(self):
        from molbuilder.core.bond_data import TORSION_BARRIERS
        entry = TORSION_BARRIERS["default"]
        self.assertAlmostEqual(entry["V3"], 1.00, places=2)

    # --- polarity thresholds ---
    def test_nonpolar_threshold(self):
        from molbuilder.core.bond_data import NONPOLAR_THRESHOLD
        self.assertAlmostEqual(NONPOLAR_THRESHOLD, 0.4, places=2)

    def test_polar_covalent_max(self):
        from molbuilder.core.bond_data import POLAR_COVALENT_MAX
        self.assertAlmostEqual(POLAR_COVALENT_MAX, 1.7, places=2)


if __name__ == "__main__":
    unittest.main()
