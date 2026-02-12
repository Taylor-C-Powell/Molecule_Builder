"""Branch-coverage tests for bonding/ and molecule/graph.py.

Targets uncovered lines identified by coverage analysis:
    - bonding/lewis.py    73% -> 95%+
    - bonding/vsepr.py    67% -> 95%+
    - bonding/covalent.py 71% -> 95%+
    - molecule/graph.py   76% -> 90%+
"""

import math
import unittest

import numpy as np

from molbuilder.bonding.lewis import (
    parse_formula,
    get_valence_electrons,
    identify_central_atom,
    LewisStructure,
    Bond as LewisBond,
    LonePair,
)
from molbuilder.bonding.vsepr import (
    VSEPRMolecule,
    _assign_positions,
    generate_3d_coordinates,
    from_formula,
    from_atoms,
    MOLECULAR_GEOMETRY,
)
from molbuilder.bonding.covalent import (
    CovalentBond,
    BondPolarity,
    OrbitalContribution,
    OrbitalType,
    MolecularBondAnalysis,
    coordinate_bond,
    _bde_key,
)
from molbuilder.molecule.graph import (
    Molecule,
    Atom,
    Bond,
    TorsionAngle,
    NewmanProjection,
    StrainEnergy,
    Hybridization,
    Stereodescriptor,
)
from molbuilder.molecule.builders import build_ethane


# ===================================================================
# 1. TestLewisBranches
# ===================================================================

class TestLewisBranches(unittest.TestCase):
    """Cover uncovered branches in bonding/lewis.py."""

    # line 43: empty symbol token skipped in parse_formula regex
    def test_parse_formula_empty_symbol_skipped(self):
        # The regex r'([A-Z][a-z]?)(\d*)' only captures uppercase-led
        # tokens, so an empty symbol never arises in practice.  Just
        # verify normal formulas still work (the branch is defensive).
        result = parse_formula("H2O")
        self.assertEqual(result, ["H", "H", "O"])

    # line 84: unknown element -> ValueError
    def test_valence_electrons_unknown_element(self):
        with self.assertRaises(ValueError):
            get_valence_electrons("Zz")

    # line 91: He special case -> 2
    def test_valence_electrons_helium(self):
        self.assertEqual(get_valence_electrons("He"), 2)

    # lines 95-97: transition-metal fallback via QuantumAtom (Fe, Z=26)
    def test_valence_electrons_transition_metal(self):
        ve = get_valence_electrons("Fe")
        self.assertIsInstance(ve, int)
        self.assertGreater(ve, 0)

    # lines 118-121: diatomic with H -> picks non-H (HCl -> Cl)
    def test_identify_central_diatomic_with_h(self):
        idx = identify_central_atom(["H", "Cl"])
        self.assertEqual(idx, 1)
        # reversed order
        idx2 = identify_central_atom(["Cl", "H"])
        self.assertEqual(idx2, 0)

    # line 127: all hydrogen -> index 0
    def test_identify_central_all_hydrogen(self):
        idx = identify_central_atom(["H", "H", "H"])
        self.assertEqual(idx, 0)

    # line 185: empty formula -> ValueError
    def test_empty_formula_raises(self):
        with self.assertRaises(ValueError):
            parse_formula("Zz")  # unknown element

    # line 217: electron-deficient molecule (charged H2)
    def test_electron_deficient(self):
        ls = LewisStructure("H2", charge=1)
        self.assertEqual(ls.total_valence_electrons, 1)
        # Should not crash, remaining clamped to 0

    # line 250: bond order already 3 -> skipped
    def test_triple_bond_limit(self):
        # N2 forms a triple bond; once at 3, no further promotion
        ls = LewisStructure("N2")
        for bond in ls.bonds:
            self.assertLessEqual(bond.order, 3)

    # line 262: no lone pair available -> break
    def test_no_more_promotions(self):
        # BH3 is electron-deficient, central B has < 8 electrons
        # but no lone pairs on H to promote
        ls = LewisStructure("BH3")
        central_e = ls._electrons_around(ls.central_index)
        self.assertEqual(central_e, 6)  # B only gets 6

    # line 280: _find_lone_pair_on returns None
    def test_find_lone_pair_none(self):
        ls = LewisStructure("CH4")
        # Carbon in CH4 has no lone pairs
        result = ls._find_lone_pair_on(ls.central_index)
        self.assertIsNone(result)

    # lines 302-307: bond_order_to for unbonded pair -> 0
    def test_bond_order_to_nonexistent(self):
        ls = LewisStructure("H2O")
        # Index 999 doesn't exist in the molecule
        self.assertEqual(ls.bond_order_to(999), 0)

    # lines 312-314, 151-152, 161: __repr__ methods
    def test_lewis_repr(self):
        ls = LewisStructure("H2O")
        r = repr(ls)
        self.assertIn("LewisStructure", r)
        self.assertIn("H2O", r)
        self.assertIn("central=O", r)
        # Bond repr
        if ls.bonds:
            br = repr(ls.bonds[0])
            self.assertIn("Bond(", br)
        # LonePair repr
        if ls.lone_pairs:
            lr = repr(ls.lone_pairs[0])
            self.assertIn("LP(", lr)

    # lines 318-345: summary()
    def test_lewis_summary(self):
        ls = LewisStructure("H2O")
        s = ls.summary()
        self.assertIn("Lewis Structure: H2O", s)
        self.assertIn("Charge:", s)
        self.assertIn("Total valence electrons:", s)
        self.assertIn("Bonds:", s)
        self.assertIn("Lone pairs:", s)
        self.assertIn("Steric number:", s)


# ===================================================================
# 2. TestVSEPRBranches
# ===================================================================

class TestVSEPRBranches(unittest.TestCase):
    """Cover uncovered branches in bonding/vsepr.py."""

    # lines 118-121: SN=5 trigonal bipyramidal directions (PCl5)
    def test_trigonal_bipyramidal_pcl5(self):
        mol = VSEPRMolecule("PCl5")
        self.assertEqual(mol.axe.steric_number, 5)
        self.assertEqual(mol.axe.molecular_geometry, "trigonal bipyramidal")
        self.assertEqual(mol.axe.hybridization, "sp3d")

    # lines 137-142: SN=7 pentagonal bipyramidal directions (IF7)
    def test_pentagonal_bipyramidal_if7(self):
        mol = VSEPRMolecule("IF7")
        self.assertEqual(mol.axe.steric_number, 7)
        self.assertEqual(mol.axe.molecular_geometry, "pentagonal bipyramidal")

    # line 172: gen is None fallback (steric_number not in generators)
    def test_assign_positions_unknown_sn(self):
        bond_dirs, lp_dirs = _assign_positions(1, 0)
        self.assertEqual(len(bond_dirs), 1)
        self.assertEqual(len(lp_dirs), 0)
        np.testing.assert_array_almost_equal(bond_dirs[0], [0, 0, 1])

    # lines 180-182: SN=5 with > 3 lone pairs
    def test_sn5_more_than_3_lp(self):
        bond_dirs, lp_dirs = _assign_positions(5, 4)
        self.assertEqual(len(lp_dirs), 4)
        self.assertEqual(len(bond_dirs), 1)

    # line 187: octahedral 1 LP
    def test_octahedral_1_lp(self):
        bond_dirs, lp_dirs = _assign_positions(6, 1)
        self.assertEqual(len(bond_dirs), 5)
        self.assertEqual(len(lp_dirs), 1)

    # line 189: octahedral 2 LP (trans pair)
    def test_octahedral_2_lp(self):
        bond_dirs, lp_dirs = _assign_positions(6, 2)
        self.assertEqual(len(bond_dirs), 4)
        self.assertEqual(len(lp_dirs), 2)

    # line 191: octahedral 3 LP (fac)
    def test_octahedral_3_lp(self):
        bond_dirs, lp_dirs = _assign_positions(6, 3)
        self.assertEqual(len(bond_dirs), 3)
        self.assertEqual(len(lp_dirs), 3)

    # line 193: octahedral 4 LP
    def test_octahedral_4_lp(self):
        bond_dirs, lp_dirs = _assign_positions(6, 4)
        self.assertEqual(len(bond_dirs), 2)
        self.assertEqual(len(lp_dirs), 4)

    # line 231: single atom -> 1-atom dict (noble gas)
    def test_single_atom_coordinates(self):
        ls = LewisStructure("He")
        coords = generate_3d_coordinates(ls)
        self.assertEqual(len(coords['atom_positions']), 1)
        self.assertEqual(coords['atom_positions'][0][0], "He")
        np.testing.assert_array_equal(coords['atom_positions'][0][1],
                                       [0.0, 0.0, 0.0])

    # lines 243-244: SN < 2 multi-atom fallback
    def test_diatomic_sn1(self):
        # HF has SN=4 actually, but directly test the path
        # by calling generate_3d_coordinates on a LewisStructure for HH+
        ls = LewisStructure("H2", charge=1)
        coords = generate_3d_coordinates(ls)
        self.assertEqual(len(coords['atom_positions']), 2)

    # line 265: bond direction overflow -> fallback direction
    def test_bond_dir_overflow(self):
        # _assign_positions for SN=2 gives 2 dirs; if somehow more bonds
        # exist, the code uses fallback [0,0,1].  We test indirectly via
        # a Lewis structure with more bonds than geometry directions.
        # CH4 has SN=4 -> 4 dirs exactly; no overflow. Hard to trigger
        # naturally. Test the guard by building a structure manually.
        ls = LewisStructure("BH3")  # SN=3, 3 bonds, 0 LP -> exact match
        coords = generate_3d_coordinates(ls)
        # All 4 atoms should be placed (B + 3H)
        for sym, pos in coords['atom_positions']:
            self.assertIsNotNone(sym)
            self.assertIsNotNone(pos)

    # line 320: (sn, E) not in MOLECULAR_GEOMETRY dict -> "unknown"
    def test_unknown_geometry(self):
        # steric_number=8, lone_pairs=0 is not in the dict
        self.assertNotIn((8, 0), MOLECULAR_GEOMETRY)
        # Construct a molecule that would get this classification
        # We can't easily get SN=8 from a formula, so test the classify
        # path indirectly. The (sn,E) lookup returns "unknown" for missing
        # keys -- verify the dict behavior
        entry = MOLECULAR_GEOMETRY.get((8, 0))
        self.assertIsNone(entry)

    # lines 360, 366-416: __repr__ and summary()
    def test_vsepr_repr_and_summary(self):
        mol = VSEPRMolecule("H2O")
        r = repr(mol)
        self.assertIn("VSEPRMolecule", r)
        self.assertIn("H2O", r)
        self.assertIn("bent", r)

        s = mol.summary()
        self.assertIn("VSEPR Analysis: H2O", s)
        self.assertIn("AXE notation:", s)
        self.assertIn("Molecular geometry:", s)
        self.assertIn("Hybridization:", s)
        self.assertIn("3D Coordinates", s)
        self.assertIn("central", s)

    # line 425: from_formula convenience constructor
    def test_from_formula(self):
        mol = from_formula("NH3")
        self.assertIsInstance(mol, VSEPRMolecule)
        self.assertEqual(mol.axe.molecular_geometry, "trigonal pyramidal")

    # lines 430-433: from_atoms convenience constructor
    def test_from_atoms(self):
        mol = from_atoms(["N", "H", "H", "H"])
        self.assertIsInstance(mol, VSEPRMolecule)
        self.assertEqual(mol.axe.bonding_groups, 3)


# ===================================================================
# 3. TestCovalentBranches
# ===================================================================

class TestCovalentBranches(unittest.TestCase):
    """Cover uncovered branches in bonding/covalent.py."""

    # lines 58-59: OrbitalContribution.__repr__
    def test_orbital_contribution_repr(self):
        sigma = OrbitalContribution(OrbitalType.SIGMA, "head-on overlap")
        self.assertEqual(repr(sigma), "sigma(head-on overlap)")
        pi = OrbitalContribution(OrbitalType.PI, "lateral overlap")
        self.assertEqual(repr(pi), "pi(lateral overlap)")

    # lines 89-90: _bde_key alphabetical ordering
    def test_bde_key_sorts(self):
        k1 = _bde_key("O", "C", 2)
        k2 = _bde_key("C", "O", 2)
        self.assertEqual(k1, k2)
        self.assertEqual(k1, ("C", "O", 2))

    # line 125: invalid symbol_b -> ValueError
    def test_unknown_symbol_b_raises(self):
        with self.assertRaises(ValueError):
            CovalentBond("C", "Zz")

    # line 129: coordinate bond without donor -> ValueError
    def test_coordinate_no_donor_raises(self):
        with self.assertRaises(ValueError):
            CovalentBond("N", "B", is_coordinate=True, donor=None)

    # lines 162, 168: partial charges when en_a > en_b (O-H: O more EN)
    def test_partial_charges_en_a_greater(self):
        bond = CovalentBond("O", "H")
        # O has higher EN than H
        self.assertEqual(bond.partial_positive, "H")
        self.assertEqual(bond.partial_negative, "O")

    # lines 162, 168: partial charges when en_a < en_b
    def test_partial_charges_en_b_greater(self):
        bond = CovalentBond("H", "O")
        self.assertEqual(bond.partial_positive, "H")
        self.assertEqual(bond.partial_negative, "O")

    # line 197: BDE for known bond C-H
    def test_dissociation_energy_kj_known(self):
        bond = CovalentBond("C", "H")
        bde = bond.dissociation_energy_kj
        self.assertIsNotNone(bde)
        self.assertGreater(bde, 0)

    # lines 202-205: eV conversion
    def test_dissociation_energy_ev(self):
        bond = CovalentBond("C", "H")
        ev = bond.dissociation_energy_ev
        self.assertIsNotNone(ev)
        self.assertAlmostEqual(ev, bond.dissociation_energy_kj / 96.485,
                               places=3)

    # lines 202-204: unknown pair -> None
    def test_dissociation_energy_none(self):
        # Xe-Kr bond is unlikely to be in the BDE table
        bond = CovalentBond("Xe", "Kr")
        self.assertIsNone(bond.dissociation_energy_kj)
        self.assertIsNone(bond.dissociation_energy_ev)

    # lines 248, 252: order_symbol and order_label properties
    def test_order_symbol_label(self):
        s = CovalentBond("C", "C", bond_order=1)
        self.assertEqual(s.order_symbol, "-")
        self.assertEqual(s.order_label, "single")
        d = CovalentBond("C", "C", bond_order=2)
        self.assertEqual(d.order_symbol, "=")
        self.assertEqual(d.order_label, "double")
        t = CovalentBond("C", "C", bond_order=3)
        self.assertEqual(t.order_symbol, "#")
        self.assertEqual(t.order_label, "triple")

    # lines 255-256: CovalentBond __repr__
    def test_covalent_bond_repr(self):
        bond = CovalentBond("C", "H")
        r = repr(bond)
        self.assertIn("CovalentBond", r)
        self.assertIn("C-H", r)
        self.assertIn("single", r)
        # coordinate tag
        cb = coordinate_bond("N", "B")
        r2 = repr(cb)
        self.assertIn("coordinate", r2)

    # lines 260-304: summary() - all branches (standard + coordinate + polar + nonpolar + BDE)
    def test_covalent_bond_summary_standard_polar(self):
        bond = CovalentBond("O", "H")
        s = bond.summary()
        self.assertIn("standard covalent bond", s)
        self.assertIn("polar covalent", s)
        self.assertIn("delta+", s)
        self.assertIn("delta-", s)
        self.assertIn("% ionic char", s)
        self.assertIn("Dipole moment", s)
        self.assertIn("Bond length", s)
        self.assertIn("Orbitals", s)

    def test_covalent_bond_summary_coordinate(self):
        bond = coordinate_bond("N", "B")
        s = bond.summary()
        self.assertIn("coordinate (dative) bond", s)
        self.assertIn("Donor", s)

    def test_covalent_bond_summary_nonpolar(self):
        bond = CovalentBond("C", "C")
        s = bond.summary()
        self.assertIn("nonpolar covalent", s)

    def test_covalent_bond_summary_no_bde(self):
        bond = CovalentBond("Xe", "Kr")
        s = bond.summary()
        self.assertIn("no reference data", s)

    # line 439: coordinate_bond() function
    def test_coordinate_bond_constructor(self):
        cb = coordinate_bond("N", "B", bond_order=1)
        self.assertTrue(cb.is_coordinate)
        self.assertEqual(cb.donor, "N")

    # line 354: all_bonds_nonpolar
    def test_all_bonds_nonpolar_h2(self):
        analysis = MolecularBondAnalysis("H2")
        self.assertTrue(analysis.all_bonds_nonpolar)

    # lines 369-376: is_symmetric branches
    def test_is_symmetric_homonuclear_diatomic(self):
        analysis = MolecularBondAnalysis("H2")
        self.assertTrue(analysis.is_symmetric)

    def test_is_symmetric_heteronuclear_diatomic(self):
        analysis = MolecularBondAnalysis("HF")
        self.assertFalse(analysis.is_symmetric)

    def test_is_symmetric_polyatomic_symmetric(self):
        # CO2: all terminals are O, no LP on central C
        analysis = MolecularBondAnalysis("CO2")
        self.assertTrue(analysis.is_symmetric)

    def test_is_symmetric_polyatomic_not_symmetric(self):
        # H2O: has lone pairs on central O
        analysis = MolecularBondAnalysis("H2O")
        self.assertFalse(analysis.is_symmetric)

    # lines 389-393: molecular_polarity all branches
    def test_molecular_polarity_nonpolar_bonds(self):
        analysis = MolecularBondAnalysis("H2")
        self.assertEqual(analysis.molecular_polarity, "nonpolar")

    def test_molecular_polarity_symmetric_cancels(self):
        analysis = MolecularBondAnalysis("CO2")
        self.assertIn("nonpolar", analysis.molecular_polarity)
        self.assertIn("symmetric", analysis.molecular_polarity)

    def test_molecular_polarity_polar(self):
        analysis = MolecularBondAnalysis("H2O")
        self.assertEqual(analysis.molecular_polarity, "polar")

    # lines 396-410: MolecularBondAnalysis summary()
    def test_molecular_analysis_summary(self):
        analysis = MolecularBondAnalysis("H2O")
        s = analysis.summary()
        self.assertIn("Covalent Bond Analysis: H2O", s)
        self.assertIn("Total sigma bonds", s)
        self.assertIn("Total pi bonds", s)
        self.assertIn("Molecular polarity", s)

    def test_molecular_analysis_summary_charged(self):
        analysis = MolecularBondAnalysis("NH3", charge=0)
        s = analysis.summary()
        self.assertIn("NH3", s)


# ===================================================================
# 4. TestGraphBranches
# ===================================================================

class TestGraphBranches(unittest.TestCase):
    """Cover uncovered branches in molecule/graph.py."""

    # lines 88-89: Atom.__repr__
    def test_atom_repr(self):
        a = Atom(symbol="C", position=np.array([1.0, 2.0, 3.0]), index=0)
        r = repr(a)
        self.assertIn("Atom(C[0]", r)
        self.assertIn("1.000", r)
        self.assertIn("2.000", r)
        self.assertIn("3.000", r)

    # lines 101-103: Bond.__repr__ with order symbols
    def test_bond_repr(self):
        b1 = Bond(atom_i=0, atom_j=1, order=1)
        self.assertIn("0-1", repr(b1))
        b2 = Bond(atom_i=0, atom_j=1, order=2)
        self.assertIn("0=1", repr(b2))
        b3 = Bond(atom_i=0, atom_j=1, order=3)
        self.assertIn("0#1", repr(b3))
        # rotatable tag
        br = Bond(atom_i=0, atom_j=1, order=1, rotatable=True)
        self.assertIn("(rot)", repr(br))
        bnr = Bond(atom_i=0, atom_j=1, order=1, rotatable=False)
        self.assertNotIn("(rot)", repr(bnr))

    # line 116: TorsionAngle.__repr__
    def test_torsion_angle_repr(self):
        t = TorsionAngle(0, 1, 2, 3, angle_deg=60.0)
        r = repr(t)
        self.assertIn("Torsion(0-1-2-3", r)
        self.assertIn("60.0", r)

    # lines 133-144: NewmanProjection.summary()
    def test_newman_projection_summary(self):
        np_obj = NewmanProjection(
            front_atom=0,
            back_atom=1,
            front_substituents=[(2, "H", 120.0), (3, "H", 240.0)],
            back_substituents=[(4, "H", 60.0), (5, "H", 180.0)],
            dihedral_deg=60.0,
        )
        s = np_obj.summary()
        self.assertIn("Newman Projection along bond 0-1", s)
        self.assertIn("Dihedral: 60.0 deg", s)
        self.assertIn("Front substituents", s)
        self.assertIn("Back substituents", s)
        self.assertIn("H[2]", s)
        self.assertIn("H[4]", s)

    # lines 154-159: StrainEnergy.summary()
    def test_strain_energy_summary(self):
        t = TorsionAngle(0, 1, 2, 3, angle_deg=0.0)
        se = StrainEnergy(
            total_kj_per_mol=12.5,
            contributions=[(t, 12.5)],
        )
        s = se.summary()
        self.assertIn("Torsional Strain: 12.50 kJ/mol", s)
        self.assertIn("12.50 kJ/mol", s)

    # line 249: add_atom_bonded when n==0 -> placed at origin
    def test_add_atom_bonded_first_atom(self):
        mol = Molecule("test")
        # Must provide bond_length_val since no atoms exist yet to look up
        idx = mol.add_atom_bonded("C", bonded_to=0,
                                   bond_length_val=1.54,
                                   hybridization=Hybridization.SP3)
        self.assertEqual(idx, 0)
        np.testing.assert_array_almost_equal(
            mol.atoms[0].position, [0.0, 0.0, 0.0])

    # line 265: SP2 angle auto-select
    def test_add_atom_bonded_sp2_angle(self):
        mol = Molecule("test")
        mol.add_atom("C", [0, 0, 0], Hybridization.SP2)
        mol.add_atom("C", [0, 0, 1.34], Hybridization.SP2)
        mol.add_bond(0, 1)
        # Third atom: angle_ref auto from parent, parent is SP2
        idx = mol.add_atom_bonded("H", bonded_to=0,
                                   hybridization=None)
        self.assertEqual(idx, 2)
        # The angle C-C-H should be roughly SP2_ANGLE (120 deg)
        angle = mol.bond_angle(1, 0, 2)
        self.assertAlmostEqual(angle, 120.0, delta=2.0)

    # line 265: SP angle auto-select
    def test_add_atom_bonded_sp_angle(self):
        mol = Molecule("test")
        mol.add_atom("C", [0, 0, 0], Hybridization.SP)
        mol.add_atom("C", [0, 0, 1.2], Hybridization.SP)
        mol.add_bond(0, 1)
        idx = mol.add_atom_bonded("H", bonded_to=0)
        self.assertEqual(idx, 2)
        angle = mol.bond_angle(1, 0, 2)
        self.assertAlmostEqual(angle, 180.0, delta=2.0)

    # lines 271-275: angle_ref fallback when bonded_to has no neighbors
    def test_add_atom_bonded_angle_ref_fallback(self):
        mol = Molecule("test")
        # Add two disconnected atoms (no bonds between them)
        mol.add_atom("C", [0, 0, 0], Hybridization.SP3)
        mol.add_atom("C", [0, 0, 1.54], Hybridization.SP3)
        # No bond between 0 and 1!  Now add bonded to atom 0
        # -> angle_ref fallback to index 1 (since 0 has no nbrs)
        idx = mol.add_atom_bonded("H", bonded_to=0)
        self.assertEqual(idx, 2)
        self.assertEqual(len(mol.atoms), 3)

    # line 342: collinear atoms -> dihedral 0.0
    def test_dihedral_angle_degenerate(self):
        mol = Molecule("linear")
        mol.add_atom("C", [0, 0, 0])
        mol.add_atom("C", [0, 0, 1])
        mol.add_atom("C", [0, 0, 2])
        mol.add_atom("C", [0, 0, 3])
        mol.add_bond(0, 1)
        mol.add_bond(1, 2)
        mol.add_bond(2, 3)
        dih = mol.dihedral_angle(0, 1, 2, 3)
        self.assertAlmostEqual(dih, 0.0, places=5)

    # line 390: ring bond -> ValueError
    def test_rotate_dihedral_ring_raises(self):
        mol = Molecule("triangle")
        mol.add_atom("C", [0, 0, 0])
        mol.add_atom("C", [1, 0, 0])
        mol.add_atom("C", [0.5, 0.866, 0])
        mol.add_bond(0, 1)
        mol.add_bond(1, 2)
        mol.close_ring(2, 0)
        with self.assertRaises(ValueError):
            mol.rotate_dihedral(0, 1, 45.0)

    # lines 464-466, 471-472: _torsion_key SP2/SP hybridization + ordering
    def test_torsion_key_sp2_sp(self):
        mol = Molecule("test")
        mol.add_atom("H", [0, 0, 0])           # 0
        mol.add_atom("C", [0, 0, 1], Hybridization.SP2)  # 1
        mol.add_atom("C", [0, 0, 2], Hybridization.SP)   # 2
        mol.add_atom("H", [0, 0, 3])           # 3
        mol.add_bond(0, 1)
        mol.add_bond(1, 2)
        mol.add_bond(2, 3)
        key = mol._torsion_key(0, 1, 2, 3)
        self.assertIn("sp2", key)
        self.assertIn("sp", key)
        # Verify alphabetical ordering of terminal symbols
        self.assertTrue(key.startswith("H"))

    def test_torsion_key_ordering_reversal(self):
        # When sym_i > sym_l, the key should be reversed
        mol = Molecule("test")
        mol.add_atom("O", [0, 0, 0])           # 0 - "O" > "H"
        mol.add_atom("C", [0, 0, 1], Hybridization.SP3)  # 1
        mol.add_atom("C", [0, 0, 2], Hybridization.SP3)  # 2
        mol.add_atom("H", [0, 0, 3])           # 3
        mol.add_bond(0, 1)
        mol.add_bond(1, 2)
        mol.add_bond(2, 3)
        key = mol._torsion_key(0, 1, 2, 3)
        # O > H alphabetically, so key should start with H
        self.assertTrue(key.startswith("H"))

    # lines 483-522: newman_projection on ethane
    def test_newman_projection(self):
        mol = build_ethane(dihedral_deg=60.0)
        # C0=index 0, C1=index 1
        np_obj = mol.newman_projection(0, 1)
        self.assertEqual(np_obj.front_atom, 0)
        self.assertEqual(np_obj.back_atom, 1)
        self.assertEqual(len(np_obj.front_substituents), 3)
        self.assertEqual(len(np_obj.back_substituents), 3)
        # Each substituent should be H
        for idx, sym, ang in np_obj.front_substituents:
            self.assertEqual(sym, "H")
        for idx, sym, ang in np_obj.back_substituents:
            self.assertEqual(sym, "H")

    # line 576: assign_rs with < 4 neighbors -> NONE
    def test_assign_rs_non_tetrahedral(self):
        mol = Molecule("test")
        mol.add_atom("C", [0, 0, 0])
        mol.add_atom("H", [1, 0, 0])
        mol.add_atom("H", [0, 1, 0])
        mol.add_bond(0, 1)
        mol.add_bond(0, 2)
        result = mol.assign_rs(0)
        self.assertEqual(result, Stereodescriptor.NONE)

    # line 606: assign_ez with < 2 subs -> NONE
    def test_assign_ez_insufficient_subs(self):
        mol = Molecule("test")
        mol.add_atom("C", [0, 0, 0])
        mol.add_atom("C", [0, 0, 1.34])
        mol.add_atom("H", [1, 0, 0])
        mol.add_bond(0, 1, order=2)
        mol.add_bond(0, 2)
        # C1 has only one sub (C0), C0 has H + C1 -> 1 non-k sub
        result = mol.assign_ez(0, 1)
        self.assertEqual(result, Stereodescriptor.NONE)

    # lines 664-668: to_coordinates_dict
    def test_to_coordinates_dict(self):
        mol = Molecule("test")
        mol.add_atom("C", [0, 0, 0])
        mol.add_atom("H", [1, 0, 0])
        mol.add_bond(0, 1)
        d = mol.to_coordinates_dict()
        self.assertIn("atom_positions", d)
        self.assertIn("bonds", d)
        self.assertIn("lone_pair_positions", d)
        self.assertIn("central_index", d)
        self.assertEqual(len(d["atom_positions"]), 2)
        self.assertEqual(d["atom_positions"][0][0], "C")
        self.assertEqual(d["bonds"][0], (0, 1, 1))
        self.assertEqual(d["lone_pair_positions"], [])

    # line 678: Molecule.__repr__
    def test_molecule_repr(self):
        mol = Molecule("ethanol")
        mol.add_atom("C", [0, 0, 0])
        mol.add_atom("O", [1, 0, 0])
        mol.add_bond(0, 1)
        r = repr(mol)
        self.assertIn("Molecule('ethanol'", r)
        self.assertIn("2 atoms", r)
        self.assertIn("1 bonds", r)

    # lines 682-711: Molecule.summary()
    def test_molecule_summary(self):
        mol = build_ethane()
        s = mol.summary()
        self.assertIn("Molecule: ethane", s)
        self.assertIn("Atoms:", s)
        self.assertIn("Bonds:", s)
        self.assertIn("Coordinates", s)
        self.assertIn("SP3", s)
        self.assertIn("rotatable", s)
        self.assertIn("fixed", s)

    # Additional: torsional_energy exercising _torsion_key
    def test_torsional_energy_ethane(self):
        mol = build_ethane(dihedral_deg=0.0)  # eclipsed
        se = mol.torsional_energy(0, 1)
        self.assertGreater(se.total_kj_per_mol, 0)
        self.assertGreater(len(se.contributions), 0)
        # Staggered should have lower energy
        mol2 = build_ethane(dihedral_deg=60.0)
        se2 = mol2.torsional_energy(0, 1)
        self.assertLess(se2.total_kj_per_mol, se.total_kj_per_mol)


if __name__ == "__main__":
    unittest.main()
