"""Tests for the molbuilder.molecule package.

Covers:
    - graph: Molecule class (add_atom, add_bond, neighbors, get_bond, distance, bond_angle)
    - builders: build_ethane, build_butane, build_cyclohexane, build_chiral_molecule, build_2_butene
    - conformations: classify_conformation, scan_torsion
    - amino_acids: all 20 AminoAcidType members, chirality, peptide bond formation
    - functional_groups: add_hydroxyl, add_amino
"""

import math
import unittest

import numpy as np

from molbuilder.molecule.graph import (
    Molecule,
    Hybridization,
    Stereodescriptor,
    ConformationType,
)
from molbuilder.molecule.builders import (
    build_ethane,
    build_butane,
    build_cyclohexane,
    build_2_butene,
    build_chiral_molecule,
)
from molbuilder.molecule.conformations import classify_conformation, scan_torsion
from molbuilder.molecule.amino_acids import (
    AminoAcidType,
    build_amino_acid,
    form_peptide_bond,
)
from molbuilder.molecule.functional_groups import add_hydroxyl, add_amino


# ===================================================================
# Molecule graph tests
# ===================================================================

class TestMoleculeAddAtomAndBond(unittest.TestCase):
    """Test basic Molecule operations: add_atom, add_bond."""

    def setUp(self):
        self.mol = Molecule("test")

    def test_add_atom_returns_index(self):
        idx = self.mol.add_atom("C", [0.0, 0.0, 0.0])
        self.assertEqual(idx, 0)

    def test_add_multiple_atoms(self):
        i0 = self.mol.add_atom("C", [0.0, 0.0, 0.0])
        i1 = self.mol.add_atom("H", [1.0, 0.0, 0.0])
        self.assertEqual(i0, 0)
        self.assertEqual(i1, 1)
        self.assertEqual(len(self.mol.atoms), 2)

    def test_atom_symbol(self):
        self.mol.add_atom("O", [0.0, 0.0, 0.0])
        self.assertEqual(self.mol.atoms[0].symbol, "O")

    def test_add_bond(self):
        self.mol.add_atom("C", [0.0, 0.0, 0.0])
        self.mol.add_atom("H", [1.0, 0.0, 0.0])
        bond = self.mol.add_bond(0, 1)
        self.assertEqual(len(self.mol.bonds), 1)
        self.assertEqual(bond.atom_i, 0)
        self.assertEqual(bond.atom_j, 1)

    def test_add_bond_with_order(self):
        self.mol.add_atom("C", [0.0, 0.0, 0.0])
        self.mol.add_atom("O", [1.2, 0.0, 0.0])
        bond = self.mol.add_bond(0, 1, order=2)
        self.assertEqual(bond.order, 2)


class TestMoleculeNeighborsAndGetBond(unittest.TestCase):
    """Test neighbors() and get_bond() methods."""

    def setUp(self):
        self.mol = Molecule("methane-like")
        self.mol.add_atom("C", [0.0, 0.0, 0.0])
        self.mol.add_atom("H", [1.0, 0.0, 0.0])
        self.mol.add_atom("H", [0.0, 1.0, 0.0])
        self.mol.add_bond(0, 1)
        self.mol.add_bond(0, 2)

    def test_neighbors(self):
        nbrs = self.mol.neighbors(0)
        self.assertIn(1, nbrs)
        self.assertIn(2, nbrs)
        self.assertEqual(len(nbrs), 2)

    def test_neighbors_leaf(self):
        nbrs = self.mol.neighbors(1)
        self.assertEqual(nbrs, [0])

    def test_get_bond_exists(self):
        bond = self.mol.get_bond(0, 1)
        self.assertIsNotNone(bond)
        self.assertEqual(bond.order, 1)

    def test_get_bond_reversed(self):
        bond = self.mol.get_bond(1, 0)
        self.assertIsNotNone(bond)

    def test_get_bond_nonexistent(self):
        bond = self.mol.get_bond(1, 2)
        self.assertIsNone(bond)


class TestMoleculeDistanceAndAngle(unittest.TestCase):
    """Test distance() and bond_angle() methods."""

    def setUp(self):
        self.mol = Molecule("angle test")
        # Triangle: right angle at vertex 1
        self.mol.add_atom("C", [0.0, 0.0, 0.0])  # 0
        self.mol.add_atom("C", [1.0, 0.0, 0.0])  # 1
        self.mol.add_atom("C", [1.0, 1.0, 0.0])  # 2
        self.mol.add_bond(0, 1)
        self.mol.add_bond(1, 2)

    def test_distance(self):
        d = self.mol.distance(0, 1)
        self.assertAlmostEqual(d, 1.0, places=5)

    def test_distance_diagonal(self):
        d = self.mol.distance(0, 2)
        self.assertAlmostEqual(d, math.sqrt(2.0), places=5)

    def test_bond_angle(self):
        angle = self.mol.bond_angle(0, 1, 2)
        self.assertAlmostEqual(angle, 90.0, places=3)


# ===================================================================
# Builder tests
# ===================================================================

class TestBuildEthane(unittest.TestCase):
    """Test build_ethane: C2H6 should have 8 atoms and 7 bonds."""

    def test_atom_count(self):
        mol = build_ethane()
        self.assertEqual(len(mol.atoms), 8)

    def test_bond_count(self):
        mol = build_ethane()
        self.assertEqual(len(mol.bonds), 7)

    def test_carbon_count(self):
        mol = build_ethane()
        c_count = sum(1 for a in mol.atoms if a.symbol == "C")
        self.assertEqual(c_count, 2)

    def test_hydrogen_count(self):
        mol = build_ethane()
        h_count = sum(1 for a in mol.atoms if a.symbol == "H")
        self.assertEqual(h_count, 6)


class TestBuildButane(unittest.TestCase):
    """Test build_butane: C4H10 should have 14 atoms."""

    def test_atom_count(self):
        mol = build_butane()
        self.assertEqual(len(mol.atoms), 14)

    def test_carbon_count(self):
        mol = build_butane()
        c_count = sum(1 for a in mol.atoms if a.symbol == "C")
        self.assertEqual(c_count, 4)

    def test_hydrogen_count(self):
        mol = build_butane()
        h_count = sum(1 for a in mol.atoms if a.symbol == "H")
        self.assertEqual(h_count, 10)


class TestBuildCyclohexane(unittest.TestCase):
    """Test build_cyclohexane: C6H12 should have 18 atoms and 18 bonds."""

    def test_atom_count(self):
        mol = build_cyclohexane()
        self.assertEqual(len(mol.atoms), 18)

    def test_bond_count(self):
        mol = build_cyclohexane()
        self.assertEqual(len(mol.bonds), 18)

    def test_carbon_count(self):
        mol = build_cyclohexane()
        c_count = sum(1 for a in mol.atoms if a.symbol == "C")
        self.assertEqual(c_count, 6)

    def test_hydrogen_count(self):
        mol = build_cyclohexane()
        h_count = sum(1 for a in mol.atoms if a.symbol == "H")
        self.assertEqual(h_count, 12)


class TestBuildChiralMolecule(unittest.TestCase):
    """Test build_chiral_molecule: CHFClBr with R configuration."""

    def test_atom_count(self):
        mol = build_chiral_molecule()
        self.assertEqual(len(mol.atoms), 5)

    def test_has_four_different_substituents(self):
        mol = build_chiral_molecule()
        symbols = {mol.atoms[i].symbol for i in mol.neighbors(0)}
        self.assertEqual(len(symbols), 4)

    def test_r_configuration(self):
        mol = build_chiral_molecule()
        config = mol.assign_rs(0)
        self.assertEqual(config, Stereodescriptor.R)


class TestBuild2Butene(unittest.TestCase):
    """Test build_2_butene with is_cis parameter."""

    def test_cis_z_configuration(self):
        mol = build_2_butene(is_cis=True)
        # Find the double bond (C1=C2)
        double_bond_found = False
        for bond in mol.bonds:
            if bond.order == 2:
                j, k = bond.atom_i, bond.atom_j
                config = mol.assign_ez(j, k)
                self.assertEqual(config, Stereodescriptor.Z)
                double_bond_found = True
                break
        self.assertTrue(double_bond_found, "No double bond found in 2-butene")

    def test_trans_e_configuration(self):
        mol = build_2_butene(is_cis=False)
        double_bond_found = False
        for bond in mol.bonds:
            if bond.order == 2:
                j, k = bond.atom_i, bond.atom_j
                config = mol.assign_ez(j, k)
                self.assertEqual(config, Stereodescriptor.E)
                double_bond_found = True
                break
        self.assertTrue(double_bond_found, "No double bond found in 2-butene")

    def test_atom_counts(self):
        mol = build_2_butene(is_cis=True)
        c_count = sum(1 for a in mol.atoms if a.symbol == "C")
        h_count = sum(1 for a in mol.atoms if a.symbol == "H")
        self.assertEqual(c_count, 4)
        self.assertEqual(h_count, 8)


# ===================================================================
# Conformations tests
# ===================================================================

class TestClassifyConformation(unittest.TestCase):
    """Test classify_conformation for standard dihedral angles."""

    def test_eclipsed_at_0(self):
        self.assertEqual(classify_conformation(0.0), ConformationType.ECLIPSED)

    def test_gauche_at_60(self):
        self.assertEqual(classify_conformation(60.0), ConformationType.GAUCHE)

    def test_eclipsed_at_120(self):
        self.assertEqual(classify_conformation(120.0), ConformationType.ECLIPSED)

    def test_anti_at_180(self):
        self.assertEqual(classify_conformation(180.0), ConformationType.ANTI)

    def test_gauche_at_minus_60(self):
        self.assertEqual(classify_conformation(-60.0), ConformationType.GAUCHE)


class TestScanTorsion(unittest.TestCase):
    """Test scan_torsion returns the correct number of points."""

    def test_scan_ethane(self):
        mol = build_ethane()
        # Find the C-C bond (atoms 0 and 1)
        c0, c1 = 0, 1
        # Find an H on each carbon for ref
        h_on_c0 = [i for i in mol.neighbors(c0) if mol.atoms[i].symbol == "H"]
        h_on_c1 = [i for i in mol.neighbors(c1) if mol.atoms[i].symbol == "H"]
        self.assertTrue(len(h_on_c0) > 0 and len(h_on_c1) > 0)

        results = scan_torsion(mol, c0, c1,
                               ref_i=h_on_c0[0], ref_l=h_on_c1[0],
                               steps=36)
        self.assertEqual(len(results), 36)

    def test_scan_returns_angle_energy_pairs(self):
        mol = build_ethane()
        c0, c1 = 0, 1
        h_on_c0 = [i for i in mol.neighbors(c0) if mol.atoms[i].symbol == "H"]
        h_on_c1 = [i for i in mol.neighbors(c1) if mol.atoms[i].symbol == "H"]

        results = scan_torsion(mol, c0, c1,
                               ref_i=h_on_c0[0], ref_l=h_on_c1[0],
                               steps=12)
        self.assertEqual(len(results), 12)
        for angle, energy in results:
            self.assertIsInstance(angle, float)
            self.assertIsInstance(energy, float)


# ===================================================================
# Amino acid tests
# ===================================================================

class TestAllAminoAcidsBuild(unittest.TestCase):
    """Test that all 20 standard amino acids build successfully."""

    def test_all_20_build(self):
        for aa_type in AminoAcidType:
            with self.subTest(aa=aa_type.name):
                mol = build_amino_acid(aa_type)
                self.assertIsNotNone(mol)
                self.assertGreater(len(mol.atoms), 0)
                self.assertTrue(hasattr(mol, "backbone"))


class TestAminoAcidChirality(unittest.TestCase):
    """All non-glycine amino acids should have S chirality at C-alpha."""

    def test_non_glycine_s_chirality(self):
        for aa_type in AminoAcidType:
            if aa_type == AminoAcidType.GLY:
                continue
            with self.subTest(aa=aa_type.name):
                mol = build_amino_acid(aa_type)
                ca = mol.backbone.CA
                config = mol.assign_rs(ca)
                self.assertEqual(
                    config, Stereodescriptor.S,
                    f"{aa_type.name}: expected S, got {config}")


class TestPeptideBond(unittest.TestCase):
    """Test that peptide bonds can be formed."""

    def test_dipeptide_builds(self):
        ala = build_amino_acid(AminoAcidType.ALA)
        gly = build_amino_acid(AminoAcidType.GLY)
        dipeptide = form_peptide_bond(ala, gly)
        self.assertIsNotNone(dipeptide)
        self.assertGreater(len(dipeptide.atoms), 0)

    def test_dipeptide_has_more_atoms_than_monomer(self):
        ala1 = build_amino_acid(AminoAcidType.ALA)
        ala2 = build_amino_acid(AminoAcidType.ALA)
        n_ala = len(ala1.atoms)
        dipeptide = form_peptide_bond(ala1, ala2)
        # Dipeptide loses H2O (3 atoms) from the two monomers
        # So combined - 3 removed atoms
        self.assertGreater(len(dipeptide.atoms), n_ala)

    def test_peptide_bond_length(self):
        ala = build_amino_acid(AminoAcidType.ALA)
        gly = build_amino_acid(AminoAcidType.GLY)
        dipeptide = form_peptide_bond(ala, gly)
        # Find the peptide bond: the C-N bond between residues
        residues = dipeptide.residues
        peptide_c = residues[0].C
        peptide_n = residues[1].N
        bl = dipeptide.distance(peptide_c, peptide_n)
        # Peptide bond length should be approximately 1.33 Angstroms
        self.assertAlmostEqual(bl, 1.33, delta=0.15)


class TestAminoAcidBackbone(unittest.TestCase):
    """Test backbone indices are correctly set."""

    def test_alanine_backbone_indices(self):
        mol = build_amino_acid(AminoAcidType.ALA)
        bb = mol.backbone
        self.assertIsNotNone(bb.N)
        self.assertIsNotNone(bb.CA)
        self.assertIsNotNone(bb.C)
        self.assertIsNotNone(bb.O)
        self.assertIsNotNone(bb.CB)

    def test_glycine_no_cb(self):
        mol = build_amino_acid(AminoAcidType.GLY)
        bb = mol.backbone
        self.assertIsNone(bb.CB)

    def test_backbone_atom_symbols(self):
        mol = build_amino_acid(AminoAcidType.ALA)
        bb = mol.backbone
        self.assertEqual(mol.atoms[bb.N].symbol, "N")
        self.assertEqual(mol.atoms[bb.CA].symbol, "C")
        self.assertEqual(mol.atoms[bb.C].symbol, "C")
        self.assertEqual(mol.atoms[bb.O].symbol, "O")


# ===================================================================
# Functional group tests
# ===================================================================

class TestAddHydroxyl(unittest.TestCase):
    """Test that add_hydroxyl attaches an O atom."""

    def test_adds_oxygen(self):
        mol = Molecule("test-OH")
        mol.add_atom("C", [0.0, 0.0, 0.0], Hybridization.SP3)
        mol.add_atom("H", [1.0, 0.0, 0.0])
        mol.add_bond(0, 1)

        initial_atoms = len(mol.atoms)
        result = add_hydroxyl(mol, 0)
        self.assertIn("O", result)
        self.assertEqual(mol.atoms[result["O"]].symbol, "O")
        # Should add 2 atoms: O and H
        self.assertEqual(len(mol.atoms), initial_atoms + 2)

    def test_returns_indices(self):
        mol = Molecule("test-OH")
        mol.add_atom("C", [0.0, 0.0, 0.0], Hybridization.SP3)
        mol.add_atom("H", [1.0, 0.0, 0.0])
        mol.add_bond(0, 1)

        result = add_hydroxyl(mol, 0)
        self.assertIn("O", result)
        self.assertIn("H", result)


class TestAddAmino(unittest.TestCase):
    """Test that add_amino (amino group) attaches an N atom."""

    def test_adds_nitrogen(self):
        mol = Molecule("test-NH2")
        mol.add_atom("C", [0.0, 0.0, 0.0], Hybridization.SP3)
        mol.add_atom("H", [1.0, 0.0, 0.0])
        mol.add_bond(0, 1)

        initial_atoms = len(mol.atoms)
        result = add_amino(mol, 0)
        self.assertIn("N", result)
        self.assertEqual(mol.atoms[result["N"]].symbol, "N")
        # Should add 3 atoms: N, H1, H2
        self.assertEqual(len(mol.atoms), initial_atoms + 3)

    def test_returns_indices(self):
        mol = Molecule("test-NH2")
        mol.add_atom("C", [0.0, 0.0, 0.0], Hybridization.SP3)
        mol.add_atom("H", [1.0, 0.0, 0.0])
        mol.add_bond(0, 1)

        result = add_amino(mol, 0)
        self.assertIn("N", result)
        self.assertIn("H1", result)
        self.assertIn("H2", result)


if __name__ == "__main__":
    unittest.main()
