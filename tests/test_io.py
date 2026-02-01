"""Unit tests for molbuilder I/O modules: XYZ, JSON, MOL/SDF, PDB."""

import os
import tempfile
import unittest

from molbuilder.molecule.builders import build_ethane, build_butane
from molbuilder.io.xyz import write_xyz, read_xyz, to_xyz_string, from_xyz_string
from molbuilder.io.json_io import write_json, read_json, to_json_string, from_json_string
from molbuilder.io.mol_sdf import write_mol, read_mol, write_sdf, read_sdf
from molbuilder.io.pdb import write_pdb, read_pdb


class TestXYZ(unittest.TestCase):
    """Tests for XYZ file round-trip."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.ethane = build_ethane()

    def test_xyz_round_trip_atom_count(self):
        """write_xyz then read_xyz should preserve the total atom count."""
        path = os.path.join(self.tmpdir, "ethane.xyz")
        write_xyz(self.ethane, path)
        mol = read_xyz(path)
        self.assertEqual(len(mol.atoms), len(self.ethane.atoms))

    def test_xyz_round_trip_elements(self):
        """write_xyz then read_xyz should preserve element symbols."""
        path = os.path.join(self.tmpdir, "ethane.xyz")
        write_xyz(self.ethane, path)
        mol = read_xyz(path)
        orig_symbols = sorted(a.symbol for a in self.ethane.atoms)
        read_symbols = sorted(a.symbol for a in mol.atoms)
        self.assertEqual(orig_symbols, read_symbols)

    def test_xyz_string_round_trip(self):
        """to_xyz_string then from_xyz_string should preserve atom count."""
        s = to_xyz_string(self.ethane)
        mol = from_xyz_string(s)
        self.assertEqual(len(mol.atoms), len(self.ethane.atoms))


class TestJSON(unittest.TestCase):
    """Tests for JSON file round-trip."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.ethane = build_ethane()
        self.butane = build_butane()

    def test_json_round_trip_atom_count(self):
        """write_json then read_json should preserve atom count."""
        path = os.path.join(self.tmpdir, "ethane.json")
        write_json(self.ethane, path)
        mol = read_json(path)
        self.assertEqual(len(mol.atoms), len(self.ethane.atoms))

    def test_json_round_trip_bond_count(self):
        """write_json then read_json should preserve bond count."""
        path = os.path.join(self.tmpdir, "ethane.json")
        write_json(self.ethane, path)
        mol = read_json(path)
        self.assertEqual(len(mol.bonds), len(self.ethane.bonds))

    def test_json_round_trip_bond_orders(self):
        """write_json then read_json should preserve bond orders."""
        path = os.path.join(self.tmpdir, "butane.json")
        write_json(self.butane, path)
        mol = read_json(path)
        orig_orders = sorted(b.order for b in self.butane.bonds)
        read_orders = sorted(b.order for b in mol.bonds)
        self.assertEqual(orig_orders, read_orders)

    def test_json_string_round_trip(self):
        """to_json_string then from_json_string should preserve atom count."""
        s = to_json_string(self.ethane)
        mol = from_json_string(s)
        self.assertEqual(len(mol.atoms), len(self.ethane.atoms))
        self.assertEqual(len(mol.bonds), len(self.ethane.bonds))


class TestMOL(unittest.TestCase):
    """Tests for MOL file round-trip."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.ethane = build_ethane()
        self.butane = build_butane()

    def test_mol_round_trip_atom_count(self):
        """write_mol then read_mol should preserve atom count."""
        path = os.path.join(self.tmpdir, "ethane.mol")
        write_mol(self.ethane, path)
        mol = read_mol(path)
        self.assertEqual(len(mol.atoms), len(self.ethane.atoms))

    def test_mol_round_trip_bond_count(self):
        """write_mol then read_mol should preserve bond count."""
        path = os.path.join(self.tmpdir, "butane.mol")
        write_mol(self.butane, path)
        mol = read_mol(path)
        self.assertEqual(len(mol.bonds), len(self.butane.bonds))


class TestSDF(unittest.TestCase):
    """Tests for SDF multi-molecule round-trip."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.ethane = build_ethane()
        self.butane = build_butane()

    def test_sdf_round_trip_molecule_count(self):
        """write_sdf with two molecules then read_sdf should return 2."""
        path = os.path.join(self.tmpdir, "multi.sdf")
        write_sdf([self.ethane, self.butane], path)
        mols = read_sdf(path)
        self.assertEqual(len(mols), 2)

    def test_sdf_preserves_atom_counts(self):
        """Each molecule read from SDF should have correct atom count."""
        path = os.path.join(self.tmpdir, "multi.sdf")
        write_sdf([self.ethane, self.butane], path)
        mols = read_sdf(path)
        self.assertEqual(len(mols[0].atoms), len(self.ethane.atoms))
        self.assertEqual(len(mols[1].atoms), len(self.butane.atoms))


class TestPDB(unittest.TestCase):
    """Tests for PDB file round-trip."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.ethane = build_ethane()

    def test_pdb_round_trip_atom_count(self):
        """write_pdb then read_pdb should preserve atom count."""
        path = os.path.join(self.tmpdir, "ethane.pdb")
        write_pdb(self.ethane, path)
        mol = read_pdb(path)
        self.assertEqual(len(mol.atoms), len(self.ethane.atoms))


if __name__ == "__main__":
    unittest.main()
