"""Tests for the coords module: distance geometry + FF optimization pipeline."""

import math

import numpy as np
import pytest

from molbuilder.smiles import parse
from molbuilder.molecule.graph import Molecule, Hybridization
from molbuilder.coords import generate_3d
from molbuilder.coords.builtin import generate_builtin


# ===================================================================
# Helpers
# ===================================================================

def _planarity(mol, indices):
    """Max deviation from best-fit plane through the given atom indices."""
    pos = np.array([mol.atoms[i].position for i in indices])
    centroid = pos.mean(axis=0)
    centered = pos - centroid
    _, S, Vt = np.linalg.svd(centered)
    normal = Vt[2]
    return float(np.abs(centered @ normal).max())


def _ring_carbon_indices(mol):
    """Return indices of SP2 carbons (ring carbons in aromatic systems)."""
    return [a.index for a in mol.atoms
            if a.symbol == "C" and a.hybridization == Hybridization.SP2]


# ===================================================================
# Benzene: planar, regular hexagon, aromatic bond lengths
# ===================================================================

class TestBenzene:
    @pytest.fixture
    def benzene(self):
        return parse("c1ccccc1")

    def test_atom_count(self, benzene):
        assert len(benzene.atoms) == 12  # 6C + 6H

    def test_planar(self, benzene):
        c_indices = _ring_carbon_indices(benzene)
        assert _planarity(benzene, c_indices) < 0.05

    def test_bond_lengths(self, benzene):
        carbons = [a for a in benzene.atoms if a.symbol == "C"]
        for i in range(6):
            j = (i + 1) % 6
            d = np.linalg.norm(carbons[j].position - carbons[i].position)
            assert abs(d - 1.40) < 0.05, f"C{i}-C{j} = {d:.3f}, expected ~1.40"

    def test_angles(self, benzene):
        carbons = [a for a in benzene.atoms if a.symbol == "C"]
        for i in range(6):
            a1 = carbons[(i - 1) % 6].index
            a2 = carbons[i].index
            a3 = carbons[(i + 1) % 6].index
            angle = benzene.bond_angle(a1, a2, a3)
            assert abs(angle - 120.0) < 2.0, f"Angle at C{i} = {angle:.1f}"


# ===================================================================
# Cyclohexane: puckered chair, SP3 bond lengths
# ===================================================================

class TestCyclohexane:
    @pytest.fixture
    def cyclohexane(self):
        return parse("C1CCCCC1")

    def test_atom_count(self, cyclohexane):
        assert len(cyclohexane.atoms) == 18  # 6C + 12H

    def test_bond_lengths(self, cyclohexane):
        carbons = [a for a in cyclohexane.atoms if a.symbol == "C"]
        for i in range(6):
            j = (i + 1) % 6
            d = np.linalg.norm(carbons[j].position - carbons[i].position)
            assert abs(d - 1.54) < 0.05, f"C{i}-C{j} = {d:.3f}, expected ~1.54"

    def test_puckered(self, cyclohexane):
        """Cyclohexane should NOT be planar (chair conformation)."""
        carbons = [a for a in cyclohexane.atoms if a.symbol == "C"]
        c_indices = [c.index for c in carbons]
        planarity = _planarity(cyclohexane, c_indices)
        assert planarity > 0.1, f"Expected puckered, got planarity={planarity:.3f}"

    def test_angles(self, cyclohexane):
        carbons = [a for a in cyclohexane.atoms if a.symbol == "C"]
        for i in range(6):
            a1 = carbons[(i - 1) % 6].index
            a2 = carbons[i].index
            a3 = carbons[(i + 1) % 6].index
            angle = cyclohexane.bond_angle(a1, a2, a3)
            assert abs(angle - 111.0) < 5.0, f"Angle at C{i} = {angle:.1f}"


# ===================================================================
# Ethanol: C-C and C-O bond lengths, C-C-O angle
# ===================================================================

class TestEthanol:
    @pytest.fixture
    def ethanol(self):
        return parse("CCO")

    def test_cc_bond(self, ethanol):
        d = ethanol.distance(0, 1)
        assert abs(d - 1.54) < 0.05, f"C-C = {d:.3f}"

    def test_co_bond(self, ethanol):
        d = ethanol.distance(1, 2)
        assert abs(d - 1.43) < 0.05, f"C-O = {d:.3f}"

    def test_cco_angle(self, ethanol):
        angle = ethanol.bond_angle(0, 1, 2)
        assert abs(angle - 109.5) < 5.0, f"C-C-O = {angle:.1f}"


# ===================================================================
# Water: O-H bond length, H-O-H angle
# ===================================================================

class TestWater:
    @pytest.fixture
    def water(self):
        mol = parse("O")
        generate_3d(mol, backend="builtin")
        return mol

    def test_oh_bonds(self, water):
        d1 = water.distance(0, 1)
        d2 = water.distance(0, 2)
        assert abs(d1 - 0.96) < 0.05, f"O-H = {d1:.3f}"
        assert abs(d2 - 0.96) < 0.05, f"O-H = {d2:.3f}"

    def test_hoh_angle(self, water):
        angle = water.bond_angle(1, 0, 2)
        # SP3 ideal angle is 109.47; water's actual is 104.5 (lone pair compression)
        # Our force field uses SP3 ideal, so accept 104-115 range
        assert 100.0 < angle < 115.0, f"H-O-H = {angle:.1f}"


# ===================================================================
# Aspirin: planar aromatic ring, no steric clashes
# ===================================================================

class TestAspirin:
    @pytest.fixture
    def aspirin(self):
        return parse("CC(=O)Oc1ccccc1C(=O)O")

    def test_no_clashes(self, aspirin):
        clashes = aspirin.check_steric_clashes(0.5)
        assert len(clashes) == 0, f"Found {len(clashes)} steric clashes"

    def test_aromatic_ring_planar(self, aspirin):
        # Find actual ring carbons: SP2 carbons with 2+ SP2-carbon neighbors
        ring_c = []
        all_sp2_c = set(_ring_carbon_indices(aspirin))
        for idx in all_sp2_c:
            sp2_c_nbrs = [n for n in aspirin.neighbors(idx)
                          if n in all_sp2_c]
            if len(sp2_c_nbrs) >= 2:
                ring_c.append(idx)
        if len(ring_c) >= 5:
            planarity = _planarity(aspirin, ring_c)
            assert planarity < 0.1


# ===================================================================
# Naphthalene: fused aromatic rings
# ===================================================================

class TestNaphthalene:
    @pytest.fixture
    def naphthalene(self):
        return parse("c1ccc2ccccc2c1")

    def test_planar(self, naphthalene):
        ring_c = _ring_carbon_indices(naphthalene)
        assert _planarity(naphthalene, ring_c) < 0.05

    def test_bond_lengths(self, naphthalene):
        """All aromatic C-C bonds should be ~1.40 A."""
        carbons = [a for a in naphthalene.atoms if a.symbol == "C"]
        for c in carbons:
            for j in naphthalene.neighbors(c.index):
                if naphthalene.atoms[j].symbol == "C" and j > c.index:
                    d = naphthalene.distance(c.index, j)
                    assert abs(d - 1.40) < 0.05, f"C-C bond {c.index}-{j} = {d:.3f}"


# ===================================================================
# Backend dispatch
# ===================================================================

class TestBackendDispatch:
    def test_builtin_backend(self):
        mol = parse("C")  # default uses auto -> builtin if no rdkit
        assert len(mol.atoms) == 5

    def test_explicit_builtin(self):
        from molbuilder.smiles.tokenizer import tokenize
        from molbuilder.smiles.parser import (
            _build_graph, _add_implicit_hydrogens, _determine_hybridization,
        )
        from molbuilder.smiles.aromatic import kekulize

        tokens = tokenize("C")
        atoms, bonds = _build_graph(tokens)
        kekulize(atoms, bonds)
        atoms, bonds = _add_implicit_hydrogens(atoms, bonds)
        mol = Molecule(name="C")
        for ai in atoms:
            hyb = _determine_hybridization(ai.index, atoms, bonds)
            mol.add_atom(symbol=ai.symbol, position=[0, 0, 0], hybridization=hyb)
        for bi in bonds:
            mol.add_bond(bi.atom_i, bi.atom_j, order=bi.order)

        generate_3d(mol, backend="builtin")
        for a in mol.atoms:
            assert np.all(np.isfinite(a.position))

    def test_invalid_backend(self):
        mol = Molecule(name="test")
        mol.add_atom("C", [0, 0, 0])
        with pytest.raises(ValueError, match="Unknown backend"):
            generate_3d(mol, backend="invalid")

    def test_rdkit_backend_import_error(self):
        """rdkit backend raises ImportError if rdkit not installed."""
        # This test is conditional -- it only tests the error path
        try:
            import rdkit  # noqa: F401
            pytest.skip("rdkit is installed, skipping import error test")
        except ImportError:
            mol = Molecule(name="test")
            mol.add_atom("C", [0, 0, 0])
            with pytest.raises(ImportError):
                generate_3d(mol, backend="rdkit")


# ===================================================================
# Drug molecule validation: all atoms have finite positions, no clashes
# ===================================================================

_DRUG_SMILES = [
    ("methane", "C"),
    ("ethane", "CC"),
    ("propane", "CCC"),
    ("butane", "CCCC"),
    ("methanol", "CO"),
    ("acetic_acid", "CC(=O)O"),
    ("acetone", "CC(=O)C"),
    ("diethyl_ether", "CCOCC"),
    ("toluene", "Cc1ccccc1"),
    ("phenol", "Oc1ccccc1"),
    ("aniline", "Nc1ccccc1"),
    ("pyridine", "c1ccncc1"),
    ("furan", "c1ccoc1"),
    ("thiophene", "c1ccsc1"),
    ("imidazole", "c1cnc[nH]1"),
    ("cyclopentane", "C1CCCC1"),
    ("caffeine", "Cn1c(=O)c2c(ncn2C)n(C)c1=O"),
]


@pytest.mark.parametrize("name,smiles", _DRUG_SMILES)
def test_drug_molecule_valid_geometry(name, smiles):
    """All atoms should have finite positions after coordinate generation."""
    mol = parse(smiles)
    generate_3d(mol, backend="builtin")
    for atom in mol.atoms:
        assert np.all(np.isfinite(atom.position)), (
            f"{name}: atom {atom.index} ({atom.symbol}) has non-finite position"
        )


@pytest.mark.parametrize("name,smiles", _DRUG_SMILES)
def test_drug_molecule_no_severe_clashes(name, smiles):
    """No atom pairs should be closer than 0.4 A (severe clash)."""
    mol = parse(smiles)
    generate_3d(mol, backend="builtin")
    clashes = mol.check_steric_clashes(min_distance=0.4)
    assert len(clashes) == 0, (
        f"{name}: {len(clashes)} severe steric clashes found"
    )


# ===================================================================
# Seed reproducibility
# ===================================================================

def test_seed_reproducibility():
    """Same seed should produce identical coordinates."""
    from molbuilder.smiles.tokenizer import tokenize
    from molbuilder.smiles.parser import (
        _build_graph, _add_implicit_hydrogens, _determine_hybridization,
    )
    from molbuilder.smiles.aromatic import kekulize

    def _make_mol():
        tokens = tokenize("CCO")
        atoms, bonds = _build_graph(tokens)
        kekulize(atoms, bonds)
        atoms, bonds = _add_implicit_hydrogens(atoms, bonds)
        mol = Molecule(name="CCO")
        for ai in atoms:
            hyb = _determine_hybridization(ai.index, atoms, bonds)
            mol.add_atom(symbol=ai.symbol, position=[0, 0, 0], hybridization=hyb)
        for bi in bonds:
            mol.add_bond(bi.atom_i, bi.atom_j, order=bi.order)
        return mol

    mol1 = _make_mol()
    mol2 = _make_mol()
    generate_3d(mol1, backend="builtin", seed=42)
    generate_3d(mol2, backend="builtin", seed=42)

    for a1, a2 in zip(mol1.atoms, mol2.atoms):
        np.testing.assert_allclose(a1.position, a2.position, atol=1e-10)
