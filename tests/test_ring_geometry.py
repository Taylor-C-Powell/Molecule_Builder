"""Tests for ring geometry correction in SMILES-parsed molecules.

Validates that ring structures have correct planarity, bond lengths,
internal angles, and that substituents are properly repositioned.
"""

import math

import numpy as np
import pytest

from molbuilder.smiles.parser import parse
from molbuilder.smiles.ring_geometry import find_all_rings, correct_ring_geometry


# ===================================================================
# Helper functions
# ===================================================================

def _ring_heavy_atoms(mol):
    """Return indices of heavy (non-H) atoms."""
    return [a.index for a in mol.atoms if a.symbol != "H"]


def _max_plane_deviation(mol, indices):
    """Max deviation of atoms from their best-fit plane (Angstroms)."""
    coords = np.array([mol.atoms[i].position for i in indices])
    centroid = coords.mean(axis=0)
    centered = coords - centroid

    # SVD to find the normal of the best-fit plane
    _, S, Vt = np.linalg.svd(centered)
    normal = Vt[-1]  # smallest singular value direction

    deviations = np.abs(centered @ normal)
    return float(deviations.max())


def _ring_bond_lengths(mol, ring):
    """Return list of bond lengths for consecutive ring atoms."""
    n = len(ring)
    lengths = []
    for k in range(n):
        i = ring[k]
        j = ring[(k + 1) % n]
        d = float(np.linalg.norm(
            mol.atoms[i].position - mol.atoms[j].position))
        lengths.append(d)
    return lengths


def _ring_internal_angles(mol, ring):
    """Return list of internal angles at each ring vertex (degrees)."""
    n = len(ring)
    angles = []
    for k in range(n):
        i = ring[(k - 1) % n]
        j = ring[k]
        l = ring[(k + 1) % n]
        angles.append(mol.bond_angle(i, j, l))
    return angles


def _find_ring_of_size(mol, size):
    """Find a ring of a specific size in the molecule."""
    rings = find_all_rings(mol)
    for ring in rings:
        if len(ring) == size:
            return ring
    return None


# ===================================================================
# Benzene tests
# ===================================================================

class TestBenzene:
    """Tests for benzene (c1ccccc1) ring geometry."""

    @pytest.fixture
    def benzene(self):
        return parse("c1ccccc1")

    def test_benzene_planarity(self, benzene):
        """All 6 carbon atoms should be coplanar (< 0.01 A deviation)."""
        carbons = [a.index for a in benzene.atoms if a.symbol == "C"]
        assert len(carbons) == 6
        deviation = _max_plane_deviation(benzene, carbons)
        assert deviation < 0.01, f"Plane deviation {deviation:.4f} A > 0.01 A"

    def test_benzene_bond_lengths(self, benzene):
        """All C-C bonds in benzene should be ~1.40 A (aromatic)."""
        ring = _find_ring_of_size(benzene, 6)
        assert ring is not None, "No 6-membered ring found"
        lengths = _ring_bond_lengths(benzene, ring)
        for bl in lengths:
            assert 1.30 < bl < 1.50, f"Bond length {bl:.3f} A out of range"

    def test_benzene_angles(self, benzene):
        """Internal angles should be ~120 degrees."""
        ring = _find_ring_of_size(benzene, 6)
        assert ring is not None
        angles = _ring_internal_angles(benzene, ring)
        for angle in angles:
            assert abs(angle - 120.0) < 5.0, (
                f"Internal angle {angle:.1f} deg, expected ~120 deg")

    def test_benzene_bond_uniformity(self, benzene):
        """All ring bond lengths should be similar (within 0.15 A)."""
        ring = _find_ring_of_size(benzene, 6)
        assert ring is not None
        lengths = _ring_bond_lengths(benzene, ring)
        assert max(lengths) - min(lengths) < 0.15, (
            f"Bond length spread {max(lengths) - min(lengths):.3f} A")


# ===================================================================
# Cyclopentane test
# ===================================================================

class TestCyclopentane:
    """Tests for cyclopentane (C1CCCC1) ring geometry."""

    @pytest.fixture
    def cyclopentane(self):
        return parse("C1CCCC1")

    def test_cyclopentane_bond_lengths(self, cyclopentane):
        """C-C bonds should be ~1.54 A."""
        ring = _find_ring_of_size(cyclopentane, 5)
        assert ring is not None, "No 5-membered ring found"
        lengths = _ring_bond_lengths(cyclopentane, ring)
        for bl in lengths:
            assert 1.35 < bl < 1.75, f"Bond length {bl:.3f} A out of range"

    def test_cyclopentane_angles(self, cyclopentane):
        """Internal angles should be ~108 degrees."""
        ring = _find_ring_of_size(cyclopentane, 5)
        assert ring is not None
        angles = _ring_internal_angles(cyclopentane, ring)
        for angle in angles:
            assert abs(angle - 108.0) < 8.0, (
                f"Internal angle {angle:.1f} deg, expected ~108 deg")


# ===================================================================
# Cyclohexane chair test
# ===================================================================

class TestCyclohexane:
    """Tests for cyclohexane (C1CCCCC1) ring geometry."""

    @pytest.fixture
    def cyclohexane(self):
        return parse("C1CCCCC1")

    def test_cyclohexane_bond_lengths(self, cyclohexane):
        """C-C bonds should be ~1.54 A."""
        ring = _find_ring_of_size(cyclohexane, 6)
        assert ring is not None, "No 6-membered ring found"
        lengths = _ring_bond_lengths(cyclohexane, ring)
        for bl in lengths:
            assert 1.40 < bl < 1.65, f"Bond length {bl:.3f} A out of range"

    def test_cyclohexane_puckered(self, cyclohexane):
        """Chair cyclohexane should NOT be planar (deviation > 0.1 A)."""
        ring = _find_ring_of_size(cyclohexane, 6)
        assert ring is not None
        deviation = _max_plane_deviation(cyclohexane, ring)
        assert deviation > 0.1, (
            f"Cyclohexane too flat: deviation {deviation:.4f} A")


# ===================================================================
# Pyridine test
# ===================================================================

class TestPyridine:
    """Tests for pyridine (c1ccncc1) ring geometry."""

    @pytest.fixture
    def pyridine(self):
        return parse("c1ccncc1")

    def test_pyridine_planarity(self, pyridine):
        """Pyridine should be planar (< 0.01 A deviation)."""
        heavy = _ring_heavy_atoms(pyridine)
        ring_atoms = [i for i in heavy if pyridine.atoms[i].symbol in ("C", "N")]
        # Find the 6 ring atoms (5C + 1N)
        ring = _find_ring_of_size(pyridine, 6)
        if ring is not None:
            deviation = _max_plane_deviation(pyridine, ring)
            assert deviation < 0.01, (
                f"Pyridine plane deviation {deviation:.4f} A > 0.01 A")

    def test_pyridine_has_nitrogen(self, pyridine):
        """Pyridine should contain one nitrogen in the ring."""
        ring = _find_ring_of_size(pyridine, 6)
        assert ring is not None
        n_count = sum(1 for i in ring if pyridine.atoms[i].symbol == "N")
        assert n_count == 1, f"Expected 1 N in ring, found {n_count}"


# ===================================================================
# Naphthalene (fused rings) test
# ===================================================================

class TestNaphthalene:
    """Tests for naphthalene (c1ccc2ccccc2c1) ring geometry."""

    @pytest.fixture
    def naphthalene(self):
        return parse("c1ccc2ccccc2c1")

    def test_naphthalene_planarity(self, naphthalene):
        """All carbons should be coplanar (< 0.05 A deviation)."""
        carbons = [a.index for a in naphthalene.atoms if a.symbol == "C"]
        assert len(carbons) == 10
        deviation = _max_plane_deviation(naphthalene, carbons)
        assert deviation < 0.05, (
            f"Naphthalene plane deviation {deviation:.4f} A > 0.05 A")

    def test_naphthalene_two_rings(self, naphthalene):
        """Should detect two 6-membered rings."""
        rings = find_all_rings(naphthalene)
        six_rings = [r for r in rings if len(r) == 6]
        assert len(six_rings) >= 2, (
            f"Expected >= 2 six-membered rings, found {len(six_rings)}")

    def test_naphthalene_bond_lengths(self, naphthalene):
        """All C-C bonds should be in aromatic range."""
        rings = find_all_rings(naphthalene)
        six_rings = [r for r in rings if len(r) == 6]
        for ring in six_rings:
            lengths = _ring_bond_lengths(naphthalene, ring)
            for bl in lengths:
                assert 1.15 < bl < 1.70, (
                    f"Bond length {bl:.3f} A out of range")


# ===================================================================
# Toluene substituent test
# ===================================================================

class TestToluene:
    """Tests for toluene (Cc1ccccc1) substituent positioning."""

    @pytest.fixture
    def toluene(self):
        return parse("Cc1ccccc1")

    def test_toluene_ring_planar(self, toluene):
        """The benzene ring in toluene should be planar."""
        ring = _find_ring_of_size(toluene, 6)
        assert ring is not None
        deviation = _max_plane_deviation(toluene, ring)
        assert deviation < 0.01, (
            f"Toluene ring deviation {deviation:.4f} A > 0.01 A")

    def test_toluene_methyl_connected(self, toluene):
        """Methyl carbon should be bonded to a ring carbon at reasonable distance."""
        # Find the sp3 carbon (methyl)
        ring = _find_ring_of_size(toluene, 6)
        assert ring is not None
        ring_set = set(ring)

        methyl = None
        for a in toluene.atoms:
            if a.symbol == "C" and a.index not in ring_set:
                methyl = a.index
                break

        assert methyl is not None, "No methyl carbon found"

        # Check that methyl is bonded to a ring atom at reasonable distance
        neighbors = toluene.neighbors(methyl)
        ring_neighbor = [n for n in neighbors if n in ring_set]
        assert len(ring_neighbor) == 1, "Methyl should bond to exactly one ring atom"

        d = float(np.linalg.norm(
            toluene.atoms[methyl].position
            - toluene.atoms[ring_neighbor[0]].position))
        assert 1.2 < d < 2.0, f"C-C(ring) distance {d:.3f} A out of range"


# ===================================================================
# Acyclic molecules unchanged test
# ===================================================================

class TestAcyclicUnchanged:
    """Verify that acyclic molecules are not affected by ring correction."""

    def test_ethane_no_rings(self):
        """Ethane should have no rings detected."""
        mol = parse("CC")
        rings = find_all_rings(mol)
        assert len(rings) == 0

    def test_butane_no_rings(self):
        """Butane should have no rings detected."""
        mol = parse("CCCC")
        rings = find_all_rings(mol)
        assert len(rings) == 0

    def test_ethanol_unchanged(self):
        """Ethanol geometry should be valid after parsing."""
        mol = parse("CCO")
        # Just verify it parsed without error and has reasonable geometry
        assert len(mol.atoms) > 0
        for a in mol.atoms:
            assert np.all(np.isfinite(a.position)), (
                f"Atom {a.index} has non-finite position")


# ===================================================================
# Ring detection tests
# ===================================================================

class TestRingDetection:
    """Tests for the find_all_rings function."""

    def test_benzene_one_ring(self):
        mol = parse("c1ccccc1")
        rings = find_all_rings(mol)
        six_rings = [r for r in rings if len(r) == 6]
        assert len(six_rings) >= 1

    def test_cyclopropane(self):
        mol = parse("C1CC1")
        rings = find_all_rings(mol)
        three_rings = [r for r in rings if len(r) == 3]
        assert len(three_rings) >= 1

    def test_no_ring_in_chain(self):
        mol = parse("CCCCC")
        rings = find_all_rings(mol)
        assert len(rings) == 0
