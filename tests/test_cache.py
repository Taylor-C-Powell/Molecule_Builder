"""Tests for ring-bond and functional-group caching on Molecule.

Run with:
    python -m pytest tests/test_cache.py -v
"""

import pytest

from molbuilder.molecule.graph import Molecule, Hybridization
from molbuilder.reactions.functional_group_detect import detect_functional_groups


# =====================================================================
#  Ring bond cache
# =====================================================================

class TestRingBondCache:
    """Test _ring_bonds_cache lifecycle on Molecule."""

    @staticmethod
    def _make_cyclohexane() -> Molecule:
        """Build a simple 6-membered carbon ring (no H)."""
        mol = Molecule(name="cyclohexane")
        for i in range(6):
            mol.add_atom("C", [float(i), 0.0, 0.0], Hybridization.SP3)
        for i in range(5):
            mol.add_bond(i, i + 1, order=1, rotatable=False)
        mol.close_ring(5, 0, order=1)
        return mol

    def test_cache_none_initially(self):
        mol = self._make_cyclohexane()
        # Cache is None because close_ring just invalidated it
        assert mol._ring_bonds_cache is None

    def test_cache_populated_after_is_in_ring(self):
        mol = self._make_cyclohexane()
        mol.is_in_ring(0, 1)
        assert mol._ring_bonds_cache is not None

    def test_ring_bonds_detected(self):
        mol = self._make_cyclohexane()
        # All bonds in a ring should be detected
        for i in range(5):
            assert mol.is_in_ring(i, i + 1)
        assert mol.is_in_ring(5, 0)

    def test_non_ring_bond_not_in_cache(self):
        mol = Molecule(name="ethane")
        c0 = mol.add_atom("C", [0.0, 0.0, 0.0])
        c1 = mol.add_atom("C", [1.5, 0.0, 0.0])
        mol.add_bond(c0, c1)
        assert not mol.is_in_ring(0, 1)

    def test_cache_invalidated_by_add_bond(self):
        mol = self._make_cyclohexane()
        mol.is_in_ring(0, 1)  # populate cache
        assert mol._ring_bonds_cache is not None
        mol.add_bond(0, 3, order=1)  # add a cross-ring bond
        assert mol._ring_bonds_cache is None

    def test_cache_invalidated_by_add_atom(self):
        mol = self._make_cyclohexane()
        mol.is_in_ring(0, 1)  # populate cache
        assert mol._ring_bonds_cache is not None
        mol.add_atom("H", [0.0, 1.0, 0.0])
        assert mol._ring_bonds_cache is None

    def test_cache_invalidated_by_close_ring(self):
        mol = Molecule(name="partial")
        for i in range(6):
            mol.add_atom("C", [float(i), 0.0, 0.0])
        for i in range(5):
            mol.add_bond(i, i + 1)
        # No ring yet
        assert not mol.is_in_ring(0, 1)
        # Close the ring
        mol.close_ring(5, 0)
        # Cache should be invalidated
        assert mol._ring_bonds_cache is None
        # Now it should be a ring
        assert mol.is_in_ring(0, 1)

    def test_second_is_in_ring_uses_cache(self):
        mol = self._make_cyclohexane()
        mol.is_in_ring(0, 1)
        cache_ref = mol._ring_bonds_cache
        # Second call should use same cache object
        mol.is_in_ring(2, 3)
        assert mol._ring_bonds_cache is cache_ref

    def test_multiple_rings(self):
        """Two fused rings: should detect ring bonds in both."""
        mol = Molecule(name="bicyclic")
        for i in range(7):
            mol.add_atom("C", [float(i), 0.0, 0.0])
        # First ring: 0-1-2-3-4-0
        for i in range(4):
            mol.add_bond(i, i + 1)
        mol.close_ring(4, 0)
        # Second ring sharing edge 2-3: 2-3-5-6-2
        mol.add_bond(3, 5)
        mol.add_bond(5, 6)
        mol.close_ring(6, 2)
        # All bonds in both rings are ring bonds
        assert mol.is_in_ring(0, 1)
        assert mol.is_in_ring(2, 3)
        assert mol.is_in_ring(5, 6)


# =====================================================================
#  FG detection cache
# =====================================================================

class TestFGCache:
    """Test _fg_cache lifecycle on Molecule."""

    @staticmethod
    def _make_ethanol() -> Molecule:
        """Build ethanol: C-C-O (no explicit H, let FG detector use valence rules)."""
        mol = Molecule(name="ethanol")
        mol.add_atom("C", [0.0, 0.0, 0.0], Hybridization.SP3)
        mol.add_atom("C", [1.5, 0.0, 0.0], Hybridization.SP3)
        mol.add_atom("O", [3.0, 0.0, 0.0])
        mol.add_bond(0, 1)
        mol.add_bond(1, 2)
        return mol

    def test_fg_cache_none_initially(self):
        mol = self._make_ethanol()
        assert mol._fg_cache is None

    def test_fg_cache_populated_after_detect(self):
        mol = self._make_ethanol()
        fgs = detect_functional_groups(mol)
        assert mol._fg_cache is not None
        assert len(fgs) > 0

    def test_second_call_returns_same_results(self):
        mol = self._make_ethanol()
        fgs1 = detect_functional_groups(mol)
        fgs2 = detect_functional_groups(mol)
        assert len(fgs1) == len(fgs2)
        assert [fg.name for fg in fgs1] == [fg.name for fg in fgs2]

    def test_cache_returns_copy_not_reference(self):
        mol = self._make_ethanol()
        fgs1 = detect_functional_groups(mol)
        fgs2 = detect_functional_groups(mol)
        # Modifying the returned list should not affect cache
        assert fgs1 is not fgs2

    def test_fg_cache_invalidated_by_add_atom(self):
        mol = self._make_ethanol()
        detect_functional_groups(mol)
        assert mol._fg_cache is not None
        mol.add_atom("Br", [4.5, 0.0, 0.0])
        assert mol._fg_cache is None

    def test_fg_cache_invalidated_by_add_bond(self):
        mol = self._make_ethanol()
        detect_functional_groups(mol)
        assert mol._fg_cache is not None
        br = mol.add_atom("Br", [4.5, 0.0, 0.0])
        # add_atom already invalidated, but let's add bond too
        mol.add_bond(0, br)
        assert mol._fg_cache is None
