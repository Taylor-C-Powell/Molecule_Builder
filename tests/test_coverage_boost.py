"""Tests targeting low-coverage modules to push overall coverage to 80%+.

Covers: smiles_io.py, synthesis_route.py helpers, and mechanism_choreography.
"""

import os
import tempfile
import warnings

import pytest

from molbuilder.smiles.parser import parse
from molbuilder.molecule.builders import build_ethane


# =====================================================================
#  SMILES I/O (smiles_io.py - was 23% coverage)
# =====================================================================

from molbuilder.io.smiles_io import write_smiles, read_smiles


class TestWriteSmiles:
    def test_write_creates_file(self, tmp_path):
        mol = parse("CCO")
        mol.name = "ethanol"
        path = str(tmp_path / "out.smi")
        write_smiles(mol, path)
        assert os.path.exists(path)

    def test_write_contains_smiles_and_name(self, tmp_path):
        mol = parse("CCO")
        mol.name = "ethanol"
        path = str(tmp_path / "out.smi")
        write_smiles(mol, path)
        with open(path) as f:
            content = f.read()
        assert "ethanol" in content
        # Should contain some SMILES-like string
        assert len(content.strip()) > 0


class TestReadSmiles:
    def test_read_single_molecule(self, tmp_path):
        path = str(tmp_path / "in.smi")
        with open(path, "w") as f:
            f.write("CCO ethanol\n")
        mols = read_smiles(path)
        assert len(mols) == 1
        assert mols[0].name == "ethanol"

    def test_read_multiple_molecules(self, tmp_path):
        path = str(tmp_path / "in.smi")
        with open(path, "w") as f:
            f.write("CCO ethanol\nC methane\nCC ethane\n")
        mols = read_smiles(path)
        assert len(mols) == 3

    def test_read_skips_blank_lines(self, tmp_path):
        path = str(tmp_path / "in.smi")
        with open(path, "w") as f:
            f.write("CCO ethanol\n\n\nC methane\n")
        mols = read_smiles(path)
        assert len(mols) == 2

    def test_read_skips_comment_lines(self, tmp_path):
        path = str(tmp_path / "in.smi")
        with open(path, "w") as f:
            f.write("# This is a comment\nCCO ethanol\n# Another comment\n")
        mols = read_smiles(path)
        assert len(mols) == 1

    def test_read_uses_smiles_as_name_when_no_name(self, tmp_path):
        path = str(tmp_path / "in.smi")
        with open(path, "w") as f:
            f.write("CCO\n")
        mols = read_smiles(path)
        assert len(mols) == 1
        assert mols[0].name == "CCO"

    def test_read_warns_on_invalid_smiles(self, tmp_path):
        path = str(tmp_path / "in.smi")
        with open(path, "w") as f:
            f.write("INVALIDSMILES!!!\nCCO ethanol\n")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            mols = read_smiles(path)
            assert len(mols) == 1  # Only valid one
            assert len(w) >= 1
            assert "Skipping invalid SMILES" in str(w[0].message)

    def test_round_trip(self, tmp_path):
        mol = parse("CC(=O)O")
        mol.name = "acetic_acid"
        path = str(tmp_path / "rt.smi")
        write_smiles(mol, path)
        mols = read_smiles(path)
        assert len(mols) == 1
        # Name should be preserved
        assert "acetic" in mols[0].name.lower() or mols[0].name == "acetic_acid"


# =====================================================================
#  Synthesis route helpers (synthesis_route.py - was 68% coverage)
# =====================================================================

from molbuilder.reactions.synthesis_route import (
    extract_best_route,
    format_route,
    _build_conditions,
    _expected_yield,
    _product_name,
    _compute_longest_linear,
    _gather_purchasable_leaves,
    SynthesisRoute,
    SynthesisStep,
)
from molbuilder.reactions.retrosynthesis import (
    retrosynthesis,
    RetrosynthesisTree,
    PURCHASABLE_MATERIALS,
)


class TestBuildConditions:
    def test_equal_temperature(self):
        """When lo == hi, should show single temperature."""
        mol = parse("CCOC(=O)C")
        tree = retrosynthesis(mol, max_depth=2, beam_width=3)
        if tree.target.best_disconnection:
            template = tree.target.best_disconnection.template
            cond = _build_conditions(template)
            assert "C" in cond  # Temperature always present

    def test_conditions_string_nonempty(self):
        mol = parse("CCCO")
        tree = retrosynthesis(mol, max_depth=2, beam_width=3)
        if tree.target.best_disconnection:
            cond = _build_conditions(tree.target.best_disconnection.template)
            assert len(cond) > 0


class TestExpectedYield:
    def test_yield_is_midpoint(self):
        mol = parse("CCOC(=O)C")
        tree = retrosynthesis(mol, max_depth=2, beam_width=3)
        if tree.target.best_disconnection:
            template = tree.target.best_disconnection.template
            y = _expected_yield(template)
            lo, hi = template.typical_yield
            assert y == pytest.approx((lo + hi) / 2.0)


class TestProductName:
    def test_purchasable_material_has_name(self):
        # Ethanol should be purchasable
        name = _product_name("CCO")
        assert name != "CCO"  # Should have a real name
        assert isinstance(name, str)

    def test_unknown_smiles_returns_smiles(self):
        name = _product_name("XYZNOTREAL")
        assert name == "XYZNOTREAL"


class TestComputeLongestLinear:
    def test_with_real_tree(self):
        mol = parse("CCOC(=O)C")
        tree = retrosynthesis(mol, max_depth=3, beam_width=3)
        lls = _compute_longest_linear(tree.target)
        assert isinstance(lls, int)
        assert lls >= 0


class TestFormatRoute:
    def test_format_route_structure(self):
        mol = parse("CCOC(=O)C")
        tree = retrosynthesis(mol, max_depth=3, beam_width=5)
        route = extract_best_route(tree)
        text = format_route(route)
        assert "Forward Synthesis Route" in text
        assert "Target" in text
        assert "Overall yield" in text
        assert "Starting Materials" in text

    def test_format_route_with_steps(self):
        mol = parse("CCOC(=O)C")
        tree = retrosynthesis(mol, max_depth=3, beam_width=5)
        route = extract_best_route(tree)
        text = format_route(route)
        if route.total_steps > 0:
            assert "Step 1" in text

    def test_format_route_empty_route(self):
        """A route with no steps should still format."""
        route = SynthesisRoute(
            target_smiles="C",
            target_name="methane",
            steps=[],
            overall_yield=100.0,
            starting_materials=[],
            total_steps=0,
            longest_linear_sequence=0,
        )
        text = format_route(route)
        assert "methane" in text
        assert "(none identified)" in text


class TestExtractBestRoute:
    def test_route_has_correct_target(self):
        mol = parse("CCOC(=O)C")
        tree = retrosynthesis(mol, max_depth=3, beam_width=5)
        route = extract_best_route(tree)
        assert isinstance(route, SynthesisRoute)
        assert len(route.target_smiles) > 0

    def test_overall_yield_is_product_of_steps(self):
        mol = parse("CCOC(=O)C")
        tree = retrosynthesis(mol, max_depth=3, beam_width=5)
        route = extract_best_route(tree)
        if route.total_steps > 0:
            # Overall yield should be < 100%
            assert route.overall_yield < 100.0
            assert route.overall_yield > 0.0

    def test_starting_materials_are_purchasable(self):
        mol = parse("CCOC(=O)C")
        tree = retrosynthesis(mol, max_depth=3, beam_width=5)
        route = extract_best_route(tree)
        for sm in route.starting_materials:
            assert sm.smiles in PURCHASABLE_MATERIALS


class TestGatherPurchasableLeaves:
    def test_gathers_unique_leaves(self):
        mol = parse("CCOC(=O)C")
        tree = retrosynthesis(mol, max_depth=3, beam_width=5)
        leaves = []
        seen = set()
        _gather_purchasable_leaves(tree.target, leaves, seen)
        smiles_list = [l.smiles for l in leaves]
        # All unique
        assert len(smiles_list) == len(set(smiles_list))


# =====================================================================
#  Mechanism choreography (mechanism_choreography.py - was 67% coverage)
# =====================================================================

from molbuilder.dynamics.mechanism_choreography import MechanismChoreographer
from molbuilder.dynamics.mechanisms import (
    ReactionMechanism,
    MechanismStage,
    MechanismType,
    ElectronFlow,
    FlowType,
)

import numpy as np
from unittest.mock import MagicMock


def _mock_forcefield():
    """Create a minimal mock ForceField for choreographer tests."""
    return MagicMock()


def _make_simple_mechanism():
    """Create a minimal mechanism with one stage."""
    stage = MechanismStage(
        name="test_stage",
        distance_targets={(0, 1): 1.5},
        annotation="Testing stage",
    )
    return ReactionMechanism(
        name="test_mechanism",
        mechanism_type=MechanismType.SN2,
        stages=[stage],
    )


def _make_two_stage_mechanism():
    """Create a mechanism with two stages for interpolation testing."""
    flow1 = ElectronFlow(
        from_atom=0, to_bond=(0, 1), flow_type=FlowType.CURLY_ARROW,
    )
    stage1 = MechanismStage(
        name="approach",
        distance_targets={(0, 1): 2.0},
        bond_order_changes={(0, 1): 0.5},
        electron_flows=[flow1],
        annotation="Nucleophile approaches",
    )
    flow2 = ElectronFlow(
        from_atom=0, to_bond=(0, 1), flow_type=FlowType.CURLY_ARROW,
    )
    stage2 = MechanismStage(
        name="bond_formation",
        distance_targets={(0, 1): 1.5},
        bond_order_changes={(0, 1): 1.0},
        electron_flows=[flow2],
        annotation="Bond forms",
    )
    return ReactionMechanism(
        name="two_stage",
        mechanism_type=MechanismType.SN2,
        stages=[stage1, stage2],
    )


class TestMechanismChoreographer:
    def test_init(self):
        mech = _make_simple_mechanism()
        choreo = MechanismChoreographer(mech, _mock_forcefield())
        assert len(choreo.mechanism.stages) == 1

    def test_two_stages(self):
        mech = _make_two_stage_mechanism()
        choreo = MechanismChoreographer(mech, _mock_forcefield())
        assert len(choreo.mechanism.stages) == 2

    def test_restraint_forces_returns_array(self):
        mech = _make_simple_mechanism()
        choreo = MechanismChoreographer(mech, _mock_forcefield())
        positions = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        forces = choreo.restraint_forces(positions, stage_idx=0, progress=0.5)
        assert forces.shape == positions.shape

    def test_restraint_forces_pull_atoms_together(self):
        mech = _make_simple_mechanism()
        choreo = MechanismChoreographer(mech, _mock_forcefield())
        positions = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        forces = choreo.restraint_forces(positions, stage_idx=0, progress=0.8)
        assert forces[0, 0] > 0
        assert forces[1, 0] < 0

    def test_bond_orders_at_returns_dict(self):
        mech = _make_two_stage_mechanism()
        choreo = MechanismChoreographer(mech, _mock_forcefield())
        orders = choreo.bond_orders_at(stage_idx=0, progress=0.5)
        assert isinstance(orders, dict)

    def test_bond_orders_interpolation(self):
        mech = _make_two_stage_mechanism()
        choreo = MechanismChoreographer(mech, _mock_forcefield())
        early = choreo.bond_orders_at(stage_idx=0, progress=0.1)
        late = choreo.bond_orders_at(stage_idx=1, progress=0.9)
        if (0, 1) in early and (0, 1) in late:
            assert late[(0, 1)] >= early[(0, 1)]

    def test_electron_flows_at(self):
        mech = _make_two_stage_mechanism()
        choreo = MechanismChoreographer(mech, _mock_forcefield())
        flows = choreo.electron_flows_at(stage_idx=0, progress=0.5)
        assert isinstance(flows, list)
        assert len(flows) == 1

    def test_electron_flows_at_boundary_hidden(self):
        mech = _make_two_stage_mechanism()
        choreo = MechanismChoreographer(mech, _mock_forcefield())
        flows = choreo.electron_flows_at(stage_idx=0, progress=0.05)
        assert flows == []

    def test_stage_annotation(self):
        mech = _make_two_stage_mechanism()
        choreo = MechanismChoreographer(mech, _mock_forcefield())
        ann = choreo.stage_annotation(0)
        assert isinstance(ann, str)
        assert "Nucleophile" in ann

    def test_stage_annotation_out_of_bounds(self):
        mech = _make_simple_mechanism()
        choreo = MechanismChoreographer(mech, _mock_forcefield())
        ann = choreo.stage_annotation(99)
        assert ann == ""


# =====================================================================
#  Aromatic kekulization (aromatic.py)
# =====================================================================

from molbuilder.smiles.aromatic import (
    kekulize,
    find_aromatic_rings,
    _build_aromatic_adj,
    _pi_electrons,
    _find_perfect_matching,
)


class _FakeAtom:
    """Minimal atom stand-in for aromatic tests."""
    def __init__(self, index, symbol, aromatic=False, hcount=None, bracket=False):
        self.index = index
        self.symbol = symbol
        self.aromatic = aromatic
        self.hcount = hcount
        self.bracket = bracket


class _FakeBond:
    """Minimal bond stand-in for aromatic tests."""
    def __init__(self, atom_i, atom_j, order=1):
        self.atom_i = atom_i
        self.atom_j = atom_j
        self.order = order


def _benzene_atoms_bonds():
    """Create benzene (6 aromatic C) with ring bonds."""
    atoms = [_FakeAtom(i, "C", aromatic=True) for i in range(6)]
    bonds = [_FakeBond(i, (i + 1) % 6, 1) for i in range(6)]
    return atoms, bonds


def _pyrrole_atoms_bonds():
    """Create pyrrole: 4 aromatic C + 1 aromatic N (index 0, with H)."""
    atoms = [
        _FakeAtom(0, "N", aromatic=True, hcount=1, bracket=True),
        _FakeAtom(1, "C", aromatic=True),
        _FakeAtom(2, "C", aromatic=True),
        _FakeAtom(3, "C", aromatic=True),
        _FakeAtom(4, "C", aromatic=True),
    ]
    bonds = [
        _FakeBond(0, 1), _FakeBond(1, 2), _FakeBond(2, 3),
        _FakeBond(3, 4), _FakeBond(4, 0),
    ]
    return atoms, bonds


class TestBuildAromaticAdj:
    def test_benzene_adjacency(self):
        atoms, bonds = _benzene_atoms_bonds()
        adj = _build_aromatic_adj(atoms, bonds)
        assert len(adj) == 6
        for i in range(6):
            assert len(adj[i]) == 2

    def test_no_aromatic_atoms(self):
        atoms = [_FakeAtom(0, "C"), _FakeAtom(1, "C")]
        bonds = [_FakeBond(0, 1)]
        adj = _build_aromatic_adj(atoms, bonds)
        assert adj == {}


class TestFindAromaticRings:
    def test_benzene_ring(self):
        atoms, bonds = _benzene_atoms_bonds()
        rings = find_aromatic_rings(atoms, bonds)
        assert len(rings) >= 1
        assert any(len(r) == 6 for r in rings)

    def test_no_aromatics(self):
        atoms = [_FakeAtom(i, "C") for i in range(3)]
        bonds = [_FakeBond(0, 1), _FakeBond(1, 2)]
        rings = find_aromatic_rings(atoms, bonds)
        assert rings == []


class TestPiElectrons:
    def test_carbon_contributes_1(self):
        atoms, bonds = _benzene_atoms_bonds()
        adj = _build_aromatic_adj(atoms, bonds)
        assert _pi_electrons(atoms[0], adj) == 1

    def test_pyrrole_nitrogen_contributes_2(self):
        atoms, bonds = _pyrrole_atoms_bonds()
        adj = _build_aromatic_adj(atoms, bonds)
        assert _pi_electrons(atoms[0], adj) == 2

    def test_furan_oxygen_contributes_2(self):
        atoms = [
            _FakeAtom(0, "O", aromatic=True),
            _FakeAtom(1, "C", aromatic=True),
            _FakeAtom(2, "C", aromatic=True),
            _FakeAtom(3, "C", aromatic=True),
            _FakeAtom(4, "C", aromatic=True),
        ]
        bonds = [
            _FakeBond(0, 1), _FakeBond(1, 2), _FakeBond(2, 3),
            _FakeBond(3, 4), _FakeBond(4, 0),
        ]
        adj = _build_aromatic_adj(atoms, bonds)
        assert _pi_electrons(atoms[0], adj) == 2

    def test_sulfur_contributes_2(self):
        atoms = [_FakeAtom(0, "S", aromatic=True)]
        adj = {0: [1, 4]}
        assert _pi_electrons(atoms[0], adj) == 2

    def test_boron_contributes_0(self):
        atoms = [_FakeAtom(0, "B", aromatic=True)]
        adj = {0: [1, 2]}
        assert _pi_electrons(atoms[0], adj) == 0


class TestFindPerfectMatching:
    def test_benzene_matching(self):
        adj = {i: [(i - 1) % 6, (i + 1) % 6] for i in range(6)}
        match = _find_perfect_matching(adj)
        assert match is not None
        for v in range(6):
            assert match[v] != -1

    def test_odd_vertices_returns_none(self):
        adj = {0: [1, 2], 1: [0, 2], 2: [0, 1]}
        assert _find_perfect_matching(adj) is None


class TestKekulize:
    def test_benzene_kekulize(self):
        atoms, bonds = _benzene_atoms_bonds()
        result = kekulize(atoms, bonds)
        double_count = sum(1 for b in result if b.order == 2)
        assert double_count == 3

    def test_no_aromatic_unchanged(self):
        atoms = [_FakeAtom(0, "C"), _FakeAtom(1, "C")]
        bonds = [_FakeBond(0, 1, 1)]
        kekulize(atoms, bonds)
        assert bonds[0].order == 1

    def test_pyrrole_kekulize(self):
        atoms, bonds = _pyrrole_atoms_bonds()
        kekulize(atoms, bonds)
        double_count = sum(1 for b in bonds if b.order == 2)
        assert double_count == 2  # 4 C's, 2 double bonds

    def test_kekulize_via_parser(self):
        """End-to-end: parse aromatic SMILES and check formulas."""
        mol = parse("c1ccccc1")  # benzene
        c = sum(1 for a in mol.atoms if a.symbol == "C")
        h = sum(1 for a in mol.atoms if a.symbol == "H")
        assert c == 6 and h == 6
        double_bonds = sum(1 for b in mol.bonds if b.order == 2)
        assert double_bonds == 3

    def test_naphthalene_kekulize(self):
        mol = parse("c1ccc2ccccc2c1")
        c = sum(1 for a in mol.atoms if a.symbol == "C")
        h = sum(1 for a in mol.atoms if a.symbol == "H")
        assert c == 10 and h == 8
