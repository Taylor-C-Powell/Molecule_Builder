"""Tests for reaction mechanism templates and choreography."""

import math
import numpy as np
import pytest

from molbuilder.molecule.graph import Molecule, Hybridization
from molbuilder.core.bond_data import bond_length


# ===================================================================
# Mechanism data model tests
# ===================================================================

class TestMechanismTypes:
    """Test mechanism enumeration and data classes."""

    def test_mechanism_type_enum(self):
        """MechanismType should have expected members."""
        from molbuilder.dynamics.mechanisms import MechanismType
        assert hasattr(MechanismType, "SN2")
        assert hasattr(MechanismType, "E2")
        assert hasattr(MechanismType, "RADICAL_SUBSTITUTION")
        assert hasattr(MechanismType, "NUCLEOPHILIC_ADDITION")

    def test_flow_type_enum(self):
        """FlowType should have CURLY_ARROW and FISHHOOK_ARROW."""
        from molbuilder.dynamics.mechanisms import FlowType
        assert hasattr(FlowType, "CURLY_ARROW")
        assert hasattr(FlowType, "FISHHOOK_ARROW")

    def test_electron_flow_creation(self):
        """ElectronFlow should be creatable with all fields."""
        from molbuilder.dynamics.mechanisms import ElectronFlow, FlowType
        flow = ElectronFlow(
            from_atom=0,
            to_bond=(0, 1),
            flow_type=FlowType.CURLY_ARROW,
            label="test")
        assert flow.from_atom == 0
        assert flow.to_bond == (0, 1)
        assert flow.label == "test"

    def test_mechanism_stage_defaults(self):
        """MechanismStage should have sensible defaults."""
        from molbuilder.dynamics.mechanisms import MechanismStage
        stage = MechanismStage(name="test")
        assert stage.name == "test"
        assert stage.distance_targets == {}
        assert stage.bond_order_changes == {}
        assert stage.electron_flows == []
        assert stage.duration_weight == 1.0


# ===================================================================
# Predefined template tests
# ===================================================================

class TestSN2Mechanism:
    """Test the SN2 mechanism template."""

    def test_sn2_has_three_stages(self):
        """SN2 mechanism should have 3 stages."""
        from molbuilder.dynamics.mechanisms import sn2_mechanism, MechanismType
        mech = sn2_mechanism(
            substrate_C=0, nucleophile=1, leaving_group=2)
        assert len(mech.stages) == 3
        assert mech.mechanism_type == MechanismType.SN2

    def test_sn2_atom_roles(self):
        """SN2 mechanism should record atom roles."""
        from molbuilder.dynamics.mechanisms import sn2_mechanism
        mech = sn2_mechanism(
            substrate_C=0, nucleophile=5, leaving_group=3)
        assert mech.atom_roles["substrate_C"] == 0
        assert mech.atom_roles["nucleophile"] == 5
        assert mech.atom_roles["leaving_group"] == 3

    def test_sn2_stage_names(self):
        """SN2 stages should have descriptive names."""
        from molbuilder.dynamics.mechanisms import sn2_mechanism
        mech = sn2_mechanism(
            substrate_C=0, nucleophile=1, leaving_group=2)
        names = [s.name for s in mech.stages]
        assert "approach" in names[0].lower()
        assert "transition" in names[1].lower()
        assert "product" in names[2].lower()

    def test_sn2_electron_flows_in_ts(self):
        """Transition state should have electron flow arrows."""
        from molbuilder.dynamics.mechanisms import sn2_mechanism
        mech = sn2_mechanism(
            substrate_C=0, nucleophile=1, leaving_group=2)
        ts_stage = mech.stages[1]
        assert len(ts_stage.electron_flows) >= 2

    def test_sn2_bond_orders_progress(self):
        """Bond orders should progress: Nu bond grows, LG bond shrinks."""
        from molbuilder.dynamics.mechanisms import sn2_mechanism
        mech = sn2_mechanism(
            substrate_C=0, nucleophile=1, leaving_group=2)
        nuc_key = (0, 1)
        lg_key = (0, 2)

        stage0_nuc = mech.stages[0].bond_order_changes.get(nuc_key, 0)
        stage2_nuc = mech.stages[2].bond_order_changes.get(nuc_key, 0)
        stage0_lg = mech.stages[0].bond_order_changes.get(lg_key, 1)
        stage2_lg = mech.stages[2].bond_order_changes.get(lg_key, 1)

        assert stage2_nuc > stage0_nuc  # Nu bond strengthens
        assert stage2_lg < stage0_lg    # LG bond weakens


class TestE2Mechanism:
    """Test the E2 elimination template."""

    def test_e2_has_three_stages(self):
        """E2 mechanism should have 3 stages."""
        from molbuilder.dynamics.mechanisms import e2_mechanism, MechanismType
        mech = e2_mechanism(
            alpha_C=0, beta_H=1, base=2, leaving_group=3)
        assert len(mech.stages) == 3
        assert mech.mechanism_type == MechanismType.E2

    def test_e2_atom_roles(self):
        """E2 mechanism should record all atom roles."""
        from molbuilder.dynamics.mechanisms import e2_mechanism
        mech = e2_mechanism(
            alpha_C=0, beta_H=1, base=2, leaving_group=3)
        assert "alpha_C" in mech.atom_roles
        assert "beta_H" in mech.atom_roles
        assert "base" in mech.atom_roles
        assert "leaving_group" in mech.atom_roles


class TestRadicalMechanism:
    """Test the radical substitution template."""

    def test_radical_has_three_stages(self):
        """Radical mechanism should have 3 stages."""
        from molbuilder.dynamics.mechanisms import (
            radical_substitution_mechanism, MechanismType,
        )
        mech = radical_substitution_mechanism(
            target_C=0, target_H=1, radical=2)
        assert len(mech.stages) == 3
        assert mech.mechanism_type == MechanismType.RADICAL_SUBSTITUTION

    def test_radical_uses_fishhook_arrows(self):
        """Radical mechanism should use fishhook arrows."""
        from molbuilder.dynamics.mechanisms import (
            radical_substitution_mechanism, FlowType,
        )
        mech = radical_substitution_mechanism(
            target_C=0, target_H=1, radical=2)
        all_flows = []
        for stage in mech.stages:
            all_flows.extend(stage.electron_flows)
        fishhooks = [f for f in all_flows
                     if f.flow_type == FlowType.FISHHOOK_ARROW]
        assert len(fishhooks) >= 2


class TestNucleophilicAddition:
    """Test the nucleophilic addition template."""

    def test_nuc_add_has_three_stages(self):
        """Nucleophilic addition should have 3 stages."""
        from molbuilder.dynamics.mechanisms import (
            nucleophilic_addition_mechanism, MechanismType,
        )
        mech = nucleophilic_addition_mechanism(
            carbonyl_C=0, carbonyl_O=1, nucleophile=2)
        assert len(mech.stages) == 3
        assert mech.mechanism_type == MechanismType.NUCLEOPHILIC_ADDITION


# ===================================================================
# Choreography tests
# ===================================================================

class TestMechanismChoreographer:
    """Test the steered-MD choreographer."""

    def _build_test_system(self):
        """Build a simple 3-atom test system."""
        from molbuilder.dynamics.forcefield import ForceField
        from molbuilder.dynamics.mechanisms import sn2_mechanism
        from molbuilder.dynamics.mechanism_choreography import MechanismChoreographer

        mol = Molecule("test")
        mol.add_atom("C", [0.0, 0.0, 0.0], Hybridization.SP3)
        mol.add_atom("O", [3.0, 0.0, 0.0])
        mol.add_atom("Cl", [-2.0, 0.0, 0.0])
        mol.add_bond(0, 2, order=1)

        ff = ForceField.from_molecule(mol)
        mech = sn2_mechanism(substrate_C=0, nucleophile=1, leaving_group=2)
        choreo = MechanismChoreographer(mech, ff, n_steps_per_stage=100)

        return mol, ff, mech, choreo

    def test_restraint_forces_shape(self):
        """Restraint forces should have correct shape."""
        mol, ff, mech, choreo = self._build_test_system()
        pos = np.array([a.position for a in mol.atoms])
        forces = choreo.restraint_forces(pos, stage_idx=0, progress=0.5)
        assert forces.shape == (len(mol.atoms), 3)

    def test_restraint_forces_finite(self):
        """Restraint forces should be finite."""
        mol, ff, mech, choreo = self._build_test_system()
        pos = np.array([a.position for a in mol.atoms])
        forces = choreo.restraint_forces(pos, stage_idx=0, progress=0.5)
        assert np.all(np.isfinite(forces))

    def test_restraint_forces_increase_with_progress(self):
        """Restraint forces should increase as progress increases (sigmoid)."""
        mol, ff, mech, choreo = self._build_test_system()
        pos = np.array([a.position for a in mol.atoms])
        f_early = choreo.restraint_forces(pos, stage_idx=0, progress=0.1)
        f_late = choreo.restraint_forces(pos, stage_idx=0, progress=0.9)
        # Late forces should generally be larger
        assert np.max(np.abs(f_late)) >= np.max(np.abs(f_early)) * 0.5

    def test_bond_orders_at_progress(self):
        """Bond orders should be computed at any progress."""
        mol, ff, mech, choreo = self._build_test_system()
        bo = choreo.bond_orders_at(stage_idx=0, progress=0.5)
        assert isinstance(bo, dict)

    def test_electron_flows_at_midstage(self):
        """Electron flows should be returned at mid-stage."""
        mol, ff, mech, choreo = self._build_test_system()
        flows = choreo.electron_flows_at(stage_idx=0, progress=0.5)
        assert isinstance(flows, list)

    def test_electron_flows_empty_at_edges(self):
        """Electron flows should be empty near stage boundaries."""
        mol, ff, mech, choreo = self._build_test_system()
        flows_start = choreo.electron_flows_at(stage_idx=0, progress=0.0)
        flows_end = choreo.electron_flows_at(stage_idx=0, progress=1.0)
        assert flows_start == []
        assert flows_end == []

    def test_stage_annotation(self):
        """Stage annotation should return a non-empty string."""
        mol, ff, mech, choreo = self._build_test_system()
        ann = choreo.stage_annotation(0)
        assert isinstance(ann, str)
        assert len(ann) > 0

    def test_stage_annotation_out_of_range(self):
        """Out-of-range stage index should return empty string."""
        mol, ff, mech, choreo = self._build_test_system()
        ann = choreo.stage_annotation(99)
        assert ann == ""
