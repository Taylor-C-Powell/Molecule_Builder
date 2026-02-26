"""Tests for the RetroCast adapter module.

Verifies conversion of MolBuilder RetrosynthesisTree structures into the
RetroCast canonical Route schema (dict-based, JSON-serializable).
"""

from __future__ import annotations

import gzip
import json
import os
import tempfile
from unittest.mock import MagicMock

import pytest

from molbuilder.reactions.retrosynthesis import (
    RetrosynthesisTree,
    RetroNode,
    Disconnection,
    Precursor,
)
from molbuilder.reactions.retrocast_adapter import (
    tree_to_retrocast_routes,
    export_retrocast_json,
)
import molbuilder


# =====================================================================
#  Test helpers
# =====================================================================

def _mock_template(name="Test Reaction", category_value="coupling", score=75.0):
    """Create a mock ReactionTemplate with the attributes the adapter reads."""
    t = MagicMock()
    t.name = name
    t.named_reaction = name
    t.category.value = category_value
    return t


def _mock_molecule(name="mock"):
    """Create a mock Molecule object (not accessed by the adapter)."""
    m = MagicMock()
    m.name = name
    return m


def _make_leaf_node(smiles, depth=1):
    """Build a purchasable leaf RetroNode."""
    return RetroNode(
        smiles=smiles,
        molecule=_mock_molecule(smiles),
        is_purchasable=True,
        depth=depth,
    )


def _make_disconnection(template_name, precursor_smiles_list, score):
    """Build a Disconnection with mock template and precursors."""
    tmpl = _mock_template(name=template_name, score=score)
    precursors = [
        Precursor(smiles=smi, molecule=None, name=f"precursor-{smi}",
                  cost_per_kg=10.0)
        for smi in precursor_smiles_list
    ]
    return Disconnection(template=tmpl, precursors=precursors, score=score)


# =====================================================================
#  1. Empty tree (purchasable target) produces empty routes list
# =====================================================================

class TestEmptyTree:
    def test_purchasable_target_returns_no_routes(self):
        root = RetroNode(
            smiles="CCO",
            molecule=_mock_molecule("ethanol"),
            is_purchasable=True,
            depth=0,
        )
        tree = RetrosynthesisTree(
            target=root, max_depth=3, beam_width=5, routes_found=0,
        )
        routes = tree_to_retrocast_routes(tree)
        assert routes == []

    def test_no_disconnections_returns_no_routes(self):
        root = RetroNode(
            smiles="C1=CC=CC=C1",
            molecule=_mock_molecule("benzene"),
            is_purchasable=False,
            depth=0,
        )
        tree = RetrosynthesisTree(
            target=root, max_depth=3, beam_width=5, routes_found=0,
        )
        routes = tree_to_retrocast_routes(tree)
        assert routes == []


# =====================================================================
#  2. 1-step route conversion
# =====================================================================

class TestOneStepRoute:
    def test_single_disconnection_to_purchasable_precursors(self):
        disc = _make_disconnection("Grignard addition", ["CBr", "CC=O"], 80.0)

        child_a = _make_leaf_node("CBr", depth=1)
        child_b = _make_leaf_node("CC=O", depth=1)

        root = RetroNode(
            smiles="CCC(O)=O",
            molecule=_mock_molecule("target"),
            is_purchasable=False,
            disconnections=[disc],
            best_disconnection=disc,
            children=[child_a, child_b],
            depth=0,
        )
        tree = RetrosynthesisTree(
            target=root, max_depth=3, beam_width=5, routes_found=1,
        )
        routes = tree_to_retrocast_routes(tree)

        assert len(routes) == 1
        route = routes[0]
        assert len(route["steps"]) == 1
        step = route["steps"][0]
        assert step["index"] == 0
        assert step["reaction_name"] == "Grignard addition"
        assert step["product"]["smiles"] == "CCC(O)=O"
        assert len(step["reactants"]) == 2


# =====================================================================
#  3. Multi-step route (depth 2)
# =====================================================================

class TestMultiStepRoute:
    def test_depth_two_route(self):
        # Inner disconnection: intermediate -> leaf precursors
        inner_disc = _make_disconnection("Reduction", ["CC=O"], 60.0)
        inner_leaf = _make_leaf_node("CC=O", depth=2)

        intermediate = RetroNode(
            smiles="CCO",
            molecule=_mock_molecule("ethanol"),
            is_purchasable=False,
            disconnections=[inner_disc],
            best_disconnection=inner_disc,
            children=[inner_leaf],
            depth=1,
        )

        # Root disconnection: target -> intermediate + purchasable
        root_disc = _make_disconnection("Esterification", ["CCO", "CC(O)=O"], 75.0)
        acetic = _make_leaf_node("CC(O)=O", depth=1)

        root = RetroNode(
            smiles="CCOC(C)=O",
            molecule=_mock_molecule("ethyl acetate"),
            is_purchasable=False,
            disconnections=[root_disc],
            best_disconnection=root_disc,
            children=[intermediate, acetic],
            depth=0,
        )
        tree = RetrosynthesisTree(
            target=root, max_depth=4, beam_width=5, routes_found=1,
        )
        routes = tree_to_retrocast_routes(tree)

        assert len(routes) == 1
        route = routes[0]
        # 2 steps: inner reduction then root esterification
        assert len(route["steps"]) == 2
        # Forward order: first step produces intermediate from leaf
        assert route["steps"][0]["reaction_name"] == "Reduction"
        assert route["steps"][1]["reaction_name"] == "Esterification"


# =====================================================================
#  4. Leaf nodes have synthesis_step=None
# =====================================================================

class TestLeafSynthesisStep:
    def test_leaf_reactants_have_none_synthesis_step(self):
        disc = _make_disconnection("Coupling", ["C=C", "CBr"], 70.0)
        child_a = _make_leaf_node("C=C", depth=1)
        child_b = _make_leaf_node("CBr", depth=1)

        root = RetroNode(
            smiles="CC=C",
            molecule=_mock_molecule("propene"),
            is_purchasable=False,
            disconnections=[disc],
            best_disconnection=disc,
            children=[child_a, child_b],
            depth=0,
        )
        tree = RetrosynthesisTree(
            target=root, max_depth=3, beam_width=5, routes_found=1,
        )
        routes = tree_to_retrocast_routes(tree)
        step = routes[0]["steps"][0]
        for reactant in step["reactants"]:
            assert reactant["synthesis_step"] is None, (
                "Leaf (purchasable) reactants must have synthesis_step=None"
            )


# =====================================================================
#  5. Rank ordering matches score ordering
# =====================================================================

class TestRankOrdering:
    def test_ranks_follow_descending_score(self):
        disc_low = _make_disconnection("Low Score Rxn", ["C"], 30.0)
        disc_mid = _make_disconnection("Mid Score Rxn", ["CC"], 55.0)
        disc_high = _make_disconnection("High Score Rxn", ["CCC"], 90.0)

        # best_disconnection is disc_high (first in sorted order)
        root = RetroNode(
            smiles="CCCC",
            molecule=_mock_molecule("butane"),
            is_purchasable=False,
            disconnections=[disc_low, disc_mid, disc_high],
            best_disconnection=disc_high,
            children=[_make_leaf_node("CCC", depth=1)],
            depth=0,
        )
        tree = RetrosynthesisTree(
            target=root, max_depth=3, beam_width=5, routes_found=1,
        )
        routes = tree_to_retrocast_routes(tree)

        assert len(routes) == 3
        assert routes[0]["rank"] == 1
        assert routes[0]["metadata"]["score"] == 90.0
        assert routes[1]["rank"] == 2
        assert routes[1]["metadata"]["score"] == 55.0
        assert routes[2]["rank"] == 3
        assert routes[2]["metadata"]["score"] == 30.0


# =====================================================================
#  6. Metadata includes source / version / score
# =====================================================================

class TestMetadata:
    def test_metadata_fields_present(self):
        disc = _make_disconnection("Some Reaction", ["C"], 65.0)
        root = RetroNode(
            smiles="CC",
            molecule=_mock_molecule("ethane"),
            is_purchasable=False,
            disconnections=[disc],
            best_disconnection=disc,
            children=[_make_leaf_node("C", depth=1)],
            depth=0,
        )
        tree = RetrosynthesisTree(
            target=root, max_depth=3, beam_width=5, routes_found=1,
        )
        routes = tree_to_retrocast_routes(tree)
        meta = routes[0]["metadata"]
        assert meta["source"] == "MolBuilder"
        assert meta["version"] == molbuilder.__version__
        assert meta["score"] == 65.0


# =====================================================================
#  7. JSON export round-trip
# =====================================================================

class TestJsonExport:
    def test_gzip_round_trip(self, tmp_path):
        disc = _make_disconnection("Aldol", ["CC=O", "CC(C)=O"], 72.0)
        root = RetroNode(
            smiles="CC(O)CC(C)=O",
            molecule=_mock_molecule("aldol product"),
            is_purchasable=False,
            disconnections=[disc],
            best_disconnection=disc,
            children=[
                _make_leaf_node("CC=O", depth=1),
                _make_leaf_node("CC(C)=O", depth=1),
            ],
            depth=0,
        )
        tree = RetrosynthesisTree(
            target=root, max_depth=3, beam_width=5, routes_found=1,
        )

        out_path = str(tmp_path / "routes.json.gz")
        export_retrocast_json(tree, out_path)

        # Read back and verify
        with gzip.open(out_path, "rt", encoding="utf-8") as fh:
            data = json.load(fh)

        assert isinstance(data, list)
        assert len(data) == 1
        route = data[0]
        assert "target" in route
        assert "steps" in route
        assert "rank" in route
        assert "metadata" in route
        assert route["metadata"]["source"] == "MolBuilder"

    def test_export_creates_valid_gzip(self, tmp_path):
        disc = _make_disconnection("Rxn", ["C"], 50.0)
        root = RetroNode(
            smiles="CC",
            molecule=_mock_molecule("ethane"),
            is_purchasable=False,
            disconnections=[disc],
            best_disconnection=disc,
            children=[_make_leaf_node("C", depth=1)],
            depth=0,
        )
        tree = RetrosynthesisTree(
            target=root, max_depth=3, beam_width=5, routes_found=1,
        )

        out_path = str(tmp_path / "test.json.gz")
        export_retrocast_json(tree, out_path)

        assert os.path.exists(out_path)
        # Should be a valid gzip file
        with gzip.open(out_path, "rb") as fh:
            raw = fh.read()
        # Must be valid JSON
        parsed = json.loads(raw)
        assert isinstance(parsed, list)


# =====================================================================
#  8. Route steps are in forward order (leaves first, target last)
# =====================================================================

class TestForwardOrder:
    def test_steps_ordered_leaves_to_target(self):
        # depth-2 tree: target -> intermediate -> leaf
        inner_disc = _make_disconnection("Step A", ["C"], 50.0)
        leaf = _make_leaf_node("C", depth=2)

        intermediate = RetroNode(
            smiles="CC",
            molecule=_mock_molecule("intermediate"),
            is_purchasable=False,
            disconnections=[inner_disc],
            best_disconnection=inner_disc,
            children=[leaf],
            depth=1,
        )

        root_disc = _make_disconnection("Step B", ["CC", "O"], 70.0)
        water = _make_leaf_node("O", depth=1)

        root = RetroNode(
            smiles="CCO",
            molecule=_mock_molecule("ethanol"),
            is_purchasable=False,
            disconnections=[root_disc],
            best_disconnection=root_disc,
            children=[intermediate, water],
            depth=0,
        )
        tree = RetrosynthesisTree(
            target=root, max_depth=4, beam_width=5, routes_found=1,
        )
        routes = tree_to_retrocast_routes(tree)
        steps = routes[0]["steps"]

        # First step should produce the intermediate (leaf -> intermediate)
        assert steps[0]["product"]["smiles"] == "CC"
        # Last step should produce the target
        assert steps[-1]["product"]["smiles"] == "CCO"


# =====================================================================
#  9. Each step has correct reactants from precursors
# =====================================================================

class TestReactantsCorrectness:
    def test_step_reactants_match_disconnection_precursors(self):
        disc = _make_disconnection(
            "Williamson Ether", ["CCBr", "CO"], 68.0,
        )
        root = RetroNode(
            smiles="CCOCC",
            molecule=_mock_molecule("diethyl ether"),
            is_purchasable=False,
            disconnections=[disc],
            best_disconnection=disc,
            children=[
                _make_leaf_node("CCBr", depth=1),
                _make_leaf_node("CO", depth=1),
            ],
            depth=0,
        )
        tree = RetrosynthesisTree(
            target=root, max_depth=3, beam_width=5, routes_found=1,
        )
        routes = tree_to_retrocast_routes(tree)
        step = routes[0]["steps"][0]
        reactant_smiles = {r["smiles"] for r in step["reactants"]}
        assert reactant_smiles == {"CCBr", "CO"}


# =====================================================================
#  10. Multiple disconnections at root produce multiple routes
# =====================================================================

class TestMultipleRootDisconnections:
    def test_each_root_disconnection_yields_a_route(self):
        disc_a = _make_disconnection("Route A Rxn", ["C", "C"], 85.0)
        disc_b = _make_disconnection("Route B Rxn", ["CC"], 60.0)
        disc_c = _make_disconnection("Route C Rxn", ["CCC"], 45.0)

        root = RetroNode(
            smiles="CCCC",
            molecule=_mock_molecule("butane"),
            is_purchasable=False,
            disconnections=[disc_a, disc_b, disc_c],
            best_disconnection=disc_a,
            children=[
                _make_leaf_node("C", depth=1),
                _make_leaf_node("C", depth=1),
            ],
            depth=0,
        )
        tree = RetrosynthesisTree(
            target=root, max_depth=3, beam_width=5, routes_found=1,
        )
        routes = tree_to_retrocast_routes(tree)
        assert len(routes) == 3

        # Each route should have at least one step
        for route in routes:
            assert len(route["steps"]) >= 1

        # Reaction names should cover all three disconnections
        rxn_names = {route["steps"][-1]["reaction_name"] for route in routes}
        assert "Route A Rxn" in rxn_names
        assert "Route B Rxn" in rxn_names
        assert "Route C Rxn" in rxn_names

    def test_alternate_disconnections_produce_single_step_routes(self):
        """Alternate (non-best) disconnections produce shallow 1-step routes."""
        disc_best = _make_disconnection("Best Rxn", ["C", "CC"], 90.0)
        disc_alt = _make_disconnection("Alt Rxn", ["CCC"], 40.0)

        root = RetroNode(
            smiles="CCCC",
            molecule=_mock_molecule("butane"),
            is_purchasable=False,
            disconnections=[disc_best, disc_alt],
            best_disconnection=disc_best,
            children=[
                _make_leaf_node("C", depth=1),
                _make_leaf_node("CC", depth=1),
            ],
            depth=0,
        )
        tree = RetrosynthesisTree(
            target=root, max_depth=3, beam_width=5, routes_found=1,
        )
        routes = tree_to_retrocast_routes(tree)
        assert len(routes) == 2

        # The alternate disconnection should produce a 1-step route
        alt_route = [r for r in routes
                     if r["steps"][-1]["reaction_name"] == "Alt Rxn"][0]
        assert len(alt_route["steps"]) == 1
