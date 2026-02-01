"""Unit tests for molbuilder reaction modules."""

import unittest

from molbuilder.smiles import parse

from molbuilder.reactions.reaction_types import ReactionCategory, ReactionTemplate
from molbuilder.reactions.reagent_data import REAGENT_DB, SOLVENT_DB
from molbuilder.reactions.knowledge_base import (
    REACTION_TEMPLATES,
    lookup_by_category,
    lookup_by_name,
)
from molbuilder.reactions.functional_group_detect import (
    detect_functional_groups,
    FunctionalGroup,
)
from molbuilder.reactions.retrosynthesis import (
    PURCHASABLE_MATERIALS,
    is_purchasable,
    retrosynthesis,
    RetrosynthesisTree,
)
from molbuilder.reactions.synthesis_route import (
    extract_best_route,
    SynthesisRoute,
)


# ===================================================================
#  Reaction types
# ===================================================================

class TestReactionCategory(unittest.TestCase):
    """Tests for ReactionCategory enum."""

    def test_at_least_10_members(self):
        """ReactionCategory should have at least 10 members."""
        members = list(ReactionCategory)
        self.assertGreaterEqual(len(members), 10)


class TestReactionTemplate(unittest.TestCase):
    """Tests for ReactionTemplate dataclass fields and helpers."""

    def _make_template(self):
        return ReactionTemplate(
            name="Test reaction",
            named_reaction="Test",
            category=ReactionCategory.SUBSTITUTION,
            reagents=["NaOH"],
            solvents=["DMSO"],
            functional_group_required=["alkyl_halide"],
            functional_group_produced=["alcohol"],
            functional_group_incompatible=["ester"],
        )

    def test_fields_exist(self):
        """ReactionTemplate should have name, category, reagents, etc."""
        t = self._make_template()
        self.assertEqual(t.name, "Test reaction")
        self.assertEqual(t.category, ReactionCategory.SUBSTITUTION)
        self.assertIsInstance(t.reagents, list)
        self.assertIsInstance(t.solvents, list)

    def test_requires(self):
        """requires() should return True for listed functional groups."""
        t = self._make_template()
        self.assertTrue(t.requires("alkyl_halide"))
        self.assertFalse(t.requires("alcohol"))

    def test_produces(self):
        """produces() should return True for listed functional groups."""
        t = self._make_template()
        self.assertTrue(t.produces("alcohol"))
        self.assertFalse(t.produces("ketone"))

    def test_is_compatible(self):
        """is_compatible() should return False if an incompatible FG is present."""
        t = self._make_template()
        self.assertTrue(t.is_compatible(["alkyl_halide", "alcohol"]))
        self.assertFalse(t.is_compatible(["ester"]))


# ===================================================================
#  Reagent and solvent databases
# ===================================================================

class TestReagentData(unittest.TestCase):
    """Tests for reagent and solvent databases."""

    def test_reagent_db_size(self):
        """REAGENT_DB should have at least 100 entries."""
        self.assertGreaterEqual(len(REAGENT_DB), 100)

    def test_solvent_db_size(self):
        """SOLVENT_DB should have at least 30 entries."""
        self.assertGreaterEqual(len(SOLVENT_DB), 30)


# ===================================================================
#  Knowledge base
# ===================================================================

class TestKnowledgeBase(unittest.TestCase):
    """Tests for the reaction knowledge base."""

    def test_template_count(self):
        """REACTION_TEMPLATES should have at least 80 entries."""
        self.assertGreaterEqual(len(REACTION_TEMPLATES), 80)

    def test_lookup_by_category_oxidation(self):
        """lookup_by_category(OXIDATION) should return non-empty list."""
        results = lookup_by_category(ReactionCategory.OXIDATION)
        self.assertGreater(len(results), 0)

    def test_lookup_by_category_reduction(self):
        """lookup_by_category(REDUCTION) should return non-empty list."""
        results = lookup_by_category(ReactionCategory.REDUCTION)
        self.assertGreater(len(results), 0)

    def test_lookup_by_category_coupling(self):
        """lookup_by_category(COUPLING) should return non-empty list."""
        results = lookup_by_category(ReactionCategory.COUPLING)
        self.assertGreater(len(results), 0)

    def test_lookup_by_name_grignard(self):
        """lookup_by_name('grignard') should return at least 2 results."""
        results = lookup_by_name("grignard")
        self.assertGreaterEqual(len(results), 2)


# ===================================================================
#  Functional group detection
# ===================================================================

class TestFunctionalGroupDetection(unittest.TestCase):
    """Tests for detect_functional_groups."""

    def _detect_names(self, smiles):
        """Helper: parse SMILES and return set of detected FG names."""
        mol = parse(smiles)
        fgs = detect_functional_groups(mol)
        return {fg.name for fg in fgs}

    def test_detect_alcohol(self):
        """CCO should be detected as containing an alcohol."""
        names = self._detect_names("CCO")
        self.assertIn("alcohol", names)

    def test_detect_aldehyde(self):
        """CC=O (acetaldehyde) should be detected as containing an aldehyde."""
        names = self._detect_names("CC=O")
        self.assertIn("aldehyde", names)

    def test_detect_carboxylic_acid(self):
        """CC(=O)O should be detected as containing a carboxylic acid."""
        names = self._detect_names("CC(=O)O")
        self.assertIn("carboxylic_acid", names)

    def test_detect_primary_amine(self):
        """CCN should be detected as containing a primary amine."""
        names = self._detect_names("CCN")
        self.assertIn("primary_amine", names)

    def test_detect_alkene(self):
        """C=C should be detected as containing an alkene."""
        names = self._detect_names("C=C")
        self.assertIn("alkene", names)

    def test_detect_alkyne(self):
        """C#C should be detected as containing an alkyne."""
        names = self._detect_names("C#C")
        self.assertIn("alkyne", names)


class TestFunctionalGroupDataclass(unittest.TestCase):
    """Tests for FunctionalGroup dataclass."""

    def test_fg_has_name_atoms_center(self):
        """FunctionalGroup should have name, atoms, and center fields."""
        fg = FunctionalGroup(name="alcohol", smarts_like="[C]-[OH]",
                             atoms=[0, 1, 2], center=1)
        self.assertEqual(fg.name, "alcohol")
        self.assertEqual(fg.atoms, [0, 1, 2])
        self.assertEqual(fg.center, 1)


# ===================================================================
#  Retrosynthesis
# ===================================================================

class TestRetrosynthesis(unittest.TestCase):
    """Tests for the retrosynthesis module."""

    def test_purchasable_materials_size(self):
        """PURCHASABLE_MATERIALS should have at least 200 entries."""
        self.assertGreaterEqual(len(PURCHASABLE_MATERIALS), 200)

    def test_is_purchasable_ethanol(self):
        """is_purchasable('CCO') should be True (ethanol is available)."""
        self.assertTrue(is_purchasable("CCO"))

    def test_retrosynthesis_returns_tree(self):
        """retrosynthesis on a simple molecule should return a tree."""
        mol = parse("CCCO")
        tree = retrosynthesis(mol, max_depth=2, beam_width=3)
        self.assertIsInstance(tree, RetrosynthesisTree)
        self.assertIsNotNone(tree.target)
        self.assertIsNotNone(tree.target.smiles)


# ===================================================================
#  Synthesis route
# ===================================================================

class TestSynthesisRoute(unittest.TestCase):
    """Tests for the synthesis route module."""

    def test_extract_best_route_from_tree_with_disconnection(self):
        """extract_best_route should produce a SynthesisRoute when the tree
        has a best_disconnection on its target node."""
        # Use a non-purchasable molecule that will generate disconnections
        mol = parse("CCCO")
        tree = retrosynthesis(mol, max_depth=3, beam_width=5)
        if tree.target.best_disconnection is not None:
            route = extract_best_route(tree)
            self.assertIsInstance(route, SynthesisRoute)
            self.assertTrue(len(route.target_smiles) > 0)
        else:
            # If no disconnection was found, the route should still work
            route = extract_best_route(tree)
            self.assertIsInstance(route, SynthesisRoute)


if __name__ == "__main__":
    unittest.main()
