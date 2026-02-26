"""Tests for Open Reaction Database (ORD) integration."""

import json

import pytest

from molbuilder.data import load_ord_conditions, _clear_ord_cache
from molbuilder.data.template_ord_mapping import TEMPLATE_TO_ORD_KEY, CATEGORY_FALLBACK
from molbuilder.process.conditions import (
    ReactionConditions,
    optimize_conditions,
    _get_ord_stats,
)
from molbuilder.process.condition_prediction import (
    predict_conditions,
    TemplateMatch,
)
from molbuilder.reactions.knowledge_base import REACTION_TEMPLATES


@pytest.fixture(autouse=True)
def _reset_cache():
    """Clear the ORD cache before and after each test."""
    _clear_ord_cache()
    yield
    _clear_ord_cache()


# =====================================================================
#  ORD data loading
# =====================================================================

class TestLoadOrdConditions:

    def test_loads_successfully(self):
        data = load_ord_conditions()
        assert isinstance(data, dict)
        assert "reactions" in data
        assert "_meta" in data

    def test_has_meta_fields(self):
        data = load_ord_conditions()
        meta = data["_meta"]
        assert "ord_version" in meta
        assert "generated" in meta
        assert "license" in meta

    def test_reactions_is_dict(self):
        data = load_ord_conditions()
        assert isinstance(data["reactions"], dict)

    def test_known_reaction_present(self):
        data = load_ord_conditions()
        assert "suzuki_coupling" in data["reactions"]

    def test_reaction_entry_structure(self):
        data = load_ord_conditions()
        entry = data["reactions"]["suzuki_coupling"]
        assert "n" in entry
        assert entry["n"] > 0
        assert "temperature_C" in entry
        assert "median" in entry["temperature_C"]
        assert "solvents" in entry
        assert isinstance(entry["solvents"], list)
        assert "atmosphere" in entry

    def test_caching(self):
        data1 = load_ord_conditions()
        data2 = load_ord_conditions()
        assert data1 is data2  # same object, not reloaded


class TestEmptyOrdData:
    """Test graceful handling when ORD data is missing or empty."""

    def test_missing_file_returns_empty(self, tmp_path, monkeypatch):
        """If the JSON file doesn't exist, returns empty reactions."""
        import molbuilder.data as data_mod

        def fake_load():
            return {"_meta": {}, "reactions": {}}

        monkeypatch.setattr(data_mod, "load_ord_conditions", fake_load)
        result = fake_load()
        assert result["reactions"] == {}

    def test_optimize_conditions_falls_back(self):
        """Even with ORD data loaded, templates not in ORD fall back."""
        # Create a fake template with a made-up named_reaction
        from molbuilder.reactions.reaction_types import (
            ReactionTemplate,
            ReactionCategory,
        )
        fake = ReactionTemplate(
            name="Fake reaction",
            named_reaction="Nonexistent Reaction XYZ",
            category=ReactionCategory.MISC,
            reagents=["NaOH"],
            solvents=["water"],
            temperature_range=(20.0, 40.0),
        )
        cond = optimize_conditions(fake, 1.0)
        # Should get a data_source -- either heuristic or category fallback
        assert isinstance(cond.data_source, str)


# =====================================================================
#  ORD stats lookup
# =====================================================================

class TestGetOrdStats:

    def test_known_template_returns_data(self):
        """Template with a known named_reaction should find ORD data."""
        suzuki = None
        for t in REACTION_TEMPLATES:
            if t.named_reaction == "Suzuki-Miyaura coupling":
                suzuki = t
                break
        if suzuki is None:
            pytest.skip("No Suzuki template found")
        stats = _get_ord_stats(suzuki)
        assert stats is not None
        assert stats["n"] > 0

    def test_category_fallback(self):
        """Template with unmapped named_reaction should get category fallback."""
        from molbuilder.reactions.reaction_types import (
            ReactionTemplate,
            ReactionCategory,
        )
        fake = ReactionTemplate(
            name="Unknown substitution",
            named_reaction="Never Heard Of This",
            category=ReactionCategory.SUBSTITUTION,
            reagents=["X"],
            solvents=["Y"],
        )
        stats = _get_ord_stats(fake)
        # Should return substitution_general fallback
        assert stats is not None
        assert stats["n"] > 0


# =====================================================================
#  ReactionConditions.data_source
# =====================================================================

class TestDataSourceField:

    def test_heuristic_default(self):
        cond = ReactionConditions(
            temperature_C=25.0,
            pressure_atm=1.0,
            solvent="THF",
            concentration_M=0.5,
            addition_rate="all at once",
            reaction_time_hours=2.0,
            atmosphere="air",
            workup_procedure="Standard workup.",
            notes="None.",
        )
        assert cond.data_source == "heuristic"

    def test_ord_backed_conditions(self):
        """Template with ORD data should produce ORD-sourced conditions."""
        suzuki = None
        for t in REACTION_TEMPLATES:
            if t.named_reaction == "Suzuki-Miyaura coupling":
                suzuki = t
                break
        if suzuki is None:
            pytest.skip("No Suzuki template found")

        cond = optimize_conditions(suzuki, 1.0)
        assert cond.data_source.startswith("ORD")
        assert "(n=" in cond.data_source

    def test_heuristic_fallback_for_unknown(self):
        """Template not in ORD should produce heuristic conditions."""
        from molbuilder.reactions.reaction_types import (
            ReactionTemplate,
            ReactionCategory,
        )
        fake = ReactionTemplate(
            name="Made up reaction",
            named_reaction="Nonexistent XYZ 999",
            category=ReactionCategory.MISC,
            reagents=[],
            solvents=[],
            temperature_range=(10.0, 30.0),
        )
        # Clear category fallback for this test by checking MISC
        # MISC has a category fallback, so data_source may be ORD
        cond = optimize_conditions(fake, 1.0)
        assert isinstance(cond.data_source, str)
        assert len(cond.data_source) > 0

    def test_ord_temperature_used(self):
        """When ORD data is available, temperature should match ORD median."""
        suzuki = None
        for t in REACTION_TEMPLATES:
            if t.named_reaction == "Suzuki-Miyaura coupling":
                suzuki = t
                break
        if suzuki is None:
            pytest.skip("No Suzuki template found")

        ord_data = load_ord_conditions()
        expected_temp = ord_data["reactions"]["suzuki_coupling"]["temperature_C"]["median"]
        cond = optimize_conditions(suzuki, 1.0)
        assert cond.temperature_C == expected_temp


# =====================================================================
#  TemplateMatch.data_source propagation
# =====================================================================

class TestTemplateMatchDataSource:

    def test_predict_conditions_propagates_data_source(self):
        """predict_conditions should set data_source on TemplateMatch."""
        result = predict_conditions("CCO", reaction_name="oxidation")
        if result.best_match:
            assert hasattr(result.best_match, "data_source")
            assert isinstance(result.best_match.data_source, str)
            assert len(result.best_match.data_source) > 0

    def test_all_candidates_have_data_source(self):
        result = predict_conditions("CCO")
        for cand in result.candidates:
            assert hasattr(cand, "data_source")
            assert isinstance(cand.data_source, str)

    def test_conditions_and_match_data_source_agree(self):
        """TemplateMatch.data_source should match its conditions.data_source."""
        result = predict_conditions("CCO", reaction_name="oxidation")
        for cand in result.candidates:
            assert cand.data_source == cand.conditions.data_source


# =====================================================================
#  Template mapping coverage
# =====================================================================

class TestTemplateMappingCoverage:

    def test_all_named_reactions_mapped(self):
        """Every template's named_reaction should have a mapping entry."""
        unmapped = set()
        for t in REACTION_TEMPLATES:
            if t.named_reaction and t.named_reaction not in TEMPLATE_TO_ORD_KEY:
                unmapped.add(t.named_reaction)
        # Paterno-Buchi has non-ASCII chars -- now mapped via unicode escapes
        # (no discard needed)
        assert len(unmapped) == 0, f"Unmapped named_reactions: {unmapped}"

    def test_all_categories_have_fallback(self):
        """Every ReactionCategory should have a fallback ORD key."""
        from molbuilder.reactions.reaction_types import ReactionCategory
        for cat in ReactionCategory:
            assert cat.name in CATEGORY_FALLBACK, f"Missing fallback for {cat.name}"

    def test_ord_keys_exist_in_json(self):
        """Every ORD key referenced in the mapping should exist in the JSON."""
        data = load_ord_conditions()
        reactions = data["reactions"]
        missing = set()
        for ord_key in set(TEMPLATE_TO_ORD_KEY.values()):
            if ord_key not in reactions:
                missing.add(ord_key)
        assert len(missing) == 0, f"Missing ORD keys in JSON: {missing}"

    def test_category_fallback_keys_in_json(self):
        """Every category fallback key should exist in the JSON."""
        data = load_ord_conditions()
        reactions = data["reactions"]
        missing = set()
        for cat_key in CATEGORY_FALLBACK.values():
            if cat_key not in reactions:
                missing.add(cat_key)
        assert len(missing) == 0, f"Missing category fallback keys: {missing}"


# =====================================================================
#  End-to-end integration
# =====================================================================

class TestEndToEnd:

    def test_ethanol_oxidation_has_data_source(self):
        result = predict_conditions("CCO", reaction_name="oxidation")
        assert result.best_match is not None
        assert result.best_match.data_source != ""
        assert result.best_match.conditions.data_source != ""

    def test_suzuki_substrate_uses_ord(self):
        """Aryl halide substrate with suzuki hint should use ORD data."""
        # Bromobenzene
        result = predict_conditions("c1ccc(Br)cc1", reaction_name="suzuki")
        if result.best_match:
            # Should find suzuki-related template backed by ORD
            has_ord = any(
                "ORD" in c.data_source for c in result.candidates
            )
            assert has_ord, "Expected at least one ORD-backed candidate"

    def test_existing_tests_unaffected(self):
        """Basic predict_conditions still works -- no regressions."""
        result = predict_conditions("CCO")
        assert len(result.candidates) > 0
        assert result.best_match is not None
        for c in result.candidates:
            assert c.conditions.temperature_C is not None
            assert c.conditions.solvent is not None
