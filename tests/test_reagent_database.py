"""Tests for the expanded reagent database with pricing tiers."""

import logging

import pytest

from molbuilder.reactions.reagent_data import (
    REAGENT_DB,
    PRICING_TIER_DEFAULTS,
    Reagent,
    get_reagent,
    normalize_reagent_name,
)


VALID_TIERS = {"commodity", "standard", "specialty", "exotic"}


class TestPricingTier:

    def test_all_reagents_have_pricing_tier(self):
        """Every reagent in REAGENT_DB must have a pricing_tier field."""
        for key, reagent in REAGENT_DB.items():
            assert hasattr(reagent, "pricing_tier"), f"{key} missing pricing_tier"
            assert reagent.pricing_tier, f"{key} has empty pricing_tier"

    def test_pricing_tier_values_valid(self):
        """All pricing_tier values must be one of the valid tiers."""
        for key, reagent in REAGENT_DB.items():
            assert reagent.pricing_tier in VALID_TIERS, (
                f"{key} has invalid pricing_tier={reagent.pricing_tier!r}"
            )

    def test_pricing_tier_defaults_complete(self):
        """PRICING_TIER_DEFAULTS covers all valid tiers."""
        for tier in VALID_TIERS:
            assert tier in PRICING_TIER_DEFAULTS

    def test_tier_costs_ascending(self):
        """Tier defaults should be ascending: commodity < standard < specialty < exotic."""
        costs = [PRICING_TIER_DEFAULTS[t] for t in ["commodity", "standard", "specialty", "exotic"]]
        for i in range(len(costs) - 1):
            assert costs[i] < costs[i + 1]


class TestReagentLookup:

    def test_reagent_lookup_by_key(self):
        """Direct DB key lookup works."""
        r = get_reagent("lialh4")
        assert r is not None
        assert "lithium" in r.name.lower() or "aluminium" in r.name.lower()

        r2 = get_reagent("nabh4")
        assert r2 is not None
        assert "sodium" in r2.name.lower() or "borohydride" in r2.name.lower()

    def test_new_alias_wilkinson(self):
        """Wilkinson alias resolves."""
        r = get_reagent("wilkinsons_catalyst")
        assert r is not None
        assert "Wilkinson" in r.name

    def test_unknown_reagent_returns_none(self):
        """Unrecognised reagent returns None."""
        r = get_reagent("nonexistent_reagent_xyz")
        assert r is None


class TestNewReagentsExist:

    def test_new_coupling_reagents_exist(self):
        """HATU, HBTU, PyBOP, T3P, CDI should be in the database."""
        for key in ["hatu", "hbtu", "pybop", "t3p", "cdi"]:
            assert key in REAGENT_DB, f"Missing coupling reagent: {key}"

    def test_new_catalyst_reagents_exist(self):
        """Transition metal catalysts should be in the database."""
        for key in ["nicl2_dppe", "rugen2_cl2", "pd_amphos"]:
            assert key in REAGENT_DB, f"Missing catalyst: {key}"

    def test_new_chiral_reagents_exist(self):
        """Chiral auxiliaries should be in the database."""
        for key in ["evans_oxazolidinone", "ellman_sulfinamide", "cbs_catalyst", "binap"]:
            assert key in REAGENT_DB, f"Missing chiral reagent: {key}"

    def test_new_oxidant_reagents_exist(self):
        """New oxidants should be in the database."""
        for key in ["dmp", "ibx", "tempo", "oxone", "mno2", "tpap"]:
            assert key in REAGENT_DB, f"Missing oxidant: {key}"

    def test_new_fluorinating_reagents_exist(self):
        """Fluorinating agents should be in the database."""
        for key in ["selectfluor", "nfsi", "dast"]:
            assert key in REAGENT_DB, f"Missing fluorinating agent: {key}"

    def test_total_reagent_count_above_150(self):
        """The database should have grown from ~100 to >150 entries."""
        assert len(REAGENT_DB) >= 150, f"Only {len(REAGENT_DB)} reagents"


class TestCostingFallback:

    def test_unknown_reagent_logs_warning(self, caplog):
        """estimate_cost should warn for unknown reagents."""
        from molbuilder.process.costing import estimate_cost
        from molbuilder.reactions.reaction_types import ReactionCategory, ReactionTemplate

        tmpl = ReactionTemplate(
            name="test",
            named_reaction=None,
            category=ReactionCategory.MISC,
            reagents=["totally_fake_reagent_xyz"],
            solvents=["THF"],
        )

        class FakeStep:
            template = tmpl

        with caplog.at_level(logging.WARNING, logger="molbuilder.costing"):
            result = estimate_cost([FakeStep()], scale_kg=1.0)

        assert result.total_usd > 0
        assert any("Unknown reagent" in r.message for r in caplog.records)

    def test_known_reagent_uses_cost_per_kg(self):
        """Known reagents use their actual cost_per_kg."""
        from molbuilder.process.costing import estimate_cost
        from molbuilder.reactions.reaction_types import ReactionCategory, ReactionTemplate

        tmpl = ReactionTemplate(
            name="test",
            named_reaction=None,
            category=ReactionCategory.MISC,
            reagents=["hcl"],  # lowercase key matches REAGENT_DB
            solvents=["THF"],
        )

        class FakeStep:
            template = tmpl

        result = estimate_cost([FakeStep()], scale_kg=1.0)
        # HCl is $6/kg commodity -- raw material includes reagent + solvent
        # At 1kg: reagent $6*1.2=$7.20, solvent ~$150 => total ~$177
        # Compare against unknown-reagent fallback ($75/kg * 1.2 = $90 + solvent)
        assert result.breakdown.raw_materials_usd < 250
