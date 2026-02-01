"""Tests for molbuilder.process -- reactor, solvent, purification,
conditions, costing, safety, and scale-up modules.

Run with:
    python -m pytest tests/test_process.py
"""

import unittest

from molbuilder.reactions.knowledge_base import lookup_by_name
from molbuilder.reactions.reaction_types import ReactionCategory, ReactionTemplate

from molbuilder.process.reactor import ReactorType, ReactorSpec, select_reactor
from molbuilder.process.solvent_systems import SolventRecommendation, select_solvent
from molbuilder.process.purification import (
    PurificationMethod,
    PurificationStep,
    recommend_purification,
)
from molbuilder.process.conditions import ReactionConditions, optimize_conditions
from molbuilder.process.costing import CostBreakdown, CostEstimate, estimate_cost
from molbuilder.process.safety import (
    GHS_PICTOGRAMS,
    GHS_HAZARD_STATEMENTS,
    HazardInfo,
    SafetyAssessment,
    assess_safety,
)
from molbuilder.process.scale_up import ScaleUpAnalysis, analyze_scale_up


# =====================================================================
#  Mock helpers for costing / safety / scale-up (duck typing)
# =====================================================================

class _MockPrecursor:
    def __init__(self, name, smiles):
        self.name = name
        self.smiles = smiles
        self.cost_per_kg = 10.0


class _MockStep:
    def __init__(self, template):
        self.template = template
        self.precursors = [_MockPrecursor("ethanol", "CCO")]
        self.step_number = 1


# =====================================================================
#  Helpers to retrieve real templates
# =====================================================================

def _get_template(name_query):
    """Return the first template matching *name_query*, or raise."""
    results = lookup_by_name(name_query)
    if not results:
        raise RuntimeError(
            "No template found for query %r -- check knowledge_base" % name_query
        )
    return results[0]


# =====================================================================
#  Test Reactor
# =====================================================================

class TestReactor(unittest.TestCase):
    """Tests for molbuilder.process.reactor."""

    def test_reactor_type_enum_members(self):
        """ReactorType should expose BATCH, CSTR, PFR, MICROREACTOR."""
        self.assertIsNotNone(ReactorType.BATCH)
        self.assertIsNotNone(ReactorType.CSTR)
        self.assertIsNotNone(ReactorType.PFR)
        self.assertIsNotNone(ReactorType.MICROREACTOR)

    def test_reactor_type_enum_has_semi_batch_and_fixed_bed(self):
        """ReactorType should also have SEMI_BATCH and FIXED_BED."""
        self.assertIsNotNone(ReactorType.SEMI_BATCH)
        self.assertIsNotNone(ReactorType.FIXED_BED)

    def test_select_reactor_returns_reactor_spec(self):
        """select_reactor should return a ReactorSpec instance."""
        template = _get_template("NaBH4 reduction")
        spec = select_reactor(template, scale_kg=1.0)
        self.assertIsInstance(spec, ReactorSpec)

    def test_select_reactor_small_scale_batch(self):
        """A standard 1 kg reaction should default to BATCH."""
        template = _get_template("NaBH4 reduction")
        spec = select_reactor(template, scale_kg=1.0)
        self.assertEqual(spec.reactor_type, ReactorType.BATCH)

    def test_reactor_spec_fields(self):
        """ReactorSpec should carry all expected fields."""
        template = _get_template("NaBH4 reduction")
        spec = select_reactor(template, scale_kg=1.0)

        self.assertIsInstance(spec.reactor_type, ReactorType)
        self.assertGreater(spec.volume_L, 0)
        self.assertIsInstance(spec.temperature_C, float)
        self.assertGreater(spec.pressure_atm, 0)
        self.assertGreater(spec.residence_time_min, 0)
        self.assertIn(spec.mixing_type, ("mechanical", "static", "none"))
        self.assertIn(spec.heat_transfer, ("jacketed", "coil", "adiabatic"))
        self.assertIn(spec.material, ("glass", "stainless steel", "hastelloy"))
        self.assertGreater(spec.estimated_cost_usd, 0)
        self.assertIsInstance(spec.notes, str)

    def test_exothermic_small_gives_microreactor(self):
        """A fast/exothermic reaction at sub-kg scale should give MICROREACTOR."""
        template = _get_template("HBr addition")  # ADDITION category
        spec = select_reactor(template, scale_kg=0.5)
        self.assertEqual(spec.reactor_type, ReactorType.MICROREACTOR)

    def test_large_scale_commodity_gives_pfr(self):
        """scale > 500 kg, non-exothermic, non-catalytic should give PFR."""
        # SN1 solvolysis: SUBSTITUTION, no catalysts, no reagents
        template = _get_template("SN1 tertiary halide")
        spec = select_reactor(template, scale_kg=600.0)
        self.assertEqual(spec.reactor_type, ReactorType.PFR)


# =====================================================================
#  Test Solvent Systems
# =====================================================================

class TestSolventSystems(unittest.TestCase):
    """Tests for molbuilder.process.solvent_systems."""

    def test_select_solvent_returns_list(self):
        """select_solvent should return a list."""
        template = _get_template("NaBH4 reduction")
        result = select_solvent(template, scale_kg=1.0)
        self.assertIsInstance(result, list)

    def test_select_solvent_nonempty(self):
        """The solvent list should not be empty for a common reaction."""
        template = _get_template("NaBH4 reduction")
        result = select_solvent(template, scale_kg=1.0)
        self.assertGreater(len(result), 0)

    def test_recommendation_instance(self):
        """Each element should be a SolventRecommendation."""
        template = _get_template("Suzuki coupling")
        recs = select_solvent(template, scale_kg=1.0)
        for rec in recs:
            self.assertIsInstance(rec, SolventRecommendation)

    def test_recommendation_fields(self):
        """SolventRecommendation should have score, solvent_name, green_score."""
        template = _get_template("Suzuki coupling")
        recs = select_solvent(template, scale_kg=1.0)
        self.assertGreater(len(recs), 0)
        rec = recs[0]
        self.assertIsInstance(rec.solvent_name, str)
        self.assertIsInstance(rec.score, float)
        self.assertIsInstance(rec.green_score, int)

    def test_scores_in_range(self):
        """Scores should be between 0 and 100 inclusive."""
        template = _get_template("NaBH4 reduction")
        recs = select_solvent(template, scale_kg=1.0)
        for rec in recs:
            self.assertGreaterEqual(rec.score, 0)
            self.assertLessEqual(rec.score, 100)

    def test_green_score_in_range(self):
        """green_score should be between 1 and 10."""
        template = _get_template("NaBH4 reduction")
        recs = select_solvent(template, scale_kg=1.0)
        for rec in recs:
            self.assertGreaterEqual(rec.green_score, 1)
            self.assertLessEqual(rec.green_score, 10)

    def test_sorted_descending(self):
        """Recommendations should be sorted by score descending."""
        template = _get_template("NaBH4 reduction")
        recs = select_solvent(template, scale_kg=1.0)
        if len(recs) > 1:
            scores = [r.score for r in recs]
            self.assertEqual(scores, sorted(scores, reverse=True))


# =====================================================================
#  Test Purification
# =====================================================================

class TestPurification(unittest.TestCase):
    """Tests for molbuilder.process.purification."""

    def test_purification_method_enum_members(self):
        """PurificationMethod should have DISTILLATION, RECRYSTALLIZATION,
        and COLUMN_CHROMATOGRAPHY."""
        self.assertIsNotNone(PurificationMethod.DISTILLATION)
        self.assertIsNotNone(PurificationMethod.RECRYSTALLIZATION)
        self.assertIsNotNone(PurificationMethod.COLUMN_CHROMATOGRAPHY)

    def test_recommend_purification_returns_list(self):
        """recommend_purification should return a list of PurificationStep."""
        template = _get_template("NaBH4 reduction")
        steps = recommend_purification(template, scale_kg=1.0)
        self.assertIsInstance(steps, list)
        self.assertGreater(len(steps), 0)

    def test_purification_step_instance(self):
        """Each step should be a PurificationStep."""
        template = _get_template("Suzuki coupling")
        steps = recommend_purification(template, scale_kg=1.0)
        for step in steps:
            self.assertIsInstance(step, PurificationStep)

    def test_purification_step_fields(self):
        """PurificationStep should have method, estimated_recovery,
        estimated_purity."""
        template = _get_template("NaBH4 reduction")
        steps = recommend_purification(template, scale_kg=1.0)
        self.assertGreater(len(steps), 0)
        step = steps[0]
        self.assertIsInstance(step.method, PurificationMethod)
        self.assertIsInstance(step.estimated_recovery, float)
        self.assertIsInstance(step.estimated_purity, float)
        self.assertGreater(step.estimated_recovery, 0)
        self.assertLessEqual(step.estimated_recovery, 100)
        self.assertGreater(step.estimated_purity, 0)
        self.assertLessEqual(step.estimated_purity, 100)

    def test_solid_product_includes_recrystallization(self):
        """Solid products (COUPLING) should include extraction +
        recrystallization."""
        # Suzuki coupling is COUPLING category, in _SOLID_PRODUCT_CATEGORIES
        template = _get_template("Suzuki coupling")
        steps = recommend_purification(template, scale_kg=1.0)
        methods = [s.method for s in steps]
        # Has catalyst -> filtration first, then extraction, then recryst
        self.assertIn(PurificationMethod.FILTRATION, methods)
        self.assertIn(PurificationMethod.RECRYSTALLIZATION, methods)

    def test_catalytic_starts_with_filtration(self):
        """Catalytic reactions should start with a filtration step."""
        template = _get_template("Suzuki coupling")
        steps = recommend_purification(template, scale_kg=1.0)
        self.assertEqual(steps[0].method, PurificationMethod.FILTRATION)


# =====================================================================
#  Test Conditions
# =====================================================================

class TestConditions(unittest.TestCase):
    """Tests for molbuilder.process.conditions."""

    def test_optimize_conditions_returns_dataclass(self):
        """optimize_conditions should return a ReactionConditions."""
        template = _get_template("NaBH4 reduction")
        cond = optimize_conditions(template, scale_kg=1.0)
        self.assertIsInstance(cond, ReactionConditions)

    def test_conditions_fields(self):
        """ReactionConditions should have temperature_C, solvent,
        atmosphere, workup_procedure."""
        template = _get_template("NaBH4 reduction")
        cond = optimize_conditions(template, scale_kg=1.0)
        self.assertIsInstance(cond.temperature_C, float)
        self.assertIsInstance(cond.solvent, str)
        self.assertIsInstance(cond.atmosphere, str)
        self.assertIsInstance(cond.workup_procedure, str)
        self.assertGreater(len(cond.workup_procedure), 0)

    def test_temperature_within_template_range(self):
        """Optimised temperature should be within the template's range."""
        template = _get_template("NaBH4 reduction")
        cond = optimize_conditions(template, scale_kg=1.0)
        lo, hi = template.temperature_range
        self.assertGreaterEqual(cond.temperature_C, lo)
        self.assertLessEqual(cond.temperature_C, hi)

    def test_inert_atmosphere_for_coupling(self):
        """Coupling reactions should require inert atmosphere (N2)."""
        template = _get_template("Suzuki coupling")
        cond = optimize_conditions(template, scale_kg=1.0)
        self.assertIn(cond.atmosphere, ("N2", "Ar"))

    def test_pressure_is_positive(self):
        """Pressure should be positive."""
        template = _get_template("Fischer esterification")
        cond = optimize_conditions(template, scale_kg=1.0)
        self.assertGreater(cond.pressure_atm, 0)


# =====================================================================
#  Test Costing
# =====================================================================

class TestCosting(unittest.TestCase):
    """Tests for molbuilder.process.costing."""

    def test_estimate_cost_returns_cost_estimate(self):
        """estimate_cost should return a CostEstimate."""
        template = _get_template("NaBH4 reduction")
        step = _MockStep(template)
        result = estimate_cost([step], scale_kg=1.0)
        self.assertIsInstance(result, CostEstimate)

    def test_total_usd_positive(self):
        """total_usd should be positive for a real step."""
        template = _get_template("NaBH4 reduction")
        step = _MockStep(template)
        result = estimate_cost([step], scale_kg=1.0)
        self.assertGreater(result.total_usd, 0)

    def test_per_kg_usd_positive(self):
        """per_kg_usd should be positive."""
        template = _get_template("NaBH4 reduction")
        step = _MockStep(template)
        result = estimate_cost([step], scale_kg=1.0)
        self.assertGreater(result.per_kg_usd, 0)

    def test_breakdown_fields(self):
        """CostBreakdown should have raw_materials_usd, labor_usd, etc."""
        template = _get_template("NaBH4 reduction")
        step = _MockStep(template)
        result = estimate_cost([step], scale_kg=1.0)
        bd = result.breakdown
        self.assertIsInstance(bd, CostBreakdown)
        self.assertIsInstance(bd.raw_materials_usd, float)
        self.assertIsInstance(bd.labor_usd, float)
        self.assertIsInstance(bd.equipment_usd, float)
        self.assertIsInstance(bd.energy_usd, float)
        self.assertIsInstance(bd.waste_disposal_usd, float)
        self.assertIsInstance(bd.overhead_usd, float)

    def test_empty_steps_returns_zero(self):
        """No steps should return a CostEstimate with total 0."""
        result = estimate_cost([], scale_kg=1.0)
        self.assertEqual(result.total_usd, 0.0)

    def test_multi_step_costing(self):
        """Multiple steps should produce a higher cost than a single step."""
        t1 = _get_template("NaBH4 reduction")
        t2 = _get_template("Suzuki coupling")
        single = estimate_cost([_MockStep(t1)], scale_kg=1.0)
        multi = estimate_cost([_MockStep(t1), _MockStep(t2)], scale_kg=1.0)
        self.assertGreater(multi.total_usd, single.total_usd)


# =====================================================================
#  Test Safety
# =====================================================================

class TestSafety(unittest.TestCase):
    """Tests for molbuilder.process.safety."""

    def test_ghs_pictograms_populated(self):
        """GHS_PICTOGRAMS should be a non-empty dict."""
        self.assertIsInstance(GHS_PICTOGRAMS, dict)
        self.assertGreater(len(GHS_PICTOGRAMS), 0)

    def test_ghs_hazard_statements_populated(self):
        """GHS_HAZARD_STATEMENTS should be a non-empty dict."""
        self.assertIsInstance(GHS_HAZARD_STATEMENTS, dict)
        self.assertGreater(len(GHS_HAZARD_STATEMENTS), 0)

    def test_assess_safety_returns_list(self):
        """assess_safety should return a list of SafetyAssessment."""
        template = _get_template("NaBH4 reduction")
        step = _MockStep(template)
        results = assess_safety([step])
        self.assertIsInstance(results, list)
        self.assertEqual(len(results), 1)

    def test_safety_assessment_instance(self):
        """Each element should be a SafetyAssessment."""
        template = _get_template("NaBH4 reduction")
        step = _MockStep(template)
        results = assess_safety([step])
        self.assertIsInstance(results[0], SafetyAssessment)

    def test_safety_assessment_fields(self):
        """SafetyAssessment should have risk_level, ppe_required, hazards."""
        template = _get_template("NaBH4 reduction")
        step = _MockStep(template)
        results = assess_safety([step])
        sa = results[0]
        self.assertIn(sa.risk_level, ("low", "medium", "high"))
        self.assertIsInstance(sa.ppe_required, list)
        self.assertGreater(len(sa.ppe_required), 0)
        self.assertIsInstance(sa.hazards, list)

    def test_hazard_info_for_each_reagent(self):
        """There should be at least one HazardInfo per reagent/catalyst."""
        template = _get_template("NaBH4 reduction")
        step = _MockStep(template)
        results = assess_safety([step])
        sa = results[0]
        expected_count = len(template.reagents) + len(template.catalysts)
        self.assertEqual(len(sa.hazards), expected_count)
        for hi in sa.hazards:
            self.assertIsInstance(hi, HazardInfo)

    def test_multi_step_safety(self):
        """assess_safety with two steps should return two assessments."""
        t1 = _get_template("NaBH4 reduction")
        t2 = _get_template("Suzuki coupling")
        results = assess_safety([_MockStep(t1), _MockStep(t2)])
        self.assertEqual(len(results), 2)
        self.assertEqual(results[0].step_number, 1)
        self.assertEqual(results[1].step_number, 2)


# =====================================================================
#  Test Scale-Up
# =====================================================================

class TestScaleUp(unittest.TestCase):
    """Tests for molbuilder.process.scale_up."""

    def test_analyze_scale_up_returns_dataclass(self):
        """analyze_scale_up should return a ScaleUpAnalysis."""
        template = _get_template("NaBH4 reduction")
        step = _MockStep(template)
        result = analyze_scale_up([step], target_annual_kg=1000.0)
        self.assertIsInstance(result, ScaleUpAnalysis)

    def test_scale_up_fields(self):
        """ScaleUpAnalysis should have recommended_mode, annual_capacity_kg,
        capital_cost_usd."""
        template = _get_template("NaBH4 reduction")
        step = _MockStep(template)
        result = analyze_scale_up([step], target_annual_kg=1000.0)
        self.assertIn(result.recommended_mode, ("batch", "continuous", "semi-continuous"))
        self.assertGreater(result.annual_capacity_kg, 0)
        self.assertGreater(result.capital_cost_usd, 0)

    def test_small_target_gives_batch(self):
        """A modest annual target should recommend batch mode."""
        template = _get_template("NaBH4 reduction")
        step = _MockStep(template)
        result = analyze_scale_up([step], target_annual_kg=500.0)
        self.assertEqual(result.recommended_mode, "batch")

    def test_cycle_time_positive(self):
        """cycle_time_hours should be positive."""
        template = _get_template("NaBH4 reduction")
        step = _MockStep(template)
        result = analyze_scale_up([step], target_annual_kg=1000.0)
        self.assertGreater(result.cycle_time_hours, 0)

    def test_operating_cost_positive(self):
        """operating_cost_annual_usd should be positive."""
        template = _get_template("NaBH4 reduction")
        step = _MockStep(template)
        result = analyze_scale_up([step], target_annual_kg=1000.0)
        self.assertGreater(result.operating_cost_annual_usd, 0)

    def test_risks_is_list(self):
        """scale_up_risks should be a list of strings."""
        template = _get_template("NaBH4 reduction")
        step = _MockStep(template)
        result = analyze_scale_up([step], target_annual_kg=1000.0)
        self.assertIsInstance(result.scale_up_risks, list)

    def test_recommendations_is_list(self):
        """recommendations should be a non-empty list."""
        template = _get_template("NaBH4 reduction")
        step = _MockStep(template)
        result = analyze_scale_up([step], target_annual_kg=1000.0)
        self.assertIsInstance(result.recommendations, list)
        self.assertGreater(len(result.recommendations), 0)


if __name__ == "__main__":
    unittest.main()
