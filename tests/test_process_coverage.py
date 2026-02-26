"""Branch-coverage tests for molbuilder.process subpackage.

Targets every uncovered branch in purification, conditions, safety, and
scale_up using directly constructed ReactionTemplate objects.  No database
lookups required -- all templates are built via ``_make_template``.

Run with:
    python -m pytest tests/test_process_coverage.py -v

Combined coverage check:
    python -m pytest tests/test_process.py tests/test_process_coverage.py \
        --cov=molbuilder/process --cov-report=term-missing
"""

import pytest

from molbuilder.reactions.reaction_types import ReactionCategory, ReactionTemplate

# -- purification ---------------------------------------------------------
from molbuilder.process.purification import (
    PurificationMethod,
    recommend_purification,
)

# -- conditions -----------------------------------------------------------
from molbuilder.process.conditions import (
    ReactionConditions,
    optimize_conditions,
    _select_atmosphere,
    _addition_rate,
    _workup_procedure,
    _select_solvent,
)

# -- safety ---------------------------------------------------------------
from molbuilder.process.safety import (
    _determine_ppe,
    _determine_engineering_controls,
    _determine_emergency_procedures,
    _determine_incompatibilities,
    _classify_waste,
    _calculate_risk_level,
)

# -- scale_up -------------------------------------------------------------
from molbuilder.process.scale_up import (
    analyze_scale_up,
    _estimate_cycle_time,
    _is_continuous_candidate,
    _identify_risks,
    _generate_recommendations,
)


# =====================================================================
#  Test helpers
# =====================================================================

def _make_template(
    category,
    reagents=None,
    catalysts=None,
    solvents=None,
    temperature_range=(20, 25),
    typical_yield=(60, 90),
    safety_notes="",
):
    """Build a ReactionTemplate with controllable fields."""
    return ReactionTemplate(
        name=f"Test {category.name}",
        named_reaction=None,
        category=category,
        reagents=reagents or [],
        solvents=solvents or [],
        catalysts=catalysts or [],
        temperature_range=temperature_range,
        typical_yield=typical_yield,
        safety_notes=safety_notes,
    )


class _MockStep:
    """Minimal duck-typed synthesis step."""

    def __init__(self, template):
        self.template = template
        self.precursors = []
        self.step_number = 1


# =====================================================================
#  1. Purification branches
# =====================================================================

class TestPurificationBranches:
    """Exercise every category routing branch + scale effects."""

    @staticmethod
    def _methods(steps):
        return [s.method for s in steps]

    # -- Liquid products (ELIMINATION, RADICAL) --

    def test_elimination_gives_extraction_distillation(self):
        t = _make_template(ReactionCategory.ELIMINATION)
        m = self._methods(recommend_purification(t, scale_kg=1.0))
        assert m == [PurificationMethod.EXTRACTION, PurificationMethod.DISTILLATION]

    def test_radical_gives_extraction_distillation(self):
        t = _make_template(ReactionCategory.RADICAL)
        m = self._methods(recommend_purification(t, scale_kg=1.0))
        assert m == [PurificationMethod.EXTRACTION, PurificationMethod.DISTILLATION]

    # -- Solid products (COUPLING no catalyst, PROTECTION, CARBONYL, PERICYCLIC) --

    def test_coupling_no_catalyst_gives_extraction_recrystallization(self):
        t = _make_template(ReactionCategory.COUPLING)
        m = self._methods(recommend_purification(t, scale_kg=1.0))
        assert m == [PurificationMethod.EXTRACTION, PurificationMethod.RECRYSTALLIZATION]

    def test_protection_gives_extraction_recrystallization(self):
        t = _make_template(ReactionCategory.PROTECTION)
        m = self._methods(recommend_purification(t, scale_kg=1.0))
        assert m == [PurificationMethod.EXTRACTION, PurificationMethod.RECRYSTALLIZATION]

    # -- Complex mixture, small scale (REARRANGEMENT, MISC) --

    def test_rearrangement_small_gives_extraction_flash(self):
        t = _make_template(ReactionCategory.REARRANGEMENT)
        m = self._methods(recommend_purification(t, scale_kg=1.0))
        assert m == [PurificationMethod.EXTRACTION, PurificationMethod.FLASH_CHROMATOGRAPHY]

    def test_misc_small_gives_extraction_flash(self):
        t = _make_template(ReactionCategory.MISC)
        m = self._methods(recommend_purification(t, scale_kg=1.0))
        assert m == [PurificationMethod.EXTRACTION, PurificationMethod.FLASH_CHROMATOGRAPHY]

    # -- Complex mixture, large scale (>5 kg) --

    def test_rearrangement_large_gives_extraction_precip_recryst(self):
        t = _make_template(ReactionCategory.REARRANGEMENT)
        m = self._methods(recommend_purification(t, scale_kg=10.0))
        assert m == [
            PurificationMethod.EXTRACTION,
            PurificationMethod.PRECIPITATION,
            PurificationMethod.RECRYSTALLIZATION,
        ]

    def test_misc_large_gives_extraction_precip_recryst(self):
        t = _make_template(ReactionCategory.MISC)
        m = self._methods(recommend_purification(t, scale_kg=10.0))
        assert m == [
            PurificationMethod.EXTRACTION,
            PurificationMethod.PRECIPITATION,
            PurificationMethod.RECRYSTALLIZATION,
        ]

    # -- Oxidation / Reduction --

    def test_oxidation_small_gives_extraction_flash(self):
        """scale < 0.5 and chromatography-appropriate -> flash."""
        t = _make_template(ReactionCategory.OXIDATION)
        m = self._methods(recommend_purification(t, scale_kg=0.3))
        assert m == [PurificationMethod.EXTRACTION, PurificationMethod.FLASH_CHROMATOGRAPHY]

    def test_oxidation_normal_gives_extraction_distillation(self):
        """scale >= 0.5 -> distillation."""
        t = _make_template(ReactionCategory.OXIDATION)
        m = self._methods(recommend_purification(t, scale_kg=1.0))
        assert m == [PurificationMethod.EXTRACTION, PurificationMethod.DISTILLATION]

    def test_reduction_normal_gives_extraction_distillation(self):
        t = _make_template(ReactionCategory.REDUCTION)
        m = self._methods(recommend_purification(t, scale_kg=1.0))
        assert m == [PurificationMethod.EXTRACTION, PurificationMethod.DISTILLATION]

    # -- Substitution / Addition --

    def test_substitution_small_gives_extraction_flash(self):
        """scale < 1.0 and chromatography-appropriate -> flash."""
        t = _make_template(ReactionCategory.SUBSTITUTION)
        m = self._methods(recommend_purification(t, scale_kg=0.5))
        assert m == [PurificationMethod.EXTRACTION, PurificationMethod.FLASH_CHROMATOGRAPHY]

    def test_substitution_normal_gives_extraction_distillation(self):
        """scale >= 1.0 -> distillation."""
        t = _make_template(ReactionCategory.SUBSTITUTION)
        m = self._methods(recommend_purification(t, scale_kg=5.0))
        assert m == [PurificationMethod.EXTRACTION, PurificationMethod.DISTILLATION]

    def test_addition_small_gives_extraction_flash(self):
        t = _make_template(ReactionCategory.ADDITION)
        m = self._methods(recommend_purification(t, scale_kg=0.5))
        assert m == [PurificationMethod.EXTRACTION, PurificationMethod.FLASH_CHROMATOGRAPHY]

    # -- Deprotection --

    def test_deprotection_gives_extraction_precipitation(self):
        t = _make_template(ReactionCategory.DEPROTECTION)
        m = self._methods(recommend_purification(t, scale_kg=1.0))
        assert m == [PurificationMethod.EXTRACTION, PurificationMethod.PRECIPITATION]

    # -- Polymerization --

    def test_polymerization_gives_precipitation_filtration(self):
        t = _make_template(ReactionCategory.POLYMERIZATION)
        m = self._methods(recommend_purification(t, scale_kg=1.0))
        assert m == [PurificationMethod.PRECIPITATION, PurificationMethod.FILTRATION]

    # -- Catalytic prefix --

    def test_catalytic_reaction_starts_with_filtration(self):
        t = _make_template(ReactionCategory.COUPLING, catalysts=["Pd(PPh3)4"])
        steps = recommend_purification(t, scale_kg=1.0)
        assert steps[0].method == PurificationMethod.FILTRATION
        # Followed by solid-product path (extraction + recrystallization)
        assert steps[1].method == PurificationMethod.EXTRACTION
        assert steps[2].method == PurificationMethod.RECRYSTALLIZATION


# =====================================================================
#  2. Conditions -- atmosphere
# =====================================================================

class TestConditionsAtmosphere:
    """Test _select_atmosphere for inert categories and sensitive reagents.

    NOTE: reagent names must normalize to keys in the sensitive_keywords set
    *without* hitting the alias table (aliases resolve to canonical names
    that are NOT in the keyword set).  Use e.g. "DIBAL" not "LiAlH4".
    """

    def test_coupling_gives_n2(self):
        t = _make_template(ReactionCategory.COUPLING)
        assert _select_atmosphere(t) == "N2"

    def test_reduction_gives_n2(self):
        t = _make_template(ReactionCategory.REDUCTION)
        assert _select_atmosphere(t) == "N2"

    def test_dibal_gives_ar(self):
        t = _make_template(ReactionCategory.SUBSTITUTION, reagents=["DIBAL"])
        assert _select_atmosphere(t) == "Ar"

    def test_memgbr_gives_ar(self):
        t = _make_template(ReactionCategory.SUBSTITUTION, reagents=["MeMgBr"])
        assert _select_atmosphere(t) == "Ar"

    def test_non_sensitive_substitution_gives_air(self):
        t = _make_template(ReactionCategory.SUBSTITUTION, reagents=["benzene"])
        assert _select_atmosphere(t) == "air"


# =====================================================================
#  3. Conditions -- addition rate
# =====================================================================

class TestConditionsAdditionRate:
    """Test _addition_rate for scale and category branches."""

    def test_non_controlled_tiny_all_at_once(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        assert _addition_rate(t, 0.05) == "all at once"

    def test_non_controlled_normal_portionwise(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        assert _addition_rate(t, 1.0) == "portion-wise over 10 min"

    def test_controlled_small_dropwise_short(self):
        t = _make_template(ReactionCategory.ADDITION)
        assert _addition_rate(t, 0.5) == "dropwise over 15-30 min"

    def test_controlled_medium_dropwise_long(self):
        t = _make_template(ReactionCategory.ADDITION)
        result = _addition_rate(t, 5.0)
        assert "dropwise over 30-60 min" in result

    def test_controlled_large_metered_1_2h(self):
        t = _make_template(ReactionCategory.ADDITION)
        result = _addition_rate(t, 50.0)
        assert "metered addition over 1-2 h" in result

    def test_huge_scale_always_metered_2_4h(self):
        # scale > 100 triggers needs_control regardless of category
        t = _make_template(ReactionCategory.SUBSTITUTION)
        result = _addition_rate(t, 200.0)
        assert "metered addition over 2-4 h" in result


# =====================================================================
#  4. Conditions -- workup procedures and solvent selection
# =====================================================================

class TestConditionsWorkupAndSolvent:
    """Test _workup_procedure per category + _select_solvent."""

    def test_oxidation_workup_quench_and_ethyl_acetate(self):
        t = _make_template(ReactionCategory.OXIDATION)
        workup = _workup_procedure(t, 1.0)
        assert "Quench" in workup
        assert "ethyl acetate" in workup

    def test_reduction_workup_quench(self):
        t = _make_template(ReactionCategory.REDUCTION)
        workup = _workup_procedure(t, 1.0)
        assert "Quench" in workup

    def test_coupling_workup_celite(self):
        t = _make_template(ReactionCategory.COUPLING)
        workup = _workup_procedure(t, 1.0)
        assert "Celite" in workup

    def test_substitution_workup_ice_water(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        workup = _workup_procedure(t, 1.0)
        assert "ice water" in workup

    def test_elimination_workup_ice_water(self):
        t = _make_template(ReactionCategory.ELIMINATION)
        workup = _workup_procedure(t, 1.0)
        assert "ice water" in workup

    def test_carbonyl_workup_nh4cl(self):
        t = _make_template(ReactionCategory.CARBONYL)
        workup = _workup_procedure(t, 1.0)
        assert "NH4Cl" in workup

    def test_protection_workup_1m_hcl(self):
        t = _make_template(ReactionCategory.PROTECTION)
        workup = _workup_procedure(t, 1.0)
        assert "1M HCl" in workup

    def test_deprotection_workup_1m_hcl(self):
        t = _make_template(ReactionCategory.DEPROTECTION)
        workup = _workup_procedure(t, 1.0)
        assert "1M HCl" in workup

    def test_polymerization_workup_anti_solvent(self):
        t = _make_template(ReactionCategory.POLYMERIZATION)
        workup = _workup_procedure(t, 1.0)
        assert "anti-solvent" in workup

    def test_pericyclic_workup_standard_aqueous(self):
        t = _make_template(ReactionCategory.PERICYCLIC)
        workup = _workup_procedure(t, 1.0)
        assert "Standard aqueous" in workup

    def test_large_scale_adds_mixer_settler(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        workup = _workup_procedure(t, 100.0)
        assert "mixer-settler" in workup

    def test_solvent_from_template(self):
        t = _make_template(ReactionCategory.COUPLING, solvents=["DMF"])
        assert _select_solvent(t) == "DMF"

    def test_solvent_default_coupling_thf(self):
        t = _make_template(ReactionCategory.COUPLING)
        assert _select_solvent(t) == "THF"

    def test_solvent_default_polymerization_toluene(self):
        t = _make_template(ReactionCategory.POLYMERIZATION)
        assert _select_solvent(t) == "toluene"


# =====================================================================
#  5. Conditions -- scale effects on concentration and pressure
# =====================================================================

class TestConditionsScaleEffects:
    """Test optimize_conditions for scale-dependent adjustments."""

    def test_scale_over_100_dilutes_30_percent(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)  # base conc 0.5
        cond = optimize_conditions(t, scale_kg=200.0)
        assert cond.concentration_M == pytest.approx(0.5 * 0.7, abs=0.01)

    def test_scale_10_to_100_dilutes_15_percent(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        cond = optimize_conditions(t, scale_kg=50.0)
        assert cond.concentration_M == pytest.approx(0.5 * 0.85, abs=0.01)

    def test_scale_under_10_keeps_base(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        cond = optimize_conditions(t, scale_kg=5.0)
        assert cond.concentration_M == pytest.approx(0.5, abs=0.01)

    def test_high_temp_sets_pressure_3_atm(self):
        t = _make_template(
            ReactionCategory.SUBSTITUTION, temperature_range=(160, 180)
        )
        cond = optimize_conditions(t, scale_kg=1.0)
        assert cond.pressure_atm == 3.0

    def test_pd_c_hydrogenation_sets_pressure_4_atm(self):
        t = _make_template(ReactionCategory.REDUCTION, reagents=["H2 Pd C"])
        cond = optimize_conditions(t, scale_kg=1.0)
        assert cond.pressure_atm == 4.0

    def test_lindlar_hydrogenation_sets_pressure_1_atm(self):
        # Lindlar catalyst: poisoned Pd, balloon pressure (1 atm)
        t = _make_template(ReactionCategory.REDUCTION, reagents=["Lindlar"])
        cond = optimize_conditions(t, scale_kg=1.0)
        assert cond.pressure_atm == 1.0

    def test_normal_conditions_pressure_1_atm(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        cond = optimize_conditions(t, scale_kg=1.0)
        assert cond.pressure_atm == 1.0


# =====================================================================
#  6. Safety -- PPE
# =====================================================================

class TestSafetyPPE:
    """Test _determine_ppe with controlled hazard/pictogram sets."""

    def test_baseline_empty_hazards(self):
        ppe = _determine_ppe(set(), set())
        assert "Safety goggles (splash-proof)" in ppe
        assert "Lab coat" in ppe
        assert "Closed-toe shoes" in ppe
        assert any("Nitrile" in p for p in ppe)

    def test_corrosive_ghs05_adds_face_shield_and_apron(self):
        ppe = _determine_ppe(set(), {"GHS05"})
        assert any("Face shield" in p for p in ppe)
        assert any("apron" in p.lower() for p in ppe)
        assert any("Chemical-resistant gloves" in p for p in ppe)

    def test_inhalation_h330_adds_respiratory(self):
        ppe = _determine_ppe({"H330"}, set())
        assert any("respiratory" in p.lower() for p in ppe)

    def test_acute_toxicity_ghs06_adds_eyewash(self):
        ppe = _determine_ppe(set(), {"GHS06"})
        assert any("eyewash" in p.lower() for p in ppe)

    def test_pyrophoric_h250_adds_fire_resistant(self):
        ppe = _determine_ppe({"H250"}, set())
        assert any("Fire-resistant" in p or "Nomex" in p for p in ppe)


# =====================================================================
#  7. Safety -- engineering controls
# =====================================================================

class TestSafetyEngineering:
    """Test _determine_engineering_controls with controlled inputs."""

    def test_baseline_fume_hood(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        controls = _determine_engineering_controls(set(), set(), t)
        assert any("fume hood" in c.lower() for c in controls)

    def test_ghs02_non_sparking_and_grounding(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        controls = _determine_engineering_controls(set(), {"GHS02"}, t)
        joined = " ".join(controls)
        assert "non-sparking" in joined
        assert "Ground" in joined

    def test_h250_schlenk_or_glovebox(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        controls = _determine_engineering_controls({"H250"}, set(), t)
        joined = " ".join(controls)
        assert "Schlenk" in joined or "glovebox" in joined

    def test_ghs03_oxidizer_segregation(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        controls = _determine_engineering_controls(set(), {"GHS03"}, t)
        joined = " ".join(controls)
        assert "oxidizer" in joined.lower() or "separated" in joined.lower()

    def test_cmr_h340_designated_area_and_hepa(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        controls = _determine_engineering_controls({"H340"}, set(), t)
        joined = " ".join(controls)
        assert "CMR" in joined or "Designated" in joined
        assert "HEPA" in joined

    def test_ghs09_secondary_containment(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        controls = _determine_engineering_controls(set(), {"GHS09"}, t)
        joined = " ".join(controls)
        assert "Secondary containment" in joined

    def test_cryogenic_cooling_note(self):
        t = _make_template(
            ReactionCategory.SUBSTITUTION, temperature_range=(-78, -70)
        )
        controls = _determine_engineering_controls(set(), set(), t)
        joined = " ".join(controls)
        assert "Cryogenic" in joined or "cryogenic" in joined

    def test_high_temp_thermal_insulation(self):
        t = _make_template(
            ReactionCategory.SUBSTITUTION, temperature_range=(160, 180)
        )
        controls = _determine_engineering_controls(set(), set(), t)
        joined = " ".join(controls)
        assert "High-temperature" in joined or "thermal" in joined.lower()


# =====================================================================
#  8. Safety -- emergency procedures, waste, risk level
# =====================================================================

class TestSafetyEmergencyAndWaste:
    """Test emergency procedures, waste classification, and risk level."""

    # -- Emergency procedures --

    def test_ghs06_ingestion_no_vomiting(self):
        procs = _determine_emergency_procedures(set(), {"GHS06"})
        joined = " ".join(procs)
        assert "Do NOT induce vomiting" in joined

    def test_h300_ingestion_no_vomiting(self):
        procs = _determine_emergency_procedures({"H300"}, set())
        joined = " ".join(procs)
        assert "Do NOT induce vomiting" in joined

    def test_h310_skin_remove_clothing(self):
        procs = _determine_emergency_procedures({"H310"}, set())
        joined = " ".join(procs)
        assert "Remove contaminated clothing" in joined

    def test_h311_skin_remove_clothing(self):
        procs = _determine_emergency_procedures({"H311"}, set())
        joined = " ".join(procs)
        assert "Remove contaminated clothing" in joined

    def test_h330_inhalation_fresh_air(self):
        procs = _determine_emergency_procedures({"H330"}, set())
        joined = " ".join(procs)
        assert "fresh air" in joined

    def test_h314_eye_flush(self):
        procs = _determine_emergency_procedures({"H314"}, set())
        joined = " ".join(procs)
        assert "Flush eyes" in joined

    def test_h318_eye_flush(self):
        procs = _determine_emergency_procedures({"H318"}, set())
        joined = " ".join(procs)
        assert "Flush eyes" in joined

    def test_ghs02_fire_dry_chemical_co2(self):
        procs = _determine_emergency_procedures(set(), {"GHS02"})
        joined = " ".join(procs)
        assert "dry chemical" in joined or "CO2" in joined

    def test_h260_water_reactive_dry_sand(self):
        procs = _determine_emergency_procedures({"H260"}, set())
        joined = " ".join(procs)
        assert "dry" in joined.lower()
        assert "sand" in joined.lower()

    def test_baseline_always_has_spill_procedure(self):
        procs = _determine_emergency_procedures(set(), set())
        joined = " ".join(procs)
        assert "spill" in joined.lower()

    # -- Waste classification --

    def test_h300_waste_category_1(self):
        assert "Category 1" in _classify_waste({"H300"})

    def test_h340_waste_category_1(self):
        assert "Category 1" in _classify_waste({"H340"})

    def test_h301_waste_category_2(self):
        assert "Category 2" in _classify_waste({"H301"})

    def test_h314_waste_category_2(self):
        assert "Category 2" in _classify_waste({"H314"})

    def test_h225_waste_flammable(self):
        assert "flammable" in _classify_waste({"H225"})

    def test_h224_waste_flammable(self):
        assert "flammable" in _classify_waste({"H224"})

    def test_empty_waste_non_hazardous(self):
        assert "Non-hazardous" in _classify_waste(set())

    # -- Risk level --

    def test_h300_risk_high(self):
        assert _calculate_risk_level({"H300"}, set()) == "high"

    def test_h250_risk_high(self):
        assert _calculate_risk_level({"H250"}, set()) == "high"

    def test_h314_risk_medium(self):
        assert _calculate_risk_level({"H314"}, set()) == "medium"

    def test_h225_risk_medium(self):
        assert _calculate_risk_level({"H225"}, set()) == "medium"

    def test_empty_risk_low(self):
        assert _calculate_risk_level(set(), set()) == "low"


# =====================================================================
#  9. Safety -- incompatibilities
# =====================================================================

class TestSafetyIncompatibilities:
    """Test _determine_incompatibilities with mock templates.

    NOTE: reagent names must normalize to keys in the internal keyword
    sets WITHOUT hitting the alias table.  E.g. use "CrO3" (not
    "KMnO4" which aliases to "potassium_permanganate" and misses the
    "kmno4" keyword set entry).
    """

    def test_oxidizer_plus_reducer_critical(self):
        # CrO3 -> "cro3" (in oxidizers), DIBAL_H -> "dibal_h" (in reducers)
        t = _make_template(ReactionCategory.MISC, reagents=["CrO3", "DIBAL_H"])
        incompat = _determine_incompatibilities(t)
        joined = " ".join(str(w) for w in incompat)
        assert "CRITICAL" in joined
        assert "Oxidizer" in joined

    def test_acid_plus_cyanide_hcn(self):
        # HNO3 -> "hno3" (in acids), NaCN -> "nacn" (in cyanides)
        t = _make_template(ReactionCategory.MISC, reagents=["HNO3", "NaCN"])
        incompat = _determine_incompatibilities(t)
        joined = " ".join(str(w) for w in incompat)
        assert "HCN" in joined

    def test_acid_plus_azide_hn3(self):
        # HNO3 -> "hno3" (in acids), NaN3 -> "nan3" (in azides)
        t = _make_template(ReactionCategory.MISC, reagents=["HNO3", "NaN3"])
        incompat = _determine_incompatibilities(t)
        joined = " ".join(str(w) for w in incompat)
        assert "HN3" in joined

    def test_water_reactive_anhydrous_warning(self):
        # SOCl2 -> "socl2" (in water_reactive)
        t = _make_template(ReactionCategory.MISC, reagents=["SOCl2"])
        incompat = _determine_incompatibilities(t)
        joined = " ".join(str(w) for w in incompat)
        assert "anhydrous" in joined

    def test_permanganate_plus_organic_fire(self):
        # NaMnO4 -> "namno4" (in permanganates), acetone -> "acetone" (in flammable_organics)
        t = _make_template(ReactionCategory.MISC, reagents=["NaMnO4", "acetone"])
        incompat = _determine_incompatibilities(t)
        joined = " ".join(str(w) for w in incompat)
        assert "fire" in joined.lower() or "explosion" in joined.lower()

    def test_acid_plus_base_exothermic(self):
        # HNO3 -> "hno3" (in conc_acids), NaOH_aq -> "naoh_aq" (in bases)
        t = _make_template(ReactionCategory.MISC, reagents=["HNO3", "NaOH_aq"])
        incompat = _determine_incompatibilities(t)
        joined = " ".join(str(w) for w in incompat)
        assert "exothermic" in joined

    def test_no_incompatible_empty(self):
        # "benzene" normalizes to "benzene" which is not in any keyword set
        t = _make_template(ReactionCategory.MISC, reagents=["benzene"])
        assert _determine_incompatibilities(t) == []

    def test_organometallic_plus_protic(self):
        # MeMgBr -> "memgbr" (in organometallics), MeOH -> "meoh" (in protic_solvents)
        t = _make_template(ReactionCategory.MISC, reagents=["MeMgBr", "MeOH"])
        incompat = _determine_incompatibilities(t)
        joined = " ".join(str(w) for w in incompat)
        assert "Organometallic" in joined or "protic" in joined


# =====================================================================
#  10. Scale-up modes
# =====================================================================

class TestScaleUpModes:
    """Test mode selection, cycle time, risks, and recommendations."""

    # -- Mode selection --

    def test_small_target_gives_batch(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        result = analyze_scale_up([_MockStep(t)], target_annual_kg=500.0)
        assert result.recommended_mode == "batch"

    def test_large_fast_gives_continuous(self):
        t = _make_template(
            ReactionCategory.SUBSTITUTION, typical_yield=(80, 95)
        )
        result = analyze_scale_up([_MockStep(t)], target_annual_kg=100_000.0)
        assert result.recommended_mode == "continuous"

    def test_mid_range_few_steps_gives_semi_continuous(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        result = analyze_scale_up([_MockStep(t)], target_annual_kg=15_000.0)
        assert result.recommended_mode == "semi-continuous"

    def test_empty_steps_gives_zeroed_analysis(self):
        result = analyze_scale_up([], target_annual_kg=1000.0)
        assert result.cycle_time_hours == 0.0
        assert result.annual_capacity_kg == 0.0
        assert result.capital_cost_usd == 0.0

    # -- Cycle time --

    def test_cycle_time_cryogenic_slower(self):
        cryo = _make_template(
            ReactionCategory.SUBSTITUTION, temperature_range=(-78, -60)
        )
        normal = _make_template(
            ReactionCategory.SUBSTITUTION, temperature_range=(20, 25)
        )
        ct_cryo = _estimate_cycle_time([_MockStep(cryo)])
        ct_normal = _estimate_cycle_time([_MockStep(normal)])
        assert ct_cryo > ct_normal

    def test_cycle_time_high_temp_faster(self):
        hot = _make_template(
            ReactionCategory.SUBSTITUTION, temperature_range=(130, 150)
        )
        normal = _make_template(
            ReactionCategory.SUBSTITUTION, temperature_range=(20, 25)
        )
        ct_hot = _estimate_cycle_time([_MockStep(hot)])
        ct_normal = _estimate_cycle_time([_MockStep(normal)])
        assert ct_hot < ct_normal

    def test_cycle_time_polymerization_longer_workup(self):
        poly = _make_template(ReactionCategory.POLYMERIZATION)
        sub = _make_template(
            ReactionCategory.SUBSTITUTION,
            typical_yield=poly.typical_yield,
            temperature_range=poly.temperature_range,
        )
        ct_poly = _estimate_cycle_time([_MockStep(poly)])
        ct_sub = _estimate_cycle_time([_MockStep(sub)])
        # POLYMERIZATION workup = 3.0h vs default 1.5h
        assert ct_poly > ct_sub

    def test_cycle_time_coupling_longer_workup(self):
        coup = _make_template(ReactionCategory.COUPLING)
        sub = _make_template(
            ReactionCategory.SUBSTITUTION,
            typical_yield=coup.typical_yield,
            temperature_range=coup.temperature_range,
        )
        ct_coup = _estimate_cycle_time([_MockStep(coup)])
        ct_sub = _estimate_cycle_time([_MockStep(sub)])
        # COUPLING workup = 2.0h vs default 1.5h
        assert ct_coup > ct_sub

    # -- Continuous candidacy --

    def test_polymerization_not_continuous_candidate(self):
        t = _make_template(ReactionCategory.POLYMERIZATION)
        assert _is_continuous_candidate([_MockStep(t)], 100_000.0) is False

    def test_below_50k_not_continuous_candidate(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        assert _is_continuous_candidate([_MockStep(t)], 40_000.0) is False

    # -- Risks --

    def test_extreme_temp_heat_transfer_risk(self):
        t = _make_template(
            ReactionCategory.SUBSTITUTION, temperature_range=(160, 180)
        )
        risks = _identify_risks([_MockStep(t)], 1000.0, "batch")
        joined = " ".join(risks)
        assert "heat transfer" in joined.lower() or "Extreme temperature" in joined

    def test_addition_exothermic_risk(self):
        t = _make_template(ReactionCategory.ADDITION)
        risks = _identify_risks([_MockStep(t)], 1000.0, "batch")
        joined = " ".join(risks)
        assert "exothermic" in joined.lower()

    def test_regulatory_risk_over_10k(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        risks = _identify_risks([_MockStep(t)], 15_000.0, "batch")
        joined = " ".join(risks)
        assert "Regulatory" in joined or "10 t/year" in joined

    def test_continuous_mode_pat_risk(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        risks = _identify_risks([_MockStep(t)], 100_000.0, "continuous")
        joined = " ".join(risks)
        assert "PAT" in joined or "analytical technology" in joined

    def test_safety_notes_propagated_to_risk(self):
        t = _make_template(
            ReactionCategory.SUBSTITUTION, safety_notes="Toxic fumes"
        )
        risks = _identify_risks([_MockStep(t)], 1000.0, "batch")
        joined = " ".join(risks)
        assert "Toxic fumes" in joined

    # -- Recommendations --

    def test_batch_recommendations_include_calorimetry(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        recs = _generate_recommendations([_MockStep(t)], 1000.0, "batch", 100.0)
        joined = " ".join(recs)
        assert "calorimetry" in joined.lower()

    def test_continuous_recommendations_include_rtd_and_pat(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        recs = _generate_recommendations([_MockStep(t)], 100_000.0, "continuous", None)
        joined = " ".join(recs)
        assert "RTD" in joined or "residence time" in joined.lower()
        assert "PAT" in joined or "process analytical" in joined.lower()

    def test_semi_continuous_recommendations_identify_steps(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        recs = _generate_recommendations(
            [_MockStep(t)], 15_000.0, "semi-continuous", 500.0
        )
        joined = " ".join(recs)
        assert "continuous operation" in joined.lower() or "batchwise" in joined

    def test_large_batch_recommends_inline_analytics(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        recs = _generate_recommendations([_MockStep(t)], 5000.0, "batch", 600.0)
        joined = " ".join(recs)
        assert "in-line" in joined.lower() or "FTIR" in joined or "Raman" in joined

    def test_high_volume_recommends_catalyst_recycling(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        recs = _generate_recommendations([_MockStep(t)], 5000.0, "batch", 500.0)
        joined = " ".join(recs)
        assert "catalyst recycling" in joined.lower() or "solvent recovery" in joined.lower()


# =====================================================================
#  11. Reactor branches
# =====================================================================

from molbuilder.process.reactor import (
    ReactorType,
    select_reactor,
    _is_cryogenic,
    _select_material,
    _estimate_reactor_cost,
)


class TestReactorBranches:
    """Cover internal helpers and decision-tree branches in reactor.py."""

    # -- Internal helpers --

    def test_is_cryogenic_true(self):
        t = _make_template(ReactionCategory.SUBSTITUTION, temperature_range=(-78, -60))
        assert _is_cryogenic(t) is True

    def test_is_cryogenic_false(self):
        t = _make_template(ReactionCategory.SUBSTITUTION, temperature_range=(20, 25))
        assert _is_cryogenic(t) is False

    def test_select_material_corrosive_hastelloy(self):
        # HNO3 normalizes to "hno3" (no alias), which is in corrosive_keywords
        t = _make_template(ReactionCategory.SUBSTITUTION, reagents=["HNO3"])
        assert _select_material(t) == "hastelloy"

    def test_select_material_high_temp_stainless(self):
        t = _make_template(ReactionCategory.SUBSTITUTION, temperature_range=(160, 180))
        assert _select_material(t) == "stainless steel"

    def test_select_material_default_glass(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        assert _select_material(t) == "glass"

    def test_estimate_reactor_cost_zero_volume(self):
        # volume_L <= 0 guard sets volume to 1.0
        cost = _estimate_reactor_cost(ReactorType.BATCH, 0.0)
        assert cost > 0

    # -- Decision tree: Branch 1 (MICROREACTOR) --

    def test_microreactor_high_temp_pressure_2(self):
        # Exothermic + small scale + mean_t > 100 -> MICROREACTOR with pressure 2.0
        t = _make_template(ReactionCategory.ADDITION, temperature_range=(110, 130))
        spec = select_reactor(t, scale_kg=0.5)
        assert spec.reactor_type == ReactorType.MICROREACTOR
        assert spec.pressure_atm == 2.0

    # -- Decision tree: Branch 2 (CSTR) --

    def test_cstr_normal_temp_pressure_1(self):
        # Exothermic + large scale + mean_t < 100 -> CSTR, pressure 1.0
        t = _make_template(ReactionCategory.ADDITION, temperature_range=(20, 30))
        spec = select_reactor(t, scale_kg=5.0)
        assert spec.reactor_type == ReactorType.CSTR
        assert spec.pressure_atm == 1.0

    def test_cstr_high_temp_pressure_3(self):
        # Exothermic + large scale + mean_t > 100 -> CSTR, pressure 3.0
        t = _make_template(ReactionCategory.RADICAL, temperature_range=(110, 130))
        spec = select_reactor(t, scale_kg=5.0)
        assert spec.reactor_type == ReactorType.CSTR
        assert spec.pressure_atm == 3.0

    # -- Decision tree: Branch 3a (FIXED_BED) --

    def test_fixed_bed_large_catalytic(self):
        # Multiphase + catalyst + scale > 100 -> FIXED_BED
        t = _make_template(
            ReactionCategory.COUPLING, catalysts=["Pd(PPh3)4"],
        )
        spec = select_reactor(t, scale_kg=200.0)
        assert spec.reactor_type == ReactorType.FIXED_BED
        assert spec.pressure_atm == 3.0

    # -- Decision tree: Branch 3b (BATCH multiphase) --

    def test_batch_multiphase_no_catalyst(self):
        # Multiphase + no catalyst + scale <= 100 -> BATCH
        t = _make_template(ReactionCategory.COUPLING)
        spec = select_reactor(t, scale_kg=10.0)
        assert spec.reactor_type == ReactorType.BATCH
        assert spec.mixing_type == "mechanical"

    # -- Decision tree: Branch 5 (SEMI_BATCH) --

    def test_semi_batch_carbonyl(self):
        # CARBONYL + scale > 10 -> SEMI_BATCH
        t = _make_template(ReactionCategory.CARBONYL)
        spec = select_reactor(t, scale_kg=20.0)
        assert spec.reactor_type == ReactorType.SEMI_BATCH

    def test_semi_batch_substitution(self):
        # SUBSTITUTION + scale > 10 -> SEMI_BATCH
        t = _make_template(ReactionCategory.SUBSTITUTION)
        spec = select_reactor(t, scale_kg=15.0)
        assert spec.reactor_type == ReactorType.SEMI_BATCH

    # -- Decision tree: Branch 6 (default BATCH) --

    def test_default_batch_small_adiabatic(self):
        # Non-exothermic, non-multiphase, small scale -> BATCH, adiabatic
        # scale_kg=1 -> volume = round(1*5/0.75, 1) = 6.7 <= 20
        t = _make_template(ReactionCategory.PERICYCLIC)
        spec = select_reactor(t, scale_kg=1.0)
        assert spec.reactor_type == ReactorType.BATCH
        assert spec.heat_transfer == "adiabatic"

    def test_default_batch_large_jacketed(self):
        # Non-exothermic, non-multiphase, moderate scale -> BATCH, jacketed
        # scale_kg=5 -> volume = round(5*5/0.75, 1) = 33.3 > 20
        t = _make_template(ReactionCategory.PERICYCLIC)
        spec = select_reactor(t, scale_kg=5.0)
        assert spec.reactor_type == ReactorType.BATCH
        assert spec.heat_transfer == "jacketed"

    # -- TypeError guard --

    def test_template_missing_temperature_range_raises(self):
        class BadTemplate:
            pass
        with pytest.raises(TypeError, match="temperature_range"):
            select_reactor(BadTemplate(), scale_kg=1.0)


# =====================================================================
#  12. Costing branches
# =====================================================================

from molbuilder.process.costing import (
    _normalise_key,
    _labor_cost_for_step,
    _waste_cost_for_step,
    estimate_cost,
    CostParameters,
)


class TestCostingBranches:
    """Cover labor/waste/notes branches in costing.py."""

    # -- _normalise_key direct call (line 80) --

    def test_normalise_key_lowercases(self):
        assert _normalise_key("HNO3") == "hno3"

    # -- Labor cost branches --

    def test_labor_baseline(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        cost = _labor_cost_for_step(t, scale_kg=1.0, params=CostParameters())
        # Base 3h * $75/h = $225
        assert cost == pytest.approx(225.0)

    def test_labor_cryogenic_adds_hour(self):
        t = _make_template(ReactionCategory.SUBSTITUTION, temperature_range=(-78, -60))
        cost = _labor_cost_for_step(t, scale_kg=1.0, params=CostParameters())
        # 3h + 1h cryo = 4h * $75 = $300
        assert cost == pytest.approx(300.0)

    def test_labor_rearrangement_adds_hour(self):
        t = _make_template(ReactionCategory.REARRANGEMENT)
        cost = _labor_cost_for_step(t, scale_kg=1.0, params=CostParameters())
        # 3h + 1h complexity = 4h * $75 = $300
        assert cost == pytest.approx(300.0)

    def test_labor_misc_adds_hour(self):
        t = _make_template(ReactionCategory.MISC)
        cost = _labor_cost_for_step(t, scale_kg=1.0, params=CostParameters())
        assert cost == pytest.approx(300.0)

    def test_labor_scale_over_50_adds_hour(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        cost = _labor_cost_for_step(t, scale_kg=60.0, params=CostParameters())
        # 3h + 1h material handling = 4h * $75 = $300
        assert cost == pytest.approx(300.0)

    def test_labor_scale_over_200_adds_two_hours(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        cost = _labor_cost_for_step(t, scale_kg=250.0, params=CostParameters())
        # 3h + 1h (>50) + 1h (>200) = 5h * $75 = $375
        assert cost == pytest.approx(375.0)

    # -- Waste cost branches --

    def test_waste_hazardous_nacn(self):
        t = _make_template(ReactionCategory.SUBSTITUTION, reagents=["NaCN"])
        cost = _waste_cost_for_step(t, scale_kg=1.0, params=CostParameters())
        # hazardous: waste_kg = 1 * 8 = 8, disposal = $2.50 * 2 = $5.00
        assert cost == pytest.approx(8.0 * 5.0)

    def test_waste_non_hazardous(self):
        t = _make_template(ReactionCategory.SUBSTITUTION, reagents=["benzene"])
        cost = _waste_cost_for_step(t, scale_kg=1.0, params=CostParameters())
        # non-hazardous: waste_kg = 1 * 5 = 5, disposal = $2.50
        assert cost == pytest.approx(5.0 * 2.50)

    # -- estimate_cost edge cases --

    def test_estimate_cost_scale_zero_guard(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        result = estimate_cost([_MockStep(t)], scale_kg=0.0)
        # scale forced to 0.001, should not raise
        assert result.scale_kg == 0.001
        assert result.total_usd > 0

    def test_estimate_cost_zero_yield_fallback(self):
        t = _make_template(ReactionCategory.SUBSTITUTION, typical_yield=(0, 0))
        result = estimate_cost([_MockStep(t)], scale_kg=1.0)
        # avg_yield forced to 0.5, should not error
        assert result.total_usd > 0

    def test_estimate_cost_lab_scale_note(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        result = estimate_cost([_MockStep(t)], scale_kg=0.05)
        assert any("Lab-scale" in n for n in result.notes)

    def test_estimate_cost_bulk_scale_note(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        result = estimate_cost([_MockStep(t)], scale_kg=200.0)
        assert any("Bulk pricing" in n for n in result.notes)

    def test_estimate_cost_expensive_reagent_note(self):
        t = _make_template(ReactionCategory.SUBSTITUTION, reagents=["OsO4"])
        result = estimate_cost([_MockStep(t)], scale_kg=1.0)
        assert any("expensive reagent" in n.lower() for n in result.notes)

    def test_estimate_cost_step_without_template_raises(self):
        class BadStep:
            pass
        with pytest.raises(TypeError, match="template"):
            estimate_cost([BadStep()], scale_kg=1.0)


# =====================================================================
#  13. CostParameters
# =====================================================================

class TestCostParameters:
    """Test configurable costing model."""

    def test_default_matches_hardcoded(self):
        p = CostParameters()
        assert p.labor_rate_usd_h == 75.0
        assert p.base_labor_hours_per_step == 3.0
        assert p.energy_cost_kwh == 0.10
        assert p.waste_disposal_per_kg == 2.50
        assert p.overhead_fraction == 0.25
        assert p.equipment_depreciation_per_batch == 0.002
        assert p.default_reagent_equiv == 1.2
        assert p.solvent_l_per_kg == 7.0

    def test_us_default_is_identical(self):
        us = CostParameters.us_default()
        default = CostParameters()
        assert us == default

    def test_eu_different_from_us(self):
        eu = CostParameters.eu_default()
        us = CostParameters.us_default()
        assert eu.labor_rate_usd_h > us.labor_rate_usd_h
        assert eu.waste_disposal_per_kg > us.waste_disposal_per_kg

    def test_india_lower_labor(self):
        india = CostParameters.india_default()
        us = CostParameters.us_default()
        assert india.labor_rate_usd_h < us.labor_rate_usd_h

    def test_china_lower_labor(self):
        china = CostParameters.china_default()
        us = CostParameters.us_default()
        assert china.labor_rate_usd_h < us.labor_rate_usd_h

    def test_region_presets_produce_different_totals(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        step = _MockStep(t)
        us_cost = estimate_cost([step], scale_kg=10.0, params=CostParameters.us_default())
        eu_cost = estimate_cost([step], scale_kg=10.0, params=CostParameters.eu_default())
        assert us_cost.total_usd != eu_cost.total_usd

    def test_custom_params_override(self):
        p = CostParameters(labor_rate_usd_h=200.0)
        t = _make_template(ReactionCategory.SUBSTITUTION)
        step = _MockStep(t)
        expensive = estimate_cost([step], scale_kg=1.0, params=p)
        default = estimate_cost([step], scale_kg=1.0)
        assert expensive.breakdown.labor_usd > default.breakdown.labor_usd

    def test_params_none_uses_default(self):
        t = _make_template(ReactionCategory.SUBSTITUTION)
        step = _MockStep(t)
        r1 = estimate_cost([step], scale_kg=1.0, params=None)
        r2 = estimate_cost([step], scale_kg=1.0)
        assert r1.total_usd == r2.total_usd


# =====================================================================
#  14. New safety incompatibility pairs
# =====================================================================

from molbuilder.process.safety import IncompatibilityWarning  # noqa: E402


class TestNewIncompatibilities:
    """Test the 10 new incompatibility pairs added in Tier 1."""

    def test_azide_plus_heavy_metal(self):
        t = _make_template(ReactionCategory.MISC, reagents=["NaN3", "AgNO3"])
        incompat = _determine_incompatibilities(t)
        joined = " ".join(str(w) for w in incompat)
        assert "detonation" in joined.lower() or "metal azide" in joined.lower()

    def test_nitric_acid_plus_alcohol(self):
        t = _make_template(ReactionCategory.MISC, reagents=["HNO3", "MeOH"])
        incompat = _determine_incompatibilities(t)
        joined = " ".join(str(w) for w in incompat)
        assert "nitrate" in joined.lower()

    def test_chloroform_plus_strong_base(self):
        t = _make_template(ReactionCategory.MISC, reagents=["CHCl3", "NaOH"])
        incompat = _determine_incompatibilities(t)
        joined = " ".join(str(w) for w in incompat)
        assert "phosgene" in joined.lower()

    def test_dmso_plus_acyl_halide(self):
        t = _make_template(ReactionCategory.MISC, reagents=["DMSO", "SOCl2"])
        incompat = _determine_incompatibilities(t)
        joined = " ".join(str(w) for w in incompat)
        assert "DMSO" in joined or "decomposition" in joined.lower()

    def test_h2o2_plus_acetone(self):
        t = _make_template(ReactionCategory.MISC, reagents=["H2O2", "acetone"])
        incompat = _determine_incompatibilities(t)
        joined = " ".join(str(w) for w in incompat)
        assert "TATP" in joined or "explosive" in joined.lower()

    def test_ether_plus_strong_oxidizer(self):
        t = _make_template(ReactionCategory.MISC, reagents=["diethyl_ether", "CrO3"])
        incompat = _determine_incompatibilities(t)
        joined = " ".join(str(w) for w in incompat)
        assert "peroxide" in joined.lower() or "ether" in joined.lower()

    def test_severity_on_warning(self):
        t = _make_template(ReactionCategory.MISC, reagents=["CrO3", "DIBAL_H"])
        incompat = _determine_incompatibilities(t)
        assert len(incompat) > 0
        assert all(isinstance(w, IncompatibilityWarning) for w in incompat)
        assert all(w.severity in ("critical", "high", "medium") for w in incompat)

    def test_mitigation_text_present(self):
        t = _make_template(ReactionCategory.MISC, reagents=["HNO3", "NaCN"])
        incompat = _determine_incompatibilities(t)
        for w in incompat:
            assert w.mitigation, f"Missing mitigation on {w}"

    def test_str_compatibility(self):
        """IncompatibilityWarning.__str__ should contain the hazard text."""
        w = IncompatibilityWarning(
            reagent_a="A", reagent_b="B",
            hazard="Test hazard", severity="critical",
            mitigation="Fix it",
        )
        assert "Test hazard" in str(w)
        assert "CRITICAL" in str(w)
