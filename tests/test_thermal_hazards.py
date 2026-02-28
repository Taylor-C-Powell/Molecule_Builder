"""Tests for thermal hazard detection in safety assessments."""

import pytest

from molbuilder.process.safety import (
    ThermalHazard,
    _assess_thermal_hazards,
    _determine_incompatibilities,
    assess_safety,
)
from molbuilder.reactions.reaction_types import ReactionCategory, ReactionTemplate
from molbuilder.reports.safety_report import generate_safety_report


def _make_template(name="Test reaction", named_reaction=None, category=ReactionCategory.MISC,
                   reagents=None, catalysts=None, temperature_range=(-10, 30)):
    """Helper to create a minimal ReactionTemplate for testing."""
    return ReactionTemplate(
        name=name,
        named_reaction=named_reaction,
        category=category,
        reagents=reagents or [],
        catalysts=catalysts or [],
        solvents=["THF"],
        temperature_range=temperature_range,
        typical_yield=(50, 80),
    )


# =====================================================================
#  Thermal hazard detection
# =====================================================================

class TestThermalHazards:

    def test_grignard_thermal_hazard(self):
        """Grignard reagent triggers critical thermal hazard."""
        tmpl = _make_template(
            named_reaction="Grignard reaction",
            reagents=["EtMgBr"],
        )
        hazards = _assess_thermal_hazards(tmpl)
        assert len(hazards) >= 1
        grignard = [h for h in hazards if "Grignard" in h.reaction_type]
        assert len(grignard) == 1
        assert grignard[0].severity == "critical"
        assert grignard[0].max_temp_c == 35.0

    def test_hydrogenation_thermal_hazard(self):
        """Catalytic hydrogenation triggers high thermal hazard."""
        tmpl = _make_template(
            named_reaction="Catalytic hydrogenation",
            reagents=["H2_Pd_C"],
        )
        hazards = _assess_thermal_hazards(tmpl)
        hydro = [h for h in hazards if "hydrogenation" in h.reaction_type.lower()]
        assert len(hydro) == 1
        assert hydro[0].severity == "high"
        assert hydro[0].max_temp_c == 50.0

    def test_diazotization_thermal_hazard(self):
        """Diazotization triggers critical thermal hazard."""
        tmpl = _make_template(
            named_reaction="Diazotization",
            reagents=["NaNO2", "HCl"],
            temperature_range=(-5, 5),
        )
        hazards = _assess_thermal_hazards(tmpl)
        diazo = [h for h in hazards if "Diazotization" in h.reaction_type]
        assert len(diazo) == 1
        assert diazo[0].severity == "critical"
        assert diazo[0].max_temp_c == 5.0

    def test_nitration_thermal_hazard(self):
        """Mixed acid nitration triggers critical thermal hazard."""
        tmpl = _make_template(
            named_reaction="Nitration",
            reagents=["HNO3", "H2SO4"],
            temperature_range=(0, 10),
        )
        hazards = _assess_thermal_hazards(tmpl)
        nitro = [h for h in hazards if "Nitration" in h.reaction_type]
        assert len(nitro) == 1
        assert nitro[0].severity == "critical"

    def test_organolithium_thermal_hazard(self):
        """n-BuLi triggers critical thermal hazard."""
        tmpl = _make_template(reagents=["n-BuLi"], temperature_range=(-78, -60))
        hazards = _assess_thermal_hazards(tmpl)
        li = [h for h in hazards if "lithium" in h.reaction_type.lower()]
        assert len(li) == 1
        assert li[0].severity == "critical"
        assert li[0].max_temp_c == -78.0

    def test_no_thermal_hazard_for_mild_reaction(self):
        """A simple coupling with mild reagents has no thermal hazards."""
        tmpl = _make_template(
            named_reaction="Suzuki coupling",
            reagents=["K2CO3"],
            catalysts=["Pd(PPh3)4"],
            category=ReactionCategory.COUPLING,
        )
        hazards = _assess_thermal_hazards(tmpl)
        assert len(hazards) == 0

    def test_thermal_hazard_in_safety_report(self):
        """End-to-end: thermal hazards appear in the ASCII safety report."""

        class FakeStep:
            def __init__(self):
                self.template = _make_template(
                    named_reaction="Grignard reaction",
                    reagents=["EtMgBr"],
                )

        assessments = assess_safety([FakeStep()])
        assert len(assessments) == 1
        assert len(assessments[0].thermal_hazards) >= 1

        report = generate_safety_report(assessments)
        assert "Thermal Hazards" in report
        assert "Grignard" in report
        assert "CRITICAL" in report


# =====================================================================
#  Solvent-reagent incompatibility checks
# =====================================================================

class TestSolventReagentIncompatibility:

    def test_nah_in_water_critical(self):
        """NaH in water/protic solvent should be critical."""
        tmpl = _make_template(reagents=["NaH"])
        warnings = _determine_incompatibilities(tmpl, solvent="water")
        solvent_warns = [w for w in warnings if "protic solvent" in w.reagent_a.lower()
                         or "water-reactive" in w.reagent_a.lower()]
        assert any(w.severity == "critical" for w in solvent_warns)

    def test_nbuli_in_thf_ok(self):
        """n-BuLi in THF should NOT trigger solvent-reagent incompatibility.
        (THF is an ethereal solvent, not protic.)"""
        tmpl = _make_template(reagents=["n-BuLi"])
        warnings = _determine_incompatibilities(tmpl, solvent="THF")
        protic_warns = [w for w in warnings if "protic solvent" in str(w)]
        assert len(protic_warns) == 0

    def test_dmso_high_temp_warning(self):
        """DMSO as solvent at high temperature should warn."""
        tmpl = _make_template(
            reagents=["K2CO3"],
            temperature_range=(160, 200),
        )
        warnings = _determine_incompatibilities(tmpl, solvent="DMSO")
        dmso_warns = [w for w in warnings if "DMSO" in w.reagent_a and "high temperature" in w.reagent_b]
        assert len(dmso_warns) == 1
        assert dmso_warns[0].severity == "high"

    def test_no_solvent_no_crash(self):
        """Calling without solvent still works (backwards-compatible)."""
        tmpl = _make_template(reagents=["NaH"])
        warnings = _determine_incompatibilities(tmpl)
        # Should not crash; NaH still triggers water-reactive warning
        assert isinstance(warnings, list)
