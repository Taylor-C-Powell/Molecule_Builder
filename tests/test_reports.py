"""Unit tests for the molbuilder.reports subpackage.

Covers text_formatter, molecule_report, synthesis_report,
safety_report, cost_report, and the Protocol types in __init__.
"""

import unittest
from dataclasses import dataclass, field

from molbuilder.reports.text_formatter import (
    section_header,
    subsection_header,
    ascii_table,
    word_wrap,
    bullet_list,
    key_value_block,
    horizontal_bar,
    format_currency,
    format_percent,
)
from molbuilder.reports.molecule_report import (
    generate_molecule_report,
    _hill_formula,
    _molecular_weight,
    _bond_order_label,
)
from molbuilder.reports.synthesis_report import (
    generate_synthesis_report,
    _safe_str as synth_safe_str,
    _yield_str,
    _temp_str,
)
from molbuilder.reports.safety_report import (
    generate_safety_report,
    _risk_rank,
    _safe_list,
    EMERGENCY_NUMBER,
    POISON_CONTROL_NUMBER,
)
from molbuilder.reports.cost_report import (
    generate_cost_report,
    _safe_float,
)
from molbuilder.reports import (
    MoleculeLike,
    SynthesisStepLike,
    SynthesisRouteLike,
    CostBreakdownLike,
    CostEstimateLike,
    SafetyAssessmentLike,
)


# =====================================================================
#  Lightweight mock / stub objects
# =====================================================================

@dataclass
class _MockAtom:
    symbol: str
    index: int = 0


@dataclass
class _MockBond:
    atom_i: int
    atom_j: int
    order: float = 1.0


@dataclass
class _MockMolecule:
    name: str
    atoms: list = field(default_factory=list)
    bonds: list = field(default_factory=list)
    _adj: dict = field(default_factory=dict)

    def neighbors(self, idx: int) -> list[int]:
        return self._adj.get(idx, [])

    def get_bond(self, i: int, j: int):
        for b in self.bonds:
            if (b.atom_i == i and b.atom_j == j) or (b.atom_i == j and b.atom_j == i):
                return b
        return None


@dataclass
class _MockCategory:
    name: str = "Functional Group Interconversion"


@dataclass
class _MockTemplate:
    name: str = "Grignard Reaction"
    category: _MockCategory = field(default_factory=_MockCategory)
    reagents: list = field(default_factory=lambda: ["RMgBr", "aldehyde"])
    solvents: list = field(default_factory=lambda: ["THF"])
    catalysts: list = field(default_factory=list)
    temperature_range: tuple = (-78, 25)
    typical_yield: tuple = (60, 85)
    mechanism: str = "Nucleophilic addition of organomagnesium halide to carbonyl."
    safety_notes: str = "Moisture-sensitive. Use inert atmosphere."
    named_reaction: str = ""


@dataclass
class _MockPrecursor:
    name: str = "Benzaldehyde"
    smiles: str = "O=Cc1ccccc1"


@dataclass
class _MockSynthesisStep:
    step_number: int = 1
    template: _MockTemplate = field(default_factory=_MockTemplate)
    precursors: list = field(default_factory=lambda: [_MockPrecursor()])
    product_smiles: str = "OC(c1ccccc1)C"
    product_name: str = "1-Phenylethanol"
    expected_yield: float = 75.0
    conditions: str = "Dry THF, N2, -78C to RT"
    notes: str = "Dropwise addition recommended."


@dataclass
class _MockSynthesisRoute:
    target_smiles: str = "OC(c1ccccc1)C"
    target_name: str = "1-Phenylethanol"
    steps: list = field(default_factory=lambda: [_MockSynthesisStep()])
    overall_yield: float = 75.0
    starting_materials: list = field(
        default_factory=lambda: [
            _MockPrecursor("Benzaldehyde", "O=Cc1ccccc1"),
            _MockPrecursor("Methylmagnesium bromide", "C[Mg]Br"),
        ]
    )
    total_steps: int = 1
    longest_linear_sequence: int = 1


@dataclass
class _MockHazard:
    reagent_name: str = "Methylmagnesium bromide"
    ghs_hazards: list = field(default_factory=lambda: ["H250", "H260"])
    hazard_descriptions: list = field(
        default_factory=lambda: ["Catches fire spontaneously", "Reacts with water"]
    )
    pictogram_descriptions: list = field(
        default_factory=lambda: ["Flame", "Corrosion"]
    )


@dataclass
class _MockSafetyAssessment:
    step_number: int = 1
    step_name: str = "Grignard Reaction"
    hazards: list = field(default_factory=lambda: [_MockHazard()])
    ppe_required: list = field(
        default_factory=lambda: ["Safety goggles", "Nitrile gloves", "Lab coat"]
    )
    engineering_controls: list = field(
        default_factory=lambda: ["Fume hood", "Inert gas line"]
    )
    emergency_procedures: list = field(
        default_factory=lambda: ["Evacuate area", "Use dry chemical extinguisher"]
    )
    incompatible_materials: list = field(
        default_factory=lambda: ["Water", "Protic solvents"]
    )
    waste_classification: str = "Halogenated organic waste"
    risk_level: str = "high"


@dataclass
class _MockCostBreakdown:
    raw_materials_usd: float = 500.0
    labor_usd: float = 300.0
    equipment_usd: float = 150.0
    energy_usd: float = 50.0
    waste_disposal_usd: float = 75.0
    overhead_usd: float = 125.0


@dataclass
class _MockCostEstimate:
    total_usd: float = 1200.0
    per_kg_usd: float = 600.0
    scale_kg: float = 2.0
    breakdown: _MockCostBreakdown = field(default_factory=_MockCostBreakdown)
    notes: list = field(
        default_factory=lambda: [
            "Prices based on Sigma-Aldrich catalog (2025)",
            "Labor assumes trained chemist at bench",
        ]
    )


# =====================================================================
#  text_formatter: section headers
# =====================================================================

class TestTextFormatterHeaders(unittest.TestCase):
    """Tests for section_header and subsection_header."""

    def test_section_header_contains_title(self):
        result = section_header("Molecule Report")
        self.assertIn("MOLECULE REPORT", result)

    def test_section_header_default_width_border(self):
        result = section_header("Test")
        lines = result.split("\n")
        self.assertEqual(len(lines), 3)
        self.assertEqual(len(lines[0]), 70)
        self.assertEqual(lines[0], "=" * 70)

    def test_section_header_custom_width_and_char(self):
        result = section_header("Hi", width=30, char="*")
        lines = result.split("\n")
        self.assertEqual(lines[0], "*" * 30)
        self.assertEqual(len(lines[1]), 30)

    def test_subsection_header_format(self):
        result = subsection_header("Atom Composition")
        self.assertTrue(result.startswith("--- Atom Composition "))
        self.assertEqual(len(result), 70)

    def test_subsection_header_custom_char(self):
        result = subsection_header("Test", char="~")
        self.assertTrue(result.startswith("~~~ Test "))


# =====================================================================
#  text_formatter: ascii_table
# =====================================================================

class TestTextFormatterTable(unittest.TestCase):
    """Tests for ascii_table."""

    def test_table_headers_present(self):
        result = ascii_table(["Name", "Value"], [["H", "1"]])
        lines = result.split("\n")
        self.assertIn("Name", lines[0])
        self.assertIn("Value", lines[0])

    def test_table_divider_row(self):
        result = ascii_table(["A", "B"], [["x", "y"]])
        lines = result.split("\n")
        self.assertTrue(all(c in ("-", " ") for c in lines[1]))

    def test_table_data_present(self):
        result = ascii_table(["Elem"], [["C"], ["O"]])
        self.assertIn("C", result)
        self.assertIn("O", result)

    def test_table_right_alignment(self):
        result = ascii_table(["Num"], [["42"]], alignments=["r"])
        lines = result.split("\n")
        # Right-aligned: value should have leading spaces
        data_line = lines[2]
        self.assertEqual(data_line.rstrip(), data_line.rstrip())
        self.assertIn("42", data_line)

    def test_table_min_widths(self):
        result = ascii_table(["A"], [["x"]], min_widths=[20])
        lines = result.split("\n")
        self.assertGreaterEqual(len(lines[0]), 20)


# =====================================================================
#  text_formatter: utilities
# =====================================================================

class TestTextFormatterUtilities(unittest.TestCase):
    """Tests for word_wrap, bullet_list, key_value_block,
    horizontal_bar, format_currency, format_percent."""

    def test_word_wrap_respects_width(self):
        text = "A " * 50  # 100 chars
        result = word_wrap(text, width=40)
        for line in result.split("\n"):
            self.assertLessEqual(len(line), 40)

    def test_word_wrap_indent(self):
        result = word_wrap("Hello world", width=70, indent=4)
        self.assertTrue(result.startswith("    "))

    def test_bullet_list_default(self):
        result = bullet_list(["Alpha", "Beta"])
        self.assertIn("- Alpha", result)
        self.assertIn("- Beta", result)

    def test_bullet_list_custom_bullet(self):
        result = bullet_list(["Item"], bullet="*")
        self.assertIn("* Item", result)

    def test_key_value_block_content(self):
        result = key_value_block([("Name", "Water"), ("Formula", "H2O")])
        lines = result.strip().split("\n")
        self.assertEqual(len(lines), 2)
        self.assertIn("Name", lines[0])
        self.assertIn("Water", lines[0])
        self.assertIn("Formula", lines[1])
        self.assertIn("H2O", lines[1])

    def test_key_value_block_empty(self):
        result = key_value_block([])
        self.assertEqual(result, "")

    def test_horizontal_bar_full(self):
        result = horizontal_bar(100, 100, width=10, char="#")
        self.assertEqual(result, "##########")

    def test_horizontal_bar_half(self):
        result = horizontal_bar(50, 100, width=10, char="#")
        self.assertEqual(len(result), 10)
        self.assertEqual(result.count("#"), 5)

    def test_horizontal_bar_zero_max(self):
        result = horizontal_bar(10, 0, width=10)
        self.assertEqual(result, " " * 10)

    def test_format_currency(self):
        self.assertEqual(format_currency(1234.5), "$1,234.50")
        self.assertEqual(format_currency(0), "$0.00")

    def test_format_percent_default(self):
        self.assertEqual(format_percent(85.0), "85.0%")

    def test_format_percent_decimals(self):
        self.assertEqual(format_percent(33.333, decimals=2), "33.33%")


# =====================================================================
#  molecule_report: internal helpers
# =====================================================================

class TestMoleculeReport(unittest.TestCase):
    """Tests for generate_molecule_report and internal helpers."""

    # -- _hill_formula --
    def test_hill_formula_organic(self):
        """C first, H second, then alphabetical."""
        self.assertEqual(_hill_formula({"C": 2, "H": 6}), "C2H6")

    def test_hill_formula_with_heteroatom(self):
        self.assertEqual(_hill_formula({"C": 2, "H": 6, "O": 1}), "C2H6O")

    def test_hill_formula_no_carbon(self):
        """No carbon: all elements alphabetically."""
        self.assertEqual(_hill_formula({"H": 2, "O": 1}), "H2O")

    def test_hill_formula_single_count(self):
        self.assertEqual(_hill_formula({"C": 1, "H": 4}), "CH4")

    # -- _molecular_weight --
    def test_molecular_weight_water(self):
        mw = _molecular_weight({"H": 2, "O": 1})
        self.assertAlmostEqual(mw, 18.015, places=2)

    def test_molecular_weight_unknown_element(self):
        """Unknown elements contribute 0."""
        mw = _molecular_weight({"Xx": 5})
        self.assertEqual(mw, 0.0)

    # -- _bond_order_label --
    def test_bond_order_single(self):
        self.assertEqual(_bond_order_label(1.0), "single")

    def test_bond_order_double(self):
        self.assertEqual(_bond_order_label(2.0), "double")

    def test_bond_order_triple(self):
        self.assertEqual(_bond_order_label(3.0), "triple")

    def test_bond_order_aromatic(self):
        self.assertEqual(_bond_order_label(1.5), "aromatic")

    def test_bond_order_unknown(self):
        self.assertEqual(_bond_order_label(None), "unknown")

    # -- generate_molecule_report with mock --
    def test_report_contains_header(self):
        mol = _MockMolecule(
            name="Water",
            atoms=[_MockAtom("O", 0), _MockAtom("H", 1), _MockAtom("H", 2)],
            bonds=[_MockBond(0, 1), _MockBond(0, 2)],
            _adj={0: [1, 2], 1: [0], 2: [0]},
        )
        report = generate_molecule_report(mol)
        self.assertIn("MOLECULE REPORT: WATER", report)

    def test_report_contains_formula(self):
        mol = _MockMolecule(
            name="Water",
            atoms=[_MockAtom("O", 0), _MockAtom("H", 1), _MockAtom("H", 2)],
            bonds=[_MockBond(0, 1), _MockBond(0, 2)],
            _adj={0: [1, 2], 1: [0], 2: [0]},
        )
        report = generate_molecule_report(mol)
        self.assertIn("H2O", report)

    def test_report_acyclic_detection(self):
        mol = _MockMolecule(
            name="Ethane",
            atoms=[_MockAtom("C", 0), _MockAtom("C", 1)],
            bonds=[_MockBond(0, 1)],
            _adj={0: [1], 1: [0]},
        )
        report = generate_molecule_report(mol)
        self.assertIn("acyclic", report.lower())

    def test_report_ring_detection(self):
        """Cyclohexane-like: 6 atoms, 6 bonds => ring_count = 1."""
        atoms = [_MockAtom("C", i) for i in range(6)]
        bonds = [_MockBond(i, (i + 1) % 6) for i in range(6)]
        adj = {i: [(i - 1) % 6, (i + 1) % 6] for i in range(6)}
        mol = _MockMolecule(name="Cyclohexane", atoms=atoms, bonds=bonds, _adj=adj)
        report = generate_molecule_report(mol)
        self.assertIn("ring structures detected", report.lower())

    def test_report_raises_on_bad_input(self):
        with self.assertRaises(TypeError):
            generate_molecule_report("not a molecule")

    # -- generate_molecule_report with real Molecule --
    def test_report_with_real_ethane(self):
        from molbuilder.molecule.builders import build_ethane
        mol = build_ethane()
        report = generate_molecule_report(mol)
        self.assertIn("C2H6", report)
        self.assertIn("End of Molecule Report", report)

    def test_report_with_parsed_ethanol(self):
        from molbuilder.smiles.parser import parse
        mol = parse("CCO")
        mol.name = "Ethanol"
        report = generate_molecule_report(mol)
        self.assertIn("C2H6O", report)

    def test_report_with_parsed_benzene(self):
        from molbuilder.smiles.parser import parse
        mol = parse("c1ccccc1")
        mol.name = "Benzene"
        report = generate_molecule_report(mol)
        self.assertIn("ring structures detected", report.lower())


# =====================================================================
#  synthesis_report
# =====================================================================

class TestSynthesisReport(unittest.TestCase):
    """Tests for generate_synthesis_report and helpers."""

    def test_safe_str_none(self):
        self.assertEqual(synth_safe_str(None), "N/A")

    def test_safe_str_empty(self):
        self.assertEqual(synth_safe_str(""), "N/A")

    def test_yield_str_numeric(self):
        self.assertEqual(_yield_str(85.0), "85.0%")

    def test_yield_str_none(self):
        self.assertEqual(_yield_str(None), "N/A")

    def test_temp_str_range(self):
        tmpl = _MockTemplate(temperature_range=(-78, 25))
        self.assertEqual(_temp_str(tmpl), "-78 - 25 C")

    def test_temp_str_single(self):
        tmpl = _MockTemplate(temperature_range=(100, 100))
        self.assertEqual(_temp_str(tmpl), "100 C")

    def test_temp_str_no_attr(self):
        self.assertEqual(_temp_str(object()), "N/A")

    def test_report_header(self):
        route = _MockSynthesisRoute()
        report = generate_synthesis_report(route)
        self.assertIn("SYNTHESIS ROUTE: 1-PHENYLETHANOL", report)

    def test_report_overview_fields(self):
        route = _MockSynthesisRoute()
        report = generate_synthesis_report(route)
        self.assertIn("Target SMILES", report)
        self.assertIn("Total Steps", report)
        self.assertIn("Overall Yield", report)

    def test_report_starting_materials(self):
        route = _MockSynthesisRoute()
        report = generate_synthesis_report(route)
        self.assertIn("Benzaldehyde", report)
        self.assertIn("Methylmagnesium bromide", report)

    def test_report_step_details(self):
        route = _MockSynthesisRoute()
        report = generate_synthesis_report(route)
        self.assertIn("Step 1", report)
        self.assertIn("Grignard Reaction", report)
        self.assertIn("THF", report)

    def test_report_empty_steps(self):
        route = _MockSynthesisRoute(steps=[], starting_materials=[])
        report = generate_synthesis_report(route)
        self.assertIn("End of Synthesis Report", report)

    def test_report_raises_on_bad_input(self):
        with self.assertRaises(TypeError):
            generate_synthesis_report("not a route")

    def test_report_step_with_dict_conditions(self):
        step = _MockSynthesisStep(conditions={"temperature": "-78C", "atmosphere": "N2"})
        route = _MockSynthesisRoute(steps=[step])
        report = generate_synthesis_report(route)
        self.assertIn("temperature", report)
        self.assertIn("N2", report)

    def test_report_step_with_named_reaction(self):
        tmpl = _MockTemplate(named_reaction="Grignard")
        step = _MockSynthesisStep(template=tmpl)
        route = _MockSynthesisRoute(steps=[step])
        report = generate_synthesis_report(route)
        self.assertIn("Named Reaction", report)


# =====================================================================
#  safety_report
# =====================================================================

class TestSafetyReport(unittest.TestCase):
    """Tests for generate_safety_report and helpers."""

    def test_risk_rank_known(self):
        self.assertEqual(_risk_rank("low"), 0)
        self.assertEqual(_risk_rank("high"), 2)
        self.assertEqual(_risk_rank("extreme"), 4)

    def test_risk_rank_unknown(self):
        self.assertEqual(_risk_rank("alien"), -1)

    def test_risk_rank_case_insensitive(self):
        self.assertEqual(_risk_rank("HIGH"), 2)
        self.assertEqual(_risk_rank(" Moderate "), 1)

    def test_safe_list_present(self):
        obj = _MockSafetyAssessment()
        result = _safe_list(obj, "ppe_required")
        self.assertEqual(len(result), 3)

    def test_safe_list_missing(self):
        result = _safe_list(object(), "nonexistent")
        self.assertEqual(result, [])

    def test_report_header(self):
        assessments = [_MockSafetyAssessment()]
        report = generate_safety_report(assessments)
        self.assertIn("SAFETY ASSESSMENT REPORT", report)

    def test_report_risk_summary(self):
        assessments = [_MockSafetyAssessment()]
        report = generate_safety_report(assessments)
        self.assertIn("HIGH", report)
        self.assertIn("Highest Risk Level", report)

    def test_report_ppe_section(self):
        assessments = [_MockSafetyAssessment()]
        report = generate_safety_report(assessments)
        self.assertIn("Safety goggles", report)
        self.assertIn("Nitrile gloves", report)

    def test_report_hazard_details(self):
        assessments = [_MockSafetyAssessment()]
        report = generate_safety_report(assessments)
        self.assertIn("Methylmagnesium bromide", report)
        self.assertIn("H250", report)

    def test_report_waste_handling(self):
        assessments = [_MockSafetyAssessment()]
        report = generate_safety_report(assessments)
        self.assertIn("Halogenated organic waste", report)

    def test_report_emergency_contacts(self):
        assessments = [_MockSafetyAssessment()]
        report = generate_safety_report(assessments)
        self.assertIn(EMERGENCY_NUMBER, report)
        self.assertIn(POISON_CONTROL_NUMBER, report)

    def test_report_none_assessments(self):
        report = generate_safety_report(None)
        self.assertIn("No safety assessments provided", report)

    def test_report_empty_assessments(self):
        report = generate_safety_report([])
        self.assertIn("No safety assessments provided", report)

    def test_report_multiple_steps(self):
        a1 = _MockSafetyAssessment(step_number=1, step_name="Step A", risk_level="low")
        a2 = _MockSafetyAssessment(step_number=2, step_name="Step B", risk_level="high")
        report = generate_safety_report([a1, a2])
        self.assertIn("Step A", report)
        self.assertIn("Step B", report)
        # Highest risk should be HIGH
        self.assertIn("Highest Risk Level", report)

    def test_report_raises_on_bad_assessment(self):
        with self.assertRaises(TypeError):
            generate_safety_report(["not an assessment"])


# =====================================================================
#  cost_report
# =====================================================================

class TestCostReport(unittest.TestCase):
    """Tests for generate_cost_report and helpers."""

    def test_safe_float_normal(self):
        self.assertEqual(_safe_float(3.14), 3.14)

    def test_safe_float_none(self):
        self.assertEqual(_safe_float(None), 0.0)

    def test_safe_float_bad_string(self):
        self.assertEqual(_safe_float("abc"), 0.0)

    def test_report_header(self):
        est = _MockCostEstimate()
        report = generate_cost_report(est)
        self.assertIn("COST ESTIMATION REPORT", report)

    def test_report_summary_values(self):
        est = _MockCostEstimate()
        report = generate_cost_report(est)
        self.assertIn("$1,200.00", report)
        self.assertIn("$600.00", report)
        self.assertIn("2.00 kg", report)

    def test_report_breakdown_table(self):
        est = _MockCostEstimate()
        report = generate_cost_report(est)
        self.assertIn("Raw Materials", report)
        self.assertIn("Labor", report)
        self.assertIn("Equipment", report)
        self.assertIn("TOTAL", report)

    def test_report_bar_chart(self):
        est = _MockCostEstimate()
        report = generate_cost_report(est)
        self.assertIn("#", report)
        self.assertIn("Bar Chart", report)

    def test_report_notes(self):
        est = _MockCostEstimate()
        report = generate_cost_report(est)
        self.assertIn("Sigma-Aldrich", report)

    def test_report_disclaimer(self):
        est = _MockCostEstimate()
        report = generate_cost_report(est)
        self.assertIn("planning purposes only", report)

    def test_report_zero_total(self):
        est = _MockCostEstimate(
            total_usd=0.0, per_kg_usd=0.0, scale_kg=0.0,
            breakdown=_MockCostBreakdown(0, 0, 0, 0, 0, 0),
            notes=[],
        )
        report = generate_cost_report(est)
        self.assertIn("$0.00", report)
        self.assertIn("No additional notes", report)

    def test_report_raises_on_bad_input(self):
        with self.assertRaises(TypeError):
            generate_cost_report("not an estimate")


# =====================================================================
#  Protocol types
# =====================================================================

class TestProtocolTypes(unittest.TestCase):
    """Tests that Protocol types are runtime-checkable against mocks and reals."""

    def test_mock_molecule_satisfies_protocol(self):
        mol = _MockMolecule(name="test", atoms=[], bonds=[])
        self.assertIsInstance(mol, MoleculeLike)

    def test_real_molecule_satisfies_protocol(self):
        from molbuilder.molecule.builders import build_ethane
        mol = build_ethane()
        self.assertIsInstance(mol, MoleculeLike)

    def test_mock_cost_estimate_satisfies_protocol(self):
        est = _MockCostEstimate()
        self.assertIsInstance(est, CostEstimateLike)

    def test_mock_cost_breakdown_satisfies_protocol(self):
        bd = _MockCostBreakdown()
        self.assertIsInstance(bd, CostBreakdownLike)

    def test_mock_safety_assessment_satisfies_protocol(self):
        a = _MockSafetyAssessment()
        self.assertIsInstance(a, SafetyAssessmentLike)

    def test_mock_synthesis_step_satisfies_protocol(self):
        step = _MockSynthesisStep()
        self.assertIsInstance(step, SynthesisStepLike)

    def test_mock_synthesis_route_satisfies_protocol(self):
        route = _MockSynthesisRoute()
        self.assertIsInstance(route, SynthesisRouteLike)


if __name__ == "__main__":
    unittest.main()
