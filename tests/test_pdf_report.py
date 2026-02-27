"""Tests for PDF report generation."""

import unittest

reportlab = __import__("pytest").importorskip("reportlab")

from molbuilder.smiles.parser import parse
from molbuilder.reports.pdf_report import (
    generate_molecule_pdf,
    generate_process_pdf,
    _ascii_to_pdf_elements,
)


class TestPdfReport(unittest.TestCase):
    """PDF report generation tests (requires reportlab)."""

    def test_molecule_pdf_bytes(self):
        """generate_molecule_pdf returns valid PDF bytes."""
        mol = parse("CCO")
        mol.name = "ethanol"
        pdf = generate_molecule_pdf(mol)
        self.assertIsInstance(pdf, bytes)
        self.assertTrue(pdf.startswith(b"%PDF-"))

    def test_molecule_pdf_nonempty(self):
        """PDF for benzene should have reasonable size."""
        mol = parse("c1ccccc1")
        mol.name = "benzene"
        pdf = generate_molecule_pdf(mol)
        self.assertGreater(len(pdf), 500)

    def test_ascii_to_elements_section_header(self):
        """Section headers (=== lines) are parsed into elements."""
        text = (
            "=" * 70 + "\n"
            + "SYNTHESIS REPORT".center(70) + "\n"
            + "=" * 70
        )
        elements = _ascii_to_pdf_elements(text)
        self.assertGreater(len(elements), 0)

    def test_ascii_to_elements_subsection(self):
        """Subsection headers (--- Title ---) are parsed."""
        text = "--- Step 1: Grignard Reaction " + "-" * 40
        elements = _ascii_to_pdf_elements(text)
        self.assertGreater(len(elements), 0)

    def test_ascii_to_elements_table(self):
        """Tables (header + divider + rows) are parsed."""
        text = (
            "Name          Value\n"
            "----          -----\n"
            "Temperature   100 C\n"
            "Pressure      1 atm\n"
        )
        elements = _ascii_to_pdf_elements(text)
        self.assertGreater(len(elements), 0)

    def test_ascii_to_elements_bullets(self):
        """Bullet lists are parsed."""
        text = "  - First item\n  - Second item\n  - Third item\n"
        elements = _ascii_to_pdf_elements(text)
        self.assertEqual(len(elements), 3)

    def test_process_pdf_with_route(self):
        """generate_process_pdf with a synthesis route produces valid PDF."""
        from molbuilder.reactions.retrosynthesis import retrosynthesis
        from molbuilder.reactions.synthesis_route import extract_best_route

        mol = parse("CCO")
        mol.name = "ethanol"
        tree = retrosynthesis(mol)
        route = extract_best_route(tree)
        if route is None:
            self.skipTest("No route found for ethanol")

        pdf = generate_process_pdf(
            mol_name="ethanol",
            route=route,
        )
        self.assertIsInstance(pdf, bytes)
        self.assertTrue(pdf.startswith(b"%PDF-"))


if __name__ == "__main__":
    unittest.main()
