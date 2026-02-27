"""PDF report generation using ReportLab.

Converts existing ASCII text reports into professional PDF documents.
Requires the optional ``reportlab`` dependency::

    pip install molbuilder[pdf]
    # or: pip install reportlab>=4.0

The public API mirrors the text report generators but returns raw PDF bytes.
"""

from __future__ import annotations

import io
from typing import Any


def _check_reportlab():
    """Raise a helpful error if reportlab is not installed."""
    try:
        import reportlab  # noqa: F401
    except ImportError:
        raise ImportError(
            "PDF report generation requires reportlab. "
            "Install it with: pip install molbuilder[pdf]  "
            "or: pip install reportlab>=4.0"
        )


def _ascii_to_pdf_elements(ascii_text: str) -> list:
    """Parse ASCII report text into ReportLab Platypus elements.

    Recognises patterns from molbuilder.reports.text_formatter:
    - Section headers: lines of === with centred title
    - Subsection headers: lines starting with --- <title>
    - Tables: header + divider (---) + data rows
    - Bullet lists: lines starting with '  - '
    - Key-value blocks: lines with ':' separator
    - Plain text paragraphs
    """
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import inch
    from reportlab.lib.enums import TA_CENTER, TA_LEFT
    from reportlab.lib.colors import HexColor
    from reportlab.platypus import Paragraph, Spacer, Table, TableStyle

    styles = getSampleStyleSheet()

    # Custom styles
    title_style = ParagraphStyle(
        "ReportTitle",
        parent=styles["Heading1"],
        alignment=TA_CENTER,
        fontSize=16,
        spaceAfter=12,
        textColor=HexColor("#1a1a2e"),
    )
    section_style = ParagraphStyle(
        "SectionHead",
        parent=styles["Heading2"],
        fontSize=13,
        spaceBefore=16,
        spaceAfter=8,
        textColor=HexColor("#16213e"),
    )
    subsection_style = ParagraphStyle(
        "SubsectionHead",
        parent=styles["Heading3"],
        fontSize=11,
        spaceBefore=10,
        spaceAfter=6,
        textColor=HexColor("#0f3460"),
    )
    body_style = ParagraphStyle(
        "BodyMono",
        parent=styles["Normal"],
        fontName="Courier",
        fontSize=8,
        leading=10,
        spaceAfter=4,
    )
    bullet_style = ParagraphStyle(
        "BulletItem",
        parent=styles["Normal"],
        fontName="Courier",
        fontSize=8,
        leading=10,
        leftIndent=20,
        spaceAfter=2,
    )

    elements: list = []
    lines = ascii_text.splitlines()
    i = 0

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # Skip empty lines
        if not stripped:
            i += 1
            continue

        # Section header: full line of === (at least 10 chars)
        if stripped == "=" * len(stripped) and len(stripped) >= 10:
            # Next line is the title, then another === line
            if i + 2 < len(lines):
                title = lines[i + 1].strip()
                elements.append(Paragraph(_escape(title), title_style))
                i += 3  # skip both border lines + title
                continue
            i += 1
            continue

        # Subsection header: starts with '--- '
        if stripped.startswith("--- "):
            title = stripped.lstrip("- ").rstrip("- ").strip()
            elements.append(Paragraph(_escape(title), subsection_style))
            i += 1
            continue

        # Table detection: check if next line is a divider (all dashes/spaces)
        if i + 1 < len(lines):
            next_stripped = lines[i + 1].strip()
            if (next_stripped and
                    all(c in "-  " for c in next_stripped) and
                    len(next_stripped) >= 5 and
                    "  " in stripped):
                # This is a table header
                table_data, rows_consumed = _parse_ascii_table(lines, i)
                if table_data:
                    t = Table(table_data)
                    t.setStyle(TableStyle([
                        ("FONTNAME", (0, 0), (-1, 0), "Courier-Bold"),
                        ("FONTSIZE", (0, 0), (-1, -1), 7),
                        ("BOTTOMPADDING", (0, 0), (-1, 0), 4),
                        ("TOPPADDING", (0, 0), (-1, -1), 2),
                        ("BOTTOMPADDING", (0, 1), (-1, -1), 2),
                        ("LINEBELOW", (0, 0), (-1, 0), 0.5, HexColor("#333333")),
                        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
                         [HexColor("#f8f8f8"), HexColor("#ffffff")]),
                    ]))
                    elements.append(Spacer(1, 6))
                    elements.append(t)
                    elements.append(Spacer(1, 6))
                    i += rows_consumed
                    continue

        # Bullet list item
        if stripped.startswith("- ") and line.startswith("  "):
            elements.append(
                Paragraph("&#8226; " + _escape(stripped[2:]), bullet_style)
            )
            i += 1
            continue

        # Footer line (all = or - at end of report)
        if (stripped == "=" * len(stripped) or stripped == "-" * len(stripped)) and len(stripped) >= 10:
            i += 1
            continue

        # Plain text line
        elements.append(Paragraph(_escape(stripped), body_style))
        i += 1

    return elements


def _escape(text: str) -> str:
    """Escape special XML characters for ReportLab Paragraphs."""
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


def _parse_ascii_table(lines: list[str], start: int) -> tuple[list[list[str]], int]:
    """Parse an ASCII table starting at header line.

    Returns (table_data, lines_consumed).
    table_data is a list of rows (each row a list of cell strings).
    """
    header_line = lines[start]
    divider_line = lines[start + 1] if start + 1 < len(lines) else ""

    # Find column boundaries from the divider
    # Split by two or more spaces
    headers = [h.strip() for h in header_line.split("  ") if h.strip()]
    if not headers:
        return [], 1

    table_data = [headers]
    i = start + 2  # skip header + divider

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        # Stop at empty line or another section marker
        if not stripped or stripped.startswith("---") or stripped == "=" * len(stripped):
            break
        cells = [c.strip() for c in line.split("  ") if c.strip()]
        if cells:
            table_data.append(cells)
        i += 1

    return table_data, i - start


def _build_pdf(elements: list) -> bytes:
    """Build a PDF from ReportLab elements, returns raw bytes."""
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.units import inch
    from reportlab.platypus import SimpleDocTemplate, Spacer

    buf = io.BytesIO()
    doc = SimpleDocTemplate(
        buf,
        pagesize=letter,
        leftMargin=0.75 * inch,
        rightMargin=0.75 * inch,
        topMargin=0.75 * inch,
        bottomMargin=0.75 * inch,
    )
    if not elements:
        from reportlab.platypus import Paragraph
        from reportlab.lib.styles import getSampleStyleSheet
        elements = [Paragraph("No data available.", getSampleStyleSheet()["Normal"])]
    doc.build(elements)
    return buf.getvalue()


def generate_molecule_pdf(mol: Any) -> bytes:
    """Generate a PDF molecule report.

    Parameters
    ----------
    mol : Molecule
        A Molecule object.

    Returns
    -------
    bytes
        Raw PDF document bytes.
    """
    _check_reportlab()
    from molbuilder.reports.molecule_report import generate_molecule_report
    text = generate_molecule_report(mol)
    elements = _ascii_to_pdf_elements(text)
    return _build_pdf(elements)


def generate_process_pdf(
    mol_name: str,
    route: Any,
    cost: Any | None = None,
    safety: Any | None = None,
) -> bytes:
    """Generate a PDF process engineering report.

    Combines synthesis route, cost breakdown, and safety assessment
    into a single professional PDF document.

    Parameters
    ----------
    mol_name : str
        Display name for the target molecule.
    route : SynthesisRouteLike
        Synthesis route object.
    cost : CostEstimateLike | None
        Cost estimate, if available.
    safety : iterable of SafetyAssessmentLike | None
        Safety assessments, if available.

    Returns
    -------
    bytes
        Raw PDF document bytes.
    """
    _check_reportlab()
    from reportlab.platypus import Spacer

    from molbuilder.reports.synthesis_report import generate_synthesis_report
    from molbuilder.reports.cost_report import generate_cost_report
    from molbuilder.reports.safety_report import generate_safety_report

    elements: list = []

    # Synthesis route section
    synth_text = generate_synthesis_report(route)
    elements.extend(_ascii_to_pdf_elements(synth_text))
    elements.append(Spacer(1, 12))

    # Cost section
    if cost is not None:
        cost_text = generate_cost_report(cost)
        elements.extend(_ascii_to_pdf_elements(cost_text))
        elements.append(Spacer(1, 12))

    # Safety section
    if safety is not None:
        safety_text = generate_safety_report(safety)
        elements.extend(_ascii_to_pdf_elements(safety_text))

    return _build_pdf(elements)
