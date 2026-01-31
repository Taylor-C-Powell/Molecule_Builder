"""Report generation: synthesis, safety, costing, molecule summaries.

Public API
----------
Text formatting utilities:
    section_header, subsection_header, ascii_table, word_wrap,
    bullet_list, key_value_block, horizontal_bar,
    format_currency, format_percent

Report generators:
    generate_molecule_report  -- comprehensive molecule analysis
    generate_synthesis_report -- multi-step synthesis route
    generate_safety_report    -- hazard and PPE assessment
    generate_cost_report      -- cost breakdown with charts
"""

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

from molbuilder.reports.molecule_report import generate_molecule_report
from molbuilder.reports.synthesis_report import generate_synthesis_report
from molbuilder.reports.safety_report import generate_safety_report
from molbuilder.reports.cost_report import generate_cost_report

__all__ = [
    # Formatters
    "section_header",
    "subsection_header",
    "ascii_table",
    "word_wrap",
    "bullet_list",
    "key_value_block",
    "horizontal_bar",
    "format_currency",
    "format_percent",
    # Report generators
    "generate_molecule_report",
    "generate_synthesis_report",
    "generate_safety_report",
    "generate_cost_report",
]
