"""Cost estimation report generator.

Produces an ASCII cost breakdown report with summary, category
breakdown table, bar chart, percentage distribution, and notes.
All output is cp1252-safe.
"""

from __future__ import annotations

from molbuilder.reports.text_formatter import (
    section_header,
    subsection_header,
    ascii_table,
    bullet_list,
    key_value_block,
    horizontal_bar,
    format_currency,
    format_percent,
    word_wrap,
)


# =====================================================================
#  Helpers
# =====================================================================

_CATEGORY_LABELS: list[tuple[str, str]] = [
    ("raw_materials_usd",   "Raw Materials"),
    ("labor_usd",           "Labor"),
    ("equipment_usd",       "Equipment"),
    ("energy_usd",          "Energy"),
    ("waste_disposal_usd",  "Waste Disposal"),
    ("overhead_usd",        "Overhead"),
]


def _safe_float(value, default: float = 0.0) -> float:
    """Coerce *value* to float, returning *default* on failure."""
    if value is None:
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


# =====================================================================
#  Public API
# =====================================================================

def generate_cost_report(estimate) -> str:
    """Generate an ASCII cost breakdown report.

    Uses duck typing -- *estimate* should expose:

    * ``.total_usd`` -- float, total estimated cost
    * ``.per_kg_usd`` -- float, cost per kilogram of product
    * ``.scale_kg`` -- float, production scale in kg
    * ``.breakdown`` -- object with attributes:
        - ``.raw_materials_usd``
        - ``.labor_usd``
        - ``.equipment_usd``
        - ``.energy_usd``
        - ``.waste_disposal_usd``
        - ``.overhead_usd``
    * ``.notes`` -- list[str], assumptions and caveats
    """
    lines: list[str] = []

    # ------------------------------------------------------------------
    # 1. Header
    # ------------------------------------------------------------------
    lines.append(section_header("Cost Estimation Report"))
    lines.append("")

    # ------------------------------------------------------------------
    # 2. Summary
    # ------------------------------------------------------------------
    lines.append(subsection_header("Summary"))

    total = _safe_float(getattr(estimate, "total_usd", None))
    per_kg = _safe_float(getattr(estimate, "per_kg_usd", None))
    scale = _safe_float(getattr(estimate, "scale_kg", None))

    summary_pairs = [
        ("Production Scale",  f"{scale:.2f} kg"),
        ("Total Cost",        format_currency(total)),
        ("Cost per kg",       format_currency(per_kg)),
    ]
    lines.append(key_value_block(summary_pairs))
    lines.append("")

    # ------------------------------------------------------------------
    # 3. Cost Breakdown -- table and bar chart
    # ------------------------------------------------------------------
    lines.append(subsection_header("Cost Breakdown"))

    breakdown = getattr(estimate, "breakdown", None)
    categories: list[tuple[str, float]] = []
    for attr, label in _CATEGORY_LABELS:
        val = _safe_float(getattr(breakdown, attr, None) if breakdown else None)
        categories.append((label, val))

    # Table
    tbl_headers = ["Category", "Amount (USD)", "% of Total"]
    tbl_rows: list[list[str]] = []
    for label, amount in categories:
        pct = (amount / total * 100.0) if total > 0 else 0.0
        tbl_rows.append([
            label,
            format_currency(amount),
            format_percent(pct),
        ])
    # Total row
    tbl_rows.append([
        "TOTAL",
        format_currency(total),
        format_percent(100.0) if total > 0 else format_percent(0.0),
    ])

    lines.append(ascii_table(
        tbl_headers, tbl_rows,
        alignments=["l", "r", "r"],
        min_widths=[16, 14, 10],
    ))
    lines.append("")

    # Bar chart
    lines.append(subsection_header("Cost Distribution (Bar Chart)"))
    max_amount = max((amt for _, amt in categories), default=0.0)
    max_label_len = max((len(label) for label, _ in categories), default=0)

    for label, amount in categories:
        bar = horizontal_bar(amount, max_amount, width=35, char="#")
        amount_str = format_currency(amount)
        lines.append(
            f"  {label:<{max_label_len}}  |{bar}| {amount_str}"
        )
    lines.append("")

    # ------------------------------------------------------------------
    # 4. Percentage Distribution
    # ------------------------------------------------------------------
    lines.append(subsection_header("Percentage Distribution"))

    pct_headers = ["Category", "Percentage", "Visual"]
    pct_rows: list[list[str]] = []
    for label, amount in categories:
        pct = (amount / total * 100.0) if total > 0 else 0.0
        bar = horizontal_bar(pct, 100.0, width=20, char="=")
        pct_rows.append([label, format_percent(pct), bar])

    lines.append(ascii_table(
        pct_headers, pct_rows,
        alignments=["l", "r", "l"],
        min_widths=[16, 10, 20],
    ))
    lines.append("")

    # ------------------------------------------------------------------
    # 5. Notes and Assumptions
    # ------------------------------------------------------------------
    lines.append(subsection_header("Notes and Assumptions"))

    notes = getattr(estimate, "notes", None)
    if notes and hasattr(notes, "__iter__") and not isinstance(notes, str):
        note_list = [str(n) for n in notes if n]
    elif notes:
        note_list = [str(notes)]
    else:
        note_list = []

    if note_list:
        lines.append(bullet_list(note_list))
    else:
        lines.append("  No additional notes.")
    lines.append("")

    # Standard disclaimer
    lines.append(subsection_header("Disclaimer"))
    disclaimer = (
        "This cost estimate is for planning purposes only.  Actual costs "
        "may vary based on supplier pricing, scale-up effects, local "
        "labor rates, and regulatory requirements.  All figures are "
        "approximate and should be validated before procurement."
    )
    lines.append(word_wrap(disclaimer, width=66, indent=2))
    lines.append("")

    # Footer
    lines.append("=" * 70)
    lines.append("  End of Cost Estimation Report")
    lines.append("=" * 70)

    return "\n".join(lines)
