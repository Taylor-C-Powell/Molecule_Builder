"""Synthesis route report generator.

Produces a detailed ASCII report of a multi-step synthesis route
including reagents, conditions, yields, and starting materials.
All output is cp1252-safe.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from molbuilder.reports.text_formatter import (
    section_header,
    subsection_header,
    ascii_table,
    bullet_list,
    key_value_block,
    word_wrap,
    format_percent,
    format_currency,
)

if TYPE_CHECKING:
    from molbuilder.reports import SynthesisRouteLike


# =====================================================================
#  Helpers
# =====================================================================

def _safe_str(value, default: str = "N/A") -> str:
    """Convert *value* to a string, returning *default* if None or empty."""
    if value is None:
        return default
    s = str(value)
    return s if s else default


def _yield_str(yld) -> str:
    """Format a yield value (float 0-100 or None) as a percentage."""
    if yld is None:
        return "N/A"
    try:
        return format_percent(float(yld))
    except (TypeError, ValueError):
        return _safe_str(yld)


def _temp_str(template) -> str:
    """Extract a human-readable temperature string from a template."""
    try:
        lo, hi = template.temperature_range
        if lo == hi:
            return f"{lo:.0f} C"
        return f"{lo:.0f} - {hi:.0f} C"
    except Exception:
        return "N/A"


# =====================================================================
#  Public API
# =====================================================================

def generate_synthesis_report(route: SynthesisRouteLike) -> str:
    """Generate an ASCII synthesis route report.

    Uses duck typing -- *route* should expose:

    * ``.target_smiles`` -- SMILES of the target molecule
    * ``.target_name`` -- human-readable target name
    * ``.steps`` -- list of step objects (see below)
    * ``.overall_yield`` -- float (percent)
    * ``.starting_materials`` -- list of material objects with ``.name``, ``.smiles``
    * ``.total_steps`` -- int
    * ``.longest_linear_sequence`` -- int

    Each **step** should expose:

    * ``.step_number`` -- int
    * ``.template`` -- ReactionTemplate (has ``.name``, ``.category``,
      ``.reagents``, ``.solvents``, ``.catalysts``, ``.temperature_range``,
      ``.typical_yield``, ``.mechanism``, ``.safety_notes``)
    * ``.precursors`` -- list with ``.name``, ``.smiles``
    * ``.product_smiles`` -- str
    * ``.product_name`` -- str
    * ``.conditions`` -- str or dict
    * ``.expected_yield`` -- float (percent)
    * ``.notes`` -- str
    """
    if not hasattr(route, 'steps'):
        raise TypeError(
            f"route must have a 'steps' attribute, "
            f"got {type(route).__name__}"
        )

    lines: list[str] = []

    # ------------------------------------------------------------------
    # 1. Header
    # ------------------------------------------------------------------
    target_name = _safe_str(getattr(route, "target_name", None), "Target")
    lines.append(section_header(f"Synthesis Route: {target_name}"))
    lines.append("")

    # ------------------------------------------------------------------
    # 2. Route Overview
    # ------------------------------------------------------------------
    lines.append(subsection_header("Route Overview"))

    target_smiles = _safe_str(getattr(route, "target_smiles", None))
    total_steps = getattr(route, "total_steps", None)
    lls = getattr(route, "longest_linear_sequence", None)
    overall_yield = getattr(route, "overall_yield", None)

    overview = [
        ("Target", target_name),
        ("Target SMILES", target_smiles),
        ("Total Steps", _safe_str(total_steps)),
        ("Longest Linear Sequence", _safe_str(lls)),
        ("Overall Yield", _yield_str(overall_yield)),
    ]
    lines.append(key_value_block(overview))
    lines.append("")

    # Starting materials list
    starting_materials = getattr(route, "starting_materials", []) or []
    if starting_materials:
        lines.append(subsection_header("Starting Materials"))
        sm_headers = ["#", "Name", "SMILES"]
        sm_rows: list[list[str]] = []
        for i, sm in enumerate(starting_materials, 1):
            sm_name = _safe_str(getattr(sm, "name", None), f"Material {i}")
            sm_smiles = _safe_str(getattr(sm, "smiles", None))
            sm_rows.append([str(i), sm_name, sm_smiles])
        lines.append(ascii_table(
            sm_headers, sm_rows,
            alignments=["r", "l", "l"],
            min_widths=[3, 25, 25],
        ))
        lines.append("")

    # ------------------------------------------------------------------
    # 3. Step-by-step Details
    # ------------------------------------------------------------------
    steps = getattr(route, "steps", []) or []
    if steps:
        lines.append(section_header("Step-by-Step Details"))
        lines.append("")

    for step in steps:
        step_num = getattr(step, "step_number", "?")
        product_name = _safe_str(getattr(step, "product_name", None))
        product_smiles = _safe_str(getattr(step, "product_smiles", None))

        lines.append(subsection_header(f"Step {step_num}: {product_name}"))

        # Template info
        template = getattr(step, "template", None)
        if template is not None:
            rxn_name = _safe_str(getattr(template, "name", None), "Unknown reaction")
            named_rxn = getattr(template, "named_reaction", None)
            category = getattr(template, "category", None)
            cat_str = category.name if category is not None else "N/A"

            step_info = [
                ("Reaction", rxn_name),
            ]
            if named_rxn:
                step_info.append(("Named Reaction", str(named_rxn)))
            step_info.append(("Category", cat_str))
            step_info.append(("Product SMILES", product_smiles))
            lines.append(key_value_block(step_info))
            lines.append("")

            # Reagents
            reagents = getattr(template, "reagents", []) or []
            if reagents:
                lines.append("  Reagents:")
                lines.append(bullet_list(reagents, indent=4))
                lines.append("")

            # Solvents
            solvents = getattr(template, "solvents", []) or []
            if solvents:
                lines.append("  Solvents:")
                lines.append(bullet_list(solvents, indent=4))
                lines.append("")

            # Catalysts
            catalysts = getattr(template, "catalysts", []) or []
            if catalysts:
                lines.append("  Catalysts:")
                lines.append(bullet_list(catalysts, indent=4))
                lines.append("")

            # Temperature
            lines.append(f"  Temperature: {_temp_str(template)}")

            # Typical yield range from template
            try:
                lo, hi = template.typical_yield
                lines.append(f"  Typical Yield Range: {lo:.0f}% - {hi:.0f}%")
            except Exception:
                pass
        else:
            lines.append(key_value_block([
                ("Product SMILES", product_smiles),
            ]))
            lines.append("")

        # Expected yield for this specific step
        expected_yield = getattr(step, "expected_yield", None)
        lines.append(f"  Expected Yield: {_yield_str(expected_yield)}")

        # Precursors
        precursors = getattr(step, "precursors", []) or []
        if precursors:
            lines.append("")
            lines.append("  Precursors:")
            for p in precursors:
                p_name = _safe_str(getattr(p, "name", None), "Unknown")
                p_smiles = _safe_str(getattr(p, "smiles", None))
                lines.append(f"    - {p_name}  ({p_smiles})")

        # Conditions (may be a string or dict)
        conditions = getattr(step, "conditions", None)
        if conditions:
            lines.append("")
            if isinstance(conditions, dict):
                cond_pairs = [(str(k), str(v)) for k, v in conditions.items()]
                lines.append("  Conditions:")
                lines.append(key_value_block(cond_pairs, indent=4))
            else:
                lines.append(word_wrap(f"  Conditions: {conditions}", indent=2))

        # Mechanism description from template
        if template is not None:
            mechanism = getattr(template, "mechanism", "")
            if mechanism:
                lines.append("")
                lines.append("  Mechanism:")
                lines.append(word_wrap(mechanism, indent=4))

            safety = getattr(template, "safety_notes", "")
            if safety:
                lines.append("")
                lines.append("  Safety Notes:")
                lines.append(word_wrap(safety, indent=4))

        # Step notes
        notes = getattr(step, "notes", "")
        if notes:
            lines.append("")
            lines.append("  Notes:")
            lines.append(word_wrap(str(notes), indent=4))

        lines.append("")

    # ------------------------------------------------------------------
    # 4. Starting Materials Summary Table
    # ------------------------------------------------------------------
    if starting_materials:
        lines.append(section_header("Starting Materials Summary"))
        lines.append("")
        sum_headers = ["Material", "SMILES"]
        sum_rows: list[list[str]] = []
        for sm in starting_materials:
            sm_name = _safe_str(getattr(sm, "name", None), "Unknown")
            sm_smiles = _safe_str(getattr(sm, "smiles", None))
            sum_rows.append([sm_name, sm_smiles])
        lines.append(ascii_table(
            sum_headers, sum_rows,
            alignments=["l", "l"],
            min_widths=[30, 30],
        ))
        lines.append("")

    # Footer
    lines.append("=" * 70)
    lines.append("  End of Synthesis Report")
    lines.append("=" * 70)

    return "\n".join(lines)
