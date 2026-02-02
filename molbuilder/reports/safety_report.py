"""Safety assessment report generator.

Produces a detailed ASCII safety report covering hazards, PPE
requirements, engineering controls, waste handling, and emergency
procedures for every step of a synthesis route.
All output is cp1252-safe.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

from molbuilder.reports.text_formatter import (
    section_header,
    subsection_header,
    ascii_table,
    bullet_list,
    key_value_block,
    word_wrap,
)

if TYPE_CHECKING:
    from molbuilder.reports import SafetyAssessmentLike


# =====================================================================
#  Helpers
# =====================================================================

# Emergency contact numbers -- override these for non-US locales
EMERGENCY_NUMBER = "911"
POISON_CONTROL_NUMBER = "1-800-222-1222"

_RISK_ORDER = {
    "low": 0,
    "moderate": 1,
    "high": 2,
    "very high": 3,
    "extreme": 4,
}


def _risk_rank(level: str) -> int:
    """Numeric rank for a risk-level string (higher = worse)."""
    return _RISK_ORDER.get(str(level).strip().lower(), -1)


def _safe_list(obj, attr: str) -> list:
    """Return getattr(obj, attr) as a list, or [] on failure."""
    val = getattr(obj, attr, None)
    if val is None:
        return []
    if isinstance(val, (list, tuple)):
        return list(val)
    return [val]


def _safe_str(value, default: str = "N/A") -> str:
    if value is None:
        return default
    s = str(value)
    return s if s else default


# =====================================================================
#  Public API
# =====================================================================

def generate_safety_report(
    assessments: Iterable[SafetyAssessmentLike] | None,
) -> str:
    """Generate an ASCII safety report.

    Uses duck typing -- *assessments* should be a list (or iterable)
    where each item exposes:

    * ``.step_number`` -- int
    * ``.step_name`` -- str
    * ``.hazards`` -- list of hazard objects, each with:
        - ``.reagent_name`` -- str
        - ``.ghs_hazards`` -- list[str]  (e.g. ``["H225", "H301"]``)
        - ``.hazard_descriptions`` -- list[str]
        - ``.pictogram_descriptions`` -- list[str]
    * ``.ppe_required`` -- list[str]
    * ``.engineering_controls`` -- list[str]
    * ``.emergency_procedures`` -- list[str]
    * ``.incompatible_materials`` -- list[str]
    * ``.waste_classification`` -- str
    * ``.risk_level`` -- str  (e.g. ``"low"``, ``"high"``)
    """
    if assessments is not None:
        for i, a in enumerate(assessments):
            if not hasattr(a, 'step_number'):
                raise TypeError(
                    f"Assessment {i} must have a 'step_number' attribute, "
                    f"got {type(a).__name__}"
                )

    lines: list[str] = []

    # ------------------------------------------------------------------
    # 1. Header
    # ------------------------------------------------------------------
    lines.append(section_header("Safety Assessment Report"))
    lines.append("")

    assessment_list = list(assessments) if assessments else []
    if not assessment_list:
        lines.append("  No safety assessments provided.")
        lines.append("")
        lines.append("=" * 70)
        return "\n".join(lines)

    # ------------------------------------------------------------------
    # 2. Overall Risk Summary
    # ------------------------------------------------------------------
    lines.append(subsection_header("Overall Risk Summary"))

    risk_levels: list[str] = []
    for a in assessment_list:
        rl = _safe_str(getattr(a, "risk_level", None), "unknown")
        risk_levels.append(rl)

    # Determine highest risk
    highest_risk = "unknown"
    highest_rank = -1
    for rl in risk_levels:
        rank = _risk_rank(rl)
        if rank > highest_rank:
            highest_rank = rank
            highest_risk = rl

    overview = [
        ("Total Steps Assessed", str(len(assessment_list))),
        ("Highest Risk Level",   highest_risk.upper()),
    ]
    lines.append(key_value_block(overview))
    lines.append("")

    # Per-step risk overview table
    risk_headers = ["Step", "Name", "Risk Level"]
    risk_rows: list[list[str]] = []
    for a in assessment_list:
        sn = _safe_str(getattr(a, "step_number", None), "?")
        nm = _safe_str(getattr(a, "step_name", None), "Unnamed")
        rl = _safe_str(getattr(a, "risk_level", None), "unknown")
        risk_rows.append([str(sn), nm, rl.upper()])
    lines.append(ascii_table(
        risk_headers, risk_rows,
        alignments=["r", "l", "c"],
        min_widths=[5, 30, 14],
    ))
    lines.append("")

    # ------------------------------------------------------------------
    # 3. Consolidated PPE Requirements
    # ------------------------------------------------------------------
    lines.append(subsection_header("PPE Requirements (All Steps)"))
    all_ppe: set[str] = set()
    for a in assessment_list:
        for item in _safe_list(a, "ppe_required"):
            all_ppe.add(str(item))

    if all_ppe:
        lines.append(bullet_list(sorted(all_ppe)))
    else:
        lines.append("  No specific PPE requirements listed.")
    lines.append("")

    # ------------------------------------------------------------------
    # 4. Per-Step Hazard Details
    # ------------------------------------------------------------------
    lines.append(section_header("Per-Step Hazard Details"))
    lines.append("")

    for a in assessment_list:
        step_num = _safe_str(getattr(a, "step_number", None), "?")
        step_name = _safe_str(getattr(a, "step_name", None), "Unnamed")
        risk_level = _safe_str(getattr(a, "risk_level", None), "unknown")

        lines.append(subsection_header(
            f"Step {step_num}: {step_name}  [Risk: {risk_level.upper()}]"
        ))

        # Reagent hazards
        hazards = _safe_list(a, "hazards")
        if hazards:
            for h in hazards:
                reagent_name = _safe_str(getattr(h, "reagent_name", None),
                                         "Unknown reagent")
                lines.append(f"  Reagent: {reagent_name}")

                ghs = _safe_list(h, "ghs_hazards")
                if ghs:
                    lines.append(f"    GHS Codes: {', '.join(str(c) for c in ghs)}")

                descs = _safe_list(h, "hazard_descriptions")
                if descs:
                    lines.append("    Hazard Descriptions:")
                    lines.append(bullet_list(
                        [str(d) for d in descs], indent=6))

                pictos = _safe_list(h, "pictogram_descriptions")
                if pictos:
                    lines.append("    Pictograms:")
                    lines.append(bullet_list(
                        [str(p) for p in pictos], indent=6))
                lines.append("")
        else:
            lines.append("  No reagent hazards listed.")
            lines.append("")

        # Engineering controls
        controls = _safe_list(a, "engineering_controls")
        if controls:
            lines.append("  Engineering Controls:")
            lines.append(bullet_list([str(c) for c in controls], indent=4))
            lines.append("")

        # Incompatible materials
        incompat = _safe_list(a, "incompatible_materials")
        if incompat:
            lines.append("  Incompatible Materials:")
            lines.append(bullet_list([str(m) for m in incompat], indent=4))
            lines.append("")

    # ------------------------------------------------------------------
    # 5. Waste Handling
    # ------------------------------------------------------------------
    lines.append(section_header("Waste Handling"))
    lines.append("")

    waste_headers = ["Step", "Name", "Waste Classification"]
    waste_rows: list[list[str]] = []
    for a in assessment_list:
        sn = _safe_str(getattr(a, "step_number", None), "?")
        nm = _safe_str(getattr(a, "step_name", None), "Unnamed")
        wc = _safe_str(getattr(a, "waste_classification", None), "Unclassified")
        waste_rows.append([str(sn), nm, wc])

    lines.append(ascii_table(
        waste_headers, waste_rows,
        alignments=["r", "l", "l"],
        min_widths=[5, 25, 25],
    ))
    lines.append("")

    lines.append("  General Waste Disposal Guidance:")
    disposal_notes = [
        "Segregate waste by classification (halogenated, aqueous, "
        "heavy-metal, etc.)",
        "Label all waste containers with contents, hazard class, "
        "and date",
        "Never mix incompatible waste streams",
        "Contact EHS for disposal of unknown or high-hazard waste",
    ]
    lines.append(bullet_list(disposal_notes, indent=4))
    lines.append("")

    # ------------------------------------------------------------------
    # 6. Emergency Procedures
    # ------------------------------------------------------------------
    lines.append(section_header("Emergency Procedures"))
    lines.append("")

    all_emergency: list[str] = []
    for a in assessment_list:
        step_num = _safe_str(getattr(a, "step_number", None), "?")
        step_name = _safe_str(getattr(a, "step_name", None), "Unnamed")
        procs = _safe_list(a, "emergency_procedures")
        if procs:
            lines.append(subsection_header(f"Step {step_num}: {step_name}"))
            lines.append(bullet_list([str(p) for p in procs], indent=4))
            lines.append("")
            all_emergency.extend(str(p) for p in procs)

    if not all_emergency:
        lines.append("  No step-specific emergency procedures listed.")
        lines.append("")

    lines.append(subsection_header("General Emergency Contacts"))
    contacts = [
        ("Emergency Services",     EMERGENCY_NUMBER),
        ("Poison Control",         POISON_CONTROL_NUMBER),
        ("Campus/Site EHS",        "See posted numbers in laboratory"),
    ]
    lines.append(key_value_block(contacts, indent=4))
    lines.append("")

    # Footer
    lines.append("=" * 70)
    lines.append("  End of Safety Assessment Report")
    lines.append("=" * 70)

    return "\n".join(lines)
