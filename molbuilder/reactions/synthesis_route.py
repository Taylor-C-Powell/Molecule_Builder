"""Forward synthesis route planning from retrosynthesis results.

This module traverses a ``RetrosynthesisTree`` (produced by the
retrosynthetic analysis engine) along its best disconnections and
reverses the order to produce a step-by-step forward synthesis route.

Key public functions
--------------------
extract_best_route(tree) -> SynthesisRoute
format_route(route) -> str
"""

from __future__ import annotations

from dataclasses import dataclass, field

from molbuilder.reactions.retrosynthesis import (
    RetrosynthesisTree,
    RetroNode,
    Precursor,
    PURCHASABLE_MATERIALS,
    is_purchasable,
)
from molbuilder.reactions.reaction_types import ReactionTemplate


# =====================================================================
#  Data structures
# =====================================================================

@dataclass
class SynthesisStep:
    """A single step in a forward synthesis route.

    Attributes
    ----------
    step_number : int
        Sequential step number (1-based).
    template : ReactionTemplate
        The reaction template used in this step.
    precursors : list[Precursor]
        Starting materials / intermediates consumed.
    product_smiles : str
        SMILES of the product formed in this step.
    product_name : str
        Human-readable name for the product.
    conditions : str
        Summary of temperature, solvent, and catalyst.
    expected_yield : float
        Estimated isolated yield (percent), midpoint of template range.
    notes : str
        Safety or practical notes.
    """
    step_number: int
    template: ReactionTemplate
    precursors: list[Precursor]
    product_smiles: str
    product_name: str
    conditions: str
    expected_yield: float
    notes: str


@dataclass
class SynthesisRoute:
    """A complete forward synthesis route.

    Attributes
    ----------
    target_smiles : str
        SMILES of the final target molecule.
    target_name : str
        Human-readable name (from molecule or SMILES).
    steps : list[SynthesisStep]
        Ordered list of synthesis steps (first step uses only
        purchasable materials).
    overall_yield : float
        Product of individual step yields (percent).
    starting_materials : list[Precursor]
        All purchasable starting materials needed.
    total_steps : int
        Total number of steps in the route.
    longest_linear_sequence : int
        Length of the longest linear chain of dependent steps.
    """
    target_smiles: str
    target_name: str
    steps: list[SynthesisStep] = field(default_factory=list)
    overall_yield: float = 0.0
    starting_materials: list[Precursor] = field(default_factory=list)
    total_steps: int = 0
    longest_linear_sequence: int = 0


# =====================================================================
#  Conditions summary builder
# =====================================================================

def _build_conditions(template: ReactionTemplate) -> str:
    """Build a one-line conditions summary from a ReactionTemplate.

    Includes temperature range, preferred solvent, and catalyst if any.
    """
    parts: list[str] = []

    # Temperature
    lo, hi = template.temperature_range
    if lo == hi:
        parts.append(f"{lo:.0f} C")
    else:
        parts.append(f"{lo:.0f} to {hi:.0f} C")

    # Solvent
    if template.solvents:
        parts.append(template.solvents[0])

    # Catalyst
    if template.catalysts:
        parts.append(f"cat. {template.catalysts[0]}")

    # Reagents (abbreviated)
    if template.reagents:
        abbreviated = template.reagents[:2]
        parts.append(" + ".join(abbreviated))

    return "; ".join(parts)


def _expected_yield(template: ReactionTemplate) -> float:
    """Midpoint of the template's typical yield range."""
    lo, hi = template.typical_yield
    return (lo + hi) / 2.0


def _product_name(smiles: str) -> str:
    """Return a human-readable name for a SMILES string if known."""
    entry = PURCHASABLE_MATERIALS.get(smiles)
    if entry is not None:
        return entry[0]
    return smiles


# =====================================================================
#  Tree traversal: collect steps in reverse (retro -> forward)
# =====================================================================

def _collect_retro_steps(
    node: RetroNode,
    steps_accumulator: list[tuple[RetroNode, ReactionTemplate, list[Precursor]]],
    visited: set[str],
) -> None:
    """Walk the retrosynthesis tree depth-first, collecting steps.

    Each non-leaf, non-purchasable node with a best_disconnection
    contributes one step.  Children are visited first (depth-first)
    so that when the list is later reversed, leaf-level reactions
    come first in the forward direction.
    """
    if node.smiles in visited:
        return
    visited.add(node.smiles)

    if node.is_purchasable:
        return

    # First recurse into children (precursors)
    for child in node.children:
        _collect_retro_steps(child, steps_accumulator, visited)

    # Then record this node's disconnection
    if node.best_disconnection is not None:
        steps_accumulator.append((
            node,
            node.best_disconnection.template,
            node.best_disconnection.precursors,
        ))


def _gather_purchasable_leaves(node: RetroNode, result: list[Precursor],
                               seen: set[str]) -> None:
    """Collect all purchasable leaf nodes as Precursor objects."""
    if node.is_purchasable:
        if node.smiles not in seen:
            seen.add(node.smiles)
            entry = PURCHASABLE_MATERIALS.get(node.smiles)
            name = entry[0] if entry else node.smiles
            cost = entry[1] if entry else 50.0
            result.append(Precursor(
                smiles=node.smiles, molecule=None,
                name=name, cost_per_kg=cost,
            ))
        return
    for child in node.children:
        _gather_purchasable_leaves(child, result, seen)


def _compute_longest_linear(node: RetroNode) -> int:
    """Compute the longest linear sequence of steps from the root.

    The longest linear sequence is the depth of the deepest non-
    purchasable node.
    """
    if node.is_purchasable or not node.children:
        return 0
    return 1 + max(_compute_longest_linear(c) for c in node.children)


# =====================================================================
#  Public API
# =====================================================================

def extract_best_route(tree: RetrosynthesisTree) -> SynthesisRoute:
    """Extract the best forward synthesis route from a retrosynthesis tree.

    Traverses the tree along ``best_disconnection`` links, reverses the
    order to produce a forward synthesis plan, and computes the overall
    yield as the product of individual step yields.

    Parameters
    ----------
    tree : RetrosynthesisTree
        The retrosynthesis result from ``retrosynthesis()``.

    Returns
    -------
    SynthesisRoute
        A forward synthesis route with ordered steps, starting from
        purchasable materials and ending at the target.
    """
    root = tree.target

    # Collect retro steps (deepest first)
    retro_steps: list[tuple[RetroNode, ReactionTemplate, list[Precursor]]] = []
    visited: set[str] = set()
    _collect_retro_steps(root, retro_steps, visited)

    # retro_steps are already in forward order because we recurse into
    # children before appending the parent.  If we had reversed, the
    # deepest transformations would come first (correct for forward
    # synthesis).  Since _collect_retro_steps already visits children
    # first, the list is naturally in forward order.

    # Build SynthesisStep objects
    steps: list[SynthesisStep] = []
    overall_yield = 100.0

    for i, (node, template, precursors) in enumerate(retro_steps):
        step_yield = _expected_yield(template)
        overall_yield *= (step_yield / 100.0)

        notes_parts: list[str] = []
        if template.safety_notes:
            notes_parts.append(template.safety_notes)
        if template.scale_notes:
            notes_parts.append(template.scale_notes)
        notes = " ".join(notes_parts) if notes_parts else ""

        step = SynthesisStep(
            step_number=i + 1,
            template=template,
            precursors=precursors,
            product_smiles=node.smiles,
            product_name=_product_name(node.smiles),
            conditions=_build_conditions(template),
            expected_yield=step_yield,
            notes=notes,
        )
        steps.append(step)

    # Gather all purchasable starting materials
    starting_materials: list[Precursor] = []
    seen_sm: set[str] = set()
    _gather_purchasable_leaves(root, starting_materials, seen_sm)

    # Longest linear sequence
    lls = _compute_longest_linear(root)

    target_name = _product_name(root.smiles)

    return SynthesisRoute(
        target_smiles=root.smiles,
        target_name=target_name,
        steps=steps,
        overall_yield=overall_yield,
        starting_materials=starting_materials,
        total_steps=len(steps),
        longest_linear_sequence=lls,
    )


def format_route(route: SynthesisRoute) -> str:
    """Format a SynthesisRoute as a readable ASCII text summary.

    Parameters
    ----------
    route : SynthesisRoute
        The synthesis route to format.

    Returns
    -------
    str
        Multi-line text summary of the route.

    Example output::

        ============================================================
        Forward Synthesis Route
        ============================================================
        Target  : CC(=O)OCC (ethyl acetate)
        Steps   : 1
        Overall yield : 67.5%
        Longest linear sequence : 1
        ============================================================

        Starting Materials:
          - CC(O)=O  (acetic acid, $1.50/kg)
          - CCO  (ethanol, $2.00/kg)

        ------------------------------------------------------------
        Step 1: Fischer esterification
        ------------------------------------------------------------
          Precursors : CC(O)=O + CCO
          Product    : CC(=O)OCC
          Conditions : 60 to 120 C; toluene (Dean-Stark); ...
          Yield      : 67.5%
          Notes      : ...
        ------------------------------------------------------------

        Overall yield: 67.5%
        ============================================================
    """
    sep = "=" * 60
    thin_sep = "-" * 60
    lines: list[str] = []

    lines.append(sep)
    lines.append("Forward Synthesis Route")
    lines.append(sep)
    lines.append(f"Target  : {route.target_smiles} ({route.target_name})")
    lines.append(f"Steps   : {route.total_steps}")
    lines.append(f"Overall yield : {route.overall_yield:.1f}%")
    lines.append(f"Longest linear sequence : {route.longest_linear_sequence}")
    lines.append(sep)

    # Starting materials
    lines.append("")
    lines.append("Starting Materials:")
    if route.starting_materials:
        for sm in route.starting_materials:
            lines.append(f"  - {sm.smiles}  ({sm.name}, ${sm.cost_per_kg:.2f}/kg)")
    else:
        lines.append("  (none identified)")

    # Steps
    for step in route.steps:
        lines.append("")
        lines.append(thin_sep)
        lines.append(f"Step {step.step_number}: {step.template.name}")
        if step.template.named_reaction:
            lines.append(f"  Named reaction : {step.template.named_reaction}")
        lines.append(thin_sep)

        precursor_str = " + ".join(p.smiles for p in step.precursors)
        lines.append(f"  Precursors : {precursor_str}")
        lines.append(f"  Product    : {step.product_smiles}"
                     f" ({step.product_name})")
        lines.append(f"  Conditions : {step.conditions}")
        lines.append(f"  Yield      : {step.expected_yield:.1f}%")
        if step.notes:
            lines.append(f"  Notes      : {step.notes}")

    lines.append("")
    lines.append(thin_sep)
    lines.append(f"Overall yield: {route.overall_yield:.1f}%")
    lines.append(sep)

    return "\n".join(lines)
