"""Reaction categories, templates, and data structures."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum, auto


class ReactionCategory(Enum):
    SUBSTITUTION = auto()
    ELIMINATION = auto()
    ADDITION = auto()
    OXIDATION = auto()
    REDUCTION = auto()
    COUPLING = auto()
    CARBONYL = auto()
    PROTECTION = auto()
    DEPROTECTION = auto()
    REARRANGEMENT = auto()
    RADICAL = auto()
    PERICYCLIC = auto()
    POLYMERIZATION = auto()
    MISC = auto()


@dataclass
class ReactionTemplate:
    """A curated reaction template for synthesis planning.

    Each template captures the essential information needed to evaluate
    whether a reaction is applicable to a given substrate and to estimate
    its likely outcome.  Templates are consumed by the retrosynthetic
    planner and by the forward-synthesis scoring engine.

    Attributes
    ----------
    name : str
        Short human-readable name for the transformation.
    named_reaction : str | None
        Conventional named-reaction label (e.g. "Suzuki coupling").
    category : ReactionCategory
        Broad classification of the reaction type.
    reagents : list[str]
        Required reagents (stoichiometric).
    solvents : list[str]
        Suitable solvents.
    catalysts : list[str]
        Catalysts or co-catalysts (may be empty).
    temperature_range : tuple[float, float]
        Typical temperature window in degrees Celsius.
    typical_yield : tuple[float, float]
        Expected isolated-yield range in percent.
    functional_group_required : list[str]
        Functional groups that must be present on the substrate.
    functional_group_produced : list[str]
        Functional groups introduced or revealed by the reaction.
    functional_group_incompatible : list[str]
        Functional groups that would be destroyed or interfere.
    mechanism : str
        One- or two-sentence description of the mechanism.
    reverse_transform : str
        Description of the retrosynthetic disconnection.
    scale_notes : str
        Practical notes on scaling the reaction.
    safety_notes : str
        Hazard or handling warnings.
    """

    name: str
    named_reaction: str | None
    category: ReactionCategory
    reagents: list[str]
    solvents: list[str]
    catalysts: list[str] = field(default_factory=list)
    temperature_range: tuple[float, float] = (20.0, 25.0)
    typical_yield: tuple[float, float] = (60.0, 90.0)
    functional_group_required: list[str] = field(default_factory=list)
    functional_group_produced: list[str] = field(default_factory=list)
    functional_group_incompatible: list[str] = field(default_factory=list)
    mechanism: str = ""
    reverse_transform: str = ""
    scale_notes: str = ""
    safety_notes: str = ""

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------

    def requires(self, fg_name: str) -> bool:
        """Return True if *fg_name* is among the required functional groups."""
        return fg_name.lower() in [f.lower() for f in self.functional_group_required]

    def produces(self, fg_name: str) -> bool:
        """Return True if *fg_name* is among the produced functional groups."""
        return fg_name.lower() in [f.lower() for f in self.functional_group_produced]

    def is_compatible(self, fg_names: list[str]) -> bool:
        """Return True if none of *fg_names* appear in the incompatible list."""
        incompat = {f.lower() for f in self.functional_group_incompatible}
        return not any(fg.lower() in incompat for fg in fg_names)

    def __repr__(self) -> str:
        return f"ReactionTemplate({self.name}, {self.category.name})"
