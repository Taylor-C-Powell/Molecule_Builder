"""Reaction mechanism data model and predefined templates.

Provides an enumeration of mechanism types, a dataclass hierarchy for
describing mechanism stages with electron flows, and factory functions
for common organic reaction mechanisms (SN2, E2, radical substitution,
nucleophilic addition).

Each mechanism is decomposed into sequential stages, where each stage
specifies target distances, bond order changes, angle targets, and
electron flow arrows for visualization.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto


# ===================================================================
# Enumerations
# ===================================================================

class MechanismType(Enum):
    """Classification of reaction mechanism types."""
    SN2 = auto()
    SN1 = auto()
    E1 = auto()
    E2 = auto()
    ELECTROPHILIC_ADDITION = auto()
    NUCLEOPHILIC_ADDITION = auto()
    RADICAL_SUBSTITUTION = auto()
    CONCERTED_PERICYCLIC = auto()


class FlowType(Enum):
    """Type of electron flow arrow."""
    CURLY_ARROW = auto()      # Two-electron (heterolytic)
    FISHHOOK_ARROW = auto()   # One-electron (homolytic / radical)


# ===================================================================
# Data classes
# ===================================================================

@dataclass
class ElectronFlow:
    """An electron flow arrow for mechanism visualization.

    Attributes
    ----------
    from_atom : int
        Atom index where the electron pair originates (lone pair or
        bond).
    to_bond : tuple[int, int]
        Atom pair where electrons flow to (new bond being formed).
    flow_type : FlowType
        CURLY_ARROW for 2-electron, FISHHOOK_ARROW for 1-electron.
    label : str
        Optional annotation for the arrow.
    """
    from_atom: int
    to_bond: tuple[int, int]
    flow_type: FlowType = FlowType.CURLY_ARROW
    label: str = ""


@dataclass
class MechanismStage:
    """A single stage of a reaction mechanism.

    Attributes
    ----------
    name : str
        Human-readable name (e.g. "Backside attack", "Leaving group
        departure").
    distance_targets : dict[tuple[int, int], float]
        Target interatomic distances in Angstroms for this stage.
        Keys are (atom_i, atom_j) pairs.
    bond_order_changes : dict[tuple[int, int], float]
        Target fractional bond orders at the end of this stage.
    angle_targets : dict[tuple[int, int, int], float]
        Target bond angles in degrees (i-j-k triples).
    electron_flows : list[ElectronFlow]
        Electron flow arrows active during this stage.
    duration_weight : float
        Relative duration of this stage (for timing).  Stages with
        higher weight take proportionally longer in the animation.
    annotation : str
        Text annotation displayed during this stage.
    """
    name: str
    distance_targets: dict[tuple[int, int], float] = field(default_factory=dict)
    bond_order_changes: dict[tuple[int, int], float] = field(default_factory=dict)
    angle_targets: dict[tuple[int, int, int], float] = field(default_factory=dict)
    electron_flows: list[ElectronFlow] = field(default_factory=list)
    duration_weight: float = 1.0
    annotation: str = ""


@dataclass
class ReactionMechanism:
    """A complete reaction mechanism composed of sequential stages.

    Attributes
    ----------
    name : str
        Mechanism name (e.g. "SN2 of CH3Cl + OH-").
    mechanism_type : MechanismType
        Classification of the mechanism.
    stages : list[MechanismStage]
        Ordered list of stages.
    atom_roles : dict[str, int]
        Named roles mapped to atom indices (e.g. "nucleophile" -> 5,
        "leaving_group" -> 3).
    """
    name: str
    mechanism_type: MechanismType
    stages: list[MechanismStage] = field(default_factory=list)
    atom_roles: dict[str, int] = field(default_factory=dict)


# ===================================================================
# Predefined mechanism template factories
# ===================================================================

def sn2_mechanism(substrate_C: int,
                  nucleophile: int,
                  leaving_group: int,
                  nuc_approach_dist: float = 3.5,
                  nuc_bond_dist: float = 1.5,
                  lg_depart_dist: float = 3.5) -> ReactionMechanism:
    """Create an SN2 mechanism template.

    The SN2 mechanism proceeds in a single concerted step:
    nucleophile attacks the electrophilic carbon from the backside
    while the leaving group departs, with inversion of configuration.

    Parameters
    ----------
    substrate_C : int
        Index of the electrophilic carbon.
    nucleophile : int
        Index of the nucleophile atom.
    leaving_group : int
        Index of the leaving group atom.
    nuc_approach_dist : float
        Initial Nu...C distance to start from.
    nuc_bond_dist : float
        Final Nu-C bond distance.
    lg_depart_dist : float
        Final C...LG distance.

    Returns
    -------
    ReactionMechanism
    """
    nuc_key = (min(nucleophile, substrate_C), max(nucleophile, substrate_C))
    lg_key = (min(leaving_group, substrate_C), max(leaving_group, substrate_C))

    stages = [
        # Stage 0: Approach
        MechanismStage(
            name="Nucleophile approach",
            distance_targets={
                (nucleophile, substrate_C): nuc_approach_dist * 0.7,
            },
            bond_order_changes={
                nuc_key: 0.2,
                lg_key: 0.9,
            },
            electron_flows=[
                ElectronFlow(
                    from_atom=nucleophile,
                    to_bond=(nucleophile, substrate_C),
                    flow_type=FlowType.CURLY_ARROW,
                    label="Nu: attacks C",
                ),
            ],
            duration_weight=1.0,
            annotation="Nucleophile approaches electrophilic carbon",
        ),
        # Stage 1: Transition state
        MechanismStage(
            name="Transition state [Nu...C...LG]",
            distance_targets={
                (nucleophile, substrate_C): (nuc_bond_dist + nuc_approach_dist * 0.7) / 2,
                (leaving_group, substrate_C): (nuc_bond_dist + lg_depart_dist) / 2,
            },
            bond_order_changes={
                nuc_key: 0.5,
                lg_key: 0.5,
            },
            electron_flows=[
                ElectronFlow(
                    from_atom=nucleophile,
                    to_bond=(nucleophile, substrate_C),
                    flow_type=FlowType.CURLY_ARROW,
                    label="Partial bond forming",
                ),
                ElectronFlow(
                    from_atom=substrate_C,
                    to_bond=(substrate_C, leaving_group),
                    flow_type=FlowType.CURLY_ARROW,
                    label="Bond breaking",
                ),
            ],
            duration_weight=1.5,
            annotation="Transition state: pentacoordinate carbon",
        ),
        # Stage 2: Product formation
        MechanismStage(
            name="Product formation",
            distance_targets={
                (nucleophile, substrate_C): nuc_bond_dist,
                (leaving_group, substrate_C): lg_depart_dist,
            },
            bond_order_changes={
                nuc_key: 1.0,
                lg_key: 0.0,
            },
            electron_flows=[
                ElectronFlow(
                    from_atom=substrate_C,
                    to_bond=(substrate_C, leaving_group),
                    flow_type=FlowType.CURLY_ARROW,
                    label="LG departs with electrons",
                ),
            ],
            duration_weight=1.0,
            annotation="Leaving group departs; Walden inversion complete",
        ),
    ]

    return ReactionMechanism(
        name="SN2 mechanism",
        mechanism_type=MechanismType.SN2,
        stages=stages,
        atom_roles={
            "substrate_C": substrate_C,
            "nucleophile": nucleophile,
            "leaving_group": leaving_group,
        },
    )


def e2_mechanism(alpha_C: int,
                 beta_H: int,
                 base: int,
                 leaving_group: int,
                 base_approach_dist: float = 3.0,
                 base_H_dist: float = 1.0,
                 lg_depart_dist: float = 3.5) -> ReactionMechanism:
    """Create an E2 elimination mechanism template.

    The E2 mechanism is concerted: the base abstracts the beta-hydrogen
    while the leaving group departs, forming a double bond.

    Parameters
    ----------
    alpha_C : int
        Index of the alpha carbon (bearing the leaving group).
    beta_H : int
        Index of the beta hydrogen.
    base : int
        Index of the base atom.
    leaving_group : int
        Index of the leaving group atom.

    Returns
    -------
    ReactionMechanism
    """
    base_H_key = (min(base, beta_H), max(base, beta_H))
    lg_key = (min(leaving_group, alpha_C), max(leaving_group, alpha_C))

    stages = [
        MechanismStage(
            name="Base approach",
            distance_targets={
                (base, beta_H): base_approach_dist * 0.6,
            },
            bond_order_changes={
                base_H_key: 0.3,
                lg_key: 0.9,
            },
            electron_flows=[
                ElectronFlow(
                    from_atom=base,
                    to_bond=(base, beta_H),
                    flow_type=FlowType.CURLY_ARROW,
                    label="Base abstracts H",
                ),
            ],
            duration_weight=1.0,
            annotation="Base approaches beta-hydrogen",
        ),
        MechanismStage(
            name="Concerted E2 transition state",
            distance_targets={
                (base, beta_H): base_H_dist * 1.3,
                (leaving_group, alpha_C): lg_depart_dist * 0.5,
            },
            bond_order_changes={
                base_H_key: 0.7,
                lg_key: 0.3,
            },
            electron_flows=[
                ElectronFlow(
                    from_atom=base,
                    to_bond=(base, beta_H),
                    flow_type=FlowType.CURLY_ARROW,
                    label="H abstraction",
                ),
                ElectronFlow(
                    from_atom=alpha_C,
                    to_bond=(alpha_C, leaving_group),
                    flow_type=FlowType.CURLY_ARROW,
                    label="LG departure",
                ),
            ],
            duration_weight=1.5,
            annotation="Transition state: anti-periplanar elimination",
        ),
        MechanismStage(
            name="Alkene product",
            distance_targets={
                (base, beta_H): base_H_dist,
                (leaving_group, alpha_C): lg_depart_dist,
            },
            bond_order_changes={
                base_H_key: 1.0,
                lg_key: 0.0,
            },
            electron_flows=[],
            duration_weight=1.0,
            annotation="Double bond formed; leaving group departed",
        ),
    ]

    return ReactionMechanism(
        name="E2 elimination",
        mechanism_type=MechanismType.E2,
        stages=stages,
        atom_roles={
            "alpha_C": alpha_C,
            "beta_H": beta_H,
            "base": base,
            "leaving_group": leaving_group,
        },
    )


def radical_substitution_mechanism(
        target_C: int,
        target_H: int,
        radical: int,
        rad_approach_dist: float = 3.0,
        rad_H_dist: float = 1.1,
        new_bond_dist: float = 1.5) -> ReactionMechanism:
    """Create a radical substitution (propagation) mechanism.

    Models the hydrogen abstraction step in radical chain reactions:
    R* + H-C -> R-H + C*

    Parameters
    ----------
    target_C : int
        Index of the carbon losing the hydrogen.
    target_H : int
        Index of the hydrogen being abstracted.
    radical : int
        Index of the radical atom.

    Returns
    -------
    ReactionMechanism
    """
    rad_H_key = (min(radical, target_H), max(radical, target_H))
    ch_key = (min(target_C, target_H), max(target_C, target_H))

    stages = [
        MechanismStage(
            name="Radical approach",
            distance_targets={
                (radical, target_H): rad_approach_dist * 0.6,
            },
            bond_order_changes={
                rad_H_key: 0.2,
                ch_key: 0.9,
            },
            electron_flows=[
                ElectronFlow(
                    from_atom=radical,
                    to_bond=(radical, target_H),
                    flow_type=FlowType.FISHHOOK_ARROW,
                    label="Radical attacks H",
                ),
            ],
            duration_weight=1.0,
            annotation="Radical approaches C-H bond",
        ),
        MechanismStage(
            name="Transition state",
            distance_targets={
                (radical, target_H): (rad_H_dist + rad_approach_dist * 0.6) / 2,
            },
            bond_order_changes={
                rad_H_key: 0.5,
                ch_key: 0.5,
            },
            electron_flows=[
                ElectronFlow(
                    from_atom=radical,
                    to_bond=(radical, target_H),
                    flow_type=FlowType.FISHHOOK_ARROW,
                    label="1e- transfer",
                ),
                ElectronFlow(
                    from_atom=target_C,
                    to_bond=(target_C, target_H),
                    flow_type=FlowType.FISHHOOK_ARROW,
                    label="1e- transfer",
                ),
            ],
            duration_weight=1.5,
            annotation="Transition state: R...H...C",
        ),
        MechanismStage(
            name="Product (new radical)",
            distance_targets={
                (radical, target_H): rad_H_dist,
            },
            bond_order_changes={
                rad_H_key: 1.0,
                ch_key: 0.0,
            },
            electron_flows=[],
            duration_weight=1.0,
            annotation="H transferred; carbon radical formed",
        ),
    ]

    return ReactionMechanism(
        name="Radical substitution (propagation)",
        mechanism_type=MechanismType.RADICAL_SUBSTITUTION,
        stages=stages,
        atom_roles={
            "target_C": target_C,
            "target_H": target_H,
            "radical": radical,
        },
    )


def nucleophilic_addition_mechanism(
        carbonyl_C: int,
        carbonyl_O: int,
        nucleophile: int,
        nuc_approach_dist: float = 3.5,
        nuc_bond_dist: float = 1.5) -> ReactionMechanism:
    """Create a nucleophilic addition to carbonyl mechanism.

    Nucleophile attacks the carbonyl carbon, converting C=O to C-O^-.

    Parameters
    ----------
    carbonyl_C : int
        Index of the carbonyl carbon.
    carbonyl_O : int
        Index of the carbonyl oxygen.
    nucleophile : int
        Index of the nucleophile atom.

    Returns
    -------
    ReactionMechanism
    """
    nuc_C_key = (min(nucleophile, carbonyl_C), max(nucleophile, carbonyl_C))
    co_key = (min(carbonyl_C, carbonyl_O), max(carbonyl_C, carbonyl_O))

    stages = [
        MechanismStage(
            name="Nucleophile approach",
            distance_targets={
                (nucleophile, carbonyl_C): nuc_approach_dist * 0.6,
            },
            bond_order_changes={
                nuc_C_key: 0.2,
                co_key: 1.8,
            },
            electron_flows=[
                ElectronFlow(
                    from_atom=nucleophile,
                    to_bond=(nucleophile, carbonyl_C),
                    flow_type=FlowType.CURLY_ARROW,
                    label="Nu: attacks C=O",
                ),
            ],
            duration_weight=1.0,
            annotation="Nucleophile approaches electrophilic carbonyl C",
        ),
        MechanismStage(
            name="Tetrahedral intermediate forming",
            distance_targets={
                (nucleophile, carbonyl_C): (nuc_bond_dist + nuc_approach_dist * 0.6) / 2,
            },
            bond_order_changes={
                nuc_C_key: 0.6,
                co_key: 1.4,
            },
            electron_flows=[
                ElectronFlow(
                    from_atom=nucleophile,
                    to_bond=(nucleophile, carbonyl_C),
                    flow_type=FlowType.CURLY_ARROW,
                    label="Bond forming",
                ),
                ElectronFlow(
                    from_atom=carbonyl_C,
                    to_bond=(carbonyl_C, carbonyl_O),
                    flow_type=FlowType.CURLY_ARROW,
                    label="Pi electrons to O",
                ),
            ],
            duration_weight=1.5,
            annotation="Rehybridization: sp2 to sp3",
        ),
        MechanismStage(
            name="Tetrahedral product",
            distance_targets={
                (nucleophile, carbonyl_C): nuc_bond_dist,
            },
            bond_order_changes={
                nuc_C_key: 1.0,
                co_key: 1.0,
            },
            electron_flows=[],
            duration_weight=1.0,
            annotation="Tetrahedral alkoxide intermediate formed",
        ),
    ]

    return ReactionMechanism(
        name="Nucleophilic addition to carbonyl",
        mechanism_type=MechanismType.NUCLEOPHILIC_ADDITION,
        stages=stages,
        atom_roles={
            "carbonyl_C": carbonyl_C,
            "carbonyl_O": carbonyl_O,
            "nucleophile": nucleophile,
        },
    )
