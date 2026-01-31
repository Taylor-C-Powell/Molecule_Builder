"""Retrosynthetic analysis engine using beam search.

Given a target molecule, this module works backwards from the product to
identify commercially available starting materials, applying known reaction
templates in reverse (disconnection approach).  A beam search explores the
most promising disconnections at each level, producing a retrosynthesis
tree that can later be converted into a forward synthesis route.

Key public function
-------------------
retrosynthesis(mol, max_depth, beam_width) -> RetrosynthesisTree

Supporting helpers
------------------
is_purchasable(smiles) -> bool
get_purchasable(smiles) -> Precursor | None
score_disconnection(template, precursors, target_mol) -> float
format_tree(tree) -> str
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from molbuilder.molecule.graph import Molecule
from molbuilder.smiles.parser import parse
from molbuilder.smiles.writer import to_smiles
from molbuilder.reactions.reaction_types import ReactionTemplate, ReactionCategory
from molbuilder.reactions.knowledge_base import (
    REACTION_TEMPLATES,
    lookup_by_functional_group,
    find_reactions_producing,
)
from molbuilder.reactions.functional_group_detect import (
    detect_functional_groups,
    FunctionalGroup,
)


# =====================================================================
#  Purchasable starting materials database  (~200 entries)
# =====================================================================

# Each entry maps a canonical SMILES to (common_name, cost_per_kg_usd).
# Organised roughly by functional-group class so the table is easy to
# extend.  Costs are representative order-of-magnitude estimates for
# bulk laboratory quantities and are NOT authoritative pricing data.

PURCHASABLE_MATERIALS: dict[str, tuple[str, float]] = {
    # --- simple hydrocarbons / gases ---
    "C": ("methane", 0.50),
    "CC": ("ethane", 1.00),
    "CCC": ("propane", 1.20),
    "CCCC": ("n-butane", 1.50),
    "C=C": ("ethylene", 1.50),
    "CC=C": ("propylene", 2.00),
    "C=CC=C": ("1,3-butadiene", 3.00),
    "C#C": ("acetylene", 2.50),
    "CC#C": ("propyne", 5.00),
    "C1CC1": ("cyclopropane", 8.00),
    "C1CCC1": ("cyclobutane", 15.00),
    "C1CCCC1": ("cyclopentane", 4.00),
    "C1CCCCC1": ("cyclohexane", 3.00),

    # --- alkyl halides ---
    "CCl": ("chloromethane", 2.00),
    "CBr": ("bromomethane", 5.00),
    "CI": ("iodomethane", 12.00),
    "CCCl": ("chloroethane", 3.00),
    "CCBr": ("bromoethane", 6.00),
    "CCI": ("iodoethane", 15.00),
    "CCCCl": ("1-chloropropane", 4.00),
    "CCCBr": ("1-bromopropane", 7.00),
    "CCCI": ("1-iodopropane", 18.00),
    "CCCCCl": ("1-chlorobutane", 5.00),
    "CCCCBr": ("1-bromobutane", 8.00),
    "CCCCI": ("1-iodobutane", 20.00),
    "CC(C)Cl": ("2-chloropropane", 5.00),
    "CC(C)Br": ("2-bromopropane", 8.00),
    "CC(C)(C)Cl": ("tert-butyl chloride", 6.00),
    "CC(C)(C)Br": ("tert-butyl bromide", 10.00),
    "C(Cl)(Cl)Cl": ("chloroform", 2.50),
    "C(Cl)Cl": ("dichloromethane", 2.00),
    "ClC=C": ("vinyl chloride", 2.00),
    "BrC=C": ("vinyl bromide", 8.00),
    "ClCC=C": ("allyl chloride", 4.00),
    "BrCC=C": ("allyl bromide", 7.00),

    # --- alcohols ---
    "CO": ("methanol", 1.00),
    "CCO": ("ethanol", 2.00),
    "CCCO": ("1-propanol", 3.00),
    "CC(C)O": ("2-propanol", 2.50),
    "CCCCO": ("1-butanol", 3.50),
    "CC(C)(C)O": ("tert-butanol", 4.00),
    "CCCCCO": ("1-pentanol", 5.00),
    "CCCCCCO": ("1-hexanol", 6.00),
    "OCC=C": ("allyl alcohol", 5.00),
    "OC1CCCCC1": ("cyclohexanol", 5.00),
    "OCCO": ("ethylene glycol", 2.00),
    "OCC(O)CO": ("glycerol", 2.50),
    "OC(C)(C)C": ("tert-butanol (alt)", 4.00),

    # --- water and simple inorganics ---
    "O": ("water", 0.01),
    "[NH3]": ("ammonia", 0.80),
    "N": ("ammonia (SMILES variant)", 0.80),
    "Cl": ("hydrochloric acid", 0.50),
    "O=C=O": ("carbon dioxide", 0.30),
    "S": ("hydrogen sulfide", 1.00),

    # --- aldehydes ---
    "C=O": ("formaldehyde", 1.50),
    "CC=O": ("acetaldehyde", 3.00),
    "CCC=O": ("propanal", 5.00),
    "CCCC=O": ("butanal", 6.00),
    "CCCCC=O": ("pentanal", 8.00),
    "O=CC=O": ("glyoxal", 5.00),

    # --- ketones ---
    "CC(C)=O": ("acetone", 1.50),
    "CCC(C)=O": ("methyl ethyl ketone", 3.00),
    "CCC(CC)=O": ("3-pentanone", 5.00),
    "CCCC(C)=O": ("2-pentanone", 5.00),
    "O=C1CCCCC1": ("cyclohexanone", 4.00),
    "C=CC(C)=O": ("methyl vinyl ketone", 6.00),

    # --- carboxylic acids ---
    "OC=O": ("formic acid", 2.00),
    "CC(O)=O": ("acetic acid", 1.50),
    "CCC(O)=O": ("propionic acid", 3.00),
    "CCCC(O)=O": ("butyric acid", 4.00),
    "CCCCC(O)=O": ("valeric acid", 6.00),
    "OC(=O)C=C": ("acrylic acid", 3.00),
    "OC(=O)CC(O)=O": ("malonic acid", 5.00),
    "OC(=O)CCC(O)=O": ("succinic acid", 4.00),
    "OC(=O)CCCCC(O)=O": ("adipic acid", 4.50),
    "OC(=O)C(O)=O": ("oxalic acid", 3.50),

    # --- esters ---
    "COC(C)=O": ("methyl acetate", 3.00),
    "CCOC(C)=O": ("ethyl acetate", 2.50),
    "CCOC(=O)CC": ("ethyl propanoate", 4.00),
    "CCOC(=O)OCC": ("diethyl carbonate", 5.00),

    # --- ethers ---
    "COC": ("dimethyl ether", 2.00),
    "CCOCC": ("diethyl ether", 3.00),
    "C1CCOC1": ("tetrahydrofuran", 4.00),
    "C1COCCO1": ("1,4-dioxane", 5.00),
    "COC=C": ("methyl vinyl ether", 6.00),
    "COCCOCCOCC": ("diglyme", 8.00),

    # --- amines ---
    "CN": ("methylamine", 3.00),
    "CCN": ("ethylamine", 4.00),
    "CCCN": ("propylamine", 5.00),
    "CCCCN": ("butylamine", 6.00),
    "CNC": ("dimethylamine", 4.00),
    "CN(C)C": ("trimethylamine", 5.00),
    "CCN(CC)CC": ("triethylamine", 6.00),
    "NCC=C": ("allylamine", 7.00),
    "NC1CCCCC1": ("cyclohexylamine", 8.00),
    "NCCN": ("ethylenediamine", 5.00),
    "NCCCN": ("1,3-diaminopropane", 7.00),
    "NCCCCN": ("1,4-diaminobutane", 8.00),

    # --- amides ---
    "NC=O": ("formamide", 3.00),
    "CC(N)=O": ("acetamide", 4.00),
    "CN(C)C=O": ("dimethylformamide", 3.50),

    # --- nitriles ---
    "C#N": ("hydrogen cyanide", 2.00),
    "CC#N": ("acetonitrile", 3.00),
    "CCC#N": ("propionitrile", 5.00),
    "CCCC#N": ("butyronitrile", 7.00),

    # --- aromatics ---
    "c1ccccc1": ("benzene", 2.50),
    "Cc1ccccc1": ("toluene", 2.50),
    "CCc1ccccc1": ("ethylbenzene", 3.50),
    "C=Cc1ccccc1": ("styrene", 4.00),
    "c1ccc(cc1)C": ("toluene (alt)", 2.50),
    "Oc1ccccc1": ("phenol", 3.00),
    "Nc1ccccc1": ("aniline", 4.00),
    "Clc1ccccc1": ("chlorobenzene", 3.50),
    "Brc1ccccc1": ("bromobenzene", 5.00),
    "Ic1ccccc1": ("iodobenzene", 10.00),
    "OC(=O)c1ccccc1": ("benzoic acid", 3.50),
    "O=Cc1ccccc1": ("benzaldehyde", 5.00),
    "CC(=O)c1ccccc1": ("acetophenone", 5.00),
    "c1ccc2ccccc2c1": ("naphthalene", 4.00),
    "c1ccncc1": ("pyridine", 4.00),
    "C1=COC=C1": ("furan", 5.00),
    "c1cc[nH]c1": ("pyrrole", 6.00),

    # --- aromatic halides ---
    "Fc1ccccc1": ("fluorobenzene", 6.00),
    "Clc1ccc(Cl)cc1": ("1,4-dichlorobenzene", 4.00),

    # --- amino acids (common L-forms, simplified SMILES) ---
    "NCC(O)=O": ("glycine", 5.00),
    "CC(N)C(O)=O": ("alanine", 8.00),
    "CC(C)C(N)C(O)=O": ("valine", 15.00),
    "CC(CC)C(N)C(O)=O": ("isoleucine", 20.00),
    "CCCC(N)C(O)=O": ("leucine (linear approx)", 18.00),
    "NC(=O)CC(N)C(O)=O": ("asparagine", 15.00),
    "OC(=O)CC(N)C(O)=O": ("aspartic acid", 12.00),
    "OC(=O)CCC(N)C(O)=O": ("glutamic acid", 12.00),
    "NCCCCC(N)C(O)=O": ("lysine", 20.00),
    "NC(N)=NCCCC(N)C(O)=O": ("arginine", 25.00),

    # --- thiols ---
    "CS": ("methanethiol", 4.00),
    "CCS": ("ethanethiol", 5.00),
    "CCCS": ("1-propanethiol", 7.00),

    # --- acid chlorides ---
    "CC(Cl)=O": ("acetyl chloride", 4.00),
    "ClC(Cl)=O": ("phosgene", 3.00),
    "CCC(Cl)=O": ("propanoyl chloride", 6.00),
    "OC(Cl)=O": ("chloroformic acid", 5.00),

    # --- acid anhydrides ---
    "CC(=O)OC(C)=O": ("acetic anhydride", 3.00),

    # --- epoxides ---
    "C1CO1": ("ethylene oxide", 3.00),
    "CC1CO1": ("propylene oxide", 4.00),

    # --- miscellaneous building blocks ---
    "C(=O)O": ("formic acid (alt)", 2.00),
    "CCCCCCCCCCCC": ("dodecane", 4.00),
    "CCCCCCCC": ("octane", 3.00),
    "CCCCCC": ("hexane", 2.50),
    "CCCCC": ("pentane", 2.00),
    "CC(C)CC": ("isopentane", 2.50),
    "CC(C)C": ("isobutane", 2.00),
    "C=CC(=O)OC": ("methyl acrylate", 4.00),
    "C=CC(=O)OCC": ("ethyl acrylate", 5.00),
    "C=C(C)C(=O)OC": ("methyl methacrylate", 5.00),
    "C(CO)O": ("ethylene glycol (alt)", 2.00),
    "OCCCCO": ("1,4-butanediol", 4.00),
    "C(F)(F)F": ("fluoroform", 3.00),
    "C(Cl)(Cl)(Cl)Cl": ("carbon tetrachloride", 3.00),
    "C(F)(F)(F)Cl": ("chlorotrifluoromethane", 5.00),
    "CC(=O)OC=C": ("vinyl acetate", 4.00),
    "ClCCCl": ("1,2-dichloroethane", 2.00),
    "BrCCBr": ("1,2-dibromoethane", 6.00),
    "CCCCCCCCCCCCCCCCCC(O)=O": ("stearic acid", 4.00),
    "CCCCCCCC(O)=O": ("octanoic acid", 5.00),

    # --- sugars / polyols ---
    "OCC(O)C(O)C(O)C(O)CO": ("D-sorbitol", 3.50),
    "OCC(O)C(O)CO": ("erythritol", 6.00),

    # --- diacids / anhydrides ---
    "O=C1OC(=O)C=C1": ("maleic anhydride", 3.50),
    "OC(=O)C=CC(O)=O": ("maleic acid", 4.00),

    # --- phosphorus / sulfur reagents (simplified) ---
    "OP(O)(O)=O": ("phosphoric acid", 1.50),
    "OS(O)(=O)=O": ("sulfuric acid", 0.50),
    "OS(=O)=O": ("sulfurous acid", 2.00),

    # --- azides and nitro compounds ---
    "CN=[N+]=[N-]": ("methyl azide", 10.00),
    "C[N+](=O)[O-]": ("nitromethane", 5.00),
    "CC[N+](=O)[O-]": ("nitroethane", 7.00),
    "[O-][N+](=O)c1ccccc1": ("nitrobenzene", 5.00),

    # --- additional alcohols ---
    "CC(O)CC": ("2-butanol", 4.00),
    "C(CO)(CO)CO": ("pentaerythritol", 6.00),
    "OC(C)C": ("2-propanol (alt)", 2.50),

    # --- additional halides ---
    "FC(F)F": ("trifluoromethane", 5.00),
    "C(F)(F)(F)C(F)(F)F": ("hexafluoroethane", 8.00),
    "ClC(Cl)=C": ("vinylidene chloride", 4.00),
    "CC(Cl)(C)C": ("neopentyl chloride", 7.00),
    "BrCCCBr": ("1,3-dibromopropane", 8.00),
    "BrCCCCBr": ("1,4-dibromobutane", 10.00),
    "ICCl": ("chloroiodomethane", 12.00),

    # --- additional aromatics ---
    "c1ccoc1": ("furan (aromatic)", 5.00),
    "c1ccsc1": ("thiophene", 5.00),
    "c1cnc2ccccc2c1": ("quinoline", 8.00),
    "c1ccc2c(c1)cccc2": ("naphthalene (alt)", 4.00),
    "OCc1ccccc1": ("benzyl alcohol", 5.00),
    "NCc1ccccc1": ("benzylamine", 7.00),
    "ClCc1ccccc1": ("benzyl chloride", 6.00),
    "BrCc1ccccc1": ("benzyl bromide", 8.00),
    "c1ccc(O)c(O)c1": ("catechol", 6.00),
    "c1cc(O)cc(O)c1": ("resorcinol", 7.00),
    "Oc1ccc(O)cc1": ("hydroquinone", 5.00),
    "CC(=O)Oc1ccccc1": ("phenyl acetate", 6.00),

    # --- heterocycles ---
    "C1CCNCC1": ("piperidine", 5.00),
    "C1CCNC1": ("pyrrolidine", 6.00),
    "C1CCOC1": ("tetrahydrofuran (ring)", 4.00),
    "C1CCOCC1": ("tetrahydropyran", 5.00),
    "C1CNCCN1": ("piperazine", 6.00),
    "c1c[nH]cn1": ("imidazole", 7.00),
    "C1CO1": ("ethylene oxide (ring)", 3.00),

    # --- additional carboxylic acid derivatives ---
    "CC(=O)NC": ("N-methylacetamide", 5.00),
    "O=C(Cl)c1ccccc1": ("benzoyl chloride", 7.00),
    "OC(=O)CCCCCC(O)=O": ("pimelic acid", 6.00),
    "OC(=O)c1ccc(C(O)=O)cc1": ("terephthalic acid", 5.00),
}

# Also accept alkyl halide generic name
_PURCHASABLE_ALIASES: dict[str, str] = {}


# =====================================================================
#  Data structures
# =====================================================================

@dataclass
class Precursor:
    """A molecule that serves as starting material for one reaction step.

    Attributes
    ----------
    smiles : str
        SMILES representation of the precursor.
    molecule : Molecule | None
        Parsed molecule object, or None if looked up from purchasable DB.
    name : str
        Human-readable name for display.
    cost_per_kg : float
        Estimated cost per kilogram in USD.
    """
    smiles: str
    molecule: Molecule | None
    name: str
    cost_per_kg: float


@dataclass
class Disconnection:
    """One possible retrosynthetic disconnection for a target node.

    Attributes
    ----------
    template : ReactionTemplate
        The reaction template applied in reverse.
    precursors : list[Precursor]
        The precursor molecules produced by this disconnection.
    score : float
        Quality score from 0 to 100 (higher is better).
    """
    template: ReactionTemplate
    precursors: list[Precursor]
    score: float


@dataclass
class RetroNode:
    """A node in the retrosynthetic search tree.

    Each node represents one molecule.  If the molecule is not purchasable,
    its ``disconnections`` list holds candidate retrosynthetic steps, and
    ``children`` holds the recursively expanded precursor nodes for the
    best disconnection.

    Attributes
    ----------
    smiles : str
        SMILES string for this molecule.
    molecule : Molecule
        The parsed Molecule object.
    functional_groups : list[FunctionalGroup]
        Functional groups detected on this molecule.
    is_purchasable : bool
        True if this molecule appears in PURCHASABLE_MATERIALS.
    disconnections : list[Disconnection]
        Candidate retrosynthetic disconnections (best first).
    best_disconnection : Disconnection | None
        The top-scoring disconnection, if any.
    children : list[RetroNode]
        Expanded child nodes (precursors of the best disconnection).
    depth : int
        Depth of this node in the search tree (root = 0).
    """
    smiles: str
    molecule: Molecule
    functional_groups: list[FunctionalGroup] = field(default_factory=list)
    is_purchasable: bool = False
    disconnections: list[Disconnection] = field(default_factory=list)
    best_disconnection: Disconnection | None = None
    children: list["RetroNode"] = field(default_factory=list)
    depth: int = 0


@dataclass
class RetrosynthesisTree:
    """Complete retrosynthesis result.

    Attributes
    ----------
    target : RetroNode
        Root of the retrosynthesis tree (the target molecule).
    max_depth : int
        Maximum search depth used.
    beam_width : int
        Beam width (number of disconnections kept per level).
    routes_found : int
        Total number of complete routes found to purchasable materials.
    """
    target: RetroNode
    max_depth: int
    beam_width: int
    routes_found: int


# =====================================================================
#  Purchasability checks
# =====================================================================

def is_purchasable(smiles: str) -> bool:
    """Return True if *smiles* matches a known purchasable material.

    The check tries the SMILES string as-is and also attempts a round-
    trip (parse then re-serialise) to handle minor notational differences.
    """
    if smiles in PURCHASABLE_MATERIALS:
        return True
    # Try canonical round-trip
    try:
        canon = to_smiles(parse(smiles))
        if canon in PURCHASABLE_MATERIALS:
            return True
    except Exception:
        pass
    return False


def get_purchasable(smiles: str) -> Precursor | None:
    """Return a Precursor for *smiles* if it is purchasable, else None."""
    entry = PURCHASABLE_MATERIALS.get(smiles)
    if entry is not None:
        name, cost = entry
        return Precursor(smiles=smiles, molecule=None, name=name,
                         cost_per_kg=cost)
    # Try canonical form
    try:
        canon = to_smiles(parse(smiles))
        entry = PURCHASABLE_MATERIALS.get(canon)
        if entry is not None:
            name, cost = entry
            return Precursor(smiles=canon, molecule=None, name=name,
                             cost_per_kg=cost)
    except Exception:
        pass
    return None


# =====================================================================
#  Scoring helpers
# =====================================================================

def _count_heavy_atoms(mol: Molecule) -> int:
    """Count non-hydrogen atoms."""
    return sum(1 for a in mol.atoms if a.symbol != "H")


def _heavy_atom_count_from_smiles(smiles: str) -> int:
    """Count heavy atoms by parsing SMILES (returns 0 on failure)."""
    try:
        mol = parse(smiles)
        return _count_heavy_atoms(mol)
    except Exception:
        return 0


def score_disconnection(
    template: ReactionTemplate,
    precursors: list[Precursor],
    target_mol: Molecule,
) -> float:
    """Score a retrosynthetic disconnection from 0 (poor) to 100 (ideal).

    The score is a weighted sum of several heuristic factors:

    1. **Yield expectation** (0--25 pts):  Higher template yield is better.
    2. **Precursor availability** (0--30 pts):  More purchasable precursors
       contribute more points.
    3. **Complexity reduction** (0--20 pts):  If precursors are significantly
       simpler (fewer heavy atoms) than the target, the score is higher.
    4. **Strategic bond preference** (0--15 pts):  Reactions that form C-C
       bonds (coupling, Grignard, aldol, etc.) score highest because C-C
       disconnections are the backbone of retrosynthetic strategy.
    5. **Template category bonus** (0--10 pts):  Coupling and carbonyl
       reactions get a small bonus as they are the most commonly used
       strategic transforms.
    """
    score = 0.0

    # --- 1. Yield expectation (0--25) ---
    lo, hi = template.typical_yield
    mid_yield = (lo + hi) / 2.0
    score += 25.0 * (mid_yield / 100.0)

    # --- 2. Precursor availability (0--30) ---
    if precursors:
        purchasable_count = sum(1 for p in precursors if is_purchasable(p.smiles))
        frac = purchasable_count / len(precursors)
        score += 30.0 * frac

    # --- 3. Complexity reduction (0--20) ---
    target_heavy = _count_heavy_atoms(target_mol)
    if target_heavy > 0 and precursors:
        max_precursor_heavy = max(
            _heavy_atom_count_from_smiles(p.smiles) for p in precursors
        )
        if max_precursor_heavy < target_heavy:
            reduction = (target_heavy - max_precursor_heavy) / target_heavy
            score += 20.0 * min(1.0, reduction * 2.0)
        # If the precursor is no simpler, no points here.

    # --- 4. Strategic bond preference (0--15) ---
    cc_keywords = ("coupling", "grignard", "aldol", "wittig", "suzuki",
                   "heck", "sonogashira", "stille", "negishi",
                   "horner", "claisen condensation", "michael",
                   "robinson")
    name_lower = template.name.lower()
    named_lower = (template.named_reaction or "").lower()
    if any(kw in name_lower or kw in named_lower for kw in cc_keywords):
        score += 15.0
    elif template.category == ReactionCategory.COUPLING:
        score += 12.0
    elif template.category in (ReactionCategory.CARBONYL,
                               ReactionCategory.ADDITION):
        score += 6.0

    # --- 5. Template category bonus (0--10) ---
    category_bonus = {
        ReactionCategory.COUPLING: 10.0,
        ReactionCategory.CARBONYL: 8.0,
        ReactionCategory.ADDITION: 6.0,
        ReactionCategory.SUBSTITUTION: 5.0,
        ReactionCategory.REDUCTION: 4.0,
        ReactionCategory.OXIDATION: 4.0,
        ReactionCategory.ELIMINATION: 3.0,
        ReactionCategory.REARRANGEMENT: 3.0,
        ReactionCategory.PROTECTION: 1.0,
        ReactionCategory.DEPROTECTION: 1.0,
    }
    score += category_bonus.get(template.category, 2.0)

    return min(100.0, max(0.0, score))


# =====================================================================
#  Reverse transform: generate precursor SMILES from a template
# =====================================================================

def _generate_precursors_for_template(
    target_smiles: str,
    target_mol: Molecule,
    template: ReactionTemplate,
    fg: FunctionalGroup,
) -> list[Precursor]:
    """Generate precursor SMILES by conceptually reversing *template*.

    The approach is a simplification: rather than performing a full
    subgraph transform, we modify the target molecule according to
    the functional group that the reaction *produces*.  The idea is to
    remove or simplify the functional group that the forward reaction
    would create, yielding one or more simpler precursor molecules.

    For multi-component reactions (e.g. Grignard, Suzuki) two precursors
    are generated by splitting the target at the bond(s) adjacent to the
    functional group centre.

    Returns a list of Precursor objects (may be empty on failure).
    """
    precursors: list[Precursor] = []

    cat = template.category
    fg_name = fg.name
    center = fg.center
    fg_atoms = fg.atoms

    # ---- Strategy: map reaction category to precursor generation ----

    # REDUCTION or OXIDATION: the precursor is the oxidised/reduced form.
    # We approximate by swapping the FG for the one the template requires.
    if cat == ReactionCategory.REDUCTION:
        # Template reduces FG_required -> FG_produced.
        # Reverse: we have the product, so precursor has the FG_required.
        # Simplification: return a variant SMILES with the bond order changed.
        precursor_smi = _modify_fg_smiles(
            target_smiles, target_mol, fg, template, direction="oxidise")
        if precursor_smi:
            precursors.append(Precursor(
                smiles=precursor_smi, molecule=None,
                name=f"precursor ({template.name})",
                cost_per_kg=_estimate_cost(precursor_smi),
            ))
        return precursors

    if cat == ReactionCategory.OXIDATION:
        precursor_smi = _modify_fg_smiles(
            target_smiles, target_mol, fg, template, direction="reduce")
        if precursor_smi:
            precursors.append(Precursor(
                smiles=precursor_smi, molecule=None,
                name=f"precursor ({template.name})",
                cost_per_kg=_estimate_cost(precursor_smi),
            ))
        return precursors

    # COUPLING / CARBONYL: split into two fragments
    if cat in (ReactionCategory.COUPLING, ReactionCategory.CARBONYL):
        frags = _split_at_fg(target_smiles, target_mol, fg, template)
        for smi in frags:
            precursors.append(Precursor(
                smiles=smi, molecule=None,
                name=f"fragment ({template.name})",
                cost_per_kg=_estimate_cost(smi),
            ))
        return precursors

    # SUBSTITUTION: replace the produced FG with the required one
    if cat == ReactionCategory.SUBSTITUTION:
        precursor_smi = _substitute_fg(
            target_smiles, target_mol, fg, template)
        if precursor_smi:
            precursors.append(Precursor(
                smiles=precursor_smi, molecule=None,
                name=f"precursor ({template.name})",
                cost_per_kg=_estimate_cost(precursor_smi),
            ))
        # Also add the reagent as a precursor if it is recognisable
        for reagent in template.reagents:
            rp = _reagent_to_precursor(reagent)
            if rp is not None:
                precursors.append(rp)
        return precursors

    # ELIMINATION / ADDITION: forward-reverse pair
    if cat == ReactionCategory.ELIMINATION:
        # The product is an alkene; precursor is an alkyl halide or alcohol.
        precursor_smi = _add_across_double_bond(
            target_smiles, target_mol, fg, template)
        if precursor_smi:
            precursors.append(Precursor(
                smiles=precursor_smi, molecule=None,
                name=f"precursor ({template.name})",
                cost_per_kg=_estimate_cost(precursor_smi),
            ))
        return precursors

    if cat == ReactionCategory.ADDITION:
        # The product has a new FG across a former double bond.
        precursor_smi = _remove_addition(
            target_smiles, target_mol, fg, template)
        if precursor_smi:
            precursors.append(Precursor(
                smiles=precursor_smi, molecule=None,
                name=f"precursor ({template.name})",
                cost_per_kg=_estimate_cost(precursor_smi),
            ))
        for reagent in template.reagents:
            rp = _reagent_to_precursor(reagent)
            if rp is not None:
                precursors.append(rp)
        return precursors

    # PROTECTION / DEPROTECTION: the core structure is essentially kept.
    if cat in (ReactionCategory.PROTECTION, ReactionCategory.DEPROTECTION):
        # Precursor is the unprotected / protected form.
        precursor_smi = _toggle_protection(
            target_smiles, target_mol, fg, template)
        if precursor_smi:
            precursors.append(Precursor(
                smiles=precursor_smi, molecule=None,
                name=f"precursor ({template.name})",
                cost_per_kg=_estimate_cost(precursor_smi),
            ))
        return precursors

    # REARRANGEMENT: return the pre-rearrangement skeleton
    if cat == ReactionCategory.REARRANGEMENT:
        precursor_smi = _reverse_rearrangement(
            target_smiles, target_mol, fg, template)
        if precursor_smi:
            precursors.append(Precursor(
                smiles=precursor_smi, molecule=None,
                name=f"precursor ({template.name})",
                cost_per_kg=_estimate_cost(precursor_smi),
            ))
        return precursors

    # Fallback: try a simple truncation
    precursor_smi = _simplify_molecule(target_smiles, target_mol, fg)
    if precursor_smi:
        precursors.append(Precursor(
            smiles=precursor_smi, molecule=None,
            name=f"simplified precursor",
            cost_per_kg=_estimate_cost(precursor_smi),
        ))
    return precursors


# =====================================================================
#  Molecular transform helpers (heuristic / simplified)
# =====================================================================

def _modify_fg_smiles(
    target_smiles: str,
    target_mol: Molecule,
    fg: FunctionalGroup,
    template: ReactionTemplate,
    direction: str,
) -> str | None:
    """Heuristically modify a functional group for redox transforms.

    For 'oxidise': alcohol -> aldehyde/ketone, aldehyde -> carboxylic acid.
    For 'reduce' : aldehyde/ketone -> alcohol, carboxylic acid -> aldehyde.

    Returns a precursor SMILES string or None on failure.
    """
    fg_name = fg.name
    try:
        if direction == "oxidise":
            # Product was reduced, so precursor is oxidised form
            if fg_name == "alcohol":
                # Precursor is the corresponding aldehyde or ketone
                return _replace_oh_with_carbonyl(target_smiles)
            if fg_name in ("aldehyde", "ketone"):
                # Precursor might be a carboxylic acid
                return target_smiles  # keep same (template applies to it)
        else:  # reduce
            if fg_name == "aldehyde":
                return _replace_carbonyl_with_oh(target_smiles)
            if fg_name == "ketone":
                return _replace_carbonyl_with_oh(target_smiles)
            if fg_name == "carboxylic_acid":
                return target_smiles
    except Exception:
        pass
    return None


def _replace_oh_with_carbonyl(smiles: str) -> str:
    """Very rough: replace first 'O' in the SMILES with '=O' pattern."""
    # This is a heuristic simplification for the retrosynthetic planner.
    # We try to find an -OH and convert it conceptually to C=O.
    if "CO" in smiles and "C=O" not in smiles:
        return smiles.replace("CO", "C=O", 1)
    return smiles


def _replace_carbonyl_with_oh(smiles: str) -> str:
    """Rough: replace C=O with C-OH (alcohol form)."""
    if "C=O" in smiles:
        return smiles.replace("C=O", "CO", 1)
    return smiles


def _split_at_fg(
    target_smiles: str,
    target_mol: Molecule,
    fg: FunctionalGroup,
    template: ReactionTemplate,
) -> list[str]:
    """Split the target into two fragment SMILES at the functional group.

    Used for coupling / carbonyl reactions where two components combine.
    The heuristic removes the functional group centre and tries to return
    the two largest remaining fragments as SMILES.
    """
    center = fg.center
    fg_atoms_set = set(fg.atoms)

    # Find bonds connecting FG atoms to the rest of the molecule
    break_bonds: list[tuple[int, int]] = []
    for a_idx in fg.atoms:
        for nb in target_mol.neighbors(a_idx):
            if nb not in fg_atoms_set:
                break_bonds.append((a_idx, nb))

    if len(break_bonds) < 2:
        # Cannot split meaningfully -- return the whole thing simplified
        simp = _simplify_molecule(target_smiles, target_mol, fg)
        return [simp] if simp else [target_smiles]

    # Build two fragment atom sets by BFS from each side of the break
    fragments: list[set[int]] = []
    all_atoms = set(range(len(target_mol.atoms)))
    excluded = fg_atoms_set

    visited_global: set[int] = set()
    for _, outside_atom in break_bonds:
        if outside_atom in visited_global:
            continue
        # BFS from outside_atom, not crossing into fg_atoms
        frag: set[int] = set()
        stack = [outside_atom]
        while stack:
            cur = stack.pop()
            if cur in frag or cur in excluded:
                continue
            frag.add(cur)
            for nb in target_mol.neighbors(cur):
                if nb not in frag and nb not in excluded:
                    stack.append(nb)
        if frag:
            visited_global |= frag
            fragments.append(frag)

    # Convert each fragment to SMILES using a simplified approach:
    # We generate a sub-SMILES by collecting the heavy-atom symbols
    # and connecting them linearly.  This is an approximation.
    result_smiles: list[str] = []
    for frag in fragments[:2]:
        smi = _fragment_to_smiles(target_mol, frag)
        if smi:
            result_smiles.append(smi)

    # If we only got one fragment, add a simple reagent as the second
    if len(result_smiles) == 1:
        for reagent in template.reagents:
            rp = _reagent_to_precursor(reagent)
            if rp is not None:
                result_smiles.append(rp.smiles)
                break
        else:
            result_smiles.append("C")  # methane fallback

    if not result_smiles:
        result_smiles = [target_smiles]

    return result_smiles


def _fragment_to_smiles(mol: Molecule, atom_indices: set[int]) -> str:
    """Build an approximate SMILES for a subset of atoms in *mol*.

    Constructs a new Molecule from the selected atoms (excluding H),
    copies the bonds between them, and serialises with to_smiles.
    """
    heavy_indices = sorted(
        idx for idx in atom_indices if mol.atoms[idx].symbol != "H"
    )
    if not heavy_indices:
        return ""

    # Build a sub-molecule
    sub = Molecule(name="fragment")
    old_to_new: dict[int, int] = {}
    for old_idx in heavy_indices:
        atom = mol.atoms[old_idx]
        new_idx = sub.add_atom(atom.symbol, atom.position.copy(),
                               atom.hybridization)
        old_to_new[old_idx] = new_idx

    # Copy bonds within the fragment
    for bond in mol.bonds:
        if bond.atom_i in old_to_new and bond.atom_j in old_to_new:
            ni = old_to_new[bond.atom_i]
            nj = old_to_new[bond.atom_j]
            # Avoid duplicate bonds
            if sub.get_bond(ni, nj) is None:
                sub.add_bond(ni, nj, order=bond.order, rotatable=bond.rotatable)

    # Add implicit hydrogens to satisfy valence (approximate)
    # We rely on the SMILES writer to handle implicit H.
    try:
        return to_smiles(sub)
    except Exception:
        # Fallback: concatenate symbols
        return "".join(mol.atoms[i].symbol for i in heavy_indices[:6])


def _substitute_fg(
    target_smiles: str,
    target_mol: Molecule,
    fg: FunctionalGroup,
    template: ReactionTemplate,
) -> str | None:
    """For substitution reactions, swap the produced FG for the required one.

    E.g. if the template produces an alcohol from an alkyl halide, the
    precursor is the alkyl halide form.
    """
    # Determine what FG the precursor should have
    required = template.functional_group_required
    produced = template.functional_group_produced

    # Map common FG swaps in SMILES
    swap_map = {
        ("alcohol", "alkyl_halide"): ("O", "Br"),
        ("ether", "alkyl_halide"): ("OC", "Br"),
        ("ether", "alcohol"): ("OC", "O"),
        ("primary_amine", "alkyl_halide"): ("N", "Br"),
        ("nitrile", "alkyl_halide"): ("C#N", "Br"),
        ("azide", "alkyl_halide"): ("N=[N+]=[N-]", "Br"),
    }

    fg_name = fg.name
    for req in required:
        key = (fg_name, req)
        if key in swap_map:
            old_frag, new_frag = swap_map[key]
            if old_frag in target_smiles:
                return target_smiles.replace(old_frag, new_frag, 1)

    # Generic fallback: just return the target with a halide substitution
    if fg_name == "alcohol" and "O" in target_smiles:
        return target_smiles.replace("O", "Br", 1)
    return None


def _add_across_double_bond(
    target_smiles: str,
    target_mol: Molecule,
    fg: FunctionalGroup,
    template: ReactionTemplate,
) -> str | None:
    """Reverse of elimination: add HX across a double bond to get precursor."""
    # If the target has an alkene, the precursor is an alkyl halide/alcohol.
    if fg.name == "alkene" and "C=C" in target_smiles:
        # Add H and Br across the double bond
        return target_smiles.replace("C=C", "CC(Br)", 1)
    return None


def _remove_addition(
    target_smiles: str,
    target_mol: Molecule,
    fg: FunctionalGroup,
    template: ReactionTemplate,
) -> str | None:
    """Reverse of addition: remove the added group to restore alkene."""
    fg_name = fg.name

    # The template required an alkene and produced the current FG
    if "alkene" in template.functional_group_required:
        # Restore the alkene by removing the added functionality
        if fg_name == "alcohol" and "CO" in target_smiles:
            return target_smiles.replace("CO", "C=C", 1)
        if fg_name.startswith("alkyl_halide"):
            for hal in ("Br", "Cl", "I"):
                if f"C{hal}" in target_smiles:
                    return target_smiles.replace(f"C{hal}", "C=C", 1)
        if fg_name == "epoxide" and "C1OC1" in target_smiles:
            return target_smiles.replace("C1OC1", "C=C", 1)
    return None


def _toggle_protection(
    target_smiles: str,
    target_mol: Molecule,
    fg: FunctionalGroup,
    template: ReactionTemplate,
) -> str | None:
    """Toggle between protected and deprotected forms.

    Simplification: for protection templates, just return the target
    since the core structure is essentially preserved.
    """
    return target_smiles


def _reverse_rearrangement(
    target_smiles: str,
    target_mol: Molecule,
    fg: FunctionalGroup,
    template: ReactionTemplate,
) -> str | None:
    """Rough approximation for reversing a rearrangement.

    Returns the target itself as a stand-in since rearrangement
    precursors are structural isomers that are hard to derive
    without full subgraph matching.
    """
    return target_smiles


def _simplify_molecule(
    target_smiles: str,
    target_mol: Molecule,
    fg: FunctionalGroup,
) -> str | None:
    """Produce a simplified precursor by removing part of the molecule.

    Heuristic: remove the functional group atoms and return the largest
    connected fragment.
    """
    fg_set = set(fg.atoms)
    remaining = set(range(len(target_mol.atoms))) - fg_set
    heavy = {i for i in remaining if target_mol.atoms[i].symbol != "H"}

    if not heavy:
        return None

    # Find largest connected component among remaining heavy atoms
    visited: set[int] = set()
    best_comp: set[int] = set()
    for start in heavy:
        if start in visited:
            continue
        comp: set[int] = set()
        stack = [start]
        while stack:
            cur = stack.pop()
            if cur in comp or cur in fg_set:
                continue
            if cur not in heavy:
                continue
            comp.add(cur)
            for nb in target_mol.neighbors(cur):
                if nb not in comp and nb in heavy:
                    stack.append(nb)
        visited |= comp
        if len(comp) > len(best_comp):
            best_comp = comp

    if not best_comp:
        return None

    return _fragment_to_smiles(target_mol, best_comp)


def _reagent_to_precursor(reagent_str: str) -> Precursor | None:
    """Try to match a reagent string to a purchasable material."""
    # Map common reagent names to SMILES
    reagent_map: dict[str, str] = {
        "NaOH": "O",
        "NaCN": "C#N",
        "NaN3": "CN=[N+]=[N-]",
        "NaOMe": "CO",
        "NaOEt": "CCO",
        "HBr": "Br",
        "HCl": "Cl",
        "H2O": "O",
        "MeOH": "CO",
        "EtOH": "CCO",
        "NaBH4": "O",
        "LiAlH4": "O",
        "H2": "O",
        "BH3*THF": "C1CCOC1",
        "mCPBA": "O",
        "Br2": "Br",
        "PCC": "O",
        "n-BuLi": "CCCC",
    }
    smi = reagent_map.get(reagent_str)
    if smi is not None:
        p = get_purchasable(smi)
        if p is not None:
            return p
    return None


def _estimate_cost(smiles: str) -> float:
    """Estimate cost per kg for a SMILES string.

    Uses purchasable DB if available, otherwise estimates based on
    molecular size.
    """
    entry = PURCHASABLE_MATERIALS.get(smiles)
    if entry is not None:
        return entry[1]
    # Rough estimate: $10/kg per heavy atom
    try:
        mol = parse(smiles)
        n_heavy = _count_heavy_atoms(mol)
        return max(5.0, n_heavy * 10.0)
    except Exception:
        return 50.0


# =====================================================================
#  Template matching: which templates apply to a given FG?
# =====================================================================

def _find_applicable_templates(
    fg: FunctionalGroup,
    all_fg_names: list[str],
) -> list[ReactionTemplate]:
    """Find reaction templates whose *produced* FG matches *fg*.

    In retrosynthesis we work backwards: we look for reactions that
    **produce** the functional group found on the target, because
    reversing such a reaction gives us the precursors.

    Also includes templates that **require** the FG (the forward
    reaction transforms it, so running it in reverse generates a
    molecule with that FG as starting material from something simpler).
    """
    results: list[ReactionTemplate] = []
    seen_names: set[str] = set()

    # Templates that produce this FG (primary retrosynthetic match)
    for tmpl in find_reactions_producing(fg.name):
        if tmpl.name not in seen_names and tmpl.is_compatible(all_fg_names):
            results.append(tmpl)
            seen_names.add(tmpl.name)

    # Also consider templates that require this FG (the forward
    # reaction uses this FG as a handle).
    for tmpl in lookup_by_functional_group(fg.name):
        if tmpl.name not in seen_names and tmpl.is_compatible(all_fg_names):
            results.append(tmpl)
            seen_names.add(tmpl.name)

    # Handle generic alkyl_halide name for specific halides
    if fg.name.startswith("alkyl_halide_"):
        for tmpl in lookup_by_functional_group("alkyl_halide"):
            if tmpl.name not in seen_names and tmpl.is_compatible(all_fg_names):
                results.append(tmpl)
                seen_names.add(tmpl.name)

    return results


# =====================================================================
#  Beam search retrosynthesis
# =====================================================================

def _build_retro_node(
    smiles: str,
    mol: Molecule,
    depth: int,
    max_depth: int,
    beam_width: int,
    visited_smiles: set[str],
) -> RetroNode:
    """Build one node of the retrosynthesis tree.

    If the molecule is purchasable, the node is a leaf.  Otherwise,
    functional groups are detected, templates are matched, disconnections
    are scored, and the top *beam_width* disconnections are kept.  The
    best disconnection's precursors are then expanded recursively.
    """
    node = RetroNode(
        smiles=smiles,
        molecule=mol,
        depth=depth,
    )

    # Check purchasability
    if is_purchasable(smiles):
        node.is_purchasable = True
        return node

    # Detect functional groups
    fgs = detect_functional_groups(mol)
    node.functional_groups = fgs

    # Depth limit
    if depth >= max_depth:
        return node

    # Prevent infinite loops
    if smiles in visited_smiles:
        return node
    visited_smiles = visited_smiles | {smiles}

    # Collect all FG names for compatibility check
    all_fg_names = [fg.name for fg in fgs]

    # Generate disconnections
    disconnections: list[Disconnection] = []
    seen_templates: set[str] = set()

    for fg in fgs:
        templates = _find_applicable_templates(fg, all_fg_names)
        for tmpl in templates:
            if tmpl.name in seen_templates:
                continue
            seen_templates.add(tmpl.name)

            precursors = _generate_precursors_for_template(
                smiles, mol, tmpl, fg)
            if not precursors:
                continue

            score = score_disconnection(tmpl, precursors, mol)
            disconnections.append(Disconnection(
                template=tmpl, precursors=precursors, score=score))

    # Sort by score (highest first) and keep top beam_width
    disconnections.sort(key=lambda d: d.score, reverse=True)
    node.disconnections = disconnections[:beam_width]

    # Select best disconnection
    if node.disconnections:
        node.best_disconnection = node.disconnections[0]

        # Recursively expand precursors of the best disconnection
        for precursor in node.best_disconnection.precursors:
            if is_purchasable(precursor.smiles):
                child_mol = _safe_parse(precursor.smiles)
                if child_mol is None:
                    continue
                child_node = RetroNode(
                    smiles=precursor.smiles,
                    molecule=child_mol,
                    is_purchasable=True,
                    depth=depth + 1,
                )
                node.children.append(child_node)
            else:
                child_mol = _safe_parse(precursor.smiles)
                if child_mol is None:
                    continue
                child_node = _build_retro_node(
                    precursor.smiles, child_mol,
                    depth + 1, max_depth, beam_width,
                    visited_smiles,
                )
                node.children.append(child_node)

    return node


def _safe_parse(smiles: str) -> Molecule | None:
    """Parse SMILES, returning None on failure."""
    try:
        return parse(smiles)
    except Exception:
        return None


def _count_routes(node: RetroNode) -> int:
    """Count complete routes (paths from root to all-purchasable leaves)."""
    if node.is_purchasable:
        return 1
    if not node.children:
        return 0
    # A route is complete when all children are resolved
    child_counts = [_count_routes(c) for c in node.children]
    if all(c > 0 for c in child_counts):
        # Multiply: each combination of child routes is a complete route
        product = 1
        for c in child_counts:
            product *= c
        return product
    return 0


def retrosynthesis(
    mol: Molecule,
    max_depth: int = 8,
    beam_width: int = 5,
) -> RetrosynthesisTree:
    """Perform retrosynthetic analysis on a target molecule.

    Starting from the target, the algorithm works backwards:
    1. Convert target to SMILES and check purchasability.
    2. Detect functional groups on the target.
    3. For each FG, look up matching reaction templates from the
       knowledge base.
    4. For each matching template, generate precursor molecules by
       conceptual reverse transform.
    5. Score each disconnection using strategic bond preference,
       atom-count reduction, FG simplification, template yield, and
       precursor availability.
    6. Keep the top *beam_width* disconnections.
    7. Recurse on non-purchasable precursors up to *max_depth*.
    8. Mark the best route through the tree.

    Parameters
    ----------
    mol : Molecule
        The target molecule to analyse.
    max_depth : int
        Maximum number of retrosynthetic steps to explore.
    beam_width : int
        Number of disconnections to keep at each level.

    Returns
    -------
    RetrosynthesisTree
        The full retrosynthesis tree with scored disconnections.
    """
    target_smiles = to_smiles(mol)
    visited: set[str] = set()

    root = _build_retro_node(
        target_smiles, mol, 0, max_depth, beam_width, visited)

    routes = _count_routes(root)

    return RetrosynthesisTree(
        target=root,
        max_depth=max_depth,
        beam_width=beam_width,
        routes_found=routes,
    )


# =====================================================================
#  Tree formatting (ASCII text)
# =====================================================================

def _format_node(node: RetroNode, indent: str, is_last: bool,
                 lines: list[str]) -> None:
    """Recursively format a node and its children as an ASCII tree."""
    connector = "`-- " if is_last else "|-- "
    status = ""
    if node.is_purchasable:
        entry = PURCHASABLE_MATERIALS.get(node.smiles)
        name = entry[0] if entry else "purchasable"
        status = f" [AVAILABLE: {name}]"
    elif node.best_disconnection:
        tmpl_name = node.best_disconnection.template.name
        score = node.best_disconnection.score
        status = f" <-- {tmpl_name} (score={score:.1f})"
    else:
        if node.depth > 0:
            status = " [no route found]"

    lines.append(f"{indent}{connector}{node.smiles}{status}")

    # Continuation indent for children
    child_indent = indent + ("    " if is_last else "|   ")

    # Show functional groups at the root
    if node.depth == 0 and node.functional_groups:
        fg_names = ", ".join(fg.name for fg in node.functional_groups)
        lines.append(f"{child_indent}FGs: {fg_names}")

    # Show alternative disconnections (briefly)
    if node.disconnections and len(node.disconnections) > 1:
        lines.append(f"{child_indent}({len(node.disconnections)} "
                     f"disconnection(s) evaluated)")

    # Recurse into children
    for i, child in enumerate(node.children):
        is_last_child = (i == len(node.children) - 1)
        _format_node(child, child_indent, is_last_child, lines)


def format_tree(tree: RetrosynthesisTree) -> str:
    """Format a RetrosynthesisTree as an ASCII text diagram.

    Parameters
    ----------
    tree : RetrosynthesisTree
        The retrosynthesis tree to format.

    Returns
    -------
    str
        Multi-line ASCII text representation of the tree.

    Example output::

        Retrosynthetic Analysis
        ==================================================
        Target: CC(=O)O
        Max depth: 8    Beam width: 5    Routes found: 2
        ==================================================
        `-- CC(=O)O <-- Fischer esterification (score=72.3)
            FGs: carboxylic_acid, alcohol
            (3 disconnection(s) evaluated)
            |-- CC(O)=O [AVAILABLE: acetic acid]
            `-- CCO [AVAILABLE: ethanol]
    """
    lines: list[str] = []
    lines.append("Retrosynthetic Analysis")
    lines.append("=" * 58)
    lines.append(f"Target: {tree.target.smiles}")
    lines.append(
        f"Max depth: {tree.max_depth}    "
        f"Beam width: {tree.beam_width}    "
        f"Routes found: {tree.routes_found}"
    )
    lines.append("=" * 58)

    _format_node(tree.target, "", True, lines)

    return "\n".join(lines)
