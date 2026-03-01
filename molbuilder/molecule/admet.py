"""ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) prediction.

Heuristic ADMET profiling using established QSAR rules and structural
alerts.  All predictions are rule-based (no ML models) and use only
properties already computed by MolBuilder: Lipinski descriptors,
functional group detection, and atom-walking on the molecular graph.

References
----------
- Lipinski, C.A. et al. (2001). Experimental and computational
  approaches to estimate solubility and permeability. Adv. Drug Deliv.
  Rev. 46, 3-26.
- Veber, D.F. et al. (2002). Molecular properties that influence the
  oral bioavailability of drug candidates. J. Med. Chem. 45, 2615-2623.
- Ertl, P. (2000). Estimation of synthetic accessibility score of
  drug-like molecules. J. Cheminf. 1, 8.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from molbuilder.molecule.graph import Molecule
from molbuilder.molecule.properties import lipinski_properties
from molbuilder.reactions.functional_group_detect import detect_functional_groups


# =====================================================================
#  Result dataclass
# =====================================================================

@dataclass
class ADMETProfile:
    """Full ADMET profile for a molecule."""
    # Absorption
    oral_bioavailability: str
    intestinal_absorption: str
    caco2_permeability: str
    pgp_substrate: bool
    # Distribution
    bbb_penetrant: bool
    plasma_protein_binding: str
    vd_class: str
    # Metabolism
    cyp_inhibition: dict[str, bool]
    metabolic_stability: str
    # Excretion
    renal_clearance: str
    half_life_class: str
    # Toxicity
    herg_risk: str
    ames_mutagenicity: bool
    hepatotoxicity_risk: str
    structural_alerts: list[str]
    # Summary
    overall_score: float
    warnings: list[str]
    flags: list[str]


# =====================================================================
#  Private helpers -- Structural analysis
# =====================================================================

def _has_basic_nitrogen(mol: Molecule) -> bool:
    """Detect basic (protonatable) nitrogen via atom-walking.

    A nitrogen is considered basic if it is sp3 (amine-like) or is in
    a non-aromatic context with at least one H or alkyl neighbor.
    """
    for atom in mol.atoms:
        if atom.symbol != "N":
            continue
        # Skip nitro, nitrile, amide patterns
        neighbors = mol.neighbors(atom.index)
        bond_orders = []
        for n_idx in neighbors:
            bond = mol.get_bond(atom.index, n_idx)
            if bond:
                bond_orders.append(bond.order)
        # Aromatic N (bond order 4) is weakly basic at best
        if any(o == 4 for o in bond_orders):
            continue
        # If N has a triple bond (nitrile) or double bond to O (nitro), skip
        if any(o == 3 for o in bond_orders):
            continue
        double_to_o = False
        for n_idx in neighbors:
            bond = mol.get_bond(atom.index, n_idx)
            if bond and bond.order == 2 and mol.atoms[n_idx].symbol == "O":
                double_to_o = True
                break
        if double_to_o:
            continue
        # Has H neighbor or is sp3 -> basic
        has_h = any(mol.atoms[n].symbol == "H" for n in neighbors)
        if has_h or len(neighbors) <= 3:
            return True
    return False


def _aromatic_ring_count(mol: Molecule) -> int:
    """Count aromatic rings from aromatic bond pairs."""
    aromatic_atoms: set[int] = set()
    for bond in mol.bonds:
        if bond.order == 4:
            aromatic_atoms.add(bond.atom_i)
            aromatic_atoms.add(bond.atom_j)
    heavy_aromatic = {i for i in aromatic_atoms if mol.atoms[i].symbol != "H"}
    if len(heavy_aromatic) < 5:
        return 0
    return max(len(heavy_aromatic) // 5, 1)


def _fg_names(mol: Molecule) -> set[str]:
    """Return set of detected functional group names."""
    try:
        fgs = detect_functional_groups(mol)
        return {fg.name for fg in fgs}
    except Exception:
        return set()


def _tpsa(props: object) -> float:
    """Extract TPSA from Lipinski properties."""
    return getattr(props, "tpsa", 0.0) or 0.0


# =====================================================================
#  Private helpers -- ADMET predictions
# =====================================================================

def _oral_bioavailability(props: object) -> str:
    """Predict oral bioavailability from Lipinski Ro5 + Veber rules."""
    violations = 0
    if props.molecular_weight > 500:
        violations += 1
    if props.logp > 5:
        violations += 1
    if props.hbd > 5:
        violations += 1
    if props.hba > 10:
        violations += 1
    # Veber criteria
    tpsa = _tpsa(props)
    if props.rotatable_bonds > 10:
        violations += 1
    if tpsa < 20 or tpsa > 130:
        violations += 1
    if violations == 0:
        return "high"
    if violations <= 2:
        return "moderate"
    return "low"


def _intestinal_absorption(props: object) -> str:
    """Predict intestinal absorption from TPSA."""
    tpsa = _tpsa(props)
    if tpsa < 140:
        return "high"
    if tpsa < 200:
        return "moderate"
    return "low"


def _caco2_permeability(props: object) -> str:
    """Predict Caco-2 permeability from logP and TPSA."""
    tpsa = _tpsa(props)
    if 1 <= props.logp <= 3 and tpsa < 80:
        return "high"
    if props.logp < 0 or tpsa > 120:
        return "low"
    return "moderate"


def _pgp_substrate(props: object) -> bool:
    """Predict P-glycoprotein substrate likelihood."""
    return (
        props.molecular_weight > 400
        and (props.hbd > 2 or props.hba > 8)
        and props.logp > 1
    )


def _bbb_penetration(props: object) -> bool:
    """Predict blood-brain barrier penetration."""
    tpsa = _tpsa(props)
    return (
        props.molecular_weight < 400
        and 1 <= props.logp <= 5
        and tpsa < 60
        and props.hbd <= 3
    )


def _plasma_protein_binding(logp: float) -> str:
    """Predict plasma protein binding from logP."""
    if logp > 4:
        return "high"
    if logp > 2:
        return "moderate"
    return "low"


def _volume_of_distribution(logp: float) -> str:
    """Predict volume of distribution class from logP."""
    if logp > 3:
        return "high"
    if logp > 1:
        return "moderate"
    return "low"


def _cyp_inhibition(mol: Molecule, fgs: set[str]) -> dict[str, bool]:
    """Predict CYP inhibition for 5 major isoforms.

    Uses structural features: azole rings (imidazole, triazole),
    aromatic amines, large lipophilic scaffolds.
    """
    result = {
        "CYP1A2": False,
        "CYP2C9": False,
        "CYP2C19": False,
        "CYP2D6": False,
        "CYP3A4": False,
    }

    # Detect azole-like N-containing heterocycles (strong CYP inhibitors)
    has_azole = "imidazole" in fgs or "triazole" in fgs or "pyrazole" in fgs
    has_aromatic_amine = "aniline" in fgs or "aromatic_amine" in fgs

    # Count aromatic rings for CYP3A4 (large substrate cavity)
    ar_count = _aromatic_ring_count(mol)

    # CYP1A2: planar polyaromatic compounds
    if ar_count >= 2:
        result["CYP1A2"] = True

    # CYP2C9: acidic + aromatic (e.g. NSAIDs)
    if ("carboxylic_acid" in fgs or "sulfonic_acid" in fgs) and ar_count >= 1:
        result["CYP2C9"] = True

    # CYP2C19: azoles are strong inhibitors
    if has_azole:
        result["CYP2C19"] = True

    # CYP2D6: basic nitrogen + aromatic
    if _has_basic_nitrogen(mol) and ar_count >= 1:
        result["CYP2D6"] = True

    # CYP3A4: large lipophilic molecules, azoles
    if has_azole or ar_count >= 3:
        result["CYP3A4"] = True

    # Aromatic amines can inhibit CYP1A2
    if has_aromatic_amine:
        result["CYP1A2"] = True

    return result


def _metabolic_stability(props: object) -> str:
    """Predict metabolic stability from logP, MW, rotatable bonds."""
    unstable_score = 0
    if props.logp > 4:
        unstable_score += 1
    if props.molecular_weight > 500:
        unstable_score += 1
    if props.rotatable_bonds > 10:
        unstable_score += 1
    if unstable_score >= 2:
        return "low"
    if unstable_score == 1:
        return "moderate"
    return "high"


def _renal_clearance(props: object) -> str:
    """Predict renal clearance from MW and logP."""
    if props.molecular_weight < 400 and props.logp < 2:
        return "high"
    if props.molecular_weight > 500 or props.logp > 3:
        return "low"
    return "moderate"


def _half_life(stability: str, ppb: str, vd: str) -> str:
    """Estimate half-life class from stability, PPB, and Vd.

    Components that increase half-life: high PPB, high Vd, low stability.
    """
    score = 0
    if stability == "low":
        score -= 1  # cleared faster
    if stability == "high":
        score += 1
    if ppb == "high":
        score += 1  # sequestered in plasma
    if vd == "high":
        score += 1  # distributed widely
    if score >= 2:
        return "long"
    if score <= -1:
        return "short"
    return "moderate"


def _herg_risk(props: object, mol: Molecule) -> str:
    """Predict hERG channel inhibition risk.

    Key risk factors: logP > 3.7, basic nitrogen, aromatic rings, MW > 350.
    """
    score = 0
    if props.logp > 3.7:
        score += 1
    if _has_basic_nitrogen(mol):
        score += 1
    if _aromatic_ring_count(mol) >= 2:
        score += 1
    if props.molecular_weight > 350:
        score += 1
    if score >= 3:
        return "high"
    if score >= 2:
        return "moderate"
    return "low"


def _ames_mutagenicity(mol: Molecule, fgs: set[str]) -> bool:
    """Predict Ames mutagenicity from structural alerts.

    Positive for nitro aromatics, aromatic amines, epoxides, aziridines.
    """
    mutagenic_fgs = {"nitro", "epoxide", "aziridine"}
    if fgs & mutagenic_fgs:
        return True
    # Aromatic amine (aniline-type)
    if "aniline" in fgs or "aromatic_amine" in fgs:
        return True
    # Nitro group on aromatic ring: check for nitro + aromatic system
    has_nitro = any(
        a.symbol == "N" for a in mol.atoms
        if any(
            mol.get_bond(a.index, n) and mol.get_bond(a.index, n).order == 2
            and mol.atoms[n].symbol == "O"
            for n in mol.neighbors(a.index)
        )
    )
    if has_nitro and _aromatic_ring_count(mol) > 0:
        return True
    return False


def _hepatotoxicity(mol: Molecule, fgs: set[str]) -> str:
    """Predict hepatotoxicity risk from reactive metabolite alerts."""
    alerts = 0
    reactive_fgs = {
        "hydrazine", "hydrazide", "epoxide", "acyl_halide",
        "anhydride", "isocyanate", "isothiocyanate",
    }
    alerts += len(fgs & reactive_fgs)
    # Quinone-like: check for para-diketone on aromatic ring
    if "quinone" in fgs:
        alerts += 1
    if alerts >= 2:
        return "high"
    if alerts == 1:
        return "moderate"
    return "low"


def _structural_alerts(mol: Molecule, fgs: set[str]) -> list[str]:
    """Detect PAINS-like structural alerts."""
    alerts = []

    # Rhodanines
    if "thiazolidinedione" in fgs or "rhodanine" in fgs:
        alerts.append("rhodanine / thiazolidinedione")

    # Catechol (ortho-dihydroxybenzene) -- check for phenol pairs
    phenol_count = sum(1 for fg in fgs if fg in ("phenol",))
    if phenol_count >= 2:
        alerts.append("catechol (ortho-diphenol)")

    # Michael acceptors (alpha,beta-unsaturated carbonyls)
    if "enone" in fgs or "michael_acceptor" in fgs:
        alerts.append("Michael acceptor")

    # Acyl hydrazides
    if "hydrazide" in fgs:
        alerts.append("acyl hydrazide")

    # Epoxide
    if "epoxide" in fgs:
        alerts.append("epoxide (alkylating agent)")

    # Nitro groups
    if "nitro" in fgs:
        alerts.append("nitro group (redox cycling)")

    # Aldehyde (reactive carbonyl)
    if "aldehyde" in fgs:
        alerts.append("aldehyde (reactive carbonyl)")

    return alerts


# =====================================================================
#  Overall score
# =====================================================================

def _overall_score(profile: ADMETProfile) -> float:
    """Compute weighted ADMET score (0-10, higher = better profile).

    Weights reflect clinical importance: toxicity > absorption > metabolism.
    """
    score = 10.0

    # Absorption (weight ~2.5)
    if profile.oral_bioavailability == "low":
        score -= 1.5
    elif profile.oral_bioavailability == "moderate":
        score -= 0.5
    if profile.intestinal_absorption == "low":
        score -= 0.5
    if profile.caco2_permeability == "low":
        score -= 0.5

    # Distribution (weight ~1.0)
    if profile.plasma_protein_binding == "high":
        score -= 0.5
    # BBB is context-dependent, not penalized

    # Metabolism (weight ~1.5)
    if profile.metabolic_stability == "low":
        score -= 1.0
    elif profile.metabolic_stability == "moderate":
        score -= 0.3
    cyp_count = sum(1 for v in profile.cyp_inhibition.values() if v)
    score -= cyp_count * 0.1

    # Excretion (weight ~0.5)
    if profile.half_life_class == "short":
        score -= 0.5

    # Toxicity (weight ~3.0)
    if profile.herg_risk == "high":
        score -= 1.5
    elif profile.herg_risk == "moderate":
        score -= 0.5
    if profile.ames_mutagenicity:
        score -= 1.0
    if profile.hepatotoxicity_risk == "high":
        score -= 1.0
    elif profile.hepatotoxicity_risk == "moderate":
        score -= 0.3
    score -= len(profile.structural_alerts) * 0.2

    return round(max(0.0, min(10.0, score)), 1)


# =====================================================================
#  Warnings and flags
# =====================================================================

def _generate_warnings(profile: ADMETProfile) -> list[str]:
    """Generate human-readable warnings for concerning ADMET properties."""
    warnings = []
    if profile.oral_bioavailability == "low":
        warnings.append("Poor predicted oral bioavailability")
    if profile.herg_risk == "high":
        warnings.append("High hERG inhibition risk -- cardiac safety concern")
    if profile.ames_mutagenicity:
        warnings.append("Positive Ames mutagenicity prediction")
    if profile.hepatotoxicity_risk == "high":
        warnings.append("High hepatotoxicity risk from reactive metabolite alerts")
    if profile.metabolic_stability == "low":
        warnings.append("Low metabolic stability -- rapid clearance expected")
    if profile.pgp_substrate:
        warnings.append("Likely P-gp substrate -- may limit oral absorption")
    return warnings


def _generate_flags(profile: ADMETProfile) -> list[str]:
    """Generate flags for notable (not necessarily bad) properties."""
    flags = []
    if profile.bbb_penetrant:
        flags.append("Predicted BBB penetrant -- consider CNS side effects")
    if profile.plasma_protein_binding == "high":
        flags.append("High plasma protein binding -- may affect drug interactions")
    if profile.structural_alerts:
        flags.append(
            f"{len(profile.structural_alerts)} structural alert(s) detected"
        )
    cyp_count = sum(1 for v in profile.cyp_inhibition.values() if v)
    if cyp_count >= 3:
        flags.append(f"Predicted inhibitor of {cyp_count}/5 CYP isoforms")
    return flags


# =====================================================================
#  Public API
# =====================================================================

def predict_admet(mol: Molecule) -> ADMETProfile:
    """Predict ADMET properties for a molecule.

    Returns an ADMETProfile with absorption, distribution, metabolism,
    excretion, and toxicity predictions plus an overall score (0-10).
    """
    props = lipinski_properties(mol)
    fgs = _fg_names(mol)

    # Absorption
    oral_bio = _oral_bioavailability(props)
    ia = _intestinal_absorption(props)
    caco2 = _caco2_permeability(props)
    pgp = _pgp_substrate(props)

    # Distribution
    bbb = _bbb_penetration(props)
    ppb = _plasma_protein_binding(props.logp)
    vd = _volume_of_distribution(props.logp)

    # Metabolism
    cyp = _cyp_inhibition(mol, fgs)
    stability = _metabolic_stability(props)

    # Excretion
    renal = _renal_clearance(props)
    hl = _half_life(stability, ppb, vd)

    # Toxicity
    herg = _herg_risk(props, mol)
    ames = _ames_mutagenicity(mol, fgs)
    hepato = _hepatotoxicity(mol, fgs)
    alerts = _structural_alerts(mol, fgs)

    profile = ADMETProfile(
        oral_bioavailability=oral_bio,
        intestinal_absorption=ia,
        caco2_permeability=caco2,
        pgp_substrate=pgp,
        bbb_penetrant=bbb,
        plasma_protein_binding=ppb,
        vd_class=vd,
        cyp_inhibition=cyp,
        metabolic_stability=stability,
        renal_clearance=renal,
        half_life_class=hl,
        herg_risk=herg,
        ames_mutagenicity=ames,
        hepatotoxicity_risk=hepato,
        structural_alerts=alerts,
        overall_score=0.0,  # placeholder
        warnings=[],
        flags=[],
    )

    profile.overall_score = _overall_score(profile)
    profile.warnings = _generate_warnings(profile)
    profile.flags = _generate_flags(profile)

    return profile
