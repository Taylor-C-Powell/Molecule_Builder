"""Safety assessment for synthesis routes.

Performs GHS hazard lookup for every reagent in every step, determines
PPE requirements, engineering controls, emergency procedures, and
produces a per-step :class:`SafetyAssessment`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from molbuilder.reactions.reaction_types import ReactionTemplate
from molbuilder.reactions.reagent_data import get_reagent, normalize_reagent_name


# =====================================================================
#  GHS reference tables
# =====================================================================

GHS_PICTOGRAMS: dict[str, str] = {
    "GHS01": "Exploding bomb -- Explosives, self-reactive, organic peroxides",
    "GHS02": "Flame -- Flammable gases/liquids/solids, pyrophoric, self-heating, "
             "emits flammable gas on water contact",
    "GHS03": "Flame over circle -- Oxidizers",
    "GHS04": "Gas cylinder -- Compressed, liquefied, or dissolved gases",
    "GHS05": "Corrosion -- Corrosive to metals, skin corrosion, serious eye damage",
    "GHS06": "Skull and crossbones -- Acute toxicity (fatal or toxic)",
    "GHS07": "Exclamation mark -- Irritant, narcotic, acute toxicity (harmful), "
             "skin sensitizer",
    "GHS08": "Health hazard -- Carcinogenicity, mutagenicity, reproductive toxicity, "
             "respiratory sensitizer, organ toxicity, aspiration hazard",
    "GHS09": "Environment -- Aquatic toxicity",
}

GHS_HAZARD_STATEMENTS: dict[str, str] = {
    # Physical hazards
    "H200": "Unstable explosive",
    "H201": "Explosive; mass explosion hazard",
    "H202": "Explosive; severe projection hazard",
    "H203": "Explosive; fire, blast or projection hazard",
    "H204": "Fire or projection hazard",
    "H205": "May mass explode in fire",
    "H220": "Extremely flammable gas",
    "H221": "Flammable gas",
    "H222": "Extremely flammable aerosol",
    "H223": "Flammable aerosol",
    "H224": "Extremely flammable liquid and vapour",
    "H225": "Highly flammable liquid and vapour",
    "H226": "Flammable liquid and vapour",
    "H227": "Combustible liquid",
    "H228": "Flammable solid",
    "H240": "Heating may cause an explosion",
    "H241": "Heating may cause a fire or explosion",
    "H242": "Heating may cause a fire",
    "H250": "Catches fire spontaneously if exposed to air",
    "H251": "Self-heating; may catch fire",
    "H252": "Self-heating in large quantities; may catch fire",
    "H260": "In contact with water releases flammable gases which may ignite spontaneously",
    "H261": "In contact with water releases flammable gas",
    "H270": "May cause or intensify fire; oxidizer",
    "H271": "May cause fire or explosion; strong oxidizer",
    "H272": "May intensify fire; oxidizer",
    "H280": "Contains gas under pressure; may explode if heated",
    "H281": "Contains refrigerated gas; may cause cryogenic burns or injury",
    "H290": "May be corrosive to metals",
    # Health hazards
    "H300": "Fatal if swallowed",
    "H301": "Toxic if swallowed",
    "H302": "Harmful if swallowed",
    "H304": "May be fatal if swallowed and enters airways",
    "H305": "May be harmful if swallowed and enters airways",
    "H310": "Fatal in contact with skin",
    "H311": "Toxic in contact with skin",
    "H312": "Harmful in contact with skin",
    "H314": "Causes severe skin burns and eye damage",
    "H315": "Causes skin irritation",
    "H317": "May cause an allergic skin reaction",
    "H318": "Causes serious eye damage",
    "H319": "Causes serious eye irritation",
    "H330": "Fatal if inhaled",
    "H331": "Toxic if inhaled",
    "H332": "Harmful if inhaled",
    "H334": "May cause allergy or asthma symptoms or breathing difficulties if inhaled",
    "H335": "May cause respiratory irritation",
    "H336": "May cause drowsiness or dizziness",
    "H340": "May cause genetic defects",
    "H341": "Suspected of causing genetic defects",
    "H350": "May cause cancer",
    "H351": "Suspected of causing cancer",
    "H360": "May damage fertility or the unborn child",
    "H361": "Suspected of damaging fertility or the unborn child",
    "H362": "May cause harm to breast-fed children",
    "H360F": "May damage fertility",
    "H360D": "May damage the unborn child",
    "H360Fd": "May damage fertility; suspected of damaging the unborn child",
    "H360Df": "May damage the unborn child; suspected of damaging fertility",
    "H370": "Causes damage to organs",
    "H371": "May cause damage to organs",
    "H372": "Causes damage to organs through prolonged or repeated exposure",
    "H373": "May cause damage to organs through prolonged or repeated exposure",
    # Environmental hazards
    "H400": "Very toxic to aquatic life",
    "H410": "Very toxic to aquatic life with long lasting effects",
    "H411": "Toxic to aquatic life with long lasting effects",
    "H412": "Harmful to aquatic life with long lasting effects",
    "H413": "May cause long lasting harmful effects to aquatic life",
    "H420": "Harms public health and the environment by destroying ozone in the upper atmosphere",
}


# =====================================================================
#  Data classes
# =====================================================================

@dataclass
class HazardInfo:
    """GHS hazard information for a single reagent."""

    reagent_name: str
    ghs_hazards: list[str]
    ghs_pictograms: list[str]
    hazard_descriptions: list[str]
    pictogram_descriptions: list[str]


@dataclass
class IncompatibilityWarning:
    """A chemical incompatibility with severity and mitigation guidance.

    Backwards-compatible: ``str(warning)`` returns the same format as the
    old plain-string incompatibility messages, so code that iterates
    ``incompatible_materials`` and calls ``str()`` keeps working.
    """
    reagent_a: str
    reagent_b: str
    hazard: str
    severity: str       # "critical", "high", "medium"
    mitigation: str

    def __str__(self) -> str:
        prefix = self.severity.upper() + ": " if self.severity == "critical" else ""
        return f"{prefix}{self.hazard}"

    def __repr__(self) -> str:
        return (f"IncompatibilityWarning({self.severity}: "
                f"{self.reagent_a} + {self.reagent_b})")


@dataclass
class ThermalHazard:
    """Thermal hazard assessment for a reaction step."""

    reaction_type: str
    severity: str          # "low", "medium", "high", "critical"
    description: str
    max_temp_c: float
    mitigation: str


@dataclass
class SafetyAssessment:
    """Complete safety assessment for one synthesis step."""

    step_number: int
    step_name: str
    hazards: list[HazardInfo]
    ppe_required: list[str]
    engineering_controls: list[str]
    emergency_procedures: list[str]
    incompatible_materials: list[IncompatibilityWarning]
    thermal_hazards: list[ThermalHazard]
    waste_classification: str
    risk_level: str             # "low", "medium", "high"

    @property
    def overall_severity(self) -> str | None:
        """Highest severity across all incompatibility warnings, or None."""
        order = {"critical": 3, "high": 2, "medium": 1}
        best = 0
        best_label = None
        for w in self.incompatible_materials:
            rank = order.get(w.severity, 0)
            if rank > best:
                best = rank
                best_label = w.severity
        return best_label


# =====================================================================
#  Internal helpers
# =====================================================================

def _build_hazard_info(reagent_name: str) -> HazardInfo:
    """Look up a reagent and compile its hazard information."""
    reagent = get_reagent(reagent_name)
    if reagent is None:
        return HazardInfo(
            reagent_name=reagent_name,
            ghs_hazards=[],
            ghs_pictograms=[],
            hazard_descriptions=["Hazard data not available in database"],
            pictogram_descriptions=[],
        )
    haz_descs = [
        GHS_HAZARD_STATEMENTS.get(h, f"Unknown hazard code {h}")
        for h in reagent.ghs_hazards
    ]
    pic_descs = [
        GHS_PICTOGRAMS.get(p, f"Unknown pictogram {p}")
        for p in reagent.ghs_pictograms
    ]
    return HazardInfo(
        reagent_name=reagent.name,
        ghs_hazards=list(reagent.ghs_hazards),
        ghs_pictograms=list(reagent.ghs_pictograms),
        hazard_descriptions=haz_descs,
        pictogram_descriptions=pic_descs,
    )


def _determine_ppe(all_hazards: set[str], all_pictograms: set[str]) -> list[str]:
    """Determine required PPE from the union of hazard codes."""
    ppe: list[str] = []
    # Always require baseline
    ppe.append("Safety goggles (splash-proof)")
    ppe.append("Lab coat")
    ppe.append("Closed-toe shoes")

    if "GHS05" in all_pictograms or "H314" in all_hazards:
        ppe.append("Chemical-resistant gloves (e.g. butyl rubber or nitrile)")
        ppe.append("Face shield")
        ppe.append("Chemical-resistant apron")
    else:
        ppe.append("Nitrile gloves (double-gloving recommended)")

    if any(h in all_hazards for h in ("H330", "H331", "H332", "H335")):
        ppe.append("Respiratory protection (fume hood minimum; respirator if outside hood)")

    if "GHS06" in all_pictograms or "H300" in all_hazards or "H310" in all_hazards:
        ppe.append("Emergency eyewash and safety shower must be accessible within 10 seconds")

    if "GHS02" in all_pictograms or "H250" in all_hazards or "H260" in all_hazards:
        ppe.append("Fire-resistant lab coat or Nomex coveralls for pyrophoric work")

    return ppe


def _determine_engineering_controls(
    all_hazards: set[str],
    all_pictograms: set[str],
    template: ReactionTemplate,
) -> list[str]:
    """Determine engineering controls."""
    controls: list[str] = []

    # Fume hood is baseline for organic chemistry
    controls.append("Perform all operations in a well-ventilated fume hood")

    if "GHS02" in all_pictograms:
        controls.append("Remove all ignition sources; use non-sparking tools")
        controls.append("Ground and bond all containers to prevent static discharge")

    if "H250" in all_hazards or "H260" in all_hazards:
        controls.append("Schlenk line or glovebox required for air/moisture-sensitive reagents")

    if "GHS03" in all_pictograms:
        controls.append("Keep oxidizers separated from fuels and organic materials")
        controls.append("Fire suppression system accessible")

    if any(h in all_hazards for h in ("H340", "H350", "H360")):
        controls.append("Designated area for CMR (carcinogenic/mutagenic/reprotoxic) substances")
        controls.append("HEPA-filtered ventilation for solid handling")

    if "GHS09" in all_pictograms:
        controls.append("Secondary containment to prevent environmental release")
        controls.append("Do not dispose of via drain; collect all waste")

    mean_t = (template.temperature_range[0] + template.temperature_range[1]) / 2.0
    if mean_t < -40:
        controls.append("Cryogenic cooling equipment; ensure adequate ventilation for cryogen vapours")
    if mean_t > 150:
        controls.append("High-temperature operation; thermal insulation and burn protection required")

    return controls


def _determine_emergency_procedures(
    all_hazards: set[str],
    all_pictograms: set[str],
) -> list[str]:
    """Emergency response procedures based on hazard profile."""
    procedures: list[str] = []

    procedures.append(
        "In case of spill: evacuate area, ventilate, absorb with "
        "inert material (vermiculite), dispose as hazardous waste."
    )

    if "GHS06" in all_pictograms or "H300" in all_hazards:
        procedures.append(
            "Ingestion: Do NOT induce vomiting.  Call Poison Control / "
            "emergency services immediately.  Rinse mouth with water."
        )

    if "H310" in all_hazards or "H311" in all_hazards:
        procedures.append(
            "Skin contact: Remove contaminated clothing immediately.  "
            "Wash skin with copious water for at least 15 minutes.  "
            "Seek medical attention."
        )

    if "H330" in all_hazards:
        procedures.append(
            "Inhalation: Move victim to fresh air.  If not breathing, "
            "administer artificial respiration.  Call emergency services."
        )

    if "H314" in all_hazards or "H318" in all_hazards:
        procedures.append(
            "Eye contact: Flush eyes with water for at least 15 minutes, "
            "lifting upper and lower eyelids.  Seek ophthalmological evaluation."
        )

    if "GHS02" in all_pictograms or "H250" in all_hazards:
        procedures.append(
            "Fire: Use dry chemical, CO2, or sand extinguisher.  "
            "Do NOT use water on pyrophoric / water-reactive materials.  "
            "Evacuate if fire cannot be controlled immediately."
        )

    if "H260" in all_hazards:
        procedures.append(
            "Water-reactive material: In case of spill, cover with dry "
            "sand or vermiculite.  Do NOT use water."
        )

    return procedures


def _determine_incompatibilities(
    template: ReactionTemplate,
    solvent: str | None = None,
) -> list[IncompatibilityWarning]:
    """List known incompatible material pairs in the step.

    When *solvent* is provided, it is included in the reagent key set
    so that all existing pair checks naturally catch solvent-reagent
    conflicts.  Additional solvent-specific checks are also applied.
    """
    warnings: list[IncompatibilityWarning] = []

    reagent_keys = {normalize_reagent_name(r) for r in template.reagents}
    if solvent:
        reagent_keys.add(normalize_reagent_name(solvent))

    def _warn(a: str, b: str, hazard: str, severity: str, mitigation: str):
        warnings.append(IncompatibilityWarning(
            reagent_a=a, reagent_b=b, hazard=hazard,
            severity=severity, mitigation=mitigation,
        ))

    # -- Original 14 pairs (unchanged logic) ---------------------------------

    # Oxidizer + reducer
    oxidizers = {"kmno4", "cro3", "h2o2", "naocl", "mcpba", "hno3", "naio4"}
    reducers = {"nabh4", "lialh4", "dibal_h", "na_nh3", "red_al", "nah"}
    if reagent_keys & oxidizers and reagent_keys & reducers:
        _warn("oxidizer", "reducer",
              "Oxidizer and reducer present -- do NOT mix directly",
              "critical", "Add reducer to substrate first, then introduce oxidizer slowly")

    # Acid + cyanide
    acids = {"hcl", "h2so4", "hno3", "tfa", "acoh", "bf3_oet2"}
    cyanides = {"nacn"}
    if reagent_keys & acids and reagent_keys & cyanides:
        _warn("acid", "cyanide",
              "Acid + cyanide can liberate HCN gas (fatal)",
              "critical", "Use cyanide only in basic or neutral pH; HCN detector required")

    # Acid + azide
    azides = {"nan3"}
    if reagent_keys & acids and reagent_keys & azides:
        _warn("acid", "azide",
              "Acid + azide can liberate HN3 (explosive/toxic)",
              "critical", "Maintain basic pH; blast shield required")

    # Water-reactive + aqueous
    water_reactive = {"lialh4", "nah", "n_buli", "memgbr", "etmgbr", "phmgbr",
                      "meli", "phli", "socl2", "accl", "ticl4", "pcl5"}
    if reagent_keys & water_reactive:
        _warn("water-reactive", "moisture",
              "Water-reactive reagent(s) present: ensure strictly anhydrous conditions",
              "high", "Use oven-dried glassware, Schlenk technique, molecular sieves")

    # Peroxides + metals
    peroxides = {"h2o2", "mcpba", "tbhp", "dtbp", "benzoyl_peroxide"}
    metals = {"ticl4", "alcl3", "zncl2", "cui", "fecl3", "fecl2"}
    if reagent_keys & peroxides and reagent_keys & metals:
        _warn("peroxide", "metal salt",
              "Peroxide + metal salt: risk of uncontrolled decomposition",
              "high", "Add peroxide slowly with temperature monitoring")

    # Hypochlorite + acids -> chlorine gas
    hypochlorites = {"naocl", "bleach", "calcium_hypochlorite"}
    if reagent_keys & hypochlorites and reagent_keys & acids:
        _warn("hypochlorite", "acid",
              "Hypochlorite + acid liberates Cl2 gas (toxic)",
              "critical", "Never acidify hypochlorite solutions; Cl2 detector required")

    # Permanganate + concentrated organics -> fire/explosion
    permanganates = {"kmno4", "namno4"}
    flammable_organics = {"acetone", "diethyl_ether", "thf", "ethanol", "methanol",
                          "toluene", "hexane", "pentane", "dcm", "dmf", "dmso"}
    if reagent_keys & permanganates and reagent_keys & flammable_organics:
        _warn("permanganate", "flammable organic",
              "Permanganate + flammable organic: fire/explosion risk",
              "high", "Use dilute permanganate in aqueous medium; keep organic solvent volume minimal")

    # Alkali metals + water -> violent reaction
    alkali_metals = {"na", "k", "li", "cs", "nah", "kh"}
    aqueous = {"h2o", "water", "naoh_aq", "hcl_aq"}
    if reagent_keys & alkali_metals and reagent_keys & aqueous:
        _warn("alkali metal/hydride", "water",
              "Alkali metal/hydride + water -> violent H2 evolution",
              "critical", "Use anhydrous solvents; quench slowly with IPA/ice if needed")

    # Chlorine/halogens + ammonia -> toxic gases
    halogens = {"cl2", "br2", "i2", "ncs"}
    ammonia = {"nh3", "nh4oh", "nh4cl"}
    if reagent_keys & halogens and reagent_keys & ammonia:
        _warn("halogen", "ammonia",
              "Halogen + ammonia -> toxic NCl3/NBr3",
              "critical", "Do not mix; use separate vessels and transfer lines")

    # Strong oxidizers + flammable solvents
    strong_oxidizers = {"kmno4", "cro3", "k2cr2o7", "hno3_conc", "h2o2_conc",
                        "naio4", "oxone", "dmp"}
    if reagent_keys & strong_oxidizers and reagent_keys & flammable_organics:
        _warn("strong oxidizer", "flammable solvent",
              "Strong oxidizer + flammable solvent: fire risk -- use compatible solvent",
              "high", "Switch to water or acetic acid as solvent where possible")

    # Nitrates + organics -> explosion risk
    nitrates = {"nano3", "kno3", "agno3", "nh4no3"}
    if reagent_keys & nitrates and reagent_keys & flammable_organics:
        _warn("nitrate salt", "organic solvent",
              "Nitrate salt + organic solvent: explosion risk at elevated temperature",
              "high", "Keep temperature below 100 C; avoid confinement")

    # Concentrated acids + concentrated bases -> exothermic
    bases = {"naoh", "koh", "lioh", "naoh_conc", "koh_conc", "nah", "naoh_aq"}
    conc_acids = {"h2so4", "hno3", "hcl_conc", "h3po4", "hf"}
    if reagent_keys & conc_acids and reagent_keys & bases:
        _warn("concentrated acid", "base",
              "Concentrated acid + base: highly exothermic neutralization -- add slowly with cooling",
              "high", "Add acid to base (never reverse); use ice bath and slow addition")

    # Grignard/organolithium + protic solvents
    organometallics = {"memgbr", "etmgbr", "phmgbr", "meli", "phli",
                       "n_buli", "t_buli", "s_buli", "znet2", "znme2"}
    protic_solvents = {"meoh", "etoh", "h2o", "acoh", "ipoh"}
    if reagent_keys & organometallics and reagent_keys & protic_solvents:
        _warn("organometallic", "protic solvent",
              "Organometallic + protic solvent: immediate quench -- use ethereal solvents only",
              "critical", "Use THF or diethyl ether; ensure strictly anhydrous conditions")

    # Perchloric acid + organics -> explosive perchlorates
    perchlorics = {"hclo4", "perchloric_acid"}
    if reagent_keys & perchlorics and reagent_keys & flammable_organics:
        _warn("perchloric acid", "organic",
              "Perchloric acid + organics -> explosive perchlorate esters",
              "critical", "Never heat perchloric acid with organics; use perchloric acid hood")

    # -- New pairs (CAMEO / Bretherick's) ------------------------------------

    # Azide + heavy metals -> detonation-sensitive azides
    heavy_metals = {"pb_oac2", "hgcl2", "hg_oac2", "pbcl2", "agno3"}
    if reagent_keys & azides and reagent_keys & heavy_metals:
        _warn("azide", "heavy metal",
              "Azide + heavy metals (Pb, Hg, Ag) -> detonation-sensitive metal azides",
              "critical", "Avoid combining; use copper or zinc catalysts instead if needed")

    # Nitric acid + alcohols -> alkyl nitrate explosives
    alcohols = {"meoh", "etoh", "ipoh", "n_buoh", "glycerol"}
    nitric_acids = {"hno3", "hno3_conc"}
    if reagent_keys & nitric_acids and reagent_keys & alcohols:
        _warn("nitric acid", "alcohol",
              "Nitric acid + alcohols -> alkyl nitrate explosives",
              "critical", "Use dilute HNO3; keep temperature below 0 C; blast shield required")

    # Chloroform + strong base -> phosgene
    chloroform = {"chcl3", "chloroform", "cdcl3"}
    strong_bases = {"naoh", "koh", "naoh_conc", "koh_conc", "n_buli", "t_buli",
                    "sodium_hydroxide", "potassium_hydroxide"}
    if reagent_keys & chloroform and reagent_keys & strong_bases:
        _warn("chloroform", "strong base",
              "Chloroform + strong base -> phosgene (COCl2, fatal)",
              "critical", "Use DCM instead; if unavoidable, work in fume hood with phosgene detector")

    # DMSO + acyl halides -> exothermic decomposition
    dmso = {"dmso"}
    acyl_halides = {"accl", "socl2", "pcl5", "pcl3", "oxalyl_chloride"}
    if reagent_keys & dmso and reagent_keys & acyl_halides:
        _warn("DMSO", "acyl halide",
              "DMSO + acyl halides -> exothermic decomposition",
              "high", "Use alternative solvent (DMF or NMP); if unavoidable, add at -78 C")

    # Hydrogen peroxide + acetone -> TATP precursor
    peroxide_tatp = {"h2o2", "h2o2_conc"}
    acetone_set = {"acetone"}
    if reagent_keys & peroxide_tatp and reagent_keys & acetone_set:
        _warn("hydrogen peroxide", "acetone",
              "Hydrogen peroxide + acetone -> TATP explosive precursor (in presence of acid)",
              "critical", "Never combine H2O2 with acetone; use alternative solvent for peroxide reactions")

    # Picric acid + metals -> shock-sensitive metal picrates
    picric = {"picric_acid", "2_4_6_trinitrophenol"}
    metal_containers = {"na", "k", "li", "zn", "pb_oac2", "cu", "fecl3"}
    if reagent_keys & picric and reagent_keys & metal_containers:
        _warn("picric acid", "metal",
              "Picric acid + metals -> metal picrates (shock-sensitive explosives)",
              "critical", "Store picric acid in plastic; avoid contact with metals")

    # Ether + strong oxidizers -> peroxide explosion
    ethers = {"diethyl_ether", "thf", "dioxane", "dme", "tbme", "cpme"}
    if reagent_keys & ethers and reagent_keys & strong_oxidizers:
        _warn("ether", "strong oxidizer",
              "Ether solvents + strong oxidizers -> peroxide formation and explosion risk",
              "high", "Test ether for peroxides before use; inhibit with BHT if stored")

    # Carbon disulfide + alkali metals -> violent reaction
    cs2 = {"cs2", "carbon_disulfide"}
    if reagent_keys & cs2 and reagent_keys & alkali_metals:
        _warn("carbon disulfide", "alkali metal",
              "Carbon disulfide + alkali metals -> violent exothermic reaction",
              "critical", "Do not combine; use alternative sulfur sources")

    # Ammonia + mercury compounds -> fulminating mercury
    mercury = {"hgcl2", "hg_oac2", "hg2cl2"}
    if reagent_keys & ammonia and reagent_keys & mercury:
        _warn("ammonia", "mercury compound",
              "Ammonia + mercury compounds -> fulminating mercury (shock-sensitive)",
              "critical", "Never combine; substitute mercury catalyst with safer alternative")

    # Acetone + chloroform + base -> chlorobutanol (unexpected side reaction)
    if reagent_keys & acetone_set and reagent_keys & chloroform and reagent_keys & bases:
        _warn("acetone + chloroform", "base",
              "Acetone + chloroform + base -> chlorobutanol side reaction",
              "medium", "Avoid mixing; use separate vessels for extraction steps")

    # -- Solvent-specific checks -----------------------------------------------

    if solvent:
        solvent_key = normalize_reagent_name(solvent)

        # NaH / water-reactive in protic or aqueous solvent
        protic_solvent_keys = {"meoh", "etoh", "ipoh", "h2o", "water",
                               "methanol", "ethanol", "isopropanol", "acoh"}
        water_reactive_all = {
            "lialh4", "lithium_aluminium_hydride",
            "nah", "sodium_hydride",
            "n_buli", "n_butyllithium", "t_buli", "s_buli",
            "memgbr", "etmgbr", "phmgbr", "meli", "phli",
        }
        if solvent_key in protic_solvent_keys and reagent_keys & water_reactive_all:
            _warn("protic solvent", "water-reactive reagent",
                  f"Water-reactive reagent in protic solvent ({solvent}): violent reaction risk",
                  "critical",
                  "Switch to ethereal solvent (THF, Et2O); ensure anhydrous conditions")

        # DMSO thermal decomposition at high temperature
        if solvent_key == "dmso":
            mean_t = (template.temperature_range[0] + template.temperature_range[1]) / 2.0
            if mean_t > 150:
                _warn("DMSO", "high temperature",
                      "DMSO decomposes above 189 C with toxic by-products (Me2S, formaldehyde)",
                      "high",
                      "Use alternative high-boiling solvent (NMP, DMI); "
                      "limit temperature to below 180 C")

        # Ether peroxide accumulation warning when used as primary solvent
        ether_solvents = {"diethyl_ether", "thf", "dioxane", "dme", "tbme", "cpme"}
        if solvent_key in ether_solvents and reagent_keys & strong_oxidizers:
            # Already caught above in ether + strong oxidizer check, but
            # reinforce the solvent-specific context
            pass  # covered by existing check

    return warnings


def _assess_thermal_hazards(template: ReactionTemplate) -> list[ThermalHazard]:
    """Identify thermal hazards based on template category and reagents."""
    hazards: list[ThermalHazard] = []
    reagent_keys = {normalize_reagent_name(r) for r in template.reagents}
    catalyst_keys = {normalize_reagent_name(c) for c in template.catalysts}
    all_keys = reagent_keys | catalyst_keys
    name_lower = (template.named_reaction or "").lower()

    # Grignard formation: organometallic + Mg
    grignard_reagents = {"memgbr", "etmgbr", "phmgbr"}
    if all_keys & grignard_reagents or "grignard" in name_lower:
        hazards.append(ThermalHazard(
            reaction_type="Grignard formation",
            severity="critical",
            description="Highly exothermic initiation; risk of runaway if reaction stalls then restarts",
            max_temp_c=35.0,
            mitigation="Use slow addition with ice bath; monitor temperature continuously; "
                       "have dry-ice/acetone bath ready for emergency cooling",
        ))

    # Catalytic hydrogenation
    h2_keys = {"h2_pd_c", "pd_c_10", "lindlar"}
    if all_keys & h2_keys or "hydrogenation" in name_lower:
        hazards.append(ThermalHazard(
            reaction_type="Catalytic hydrogenation",
            severity="high",
            description="Exothermic hydrogen uptake; catalyst can ignite H2/air mixtures",
            max_temp_c=50.0,
            mitigation="Purge vessel with inert gas before introducing H2; "
                       "control H2 pressure; keep catalyst wet when filtering",
        ))

    # Diazotization
    if "diazotization" in name_lower or "diazo" in name_lower:
        hazards.append(ThermalHazard(
            reaction_type="Diazotization",
            severity="critical",
            description="Diazonium salts are thermally unstable and may decompose explosively above 5 C",
            max_temp_c=5.0,
            mitigation="Maintain temperature at 0-5 C with ice-salt bath; "
                       "prepare fresh and use immediately; never isolate diazonium salts",
        ))

    # Friedel-Crafts
    fc_catalysts = {"alcl3", "bf3_oet2", "ticl4"}
    if all_keys & fc_catalysts and (
        "friedel" in name_lower or "acylation" in name_lower or "alkylation" in name_lower
    ):
        hazards.append(ThermalHazard(
            reaction_type="Friedel-Crafts reaction",
            severity="high",
            description="Exothermic complexation with Lewis acid; HCl gas evolution",
            max_temp_c=25.0,
            mitigation="Add Lewis acid catalyst slowly; maintain ice-bath cooling; "
                       "vent HCl through scrubber",
        ))

    # Nitration
    if "nitration" in name_lower or ("hno3" in all_keys and "h2so4" in all_keys):
        hazards.append(ThermalHazard(
            reaction_type="Nitration",
            severity="critical",
            description="Mixed acid nitration is highly exothermic with risk of runaway and explosion",
            max_temp_c=10.0,
            mitigation="Add nitrating agent slowly below 10 C; use blast shield; "
                       "have quench water available",
        ))

    # Sulfonation
    if "sulfonation" in name_lower:
        hazards.append(ThermalHazard(
            reaction_type="Sulfonation",
            severity="high",
            description="Exothermic reaction with fuming sulfuric acid or SO3",
            max_temp_c=40.0,
            mitigation="Add sulfonating agent slowly with cooling; "
                       "use temperature-controlled addition funnel",
        ))

    # Ozonolysis
    if "ozonolysis" in name_lower or "o3" in all_keys:
        hazards.append(ThermalHazard(
            reaction_type="Ozonolysis",
            severity="high",
            description="Ozonides are explosive intermediates; ozone is a strong oxidant",
            max_temp_c=-78.0,
            mitigation="Perform at -78 C in DCM; quench ozonides immediately with Me2S or PPh3; "
                       "never concentrate ozonide solutions",
        ))

    # Metalation (n-BuLi, t-BuLi, s-BuLi)
    lithium_bases = {"n_buli", "n_butyllithium", "t_buli", "s_buli",
                     "meli", "phli"}
    if all_keys & lithium_bases:
        hazards.append(ThermalHazard(
            reaction_type="Organolithium metalation",
            severity="critical",
            description="Pyrophoric reagent; exothermic metalation; "
                       "t-BuLi ignites on air contact",
            max_temp_c=-78.0,
            mitigation="Handle under inert atmosphere (Schlenk or glovebox); "
                       "cool to -78 C before addition; use syringe pump for controlled addition",
        ))

    # Wolff-Kishner (hydrazine + base + heat)
    if "wolff" in name_lower or ("n2h4" in all_keys or "nh2nh2" in all_keys):
        hazards.append(ThermalHazard(
            reaction_type="Wolff-Kishner reduction",
            severity="high",
            description="Exothermic decomposition of hydrazones; N2 gas evolution at reflux",
            max_temp_c=200.0,
            mitigation="Use well-ventilated reflux condenser; add hydrazine slowly; "
                       "ensure outlet for N2 gas to prevent pressure buildup",
        ))

    return hazards


def _classify_waste(all_hazards: set[str]) -> str:
    """Classify waste stream based on worst hazard codes."""
    if any(h in all_hazards for h in ("H300", "H310", "H330", "H340", "H350", "H360")):
        return "Hazardous waste -- Category 1 (acute/CMR): requires licensed disposal contractor"
    if any(h in all_hazards for h in ("H301", "H311", "H314", "H331", "H400", "H410")):
        return "Hazardous waste -- Category 2 (toxic/corrosive/ecotoxic): segregated collection"
    if any(h in all_hazards for h in ("H225", "H224", "H220", "H228")):
        return "Hazardous waste -- flammable: store in approved flammable-waste containers"
    return "Non-hazardous chemical waste: collect in appropriate waste stream"


def _calculate_risk_level(all_hazards: set[str], all_pictograms: set[str]) -> str:
    """Assign overall risk level for the step."""
    # High risk: acutely fatal, CMR, pyrophoric, explosive
    high_codes = {"H200", "H201", "H240", "H250", "H300", "H310", "H330",
                  "H340", "H350", "H360"}
    if all_hazards & high_codes:
        return "high"

    # Medium risk: toxic, corrosive, flammable, sensitizer
    medium_codes = {"H301", "H311", "H314", "H331", "H225", "H224",
                    "H260", "H334", "H317", "H370", "H372"}
    if all_hazards & medium_codes:
        return "medium"

    return "low"


# =====================================================================
#  Public API
# =====================================================================

def assess_safety(steps: list[Any]) -> list[SafetyAssessment]:
    """Produce a :class:`SafetyAssessment` for every step in *steps*.

    Parameters
    ----------
    steps : list
        Each element must have a ``.template`` attribute
        (:class:`ReactionTemplate`) and a ``.precursors`` attribute.
        Duck typing is used.

    Returns
    -------
    list[SafetyAssessment]
        One assessment per step, in order.
    """
    if not steps:
        return []
    for i, step in enumerate(steps):
        if not hasattr(step, 'template'):
            raise TypeError(
                f"Step {i} must have a 'template' attribute, "
                f"got {type(step).__name__}"
            )

    assessments: list[SafetyAssessment] = []

    for idx, step in enumerate(steps):
        template: ReactionTemplate = step.template

        # Collect hazard info for all reagents + catalysts
        all_reagent_names = list(template.reagents) + list(template.catalysts)
        hazard_infos: list[HazardInfo] = [
            _build_hazard_info(rname) for rname in all_reagent_names
        ]

        # Aggregate hazard codes and pictograms across all reagents in this step
        all_hazards: set[str] = set()
        all_pictograms: set[str] = set()
        for hi in hazard_infos:
            all_hazards.update(hi.ghs_hazards)
            all_pictograms.update(hi.ghs_pictograms)

        # Determine solvent from conditions if available
        solvent = None
        if hasattr(step, "conditions") and step.conditions is not None:
            solvent = getattr(step.conditions, "solvent", None)

        assessments.append(SafetyAssessment(
            step_number=idx + 1,
            step_name=template.name,
            hazards=hazard_infos,
            ppe_required=_determine_ppe(all_hazards, all_pictograms),
            engineering_controls=_determine_engineering_controls(
                all_hazards, all_pictograms, template,
            ),
            emergency_procedures=_determine_emergency_procedures(
                all_hazards, all_pictograms,
            ),
            incompatible_materials=_determine_incompatibilities(
                template, solvent=solvent,
            ),
            thermal_hazards=_assess_thermal_hazards(template),
            waste_classification=_classify_waste(all_hazards),
            risk_level=_calculate_risk_level(all_hazards, all_pictograms),
        ))

    return assessments
