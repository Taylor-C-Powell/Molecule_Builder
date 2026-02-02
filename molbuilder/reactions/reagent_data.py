"""Reagent and solvent databases with GHS hazard information.

This module provides two master dictionaries -- ``REAGENT_DB`` and
``SOLVENT_DB`` -- that catalogue the most common reagents and solvents
encountered in organic synthesis.  Each entry carries physical-property
data, approximate cost, GHS hazard classifications, and storage guidance.

Lookup helpers (``get_reagent``, ``get_solvent``, ``search_reagents``,
``search_solvents``) allow fuzzy retrieval by name or formula fragment.
"""

from __future__ import annotations
from dataclasses import dataclass, field


# =====================================================================
#  Reagent name normalisation
# =====================================================================

# Alias dictionary for common alternative names
_REAGENT_ALIASES: dict[str, str] = {
    "lialh4": "lithium_aluminium_hydride",
    "nabh4": "sodium_borohydride",
    "pcc": "pyridinium_chlorochromate",
    "dcc": "dicyclohexylcarbodiimide",
    "mcpba": "meta_chloroperoxybenzoic_acid",
    "tbscl": "tert_butyldimethylsilyl_chloride",
    "tmscl": "trimethylsilyl_chloride",
    "boc2o": "di_tert_butyl_dicarbonate",
    "fmoc_cl": "fluorenylmethyloxycarbonyl_chloride",
    "pd_c": "palladium_on_carbon",
    "pd/c": "palladium_on_carbon",
    "h2so4": "sulfuric_acid",
    "hcl": "hydrochloric_acid",
    "naoh": "sodium_hydroxide",
    "koh": "potassium_hydroxide",
    "kmno4": "potassium_permanganate",
    "k2cr2o7": "potassium_dichromate",
    "lioh": "lithium_hydroxide",
    "nah": "sodium_hydride",
    "n_buli": "n_butyllithium",
    "n-buli": "n_butyllithium",
    "t_buok": "potassium_tert_butoxide",
    "t-buok": "potassium_tert_butoxide",
    "et3n": "triethylamine",
    "dbu": "1_8_diazabicyclo_5_4_0_undec_7_ene",
    "dmap": "4_dimethylaminopyridine",
    "tfa": "trifluoroacetic_acid",
}


def normalize_reagent_name(name: str) -> str:
    """Normalize a reagent name to its canonical key.

    Lowercases, replaces whitespace/hyphens with underscores, and checks
    the alias dictionary for common abbreviations.
    """
    key = name.strip().lower().replace(" ", "_").replace("-", "_")
    return _REAGENT_ALIASES.get(key, key)


# =====================================================================
#  Data classes
# =====================================================================

@dataclass
class Reagent:
    name: str
    formula: str
    cas: str = ""
    mw: float = 0.0
    cost_per_kg: float = 0.0       # USD approximate
    ghs_hazards: list[str] = field(default_factory=list)
    ghs_pictograms: list[str] = field(default_factory=list)
    storage: str = ""
    physical_state: str = ""        # solid, liquid, gas


@dataclass
class Solvent:
    name: str
    formula: str
    bp: float = 0.0                # boiling point degC
    mp: float = 0.0                # melting point degC
    density: float = 1.0           # g/mL
    dielectric: float = 1.0
    polarity_index: float = 0.0
    cas: str = ""
    cost_per_L: float = 0.0
    green_score: int = 5           # 1=very green, 10=avoid
    ghs_hazards: list[str] = field(default_factory=list)
    miscible_with_water: bool = False


# =====================================================================
#  REAGENT_DB  (~100 entries)
# =====================================================================

REAGENT_DB: dict[str, Reagent] = {

    # -----------------------------------------------------------------
    #  Oxidizing agents
    # -----------------------------------------------------------------
    "kmno4": Reagent(
        name="Potassium permanganate", formula="KMnO4", cas="7722-64-7",
        mw=158.03, cost_per_kg=35,
        ghs_hazards=["H272", "H302", "H314", "H410"],
        ghs_pictograms=["GHS03", "GHS05", "GHS07", "GHS09"],
        storage="Cool dry place, away from organics", physical_state="solid",
    ),
    "cro3": Reagent(
        name="Chromium trioxide", formula="CrO3", cas="1333-82-0",
        mw=99.99, cost_per_kg=120,
        ghs_hazards=["H271", "H301", "H311", "H314", "H317", "H330", "H340", "H350", "H361", "H372", "H410"],
        ghs_pictograms=["GHS03", "GHS05", "GHS06", "GHS08", "GHS09"],
        storage="Locked oxidizer cabinet, separate from organics", physical_state="solid",
    ),
    "pcc": Reagent(
        name="Pyridinium chlorochromate", formula="C5H5NHCrO3Cl", cas="26299-14-9",
        mw=215.56, cost_per_kg=280,
        ghs_hazards=["H301", "H312", "H314", "H317", "H340", "H350"],
        ghs_pictograms=["GHS05", "GHS06", "GHS08"],
        storage="Desiccator, away from organics", physical_state="solid",
    ),
    "pdc": Reagent(
        name="Pyridinium dichromate", formula="(C5H5NH)2Cr2O7", cas="20039-37-6",
        mw=376.21, cost_per_kg=350,
        ghs_hazards=["H301", "H312", "H314", "H317", "H340", "H350"],
        ghs_pictograms=["GHS05", "GHS06", "GHS08"],
        storage="Desiccator, away from organics", physical_state="solid",
    ),
    "oxalyl_chloride": Reagent(
        name="Oxalyl chloride", formula="(COCl)2", cas="79-37-8",
        mw=126.93, cost_per_kg=75,
        ghs_hazards=["H301", "H314", "H330"],
        ghs_pictograms=["GHS05", "GHS06"],
        storage="Fridge under inert gas, moisture-sensitive", physical_state="liquid",
    ),
    "dmso": Reagent(
        name="Dimethyl sulfoxide", formula="(CH3)2SO", cas="67-68-5",
        mw=78.13, cost_per_kg=15,
        ghs_hazards=[],
        ghs_pictograms=[],
        storage="Room temperature, tightly sealed", physical_state="liquid",
    ),
    "jones_reagent": Reagent(
        name="Jones reagent (CrO3/H2SO4)", formula="CrO3/H2SO4",
        mw=0.0, cost_per_kg=0,
        ghs_hazards=["H271", "H301", "H314", "H340", "H350"],
        ghs_pictograms=["GHS03", "GHS05", "GHS06", "GHS08"],
        storage="Prepared fresh in acetone", physical_state="liquid",
    ),
    "mcpba": Reagent(
        name="meta-Chloroperoxybenzoic acid", formula="ClC6H4CO3H", cas="937-14-4",
        mw=172.57, cost_per_kg=220,
        ghs_hazards=["H242", "H302", "H315", "H318", "H335"],
        ghs_pictograms=["GHS02", "GHS05", "GHS07"],
        storage="Freezer, away from organics and metals", physical_state="solid",
    ),
    "oso4": Reagent(
        name="Osmium tetroxide", formula="OsO4", cas="20816-12-0",
        mw=254.23, cost_per_kg=25000,
        ghs_hazards=["H300", "H310", "H314", "H330"],
        ghs_pictograms=["GHS05", "GHS06"],
        storage="Double-sealed container, fume hood only", physical_state="solid",
    ),
    "naio4": Reagent(
        name="Sodium periodate", formula="NaIO4", cas="7790-28-5",
        mw=213.89, cost_per_kg=110,
        ghs_hazards=["H272", "H302", "H315", "H319"],
        ghs_pictograms=["GHS03", "GHS07"],
        storage="Cool dry place", physical_state="solid",
    ),
    "h2o2": Reagent(
        name="Hydrogen peroxide 30%", formula="H2O2", cas="7722-84-1",
        mw=34.01, cost_per_kg=12,
        ghs_hazards=["H271", "H302", "H314", "H318", "H335"],
        ghs_pictograms=["GHS03", "GHS05", "GHS07"],
        storage="Fridge, away from metals and organics", physical_state="liquid",
    ),
    "naocl": Reagent(
        name="Sodium hypochlorite (bleach)", formula="NaOCl", cas="7681-52-9",
        mw=74.44, cost_per_kg=5,
        ghs_hazards=["H314", "H400", "H410"],
        ghs_pictograms=["GHS05", "GHS09"],
        storage="Cool dark place", physical_state="liquid",
    ),
    "dess_martin": Reagent(
        name="Dess-Martin periodinane", formula="C13H13IO8", cas="87413-09-0",
        mw=424.14, cost_per_kg=3500,
        ghs_hazards=["H228", "H302", "H315", "H319"],
        ghs_pictograms=["GHS02", "GHS07"],
        storage="Freezer, moisture-sensitive", physical_state="solid",
    ),

    # -----------------------------------------------------------------
    #  Reducing agents
    # -----------------------------------------------------------------
    "nabh4": Reagent(
        name="Sodium borohydride", formula="NaBH4", cas="16940-66-2",
        mw=37.83, cost_per_kg=80,
        ghs_hazards=["H260", "H301", "H314"],
        ghs_pictograms=["GHS02", "GHS05", "GHS06"],
        storage="Dry, inert atmosphere, away from water", physical_state="solid",
    ),
    "lialh4": Reagent(
        name="Lithium aluminium hydride", formula="LiAlH4", cas="16853-85-3",
        mw=37.95, cost_per_kg=250,
        ghs_hazards=["H260", "H314"],
        ghs_pictograms=["GHS02", "GHS05"],
        storage="Glovebox or sure-seal, strict inert atmosphere", physical_state="solid",
    ),
    "dibal_h": Reagent(
        name="Diisobutylaluminium hydride", formula="(i-Bu)2AlH", cas="1191-15-7",
        mw=142.22, cost_per_kg=320,
        ghs_hazards=["H250", "H260", "H314"],
        ghs_pictograms=["GHS02", "GHS05"],
        storage="Sure-seal, strict inert atmosphere, -20 C", physical_state="liquid",
    ),
    "h2_pd_c": Reagent(
        name="Hydrogen / Palladium on carbon", formula="H2 + Pd/C",
        mw=0.0, cost_per_kg=800,
        ghs_hazards=["H220", "H315", "H319", "H335"],
        ghs_pictograms=["GHS02", "GHS07"],
        storage="Pd/C: keep wet to avoid pyrophoricity; H2: gas cylinder", physical_state="solid",
    ),
    "zn_hcl": Reagent(
        name="Zinc / Hydrochloric acid", formula="Zn/HCl",
        mw=0.0, cost_per_kg=15,
        ghs_hazards=["H290", "H314", "H335", "H410"],
        ghs_pictograms=["GHS05", "GHS07", "GHS09"],
        storage="Zinc: dry; HCl: corrosive cabinet", physical_state="solid",
    ),
    "na_nh3": Reagent(
        name="Sodium in liquid ammonia (Birch)", formula="Na/NH3(l)",
        mw=0.0, cost_per_kg=50,
        ghs_hazards=["H250", "H260", "H314", "H331"],
        ghs_pictograms=["GHS02", "GHS05", "GHS06"],
        storage="Na under oil; NH3 in cylinder", physical_state="liquid",
    ),
    "red_al": Reagent(
        name="Red-Al (sodium bis(2-methoxyethoxy)aluminium hydride)",
        formula="NaAlH2(OCH2CH2OCH3)2", cas="22722-98-1",
        mw=202.16, cost_per_kg=400,
        ghs_hazards=["H250", "H260", "H314"],
        ghs_pictograms=["GHS02", "GHS05"],
        storage="Sure-seal, inert atmosphere", physical_state="liquid",
    ),
    "cecl3": Reagent(
        name="Cerium(III) chloride heptahydrate", formula="CeCl3*7H2O", cas="18618-55-8",
        mw=372.58, cost_per_kg=90,
        ghs_hazards=["H302", "H315", "H319"],
        ghs_pictograms=["GHS07"],
        storage="Cool dry place", physical_state="solid",
    ),

    # -----------------------------------------------------------------
    #  Nucleophiles / Bases (inorganic)
    # -----------------------------------------------------------------
    "naoh": Reagent(
        name="Sodium hydroxide", formula="NaOH", cas="1310-73-2",
        mw=40.00, cost_per_kg=8,
        ghs_hazards=["H290", "H314"],
        ghs_pictograms=["GHS05"],
        storage="Tightly sealed, away from acids", physical_state="solid",
    ),
    "koh": Reagent(
        name="Potassium hydroxide", formula="KOH", cas="1310-58-3",
        mw=56.11, cost_per_kg=12,
        ghs_hazards=["H290", "H302", "H314"],
        ghs_pictograms=["GHS05", "GHS07"],
        storage="Tightly sealed, away from acids", physical_state="solid",
    ),
    "nacn": Reagent(
        name="Sodium cyanide", formula="NaCN", cas="143-33-9",
        mw=49.01, cost_per_kg=60,
        ghs_hazards=["H290", "H300", "H310", "H330", "H410"],
        ghs_pictograms=["GHS06", "GHS09"],
        storage="Locked poison cabinet, separate from acids", physical_state="solid",
    ),
    "nan3": Reagent(
        name="Sodium azide", formula="NaN3", cas="26628-22-8",
        mw=65.01, cost_per_kg=55,
        ghs_hazards=["H300", "H310", "H400", "H410"],
        ghs_pictograms=["GHS06", "GHS09"],
        storage="Poison cabinet, never contact heavy metals or acid", physical_state="solid",
    ),
    "naome": Reagent(
        name="Sodium methoxide", formula="NaOCH3", cas="124-41-4",
        mw=54.02, cost_per_kg=40,
        ghs_hazards=["H228", "H251", "H301", "H311", "H314", "H331"],
        ghs_pictograms=["GHS02", "GHS05", "GHS06"],
        storage="Inert atmosphere, away from moisture", physical_state="solid",
    ),
    "naoet": Reagent(
        name="Sodium ethoxide", formula="NaOC2H5", cas="141-52-6",
        mw=68.05, cost_per_kg=45,
        ghs_hazards=["H228", "H251", "H301", "H311", "H314", "H331"],
        ghs_pictograms=["GHS02", "GHS05", "GHS06"],
        storage="Inert atmosphere, away from moisture", physical_state="solid",
    ),
    "kotbu": Reagent(
        name="Potassium tert-butoxide", formula="KOC(CH3)3", cas="865-47-4",
        mw=112.21, cost_per_kg=90,
        ghs_hazards=["H228", "H251", "H301", "H314"],
        ghs_pictograms=["GHS02", "GHS05", "GHS06"],
        storage="Glovebox or inert atmosphere", physical_state="solid",
    ),

    # -----------------------------------------------------------------
    #  Electrophiles (alkyl/acyl halides)
    # -----------------------------------------------------------------
    "mei": Reagent(
        name="Methyl iodide", formula="CH3I", cas="74-88-4",
        mw=141.94, cost_per_kg=120,
        ghs_hazards=["H301", "H312", "H315", "H331", "H335", "H351"],
        ghs_pictograms=["GHS06", "GHS08"],
        storage="Fridge, amber bottle, inert atmosphere", physical_state="liquid",
    ),
    "etbr": Reagent(
        name="Ethyl bromide", formula="C2H5Br", cas="74-96-4",
        mw=108.97, cost_per_kg=60,
        ghs_hazards=["H225", "H302", "H332", "H351"],
        ghs_pictograms=["GHS02", "GHS07", "GHS08"],
        storage="Fridge, away from ignition", physical_state="liquid",
    ),
    "bnbr": Reagent(
        name="Benzyl bromide", formula="C6H5CH2Br", cas="100-39-0",
        mw=171.04, cost_per_kg=90,
        ghs_hazards=["H301", "H315", "H318", "H335"],
        ghs_pictograms=["GHS05", "GHS06"],
        storage="Fridge, tightly sealed", physical_state="liquid",
    ),
    "allyl_bromide": Reagent(
        name="Allyl bromide", formula="CH2=CHCH2Br", cas="106-95-6",
        mw=120.98, cost_per_kg=55,
        ghs_hazards=["H225", "H301", "H311", "H314", "H331"],
        ghs_pictograms=["GHS02", "GHS05", "GHS06"],
        storage="Fridge, inert atmosphere", physical_state="liquid",
    ),
    "accl": Reagent(
        name="Acetyl chloride", formula="CH3COCl", cas="75-36-5",
        mw=78.50, cost_per_kg=30,
        ghs_hazards=["H225", "H314", "H331"],
        ghs_pictograms=["GHS02", "GHS05", "GHS06"],
        storage="Inert atmosphere, moisture-sensitive", physical_state="liquid",
    ),
    "bzcl": Reagent(
        name="Benzoyl chloride", formula="C6H5COCl", cas="98-88-4",
        mw=140.57, cost_per_kg=40,
        ghs_hazards=["H302", "H314", "H332"],
        ghs_pictograms=["GHS05", "GHS07"],
        storage="Inert atmosphere, moisture-sensitive", physical_state="liquid",
    ),

    # -----------------------------------------------------------------
    #  Organometallics
    # -----------------------------------------------------------------
    "memgbr": Reagent(
        name="Methylmagnesium bromide (3 M in Et2O)", formula="CH3MgBr",
        cas="75-16-1", mw=119.24, cost_per_kg=600,
        ghs_hazards=["H225", "H250", "H260", "H314"],
        ghs_pictograms=["GHS02", "GHS05"],
        storage="Sure-seal, strict inert atmosphere, fridge", physical_state="liquid",
    ),
    "etmgbr": Reagent(
        name="Ethylmagnesium bromide (3 M in Et2O)", formula="C2H5MgBr",
        cas="925-90-6", mw=133.27, cost_per_kg=550,
        ghs_hazards=["H225", "H250", "H260", "H314"],
        ghs_pictograms=["GHS02", "GHS05"],
        storage="Sure-seal, strict inert atmosphere, fridge", physical_state="liquid",
    ),
    "phmgbr": Reagent(
        name="Phenylmagnesium bromide (3 M in THF)", formula="C6H5MgBr",
        cas="100-58-3", mw=181.31, cost_per_kg=500,
        ghs_hazards=["H225", "H250", "H260", "H314"],
        ghs_pictograms=["GHS02", "GHS05"],
        storage="Sure-seal, strict inert atmosphere, fridge", physical_state="liquid",
    ),
    "n_buli": Reagent(
        name="n-Butyllithium (2.5 M in hexanes)", formula="C4H9Li",
        cas="109-72-8", mw=64.06, cost_per_kg=450,
        ghs_hazards=["H225", "H250", "H260", "H304", "H314"],
        ghs_pictograms=["GHS02", "GHS05"],
        storage="Sure-seal, glovebox, -20 C", physical_state="liquid",
    ),
    "phli": Reagent(
        name="Phenyllithium (1.9 M in Bu2O)", formula="C6H5Li",
        cas="591-51-5", mw=84.04, cost_per_kg=700,
        ghs_hazards=["H225", "H250", "H260", "H314"],
        ghs_pictograms=["GHS02", "GHS05"],
        storage="Sure-seal, glovebox", physical_state="liquid",
    ),
    "meli": Reagent(
        name="Methyllithium (1.6 M in Et2O)", formula="CH3Li",
        cas="917-54-4", mw=21.98, cost_per_kg=800,
        ghs_hazards=["H225", "H250", "H260", "H314"],
        ghs_pictograms=["GHS02", "GHS05"],
        storage="Sure-seal, glovebox", physical_state="liquid",
    ),

    # -----------------------------------------------------------------
    #  Coupling catalysts / reagents
    # -----------------------------------------------------------------
    "pd_pph3_4": Reagent(
        name="Tetrakis(triphenylphosphine)palladium(0)", formula="Pd(PPh3)4",
        cas="14221-01-3", mw=1155.56, cost_per_kg=12000,
        ghs_hazards=["H315", "H319", "H335"],
        ghs_pictograms=["GHS07"],
        storage="Glovebox or Schlenk line, -20 C, light-sensitive", physical_state="solid",
    ),
    "pdcl2_dppf": Reagent(
        name="[1,1'-Bis(diphenylphosphino)ferrocene]dichloropalladium(II)",
        formula="PdCl2(dppf)", cas="72287-26-4", mw=731.70, cost_per_kg=15000,
        ghs_hazards=["H302", "H315", "H319", "H335"],
        ghs_pictograms=["GHS07"],
        storage="Inert atmosphere, light-sensitive", physical_state="solid",
    ),
    "pd_oac_2": Reagent(
        name="Palladium(II) acetate", formula="Pd(OAc)2", cas="3375-31-3",
        mw=224.51, cost_per_kg=10000,
        ghs_hazards=["H302", "H315", "H319"],
        ghs_pictograms=["GHS07"],
        storage="Cool dry place, light-sensitive", physical_state="solid",
    ),
    "cui": Reagent(
        name="Copper(I) iodide", formula="CuI", cas="7681-65-4",
        mw=190.45, cost_per_kg=120,
        ghs_hazards=["H302", "H315", "H319", "H335", "H410"],
        ghs_pictograms=["GHS07", "GHS09"],
        storage="Cool dry place, light-sensitive", physical_state="solid",
    ),
    "ni_cod_2": Reagent(
        name="Bis(1,5-cyclooctadiene)nickel(0)", formula="Ni(cod)2",
        cas="1295-35-8", mw=275.04, cost_per_kg=18000,
        ghs_hazards=["H250", "H302", "H315", "H317", "H319", "H334", "H341", "H350", "H360", "H372"],
        ghs_pictograms=["GHS02", "GHS07", "GHS08"],
        storage="Glovebox, -20 C, extremely air-sensitive", physical_state="solid",
    ),

    # -----------------------------------------------------------------
    #  Bases (organic)
    # -----------------------------------------------------------------
    "et3n": Reagent(
        name="Triethylamine", formula="(C2H5)3N", cas="121-44-8",
        mw=101.19, cost_per_kg=22,
        ghs_hazards=["H225", "H302", "H311", "H314", "H331"],
        ghs_pictograms=["GHS02", "GHS05", "GHS06"],
        storage="Flammable cabinet, away from acids", physical_state="liquid",
    ),
    "dmap": Reagent(
        name="4-Dimethylaminopyridine", formula="C7H10N2", cas="1122-58-3",
        mw=122.17, cost_per_kg=250,
        ghs_hazards=["H301", "H311", "H331"],
        ghs_pictograms=["GHS06"],
        storage="Cool dry place", physical_state="solid",
    ),
    "dbu": Reagent(
        name="1,8-Diazabicyclo[5.4.0]undec-7-ene", formula="C9H16N2",
        cas="6674-22-2", mw=152.24, cost_per_kg=110,
        ghs_hazards=["H302", "H311", "H314"],
        ghs_pictograms=["GHS05", "GHS07"],
        storage="Inert atmosphere", physical_state="liquid",
    ),
    "nah": Reagent(
        name="Sodium hydride (60% in mineral oil)", formula="NaH",
        cas="7646-69-7", mw=24.00, cost_per_kg=65,
        ghs_hazards=["H260", "H314"],
        ghs_pictograms=["GHS02", "GHS05"],
        storage="Glovebox or Schlenk line", physical_state="solid",
    ),
    "k2co3": Reagent(
        name="Potassium carbonate", formula="K2CO3", cas="584-08-7",
        mw=138.21, cost_per_kg=10,
        ghs_hazards=["H315", "H319", "H335"],
        ghs_pictograms=["GHS07"],
        storage="Dry, tightly sealed", physical_state="solid",
    ),
    "cs2co3": Reagent(
        name="Cesium carbonate", formula="Cs2CO3", cas="534-17-8",
        mw=325.82, cost_per_kg=400,
        ghs_hazards=["H315", "H319", "H335"],
        ghs_pictograms=["GHS07"],
        storage="Dry, tightly sealed, hygroscopic", physical_state="solid",
    ),
    "pyridine": Reagent(
        name="Pyridine", formula="C5H5N", cas="110-86-1",
        mw=79.10, cost_per_kg=25,
        ghs_hazards=["H225", "H302", "H312", "H332"],
        ghs_pictograms=["GHS02", "GHS07"],
        storage="Flammable cabinet", physical_state="liquid",
    ),
    "dipea": Reagent(
        name="N,N-Diisopropylethylamine (Hunig's base)", formula="(iPr)2NEt",
        cas="7087-68-5", mw=129.25, cost_per_kg=50,
        ghs_hazards=["H225", "H302", "H311", "H314", "H331"],
        ghs_pictograms=["GHS02", "GHS05", "GHS06"],
        storage="Flammable cabinet, inert atmosphere", physical_state="liquid",
    ),

    # -----------------------------------------------------------------
    #  Acids
    # -----------------------------------------------------------------
    "hcl": Reagent(
        name="Hydrochloric acid (conc. 37%)", formula="HCl", cas="7647-01-0",
        mw=36.46, cost_per_kg=6,
        ghs_hazards=["H290", "H314", "H335"],
        ghs_pictograms=["GHS05", "GHS07"],
        storage="Acid cabinet", physical_state="liquid",
    ),
    "h2so4": Reagent(
        name="Sulfuric acid (conc.)", formula="H2SO4", cas="7664-93-9",
        mw=98.08, cost_per_kg=8,
        ghs_hazards=["H290", "H314"],
        ghs_pictograms=["GHS05"],
        storage="Acid cabinet, separate from organics", physical_state="liquid",
    ),
    "hno3": Reagent(
        name="Nitric acid (conc. 70%)", formula="HNO3", cas="7697-37-2",
        mw=63.01, cost_per_kg=10,
        ghs_hazards=["H272", "H290", "H314"],
        ghs_pictograms=["GHS03", "GHS05"],
        storage="Acid cabinet, separate from organics", physical_state="liquid",
    ),
    "tfa": Reagent(
        name="Trifluoroacetic acid", formula="CF3COOH", cas="76-05-1",
        mw=114.02, cost_per_kg=55,
        ghs_hazards=["H290", "H314", "H332"],
        ghs_pictograms=["GHS05", "GHS07"],
        storage="Acid cabinet, tightly sealed", physical_state="liquid",
    ),
    "acoh": Reagent(
        name="Acetic acid (glacial)", formula="CH3COOH", cas="64-19-7",
        mw=60.05, cost_per_kg=12,
        ghs_hazards=["H226", "H314"],
        ghs_pictograms=["GHS02", "GHS05"],
        storage="Acid/flammable cabinet", physical_state="liquid",
    ),
    "p_tsoh": Reagent(
        name="para-Toluenesulfonic acid monohydrate", formula="CH3C6H4SO3H*H2O",
        cas="6192-52-5", mw=190.22, cost_per_kg=25,
        ghs_hazards=["H302", "H314", "H318"],
        ghs_pictograms=["GHS05", "GHS07"],
        storage="Dry, tightly sealed", physical_state="solid",
    ),
    "bf3_oet2": Reagent(
        name="Boron trifluoride diethyl etherate", formula="BF3*OEt2",
        cas="109-63-7", mw=141.93, cost_per_kg=45,
        ghs_hazards=["H226", "H302", "H314", "H331"],
        ghs_pictograms=["GHS02", "GHS05", "GHS06"],
        storage="Inert atmosphere, moisture-sensitive", physical_state="liquid",
    ),
    "ticl4": Reagent(
        name="Titanium(IV) chloride", formula="TiCl4", cas="7550-45-0",
        mw=189.68, cost_per_kg=60,
        ghs_hazards=["H290", "H314", "H330", "H335"],
        ghs_pictograms=["GHS05", "GHS06"],
        storage="Inert atmosphere, fume hood, moisture-sensitive", physical_state="liquid",
    ),
    "alcl3": Reagent(
        name="Aluminium chloride", formula="AlCl3", cas="7446-70-0",
        mw=133.34, cost_per_kg=20,
        ghs_hazards=["H290", "H302", "H314"],
        ghs_pictograms=["GHS05", "GHS07"],
        storage="Inert atmosphere, moisture-sensitive", physical_state="solid",
    ),

    # -----------------------------------------------------------------
    #  Protection reagents
    # -----------------------------------------------------------------
    "boc2o": Reagent(
        name="Di-tert-butyl dicarbonate", formula="(Boc)2O", cas="24424-99-5",
        mw=218.25, cost_per_kg=100,
        ghs_hazards=["H228", "H302"],
        ghs_pictograms=["GHS02", "GHS07"],
        storage="Fridge", physical_state="solid",
    ),
    "fmoc_cl": Reagent(
        name="9-Fluorenylmethyl chloroformate", formula="FmocCl",
        cas="28920-43-6", mw=258.70, cost_per_kg=350,
        ghs_hazards=["H302", "H315", "H319"],
        ghs_pictograms=["GHS07"],
        storage="Fridge, moisture-sensitive", physical_state="solid",
    ),
    "tbscl": Reagent(
        name="tert-Butyldimethylsilyl chloride", formula="TBSCl",
        cas="18162-48-6", mw=150.72, cost_per_kg=120,
        ghs_hazards=["H228", "H302", "H314"],
        ghs_pictograms=["GHS02", "GHS05", "GHS07"],
        storage="Dry, inert atmosphere", physical_state="solid",
    ),
    "tbdpscl": Reagent(
        name="tert-Butyldiphenylsilyl chloride", formula="TBDPSCl",
        cas="58479-61-1", mw=274.87, cost_per_kg=300,
        ghs_hazards=["H302", "H314"],
        ghs_pictograms=["GHS05", "GHS07"],
        storage="Inert atmosphere, moisture-sensitive", physical_state="liquid",
    ),
    "pmbcl": Reagent(
        name="para-Methoxybenzyl chloride", formula="PMBCl", cas="824-94-2",
        mw=156.61, cost_per_kg=80,
        ghs_hazards=["H302", "H315", "H319"],
        ghs_pictograms=["GHS07"],
        storage="Fridge, moisture-sensitive", physical_state="liquid",
    ),
    "ac2o": Reagent(
        name="Acetic anhydride", formula="(CH3CO)2O", cas="108-24-7",
        mw=102.09, cost_per_kg=15,
        ghs_hazards=["H226", "H302", "H314", "H332"],
        ghs_pictograms=["GHS02", "GHS05", "GHS07"],
        storage="Flammable cabinet, moisture-sensitive", physical_state="liquid",
    ),

    # -----------------------------------------------------------------
    #  Deprotection reagents
    # -----------------------------------------------------------------
    "tbaf": Reagent(
        name="Tetrabutylammonium fluoride (1 M in THF)", formula="TBAF",
        cas="429-41-4", mw=261.46, cost_per_kg=500,
        ghs_hazards=["H225", "H302", "H314"],
        ghs_pictograms=["GHS02", "GHS05", "GHS07"],
        storage="Sure-seal, fridge", physical_state="liquid",
    ),
    "hf_pyridine": Reagent(
        name="Hydrogen fluoride / pyridine (Olah's reagent)", formula="HF*pyr",
        cas="32001-55-1", mw=0.0, cost_per_kg=200,
        ghs_hazards=["H300", "H310", "H314", "H330"],
        ghs_pictograms=["GHS05", "GHS06"],
        storage="Polyethylene container, CaF2 antidote nearby", physical_state="liquid",
    ),
    "ddq": Reagent(
        name="2,3-Dichloro-5,6-dicyano-1,4-benzoquinone", formula="DDQ",
        cas="84-58-2", mw=227.00, cost_per_kg=500,
        ghs_hazards=["H301", "H311", "H319", "H331"],
        ghs_pictograms=["GHS06"],
        storage="Cool dry place, light-sensitive", physical_state="solid",
    ),

    # -----------------------------------------------------------------
    #  Miscellaneous
    # -----------------------------------------------------------------
    "pph3": Reagent(
        name="Triphenylphosphine", formula="PPh3", cas="603-35-0",
        mw=262.29, cost_per_kg=50,
        ghs_hazards=["H302", "H315", "H319"],
        ghs_pictograms=["GHS07"],
        storage="Room temperature", physical_state="solid",
    ),
    "grubbs_1": Reagent(
        name="Grubbs catalyst 1st generation", formula="Cl2(PCy3)2Ru=CHPh",
        cas="172222-30-9", mw=822.96, cost_per_kg=25000,
        ghs_hazards=["H315", "H319", "H335"],
        ghs_pictograms=["GHS07"],
        storage="Glovebox, -20 C, inert atmosphere", physical_state="solid",
    ),
    "grubbs_2": Reagent(
        name="Grubbs catalyst 2nd generation",
        formula="Cl2(IMes)(PCy3)Ru=CHPh", cas="246047-72-3",
        mw=848.97, cost_per_kg=28000,
        ghs_hazards=["H315", "H319", "H335"],
        ghs_pictograms=["GHS07"],
        storage="Glovebox, -20 C, inert atmosphere", physical_state="solid",
    ),
    "socl2": Reagent(
        name="Thionyl chloride", formula="SOCl2", cas="7719-09-7",
        mw=118.97, cost_per_kg=25,
        ghs_hazards=["H302", "H314", "H331"],
        ghs_pictograms=["GHS05", "GHS06"],
        storage="Inert atmosphere, moisture-sensitive", physical_state="liquid",
    ),
    "pcl5": Reagent(
        name="Phosphorus pentachloride", formula="PCl5", cas="10026-13-8",
        mw=208.24, cost_per_kg=40,
        ghs_hazards=["H290", "H302", "H314"],
        ghs_pictograms=["GHS05", "GHS07"],
        storage="Inert atmosphere, moisture-sensitive", physical_state="solid",
    ),
    "p2o5": Reagent(
        name="Phosphorus pentoxide", formula="P2O5", cas="1314-56-3",
        mw=141.94, cost_per_kg=20,
        ghs_hazards=["H290", "H314"],
        ghs_pictograms=["GHS05"],
        storage="Tightly sealed, extremely hygroscopic", physical_state="solid",
    ),
    "diad": Reagent(
        name="Diisopropyl azodicarboxylate", formula="iPrO2CN=NCO2iPr",
        cas="2446-83-5", mw=202.21, cost_per_kg=180,
        ghs_hazards=["H225", "H302", "H312", "H315", "H319", "H332"],
        ghs_pictograms=["GHS02", "GHS07"],
        storage="Fridge, flammable cabinet", physical_state="liquid",
    ),
    "dead": Reagent(
        name="Diethyl azodicarboxylate", formula="EtO2CN=NCO2Et",
        cas="1972-28-7", mw=174.15, cost_per_kg=160,
        ghs_hazards=["H225", "H301", "H311", "H315", "H319", "H331", "H335", "H351"],
        ghs_pictograms=["GHS02", "GHS06", "GHS08"],
        storage="Fridge, flammable cabinet", physical_state="liquid",
    ),
    "imidazole": Reagent(
        name="Imidazole", formula="C3H4N2", cas="288-32-4",
        mw=68.08, cost_per_kg=30,
        ghs_hazards=["H302", "H314"],
        ghs_pictograms=["GHS05", "GHS07"],
        storage="Room temperature", physical_state="solid",
    ),
    "n2h4": Reagent(
        name="Hydrazine monohydrate", formula="N2H4*H2O", cas="7803-57-8",
        mw=50.06, cost_per_kg=40,
        ghs_hazards=["H226", "H301", "H311", "H314", "H317", "H331", "H350", "H410"],
        ghs_pictograms=["GHS02", "GHS05", "GHS06", "GHS08", "GHS09"],
        storage="Flammable cabinet, poison area", physical_state="liquid",
    ),
    "zncl2": Reagent(
        name="Zinc chloride", formula="ZnCl2", cas="7646-85-7",
        mw=136.30, cost_per_kg=15,
        ghs_hazards=["H290", "H302", "H314", "H410"],
        ghs_pictograms=["GHS05", "GHS07", "GHS09"],
        storage="Dry, tightly sealed", physical_state="solid",
    ),
    "znbr2": Reagent(
        name="Zinc bromide", formula="ZnBr2", cas="7699-45-8",
        mw=225.20, cost_per_kg=50,
        ghs_hazards=["H302", "H314", "H410"],
        ghs_pictograms=["GHS05", "GHS07", "GHS09"],
        storage="Dry, tightly sealed, hygroscopic", physical_state="solid",
    ),
    "cul": Reagent(
        name="Copper(I) chloride", formula="CuCl", cas="7758-89-6",
        mw=98.99, cost_per_kg=30,
        ghs_hazards=["H302", "H400", "H410"],
        ghs_pictograms=["GHS07", "GHS09"],
        storage="Inert atmosphere", physical_state="solid",
    ),
    "mgso4": Reagent(
        name="Magnesium sulfate (anhydrous, drying agent)", formula="MgSO4",
        cas="7487-88-9", mw=120.37, cost_per_kg=8,
        ghs_hazards=[], ghs_pictograms=[],
        storage="Tightly sealed", physical_state="solid",
    ),
    "na2so4": Reagent(
        name="Sodium sulfate (anhydrous, drying agent)", formula="Na2SO4",
        cas="7757-82-6", mw=142.04, cost_per_kg=5,
        ghs_hazards=[], ghs_pictograms=[],
        storage="Tightly sealed", physical_state="solid",
    ),
    "nahco3": Reagent(
        name="Sodium bicarbonate", formula="NaHCO3", cas="144-55-8",
        mw=84.01, cost_per_kg=4,
        ghs_hazards=[], ghs_pictograms=[],
        storage="Room temperature", physical_state="solid",
    ),
    "nh4cl": Reagent(
        name="Ammonium chloride", formula="NH4Cl", cas="12125-02-9",
        mw=53.49, cost_per_kg=6,
        ghs_hazards=["H302", "H319"],
        ghs_pictograms=["GHS07"],
        storage="Room temperature", physical_state="solid",
    ),
    "lda": Reagent(
        name="Lithium diisopropylamide (2 M in THF/heptane)", formula="LDA",
        cas="4111-54-0", mw=107.12, cost_per_kg=500,
        ghs_hazards=["H225", "H250", "H260", "H314"],
        ghs_pictograms=["GHS02", "GHS05"],
        storage="Sure-seal, -20 C", physical_state="liquid",
    ),
    "lhmds": Reagent(
        name="Lithium bis(trimethylsilyl)amide (1 M in THF)", formula="LiHMDS",
        cas="4039-32-1", mw=167.33, cost_per_kg=450,
        ghs_hazards=["H225", "H250", "H260", "H314"],
        ghs_pictograms=["GHS02", "GHS05"],
        storage="Sure-seal, -20 C", physical_state="liquid",
    ),
    "nahmds": Reagent(
        name="Sodium bis(trimethylsilyl)amide (1 M in THF)", formula="NaHMDS",
        cas="1070-89-9", mw=183.37, cost_per_kg=400,
        ghs_hazards=["H225", "H250", "H260", "H314"],
        ghs_pictograms=["GHS02", "GHS05"],
        storage="Sure-seal, -20 C", physical_state="liquid",
    ),
    "khmds": Reagent(
        name="Potassium bis(trimethylsilyl)amide (0.5 M in toluene)",
        formula="KHMDS", cas="40949-94-8", mw=199.47, cost_per_kg=500,
        ghs_hazards=["H225", "H250", "H260", "H314"],
        ghs_pictograms=["GHS02", "GHS05"],
        storage="Sure-seal, -20 C", physical_state="liquid",
    ),
    "o3": Reagent(
        name="Ozone (generated in situ)", formula="O3",
        mw=48.00, cost_per_kg=0,
        ghs_hazards=["H270", "H314", "H330"],
        ghs_pictograms=["GHS03", "GHS05", "GHS06"],
        storage="Generated from O2 via ozonizer", physical_state="gas",
    ),
    "me2s": Reagent(
        name="Dimethyl sulfide", formula="(CH3)2S", cas="75-18-3",
        mw=62.13, cost_per_kg=15,
        ghs_hazards=["H225", "H302", "H332"],
        ghs_pictograms=["GHS02", "GHS07"],
        storage="Flammable cabinet, fume hood (stench)", physical_state="liquid",
    ),
    "ph3p_chco2et": Reagent(
        name="Ethyl (triphenylphosphoranylidene)acetate (stabilized Wittig)",
        formula="Ph3P=CHCO2Et", cas="1099-45-2", mw=348.39, cost_per_kg=300,
        ghs_hazards=["H302", "H315", "H319"],
        ghs_pictograms=["GHS07"],
        storage="Room temperature", physical_state="solid",
    ),
    "ch2i2": Reagent(
        name="Diiodomethane", formula="CH2I2", cas="75-11-6",
        mw=267.84, cost_per_kg=80,
        ghs_hazards=["H302", "H315", "H319", "H335"],
        ghs_pictograms=["GHS07"],
        storage="Amber bottle, fridge", physical_state="liquid",
    ),
    "znet2": Reagent(
        name="Diethylzinc (1 M in hexanes)", formula="ZnEt2",
        cas="557-20-0", mw=123.53, cost_per_kg=350,
        ghs_hazards=["H250", "H260", "H302", "H314"],
        ghs_pictograms=["GHS02", "GHS05", "GHS07"],
        storage="Sure-seal, glovebox", physical_state="liquid",
    ),
    "bu3snh": Reagent(
        name="Tributyltin hydride", formula="Bu3SnH", cas="688-73-3",
        mw=291.06, cost_per_kg=200,
        ghs_hazards=["H301", "H311", "H315", "H319", "H331", "H360", "H372", "H410"],
        ghs_pictograms=["GHS06", "GHS08", "GHS09"],
        storage="Fridge, inert atmosphere", physical_state="liquid",
    ),
    "aibn": Reagent(
        name="Azobisisobutyronitrile", formula="(CH3)2C(CN)N=NC(CN)(CH3)2",
        cas="78-67-1", mw=164.21, cost_per_kg=50,
        ghs_hazards=["H242", "H302", "H412"],
        ghs_pictograms=["GHS02", "GHS07"],
        storage="Freezer, shock-sensitive", physical_state="solid",
    ),
    "nbs": Reagent(
        name="N-Bromosuccinimide", formula="C4H4BrNO2", cas="128-08-5",
        mw=177.98, cost_per_kg=45,
        ghs_hazards=["H302", "H315", "H319"],
        ghs_pictograms=["GHS07"],
        storage="Cool dry place, light-sensitive", physical_state="solid",
    ),
    "ncs": Reagent(
        name="N-Chlorosuccinimide", formula="C4H4ClNO2", cas="128-09-6",
        mw=133.53, cost_per_kg=40,
        ghs_hazards=["H302", "H314"],
        ghs_pictograms=["GHS05", "GHS07"],
        storage="Cool dry place", physical_state="solid",
    ),
    "br2": Reagent(
        name="Bromine", formula="Br2", cas="7726-95-6",
        mw=159.81, cost_per_kg=30,
        ghs_hazards=["H314", "H330", "H400"],
        ghs_pictograms=["GHS05", "GHS06", "GHS09"],
        storage="Corrosive cabinet, fume hood", physical_state="liquid",
    ),
    "i2": Reagent(
        name="Iodine", formula="I2", cas="7553-56-2",
        mw=253.81, cost_per_kg=60,
        ghs_hazards=["H302", "H312", "H332", "H400"],
        ghs_pictograms=["GHS07", "GHS09"],
        storage="Cool dry place", physical_state="solid",
    ),
    "9_bbn": Reagent(
        name="9-Borabicyclo[3.3.1]nonane (0.5 M in THF)", formula="9-BBN",
        cas="280-64-8", mw=122.02, cost_per_kg=400,
        ghs_hazards=["H225", "H250", "H260", "H314"],
        ghs_pictograms=["GHS02", "GHS05"],
        storage="Sure-seal, inert atmosphere", physical_state="liquid",
    ),
    "bh3_thf": Reagent(
        name="Borane-THF complex (1 M in THF)", formula="BH3*THF",
        cas="14044-65-6", mw=85.94, cost_per_kg=250,
        ghs_hazards=["H225", "H250", "H260", "H314"],
        ghs_pictograms=["GHS02", "GHS05"],
        storage="Sure-seal, strict inert atmosphere, fridge", physical_state="liquid",
    ),
    "b2h6": Reagent(
        name="Diborane (gas)", formula="B2H6",
        cas="19287-45-7", mw=27.67, cost_per_kg=0,
        ghs_hazards=["H220", "H250", "H330", "H314"],
        ghs_pictograms=["GHS02", "GHS05", "GHS06"],
        storage="Gas cylinder, pyrophoric", physical_state="gas",
    ),
    "phthalimide_k": Reagent(
        name="Potassium phthalimide", formula="C8H4KNO2", cas="1074-82-4",
        mw=185.22, cost_per_kg=30,
        ghs_hazards=["H302", "H315", "H319"],
        ghs_pictograms=["GHS07"],
        storage="Room temperature", physical_state="solid",
    ),
    "nh2nh2": Reagent(
        name="Hydrazine (anhydrous)", formula="N2H4", cas="302-01-2",
        mw=32.05, cost_per_kg=50,
        ghs_hazards=["H226", "H301", "H311", "H314", "H317", "H331", "H350", "H410"],
        ghs_pictograms=["GHS02", "GHS05", "GHS06", "GHS08", "GHS09"],
        storage="Flammable, poison cabinet", physical_state="liquid",
    ),
    "crcl2": Reagent(
        name="Chromium(II) chloride", formula="CrCl2", cas="10049-05-5",
        mw=122.90, cost_per_kg=200,
        ghs_hazards=["H302", "H315", "H319", "H335"],
        ghs_pictograms=["GHS07"],
        storage="Glovebox, air-sensitive", physical_state="solid",
    ),
    "pd_c_10": Reagent(
        name="Palladium on carbon (10 wt%)", formula="Pd/C",
        cas="7440-05-3", mw=106.42, cost_per_kg=900,
        ghs_hazards=["H228"],
        ghs_pictograms=["GHS02"],
        storage="Keep wetted with water to prevent ignition", physical_state="solid",
    ),
    "lindlar": Reagent(
        name="Lindlar catalyst (Pd/CaCO3/Pb)", formula="Pd/CaCO3/Pb(OAc)2",
        mw=0.0, cost_per_kg=2000,
        ghs_hazards=["H302", "H332", "H360", "H373"],
        ghs_pictograms=["GHS07", "GHS08"],
        storage="Room temperature, keep dry", physical_state="solid",
    ),
    "sncl2": Reagent(
        name="Tin(II) chloride dihydrate", formula="SnCl2*2H2O",
        cas="10025-69-1", mw=225.65, cost_per_kg=25,
        ghs_hazards=["H302", "H315", "H319"],
        ghs_pictograms=["GHS07"],
        storage="Cool dry place", physical_state="solid",
    ),
    "samarium_ii_iodide": Reagent(
        name="Samarium(II) iodide (0.1 M in THF)", formula="SmI2",
        cas="32248-43-4", mw=404.16, cost_per_kg=5000,
        ghs_hazards=["H225", "H315", "H319"],
        ghs_pictograms=["GHS02", "GHS07"],
        storage="Sure-seal, inert atmosphere, air-sensitive", physical_state="liquid",
    ),
    "ticl3": Reagent(
        name="Titanium(III) chloride", formula="TiCl3", cas="7705-07-9",
        mw=154.23, cost_per_kg=80,
        ghs_hazards=["H302", "H314"],
        ghs_pictograms=["GHS05", "GHS07"],
        storage="Inert atmosphere", physical_state="solid",
    ),
    "pd2_dba_3": Reagent(
        name="Tris(dibenzylideneacetone)dipalladium(0)",
        formula="Pd2(dba)3", cas="51364-51-3", mw=915.72, cost_per_kg=14000,
        ghs_hazards=["H315", "H319", "H335"],
        ghs_pictograms=["GHS07"],
        storage="Glovebox, inert atmosphere, light-sensitive", physical_state="solid",
    ),
    "sphos": Reagent(
        name="2-Dicyclohexylphosphino-2',6'-dimethoxybiphenyl (SPhos)",
        formula="SPhos", cas="657408-07-6", mw=414.54, cost_per_kg=8000,
        ghs_hazards=["H302", "H315", "H319"],
        ghs_pictograms=["GHS07"],
        storage="Inert atmosphere, room temperature", physical_state="solid",
    ),
    "xphos": Reagent(
        name="2-Dicyclohexylphosphino-2',4',6'-triisopropylbiphenyl (XPhos)",
        formula="XPhos", cas="564483-18-7", mw=476.67, cost_per_kg=9000,
        ghs_hazards=["H302", "H315", "H319"],
        ghs_pictograms=["GHS07"],
        storage="Inert atmosphere, room temperature", physical_state="solid",
    ),
    "bu4nf": Reagent(
        name="Tetrabutylammonium fluoride trihydrate", formula="Bu4NF*3H2O",
        cas="22206-57-1", mw=315.51, cost_per_kg=450,
        ghs_hazards=["H302", "H314"],
        ghs_pictograms=["GHS05", "GHS07"],
        storage="Fridge, hygroscopic", physical_state="solid",
    ),
    "hbf4": Reagent(
        name="Tetrafluoroboric acid (48% in water)", formula="HBF4",
        cas="16872-11-0", mw=87.81, cost_per_kg=35,
        ghs_hazards=["H290", "H314"],
        ghs_pictograms=["GHS05"],
        storage="Acid cabinet", physical_state="liquid",
    ),
    "bpin2": Reagent(
        name="Bis(pinacolato)diboron", formula="B2pin2", cas="73183-34-3",
        mw=253.94, cost_per_kg=250,
        ghs_hazards=["H302", "H315", "H319"],
        ghs_pictograms=["GHS07"],
        storage="Room temperature, moisture-sensitive", physical_state="solid",
    ),
    "me3sn_snme3": Reagent(
        name="Hexamethylditin", formula="Me3SnSnMe3", cas="661-69-8",
        mw=327.59, cost_per_kg=1500,
        ghs_hazards=["H301", "H311", "H315", "H319", "H331", "H360", "H372", "H410"],
        ghs_pictograms=["GHS06", "GHS08", "GHS09"],
        storage="Fridge, inert atmosphere", physical_state="liquid",
    ),
}


# =====================================================================
#  SOLVENT_DB  (~30 entries)
# =====================================================================

SOLVENT_DB: dict[str, Solvent] = {
    "water": Solvent(
        name="Water", formula="H2O", bp=100.0, mp=0.0,
        density=1.00, dielectric=80.1, polarity_index=10.2,
        cas="7732-18-5", cost_per_L=0.05, green_score=1,
        ghs_hazards=[], miscible_with_water=True,
    ),
    "meoh": Solvent(
        name="Methanol", formula="CH3OH", bp=64.7, mp=-97.6,
        density=0.791, dielectric=32.7, polarity_index=5.1,
        cas="67-56-1", cost_per_L=8, green_score=3,
        ghs_hazards=["H225", "H301", "H311", "H331"],
        miscible_with_water=True,
    ),
    "etoh": Solvent(
        name="Ethanol", formula="C2H5OH", bp=78.4, mp=-114.1,
        density=0.789, dielectric=24.5, polarity_index=5.2,
        cas="64-17-5", cost_per_L=10, green_score=2,
        ghs_hazards=["H225"], miscible_with_water=True,
    ),
    "iproh": Solvent(
        name="Isopropanol", formula="(CH3)2CHOH", bp=82.6, mp=-89.5,
        density=0.786, dielectric=17.9, polarity_index=3.9,
        cas="67-63-0", cost_per_L=9, green_score=2,
        ghs_hazards=["H225", "H319", "H336"],
        miscible_with_water=True,
    ),
    "nbuoh": Solvent(
        name="n-Butanol", formula="C4H9OH", bp=117.7, mp=-89.8,
        density=0.810, dielectric=17.8, polarity_index=3.9,
        cas="71-36-3", cost_per_L=12, green_score=3,
        ghs_hazards=["H226", "H302", "H315", "H318", "H335", "H336"],
        miscible_with_water=False,
    ),
    "dcm": Solvent(
        name="Dichloromethane", formula="CH2Cl2", bp=39.6, mp=-96.7,
        density=1.325, dielectric=8.9, polarity_index=3.1,
        cas="75-09-2", cost_per_L=12, green_score=7,
        ghs_hazards=["H302", "H315", "H336", "H351"],
        miscible_with_water=False,
    ),
    "chcl3": Solvent(
        name="Chloroform", formula="CHCl3", bp=61.2, mp=-63.5,
        density=1.489, dielectric=4.8, polarity_index=4.1,
        cas="67-66-3", cost_per_L=14, green_score=8,
        ghs_hazards=["H302", "H315", "H319", "H331", "H336", "H351", "H361", "H372"],
        miscible_with_water=False,
    ),
    "ccl4": Solvent(
        name="Carbon tetrachloride", formula="CCl4", bp=76.7, mp=-22.9,
        density=1.594, dielectric=2.2, polarity_index=1.6,
        cas="56-23-5", cost_per_L=20, green_score=10,
        ghs_hazards=["H301", "H311", "H331", "H351", "H372", "H412", "H420"],
        miscible_with_water=False,
    ),
    "et2o": Solvent(
        name="Diethyl ether", formula="(C2H5)2O", bp=34.6, mp=-116.3,
        density=0.713, dielectric=4.3, polarity_index=2.8,
        cas="60-29-7", cost_per_L=15, green_score=5,
        ghs_hazards=["H224", "H302", "H336"],
        miscible_with_water=False,
    ),
    "thf": Solvent(
        name="Tetrahydrofuran", formula="C4H8O", bp=66.0, mp=-108.4,
        density=0.889, dielectric=7.5, polarity_index=4.0,
        cas="109-99-9", cost_per_L=18, green_score=5,
        ghs_hazards=["H225", "H302", "H319", "H335"],
        miscible_with_water=True,
    ),
    "1,4_dioxane": Solvent(
        name="1,4-Dioxane", formula="C4H8O2", bp=101.1, mp=11.8,
        density=1.033, dielectric=2.2, polarity_index=4.8,
        cas="123-91-1", cost_per_L=20, green_score=8,
        ghs_hazards=["H225", "H302", "H319", "H335", "H351"],
        miscible_with_water=True,
    ),
    "mtbe": Solvent(
        name="Methyl tert-butyl ether", formula="(CH3)3COCH3", bp=55.2, mp=-108.6,
        density=0.740, dielectric=2.6, polarity_index=2.5,
        cas="1634-04-4", cost_per_L=15, green_score=4,
        ghs_hazards=["H225", "H335"],
        miscible_with_water=False,
    ),
    "2_methf": Solvent(
        name="2-Methyltetrahydrofuran", formula="C5H10O", bp=80.0, mp=-136.0,
        density=0.854, dielectric=6.2, polarity_index=3.5,
        cas="96-47-9", cost_per_L=25, green_score=3,
        ghs_hazards=["H225", "H335"],
        miscible_with_water=False,
    ),
    "acetone": Solvent(
        name="Acetone", formula="(CH3)2CO", bp=56.1, mp=-94.7,
        density=0.791, dielectric=20.7, polarity_index=5.1,
        cas="67-64-1", cost_per_L=8, green_score=3,
        ghs_hazards=["H225", "H319", "H336"],
        miscible_with_water=True,
    ),
    "mek": Solvent(
        name="Methyl ethyl ketone (2-butanone)", formula="CH3COC2H5",
        bp=79.6, mp=-86.7,
        density=0.805, dielectric=18.5, polarity_index=4.7,
        cas="78-93-3", cost_per_L=10, green_score=4,
        ghs_hazards=["H225", "H319", "H336"],
        miscible_with_water=True,
    ),
    "etoac": Solvent(
        name="Ethyl acetate", formula="CH3COOC2H5", bp=77.1, mp=-83.6,
        density=0.902, dielectric=6.0, polarity_index=4.4,
        cas="141-78-6", cost_per_L=10, green_score=3,
        ghs_hazards=["H225", "H319", "H336"],
        miscible_with_water=False,
    ),
    "toluene": Solvent(
        name="Toluene", formula="C6H5CH3", bp=110.6, mp=-95.0,
        density=0.867, dielectric=2.4, polarity_index=2.4,
        cas="108-88-3", cost_per_L=10, green_score=5,
        ghs_hazards=["H225", "H304", "H315", "H336", "H361", "H373"],
        miscible_with_water=False,
    ),
    "benzene": Solvent(
        name="Benzene", formula="C6H6", bp=80.1, mp=5.5,
        density=0.879, dielectric=2.3, polarity_index=2.7,
        cas="71-43-2", cost_per_L=12, green_score=10,
        ghs_hazards=["H225", "H304", "H315", "H319", "H340", "H350", "H372"],
        miscible_with_water=False,
    ),
    "xylenes": Solvent(
        name="Xylenes (mixed isomers)", formula="C6H4(CH3)2", bp=139.0, mp=-47.0,
        density=0.860, dielectric=2.4, polarity_index=2.5,
        cas="1330-20-7", cost_per_L=10, green_score=6,
        ghs_hazards=["H226", "H304", "H312", "H315", "H319", "H332", "H335"],
        miscible_with_water=False,
    ),
    "hexanes": Solvent(
        name="Hexanes", formula="C6H14", bp=69.0, mp=-95.0,
        density=0.659, dielectric=1.9, polarity_index=0.1,
        cas="110-54-3", cost_per_L=10, green_score=5,
        ghs_hazards=["H225", "H304", "H315", "H336", "H361", "H373", "H411"],
        miscible_with_water=False,
    ),
    "pentane": Solvent(
        name="Pentane", formula="C5H12", bp=36.1, mp=-129.7,
        density=0.626, dielectric=1.8, polarity_index=0.0,
        cas="109-66-0", cost_per_L=12, green_score=5,
        ghs_hazards=["H225", "H304", "H336", "H411"],
        miscible_with_water=False,
    ),
    "heptane": Solvent(
        name="Heptane", formula="C7H16", bp=98.4, mp=-90.6,
        density=0.684, dielectric=1.9, polarity_index=0.1,
        cas="142-82-5", cost_per_L=12, green_score=4,
        ghs_hazards=["H225", "H304", "H315", "H336", "H411"],
        miscible_with_water=False,
    ),
    "petroleum_ether": Solvent(
        name="Petroleum ether (40-60)", formula="mixture", bp=50.0, mp=-70.0,
        density=0.640, dielectric=1.9, polarity_index=0.1,
        cas="8032-32-4", cost_per_L=8, green_score=5,
        ghs_hazards=["H225", "H304", "H336", "H411"],
        miscible_with_water=False,
    ),
    "dmf": Solvent(
        name="N,N-Dimethylformamide", formula="(CH3)2NCHO", bp=153.0, mp=-60.5,
        density=0.944, dielectric=36.7, polarity_index=6.4,
        cas="68-12-2", cost_per_L=15, green_score=7,
        ghs_hazards=["H226", "H312", "H319", "H332", "H360"],
        miscible_with_water=True,
    ),
    "dmso_solvent": Solvent(
        name="Dimethyl sulfoxide", formula="(CH3)2SO", bp=189.0, mp=19.0,
        density=1.100, dielectric=46.7, polarity_index=7.2,
        cas="67-68-5", cost_per_L=12, green_score=4,
        ghs_hazards=[], miscible_with_water=True,
    ),
    "nmp": Solvent(
        name="N-Methyl-2-pyrrolidone", formula="C5H9NO", bp=202.0, mp=-24.4,
        density=1.028, dielectric=32.0, polarity_index=6.7,
        cas="872-50-4", cost_per_L=20, green_score=7,
        ghs_hazards=["H315", "H319", "H335", "H360"],
        miscible_with_water=True,
    ),
    "dma": Solvent(
        name="N,N-Dimethylacetamide", formula="CH3CON(CH3)2", bp=165.0, mp=-20.0,
        density=0.937, dielectric=37.8, polarity_index=6.5,
        cas="127-19-5", cost_per_L=18, green_score=7,
        ghs_hazards=["H312", "H332", "H360"],
        miscible_with_water=True,
    ),
    "mecn": Solvent(
        name="Acetonitrile", formula="CH3CN", bp=81.6, mp=-43.8,
        density=0.786, dielectric=37.5, polarity_index=5.8,
        cas="75-05-8", cost_per_L=14, green_score=5,
        ghs_hazards=["H225", "H302", "H312", "H319", "H332"],
        miscible_with_water=True,
    ),
    "pyridine_solvent": Solvent(
        name="Pyridine", formula="C5H5N", bp=115.2, mp=-41.6,
        density=0.982, dielectric=12.4, polarity_index=5.3,
        cas="110-86-1", cost_per_L=25, green_score=7,
        ghs_hazards=["H225", "H302", "H312", "H332"],
        miscible_with_water=True,
    ),
    "et3n_solvent": Solvent(
        name="Triethylamine", formula="(C2H5)3N", bp=88.9, mp=-114.7,
        density=0.726, dielectric=2.4, polarity_index=1.8,
        cas="121-44-8", cost_per_L=22, green_score=5,
        ghs_hazards=["H225", "H302", "H311", "H314", "H331"],
        miscible_with_water=False,
    ),
    "acoh_solvent": Solvent(
        name="Acetic acid (glacial)", formula="CH3COOH", bp=117.9, mp=16.6,
        density=1.049, dielectric=6.2, polarity_index=6.0,
        cas="64-19-7", cost_per_L=12, green_score=4,
        ghs_hazards=["H226", "H314"], miscible_with_water=True,
    ),
    "tfa_solvent": Solvent(
        name="Trifluoroacetic acid", formula="CF3COOH", bp=72.4, mp=-15.4,
        density=1.489, dielectric=8.6, polarity_index=5.5,
        cas="76-05-1", cost_per_L=55, green_score=8,
        ghs_hazards=["H290", "H314", "H332"],
        miscible_with_water=True,
    ),
}


# =====================================================================
#  Lookup helpers
# =====================================================================

def get_reagent(name: str) -> Reagent | None:
    """Look up a reagent by name or normalised key."""
    key = normalize_reagent_name(name)
    return REAGENT_DB.get(key) or REAGENT_DB.get(name)


def get_solvent(name: str) -> Solvent | None:
    """Look up a solvent by name or normalised key."""
    key = normalize_reagent_name(name)
    return SOLVENT_DB.get(key) or SOLVENT_DB.get(name)


def search_reagents(query: str) -> list[Reagent]:
    """Return all reagents whose name or formula contains *query*."""
    q = query.lower()
    return [r for r in REAGENT_DB.values()
            if q in r.name.lower() or q in r.formula.lower()]


def search_solvents(query: str) -> list[Solvent]:
    """Return all solvents whose name contains *query*."""
    q = query.lower()
    return [s for s in SOLVENT_DB.values() if q in s.name.lower()]
