"""Curated reaction knowledge base with ~80 reaction templates.

Every template is a ``ReactionTemplate`` instance that captures the
reagents, conditions, selectivity cues, functional-group requirements,
and retrosynthetic disconnections needed by the synthesis planner.

Use the lookup helpers at the bottom of the module to query the
database by functional group, category, or name.
"""

from molbuilder.reactions.reaction_types import ReactionTemplate, ReactionCategory

# =====================================================================
#  Master list of reaction templates
# =====================================================================

REACTION_TEMPLATES: list[ReactionTemplate] = []

# Helper to build and register in one step
def _t(**kw) -> ReactionTemplate:
    tmpl = ReactionTemplate(**kw)
    REACTION_TEMPLATES.append(tmpl)
    return tmpl


# -----------------------------------------------------------------
#  SUBSTITUTION  (~10)
# -----------------------------------------------------------------

_t(
    name="SN2 hydroxide on primary alkyl halide",
    named_reaction="SN2",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["NaOH"],
    solvents=["DMSO", "DMF"],
    temperature_range=(25.0, 80.0),
    typical_yield=(70.0, 95.0),
    functional_group_required=["alkyl_halide"],
    functional_group_produced=["alcohol"],
    functional_group_incompatible=["ester"],
    mechanism="Back-side nucleophilic attack by hydroxide on the electrophilic carbon with simultaneous departure of the leaving group, proceeding through a single pentacoordinate transition state.",
    reverse_transform="Disconnect C-OH bond; re-attach halide to carbon.",
    safety_notes="NaOH is corrosive (H314).",
)

_t(
    name="SN2 cyanide on primary alkyl halide",
    named_reaction="SN2 cyanide displacement",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["NaCN"],
    solvents=["DMSO", "DMF"],
    temperature_range=(25.0, 80.0),
    typical_yield=(65.0, 90.0),
    functional_group_required=["alkyl_halide"],
    functional_group_produced=["nitrile"],
    functional_group_incompatible=[],
    mechanism="Cyanide ion attacks electrophilic carbon in an SN2 fashion, extending the carbon chain by one.",
    reverse_transform="Disconnect C-CN bond; replace CN with halide.",
    safety_notes="NaCN is acutely toxic (H300/H310/H330). Handle in fume hood with cyanide antidote kit available.",
)

_t(
    name="SN2 azide on primary alkyl halide",
    named_reaction="SN2 azide displacement",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["NaN3"],
    solvents=["DMF", "DMSO"],
    temperature_range=(25.0, 60.0),
    typical_yield=(75.0, 95.0),
    functional_group_required=["alkyl_halide"],
    functional_group_produced=["azide"],
    functional_group_incompatible=[],
    mechanism="Azide ion displaces halide in a concerted SN2 mechanism.",
    reverse_transform="Disconnect C-N3 bond; replace azide with halide.",
    safety_notes="NaN3 is highly toxic and forms explosive heavy-metal azides with Cu/Pb. Never acidify waste.",
)

_t(
    name="SN1 tertiary halide solvolysis",
    named_reaction="SN1",
    category=ReactionCategory.SUBSTITUTION,
    reagents=[],
    solvents=["water", "EtOH", "AcOH"],
    temperature_range=(25.0, 80.0),
    typical_yield=(50.0, 85.0),
    functional_group_required=["alkyl_halide"],
    functional_group_produced=["alcohol"],
    functional_group_incompatible=[],
    mechanism="Ionisation of the tertiary halide generates a carbocation, which is captured by the solvent nucleophile.",
    reverse_transform="Disconnect C-OH; reattach halide on the tertiary carbon.",
)

_t(
    name="Gabriel synthesis",
    named_reaction="Gabriel synthesis",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["potassium phthalimide", "N2H4"],
    solvents=["DMF", "EtOH"],
    temperature_range=(25.0, 100.0),
    typical_yield=(60.0, 85.0),
    functional_group_required=["alkyl_halide"],
    functional_group_produced=["primary_amine"],
    functional_group_incompatible=[],
    mechanism="Phthalimide anion displaces halide (SN2); subsequent hydrazinolysis liberates the primary amine.",
    reverse_transform="Disconnect C-NH2; replace amine with halide.",
)

_t(
    name="Finkelstein reaction",
    named_reaction="Finkelstein reaction",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["NaI"],
    solvents=["acetone"],
    temperature_range=(25.0, 60.0),
    typical_yield=(75.0, 95.0),
    functional_group_required=["alkyl_halide_cl", "alkyl_halide_br"],
    functional_group_produced=["alkyl_halide_i"],
    functional_group_incompatible=[],
    mechanism="Iodide displaces chloride or bromide in acetone; NaCl or NaBr precipitates, driving equilibrium.",
    reverse_transform="Replace iodide with the original halide.",
)

_t(
    name="Mitsunobu reaction",
    named_reaction="Mitsunobu reaction",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["PPh3", "DIAD", "nucleophile (phthalimide, carboxylic acid, etc.)"],
    solvents=["THF"],
    temperature_range=(0.0, 25.0),
    typical_yield=(60.0, 90.0),
    functional_group_required=["alcohol"],
    functional_group_produced=["ester", "ether", "azide"],
    functional_group_incompatible=[],
    mechanism="PPh3/DIAD activate the alcohol with inversion; a pronucleophile (pKa < 11) displaces via SN2.",
    reverse_transform="Disconnect the new C-nucleophile bond; restore the alcohol with inverted stereochemistry.",
    safety_notes="DIAD/DEAD are toxic and flammable. Triphenylphosphine oxide byproduct complicates purification.",
)

_t(
    name="Williamson ether synthesis",
    named_reaction="Williamson ether synthesis",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["NaH", "alkyl halide (primary)"],
    solvents=["THF", "DMF"],
    temperature_range=(0.0, 70.0),
    typical_yield=(70.0, 95.0),
    functional_group_required=["alcohol"],
    functional_group_produced=["ether"],
    functional_group_incompatible=["carboxylic_acid", "primary_amine"],
    mechanism="NaH deprotonates the alcohol to form an alkoxide, which then displaces a primary halide via SN2.",
    reverse_transform="Disconnect C-O-C; one fragment becomes alcohol, other becomes alkyl halide.",
    safety_notes="NaH is pyrophoric when dry. Generates H2 gas.",
)

_t(
    name="SN2 with KOtBu on primary halide",
    named_reaction="SN2 (bulky base can give elimination side products)",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["KOtBu"],
    solvents=["DMSO", "THF"],
    temperature_range=(25.0, 60.0),
    typical_yield=(40.0, 75.0),
    functional_group_required=["alkyl_halide"],
    functional_group_produced=["ether"],
    functional_group_incompatible=[],
    mechanism="tert-Butoxide can act as nucleophile on unhindered substrates, but elimination (E2) competes.",
    reverse_transform="Disconnect C-O-tBu; replace with halide.",
    safety_notes="KOtBu is pyrophoric and moisture-sensitive.",
)

# -----------------------------------------------------------------
#  ELIMINATION  (~6)
# -----------------------------------------------------------------

_t(
    name="E2 elimination with KOtBu",
    named_reaction="E2",
    category=ReactionCategory.ELIMINATION,
    reagents=["KOtBu"],
    solvents=["THF", "t-BuOH"],
    temperature_range=(25.0, 80.0),
    typical_yield=(65.0, 90.0),
    functional_group_required=["alkyl_halide"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=[],
    mechanism="Strong bulky base abstracts a beta-proton anti-periplanar to the leaving group in a concerted E2 mechanism.",
    reverse_transform="Add H-X across the double bond.",
    safety_notes="KOtBu is pyrophoric.",
)

_t(
    name="E1 acid-catalyzed dehydration of alcohol",
    named_reaction="E1 dehydration",
    category=ReactionCategory.ELIMINATION,
    reagents=["H2SO4"],
    solvents=["neat", "toluene"],
    temperature_range=(80.0, 180.0),
    typical_yield=(50.0, 85.0),
    functional_group_required=["alcohol"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=["epoxide"],
    mechanism="Protonation of the hydroxyl followed by loss of water generates a carbocation; loss of a beta-proton gives the alkene.",
    reverse_transform="Add water across the double bond (acid-catalysed hydration).",
    safety_notes="Concentrated H2SO4 is extremely corrosive.",
)

_t(
    name="Hofmann elimination",
    named_reaction="Hofmann elimination",
    category=ReactionCategory.ELIMINATION,
    reagents=["MeI (excess)", "Ag2O", "heat"],
    solvents=["water"],
    temperature_range=(100.0, 200.0),
    typical_yield=(40.0, 70.0),
    functional_group_required=["primary_amine", "secondary_amine", "tertiary_amine"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=[],
    mechanism="Exhaustive methylation of the amine followed by hydroxide-induced E2 elimination gives the less substituted (Hofmann) alkene.",
    reverse_transform="Amino-alkylation of the less-substituted end of the alkene.",
)

_t(
    name="Cope elimination",
    named_reaction="Cope elimination",
    category=ReactionCategory.ELIMINATION,
    reagents=["mCPBA (to form N-oxide)", "heat"],
    solvents=["DCM", "toluene"],
    temperature_range=(80.0, 150.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["tertiary_amine"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=[],
    mechanism="Oxidation of a tertiary amine to the N-oxide, which undergoes syn-periplanar intramolecular elimination via a five-membered cyclic transition state.",
    reverse_transform="Add amine across the double bond.",
)

_t(
    name="Dehydrohalogenation with NaOEt",
    named_reaction="E2 dehydrohalogenation",
    category=ReactionCategory.ELIMINATION,
    reagents=["NaOEt"],
    solvents=["EtOH"],
    temperature_range=(50.0, 80.0),
    typical_yield=(60.0, 85.0),
    functional_group_required=["alkyl_halide"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=["ester"],
    mechanism="Ethoxide base removes a beta-proton anti to the halide in a concerted E2 mechanism, favouring the Zaitsev product.",
    reverse_transform="Add HX across the alkene.",
)

_t(
    name="Zaitsev vs Hofmann selectivity (E2, size of base)",
    named_reaction="E2 selectivity",
    category=ReactionCategory.ELIMINATION,
    reagents=["KOtBu (Hofmann)", "NaOEt (Zaitsev)"],
    solvents=["THF", "EtOH"],
    temperature_range=(25.0, 80.0),
    typical_yield=(60.0, 85.0),
    functional_group_required=["alkyl_halide"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=[],
    mechanism="Bulky bases (KOtBu) give the less-substituted (Hofmann) alkene; smaller bases (NaOEt) give the more-substituted (Zaitsev) alkene.",
    reverse_transform="Add HX to the appropriate regiochemistry.",
)

# -----------------------------------------------------------------
#  ADDITION  (~10)
# -----------------------------------------------------------------

_t(
    name="HBr addition to alkene (Markovnikov)",
    named_reaction="Electrophilic addition",
    category=ReactionCategory.ADDITION,
    reagents=["HBr"],
    solvents=["DCM"],
    temperature_range=(-10.0, 25.0),
    typical_yield=(70.0, 90.0),
    functional_group_required=["alkene"],
    functional_group_produced=["alkyl_halide_br"],
    functional_group_incompatible=["alcohol"],
    mechanism="Protonation of the alkene forms the more stable carbocation (Markovnikov); bromide captures the cation.",
    reverse_transform="Remove HBr from the alkyl bromide via E2.",
)

_t(
    name="Anti-Markovnikov HBr addition (peroxide)",
    named_reaction="Radical addition (anti-Markovnikov)",
    category=ReactionCategory.ADDITION,
    reagents=["HBr", "ROOR (peroxide initiator)"],
    solvents=["CCl4", "pentane"],
    temperature_range=(25.0, 80.0),
    typical_yield=(60.0, 85.0),
    functional_group_required=["alkene"],
    functional_group_produced=["alkyl_halide_br"],
    functional_group_incompatible=[],
    mechanism="Peroxide generates Br radical; Br adds to the less substituted carbon giving the more stable radical, then H-atom transfer completes addition.",
    reverse_transform="Remove HBr from anti-Markovnikov position.",
    safety_notes="Peroxides are shock-sensitive. Use freshly prepared.",
)

_t(
    name="Hydroboration-oxidation",
    named_reaction="Hydroboration-oxidation",
    category=ReactionCategory.ADDITION,
    reagents=["BH3*THF (or 9-BBN)", "H2O2", "NaOH"],
    solvents=["THF"],
    temperature_range=(0.0, 25.0),
    typical_yield=(70.0, 90.0),
    functional_group_required=["alkene"],
    functional_group_produced=["alcohol"],
    functional_group_incompatible=["aldehyde", "ketone"],
    mechanism="Syn addition of B-H across the double bond with anti-Markovnikov regioselectivity; oxidative workup replaces B with OH with retention.",
    reverse_transform="Disconnect C-OH at the anti-Markovnikov position; restore alkene.",
    safety_notes="BH3 is pyrophoric; handle under inert atmosphere.",
)

_t(
    name="Epoxidation with mCPBA",
    named_reaction="Prilezhaev epoxidation",
    category=ReactionCategory.ADDITION,
    reagents=["mCPBA"],
    solvents=["DCM"],
    temperature_range=(0.0, 25.0),
    typical_yield=(70.0, 95.0),
    functional_group_required=["alkene"],
    functional_group_produced=["epoxide"],
    functional_group_incompatible=["primary_amine", "thiol"],
    mechanism="Concerted delivery of oxygen from the peracid to the alkene via a butterfly transition state; syn addition.",
    reverse_transform="Open the epoxide to reveal the original alkene.",
)

_t(
    name="Dihydroxylation with OsO4",
    named_reaction="Upjohn dihydroxylation",
    category=ReactionCategory.ADDITION,
    reagents=["OsO4 (cat.)", "NMO"],
    solvents=["acetone/water", "THF/water"],
    catalysts=["OsO4"],
    temperature_range=(0.0, 25.0),
    typical_yield=(75.0, 95.0),
    functional_group_required=["alkene"],
    functional_group_produced=["diol"],
    functional_group_incompatible=[],
    mechanism="OsO4 adds across the alkene in a [3+2] cycloaddition to form an osmate ester; NMO re-oxidises Os(VI) to Os(VIII) catalytically. Syn diol is produced.",
    reverse_transform="Convert 1,2-diol back to alkene (via periodate cleavage or Corey-Winter).",
    safety_notes="OsO4 is extremely toxic and volatile. Use catalytic amounts with NMO co-oxidant.",
)

_t(
    name="Catalytic hydrogenation of alkene",
    named_reaction="Catalytic hydrogenation",
    category=ReactionCategory.ADDITION,
    reagents=["H2"],
    solvents=["EtOAc", "MeOH", "EtOH"],
    catalysts=["Pd/C", "PtO2"],
    temperature_range=(20.0, 50.0),
    typical_yield=(85.0, 99.0),
    functional_group_required=["alkene"],
    functional_group_produced=[],
    functional_group_incompatible=["alkyl_halide_br", "azide"],
    mechanism="Syn addition of H2 across the double bond on the catalyst surface.",
    reverse_transform="Re-introduce unsaturation by dehydrogenation or elimination.",
    safety_notes="H2 is extremely flammable. Pd/C can be pyrophoric when dry.",
)

_t(
    name="Halogenation of alkene (Br2/CCl4)",
    named_reaction="Halogenation",
    category=ReactionCategory.ADDITION,
    reagents=["Br2"],
    solvents=["CCl4", "DCM"],
    temperature_range=(-10.0, 25.0),
    typical_yield=(75.0, 95.0),
    functional_group_required=["alkene"],
    functional_group_produced=["alkyl_halide_br"],
    functional_group_incompatible=[],
    mechanism="Bromine adds anti across the double bond via a bromonium ion intermediate. Both carbons receive Br (vicinal dibromide).",
    reverse_transform="Remove two Br atoms by Zn in AcOH to restore the alkene.",
    safety_notes="Br2 is highly corrosive and toxic. Work in fume hood.",
)

_t(
    name="Acid-catalysed hydration of alkene",
    named_reaction="Acid-catalysed hydration",
    category=ReactionCategory.ADDITION,
    reagents=["H2SO4 (cat.)", "H2O"],
    solvents=["water/THF"],
    temperature_range=(25.0, 80.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["alkene"],
    functional_group_produced=["alcohol"],
    functional_group_incompatible=["epoxide"],
    mechanism="Protonation of the alkene gives the Markovnikov carbocation; water captures it to yield the alcohol after deprotonation.",
    reverse_transform="Dehydrate the alcohol to the alkene.",
)

_t(
    name="Ozonolysis",
    named_reaction="Ozonolysis",
    category=ReactionCategory.ADDITION,
    reagents=["O3", "Me2S (or PPh3)"],
    solvents=["DCM", "MeOH"],
    temperature_range=(-78.0, -40.0),
    typical_yield=(65.0, 90.0),
    functional_group_required=["alkene"],
    functional_group_produced=["aldehyde", "ketone"],
    functional_group_incompatible=["thiol"],
    mechanism="Ozone adds across the double bond via a [3+2]/retro-[3+2]/[3+2] sequence forming an ozonide; reductive workup (Me2S or PPh3) cleaves it to carbonyls.",
    reverse_transform="Join two carbonyl fragments with a Wittig or McMurry coupling to reform the alkene.",
    safety_notes="Ozone is a powerful oxidizer and respiratory hazard. Ozonides are explosive; always quench completely before workup.",
)

_t(
    name="Simmons-Smith cyclopropanation",
    named_reaction="Simmons-Smith",
    category=ReactionCategory.ADDITION,
    reagents=["CH2I2", "Zn-Cu couple (or ZnEt2)"],
    solvents=["Et2O", "DCM"],
    temperature_range=(0.0, 40.0),
    typical_yield=(55.0, 85.0),
    functional_group_required=["alkene"],
    functional_group_produced=["cyclopropane"],
    functional_group_incompatible=[],
    mechanism="Zinc carbenoid (IZnCH2I) delivers CH2 to the less hindered face of the alkene in a concerted [2+1] cycloaddition.",
    reverse_transform="Open the cyclopropane to reveal alkene + CH2 unit.",
    safety_notes="CH2I2 is a lachrymator; ZnEt2 is pyrophoric.",
)

# -----------------------------------------------------------------
#  OXIDATION  (~8)
# -----------------------------------------------------------------

_t(
    name="PCC oxidation of primary alcohol to aldehyde",
    named_reaction="PCC oxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["PCC"],
    solvents=["DCM"],
    temperature_range=(20.0, 25.0),
    typical_yield=(70.0, 90.0),
    functional_group_required=["alcohol"],
    functional_group_produced=["aldehyde"],
    functional_group_incompatible=["primary_amine"],
    mechanism="PCC oxidises the primary alcohol to the aldehyde stage and stops (no over-oxidation to carboxylic acid) because the mild conditions and anhydrous solvent prevent further oxidation.",
    reverse_transform="Reduce aldehyde to primary alcohol (NaBH4).",
    safety_notes="PCC contains Cr(VI) -- carcinogenic. Use minimum amount and dispose of Cr waste properly.",
)

_t(
    name="Jones oxidation of alcohol to carboxylic acid / ketone",
    named_reaction="Jones oxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["CrO3", "H2SO4"],
    solvents=["acetone"],
    temperature_range=(0.0, 25.0),
    typical_yield=(70.0, 90.0),
    functional_group_required=["alcohol"],
    functional_group_produced=["carboxylic_acid", "ketone"],
    functional_group_incompatible=["alkene"],
    mechanism="Cr(VI) oxidises primary alcohols all the way to carboxylic acids and secondary alcohols to ketones.",
    reverse_transform="Reduce carboxylic acid/ketone back to alcohol.",
    safety_notes="Cr(VI) is carcinogenic. Generate minimum chromium waste.",
)

_t(
    name="Swern oxidation",
    named_reaction="Swern oxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["oxalyl chloride", "DMSO", "Et3N"],
    solvents=["DCM"],
    temperature_range=(-78.0, -40.0),
    typical_yield=(75.0, 95.0),
    functional_group_required=["alcohol"],
    functional_group_produced=["aldehyde", "ketone"],
    functional_group_incompatible=[],
    mechanism="DMSO is activated by oxalyl chloride; the activated species oxidises the alcohol via an alkoxysulfonium intermediate. Triethylamine triggers intramolecular elimination to yield the carbonyl.",
    reverse_transform="Reduce aldehyde/ketone back to alcohol.",
    safety_notes="Must keep at -78 C to avoid Pummerer side reactions. DMSO/oxalyl chloride generates CO, CO2, and (CH3)2S (foul odour).",
)

_t(
    name="KMnO4 oxidation of alkene",
    named_reaction="Permanganate oxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["KMnO4"],
    solvents=["water/acetone", "water"],
    temperature_range=(0.0, 100.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["alkene"],
    functional_group_produced=["diol", "carboxylic_acid", "ketone"],
    functional_group_incompatible=["aldehyde"],
    mechanism="Cold dilute KMnO4 gives syn-dihydroxylation; hot concentrated KMnO4 cleaves the double bond to carboxylic acids/ketones.",
    reverse_transform="Recombine fragments or reduce products.",
    safety_notes="Strong oxidizer (H272). Purple waste can be quenched with sodium bisulfite.",
)

_t(
    name="Dess-Martin periodinane oxidation",
    named_reaction="Dess-Martin oxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["Dess-Martin periodinane"],
    solvents=["DCM"],
    temperature_range=(20.0, 25.0),
    typical_yield=(80.0, 98.0),
    functional_group_required=["alcohol"],
    functional_group_produced=["aldehyde", "ketone"],
    functional_group_incompatible=[],
    mechanism="Hypervalent iodine(V) reagent oxidises alcohols to carbonyls under very mild, neutral conditions with high functional-group tolerance.",
    reverse_transform="Reduce carbonyl back to alcohol.",
    safety_notes="Shock-sensitive when dry. Keep as solution.",
)

_t(
    name="Baeyer-Villiger oxidation",
    named_reaction="Baeyer-Villiger oxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["mCPBA"],
    solvents=["DCM"],
    temperature_range=(0.0, 40.0),
    typical_yield=(60.0, 90.0),
    functional_group_required=["ketone"],
    functional_group_produced=["ester"],
    functional_group_incompatible=["primary_amine"],
    mechanism="Peracid attacks the ketone to form a Criegee intermediate; 1,2-alkyl shift inserts an oxygen between the carbonyl carbon and the migrating group.",
    reverse_transform="Disconnect the ester C-O bond; recombine as ketone.",
)

_t(
    name="Sharpless asymmetric epoxidation",
    named_reaction="Sharpless epoxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["Ti(OiPr)4", "tBuOOH", "L-(+)-DIPT or D-(-)-DIPT"],
    solvents=["DCM"],
    catalysts=["Ti(OiPr)4"],
    temperature_range=(-20.0, 0.0),
    typical_yield=(70.0, 90.0),
    functional_group_required=["alkene", "alcohol"],
    functional_group_produced=["epoxide"],
    functional_group_incompatible=["thiol"],
    mechanism="Titanium-tartrate complex directs epoxidation of allylic alcohols with high enantioselectivity using tBuOOH as terminal oxidant.",
    reverse_transform="Open epoxide to regenerate allylic alcohol.",
    safety_notes="tBuOOH is an organic peroxide -- shock sensitive above certain concentrations.",
)

_t(
    name="Wacker oxidation",
    named_reaction="Wacker oxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["PdCl2 (cat.)", "CuCl2 (co-cat.)", "O2"],
    solvents=["DMF/water"],
    catalysts=["PdCl2", "CuCl2"],
    temperature_range=(25.0, 80.0),
    typical_yield=(60.0, 85.0),
    functional_group_required=["alkene"],
    functional_group_produced=["ketone"],
    functional_group_incompatible=[],
    mechanism="Pd(II)-catalysed nucleophilic attack of water on the coordinated alkene (Markovnikov selectivity) gives a methyl ketone. CuCl2/O2 re-oxidises Pd(0) to Pd(II).",
    reverse_transform="Disconnect CH3-C(=O)-R; rebuild terminal alkene.",
)

# -----------------------------------------------------------------
#  REDUCTION  (~8)
# -----------------------------------------------------------------

_t(
    name="NaBH4 reduction of aldehyde/ketone to alcohol",
    named_reaction="NaBH4 reduction",
    category=ReactionCategory.REDUCTION,
    reagents=["NaBH4"],
    solvents=["MeOH", "EtOH"],
    temperature_range=(0.0, 25.0),
    typical_yield=(80.0, 98.0),
    functional_group_required=["aldehyde", "ketone"],
    functional_group_produced=["alcohol"],
    functional_group_incompatible=[],
    mechanism="Hydride delivery from borohydride to the electrophilic carbonyl carbon gives an alkoxide, protonated on workup.",
    reverse_transform="Oxidise alcohol to aldehyde/ketone.",
)

_t(
    name="LiAlH4 reduction of ester to primary alcohol",
    named_reaction="LiAlH4 reduction",
    category=ReactionCategory.REDUCTION,
    reagents=["LiAlH4"],
    solvents=["THF", "Et2O"],
    temperature_range=(0.0, 70.0),
    typical_yield=(75.0, 95.0),
    functional_group_required=["ester"],
    functional_group_produced=["alcohol"],
    functional_group_incompatible=["alkyl_halide", "nitrile"],
    mechanism="Two equivalents of hydride reduce the ester first to an aldehyde (hemiacetal collapses) then to the primary alcohol.",
    reverse_transform="Oxidise alcohol and re-esterify.",
    safety_notes="LiAlH4 reacts violently with water and protic solvents. Quench with careful Fieser workup (H2O, NaOH, H2O).",
)

_t(
    name="LiAlH4 reduction of amide to amine",
    named_reaction="LiAlH4 reduction",
    category=ReactionCategory.REDUCTION,
    reagents=["LiAlH4"],
    solvents=["THF"],
    temperature_range=(0.0, 70.0),
    typical_yield=(65.0, 90.0),
    functional_group_required=["amide"],
    functional_group_produced=["primary_amine", "secondary_amine", "tertiary_amine"],
    functional_group_incompatible=["alkyl_halide"],
    mechanism="LiAlH4 reduces the amide C=O completely, yielding the corresponding amine.",
    reverse_transform="Oxidise amine back to amide (difficult) or form amide from amine + acyl chloride.",
    safety_notes="LiAlH4 is pyrophoric. Use strict anhydrous conditions.",
)

_t(
    name="DIBAL-H reduction of ester to aldehyde at -78 C",
    named_reaction="DIBAL-H partial reduction",
    category=ReactionCategory.REDUCTION,
    reagents=["DIBAL-H"],
    solvents=["DCM", "toluene"],
    temperature_range=(-78.0, -60.0),
    typical_yield=(60.0, 85.0),
    functional_group_required=["ester"],
    functional_group_produced=["aldehyde"],
    functional_group_incompatible=[],
    mechanism="At -78 C, DIBAL-H delivers one equivalent of hydride to the ester, stopping at the aluminium hemiacetal stage. Aqueous workup liberates the aldehyde.",
    reverse_transform="Re-esterify the aldehyde fragment.",
    safety_notes="DIBAL-H is pyrophoric. Strict -78 C temperature control is essential to prevent over-reduction.",
)

_t(
    name="Catalytic hydrogenation (general)",
    named_reaction="Catalytic hydrogenation",
    category=ReactionCategory.REDUCTION,
    reagents=["H2"],
    solvents=["EtOAc", "MeOH", "EtOH"],
    catalysts=["Pd/C", "PtO2", "Rh/C"],
    temperature_range=(20.0, 50.0),
    typical_yield=(85.0, 99.0),
    functional_group_required=["alkene", "alkyne", "nitrile", "nitro"],
    functional_group_produced=["alcohol", "primary_amine"],
    functional_group_incompatible=[],
    mechanism="Substrate adsorbs on metal surface; H2 dissociates and adds syn to the unsaturated bond.",
    reverse_transform="Re-oxidise or dehydrogenate.",
    safety_notes="H2 is highly flammable (H220). Pd/C can be pyrophoric when dry.",
)

_t(
    name="Birch reduction",
    named_reaction="Birch reduction",
    category=ReactionCategory.REDUCTION,
    reagents=["Na (or Li)", "NH3 (liquid)", "t-BuOH"],
    solvents=["liquid NH3"],
    temperature_range=(-78.0, -33.0),
    typical_yield=(60.0, 85.0),
    functional_group_required=["aromatic_ring"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=["alkyl_halide_br"],
    mechanism="Dissolving metal (Na/Li in liq. NH3) delivers electrons to the aromatic ring; proton source (t-BuOH) quenches the radical anion to give 1,4-cyclohexadienes.",
    reverse_transform="Re-aromatise by oxidation (DDQ or Pd/C).",
    safety_notes="Liquid NH3 is extremely hazardous (bp -33 C). Na metal is pyrophoric. Requires cold-finger condenser.",
)

_t(
    name="Wolff-Kishner reduction",
    named_reaction="Wolff-Kishner reduction",
    category=ReactionCategory.REDUCTION,
    reagents=["N2H4", "KOH"],
    solvents=["diethylene glycol", "ethylene glycol"],
    temperature_range=(180.0, 230.0),
    typical_yield=(60.0, 85.0),
    functional_group_required=["aldehyde", "ketone"],
    functional_group_produced=[],
    functional_group_incompatible=["ester", "alkyl_halide"],
    mechanism="Hydrazone formation followed by base-mediated decomposition with loss of N2 converts C=O to CH2.",
    reverse_transform="Oxidise the methylene back to carbonyl.",
    safety_notes="Hydrazine is toxic and carcinogenic (H350). Very high temperatures required.",
)

_t(
    name="Clemmensen reduction",
    named_reaction="Clemmensen reduction",
    category=ReactionCategory.REDUCTION,
    reagents=["Zn(Hg)", "HCl (conc.)"],
    solvents=["water/HCl"],
    temperature_range=(80.0, 110.0),
    typical_yield=(55.0, 80.0),
    functional_group_required=["aldehyde", "ketone"],
    functional_group_produced=[],
    functional_group_incompatible=["alcohol", "primary_amine"],
    mechanism="Zinc amalgam in concentrated HCl reduces the carbonyl to a methylene group under acidic conditions (complement to Wolff-Kishner which uses base).",
    reverse_transform="Oxidise the methylene back to carbonyl.",
    safety_notes="Uses mercury amalgam (environmental hazard) and concentrated HCl.",
)

_t(
    name="Luche reduction (NaBH4/CeCl3)",
    named_reaction="Luche reduction",
    category=ReactionCategory.REDUCTION,
    reagents=["NaBH4", "CeCl3*7H2O"],
    solvents=["MeOH"],
    temperature_range=(-10.0, 0.0),
    typical_yield=(75.0, 95.0),
    functional_group_required=["ketone"],
    functional_group_produced=["alcohol"],
    functional_group_incompatible=[],
    mechanism="CeCl3 activates the carbonyl as a hard electrophile, directing NaBH4 to achieve selective 1,2-reduction of enones (leaving conjugated alkene intact).",
    reverse_transform="Oxidise allylic alcohol back to enone.",
)

# -----------------------------------------------------------------
#  COUPLING  (~10)
# -----------------------------------------------------------------

_t(
    name="Grignard addition to aldehyde",
    named_reaction="Grignard reaction",
    category=ReactionCategory.COUPLING,
    reagents=["RMgBr"],
    solvents=["THF", "Et2O"],
    temperature_range=(-20.0, 25.0),
    typical_yield=(70.0, 95.0),
    functional_group_required=["aldehyde"],
    functional_group_produced=["alcohol"],
    functional_group_incompatible=["alcohol", "carboxylic_acid", "primary_amine"],
    mechanism="Grignard reagent adds to the aldehyde carbonyl; aqueous acid workup yields a secondary alcohol.",
    reverse_transform="Disconnect C-C bond alpha to the OH; one fragment becomes RMgBr, the other becomes aldehyde.",
    safety_notes="Grignard reagents are pyrophoric in air. Strict anhydrous conditions.",
)

_t(
    name="Grignard addition to ketone",
    named_reaction="Grignard reaction",
    category=ReactionCategory.COUPLING,
    reagents=["RMgBr"],
    solvents=["THF", "Et2O"],
    temperature_range=(-20.0, 25.0),
    typical_yield=(65.0, 90.0),
    functional_group_required=["ketone"],
    functional_group_produced=["alcohol"],
    functional_group_incompatible=["alcohol", "carboxylic_acid", "primary_amine"],
    mechanism="Grignard reagent adds to the ketone carbonyl; aqueous acid workup gives a tertiary alcohol.",
    reverse_transform="Disconnect one of the three C-C bonds at the tertiary alcohol carbon.",
)

_t(
    name="Grignard addition to ester (gives tertiary alcohol)",
    named_reaction="Grignard reaction (2 equiv on ester)",
    category=ReactionCategory.COUPLING,
    reagents=["RMgBr (2 equiv)"],
    solvents=["THF"],
    temperature_range=(-20.0, 25.0),
    typical_yield=(60.0, 85.0),
    functional_group_required=["ester"],
    functional_group_produced=["alcohol"],
    functional_group_incompatible=["alcohol", "primary_amine"],
    mechanism="First equivalent of RMgBr attacks the ester; the tetrahedral intermediate collapses to a ketone + alkoxide. Second equivalent adds to the ketone to give a tertiary alcohol.",
    reverse_transform="Disconnect two identical R groups from the tertiary alcohol; one fragment is RMgBr, the other is the ester.",
)

_t(
    name="Grignard with CO2 (carboxylation)",
    named_reaction="Grignard carboxylation",
    category=ReactionCategory.COUPLING,
    reagents=["RMgBr", "CO2 (dry ice)"],
    solvents=["THF", "Et2O"],
    temperature_range=(-78.0, 0.0),
    typical_yield=(65.0, 90.0),
    functional_group_required=[],
    functional_group_produced=["carboxylic_acid"],
    functional_group_incompatible=["alcohol", "primary_amine", "carboxylic_acid"],
    mechanism="Grignard nucleophile attacks CO2 to form a carboxylate; acidic workup gives the carboxylic acid with one extra carbon.",
    reverse_transform="Disconnect the C-COOH bond; carbon fragment becomes RMgBr.",
)

_t(
    name="Suzuki coupling",
    named_reaction="Suzuki-Miyaura coupling",
    category=ReactionCategory.COUPLING,
    reagents=["ArB(OH)2 (boronic acid)", "base (K2CO3, Cs2CO3, K3PO4)"],
    solvents=["THF/water", "dioxane/water", "DMF"],
    catalysts=["Pd(PPh3)4", "PdCl2(dppf)", "Pd(OAc)2/SPhos"],
    temperature_range=(60.0, 100.0),
    typical_yield=(70.0, 95.0),
    functional_group_required=["alkyl_halide_br", "alkyl_halide_i"],
    functional_group_produced=[],
    functional_group_incompatible=[],
    mechanism="Pd(0) oxidatively adds to the aryl/vinyl halide; transmetalation with the boronate ester transfers the second partner; reductive elimination forms the new C-C bond.",
    reverse_transform="Disconnect the biaryl C-C bond; one fragment gets halide, the other gets B(OH)2.",
    safety_notes="Pd catalysts are expensive. Boronic acids are generally low toxicity.",
)

_t(
    name="Heck reaction",
    named_reaction="Heck reaction",
    category=ReactionCategory.COUPLING,
    reagents=["aryl/vinyl halide", "alkene", "Et3N (or DIPEA)"],
    solvents=["DMF", "MeCN", "toluene"],
    catalysts=["Pd(OAc)2", "Pd(PPh3)4"],
    temperature_range=(80.0, 140.0),
    typical_yield=(60.0, 90.0),
    functional_group_required=["alkene", "alkyl_halide_br"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=[],
    mechanism="Pd(0) oxidatively adds to the halide; syn migratory insertion into the alkene followed by beta-hydride elimination gives a new substituted alkene with retention of the double bond.",
    reverse_transform="Disconnect the vinylic C-C bond; one fragment gets halide.",
)

_t(
    name="Sonogashira coupling",
    named_reaction="Sonogashira coupling",
    category=ReactionCategory.COUPLING,
    reagents=["terminal alkyne", "aryl/vinyl halide", "Et3N (or DIPEA)"],
    solvents=["THF", "DMF", "Et3N (as solvent)"],
    catalysts=["Pd(PPh3)4 (or PdCl2(PPh3)2)", "CuI"],
    temperature_range=(25.0, 80.0),
    typical_yield=(65.0, 90.0),
    functional_group_required=["alkyne", "alkyl_halide_br"],
    functional_group_produced=["alkyne"],
    functional_group_incompatible=[],
    mechanism="Cu(I) generates a copper acetylide from the terminal alkyne; Pd(0) oxidatively adds to the halide; transmetalation and reductive elimination form the new C(sp)-C(sp2) bond.",
    reverse_transform="Disconnect the C(sp)-C(sp2) bond; one fragment gets halide, the other is a terminal alkyne.",
)

_t(
    name="Wittig reaction",
    named_reaction="Wittig reaction",
    category=ReactionCategory.COUPLING,
    reagents=["Ph3P=CHR (Wittig ylide)", "n-BuLi (for ylide generation)"],
    solvents=["THF"],
    temperature_range=(-78.0, 25.0),
    typical_yield=(55.0, 85.0),
    functional_group_required=["aldehyde", "ketone"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=[],
    mechanism="Phosphorus ylide attacks the carbonyl; a betaine forms and closes to an oxaphosphetane, which undergoes retro-[2+2] to give the alkene and Ph3P=O.",
    reverse_transform="Disconnect the C=C bond; one fragment becomes an aldehyde/ketone, the other becomes a phosphonium salt.",
    safety_notes="n-BuLi is pyrophoric. Ph3P=O is a troublesome byproduct (hard to remove).",
)

_t(
    name="Horner-Wadsworth-Emmons reaction",
    named_reaction="Horner-Wadsworth-Emmons",
    category=ReactionCategory.COUPLING,
    reagents=["(EtO)2P(O)CHR'CO2Et", "NaH (or LDA)"],
    solvents=["THF"],
    temperature_range=(-78.0, 25.0),
    typical_yield=(60.0, 90.0),
    functional_group_required=["aldehyde", "ketone"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=[],
    mechanism="Stabilised phosphonate carbanion attacks the carbonyl; the resulting betaine cyclises and fragments (retro-[2+2]) to give predominantly the E-alkene and water-soluble diethylphosphate.",
    reverse_transform="Disconnect the C=C; one fragment becomes aldehyde/ketone, the other a phosphonate ester.",
)

_t(
    name="Aldol reaction",
    named_reaction="Aldol reaction",
    category=ReactionCategory.COUPLING,
    reagents=["LDA (for enolate formation)"],
    solvents=["THF"],
    temperature_range=(-78.0, 0.0),
    typical_yield=(60.0, 85.0),
    functional_group_required=["aldehyde", "ketone"],
    functional_group_produced=["alcohol", "ketone"],
    functional_group_incompatible=[],
    mechanism="Base generates an enolate from one carbonyl; the enolate attacks a second carbonyl to form a beta-hydroxy carbonyl (aldol product). Heating can drive dehydration to the enone.",
    reverse_transform="Disconnect the C-C bond between the alpha and beta carbons of the beta-hydroxy carbonyl; both fragments are carbonyls.",
)

_t(
    name="Claisen condensation",
    named_reaction="Claisen condensation",
    category=ReactionCategory.COUPLING,
    reagents=["NaOEt"],
    solvents=["EtOH"],
    temperature_range=(25.0, 80.0),
    typical_yield=(50.0, 75.0),
    functional_group_required=["ester"],
    functional_group_produced=["ester", "ketone"],
    functional_group_incompatible=[],
    mechanism="Ethoxide deprotonates the ester alpha-position; the resulting enolate attacks a second ester molecule. Loss of ethoxide gives a beta-keto ester.",
    reverse_transform="Disconnect the C-C bond alpha to the ketone; both fragments become esters.",
)

_t(
    name="Stille coupling",
    named_reaction="Stille coupling",
    category=ReactionCategory.COUPLING,
    reagents=["R-SnBu3 (organostannane)", "aryl/vinyl halide"],
    solvents=["DMF", "toluene", "THF"],
    catalysts=["Pd(PPh3)4", "PdCl2(PPh3)2"],
    temperature_range=(60.0, 110.0),
    typical_yield=(60.0, 90.0),
    functional_group_required=["alkyl_halide_br", "alkyl_halide_i"],
    functional_group_produced=[],
    functional_group_incompatible=[],
    mechanism="Pd(0) oxidatively adds to the halide; transmetalation with the stannane transfers R; reductive elimination forms the new C-C bond.",
    reverse_transform="Disconnect the newly formed C-C bond; one fragment gets halide, the other gets SnBu3.",
    safety_notes="Organotin compounds are highly toxic (H301/H311/H331, H360, H372). Use with extreme care and proper waste disposal.",
)

_t(
    name="Negishi coupling",
    named_reaction="Negishi coupling",
    category=ReactionCategory.COUPLING,
    reagents=["RZnX (organozinc)", "aryl/vinyl halide"],
    solvents=["THF"],
    catalysts=["Pd(PPh3)4", "PdCl2(dppf)"],
    temperature_range=(0.0, 60.0),
    typical_yield=(65.0, 90.0),
    functional_group_required=["alkyl_halide_br", "alkyl_halide_i"],
    functional_group_produced=[],
    functional_group_incompatible=[],
    mechanism="Pd(0) oxidatively adds to the halide; transmetalation with the organozinc transfers R; reductive elimination forms the C-C bond. High functional group tolerance.",
    reverse_transform="Disconnect the C-C bond; one fragment gets halide, the other gets ZnX.",
)

# -----------------------------------------------------------------
#  CARBONYL CHEMISTRY  (~10)
# -----------------------------------------------------------------

_t(
    name="Fischer esterification",
    named_reaction="Fischer esterification",
    category=ReactionCategory.CARBONYL,
    reagents=["alcohol", "H2SO4 (cat.)"],
    solvents=["neat", "toluene (Dean-Stark)"],
    temperature_range=(60.0, 120.0),
    typical_yield=(50.0, 85.0),
    functional_group_required=["carboxylic_acid", "alcohol"],
    functional_group_produced=["ester"],
    functional_group_incompatible=["primary_amine"],
    mechanism="Acid-catalysed nucleophilic addition-elimination: protonation of C=O, attack by alcohol, loss of water. Equilibrium driven by excess alcohol or water removal.",
    reverse_transform="Disconnect the ester bond; fragments are carboxylic acid + alcohol.",
)

_t(
    name="Amide formation from acyl chloride and amine",
    named_reaction="Schotten-Baumann acylation",
    category=ReactionCategory.CARBONYL,
    reagents=["RCOCl (acyl chloride)", "Et3N (base)"],
    solvents=["DCM", "THF"],
    temperature_range=(0.0, 25.0),
    typical_yield=(75.0, 95.0),
    functional_group_required=["primary_amine", "secondary_amine"],
    functional_group_produced=["amide"],
    functional_group_incompatible=["alcohol"],
    mechanism="Nucleophilic amine attacks the acyl chloride; loss of HCl (scavenged by Et3N) gives the amide.",
    reverse_transform="Disconnect the amide C-N bond; one fragment becomes acyl chloride, the other the amine.",
)

_t(
    name="Claisen rearrangement",
    named_reaction="Claisen rearrangement",
    category=ReactionCategory.REARRANGEMENT,
    reagents=[],
    solvents=["toluene", "decalin", "neat"],
    temperature_range=(150.0, 250.0),
    typical_yield=(60.0, 90.0),
    functional_group_required=["ether", "alkene"],
    functional_group_produced=["aldehyde", "ketone"],
    functional_group_incompatible=[],
    mechanism="Allyl vinyl ether undergoes a concerted [3,3]-sigmatropic rearrangement through a chair-like transition state to give a gamma,delta-unsaturated carbonyl.",
    reverse_transform="Reconnect the sigma bond as the allyl vinyl ether.",
)

_t(
    name="Diels-Alder reaction",
    named_reaction="Diels-Alder",
    category=ReactionCategory.PERICYCLIC,
    reagents=["diene", "dienophile"],
    solvents=["toluene", "DCM", "neat"],
    temperature_range=(25.0, 200.0),
    typical_yield=(60.0, 95.0),
    functional_group_required=["alkene"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=[],
    mechanism="Concerted [4+2] cycloaddition between a conjugated diene (s-cis conformation) and a dienophile. Suprafacial on both components; endo rule favours kinetic product.",
    reverse_transform="Retro-Diels-Alder: disconnect the cyclohexene ring into diene + dienophile.",
)

_t(
    name="Beckmann rearrangement",
    named_reaction="Beckmann rearrangement",
    category=ReactionCategory.REARRANGEMENT,
    reagents=["hydroxylamine (to form oxime)", "H2SO4 (or PCl5, SOCl2)"],
    solvents=["AcOH", "DCM"],
    temperature_range=(25.0, 130.0),
    typical_yield=(55.0, 85.0),
    functional_group_required=["ketone"],
    functional_group_produced=["amide"],
    functional_group_incompatible=["primary_amine"],
    mechanism="Ketoxime undergoes acid-catalysed 1,2-shift of the group anti to the departing hydroxyl; nitrilium ion is captured by water to give the amide (lactam for cyclic ketones).",
    reverse_transform="Disconnect the amide C-N bond; recombine as ketone.",
)

_t(
    name="Cannizzaro reaction",
    named_reaction="Cannizzaro reaction",
    category=ReactionCategory.CARBONYL,
    reagents=["NaOH (conc.)"],
    solvents=["water"],
    temperature_range=(25.0, 80.0),
    typical_yield=(40.0, 50.0),
    functional_group_required=["aldehyde"],
    functional_group_produced=["alcohol", "carboxylic_acid"],
    functional_group_incompatible=[],
    mechanism="In concentrated base, an aldehyde without alpha-hydrogens disproportionates: one molecule is oxidised to the carboxylate, the other is reduced to the alcohol.",
    reverse_transform="Combine a primary alcohol and a carboxylic acid derived from the same aldehyde.",
)

_t(
    name="Tischenko reaction",
    named_reaction="Tischenko reaction",
    category=ReactionCategory.CARBONYL,
    reagents=["Al(OEt)3 (cat.)"],
    solvents=["neat", "toluene"],
    catalysts=["Al(OEt)3"],
    temperature_range=(0.0, 80.0),
    typical_yield=(60.0, 90.0),
    functional_group_required=["aldehyde"],
    functional_group_produced=["ester"],
    functional_group_incompatible=[],
    mechanism="Lewis acid-catalysed disproportionation of two aldehyde molecules: one is oxidised and one is reduced, joined as an ester. Related to Cannizzaro but gives the ester directly.",
    reverse_transform="Hydrolyse the ester; both fragments derive from the same aldehyde.",
)

_t(
    name="Michael addition",
    named_reaction="Michael addition",
    category=ReactionCategory.CARBONYL,
    reagents=["NaOEt (or LDA)"],
    solvents=["EtOH", "THF"],
    temperature_range=(-78.0, 25.0),
    typical_yield=(60.0, 90.0),
    functional_group_required=["alkene", "ketone"],
    functional_group_produced=["ketone"],
    functional_group_incompatible=[],
    mechanism="A stabilised nucleophile (enolate, malonate, etc.) adds in a 1,4-fashion to an alpha,beta-unsaturated carbonyl (Michael acceptor).",
    reverse_transform="Disconnect the C-C bond beta to the carbonyl; one fragment is the enolate donor, the other the enone.",
)

_t(
    name="Robinson annulation",
    named_reaction="Robinson annulation",
    category=ReactionCategory.CARBONYL,
    reagents=["NaOEt", "methyl vinyl ketone"],
    solvents=["EtOH", "THF"],
    temperature_range=(25.0, 80.0),
    typical_yield=(45.0, 75.0),
    functional_group_required=["ketone"],
    functional_group_produced=["ketone", "alkene"],
    functional_group_incompatible=[],
    mechanism="Michael addition of a ketone enolate to methyl vinyl ketone, followed by intramolecular aldol condensation and dehydration to form a cyclohexenone ring.",
    reverse_transform="Disconnect the cyclohexenone ring into a 1,5-diketone, then further into the ketone + MVK.",
)


# Baeyer-Villiger duplicate removed: the canonical entry under OXIDATION
# (line ~514) is the correct categorization.

# -----------------------------------------------------------------
#  PROTECTION  (~8)
# -----------------------------------------------------------------

_t(
    name="Boc protection of amine",
    named_reaction="Boc protection",
    category=ReactionCategory.PROTECTION,
    reagents=["Boc2O", "Et3N (or NaOH)"],
    solvents=["DCM", "THF", "water/dioxane"],
    temperature_range=(0.0, 25.0),
    typical_yield=(85.0, 98.0),
    functional_group_required=["primary_amine", "secondary_amine"],
    functional_group_produced=["carbamate"],
    functional_group_incompatible=[],
    mechanism="Nucleophilic amine attacks the electrophilic carbonyl of Boc2O; loss of tert-butoxide/CO2 gives the Boc-carbamate.",
    reverse_transform="Remove Boc with TFA or HCl in dioxane.",
)

_t(
    name="Boc deprotection",
    named_reaction="Boc deprotection",
    category=ReactionCategory.DEPROTECTION,
    reagents=["TFA (or HCl in dioxane)"],
    solvents=["DCM"],
    temperature_range=(0.0, 25.0),
    typical_yield=(90.0, 99.0),
    functional_group_required=["carbamate"],
    functional_group_produced=["primary_amine", "secondary_amine"],
    functional_group_incompatible=[],
    mechanism="Acid protonates the carbamate oxygen; tert-butyl cation departs, generating CO2 and the free amine.",
    reverse_transform="Protect amine with Boc2O.",
)

_t(
    name="Fmoc protection of amine",
    named_reaction="Fmoc protection",
    category=ReactionCategory.PROTECTION,
    reagents=["FmocCl (or Fmoc-OSu)", "Na2CO3 (or NaHCO3)"],
    solvents=["dioxane/water", "THF"],
    temperature_range=(0.0, 25.0),
    typical_yield=(80.0, 95.0),
    functional_group_required=["primary_amine", "secondary_amine"],
    functional_group_produced=["carbamate"],
    functional_group_incompatible=[],
    mechanism="Amine attacks FmocCl, displacing chloride; base scavenges HCl.",
    reverse_transform="Remove Fmoc with piperidine (20% in DMF).",
)

_t(
    name="Fmoc deprotection",
    named_reaction="Fmoc deprotection",
    category=ReactionCategory.DEPROTECTION,
    reagents=["piperidine (20% in DMF)"],
    solvents=["DMF"],
    temperature_range=(20.0, 25.0),
    typical_yield=(90.0, 99.0),
    functional_group_required=["carbamate"],
    functional_group_produced=["primary_amine", "secondary_amine"],
    functional_group_incompatible=[],
    mechanism="Piperidine abstracts the acidic fluorenyl proton; beta-elimination releases dibenzofulvene and CO2, liberating the free amine. Piperidine scavenges the dibenzofulvene.",
    reverse_transform="Protect amine with FmocCl.",
)

_t(
    name="TBS protection of alcohol",
    named_reaction="TBS silylation",
    category=ReactionCategory.PROTECTION,
    reagents=["TBSCl", "imidazole"],
    solvents=["DMF", "DCM"],
    temperature_range=(0.0, 25.0),
    typical_yield=(85.0, 98.0),
    functional_group_required=["alcohol"],
    functional_group_produced=["silyl_ether"],
    functional_group_incompatible=[],
    mechanism="Imidazole acts as nucleophilic catalyst, activating TBSCl; the alcohol then displaces imidazole to form the TBS ether.",
    reverse_transform="Remove TBS with TBAF in THF.",
)

_t(
    name="TBS deprotection",
    named_reaction="TBS removal",
    category=ReactionCategory.DEPROTECTION,
    reagents=["TBAF"],
    solvents=["THF"],
    temperature_range=(0.0, 25.0),
    typical_yield=(90.0, 99.0),
    functional_group_required=["silyl_ether"],
    functional_group_produced=["alcohol"],
    functional_group_incompatible=[],
    mechanism="Fluoride ion attacks silicon (very strong Si-F bond, 140 kcal/mol), displacing the alkoxide which is protonated on workup.",
    reverse_transform="Protect alcohol with TBSCl/imidazole.",
)

_t(
    name="Benzyl protection of alcohol",
    named_reaction="Benzyl protection",
    category=ReactionCategory.PROTECTION,
    reagents=["BnBr", "NaH"],
    solvents=["DMF", "THF"],
    temperature_range=(0.0, 25.0),
    typical_yield=(80.0, 95.0),
    functional_group_required=["alcohol"],
    functional_group_produced=["ether"],
    functional_group_incompatible=["alkyl_halide"],
    mechanism="NaH deprotonates the alcohol; the alkoxide displaces bromide from BnBr in an SN2 reaction.",
    reverse_transform="Remove Bn by hydrogenolysis (H2/Pd-C).",
)

_t(
    name="Benzyl deprotection (hydrogenolysis)",
    named_reaction="Hydrogenolysis",
    category=ReactionCategory.DEPROTECTION,
    reagents=["H2", "Pd/C"],
    solvents=["EtOAc", "MeOH", "EtOH"],
    catalysts=["Pd/C"],
    temperature_range=(20.0, 40.0),
    typical_yield=(85.0, 99.0),
    functional_group_required=["ether"],
    functional_group_produced=["alcohol"],
    functional_group_incompatible=["alkene", "alkyne"],
    mechanism="Pd-catalysed hydrogenolysis cleaves the benzylic C-O bond. Toluene is released as byproduct.",
    reverse_transform="Protect alcohol with BnBr/NaH.",
    safety_notes="H2/Pd-C: flammable gas + pyrophoric catalyst. Filter Pd/C through Celite before discarding.",
)

_t(
    name="Acetyl protection of alcohol",
    named_reaction="Acetylation",
    category=ReactionCategory.PROTECTION,
    reagents=["Ac2O", "Et3N (or pyridine)", "DMAP (cat.)"],
    solvents=["DCM"],
    catalysts=["DMAP"],
    temperature_range=(0.0, 25.0),
    typical_yield=(85.0, 98.0),
    functional_group_required=["alcohol"],
    functional_group_produced=["ester"],
    functional_group_incompatible=["primary_amine"],
    mechanism="DMAP-catalysed nucleophilic acyl substitution: DMAP activates acetic anhydride, the alcohol attacks, and acetic acid departs.",
    reverse_transform="Remove acetyl with K2CO3/MeOH or LiOH/THF-water.",
)

_t(
    name="PMB protection of alcohol",
    named_reaction="PMB protection",
    category=ReactionCategory.PROTECTION,
    reagents=["PMBCl", "NaH"],
    solvents=["DMF", "THF"],
    temperature_range=(0.0, 25.0),
    typical_yield=(75.0, 90.0),
    functional_group_required=["alcohol"],
    functional_group_produced=["ether"],
    functional_group_incompatible=[],
    mechanism="NaH deprotonates alcohol; the alkoxide displaces chloride from PMBCl in an SN2 reaction.",
    reverse_transform="Remove PMB with DDQ (oxidative) or TFA/anisole (acidic).",
)

_t(
    name="PMB deprotection with DDQ",
    named_reaction="PMB deprotection",
    category=ReactionCategory.DEPROTECTION,
    reagents=["DDQ"],
    solvents=["DCM/water", "DCM/pH 7 buffer"],
    temperature_range=(0.0, 25.0),
    typical_yield=(75.0, 95.0),
    functional_group_required=["ether"],
    functional_group_produced=["alcohol"],
    functional_group_incompatible=["alkene"],
    mechanism="DDQ oxidises the electron-rich PMB group to a quinone methide/oxocarbenium species; water (or buffer) traps it, releasing the free alcohol.",
    reverse_transform="Protect with PMBCl/NaH.",
)

_t(
    name="Acetonide protection of 1,2-diol",
    named_reaction="Acetonide protection",
    category=ReactionCategory.PROTECTION,
    reagents=["2,2-dimethoxypropane", "p-TsOH (cat.)"],
    solvents=["acetone", "DMF"],
    catalysts=["p-TsOH"],
    temperature_range=(20.0, 50.0),
    typical_yield=(80.0, 95.0),
    functional_group_required=["diol"],
    functional_group_produced=["acetal"],
    functional_group_incompatible=[],
    mechanism="Acid-catalysed transacetalisation: the diol displaces methanol from 2,2-dimethoxypropane to form a cyclic acetonide.",
    reverse_transform="Remove acetonide with aqueous acid (pTsOH/MeOH or AcOH/water).",
)

_t(
    name="THP protection of alcohol",
    named_reaction="THP protection",
    category=ReactionCategory.PROTECTION,
    reagents=["DHP (3,4-dihydro-2H-pyran)", "p-TsOH (cat.)"],
    solvents=["DCM"],
    catalysts=["p-TsOH"],
    temperature_range=(0.0, 25.0),
    typical_yield=(80.0, 95.0),
    functional_group_required=["alcohol"],
    functional_group_produced=["acetal"],
    functional_group_incompatible=[],
    mechanism="Acid-catalysed addition of the alcohol to the vinyl ether of DHP forms a mixed THP acetal. Note: creates a new stereocenter.",
    reverse_transform="Remove THP with p-TsOH in MeOH or PPTS in EtOH.",
)

# -----------------------------------------------------------------
#  REARRANGEMENT  (~6)
# -----------------------------------------------------------------

_t(
    name="Cope rearrangement",
    named_reaction="Cope rearrangement",
    category=ReactionCategory.REARRANGEMENT,
    reagents=[],
    solvents=["toluene", "decalin"],
    temperature_range=(150.0, 300.0),
    typical_yield=(50.0, 85.0),
    functional_group_required=["alkene"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=[],
    mechanism="[3,3]-Sigmatropic rearrangement of a 1,5-diene through a chair-like transition state. Oxy-Cope variant (3-hydroxy-1,5-diene) is driven by tautomerisation to a ketone.",
    reverse_transform="Reverse [3,3] shift.",
)

_t(
    name="Pinacol rearrangement",
    named_reaction="Pinacol rearrangement",
    category=ReactionCategory.REARRANGEMENT,
    reagents=["H2SO4"],
    solvents=["water"],
    temperature_range=(25.0, 100.0),
    typical_yield=(55.0, 80.0),
    functional_group_required=["diol"],
    functional_group_produced=["ketone"],
    functional_group_incompatible=[],
    mechanism="Acid protonates one hydroxyl; water departs to give a carbocation; 1,2-alkyl shift to the cationic centre, followed by loss of a proton gives pinacolone (a ketone).",
    reverse_transform="Add two OH groups to the ketone alpha-carbon.",
)

_t(
    name="Curtius rearrangement",
    named_reaction="Curtius rearrangement",
    category=ReactionCategory.REARRANGEMENT,
    reagents=["DPPA (diphenylphosphoryl azide)", "Et3N"],
    solvents=["toluene", "THF"],
    temperature_range=(60.0, 110.0),
    typical_yield=(60.0, 85.0),
    functional_group_required=["carboxylic_acid"],
    functional_group_produced=["primary_amine", "carbamate"],
    functional_group_incompatible=[],
    mechanism="Carboxylic acid is converted to the acyl azide; thermolysis induces 1,2-shift with loss of N2 to give an isocyanate. Trapping with water gives the amine (one fewer carbon); trapping with alcohol gives a carbamate.",
    reverse_transform="Amine + CO2 -> carboxylic acid (formal).",
    safety_notes="Acyl azides are potentially explosive. DPPA generates HN3 (toxic).",
)

_t(
    name="Hofmann rearrangement",
    named_reaction="Hofmann rearrangement",
    category=ReactionCategory.REARRANGEMENT,
    reagents=["NaOH", "Br2"],
    solvents=["water"],
    temperature_range=(0.0, 70.0),
    typical_yield=(55.0, 80.0),
    functional_group_required=["amide"],
    functional_group_produced=["primary_amine"],
    functional_group_incompatible=[],
    mechanism="Bromination of the amide nitrogen, base-induced 1,2-shift with loss of bromide gives an isocyanate, which is hydrolysed to the amine (one fewer carbon).",
    reverse_transform="Amine + CO -> amide (formal).",
)


# -----------------------------------------------------------------
#  MISCELLANEOUS  (~4)
# -----------------------------------------------------------------

_t(
    name="Appel reaction (alcohol to alkyl halide)",
    named_reaction="Appel reaction",
    category=ReactionCategory.MISC,
    reagents=["PPh3", "CCl4 (or CBr4, or I2)"],
    solvents=["DCM"],
    temperature_range=(0.0, 25.0),
    typical_yield=(70.0, 95.0),
    functional_group_required=["alcohol"],
    functional_group_produced=["alkyl_halide"],
    functional_group_incompatible=["carboxylic_acid"],
    mechanism=(
        "PPh3 reacts with CCl4 to form a phosphonium salt; the alcohol oxygen "
        "attacks phosphorus, forming an oxyphosphonium intermediate. Chloride "
        "displaces the activated oxygen via SN2 with inversion."
    ),
    reverse_transform="Disconnect C-X; replace halide with hydroxyl.",
    safety_notes="CCl4 is toxic and an ozone depleter. CBr4 is preferred where possible.",
)

_t(
    name="Staudinger reduction (azide to amine)",
    named_reaction="Staudinger reduction",
    category=ReactionCategory.MISC,
    reagents=["PPh3", "H2O"],
    solvents=["THF", "THF/water"],
    temperature_range=(20.0, 60.0),
    typical_yield=(75.0, 95.0),
    functional_group_required=["azide"],
    functional_group_produced=["primary_amine"],
    functional_group_incompatible=[],
    mechanism=(
        "PPh3 attacks the terminal nitrogen of the azide, releasing N2 and "
        "forming an iminophosphorane (aza-ylide). Hydrolysis gives the "
        "primary amine and Ph3P=O."
    ),
    reverse_transform="Convert amine back to azide via diazotransfer.",
    safety_notes="Organic azides can be shock-sensitive. Ph3P=O byproduct.",
)

_t(
    name="Corey-Chaykovsky epoxidation / cyclopropanation",
    named_reaction="Corey-Chaykovsky reaction",
    category=ReactionCategory.MISC,
    reagents=["trimethylsulfonium iodide (or trimethylsulfoxonium iodide)", "NaH"],
    solvents=["DMSO", "THF"],
    temperature_range=(0.0, 25.0),
    typical_yield=(55.0, 80.0),
    functional_group_required=["aldehyde", "ketone"],
    functional_group_produced=["epoxide", "cyclopropane"],
    functional_group_incompatible=[],
    mechanism=(
        "Deprotonation of the sulfonium salt gives the sulfur ylide, which "
        "attacks the carbonyl. With the sulfonium ylide, the betaine closes "
        "to an epoxide; with the sulfoxonium ylide, cyclopropanation of enones "
        "occurs via conjugate addition then ring closure."
    ),
    reverse_transform="Open the epoxide to retrieve the carbonyl; or open the cyclopropane.",
    safety_notes="NaH is pyrophoric. Handle under inert atmosphere.",
)

_t(
    name="Olefin metathesis (ring-closing)",
    named_reaction="Ring-closing metathesis (RCM)",
    category=ReactionCategory.MISC,
    reagents=["Grubbs 2nd generation catalyst"],
    solvents=["DCM", "toluene"],
    catalysts=["Grubbs 2nd generation"],
    temperature_range=(25.0, 40.0),
    typical_yield=(60.0, 90.0),
    functional_group_required=["alkene"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=["thiol", "primary_amine"],
    mechanism=(
        "Ruthenium carbene catalyst initiates by [2+2] cycloaddition with one "
        "terminal alkene; metallacyclobutane fragmentation releases ethylene. "
        "Intramolecular [2+2] with the second alkene closes the ring."
    ),
    reverse_transform="Open the ring by cross-metathesis; restore two terminal alkenes.",
    safety_notes="Grubbs catalysts are air-sensitive. Ethylene gas is evolved.",
)

# -----------------------------------------------------------------
#  METATHESIS / C-H ACTIVATION / ASYMMETRIC  (~5)
# -----------------------------------------------------------------

_t(
    name="Cross-metathesis (Grubbs catalyst)",
    named_reaction="Cross-metathesis",
    category=ReactionCategory.MISC,
    reagents=["Grubbs 2nd generation catalyst", "terminal alkene partner"],
    solvents=["DCM", "toluene"],
    catalysts=["Grubbs 2nd generation"],
    temperature_range=(25.0, 45.0),
    typical_yield=(40.0, 80.0),
    functional_group_required=["alkene"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=["thiol", "primary_amine"],
    mechanism=(
        "Ruthenium carbene initiates by [2+2] cycloaddition with one "
        "terminal olefin; metallacyclobutane fragmentation and re-entry "
        "with the cross-partner olefin produces the hetero-disubstituted "
        "alkene. Statistical mixtures are minimised by choosing partners "
        "of different reactivity types (Type I / Type II / Type III)."
    ),
    reverse_transform=(
        "Disconnect the internal alkene; each fragment becomes a "
        "terminal olefin."
    ),
    scale_notes=(
        "Catalyst loading typically 1-5 mol%. Ethylene removal by "
        "sparging or reduced pressure shifts equilibrium toward product."
    ),
    safety_notes=(
        "Grubbs catalysts are air-sensitive. Ethylene gas is evolved; "
        "ensure adequate ventilation."
    ),
)

_t(
    name="Pd-catalysed C-H activation (arene functionalisation)",
    named_reaction="Pd-catalysed C-H activation",
    category=ReactionCategory.COUPLING,
    reagents=["aryl halide (or oxidant)", "PivOH (or AcOH)"],
    solvents=["DMA", "toluene", "DMF"],
    catalysts=["Pd(OAc)2"],
    temperature_range=(80.0, 140.0),
    typical_yield=(30.0, 75.0),
    functional_group_required=["aromatic_ring"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=["thiol"],
    mechanism=(
        "Pd(II) coordinates to the arene and undergoes concerted "
        "metalation-deprotonation (CMD) at the least hindered C-H bond, "
        "assisted by carboxylate base. The resulting aryl-Pd species "
        "either undergoes oxidative Heck-type coupling or is oxidised "
        "to Pd(IV) before reductive elimination to form the new C-C bond."
    ),
    reverse_transform=(
        "Disconnect the newly formed aryl-aryl or aryl-R bond; "
        "one fragment keeps H, the other gets a halide or equivalent."
    ),
    scale_notes=(
        "Regioselectivity is substrate-dependent; directing groups "
        "(pyridine, amide, oxazoline) greatly improve selectivity. "
        "Typical Pd loading 5-10 mol%."
    ),
    safety_notes=(
        "Pd(OAc)2 is a skin sensitiser. DMA is a reproductive toxin "
        "(H360). High temperatures require pressure-rated equipment."
    ),
)

_t(
    name="Asymmetric hydrogenation (Rh or Ru catalysed)",
    named_reaction="Asymmetric hydrogenation",
    category=ReactionCategory.REDUCTION,
    reagents=["H2"],
    solvents=["MeOH", "iPrOH", "DCM"],
    catalysts=["Rh(cod)2BF4 / chiral bisphosphine", "Ru(OAc)2(BINAP)"],
    temperature_range=(20.0, 50.0),
    typical_yield=(85.0, 99.0),
    functional_group_required=["alkene"],
    functional_group_produced=["alkane"],
    functional_group_incompatible=["thiol"],
    mechanism=(
        "Chiral transition-metal complex (Rh or Ru with enantiomerically "
        "pure bisphosphine ligand) coordinates the prochiral olefin. "
        "Oxidative addition of H2 followed by migratory insertion and "
        "reductive elimination delivers H2 across the double bond with "
        "enantioselectivities typically >95% ee."
    ),
    reverse_transform=(
        "Disconnect two C-H bonds on adjacent carbons; restore the "
        "alkene with defined geometry."
    ),
    scale_notes=(
        "Highly scalable; used industrially for L-DOPA, menthol, etc. "
        "Catalyst loadings as low as 0.01 mol% on scale. H2 pressures "
        "of 1-50 atm are common."
    ),
    safety_notes=(
        "H2 is extremely flammable (H220). Pressure equipment must "
        "be rated. Chiral ligands are expensive."
    ),
)

_t(
    name="Ring-opening metathesis polymerisation (ROMP)",
    named_reaction="ROMP",
    category=ReactionCategory.POLYMERIZATION,
    reagents=["Grubbs 3rd generation catalyst (or Grubbs 2nd gen)"],
    solvents=["DCM", "THF", "toluene"],
    catalysts=["Grubbs 3rd generation"],
    temperature_range=(20.0, 40.0),
    typical_yield=(80.0, 99.0),
    functional_group_required=["alkene"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=["thiol", "primary_amine"],
    mechanism=(
        "Ruthenium alkylidene initiates by [2+2] cycloaddition with the "
        "strained cyclic olefin (norbornene, cyclooctene, etc.); "
        "metallacyclobutane retro-[2+2] opens the ring and regenerates "
        "the active alkylidene. Repetitive ring-opening of strained "
        "monomers produces a linear polymer with backbone unsaturation. "
        "Thermodynamic driving force is the release of ring strain."
    ),
    reverse_transform=(
        "Depolymerise the unsaturated polymer back to cyclic monomer "
        "under dilute conditions (ring-closing metathesis)."
    ),
    scale_notes=(
        "Living polymerisation: PDI close to 1.0. Control molecular "
        "weight via monomer-to-initiator ratio. Quench with ethyl "
        "vinyl ether."
    ),
    safety_notes=(
        "Strained monomers (e.g. dicyclopentadiene) are flammable and "
        "may self-polymerise. Ru catalysts are air-sensitive."
    ),
)

_t(
    name="Olefin cross-metathesis (selective, with acrylate partner)",
    named_reaction="Olefin cross-metathesis",
    category=ReactionCategory.MISC,
    reagents=["Hoveyda-Grubbs 2nd generation catalyst", "methyl acrylate"],
    solvents=["DCM", "toluene"],
    catalysts=["Hoveyda-Grubbs 2nd generation"],
    temperature_range=(25.0, 50.0),
    typical_yield=(50.0, 85.0),
    functional_group_required=["alkene"],
    functional_group_produced=["alkene", "ester"],
    functional_group_incompatible=["thiol", "primary_amine"],
    mechanism=(
        "Hoveyda-Grubbs catalyst releases the chelating isopropoxy "
        "styrenyl ligand to generate the active 14-electron Ru "
        "alkylidene. [2+2] cycloaddition with the Type I olefin "
        "(terminal alkene) followed by cycloreversion, then a second "
        "[2+2]/retro-[2+2] cycle with the Type II partner (acrylate) "
        "gives the E-cross product selectively. Ethylene by-product "
        "is removed under reduced pressure."
    ),
    reverse_transform=(
        "Disconnect the internal alkene adjacent to the ester; "
        "one fragment is a terminal olefin, the other is the acrylate."
    ),
    scale_notes=(
        "Hoveyda-Grubbs catalyst is more robust than Grubbs 2nd gen "
        "for electron-poor partners. Typical loading 2-5 mol%."
    ),
    safety_notes=(
        "Methyl acrylate is a skin sensitiser and respiratory irritant "
        "(H312/H332). Ru catalysts are air-sensitive."
    ),
)


# -----------------------------------------------------------------
#  HETEROCYCLE FUNCTIONALIZATION  (~4)
# -----------------------------------------------------------------

_t(
    name="N-methylation of secondary amine with methyl iodide",
    named_reaction="N-Alkylation",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["MeI", "K2CO3"],
    solvents=["DMF", "DMSO", "acetone"],
    temperature_range=(25.0, 80.0),
    typical_yield=(70.0, 95.0),
    functional_group_required=["secondary_amine"],
    functional_group_produced=["tertiary_amine"],
    functional_group_incompatible=[],
    mechanism=(
        "Deprotonation of the N-H by K2CO3 generates the nucleophilic "
        "nitrogen, which attacks the electrophilic carbon of methyl iodide "
        "in an SN2 mechanism, forming the N-CH3 bond with loss of iodide."
    ),
    reverse_transform=(
        "Disconnect N-CH3 bond on a tertiary amine to reveal a secondary "
        "amine precursor plus methyl iodide."
    ),
    scale_notes=(
        "MeI is volatile (bp 42 C); use excess K2CO3 to neutralise HI. "
        "Monitor by TLC or NMR for over-alkylation on polyamines."
    ),
    safety_notes=(
        "Methyl iodide is highly toxic and a suspected carcinogen "
        "(H301/H311/H331/H351). Must be handled in a well-ventilated "
        "fume hood with appropriate PPE."
    ),
)

_t(
    name="N-methylation of secondary amine with dimethyl sulfate",
    named_reaction="N-Alkylation (DMS)",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["(MeO)2SO2", "NaOH"],
    solvents=["water", "DCM"],
    temperature_range=(0.0, 25.0),
    typical_yield=(75.0, 92.0),
    functional_group_required=["secondary_amine"],
    functional_group_produced=["tertiary_amine"],
    functional_group_incompatible=["alcohol"],
    mechanism=(
        "Dimethyl sulfate methylates the nitrogen nucleophile. NaOH "
        "neutralises the acidic by-product and destroys excess DMS."
    ),
    reverse_transform=(
        "Disconnect N-CH3 bond on a tertiary amine to reveal a secondary "
        "amine precursor."
    ),
    safety_notes=(
        "Dimethyl sulfate is extremely toxic and a known carcinogen "
        "(H300/H310/H330/H350). Requires maximum PPE and fume hood."
    ),
)

_t(
    name="Amide bond formation from amine and acid chloride",
    named_reaction="Schotten-Baumann acylation",
    category=ReactionCategory.COUPLING,
    reagents=["acid chloride", "Et3N"],
    solvents=["DCM", "THF"],
    temperature_range=(0.0, 25.0),
    typical_yield=(75.0, 95.0),
    functional_group_required=["primary_amine", "secondary_amine"],
    functional_group_produced=["amide"],
    functional_group_incompatible=[],
    mechanism=(
        "Nucleophilic addition-elimination: the amine nitrogen attacks "
        "the electrophilic carbonyl carbon of the acid chloride, "
        "displacing chloride. Triethylamine scavenges HCl."
    ),
    reverse_transform=(
        "Disconnect the amide C-N bond to give an amine and an acid "
        "chloride (or carboxylic acid) as precursors."
    ),
    safety_notes="Acid chlorides are corrosive and moisture-sensitive.",
)

_t(
    name="Amide bond formation with coupling reagent",
    named_reaction="EDC/HOBt amide coupling",
    category=ReactionCategory.COUPLING,
    reagents=["EDC", "HOBt", "DIPEA"],
    solvents=["DMF", "DCM"],
    temperature_range=(0.0, 25.0),
    typical_yield=(70.0, 92.0),
    functional_group_required=["primary_amine", "secondary_amine"],
    functional_group_produced=["amide"],
    functional_group_incompatible=[],
    mechanism=(
        "EDC activates the carboxylic acid as an O-acylisourea, HOBt "
        "converts it to an active ester, and the amine displaces HOBt "
        "to form the amide bond."
    ),
    reverse_transform=(
        "Disconnect the amide C-N bond to give an amine and a carboxylic "
        "acid as precursors."
    ),
    safety_notes="EDC is an irritant. HOBt is shock-sensitive when dry.",
)


# -----------------------------------------------------------------
#  CROSS-COUPLING (additional)  (~13)
# -----------------------------------------------------------------

_t(
    name="Buchwald-Hartwig amination",
    named_reaction="Buchwald-Hartwig amination",
    category=ReactionCategory.COUPLING,
    reagents=["aryl halide", "amine", "NaOtBu"],
    solvents=["toluene", "dioxane"],
    catalysts=["Pd2(dba)3/XPhos", "Pd(OAc)2/BINAP"],
    temperature_range=(80.0, 110.0),
    typical_yield=(70.0, 95.0),
    functional_group_required=["alkyl_halide_br", "primary_amine"],
    functional_group_produced=["secondary_amine"],
    functional_group_incompatible=[],
    mechanism="Pd(0) oxidatively adds to the aryl halide; coordination and deprotonation of the amine gives a Pd-amido species; reductive elimination forms the new C-N bond.",
    reverse_transform="Disconnect the aryl C-N bond; one fragment gets halide, the other is the amine.",
)

_t(
    name="Chan-Lam coupling",
    named_reaction="Chan-Lam coupling",
    category=ReactionCategory.COUPLING,
    reagents=["ArB(OH)2", "amine (or alcohol)", "Et3N"],
    solvents=["DCM", "DMF"],
    catalysts=["Cu(OAc)2"],
    temperature_range=(20.0, 40.0),
    typical_yield=(40.0, 80.0),
    functional_group_required=["boronic_acid", "primary_amine"],
    functional_group_produced=["secondary_amine"],
    functional_group_incompatible=[],
    mechanism="Cu(II)-mediated oxidative coupling of a boronic acid with an N-H or O-H nucleophile under aerobic conditions. Air serves as terminal oxidant.",
    reverse_transform="Disconnect aryl C-N (or C-O) bond; one fragment gets B(OH)2, the other retains N-H.",
)

_t(
    name="Ullmann coupling",
    named_reaction="Ullmann coupling",
    category=ReactionCategory.COUPLING,
    reagents=["aryl halide", "Cu powder (or CuI)", "base"],
    solvents=["DMF", "DMSO"],
    catalysts=["CuI", "1,10-phenanthroline"],
    temperature_range=(80.0, 150.0),
    typical_yield=(50.0, 85.0),
    functional_group_required=["alkyl_halide_br", "alkyl_halide_i"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=[],
    mechanism="Cu(I) undergoes oxidative addition to one aryl halide, then a second aryl halide undergoes transmetalation or sigma-bond metathesis; reductive elimination gives the biaryl.",
    reverse_transform="Disconnect biaryl C-C bond; both fragments get halides.",
)

_t(
    name="Kumada coupling",
    named_reaction="Kumada coupling",
    category=ReactionCategory.COUPLING,
    reagents=["RMgX (Grignard reagent)", "aryl/vinyl halide"],
    solvents=["THF", "Et2O"],
    catalysts=["NiCl2(dppf)", "PdCl2(dppf)"],
    temperature_range=(0.0, 60.0),
    typical_yield=(65.0, 90.0),
    functional_group_required=["alkyl_halide_br", "alkyl_halide_cl"],
    functional_group_produced=[],
    functional_group_incompatible=["aldehyde", "ketone", "ester"],
    mechanism="Ni(0) or Pd(0) oxidatively adds to the halide; transmetalation with the Grignard reagent; reductive elimination forms C-C bond. Limited FG tolerance due to Grignard reactivity.",
    reverse_transform="Disconnect C-C bond; one fragment gets halide, the other becomes RMgX.",
)

_t(
    name="Hiyama coupling",
    named_reaction="Hiyama coupling",
    category=ReactionCategory.COUPLING,
    reagents=["R-Si(OMe)3 (organosilane)", "aryl halide", "TBAF (activator)"],
    solvents=["THF", "DMF"],
    catalysts=["Pd(PPh3)4"],
    temperature_range=(50.0, 100.0),
    typical_yield=(55.0, 85.0),
    functional_group_required=["alkyl_halide_br", "alkyl_halide_i"],
    functional_group_produced=[],
    functional_group_incompatible=[],
    mechanism="Fluoride activates the organosilane to form a pentacoordinate silicate; Pd(0) oxidatively adds to the halide; transmetalation with the activated silane; reductive elimination gives C-C bond.",
    reverse_transform="Disconnect C-C bond; one fragment gets halide, the other gets Si(OMe)3.",
)

_t(
    name="Miyaura borylation",
    named_reaction="Miyaura borylation",
    category=ReactionCategory.COUPLING,
    reagents=["B2pin2 (bis(pinacolato)diboron)", "KOAc"],
    solvents=["dioxane", "DMSO"],
    catalysts=["PdCl2(dppf)"],
    temperature_range=(80.0, 100.0),
    typical_yield=(70.0, 90.0),
    functional_group_required=["alkyl_halide_br", "alkyl_halide_i"],
    functional_group_produced=["boronic_acid"],
    functional_group_incompatible=[],
    mechanism="Pd(0) oxidatively adds to aryl halide; sigma-bond metathesis with B2pin2 installs the boronate ester; reductive elimination gives Ar-Bpin.",
    reverse_transform="Disconnect C-B bond; restore the halide.",
)

_t(
    name="C-H borylation (Ir-catalysed)",
    named_reaction="Ir-catalysed C-H borylation",
    category=ReactionCategory.COUPLING,
    reagents=["B2pin2", "dtbpy (4,4'-di-tert-butyl-2,2'-bipyridine)"],
    solvents=["hexane", "THF"],
    catalysts=["[Ir(cod)OMe]2"],
    temperature_range=(25.0, 80.0),
    typical_yield=(60.0, 90.0),
    functional_group_required=["aromatic_ring"],
    functional_group_produced=["boronic_acid"],
    functional_group_incompatible=[],
    mechanism="Ir(III)-trisboryl complex activates the arene C-H bond via oxidative addition; reductive elimination installs the Bpin group. Regioselectivity is governed by sterics (meta/para to substituents).",
    reverse_transform="Disconnect C-B bond; restore C-H.",
)

_t(
    name="Suzuki coupling (alkyl, Ni-catalysed)",
    named_reaction="Suzuki coupling (alkyl)",
    category=ReactionCategory.COUPLING,
    reagents=["alkyl-9-BBN", "aryl halide", "K3PO4"],
    solvents=["THF/water"],
    catalysts=["NiCl2(dme)/bathophenanthroline"],
    temperature_range=(25.0, 60.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["alkyl_halide_br"],
    functional_group_produced=[],
    functional_group_incompatible=[],
    mechanism="Ni(0) activates the alkyl halide (which resists Pd-catalysed oxidative addition); transmetalation with the alkyl borane; reductive elimination forms C(sp3)-C(sp2) bond.",
    reverse_transform="Disconnect the C(sp3)-C(sp2) bond; one fragment gets halide, the other gets 9-BBN.",
)

_t(
    name="Tsuji-Trost allylation",
    named_reaction="Tsuji-Trost reaction",
    category=ReactionCategory.COUPLING,
    reagents=["allyl acetate (or allyl carbonate)", "nucleophile", "base"],
    solvents=["THF", "DCM"],
    catalysts=["Pd(PPh3)4", "Pd2(dba)3/dppe"],
    temperature_range=(0.0, 60.0),
    typical_yield=(60.0, 90.0),
    functional_group_required=["alkene"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=[],
    mechanism="Pd(0) oxidatively adds to the allyl electrophile forming a pi-allyl-Pd complex; nucleophilic attack on the allyl terminus (soft nucleophiles attack carbon, hard nucleophiles attack Pd) forms the product.",
    reverse_transform="Disconnect the allyl C-nucleophile bond; restore the allyl leaving group.",
)

_t(
    name="Photoredox/Ni dual catalysis C-C coupling",
    named_reaction="Metallaphotoredox coupling",
    category=ReactionCategory.COUPLING,
    reagents=["alkyl radical precursor", "aryl halide"],
    solvents=["DMF", "DMA"],
    catalysts=["Ir(ppy)3", "NiCl2*glyme/dtbbpy"],
    temperature_range=(20.0, 30.0),
    typical_yield=(40.0, 75.0),
    functional_group_required=["alkyl_halide_br"],
    functional_group_produced=[],
    functional_group_incompatible=[],
    mechanism="Ir photocatalyst generates an alkyl radical from a carboxylic acid or alkyl halide via single-electron transfer; Ni(0) oxidatively adds to the aryl halide; radical capture by Ni(II) and reductive elimination form the C(sp3)-C(sp2) bond.",
    reverse_transform="Disconnect the C(sp3)-C(sp2) bond; one fragment becomes a radical precursor.",
)

_t(
    name="Liebeskind-Srogl coupling",
    named_reaction="Liebeskind-Srogl coupling",
    category=ReactionCategory.COUPLING,
    reagents=["thioester", "boronic acid", "CuTC (copper(I) thiophene-2-carboxylate)"],
    solvents=["THF"],
    catalysts=["Pd(PPh3)4"],
    temperature_range=(40.0, 60.0),
    typical_yield=(55.0, 85.0),
    functional_group_required=["ester"],
    functional_group_produced=["ketone"],
    functional_group_incompatible=[],
    mechanism="Pd(0) oxidatively adds to the thioester C-S bond; CuTC mediates transmetalation with the boronic acid; reductive elimination gives the ketone.",
    reverse_transform="Disconnect the ketone C-C bond; one fragment becomes thioester, the other boronic acid.",
)

_t(
    name="Cadiot-Chodkiewicz coupling",
    named_reaction="Cadiot-Chodkiewicz coupling",
    category=ReactionCategory.COUPLING,
    reagents=["terminal alkyne", "1-bromoalkyne", "CuCl", "NH2OH*HCl"],
    solvents=["MeOH/Et2O", "pyridine"],
    catalysts=["CuCl"],
    temperature_range=(0.0, 25.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["alkyne"],
    functional_group_produced=["alkyne"],
    functional_group_incompatible=[],
    mechanism="Cu(I) forms a copper acetylide from the terminal alkyne; oxidative coupling with the 1-bromoalkyne gives an unsymmetrical 1,3-diyne.",
    reverse_transform="Disconnect the diyne central bond; one fragment is a terminal alkyne, the other a 1-bromoalkyne.",
)

# -----------------------------------------------------------------
#  OLEFINATION  (~6)
# -----------------------------------------------------------------

_t(
    name="Julia-Kocienski olefination",
    named_reaction="Julia-Kocienski olefination",
    category=ReactionCategory.CARBONYL,
    reagents=["1-phenyl-1H-tetrazol-5-yl sulfone", "KHMDS"],
    solvents=["THF", "DMF"],
    temperature_range=(-78.0, 0.0),
    typical_yield=(55.0, 85.0),
    functional_group_required=["aldehyde"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=[],
    mechanism="Base deprotonates the sulfone; the carbanion attacks the aldehyde; beta-elimination of SO2 and the heterocycle gives predominantly the E-alkene.",
    reverse_transform="Disconnect C=C; one fragment becomes aldehyde, the other a sulfone.",
)

_t(
    name="Peterson olefination",
    named_reaction="Peterson olefination",
    category=ReactionCategory.CARBONYL,
    reagents=["alpha-silyl carbanion (e.g. TMSCH2Li)", "BF3*OEt2 (or base)"],
    solvents=["THF"],
    temperature_range=(-78.0, 25.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["aldehyde", "ketone"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=[],
    mechanism="Alpha-silyl carbanion adds to the carbonyl; the beta-hydroxy silane eliminates to give the alkene. Acid conditions give Z-alkene; base conditions give E-alkene.",
    reverse_transform="Disconnect C=C; one fragment becomes carbonyl, the other an alpha-silyl halide.",
)

_t(
    name="Tebbe olefination",
    named_reaction="Tebbe olefination",
    category=ReactionCategory.CARBONYL,
    reagents=["Tebbe reagent (Cp2Ti(mu-Cl)(mu-CH2)AlMe2)", "pyridine"],
    solvents=["THF", "toluene"],
    temperature_range=(-40.0, 25.0),
    typical_yield=(60.0, 85.0),
    functional_group_required=["aldehyde", "ketone", "ester"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=[],
    mechanism="Pyridine triggers release of the Schrock carbene Cp2Ti=CH2 from the Tebbe reagent; [2+2] cycloaddition with the carbonyl forms a titanacyclobutane; retro-[2+2] gives the alkene and Cp2Ti=O.",
    reverse_transform="Disconnect C=CH2; restore the carbonyl.",
)

_t(
    name="Aza-Wittig reaction",
    named_reaction="Aza-Wittig reaction",
    category=ReactionCategory.CARBONYL,
    reagents=["iminophosphorane (R3P=NR')", "carbonyl compound"],
    solvents=["THF", "toluene"],
    temperature_range=(25.0, 80.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["aldehyde", "ketone"],
    functional_group_produced=["imine"],
    functional_group_incompatible=[],
    mechanism="Iminophosphorane attacks the carbonyl; a betaine closes to a four-membered aza-oxaphosphetane; retro-[2+2] gives the imine and Ph3P=O, analogous to the Wittig reaction.",
    reverse_transform="Disconnect C=N; one fragment becomes carbonyl, the other an amine precursor.",
)

_t(
    name="Shapiro reaction",
    named_reaction="Shapiro reaction",
    category=ReactionCategory.ELIMINATION,
    reagents=["TsNHNH2 (tosylhydrazide)", "n-BuLi (2 equiv)"],
    solvents=["THF", "TMEDA"],
    temperature_range=(-78.0, 25.0),
    typical_yield=(50.0, 75.0),
    functional_group_required=["ketone"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=["ester"],
    mechanism="Tosylhydrazone is formed from the ketone; double deprotonation by BuLi generates a dianion which loses N2 and tosyl anion to produce a vinyl lithium species, protonated to give the less-substituted alkene.",
    reverse_transform="Disconnect the alkene; restore the ketone via hydration.",
)

_t(
    name="Bamford-Stevens reaction",
    named_reaction="Bamford-Stevens reaction",
    category=ReactionCategory.ELIMINATION,
    reagents=["TsNHNH2", "NaOMe (or NaH)"],
    solvents=["diglyme", "ethylene glycol"],
    temperature_range=(150.0, 200.0),
    typical_yield=(45.0, 70.0),
    functional_group_required=["ketone"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=[],
    mechanism="Tosylhydrazone decomposes thermally under basic conditions via a diazo intermediate and a carbene; 1,2-hydride shift gives the alkene. In protic solvents, carbocation intermediates give Zaitsev products.",
    reverse_transform="Disconnect the alkene; restore the ketone.",
)

# -----------------------------------------------------------------
#  HETEROCYCLE SYNTHESIS  (~14)
# -----------------------------------------------------------------

_t(
    name="Hantzsch pyridine synthesis",
    named_reaction="Hantzsch pyridine synthesis",
    category=ReactionCategory.MISC,
    reagents=["aldehyde", "2 equiv beta-keto ester", "NH4OAc"],
    solvents=["EtOH", "AcOH"],
    temperature_range=(60.0, 80.0),
    typical_yield=(40.0, 70.0),
    functional_group_required=["aldehyde"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=[],
    mechanism="Knoevenagel condensation of aldehyde with one beta-keto ester; Michael addition of the enamine (from the second keto ester + ammonia); cyclisation and aromatisation give the 1,4-dihydropyridine, oxidised to the pyridine.",
    reverse_transform="Disconnect the pyridine ring into aldehyde + two beta-keto esters + ammonia.",
)

_t(
    name="Skraup quinoline synthesis",
    named_reaction="Skraup synthesis",
    category=ReactionCategory.MISC,
    reagents=["glycerol", "aniline", "H2SO4", "oxidant (nitrobenzene)"],
    solvents=["neat"],
    temperature_range=(100.0, 180.0),
    typical_yield=(30.0, 60.0),
    functional_group_required=["primary_amine", "aromatic_ring"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=[],
    mechanism="Glycerol is dehydrated to acrolein by sulfuric acid; conjugate addition of aniline, followed by cyclisation and oxidative aromatisation gives quinoline.",
    reverse_transform="Disconnect quinoline ring at C2-C3; fragments are aniline and acrolein equivalent.",
)

_t(
    name="Doebner-Miller quinoline synthesis",
    named_reaction="Doebner-Miller synthesis",
    category=ReactionCategory.MISC,
    reagents=["aniline", "alpha,beta-unsaturated aldehyde", "acid catalyst"],
    solvents=["AcOH", "EtOH"],
    temperature_range=(80.0, 130.0),
    typical_yield=(40.0, 70.0),
    functional_group_required=["primary_amine", "aromatic_ring"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=[],
    mechanism="Aniline condenses with the alpha,beta-unsaturated aldehyde; conjugate addition, cyclisation, and dehydration give 2-substituted quinoline.",
    reverse_transform="Disconnect quinoline; fragments are aniline and alpha,beta-unsaturated aldehyde.",
)

_t(
    name="Friedlander quinoline synthesis",
    named_reaction="Friedlander synthesis",
    category=ReactionCategory.MISC,
    reagents=["2-aminobenzaldehyde", "ketone", "base (or acid)"],
    solvents=["EtOH", "toluene"],
    temperature_range=(60.0, 120.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["aldehyde", "ketone", "primary_amine"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=[],
    mechanism="Aldol-type condensation of the ketone with the 2-aminobenzaldehyde followed by intramolecular cyclodehydration gives the quinoline.",
    reverse_transform="Disconnect quinoline at C2-C3; fragments are 2-aminobenzaldehyde and ketone.",
)

_t(
    name="Paal-Knorr pyrrole synthesis",
    named_reaction="Paal-Knorr pyrrole synthesis",
    category=ReactionCategory.MISC,
    reagents=["1,4-dicarbonyl compound", "primary amine (or NH4OAc)"],
    solvents=["AcOH", "toluene"],
    temperature_range=(25.0, 110.0),
    typical_yield=(50.0, 85.0),
    functional_group_required=["ketone"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=[],
    mechanism="The amine condenses with both carbonyls of the 1,4-diketone; double cyclodehydration gives the pyrrole ring.",
    reverse_transform="Disconnect pyrrole ring; fragments are 1,4-diketone and amine.",
)

_t(
    name="Paal-Knorr furan synthesis",
    named_reaction="Paal-Knorr furan synthesis",
    category=ReactionCategory.MISC,
    reagents=["1,4-dicarbonyl compound", "acid catalyst (H2SO4 or P2O5)"],
    solvents=["toluene", "neat"],
    temperature_range=(80.0, 150.0),
    typical_yield=(40.0, 75.0),
    functional_group_required=["ketone"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=[],
    mechanism="Acid-catalysed intramolecular cyclodehydration of a 1,4-diketone; loss of two equivalents of water gives the furan ring.",
    reverse_transform="Disconnect furan ring into 1,4-diketone.",
)

_t(
    name="Fischer indole synthesis",
    named_reaction="Fischer indole synthesis",
    category=ReactionCategory.MISC,
    reagents=["aryl hydrazine", "aldehyde/ketone", "acid catalyst (ZnCl2 or BF3)"],
    solvents=["AcOH", "EtOH"],
    temperature_range=(80.0, 180.0),
    typical_yield=(40.0, 75.0),
    functional_group_required=["ketone", "aromatic_ring"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=[],
    mechanism="Phenylhydrazine condenses with ketone to form phenylhydrazone; acid-catalysed [3,3]-sigmatropic rearrangement, re-aromatisation, and loss of NH3 give the indole.",
    reverse_transform="Disconnect indole at C2-C3; fragments are aryl hydrazine and ketone.",
)

_t(
    name="Madelung indole synthesis",
    named_reaction="Madelung indole synthesis",
    category=ReactionCategory.MISC,
    reagents=["N-acyl-o-toluidine", "strong base (NaNH2 or BuLi)"],
    solvents=["THF", "toluene"],
    temperature_range=(200.0, 350.0),
    typical_yield=(30.0, 60.0),
    functional_group_required=["amide", "aromatic_ring"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=[],
    mechanism="Strong base deprotonates the benzylic methyl of the N-acyl-o-toluidine; intramolecular cyclisation onto the amide carbonyl forms the indole after dehydration.",
    reverse_transform="Disconnect indole ring; fragments are o-toluidine and acyl group.",
)

_t(
    name="Leimgruber-Batcho indole synthesis",
    named_reaction="Leimgruber-Batcho indole synthesis",
    category=ReactionCategory.MISC,
    reagents=["2-nitrotoluene", "DMF-DMA (dimethylformamide dimethyl acetal)"],
    solvents=["DMF", "toluene"],
    temperature_range=(100.0, 150.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["aromatic_ring", "nitro"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=[],
    mechanism="DMF-DMA condenses with the methyl group of the nitrotoluene to form an enamine; reductive cyclisation (Fe/AcOH or catalytic hydrogenation) gives the indole.",
    reverse_transform="Disconnect indole; fragments are 2-nitrotoluene and one-carbon unit.",
)

_t(
    name="Chichibabin pyridine synthesis",
    named_reaction="Chichibabin synthesis",
    category=ReactionCategory.MISC,
    reagents=["aldehyde (3 equiv)", "NH3"],
    solvents=["neat", "AcOH"],
    temperature_range=(200.0, 350.0),
    typical_yield=(20.0, 50.0),
    functional_group_required=["aldehyde"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=[],
    mechanism="Three molecules of aldehyde condense with ammonia via sequential aldol-type reactions and cyclodehydration to give a 2,4,6-trisubstituted pyridine.",
    reverse_transform="Disconnect pyridine into three aldehyde fragments and ammonia.",
)

_t(
    name="Knorr pyrazole synthesis",
    named_reaction="Knorr pyrazole synthesis",
    category=ReactionCategory.MISC,
    reagents=["1,3-diketone", "hydrazine"],
    solvents=["EtOH", "AcOH"],
    temperature_range=(25.0, 80.0),
    typical_yield=(60.0, 90.0),
    functional_group_required=["ketone"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=[],
    mechanism="Hydrazine condenses with both carbonyls of the 1,3-diketone; double cyclodehydration gives the pyrazole ring.",
    reverse_transform="Disconnect pyrazole; fragments are 1,3-diketone and hydrazine.",
)

_t(
    name="Gewald reaction (2-aminothiophene synthesis)",
    named_reaction="Gewald reaction",
    category=ReactionCategory.MISC,
    reagents=["ketone", "alpha-cyano ester (or malononitrile)", "elemental sulfur", "Et3N"],
    solvents=["EtOH", "DMF"],
    temperature_range=(50.0, 80.0),
    typical_yield=(40.0, 75.0),
    functional_group_required=["ketone"],
    functional_group_produced=["aromatic_ring", "primary_amine"],
    functional_group_incompatible=[],
    mechanism="Knoevenagel condensation followed by sulfur insertion via Gewald cyclisation; the 2-aminothiophene is obtained after base-mediated ring closure.",
    reverse_transform="Disconnect the thiophene ring; fragments are ketone, active methylene compound, and sulfur.",
)

_t(
    name="Combes quinoline synthesis",
    named_reaction="Combes synthesis",
    category=ReactionCategory.MISC,
    reagents=["aniline", "1,3-diketone", "acid catalyst"],
    solvents=["AcOH", "neat"],
    temperature_range=(100.0, 150.0),
    typical_yield=(40.0, 70.0),
    functional_group_required=["primary_amine", "aromatic_ring", "ketone"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=[],
    mechanism="Aniline condenses with the 1,3-diketone to form a Schiff base; acid-catalysed electrophilic cyclisation and dehydration give the quinoline.",
    reverse_transform="Disconnect quinoline; fragments are aniline and 1,3-diketone.",
)

# -----------------------------------------------------------------
#  AROMATIC FUNCTIONALIZATION  (~13)
# -----------------------------------------------------------------

_t(
    name="Friedel-Crafts alkylation",
    named_reaction="Friedel-Crafts alkylation",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["alkyl halide", "AlCl3 (or FeCl3)"],
    solvents=["DCM", "CS2", "nitrobenzene"],
    catalysts=["AlCl3"],
    temperature_range=(0.0, 80.0),
    typical_yield=(40.0, 80.0),
    functional_group_required=["aromatic_ring", "alkyl_halide"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=["nitro", "nitrile"],
    mechanism="AlCl3 generates carbocation from alkyl halide; electrophilic aromatic substitution installs the alkyl group. Rearrangement-prone; over-alkylation is common.",
    reverse_transform="Disconnect the aryl-alkyl C-C bond; restore halide on the alkyl fragment.",
)

_t(
    name="Friedel-Crafts acylation",
    named_reaction="Friedel-Crafts acylation",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["acyl chloride (or anhydride)", "AlCl3 (>1 equiv)"],
    solvents=["DCM", "CS2", "nitrobenzene"],
    catalysts=["AlCl3"],
    temperature_range=(0.0, 80.0),
    typical_yield=(60.0, 90.0),
    functional_group_required=["aromatic_ring"],
    functional_group_produced=["ketone"],
    functional_group_incompatible=["nitro", "primary_amine"],
    mechanism="AlCl3 activates the acyl chloride to form an acylium cation; electrophilic aromatic substitution installs the C=O group. No rearrangement (unlike alkylation). Requires >1 equiv AlCl3 (coordinates to product).",
    reverse_transform="Disconnect aryl-C(=O) bond; fragments are arene and acyl chloride.",
)

_t(
    name="Vilsmeier-Haack formylation",
    named_reaction="Vilsmeier-Haack reaction",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["DMF", "POCl3"],
    solvents=["DCM", "DMF (excess as solvent)"],
    temperature_range=(0.0, 80.0),
    typical_yield=(50.0, 85.0),
    functional_group_required=["aromatic_ring"],
    functional_group_produced=["aldehyde"],
    functional_group_incompatible=["nitro"],
    mechanism="POCl3 activates DMF to form the Vilsmeier reagent (a chloroiminium ion); electrophilic attack on the electron-rich arene, followed by hydrolysis gives the aryl aldehyde.",
    reverse_transform="Disconnect aryl-CHO bond; restore the arene.",
)

_t(
    name="Sandmeyer reaction",
    named_reaction="Sandmeyer reaction",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["NaNO2", "HX (HCl, HBr)", "CuX (CuCl, CuBr, CuCN)"],
    solvents=["water", "aqueous HCl"],
    temperature_range=(0.0, 5.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["primary_amine", "aromatic_ring"],
    functional_group_produced=["alkyl_halide_cl", "alkyl_halide_br", "nitrile"],
    functional_group_incompatible=[],
    mechanism="Diazotisation of aniline with NaNO2/HX forms the diazonium salt; CuX mediates radical substitution of N2 by halide or CN.",
    reverse_transform="Disconnect aryl-X; restore the amine (via diazotisation reverse).",
    safety_notes="Diazonium salts are explosive when dry. Keep at 0-5 C in solution.",
)

_t(
    name="Kolbe-Schmitt carboxylation",
    named_reaction="Kolbe-Schmitt reaction",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["CO2 (high pressure)", "NaOH"],
    solvents=["neat (sodium phenoxide melt)"],
    temperature_range=(100.0, 200.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["aromatic_ring", "alcohol"],
    functional_group_produced=["carboxylic_acid"],
    functional_group_incompatible=[],
    mechanism="Sodium phenoxide reacts with CO2 under pressure; electrophilic carboxylation at the ortho position (Na) or para position (K); acidification liberates salicylic acid (or p-hydroxybenzoic acid).",
    reverse_transform="Disconnect aryl-COOH; restore phenol.",
)

_t(
    name="Minisci reaction",
    named_reaction="Minisci reaction",
    category=ReactionCategory.RADICAL,
    reagents=["carboxylic acid (or alkyl halide)", "AgNO3 (cat.)", "(NH4)2S2O8"],
    solvents=["water/TFA", "DCM/TFA"],
    temperature_range=(25.0, 80.0),
    typical_yield=(30.0, 70.0),
    functional_group_required=["aromatic_ring"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=[],
    mechanism="Silver-catalysed oxidative decarboxylation generates an alkyl radical; the radical adds to the protonated heteroarene (pyridine, quinoline) in a Minisci-type substitution.",
    reverse_transform="Disconnect the alkyl-heteroaryl C-C bond.",
)

_t(
    name="Duff reaction (formylation with hexamine)",
    named_reaction="Duff reaction",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["hexamethylenetetramine (hexamine)", "AcOH (or TFA)"],
    solvents=["AcOH", "TFA"],
    temperature_range=(80.0, 120.0),
    typical_yield=(30.0, 60.0),
    functional_group_required=["aromatic_ring", "alcohol"],
    functional_group_produced=["aldehyde"],
    functional_group_incompatible=["nitro"],
    mechanism="Hexamine is hydrolysed under acidic conditions to generate formaldehyde equivalents; electrophilic aromatic formylation of activated (phenolic) arenes gives the salicylaldehyde.",
    reverse_transform="Disconnect aryl-CHO; restore the phenol.",
)

_t(
    name="Gattermann-Koch formylation",
    named_reaction="Gattermann-Koch reaction",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["CO", "HCl", "AlCl3", "CuCl (cat.)"],
    solvents=["neat"],
    temperature_range=(0.0, 25.0),
    typical_yield=(40.0, 70.0),
    functional_group_required=["aromatic_ring"],
    functional_group_produced=["aldehyde"],
    functional_group_incompatible=["alcohol", "primary_amine"],
    mechanism="CO and HCl form formyl chloride (HCOCl) in situ; AlCl3-mediated Friedel-Crafts acylation installs the formyl group.",
    reverse_transform="Disconnect aryl-CHO; restore the arene.",
    safety_notes="CO is extremely toxic; requires pressure equipment.",
)

_t(
    name="Birch alkylation",
    named_reaction="Birch alkylation",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["Na (or Li)", "NH3 (liquid)", "t-BuOH", "alkyl halide"],
    solvents=["liquid NH3"],
    temperature_range=(-78.0, -33.0),
    typical_yield=(40.0, 70.0),
    functional_group_required=["aromatic_ring", "alkyl_halide"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=[],
    mechanism="Birch reduction generates a dienolate intermediate which is trapped by the alkyl halide electrophile to give an alkylated 1,4-cyclohexadiene.",
    reverse_transform="Disconnect alkyl group from the cyclohexadienyl system.",
    safety_notes="Liquid NH3 and alkali metals: extreme hazard.",
)

_t(
    name="Balz-Schiemann fluorination",
    named_reaction="Balz-Schiemann reaction",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["NaNO2", "HBF4"],
    solvents=["water", "aqueous HBF4"],
    temperature_range=(0.0, 5.0),
    typical_yield=(40.0, 70.0),
    functional_group_required=["primary_amine", "aromatic_ring"],
    functional_group_produced=["alkyl_halide"],
    functional_group_incompatible=[],
    mechanism="Diazotisation of the aniline with NaNO2/HBF4 gives aryl diazonium tetrafluoroborate; thermal decomposition releases N2 and BF3, installing fluorine on the ring.",
    reverse_transform="Disconnect aryl-F; restore amine.",
    safety_notes="Diazonium tetrafluoroborates are shock-sensitive solids.",
)

_t(
    name="Reimer-Tiemann formylation",
    named_reaction="Reimer-Tiemann reaction",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["CHCl3", "NaOH (aq.)"],
    solvents=["water/CHCl3"],
    temperature_range=(60.0, 70.0),
    typical_yield=(20.0, 50.0),
    functional_group_required=["aromatic_ring", "alcohol"],
    functional_group_produced=["aldehyde"],
    functional_group_incompatible=[],
    mechanism="NaOH deprotonates CHCl3 to generate dichlorocarbene (:CCl2); electrophilic attack on the phenoxide, followed by hydrolysis of the gem-dichloride gives the ortho-hydroxybenzaldehyde.",
    reverse_transform="Disconnect aryl-CHO; restore phenol.",
    safety_notes="CHCl3 is toxic; dichlorocarbene is very reactive.",
)

_t(
    name="Schiemann fluorination via Balz-Schiemann",
    named_reaction="Schiemann reaction",
    category=ReactionCategory.SUBSTITUTION,
    reagents=["NaNO2", "HBF4"],
    solvents=["water"],
    temperature_range=(0.0, 5.0),
    typical_yield=(35.0, 65.0),
    functional_group_required=["primary_amine", "aromatic_ring"],
    functional_group_produced=["alkyl_halide"],
    functional_group_incompatible=[],
    mechanism="Diazonium tetrafluoroborate decomposes thermally to aryl fluoride, N2, and BF3.",
    reverse_transform="Disconnect Ar-F; restore ArNH2.",
)

# -----------------------------------------------------------------
#  MULTICOMPONENT REACTIONS  (~6)
# -----------------------------------------------------------------

_t(
    name="Strecker amino acid synthesis",
    named_reaction="Strecker synthesis",
    category=ReactionCategory.MISC,
    reagents=["aldehyde", "NH4Cl", "NaCN"],
    solvents=["water", "MeOH/water"],
    temperature_range=(0.0, 25.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["aldehyde"],
    functional_group_produced=["primary_amine", "nitrile"],
    functional_group_incompatible=[],
    mechanism="Ammonia condenses with aldehyde to form an imine; cyanide adds to the imine giving an alpha-aminonitrile; hydrolysis of the nitrile gives the alpha-amino acid.",
    reverse_transform="Disconnect alpha-amino acid into aldehyde + ammonia + HCN.",
)

_t(
    name="Mannich reaction",
    named_reaction="Mannich reaction",
    category=ReactionCategory.MISC,
    reagents=["formaldehyde", "secondary amine (or NH4+)", "ketone"],
    solvents=["water", "EtOH", "AcOH"],
    temperature_range=(20.0, 80.0),
    typical_yield=(40.0, 75.0),
    functional_group_required=["ketone"],
    functional_group_produced=["tertiary_amine"],
    functional_group_incompatible=[],
    mechanism="Iminium ion (from amine + formaldehyde) is attacked by the enol form of the ketone; the beta-aminoketone (Mannich base) is formed.",
    reverse_transform="Disconnect C-C bond beta to the amine; fragments are ketone, formaldehyde, and amine.",
)

_t(
    name="Passerini reaction",
    named_reaction="Passerini reaction",
    category=ReactionCategory.MISC,
    reagents=["carboxylic acid", "aldehyde", "isocyanide"],
    solvents=["DCM", "MeOH"],
    temperature_range=(20.0, 40.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["aldehyde", "carboxylic_acid"],
    functional_group_produced=["ester", "amide"],
    functional_group_incompatible=[],
    mechanism="Three-component reaction: aldehyde activates with carboxylic acid; isocyanide alpha-addition to the activated carbonyl; acyl transfer gives an alpha-acyloxyamide.",
    reverse_transform="Disconnect into aldehyde, carboxylic acid, and isocyanide.",
)

_t(
    name="Ugi four-component reaction",
    named_reaction="Ugi reaction",
    category=ReactionCategory.MISC,
    reagents=["aldehyde", "amine", "carboxylic acid", "isocyanide"],
    solvents=["MeOH", "DCM"],
    temperature_range=(20.0, 40.0),
    typical_yield=(40.0, 80.0),
    functional_group_required=["aldehyde", "carboxylic_acid"],
    functional_group_produced=["amide"],
    functional_group_incompatible=[],
    mechanism="Imine (from aldehyde + amine) is activated by the carboxylic acid; isocyanide attacks the iminium ion; intramolecular Mumm rearrangement gives the bis-amide product.",
    reverse_transform="Disconnect into aldehyde, amine, carboxylic acid, and isocyanide.",
)

_t(
    name="Baylis-Hillman reaction",
    named_reaction="Baylis-Hillman reaction",
    category=ReactionCategory.MISC,
    reagents=["aldehyde", "activated alkene (acrylate, MVK)", "DABCO (cat.)"],
    solvents=["DMF", "neat", "MeOH"],
    catalysts=["DABCO"],
    temperature_range=(20.0, 50.0),
    typical_yield=(40.0, 80.0),
    functional_group_required=["aldehyde", "alkene"],
    functional_group_produced=["alcohol", "alkene"],
    functional_group_incompatible=[],
    mechanism="DABCO adds to the activated alkene in 1,4-fashion; the resulting zwitterion attacks the aldehyde; proton transfer and elimination of DABCO gives the alpha-methylene-beta-hydroxy product.",
    reverse_transform="Disconnect C-C bond between the aldehyde-derived carbon and the acrylate.",
)

_t(
    name="Henry (nitroaldol) reaction",
    named_reaction="Henry reaction",
    category=ReactionCategory.MISC,
    reagents=["nitroalkane", "aldehyde", "base (Et3N or NaOH)"],
    solvents=["THF", "MeOH", "water"],
    temperature_range=(-20.0, 25.0),
    typical_yield=(50.0, 85.0),
    functional_group_required=["aldehyde", "nitro"],
    functional_group_produced=["alcohol", "nitro"],
    functional_group_incompatible=[],
    mechanism="Base deprotonates the nitroalkane to form a nitronate; nucleophilic addition to the aldehyde gives a beta-nitro alcohol. Dehydration gives the nitroalkene (Henry-Nef).",
    reverse_transform="Disconnect C-C bond between the alcohol carbon and the nitro-bearing carbon.",
)

# -----------------------------------------------------------------
#  CONDENSATION (additional)  (~3)
# -----------------------------------------------------------------

_t(
    name="Knoevenagel condensation",
    named_reaction="Knoevenagel condensation",
    category=ReactionCategory.CARBONYL,
    reagents=["active methylene compound (malonate, cyanoacetate)", "amine catalyst (piperidine)"],
    solvents=["EtOH", "toluene"],
    catalysts=["piperidine"],
    temperature_range=(25.0, 80.0),
    typical_yield=(50.0, 85.0),
    functional_group_required=["aldehyde"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=[],
    mechanism="Piperidine catalyses formation of an iminium ion from the aldehyde; active methylene compound attacks; elimination of water gives the alpha,beta-unsaturated product.",
    reverse_transform="Disconnect the C=C; fragments are aldehyde and active methylene compound.",
)

_t(
    name="Weinreb amide formation and ketone synthesis",
    named_reaction="Weinreb amide",
    category=ReactionCategory.CARBONYL,
    reagents=["N,O-dimethylhydroxylamine (HCl salt)", "EDC/HOBt (or acyl chloride)"],
    solvents=["DCM", "THF"],
    temperature_range=(0.0, 25.0),
    typical_yield=(70.0, 90.0),
    functional_group_required=["carboxylic_acid"],
    functional_group_produced=["amide"],
    functional_group_incompatible=[],
    mechanism="Carboxylic acid is activated and coupled with N,O-dimethylhydroxylamine to give the Weinreb amide. Treatment with RMgBr or RLi gives the ketone (the tetrahedral intermediate is stabilised by chelation).",
    reverse_transform="Disconnect the ketone C-C bond; fragments are Weinreb amide and organometallic.",
)

_t(
    name="Reformatsky reaction",
    named_reaction="Reformatsky reaction",
    category=ReactionCategory.CARBONYL,
    reagents=["alpha-bromo ester", "Zn dust"],
    solvents=["THF", "Et2O"],
    temperature_range=(25.0, 80.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["aldehyde", "ketone"],
    functional_group_produced=["alcohol", "ester"],
    functional_group_incompatible=[],
    mechanism="Zinc inserts into the C-Br bond of the alpha-bromo ester to form a zinc enolate (Reformatsky reagent); addition to the aldehyde or ketone gives a beta-hydroxy ester.",
    reverse_transform="Disconnect C-C bond; fragments are alpha-bromo ester (Zn) and aldehyde/ketone.",
)

# -----------------------------------------------------------------
#  REARRANGEMENT (additional)  (~5)
# -----------------------------------------------------------------

_t(
    name="Wolff rearrangement",
    named_reaction="Wolff rearrangement",
    category=ReactionCategory.REARRANGEMENT,
    reagents=["alpha-diazo ketone", "hv (or Ag2O/heat)"],
    solvents=["MeOH", "water", "THF"],
    temperature_range=(20.0, 80.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["ketone"],
    functional_group_produced=["carboxylic_acid", "ester"],
    functional_group_incompatible=[],
    mechanism="Photolytic or thermal loss of N2 from the alpha-diazo ketone generates a ketocarbene; 1,2-shift gives a ketene which is trapped by MeOH (ester) or H2O (carboxylic acid). Part of the Arndt-Eistert homologation.",
    reverse_transform="Disconnect the ester/acid; precursor is an alpha-diazo ketone.",
)

_t(
    name="Favorskii rearrangement",
    named_reaction="Favorskii rearrangement",
    category=ReactionCategory.REARRANGEMENT,
    reagents=["alpha-halo ketone", "NaOMe (or NaOH)"],
    solvents=["MeOH", "water"],
    temperature_range=(0.0, 60.0),
    typical_yield=(45.0, 75.0),
    functional_group_required=["ketone", "alkyl_halide"],
    functional_group_produced=["ester", "carboxylic_acid"],
    functional_group_incompatible=[],
    mechanism="Base abstracts alpha-proton; enolate displaces halide intramolecularly to form a cyclopropanone; nucleophilic opening by methoxide gives a ring-contracted ester.",
    reverse_transform="Disconnect the ester; precursor is an alpha-halo ketone.",
)

_t(
    name="Fries rearrangement",
    named_reaction="Fries rearrangement",
    category=ReactionCategory.REARRANGEMENT,
    reagents=["phenyl ester", "AlCl3"],
    solvents=["nitrobenzene", "CS2"],
    catalysts=["AlCl3"],
    temperature_range=(25.0, 160.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["ester", "aromatic_ring"],
    functional_group_produced=["ketone", "alcohol"],
    functional_group_incompatible=[],
    mechanism="Lewis acid cleaves the ester bond; the acylium fragment undergoes intramolecular Friedel-Crafts acylation on the phenol ring. Low temperature gives para product; high temperature gives ortho product.",
    reverse_transform="Disconnect the aryl ketone; precursor is the phenyl ester.",
)

_t(
    name="Schmidt reaction",
    named_reaction="Schmidt reaction",
    category=ReactionCategory.REARRANGEMENT,
    reagents=["HN3 (hydrazoic acid)", "H2SO4"],
    solvents=["CHCl3", "DCM"],
    temperature_range=(0.0, 40.0),
    typical_yield=(40.0, 75.0),
    functional_group_required=["ketone"],
    functional_group_produced=["amide"],
    functional_group_incompatible=["primary_amine"],
    mechanism="HN3 adds to the protonated ketone; 1,2-alkyl shift (analogous to Beckmann) with loss of N2 gives the amide (lactam for cyclic ketones).",
    reverse_transform="Disconnect the amide C-N bond; precursor is the ketone.",
    safety_notes="HN3 is extremely toxic and explosive. Use DPPA (safer) when possible.",
)

_t(
    name="Overman rearrangement",
    named_reaction="Overman rearrangement",
    category=ReactionCategory.REARRANGEMENT,
    reagents=["allylic trichloroacetimidate", "heat (or Hg(II) cat.)"],
    solvents=["toluene", "xylene"],
    temperature_range=(100.0, 160.0),
    typical_yield=(60.0, 85.0),
    functional_group_required=["alkene", "alcohol"],
    functional_group_produced=["primary_amine"],
    functional_group_incompatible=[],
    mechanism="[3,3]-Sigmatropic rearrangement of the allylic trichloroacetimidate gives an allylic trichloroacetamide with 1,3-transposition of the nitrogen. Suprafacial chirality transfer.",
    reverse_transform="Disconnect the C-N bond; restore the allylic alcohol and trichloroacetimidate.",
)

# -----------------------------------------------------------------
#  OXIDATION (additional)  (~8)
# -----------------------------------------------------------------

_t(
    name="TPAP (Ley-Griffith) oxidation",
    named_reaction="TPAP oxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["TPAP (Pr4NRuO4, cat.)", "NMO (co-oxidant)", "4A molecular sieves"],
    solvents=["DCM"],
    catalysts=["TPAP"],
    temperature_range=(20.0, 25.0),
    typical_yield=(80.0, 95.0),
    functional_group_required=["alcohol"],
    functional_group_produced=["aldehyde", "ketone"],
    functional_group_incompatible=["primary_amine"],
    mechanism="Catalytic TPAP oxidises the alcohol to aldehyde/ketone; NMO re-oxidises Ru(V) to Ru(VII). 4A molecular sieves remove water. Mild, selective, no over-oxidation of primary alcohols.",
    reverse_transform="Reduce aldehyde/ketone to alcohol.",
)

_t(
    name="IBX oxidation",
    named_reaction="IBX oxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["IBX (2-iodoxybenzoic acid)"],
    solvents=["DMSO"],
    temperature_range=(20.0, 40.0),
    typical_yield=(80.0, 95.0),
    functional_group_required=["alcohol"],
    functional_group_produced=["aldehyde", "ketone"],
    functional_group_incompatible=[],
    mechanism="Hypervalent iodine(V) oxidises alcohols to carbonyls in DMSO. Selective for primary alcohols to aldehydes (no over-oxidation). Can also oxidise alpha-position of carbonyls to enones.",
    reverse_transform="Reduce carbonyl to alcohol.",
    safety_notes="IBX is shock-sensitive as a dry solid. Keep in solution.",
)

_t(
    name="Pinnick (Lindgren) oxidation",
    named_reaction="Pinnick oxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["NaClO2", "NaH2PO4 (buffer)", "2-methyl-2-butene (scavenger)"],
    solvents=["t-BuOH/water", "THF/water"],
    temperature_range=(0.0, 25.0),
    typical_yield=(80.0, 95.0),
    functional_group_required=["aldehyde"],
    functional_group_produced=["carboxylic_acid"],
    functional_group_incompatible=["alkene"],
    mechanism="NaClO2 selectively oxidises aldehydes to carboxylic acids; 2-methyl-2-butene scavenges the HOCl byproduct to prevent side reactions. Highly chemoselective.",
    reverse_transform="Reduce carboxylic acid to aldehyde (e.g. DIBAL-H).",
)

_t(
    name="Riley (SeO2) oxidation",
    named_reaction="Riley oxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["SeO2"],
    solvents=["dioxane", "t-BuOH"],
    temperature_range=(60.0, 100.0),
    typical_yield=(40.0, 70.0),
    functional_group_required=["alkene"],
    functional_group_produced=["alcohol", "aldehyde"],
    functional_group_incompatible=[],
    mechanism="SeO2 performs an ene reaction with the alkene, followed by [2,3]-sigmatropic rearrangement of the allylseleninic acid; hydrolysis gives the allylic alcohol. Further oxidation can give the enone.",
    reverse_transform="Disconnect the allylic C-OH; restore C-H.",
    safety_notes="Selenium compounds are toxic; generate selenium waste carefully.",
)

_t(
    name="Dakin oxidation",
    named_reaction="Dakin oxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["H2O2", "NaOH (or Na2CO3)"],
    solvents=["water", "MeOH/water"],
    temperature_range=(0.0, 25.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["aldehyde", "aromatic_ring"],
    functional_group_produced=["alcohol"],
    functional_group_incompatible=[],
    mechanism="H2O2 attacks the electron-poor carbonyl of an aryl aldehyde; 1,2-aryl shift (Baeyer-Villiger-like) gives a formate ester; hydrolysis under basic conditions gives the phenol.",
    reverse_transform="Disconnect the phenol OH; restore the aldehyde.",
)

_t(
    name="TEMPO-mediated oxidation",
    named_reaction="TEMPO oxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["TEMPO (cat.)", "NaOCl (bleach, terminal oxidant)", "NaBr (co-cat.)"],
    solvents=["DCM/water (biphasic)"],
    catalysts=["TEMPO"],
    temperature_range=(0.0, 10.0),
    typical_yield=(85.0, 98.0),
    functional_group_required=["alcohol"],
    functional_group_produced=["aldehyde"],
    functional_group_incompatible=["primary_amine"],
    mechanism="TEMPO (a stable nitroxide radical) is oxidised to the oxoammonium cation, which selectively oxidises primary alcohols to aldehydes (Anelli protocol). NaBr mediates the oxidant cycle.",
    reverse_transform="Reduce aldehyde to primary alcohol.",
)

_t(
    name="Oppenauer oxidation",
    named_reaction="Oppenauer oxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["Al(OiPr)3 (or Al(OtBu)3)", "acetone (hydrogen acceptor)"],
    solvents=["toluene", "acetone (excess)"],
    catalysts=["Al(OiPr)3"],
    temperature_range=(25.0, 80.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["alcohol"],
    functional_group_produced=["ketone"],
    functional_group_incompatible=[],
    mechanism="Aluminium alkoxide transfers a hydride from the substrate alcohol to acetone via a six-membered transition state (Meerwein-Ponndorf-Verley reverse). Equilibrium driven by excess acetone.",
    reverse_transform="Reduce the ketone back to alcohol (MPV reduction).",
)

_t(
    name="Lemieux-Johnson oxidation (OsO4/NaIO4)",
    named_reaction="Lemieux-Johnson oxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["OsO4 (cat.)", "NaIO4"],
    solvents=["dioxane/water", "THF/water"],
    catalysts=["OsO4"],
    temperature_range=(0.0, 25.0),
    typical_yield=(70.0, 90.0),
    functional_group_required=["alkene"],
    functional_group_produced=["aldehyde", "ketone"],
    functional_group_incompatible=[],
    mechanism="OsO4 dihydroxylates the alkene (syn); NaIO4 cleaves the 1,2-diol to give two carbonyl fragments. One-pot oxidative cleavage.",
    reverse_transform="Join two carbonyl fragments via Wittig or metathesis.",
    safety_notes="OsO4 is extremely toxic. Use catalytic amounts.",
)

# -----------------------------------------------------------------
#  REDUCTION (additional)  (~5)
# -----------------------------------------------------------------

_t(
    name="Rosenmund reduction",
    named_reaction="Rosenmund reduction",
    category=ReactionCategory.REDUCTION,
    reagents=["H2", "BaSO4-poisoned Pd (Rosenmund catalyst)"],
    solvents=["toluene", "xylene"],
    catalysts=["Pd/BaSO4 + quinoline-S poison"],
    temperature_range=(25.0, 60.0),
    typical_yield=(55.0, 80.0),
    functional_group_required=["acid_chloride"],
    functional_group_produced=["aldehyde"],
    functional_group_incompatible=["alkene"],
    mechanism="Selective reduction of acyl chloride to aldehyde using a poisoned Pd catalyst that prevents over-reduction to the alcohol.",
    reverse_transform="Oxidise the aldehyde to acyl chloride.",
)

_t(
    name="Lindlar hydrogenation",
    named_reaction="Lindlar hydrogenation",
    category=ReactionCategory.REDUCTION,
    reagents=["H2"],
    solvents=["EtOAc", "hexane", "MeOH"],
    catalysts=["Pd/CaCO3 + Pb(OAc)2 + quinoline (Lindlar catalyst)"],
    temperature_range=(20.0, 30.0),
    typical_yield=(75.0, 95.0),
    functional_group_required=["alkyne"],
    functional_group_produced=["alkene"],
    functional_group_incompatible=[],
    mechanism="Poisoned Pd catalyst selectively reduces alkynes to cis-alkenes by syn addition of H2. The catalyst poison prevents further reduction to the alkane.",
    reverse_transform="Dehydrogenate the cis-alkene to alkyne.",
)

_t(
    name="CBS asymmetric reduction",
    named_reaction="CBS reduction",
    category=ReactionCategory.REDUCTION,
    reagents=["BH3*THF", "(S)-CBS catalyst (Corey-Bakshi-Shibata oxazaborolidine)"],
    solvents=["THF", "toluene"],
    catalysts=["(S)-CBS oxazaborolidine"],
    temperature_range=(-40.0, 0.0),
    typical_yield=(80.0, 95.0),
    functional_group_required=["ketone"],
    functional_group_produced=["alcohol"],
    functional_group_incompatible=[],
    mechanism="Chiral oxazaborolidine pre-coordinates the ketone and activates BH3; face-selective hydride delivery via a six-membered transition state gives the alcohol with >95% ee.",
    reverse_transform="Oxidise the chiral alcohol to the ketone.",
)

_t(
    name="Noyori asymmetric hydrogenation (Ru-BINAP)",
    named_reaction="Noyori asymmetric hydrogenation",
    category=ReactionCategory.REDUCTION,
    reagents=["H2"],
    solvents=["MeOH", "iPrOH"],
    catalysts=["RuCl2(BINAP)(dmf)n"],
    temperature_range=(20.0, 50.0),
    typical_yield=(85.0, 99.0),
    functional_group_required=["ketone"],
    functional_group_produced=["alcohol"],
    functional_group_incompatible=[],
    mechanism="Ru-BINAP complex catalyses enantioselective transfer of H2 to the ketone via a six-membered pericyclic transition state involving the N-H of the diamine co-ligand and Ru-H.",
    reverse_transform="Oxidise the chiral alcohol to the ketone.",
)

_t(
    name="Noyori asymmetric transfer hydrogenation",
    named_reaction="Noyori transfer hydrogenation",
    category=ReactionCategory.REDUCTION,
    reagents=["HCO2H/Et3N (5:2 azeotrope)", "or iPrOH/KOtBu"],
    solvents=["iPrOH", "DMF"],
    catalysts=["RuCl(TsDPEN)(p-cymene)"],
    temperature_range=(20.0, 40.0),
    typical_yield=(80.0, 95.0),
    functional_group_required=["ketone", "imine"],
    functional_group_produced=["alcohol", "primary_amine"],
    functional_group_incompatible=[],
    mechanism="Ru-TsDPEN catalyst transfers hydrogen from isopropanol or formic acid to the substrate via an outer-sphere bifunctional mechanism (no direct Ru-substrate bond). Enantioselectivity >95% ee for aryl ketones.",
    reverse_transform="Oxidise the product back to the ketone/imine.",
)

# -----------------------------------------------------------------
#  RADICAL / PHOTOCHEMISTRY  (~5)
# -----------------------------------------------------------------

_t(
    name="Barton decarboxylation",
    named_reaction="Barton decarboxylation",
    category=ReactionCategory.RADICAL,
    reagents=["Barton ester (thiohydroxamic ester)", "hv or heat", "radical trap (Bu3SnH or thiol)"],
    solvents=["benzene", "CCl4", "toluene"],
    temperature_range=(25.0, 80.0),
    typical_yield=(50.0, 80.0),
    functional_group_required=["carboxylic_acid"],
    functional_group_produced=[],
    functional_group_incompatible=[],
    mechanism="Photolysis of the Barton (thiohydroxamic) ester generates an alkyl radical by decarboxylation; the radical is trapped by Bu3SnH (reduction) or another radical acceptor.",
    reverse_transform="Disconnect the C-H (or C-X) bond; restore the carboxylic acid.",
)

_t(
    name="Norrish Type I photocleavage",
    named_reaction="Norrish Type I",
    category=ReactionCategory.RADICAL,
    reagents=["hv (UV light)"],
    solvents=["acetone", "benzene"],
    temperature_range=(20.0, 25.0),
    typical_yield=(30.0, 60.0),
    functional_group_required=["ketone"],
    functional_group_produced=["aldehyde", "alkene"],
    functional_group_incompatible=[],
    mechanism="UV irradiation excites the ketone to the n-pi* state; alpha-cleavage (homolysis of the C-C bond alpha to C=O) generates an acyl radical and an alkyl radical. Decarbonylation or recombination follows.",
    reverse_transform="Recombine the radical fragments as a ketone.",
)

_t(
    name="Norrish Type II photocleavage",
    named_reaction="Norrish Type II",
    category=ReactionCategory.RADICAL,
    reagents=["hv (UV light)"],
    solvents=["acetone", "benzene"],
    temperature_range=(20.0, 25.0),
    typical_yield=(30.0, 60.0),
    functional_group_required=["ketone"],
    functional_group_produced=["alkene", "ketone"],
    functional_group_incompatible=[],
    mechanism="UV-excited ketone undergoes intramolecular gamma-hydrogen abstraction via a six-membered transition state; the resulting 1,4-biradical undergoes retro-[2+2] fragmentation to give an enol (tautomerises to ketone) and an alkene.",
    reverse_transform="Join an alkene and a methyl ketone via [2+2] cyclisation.",
)

_t(
    name="Wohl-Ziegler bromination",
    named_reaction="Wohl-Ziegler bromination",
    category=ReactionCategory.RADICAL,
    reagents=["NBS (N-bromosuccinimide)", "AIBN (or hv, or BPO)"],
    solvents=["CCl4", "DCM"],
    temperature_range=(25.0, 80.0),
    typical_yield=(55.0, 80.0),
    functional_group_required=["alkene"],
    functional_group_produced=["alkyl_halide_br"],
    functional_group_incompatible=[],
    mechanism="Radical chain process: AIBN generates initiating radicals; hydrogen abstraction from the allylic position by Br radical; NBS maintains low [Br2] for selective allylic (not addition) bromination.",
    reverse_transform="Remove the allylic bromide; restore C-H.",
)

_t(
    name="Paternò-Büchi [2+2] photocycloaddition",
    named_reaction="Paternò-Büchi reaction",
    category=ReactionCategory.RADICAL,
    reagents=["hv (UV light)", "carbonyl compound", "alkene"],
    solvents=["acetone", "benzene"],
    temperature_range=(20.0, 25.0),
    typical_yield=(30.0, 65.0),
    functional_group_required=["ketone", "alkene"],
    functional_group_produced=["ether"],
    functional_group_incompatible=[],
    mechanism="Photoexcited carbonyl (n,pi* state) undergoes [2+2] cycloaddition with an alkene to form an oxetane. Proceeds via either a singlet or triplet biradical intermediate.",
    reverse_transform="Retro-[2+2] fragmentation of the oxetane gives carbonyl + alkene.",
)

# -----------------------------------------------------------------
#  ASYMMETRIC CATALYSIS  (~7)
# -----------------------------------------------------------------

_t(
    name="Sharpless asymmetric dihydroxylation",
    named_reaction="Sharpless AD",
    category=ReactionCategory.ADDITION,
    reagents=["OsO4 (cat.)", "K2OsO4*2H2O", "K3Fe(CN)6", "K2CO3", "chiral ligand (DHQD)2PHAL (or (DHQ)2PHAL)"],
    solvents=["t-BuOH/water (1:1)"],
    catalysts=["OsO4", "(DHQD)2PHAL"],
    temperature_range=(0.0, 25.0),
    typical_yield=(75.0, 95.0),
    functional_group_required=["alkene"],
    functional_group_produced=["diol"],
    functional_group_incompatible=[],
    mechanism="Chiral alkaloid ligand wraps around OsO4 creating a chiral pocket; [3+2] cycloaddition to the alkene face gives an osmate ester with high ee; K3Fe(CN)6 re-oxidises Os(VI) and releases the diol.",
    reverse_transform="Convert diol to alkene via periodate cleavage.",
)

_t(
    name="Jacobsen epoxidation (Mn-salen)",
    named_reaction="Jacobsen epoxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["NaOCl (or mCPBA)", "Mn(salen) catalyst"],
    solvents=["DCM"],
    catalysts=["Jacobsen Mn(III)-salen catalyst"],
    temperature_range=(0.0, 25.0),
    typical_yield=(60.0, 90.0),
    functional_group_required=["alkene"],
    functional_group_produced=["epoxide"],
    functional_group_incompatible=[],
    mechanism="Mn(III)-salen catalyst is oxidised to Mn(V)=O by NaOCl; the oxo species transfers oxygen to the cis-disubstituted alkene enantioselectively via a side-on approach dictated by the salen ligand.",
    reverse_transform="Open the epoxide to restore the alkene.",
)

_t(
    name="Shi asymmetric epoxidation",
    named_reaction="Shi epoxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["Oxone (KHSO5)", "fructose-derived ketone catalyst"],
    solvents=["MeCN/DMM/water"],
    catalysts=["Shi fructose ketone catalyst"],
    temperature_range=(0.0, 5.0),
    typical_yield=(60.0, 85.0),
    functional_group_required=["alkene"],
    functional_group_produced=["epoxide"],
    functional_group_incompatible=[],
    mechanism="Oxone oxidises the chiral ketone to a dioxirane in situ; the dioxirane epoxidises the trans-disubstituted alkene enantioselectively; the ketone is regenerated catalytically.",
    reverse_transform="Open the epoxide to restore the alkene.",
)

_t(
    name="Evans aldol reaction",
    named_reaction="Evans aldol",
    category=ReactionCategory.COUPLING,
    reagents=["oxazolidinone auxiliary", "Bu2BOTf", "DIPEA", "aldehyde"],
    solvents=["DCM"],
    temperature_range=(-78.0, 0.0),
    typical_yield=(70.0, 90.0),
    functional_group_required=["aldehyde"],
    functional_group_produced=["alcohol", "amide"],
    functional_group_incompatible=[],
    mechanism="Bu2BOTf generates the Z-boron enolate of the Evans N-acyloxazolidinone; Zimmerman-Traxler transition state controls facial selectivity in the aldol addition to the aldehyde, giving the syn-aldol with >95% de.",
    reverse_transform="Disconnect the beta-hydroxy carbonyl C-C bond; fragments are aldehyde and the oxazolidinone.",
)

_t(
    name="Enzymatic kinetic resolution",
    named_reaction="Enzymatic resolution",
    category=ReactionCategory.MISC,
    reagents=["lipase (CAL-B, PPL, etc.)", "vinyl acetate (acyl donor)"],
    solvents=["MTBE", "hexane", "toluene"],
    catalysts=["Candida antarctica lipase B (CAL-B)"],
    temperature_range=(25.0, 40.0),
    typical_yield=(40.0, 50.0),
    functional_group_required=["alcohol"],
    functional_group_produced=["ester"],
    functional_group_incompatible=[],
    mechanism="Lipase selectively acylates one enantiomer of a racemic alcohol; the unreacted enantiomer and the acylated product are separated. Maximum theoretical yield per enantiomer is 50%.",
    reverse_transform="Hydrolyse the ester to recover the alcohol.",
)

_t(
    name="Sharpless click reaction (CuAAC)",
    named_reaction="CuAAC (click chemistry)",
    category=ReactionCategory.COUPLING,
    reagents=["terminal alkyne", "organic azide", "CuSO4", "sodium ascorbate"],
    solvents=["t-BuOH/water", "DMF/water"],
    catalysts=["CuSO4/sodium ascorbate"],
    temperature_range=(20.0, 60.0),
    typical_yield=(80.0, 99.0),
    functional_group_required=["alkyne"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=[],
    mechanism="Cu(I) forms a copper acetylide from the terminal alkyne; [3+2] cycloaddition with the azide gives regioselectively the 1,4-disubstituted 1,2,3-triazole. Highly reliable and orthogonal.",
    reverse_transform="Disconnect the triazole ring; fragments are terminal alkyne and azide.",
)

_t(
    name="Electrochemical oxidation (anodic)",
    named_reaction="Anodic oxidation",
    category=ReactionCategory.OXIDATION,
    reagents=["electricity (constant current)", "supporting electrolyte (nBu4NBF4)"],
    solvents=["MeCN", "DMF"],
    temperature_range=(20.0, 25.0),
    typical_yield=(40.0, 80.0),
    functional_group_required=["aromatic_ring"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=[],
    mechanism="Anodic single-electron oxidation generates a radical cation; depending on conditions, coupling, methoxylation, or other functionalization occurs. No chemical oxidant required.",
    reverse_transform="Reverse the electrochemical transformation (reduction).",
)

# -----------------------------------------------------------------
#  MODERN / GREEN CHEMISTRY  (~4)
# -----------------------------------------------------------------

_t(
    name="Photocatalytic oxidative coupling",
    named_reaction="Photocatalytic C-H coupling",
    category=ReactionCategory.COUPLING,
    reagents=["photocatalyst (Ru(bpy)3Cl2 or Ir(ppy)3)", "oxidant or O2"],
    solvents=["DMF", "MeCN"],
    catalysts=["Ru(bpy)3Cl2"],
    temperature_range=(20.0, 30.0),
    typical_yield=(30.0, 70.0),
    functional_group_required=["aromatic_ring"],
    functional_group_produced=["aromatic_ring"],
    functional_group_incompatible=[],
    mechanism="Visible light excites the photocatalyst; single-electron transfer generates reactive radical intermediates from substrates; radical coupling or radical addition to heteroarenes gives C-C or C-heteroatom bonds.",
    reverse_transform="Disconnect the newly formed C-C bond; restore the starting materials.",
)

_t(
    name="Flow hydrogenation (H-Cube)",
    named_reaction="Flow hydrogenation",
    category=ReactionCategory.REDUCTION,
    reagents=["H2 (generated in situ from water electrolysis)"],
    solvents=["EtOAc", "MeOH"],
    catalysts=["Pd/C cartridge", "Pt/C cartridge", "Raney Ni cartridge"],
    temperature_range=(25.0, 80.0),
    typical_yield=(85.0, 99.0),
    functional_group_required=["alkene", "nitro", "nitrile"],
    functional_group_produced=["alcohol", "primary_amine"],
    functional_group_incompatible=[],
    mechanism="Continuous-flow reactor generates H2 from water electrolysis and passes substrate solution over a fixed-bed heterogeneous catalyst. Enhanced mass transfer and safety vs. batch hydrogenation.",
    reverse_transform="Dehydrogenation or oxidation of the reduced product.",
)

_t(
    name="Nickel-catalysed reductive cross-coupling",
    named_reaction="Ni reductive cross-coupling",
    category=ReactionCategory.COUPLING,
    reagents=["aryl halide", "alkyl halide", "Zn (or Mn) reductant"],
    solvents=["DMA", "DMF"],
    catalysts=["NiCl2*glyme/dtbbpy"],
    temperature_range=(20.0, 60.0),
    typical_yield=(40.0, 75.0),
    functional_group_required=["alkyl_halide_br"],
    functional_group_produced=[],
    functional_group_incompatible=[],
    mechanism="Ni(0) oxidatively adds to one halide; the second halide undergoes radical capture at the Ni centre; reductive elimination forms the new C-C bond. Mn or Zn reduces Ni(II) back to Ni(0).",
    reverse_transform="Disconnect C(sp3)-C(sp2) bond; both fragments get halides.",
)

_t(
    name="Biocatalytic ketone reduction (KRED)",
    named_reaction="Ketoreductase reduction",
    category=ReactionCategory.REDUCTION,
    reagents=["ketoreductase enzyme (KRED)", "NADPH cofactor", "glucose/GDH (cofactor recycling)"],
    solvents=["phosphate buffer (pH 7)", "iPrOH co-solvent"],
    catalysts=["KRED engineered enzyme"],
    temperature_range=(25.0, 37.0),
    typical_yield=(80.0, 99.0),
    functional_group_required=["ketone"],
    functional_group_produced=["alcohol"],
    functional_group_incompatible=[],
    mechanism="Enzyme binds the ketone substrate in its active site; NADPH delivers a hydride enantioselectively; glucose dehydrogenase (GDH) regenerates NADPH from NADP+ using glucose. >99% ee is routine with directed-evolution-optimised KREDs.",
    reverse_transform="Oxidise the chiral alcohol to the ketone.",
)


# =====================================================================
#  Lookup helpers
# =====================================================================

def lookup_by_category(cat: ReactionCategory) -> list[ReactionTemplate]:
    """Return all templates belonging to *cat*.

    Parameters
    ----------
    cat : ReactionCategory
        The reaction category to filter on.

    Returns
    -------
    list[ReactionTemplate]
        All templates whose ``category`` matches *cat*.
    """
    return [r for r in REACTION_TEMPLATES if r.category == cat]


def lookup_by_name(name: str) -> list[ReactionTemplate]:
    """Case-insensitive substring search on template name or named_reaction.

    Parameters
    ----------
    name : str
        Search substring (case-insensitive).

    Returns
    -------
    list[ReactionTemplate]
        Templates whose ``name`` or ``named_reaction`` contains *name*.
    """
    q = name.lower()
    return [r for r in REACTION_TEMPLATES
            if q in r.name.lower()
            or (r.named_reaction and q in r.named_reaction.lower())]


def lookup_by_functional_group(fg: str) -> list[ReactionTemplate]:
    """Find reactions that require a given functional group on the substrate.

    Parameters
    ----------
    fg : str
        Functional-group label (case-insensitive), e.g. ``"alcohol"``,
        ``"alkene"``, ``"ketone"``.

    Returns
    -------
    list[ReactionTemplate]
        Templates whose ``functional_group_required`` list contains *fg*.
    """
    q = fg.lower()
    return [r for r in REACTION_TEMPLATES
            if q in [f.lower() for f in r.functional_group_required]]


def lookup_by_reagent(reagent: str) -> list[ReactionTemplate]:
    """Find reactions that use a given reagent.

    Performs a case-insensitive substring search across each template's
    ``reagents`` list.

    Parameters
    ----------
    reagent : str
        Reagent name or fragment (case-insensitive), e.g. ``"NaBH4"``,
        ``"Pd"``, ``"mCPBA"``.

    Returns
    -------
    list[ReactionTemplate]
        Templates that reference *reagent* in their reagent list.
    """
    q = reagent.lower()
    return [r for r in REACTION_TEMPLATES
            if any(q in rgt.lower() for rgt in r.reagents)]


# ------------------------------------------------------------------
#  Legacy aliases (backwards compatibility)
# ------------------------------------------------------------------

def find_reactions_for_fg(fg_name: str) -> list[ReactionTemplate]:
    """Find reactions that can transform a given functional group.

    .. deprecated:: Use :func:`lookup_by_functional_group` instead.
    """
    return lookup_by_functional_group(fg_name)


def find_reactions_producing(fg_name: str) -> list[ReactionTemplate]:
    """Find reactions that produce a given functional group."""
    q = fg_name.lower()
    return [r for r in REACTION_TEMPLATES
            if q in [f.lower() for f in r.functional_group_produced]]


def find_reactions_by_category(category: ReactionCategory) -> list[ReactionTemplate]:
    """Return all templates that belong to *category*.

    .. deprecated:: Use :func:`lookup_by_category` instead.
    """
    return lookup_by_category(category)


def find_reactions_by_name(query: str) -> list[ReactionTemplate]:
    """Search templates by name or named-reaction string.

    .. deprecated:: Use :func:`lookup_by_name` instead.
    """
    return lookup_by_name(query)
