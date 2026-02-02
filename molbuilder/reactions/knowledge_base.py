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
    category=ReactionCategory.CARBONYL,
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
    category=ReactionCategory.CARBONYL,
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
    category=ReactionCategory.CARBONYL,
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

_t(
    name="Baeyer-Villiger oxidation (carbonyl section)",
    named_reaction="Baeyer-Villiger oxidation",
    category=ReactionCategory.CARBONYL,
    reagents=["mCPBA"],
    solvents=["DCM"],
    temperature_range=(0.0, 40.0),
    typical_yield=(60.0, 90.0),
    functional_group_required=["ketone"],
    functional_group_produced=["ester"],
    functional_group_incompatible=["primary_amine"],
    mechanism="Peracid attacks the ketone; Criegee intermediate undergoes 1,2-alkyl migration (migratory aptitude: tert > sec > aryl > pri > methyl) to insert oxygen.",
    reverse_transform="Disconnect the ester oxygen; recombine as ketone.",
)

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
