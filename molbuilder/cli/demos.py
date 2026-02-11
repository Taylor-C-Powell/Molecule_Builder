"""Demo functions for each Molecule Builder module.

Each ``demo_*`` function exercises one logical area of the package and
prints a self-contained text report (plus optional interactive
visualisation for the two matplotlib-based demos).

Migrated from legacy/main.py.
"""


# ===================================================================
# Module demos
# ===================================================================

def demo_bohr_model():
    """Bohr atomic model: orbital radii, energies, spectral transitions."""
    from molbuilder.atomic.bohr import from_symbol

    print("=" * 60)
    print("  BOHR ATOMIC MODEL")
    print("=" * 60)
    print()

    # Hydrogen
    H = from_symbol("H")
    print(H.summary())
    print()

    # Balmer series (visible transitions)
    print("  Hydrogen Balmer Series (visible spectrum):")
    print(f"  {'Transition':<12} {'Wavelength':>12} {'Energy':>12}")
    print(f"  {'-' * 40}")
    for ni in range(3, 7):
        wl = H.transition_wavelength_nm(ni, 2)
        e = H.transition_energy(ni, 2)
        print(f"  n={ni} -> n=2   {wl:>10.1f} nm  {e:>10.3f} eV")
    print()

    # Multi-electron atoms
    for sym in ["C", "Na", "Fe"]:
        atom = from_symbol(sym)
        print(f"  {atom.name} (Z={atom.atomic_number}):")
        print(f"    Ground state energy: {atom.energy_level(1):.3f} eV")
        print(f"    Ionisation energy:   {atom.ionization_energy():.3f} eV")
        print(f"    1st orbital radius:  {atom.orbital_radius_pm(1):.1f} pm")
        print()


def demo_quantum_model():
    """Quantum mechanical atom: electron configurations, Slater's rules."""
    from molbuilder.atomic.quantum_atom import from_symbol, aufbau_order

    print("=" * 60)
    print("  QUANTUM MECHANICAL ATOM")
    print("=" * 60)
    print()

    # Aufbau order (first 10 subshells)
    order = aufbau_order()
    labels = {0: "s", 1: "p", 2: "d", 3: "f"}
    print("  Aufbau filling order:")
    print("   ", " -> ".join(
        f"{n}{labels[l]}" for n, l in order[:10]))
    print()

    # Example atoms including Aufbau exceptions
    atoms = [
        ("H",  "simplest atom"),
        ("C",  "basis of organic chemistry"),
        ("Fe", "transition metal"),
        ("Cr", "Aufbau exception: half-filled 3d"),
        ("Cu", "Aufbau exception: filled 3d"),
    ]

    print(f"  {'Symbol':<6} {'Config':<28} {'Spin mult':>10}")
    print(f"  {'-' * 50}")
    for sym, note in atoms:
        atom = from_symbol(sym)
        config = atom.electron_configuration_string()
        mult = atom.spin_multiplicity()
        print(f"  {sym:<6} {config:<28} {mult:>10}")
    print()

    # Detailed view for carbon
    C = from_symbol("C")
    print(C.summary())
    print()


def demo_element_data():
    """Element property data: electronegativity, radii, CPK colours."""
    from molbuilder.core.element_properties import (
        electronegativity, covalent_radius_pm, cpk_color,
        estimated_bond_length_pm, period, can_expand_octet,
    )

    print("=" * 60)
    print("  ELEMENT DATA")
    print("=" * 60)
    print()

    elements = ["H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I", "Xe"]
    print(f"  {'Sym':<4} {'EN':>5} {'Radius':>8} {'Period':>7} {'Expand':>7} {'Color'}")
    print(f"  {'-' * 48}")
    for sym in elements:
        en = electronegativity(sym)
        r = covalent_radius_pm(sym)
        p = period(sym)
        exp = "yes" if can_expand_octet(sym) else "no"
        col = cpk_color(sym)
        print(f"  {sym:<4} {en:>5.2f} {r:>6.0f} pm {p:>6} {exp:>7}  {col}")
    print()

    # Bond lengths
    print("  Estimated bond lengths:")
    bonds = [("C", "H", 1), ("C", "C", 1), ("C", "C", 2),
             ("C", "O", 2), ("N", "N", 3)]
    for a, b, order in bonds:
        sym = {1: "-", 2: "=", 3: "#"}[order]
        bl = estimated_bond_length_pm(a, b, order)
        print(f"    {a}{sym}{b}: {bl:.0f} pm")
    print()


def demo_lewis_structure():
    """Lewis structures: bonding pairs, lone pairs, octets."""
    from molbuilder.bonding.lewis import LewisStructure

    print("=" * 60)
    print("  LEWIS STRUCTURES")
    print("=" * 60)
    print()

    molecules = [
        ("H2O",  0, "bent, 2 lone pairs"),
        ("CO2",  0, "linear, double bonds"),
        ("NH3",  0, "trigonal pyramidal"),
        ("CH4",  0, "tetrahedral"),
        ("BF3",  0, "incomplete octet"),
        ("SF6",  0, "expanded octet"),
        ("PCl5", 0, "expanded octet"),
        ("NH4",  1, "ammonium ion"),
    ]

    for formula, charge, note in molecules:
        lewis = LewisStructure(formula, charge)
        charge_str = f" (charge {charge:+d})" if charge else ""
        print(f"  {formula}{charge_str} -- {note}")
        print(f"    Valence e-: {lewis.total_valence_electrons}")
        print(f"    Central atom: {lewis.central_symbol}")
        print(f"    Bonding pairs: {lewis.bonding_pairs_on_central()}")
        print(f"    Lone pairs: {lewis.lone_pairs_on_central()}")
        print(f"    Steric number: {lewis.steric_number()}")

        bond_strs = []
        for bond in lewis.bonds:
            sa = lewis.atoms[bond.atom_a]
            sb = lewis.atoms[bond.atom_b]
            sym = {1: "-", 2: "=", 3: "#"}.get(bond.order, "?")
            bond_strs.append(f"{sa}{sym}{sb}")
        print(f"    Bonds: {', '.join(bond_strs)}")
        print()


def demo_vsepr_model():
    """VSEPR molecular geometry: AXnEm classification, 3D coordinates."""
    from molbuilder.bonding.vsepr import VSEPRMolecule

    print("=" * 60)
    print("  VSEPR MOLECULAR GEOMETRY")
    print("=" * 60)
    print()

    molecules = [
        ("CO2",  0),   # linear
        ("BF3",  0),   # trigonal planar
        ("CH4",  0),   # tetrahedral
        ("NH3",  0),   # trigonal pyramidal
        ("H2O",  0),   # bent
        ("PCl5", 0),   # trigonal bipyramidal
        ("SF4",  0),   # seesaw
        ("SF6",  0),   # octahedral
        ("XeF4", 0),   # square planar
        ("ClF3", 0),   # T-shaped
        ("IF5",  0),   # square pyramidal
    ]

    print(f"  {'Formula':<8} {'AXE':>6} {'Electron Geom':<24} {'Molecular Geom':<24}")
    print(f"  {'-' * 64}")
    for formula, charge in molecules:
        mol = VSEPRMolecule(formula, charge)
        axe = mol.axe
        print(f"  {formula:<8} {axe.axe_notation:>6} "
              f"{axe.electron_geometry:<24} {axe.molecular_geometry:<24}")
    print()

    # Detailed view for water
    water = VSEPRMolecule("H2O")
    print(water.summary())
    print()


def demo_covalent_bonds():
    """Covalent bond analysis: polarity, orbital composition, BDE."""
    from molbuilder.bonding.covalent import (
        single_bond, double_bond, triple_bond,
        MolecularBondAnalysis,
    )

    print("=" * 60)
    print("  COVALENT BONDS")
    print("=" * 60)
    print()

    # Individual bonds
    bonds = [
        single_bond("H", "H"),
        single_bond("H", "Cl"),
        double_bond("C", "O"),
        triple_bond("N", "N"),
        single_bond("C", "C"),
    ]

    print(f"  {'Bond':<8} {'Order':>5} {'dEN':>6} {'Polarity':<20}"
          f" {'BDE (kJ)':>9} {'Length':>8}")
    print(f"  {'-' * 62}")
    for b in bonds:
        sym = b.order_symbol
        label = f"{b.symbol_a}{sym}{b.symbol_b}"
        bde = f"{b.dissociation_energy_kj:.0f}" if b.dissociation_energy_kj else "---"
        bl = f"{b.bond_length_pm:.0f} pm"
        print(f"  {label:<8} {b.bond_order:>5} {b.delta_en:>6.2f}"
              f" {b.polarity.name:<20} {bde:>9} {bl:>8}")
    print()

    # Sigma/pi composition
    print("  Orbital composition:")
    for b in bonds:
        label = f"{b.symbol_a}{b.order_symbol}{b.symbol_b}"
        parts = ", ".join(f"{oc.orbital_type.name}" for oc in b.orbital_contributions)
        print(f"    {label}: {parts}")
    print()

    # Molecular polarity
    print("  Molecular polarity analysis:")
    formulas = ["H2", "HCl", "CO2", "H2O", "CH4", "NH3", "BF3", "CCl4"]
    print(f"  {'Formula':<10} {'Polarity':<20} {'sigma':>5} {'pi':>5}")
    print(f"  {'-' * 44}")
    for f in formulas:
        ma = MolecularBondAnalysis(f)
        print(f"  {f:<10} {ma.molecular_polarity:<20}"
              f" {ma.total_sigma_bonds:>5} {ma.total_pi_bonds:>5}")
    print()


def demo_molecular_conformations():
    """Molecular conformations: ethane, butane, cyclohexane, stereochemistry."""
    from molbuilder.molecule.builders import (
        build_ethane, build_butane, build_cyclohexane, build_2_butene,
        build_chiral_molecule,
    )
    from molbuilder.molecule.conformations import classify_conformation, scan_torsion
    from molbuilder.molecule.graph import RingConformation

    print("=" * 60)
    print("  MOLECULAR CONFORMATIONS & STEREOCHEMISTRY")
    print("=" * 60)
    print()

    # Ethane conformations
    print("  --- Ethane: Staggered vs Eclipsed ---")
    staggered = build_ethane(60.0)
    eclipsed = build_ethane(0.0)
    e_s = staggered.torsional_energy(0, 1).total_kj_per_mol
    e_e = eclipsed.torsional_energy(0, 1).total_kj_per_mol
    print(f"    Staggered energy: {e_s:.2f} kJ/mol")
    print(f"    Eclipsed energy:  {e_e:.2f} kJ/mol")
    print(f"    Barrier:          {e_e - e_s:.2f} kJ/mol")
    print()

    # Newman projections
    print("  Newman projection (staggered):")
    np_s = staggered.newman_projection(0, 1)
    for idx, sym, ang in np_s.front_substituents:
        print(f"    Front: {sym}[{idx}] at {ang:.1f} deg")
    for idx, sym, ang in np_s.back_substituents:
        print(f"    Back:  {sym}[{idx}] at {ang:.1f} deg")
    print()

    # Butane conformations
    print("  --- Butane conformations ---")
    print(f"  {'Dihedral':>8} {'Energy':>12} {'Type'}")
    print(f"  {'-' * 36}")
    for dih in [180, 60, 0, 120]:
        mol = build_butane(dih)
        energy = mol.torsional_energy(1, 2).total_kj_per_mol
        conf = classify_conformation(dih)
        print(f"  {dih:>6} deg {energy:>8.2f} kJ/mol  {conf.name}")
    print()

    # Butane torsion scan (condensed)
    print("  --- Butane torsion scan (every 30 deg) ---")
    butane = build_butane(180.0)
    scan = scan_torsion(butane, 1, 2, 0, 3, steps=12)
    print(f"  {'Angle':>7} {'kJ/mol':>8}  Bar")
    for angle, energy in scan:
        bar = "#" * max(0, int(energy / 2))
        print(f"  {angle:>7.0f} {energy:>8.2f}  {bar}")
    print()

    # Cyclohexane
    print("  --- Cyclohexane ---")
    chair = build_cyclohexane(RingConformation.CHAIR)
    boat = build_cyclohexane(RingConformation.BOAT)
    dists_chair = [chair.distance(i, (i+1) % 6) for i in range(6)]
    dists_boat = [boat.distance(i, (i+1) % 6) for i in range(6)]
    print(f"    Chair C-C: {min(dists_chair):.3f} - {max(dists_chair):.3f} A")
    print(f"    Boat  C-C: {min(dists_boat):.3f} - {max(dists_boat):.3f} A")
    print()

    # E/Z isomerism
    print("  --- E/Z Isomerism (2-butene) ---")
    cis = build_2_butene(is_cis=True)
    trans = build_2_butene(is_cis=False)
    print(f"    cis-2-butene:   {cis.assign_ez(1, 2).name}")
    print(f"    trans-2-butene: {trans.assign_ez(1, 2).name}")
    print()

    # R/S chirality
    print("  --- R/S Chirality (CHFClBr) ---")
    chiral = build_chiral_molecule(["H", "F", "Cl", "Br"])
    print(f"    Configuration: {chiral.assign_rs(0).name}")
    print()


def demo_amino_acids():
    """Amino acids: all 20, functional groups, peptide bonds, chirality."""
    from molbuilder.molecule.amino_acids import (
        build_amino_acid, form_peptide_bond, build_peptide,
        set_secondary_structure, add_hydroxyl, add_amino,
        add_thiol, AminoAcidType, SecondaryStructure, AMINO_ACID_DATA,
    )
    from molbuilder.molecule.graph import Molecule, Hybridization
    from molbuilder.core.bond_data import SP3_ANGLE

    print("=" * 60)
    print("  AMINO ACIDS & FUNCTIONAL GROUPS")
    print("=" * 60)
    print()

    # Functional groups
    print("  --- Functional Group Builders ---")
    demo_mol = Molecule("demo chain")
    c0 = demo_mol.add_atom("C", [0.0, 0.0, 0.0], Hybridization.SP3)
    c1 = demo_mol.add_atom_bonded("C", c0, hybridization=Hybridization.SP3)
    c2 = demo_mol.add_atom_bonded("C", c1, angle_ref=c0,
                                   bond_angle_deg=SP3_ANGLE,
                                   hybridization=Hybridization.SP3)

    oh = add_hydroxyl(demo_mol, c0)
    nh2 = add_amino(demo_mol, c1, angle_ref=c0, dihedral_deg=60.0)
    sh = add_thiol(demo_mol, c2, angle_ref=c1, dihedral_deg=60.0)

    print(f"    Carbon chain with -OH, -NH2, -SH attached")
    print(f"    Total atoms: {len(demo_mol.atoms)}, bonds: {len(demo_mol.bonds)}")
    print()

    # All 20 amino acids
    print("  --- All 20 Standard Amino Acids ---")
    print(f"  {'Name':<16} {'1-letter':>8} {'3-letter':>8} {'Atoms':>6} {'Bonds':>6}")
    print(f"  {'-' * 50}")

    all_aa = {}
    for aa_type in AminoAcidType:
        data = AMINO_ACID_DATA[aa_type]
        mol = build_amino_acid(aa_type)
        all_aa[aa_type] = mol
        print(f"  {data['name']:<16} {data['code1']:>8} {data['code3']:>8}"
              f" {len(mol.atoms):>6} {len(mol.bonds):>6}")
    print()

    # Chirality check
    print("  --- L-Chirality Verification ---")
    print(f"  {'Name':<16} {'CIP config':>10}")
    print(f"  {'-' * 28}")
    for aa_type, mol in all_aa.items():
        bb = mol.backbone
        data = AMINO_ACID_DATA[aa_type]
        if aa_type == AminoAcidType.GLY:
            config = "achiral"
        else:
            config = mol.assign_rs(bb.CA).name
        print(f"  {data['name']:<16} {config:>10}")
    print()

    # Dipeptide
    print("  --- Dipeptide: Ala-Gly ---")
    ala = build_amino_acid(AminoAcidType.ALA)
    gly = build_amino_acid(AminoAcidType.GLY)
    dipeptide = form_peptide_bond(ala, gly)
    print(f"    Atoms: {len(dipeptide.atoms)}, Bonds: {len(dipeptide.bonds)}")
    if hasattr(dipeptide, 'residues') and len(dipeptide.residues) >= 2:
        c_i = dipeptide.residues[0].C
        n_i = dipeptide.residues[1].N
        print(f"    Peptide bond: {dipeptide.distance(c_i, n_i):.3f} A")
        ca1 = dipeptide.residues[0].CA
        ca2 = dipeptide.residues[1].CA
        omega = dipeptide.dihedral_angle(ca1, c_i, n_i, ca2)
        print(f"    Omega angle:  {omega:.1f} deg (trans)")
    print()

    # Tripeptide with secondary structure
    print("  --- Tripeptide: Ala-Ala-Ala (alpha helix) ---")
    tripeptide = build_peptide(
        [AminoAcidType.ALA, AminoAcidType.ALA, AminoAcidType.ALA])
    if hasattr(tripeptide, 'residues'):
        set_secondary_structure(tripeptide, tripeptide.residues,
                                SecondaryStructure.ALPHA_HELIX)
    print(f"    Atoms: {len(tripeptide.atoms)}, Bonds: {len(tripeptide.bonds)}")
    print(f"    Residues: {len(getattr(tripeptide, 'residues', []))}")
    print(f"    Phi/psi: alpha helix (-57, -47)")
    print()


def demo_visualization():
    """3D molecule visualization (requires matplotlib)."""
    try:
        import matplotlib
    except ImportError:
        print("  matplotlib is not installed. Skipping visualization.")
        return

    from molbuilder.bonding.vsepr import VSEPRMolecule
    from molbuilder.visualization.molecule_viz import visualize_molecule, visualize_gallery

    print("=" * 60)
    print("  3D MOLECULE VISUALIZATION")
    print("=" * 60)
    print()

    choice = input("  Show: [1] Single molecule  [2] Gallery  [3] Skip: ").strip()

    if choice == "1":
        formula = input("  Enter formula (e.g. H2O, CH4, SF6): ").strip()
        if not formula:
            formula = "H2O"
        try:
            mol = VSEPRMolecule(formula)
            print(f"  Rendering {formula}...")
            visualize_molecule(mol)
        except Exception as e:
            print(f"  Error: {e}")

    elif choice == "2":
        formulas = ["H2O", "NH3", "CH4", "CO2", "BF3",
                    "SF6", "PCl5", "XeF4", "ClF3"]
        mols = []
        for f in formulas:
            try:
                mols.append(VSEPRMolecule(f))
            except Exception as e:
                print(f"  Skipping {f}: {e}")
        if mols:
            print(f"  Rendering gallery of {len(mols)} molecules...")
            visualize_gallery(mols)
    else:
        print("  Skipped visualization.")
    print()


def demo_quantum_visualization():
    """Quantum orbital visualization (requires matplotlib)."""
    try:
        import matplotlib
    except ImportError:
        print("  matplotlib is not installed. Skipping visualization.")
        return

    from molbuilder.atomic.quantum_atom import from_symbol as q_from_symbol
    from molbuilder.visualization.quantum_viz import (
        plot_orbital_3d,
        visualize_atom,
    )

    print("=" * 60)
    print("  QUANTUM ORBITAL VISUALIZATION")
    print("=" * 60)
    print()

    choice = input("  Show: [1] Hydrogen 1s  [2] Carbon 2p  "
                   "[3] Full atom overview  [4] Skip: ").strip()

    if choice == "1":
        print("  Rendering hydrogen 1s orbital...")
        plot_orbital_3d(1, 0, 0, Z=1)

    elif choice == "2":
        print("  Rendering carbon 2p orbital...")
        plot_orbital_3d(2, 1, 0, Z=6)

    elif choice == "3":
        element = input("  Enter element symbol (e.g. C, Fe, Na): ").strip()
        if not element:
            element = "C"
        try:
            atom = q_from_symbol(element)
            print(f"  Rendering full overview for {element}...")
            visualize_atom(atom)
        except Exception as e:
            print(f"  Error: {e}")
    else:
        print("  Skipped.")
    print()


def demo_slowmo_interaction():
    """Extreme slow-motion atomic interaction visualization."""
    try:
        import matplotlib
    except ImportError:
        print("  matplotlib is not installed. Skipping visualization.")
        return

    print("=" * 60)
    print("  EXTREME SLOW-MOTION ATOMIC INTERACTIONS")
    print("=" * 60)
    print()

    print("  This demo visualizes atomic interactions at sub-femtosecond")
    print("  time resolution using molecular dynamics simulation.")
    print()

    choice = input(
        "  Show: [1] Ethane MD vibration  [2] SN2 mechanism\n"
        "        [3] Bond formation        [4] Skip: ").strip()

    if choice == "1":
        from molbuilder.molecule.builders import build_ethane
        from molbuilder.visualization.interaction_viz import (
            visualize_md_trajectory, PlaybackConfig,
        )
        print("  Building ethane and running MD at 300 K...")
        mol = build_ethane(60.0)
        config = PlaybackConfig(
            show_electron_density=False,
            slowmo_factor=1e15,
        )
        visualize_md_trajectory(mol, n_steps=500, config=config)

    elif choice == "2":
        from molbuilder.molecule.builders import build_ethane
        from molbuilder.dynamics.mechanisms import sn2_mechanism
        from molbuilder.visualization.interaction_viz import (
            visualize_reaction, PlaybackConfig,
        )
        print("  Setting up SN2 mechanism demonstration...")
        mol = build_ethane(60.0)
        mechanism = sn2_mechanism(
            substrate_C=0, nucleophile=1, leaving_group=2)
        config = PlaybackConfig(
            show_electron_density=True,
            show_electron_flows=True,
        )
        visualize_reaction(
            mol, mechanism,
            n_steps_per_stage=150, config=config)

    elif choice == "3":
        from molbuilder.molecule.graph import Molecule, Hybridization
        from molbuilder.visualization.interaction_viz import (
            visualize_bond_formation, PlaybackConfig,
        )
        print("  Setting up bond formation between two carbon atoms...")
        mol = Molecule("C...C bond formation")
        c0 = mol.add_atom("C", [0.0, 0.0, 0.0], Hybridization.SP3)
        c1 = mol.add_atom("C", [3.0, 0.0, 0.0], Hybridization.SP3)
        config = PlaybackConfig(
            show_electron_density=True,
            show_electron_flows=True,
        )
        visualize_bond_formation(mol, c0, c1, config=config)

    else:
        print("  Skipped.")
    print()
