"""Interactive step-by-step molecule building wizard.

Provides a guided CLI workflow for constructing molecules via five
different approaches: SMILES input, molecular formula (VSEPR), atom-by-
atom assembly, preset molecules, and peptide/amino-acid sequences.

All output is ASCII-safe (Windows cp1252 compatible).
"""

from __future__ import annotations


# ===================================================================
# Amino acid lookup tables
# ===================================================================

# One-letter code -> three-letter code
_ONE_TO_THREE = {
    "G": "GLY", "A": "ALA", "V": "VAL", "L": "LEU", "I": "ILE",
    "P": "PRO", "F": "PHE", "W": "TRP", "M": "MET", "S": "SER",
    "T": "THR", "C": "CYS", "Y": "TYR", "N": "ASN", "Q": "GLN",
    "D": "ASP", "E": "GLU", "K": "LYS", "R": "ARG", "H": "HIS",
}

# Full name -> three-letter code
_NAME_TO_THREE = {
    "GLYCINE": "GLY", "ALANINE": "ALA", "VALINE": "VAL",
    "LEUCINE": "LEU", "ISOLEUCINE": "ILE", "PROLINE": "PRO",
    "PHENYLALANINE": "PHE", "TRYPTOPHAN": "TRP", "METHIONINE": "MET",
    "SERINE": "SER", "THREONINE": "THR", "CYSTEINE": "CYS",
    "TYROSINE": "TYR", "ASPARAGINE": "ASN", "GLUTAMINE": "GLN",
    "ASPARTATE": "ASP", "GLUTAMATE": "GLU", "LYSINE": "LYS",
    "ARGININE": "ARG", "HISTIDINE": "HIS",
}


# ===================================================================
# Helper utilities
# ===================================================================

def _prompt(message: str, default: str = "") -> str:
    """Prompt the user for input.  Returns stripped text."""
    try:
        val = input(message).strip()
    except (EOFError, KeyboardInterrupt):
        print()
        return default
    return val if val else default


def _prompt_int(message: str, lo: int, hi: int) -> int | None:
    """Prompt for an integer in [lo, hi].  Returns None on bad input."""
    raw = _prompt(message)
    try:
        val = int(raw)
        if lo <= val <= hi:
            return val
        print(f"  Please enter a number between {lo} and {hi}.")
    except ValueError:
        if raw:
            print(f"  Invalid number: {raw}")
    return None


def _molecular_formula(mol) -> str:
    """Compute a simple molecular formula string from a Molecule."""
    from collections import Counter
    counts = Counter(atom.symbol for atom in mol.atoms)
    # Hill system: C first, H second, then alphabetical
    parts = []
    for sym in ("C", "H"):
        if sym in counts:
            parts.append(sym if counts[sym] == 1 else f"{sym}{counts[sym]}")
            del counts[sym]
    for sym in sorted(counts):
        parts.append(sym if counts[sym] == 1 else f"{sym}{counts[sym]}")
    return "".join(parts)


def _molecular_weight(mol) -> float:
    """Compute the molecular weight in g/mol."""
    from molbuilder.core.elements import SYMBOL_TO_Z, ELEMENTS
    total = 0.0
    for atom in mol.atoms:
        z = SYMBOL_TO_Z.get(atom.symbol)
        if z is not None:
            total += ELEMENTS[z][2]
    return total


def _validate_element(symbol: str) -> str | None:
    """Validate and normalise an element symbol.  Returns None if invalid."""
    from molbuilder.core.elements import SYMBOL_TO_Z
    sym = symbol.strip()
    if not sym:
        return None
    if len(sym) == 1:
        sym = sym.upper()
    elif len(sym) == 2:
        sym = sym[0].upper() + sym[1].lower()
    else:
        return None
    if sym not in SYMBOL_TO_Z:
        return None
    return sym


# ===================================================================
# Flow 1: Build from SMILES
# ===================================================================

def _flow_smiles():
    """Prompt for a SMILES string, parse it, and return the Molecule."""
    print()
    print("  Enter a SMILES string (e.g. CCO, c1ccccc1, CC(=O)O).")
    smiles_str = _prompt("  SMILES: ")
    if not smiles_str:
        print("  No input provided.")
        return None

    try:
        from molbuilder.smiles import parse
        mol = parse(smiles_str)
        mol.name = mol.name if mol.name else smiles_str
        print()
        print(f"  Successfully parsed: {smiles_str}")
        print(f"  Atoms: {len(mol.atoms)}    Bonds: {len(mol.bonds)}")
        return mol
    except Exception as exc:
        print(f"  Error parsing SMILES: {exc}")
        return None


# ===================================================================
# Flow 2: Build from molecular formula (VSEPR)
# ===================================================================

def _flow_formula():
    """Prompt for a molecular formula and build via VSEPR."""
    print()
    print("  Enter a molecular formula for a simple molecule")
    print("  (e.g. H2O, CH4, NH3, BF3, SF6).")
    formula = _prompt("  Formula: ")
    if not formula:
        print("  No input provided.")
        return None

    try:
        from molbuilder.bonding.vsepr import VSEPRMolecule
        vsepr = VSEPRMolecule(formula)
        print()
        print(vsepr.summary())

        # Convert to a Molecule graph object for the analysis menu
        coords = vsepr.coordinates
        from molbuilder.molecule.graph import Molecule
        mol = Molecule(formula)
        for sym, pos in coords["atom_positions"]:
            mol.add_atom(sym, pos)
        for i, j, order in coords["bonds"]:
            mol.add_bond(i, j, order=order)
        mol.name = formula
        return mol
    except Exception as exc:
        print(f"  Error building molecule: {exc}")
        return None


# ===================================================================
# Flow 3: Step-by-step atom-by-atom builder
# ===================================================================

def _flow_step_by_step():
    """Interactively build a molecule atom by atom."""
    from molbuilder.molecule.graph import Molecule

    print()
    print("  Step-by-step Molecule Builder")
    print("  " + "-" * 40)
    print()

    # First atom
    sym = _prompt("  Enter element symbol for the first atom (e.g. C): ")
    sym = _validate_element(sym) if sym else None
    if sym is None:
        print("  Invalid element symbol.")
        return None

    mol = Molecule("custom molecule")
    mol.add_atom(sym, [0.0, 0.0, 0.0])
    print(f"  Added {sym} as atom 0.")

    while True:
        print()
        n_atoms = len(mol.atoms)
        n_bonds = len(mol.bonds)
        print(f"  Current molecule: {n_atoms} atom(s), {n_bonds} bond(s)")

        # Show atom list compactly
        atom_list = ", ".join(
            f"{a.symbol}[{a.index}]" for a in mol.atoms
        )
        if len(atom_list) > 60:
            atom_list = atom_list[:57] + "..."
        print(f"  Atoms: {atom_list}")
        print()
        print("  Options:")
        print("    [1] Add atom bonded to existing atom")
        print("    [2] Add free (unconnected) atom")
        print("    [3] Remove last atom")
        print("    [4] Show current structure")
        print("    [5] Done -- go to analysis")
        print("    [0] Cancel -- back to wizard menu")
        print()

        choice = _prompt("  Choice: ")

        if choice == "1":
            _step_add_bonded(mol)
        elif choice == "2":
            _step_add_free(mol)
        elif choice == "3":
            _step_remove_last(mol)
        elif choice == "4":
            print()
            print(mol.summary())
        elif choice == "5":
            if len(mol.atoms) == 0:
                print("  Molecule is empty. Add at least one atom first.")
                continue
            return mol
        elif choice == "0":
            return None
        else:
            print("  Invalid choice. Please enter 0-5.")


def _step_add_bonded(mol):
    """Add an atom bonded to an existing atom."""
    if len(mol.atoms) == 0:
        print("  No atoms yet. Add a free atom first.")
        return

    # Pick parent atom
    n = len(mol.atoms)
    if n == 1:
        parent_idx = 0
        print(f"  Bonding to the only atom: {mol.atoms[0].symbol}[0]")
    else:
        print(f"  Available atoms: 0 - {n - 1}")
        parent_val = _prompt_int(f"  Bond to which atom index? [0-{n-1}]: ", 0, n - 1)
        if parent_val is None:
            return
        parent_idx = parent_val

    # Element
    sym = _prompt("  Element symbol for new atom (e.g. H, C, O): ")
    sym = _validate_element(sym) if sym else None
    if sym is None:
        print("  Invalid element symbol.")
        return

    # Bond order
    order_raw = _prompt("  Bond order (1=single, 2=double, 3=triple) [1]: ", "1")
    try:
        order = int(order_raw)
        if order not in (1, 2, 3):
            print("  Invalid bond order. Using single bond.")
            order = 1
    except ValueError:
        order = 1

    try:
        idx = mol.add_atom_bonded(sym, parent_idx, bond_order=order)
        print(f"  Added {sym}[{idx}] bonded to "
              f"{mol.atoms[parent_idx].symbol}[{parent_idx}] "
              f"(order={order}).")
    except Exception as exc:
        print(f"  Error: {exc}")


def _step_add_free(mol):
    """Add a free (unconnected) atom."""
    sym = _prompt("  Element symbol (e.g. C, N, O): ")
    sym = _validate_element(sym) if sym else None
    if sym is None:
        print("  Invalid element symbol.")
        return

    # Place slightly offset from existing atoms
    offset = len(mol.atoms) * 2.0
    idx = mol.add_atom(sym, [offset, 0.0, 0.0])
    print(f"  Added free atom {sym}[{idx}].")


def _step_remove_last(mol):
    """Remove the last atom (and any bonds to it)."""
    if len(mol.atoms) == 0:
        print("  No atoms to remove.")
        return

    last_idx = len(mol.atoms) - 1
    last_sym = mol.atoms[last_idx].symbol

    # Remove bonds involving the last atom
    mol.bonds = [b for b in mol.bonds
                 if b.atom_i != last_idx and b.atom_j != last_idx]
    # Update adjacency
    for nbr in list(mol._adj.get(last_idx, [])):
        mol._adj[nbr] = [x for x in mol._adj[nbr] if x != last_idx]
    if last_idx in mol._adj:
        del mol._adj[last_idx]
    # Remove the atom
    mol.atoms.pop()
    print(f"  Removed {last_sym}[{last_idx}].")


# ===================================================================
# Flow 4: Preset molecules
# ===================================================================

_PRESETS = [
    ("Ethane (staggered)",    "builder",  "ethane"),
    ("Butane (anti)",         "builder",  "butane"),
    ("Cyclohexane (chair)",   "builder",  "cyclohexane"),
    ("2-Butene (cis/Z)",     "builder",  "2-butene"),
    ("Chiral molecule (CHFClBr)", "builder", "chiral"),
    ("Ethanol",               "smiles",   "CCO"),
    ("Acetic acid",           "smiles",   "CC(=O)O"),
    ("Benzene",               "smiles",   "c1ccccc1"),
    ("Aspirin",               "smiles",   "CC(=O)Oc1ccccc1C(=O)O"),
    ("Caffeine",              "smiles",   "Cn1cnc2c1c(=O)n(c(=O)n2C)C"),
]


def _flow_presets():
    """Let the user pick a preset molecule."""
    print()
    print("  Available preset molecules:")
    print()
    for i, (name, _, _) in enumerate(_PRESETS, 1):
        print(f"    [{i:>2}] {name}")
    print(f"    [ 0] Cancel")
    print()

    choice = _prompt_int("  Select preset: ", 0, len(_PRESETS))
    if choice is None or choice == 0:
        return None

    name, kind, key = _PRESETS[choice - 1]

    try:
        if kind == "builder":
            mol = _build_preset(key)
        else:
            from molbuilder.smiles import parse
            mol = parse(key)
            mol.name = name
        print(f"  Built: {mol.name} ({len(mol.atoms)} atoms, {len(mol.bonds)} bonds)")
        return mol
    except Exception as exc:
        print(f"  Error building {name}: {exc}")
        return None


def _build_preset(key: str):
    """Build a preset molecule by key."""
    from molbuilder.molecule.builders import (
        build_ethane, build_butane, build_cyclohexane,
        build_2_butene, build_chiral_molecule,
    )
    builders = {
        "ethane":       build_ethane,
        "butane":       build_butane,
        "cyclohexane":  build_cyclohexane,
        "2-butene":     build_2_butene,
        "chiral":       build_chiral_molecule,
    }
    return builders[key]()


# ===================================================================
# Flow 5: Peptide builder
# ===================================================================

def _parse_amino_acid_token(token: str):
    """Parse a single amino acid token.  Returns an AminoAcidType or None."""
    from molbuilder.molecule.amino_acids import AminoAcidType

    upper = token.upper()

    # Three-letter code
    try:
        return AminoAcidType[upper]
    except (KeyError, ValueError):
        pass

    # One-letter code
    if len(upper) == 1 and upper in _ONE_TO_THREE:
        return AminoAcidType[_ONE_TO_THREE[upper]]

    # Full name
    if upper in _NAME_TO_THREE:
        return AminoAcidType[_NAME_TO_THREE[upper]]

    return None


def _flow_peptide():
    """Build a peptide from an amino acid sequence."""
    print()
    print("  Peptide Builder")
    print("  " + "-" * 40)
    print()
    print("  Enter amino acid sequence using any of these formats:")
    print("    - Three-letter codes: ALA GLY PHE")
    print("    - One-letter codes:   A G F")
    print("    - Full names:         Alanine Glycine Phenylalanine")
    print("    - Or a mix:           ALA G Phenylalanine")
    print()

    raw = _prompt("  Sequence: ")
    if not raw:
        print("  No input provided.")
        return None

    tokens = raw.split()
    sequence = []
    for tok in tokens:
        aa = _parse_amino_acid_token(tok)
        if aa is None:
            print(f"  Unrecognised amino acid: '{tok}'")
            print("  Valid three-letter codes: " + ", ".join(
                t.name for t in __import__(
                    "molbuilder.molecule.amino_acids",
                    fromlist=["AminoAcidType"]).AminoAcidType))
            return None
        sequence.append(aa)

    if not sequence:
        print("  Empty sequence.")
        return None

    names = [aa.name for aa in sequence]
    print(f"  Parsed sequence: {' - '.join(names)} ({len(sequence)} residues)")

    # Optional secondary structure
    print()
    print("  Set secondary structure (optional):")
    print("    [1] Alpha helix")
    print("    [2] Beta sheet")
    print("    [3] Extended")
    print("    [0] None (default backbone angles)")
    print()
    ss_choice = _prompt("  Secondary structure [0]: ", "0")

    phi_psi = None
    ss_label = ""
    if ss_choice == "1":
        phi_psi = [(-57.0, -47.0)] * len(sequence)
        ss_label = " (alpha helix)"
    elif ss_choice == "2":
        phi_psi = [(-135.0, 135.0)] * len(sequence)
        ss_label = " (beta sheet)"
    elif ss_choice == "3":
        phi_psi = [(-180.0, 180.0)] * len(sequence)
        ss_label = " (extended)"

    try:
        from molbuilder.molecule.amino_acids import build_peptide
        mol = build_peptide(sequence, phi_psi=phi_psi)
        mol.name = "-".join(names) + ss_label
        print(f"  Built peptide: {mol.name}")
        print(f"  Atoms: {len(mol.atoms)}    Bonds: {len(mol.bonds)}")
        return mol
    except Exception as exc:
        print(f"  Error building peptide: {exc}")
        return None


# ===================================================================
# Analysis menu
# ===================================================================

def _analysis_menu(mol):
    """Post-build analysis options for a molecule."""
    while True:
        print()
        print("  " + "=" * 56)
        print("  Analysis Options")
        print("  " + "=" * 56)
        formula = _molecular_formula(mol)
        print(f"  Current: {mol.name} ({len(mol.atoms)} atoms, "
              f"{len(mol.bonds)} bonds)")
        print()
        print("    [ 1] Show molecule summary")
        print("    [ 2] Show molecular formula and weight")
        print("    [ 3] Detect functional groups")
        print("    [ 4] Bond analysis")
        print("    [ 5] Generate SMILES")
        print("    [ 6] Export to file (XYZ/MOL/PDB/JSON)")
        print("    [ 7] Retrosynthetic analysis")
        print("    [ 8] Process engineering analysis")
        print("    [ 9] Visualize 3D (requires display)")
        print("    [10] Build another molecule")
        print("    [ 0] Back to main menu")
        print()

        choice = _prompt("  Select option: ")

        if choice == "0":
            return False  # back to main menu
        elif choice == "1":
            _analysis_summary(mol)
        elif choice == "2":
            _analysis_formula_weight(mol)
        elif choice == "3":
            _analysis_functional_groups(mol)
        elif choice == "4":
            _analysis_bonds(mol)
        elif choice == "5":
            _analysis_smiles(mol)
        elif choice == "6":
            _analysis_export(mol)
        elif choice == "7":
            _analysis_retrosynthesis(mol)
        elif choice == "8":
            _analysis_process(mol)
        elif choice == "9":
            _analysis_visualize(mol)
        elif choice == "10":
            return True  # build another
        else:
            print("  Invalid choice. Please enter 0-10.")


def _analysis_summary(mol):
    """Print the full molecule summary."""
    print()
    print(mol.summary())


def _analysis_formula_weight(mol):
    """Print molecular formula and weight."""
    print()
    formula = _molecular_formula(mol)
    weight = _molecular_weight(mol)
    print(f"  Molecular formula: {formula}")
    print(f"  Molecular weight:  {weight:.3f} g/mol")
    print(f"  Atom count:        {len(mol.atoms)}")
    print(f"  Bond count:        {len(mol.bonds)}")


def _analysis_functional_groups(mol):
    """Detect and display functional groups."""
    print()
    try:
        from molbuilder.reactions.functional_group_detect import detect_functional_groups
        groups = detect_functional_groups(mol)
        if groups:
            print(f"  Functional groups found ({len(groups)}):")
            for fg in groups:
                atoms_str = ", ".join(str(a) for a in fg.atoms)
                print(f"    - {fg.name:<20s}  atoms: [{atoms_str}]")
        else:
            print("  No standard functional groups detected.")
    except Exception as exc:
        print(f"  Error detecting functional groups: {exc}")


def _analysis_bonds(mol):
    """Display bond analysis."""
    print()
    if not mol.bonds:
        print("  No bonds in molecule.")
        return

    print(f"  Bond Analysis ({len(mol.bonds)} bonds):")
    print(f"  {'Atoms':<16} {'Order':>6} {'Length (A)':>11} {'Rotatable':>10}")
    print(f"  {'-' * 46}")
    for bond in mol.bonds:
        sa = mol.atoms[bond.atom_i].symbol
        sb = mol.atoms[bond.atom_j].symbol
        sym = {1: "single", 2: "double", 3: "triple"}.get(bond.order, "?")
        dist = mol.distance(bond.atom_i, bond.atom_j)
        rot = "yes" if bond.rotatable else "no"
        label = f"{sa}[{bond.atom_i}]-{sb}[{bond.atom_j}]"
        print(f"  {label:<16} {sym:>6} {dist:>11.3f} {rot:>10}")

    # Summary counts
    from collections import Counter
    order_counts = Counter(b.order for b in mol.bonds)
    print()
    print("  Bond order summary:")
    for order in sorted(order_counts):
        label = {1: "Single", 2: "Double", 3: "Triple"}.get(order, f"Order {order}")
        print(f"    {label}: {order_counts[order]}")
    rot_count = sum(1 for b in mol.bonds if b.rotatable)
    print(f"    Rotatable bonds: {rot_count}")


def _analysis_smiles(mol):
    """Generate and display a SMILES string."""
    print()
    try:
        from molbuilder.smiles import to_smiles
        smi = to_smiles(mol)
        print(f"  SMILES: {smi}")
    except Exception as exc:
        print(f"  Error generating SMILES: {exc}")


def _analysis_export(mol):
    """Export molecule to a file in the chosen format."""
    print()
    print("  Export formats:")
    print("    [1] XYZ")
    print("    [2] MOL (V2000)")
    print("    [3] PDB")
    print("    [4] JSON")
    print("    [0] Cancel")
    print()

    fmt_choice = _prompt("  Format: ")
    if fmt_choice == "0" or not fmt_choice:
        return

    fmt_map = {
        "1": ("xyz", "xyz"),
        "2": ("mol", "mol"),
        "3": ("pdb", "pdb"),
        "4": ("json", "json"),
    }
    if fmt_choice not in fmt_map:
        print("  Invalid format choice.")
        return

    fmt_name, ext = fmt_map[fmt_choice]

    # Suggest a default filename
    safe_name = mol.name.replace(" ", "_").replace("/", "-")
    safe_name = "".join(c for c in safe_name if c.isalnum() or c in "_-.()")
    if not safe_name:
        safe_name = "molecule"
    default_fn = f"{safe_name}.{ext}"

    filepath = _prompt(f"  Filename [{default_fn}]: ", default_fn)

    try:
        from molbuilder.io import write_xyz, write_mol, write_pdb, write_json
        writers = {
            "xyz":  write_xyz,
            "mol":  write_mol,
            "pdb":  write_pdb,
            "json": write_json,
        }
        writers[fmt_name](mol, filepath)
        print(f"  Exported to: {filepath}")
    except Exception as exc:
        print(f"  Error exporting: {exc}")


def _analysis_retrosynthesis(mol):
    """Run retrosynthetic analysis (if module available)."""
    print()
    try:
        from molbuilder.reactions import (
            detect_functional_groups, lookup_by_functional_group,
        )
        groups = detect_functional_groups(mol)
        if not groups:
            print("  No functional groups detected for retrosynthetic analysis.")
            return

        print("  Retrosynthetic Analysis")
        print("  " + "-" * 40)
        print(f"  Detected functional groups: {len(groups)}")
        for fg in groups:
            print(f"    - {fg.name}")

        print()
        print("  Suggested disconnections / reactions:")
        found_any = False
        for fg in groups:
            templates = lookup_by_functional_group(fg.name)
            for tmpl in templates:
                print(f"    [{fg.name}] {tmpl.name}")
                if tmpl.description:
                    print(f"      {tmpl.description}")
                found_any = True

        if not found_any:
            print("    No reaction templates found for detected groups.")

    except ImportError:
        print("  Retrosynthetic analysis module is not yet available.")
        print("  This feature requires molbuilder.reactions to be fully")
        print("  implemented with retrosynthesis planning capabilities.")
    except Exception as exc:
        print(f"  Error in retrosynthetic analysis: {exc}")


def _analysis_process(mol):
    """Run process engineering analysis (if module available)."""
    print()
    try:
        from molbuilder.process.reactor import ReactorType
        from molbuilder.reactions import (
            detect_functional_groups, lookup_by_functional_group,
        )

        groups = detect_functional_groups(mol)
        print("  Process Engineering Analysis")
        print("  " + "-" * 40)
        formula = _molecular_formula(mol)
        weight = _molecular_weight(mol)
        print(f"  Molecule:   {mol.name}")
        print(f"  Formula:    {formula}")
        print(f"  Mol weight: {weight:.3f} g/mol")
        print(f"  Atoms:      {len(mol.atoms)}")
        print(f"  Bonds:      {len(mol.bonds)}")
        print()

        if groups:
            print("  Reactive functional groups:")
            for fg in groups:
                print(f"    - {fg.name}")
            print()
            print("  Suggested reactor types for synthesis:")
            for rt in ReactorType:
                print(f"    - {rt.name}")
        else:
            print("  No reactive functional groups detected.")
            print("  Process analysis requires identifiable reaction sites.")

    except ImportError:
        print("  Process engineering module is not yet fully available.")
        print("  This feature requires molbuilder.process to be fully")
        print("  implemented with reactor selection and costing.")
    except Exception as exc:
        print(f"  Error in process analysis: {exc}")


def _analysis_visualize(mol):
    """Launch 3D visualisation (requires matplotlib and display)."""
    print()
    try:
        from molbuilder.bonding.vsepr import VSEPRMolecule
        from molbuilder.visualization.molecule_viz import visualize_molecule

        # The visualizer expects a VSEPRMolecule, but our mol has
        # to_coordinates_dict().  Build a thin wrapper.
        coords = mol.to_coordinates_dict()

        class _MolWrapper:
            """Minimal wrapper to satisfy visualize_molecule signature."""
            def __init__(self, name, coordinates):
                self.formula = name
                self.coordinates = coordinates
                # Provide a minimal axe attribute
                self.axe = type("AXE", (), {
                    "molecular_geometry": "custom",
                    "electron_geometry": "custom",
                    "hybridisation": "---",
                    "notation": "---",
                    "ideal_bond_angles": [],
                })()
            def computed_bond_angles(self):
                return []
            def summary(self):
                return f"  {self.formula}"

        wrapper = _MolWrapper(mol.name, coords)
        print("  Launching 3D visualisation...")
        print("  (Close the figure window to return to the menu.)")
        visualize_molecule(wrapper)
    except ImportError as exc:
        print(f"  Visualisation requires matplotlib: {exc}")
        print("  Install with: pip install matplotlib")
    except Exception as exc:
        print(f"  Error during visualisation: {exc}")


# ===================================================================
# Main wizard entry point
# ===================================================================

def wizard_main():
    """Top-level interactive molecule building wizard."""
    while True:
        print()
        print("  " + "=" * 56)
        print("  Molecule Builder Wizard")
        print("  " + "=" * 56)
        print()
        print("    [1] Build from SMILES string")
        print("    [2] Build from molecular formula (simple molecules)")
        print("    [3] Build step-by-step (atom by atom)")
        print("    [4] Choose from preset molecules")
        print("    [5] Build peptide from amino acid sequence")
        print("    [6] Back to main menu")
        print()

        choice = _prompt("  Select option: ")

        mol = None

        if choice == "1":
            mol = _flow_smiles()
        elif choice == "2":
            mol = _flow_formula()
        elif choice == "3":
            mol = _flow_step_by_step()
        elif choice == "4":
            mol = _flow_presets()
        elif choice == "5":
            mol = _flow_peptide()
        elif choice == "6":
            return
        else:
            print("  Invalid choice. Please enter 1-6.")
            continue

        if mol is None:
            continue

        # Enter analysis menu
        build_another = _analysis_menu(mol)
        if not build_another:
            return  # back to main menu
        # else: loop back to wizard menu
