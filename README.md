# MolBuilder

A professional-grade molecular engineering toolkit built from scratch in pure Python. MolBuilder spans the full pipeline from atomic theory and molecular modeling through retrosynthetic analysis, process engineering, and industrial scale-up -- without depending on RDKit, OpenBabel, or any external chemistry library.

## What It Does

MolBuilder is a self-contained chemistry platform that covers seven layers of molecular science:

| Layer | Capabilities |
|-------|-------------|
| **Atomic Physics** | Bohr model, quantum numbers, electron configurations (with Aufbau exceptions), Slater's rules for effective nuclear charge, hydrogen-like wavefunctions, orbital probability densities |
| **Chemical Bonding** | Lewis structures with octet/expanded octet support, VSEPR geometry prediction (12+ molecular shapes), covalent bond analysis (polarity, dipole moments, BDE, sigma/pi orbitals) |
| **Molecular Modeling** | Full molecular graph with 3D coordinate generation, conformational analysis (eclipsed/staggered/gauche/anti), torsional energy profiles, Newman projections, cyclohexane chair/boat conformations, R/S and E/Z stereochemistry via CIP priority rules |
| **Biochemistry** | All 20 standard amino acids with L-chirality, peptide bond formation, phi/psi backbone angles, secondary structure templates (alpha helix, beta sheet) |
| **Cheminformatics** | SMILES parser and writer with chirality (`@`/`@@`), E/Z bond stereochemistry (`/`/`\`), bracket atoms (isotopes, charges, H counts), aromatic perception, Morgan canonical ordering |
| **Retrosynthesis** | Beam-search retrosynthetic engine with 91 curated reaction templates, 200+ purchasable starting materials database, scored disconnections, and forward route extraction |
| **Process Engineering** | Reactor selection (batch/CSTR/PFR/microreactor), condition optimization, purification strategy, GHS safety assessment (69 hazard codes, chemical incompatibility detection), cost estimation, and batch-vs-continuous scale-up analysis |

## Installation

**Requirements:** Python >= 3.11, numpy, scipy, matplotlib

```bash
# From the project directory
pip install -e .

# Or install dependencies manually
pip install numpy scipy matplotlib
```

The GUI uses `tkinter`, which is included with most Python distributions. On some Linux systems you may need to install it separately (`sudo apt install python3-tk`).

## Quick Start

### Interactive CLI

```bash
python -m molbuilder
```

This launches an interactive menu with 11 options:

```
[  1] Bohr Atomic Model
[  2] Quantum Mechanical Atom
[  3] Element Data
[  4] Lewis Structures
[  5] VSEPR Molecular Geometry
[  6] Covalent Bonds
[  7] Molecular Conformations
[  8] Amino Acids & Functional Groups
[  9] 3D Molecule Visualization
[ 10] Quantum Orbital Visualization
[ 11] Molecule Builder Wizard
[  a] Run all text demos (1-8)
[  q] Quit
```

Options 1-8 run educational demos covering atomic models through peptide chemistry. Option 9 renders interactive 3D ball-and-stick models. Option 10 visualizes quantum orbitals. Option 11 launches the interactive molecule builder wizard.

You can also run a specific demo directly: `python -m molbuilder 5` runs the VSEPR demo.

### Molecule Builder Wizard

The wizard (option 11) provides five ways to build molecules:

1. **SMILES input** -- Parse any SMILES string into a 3D molecule
2. **Molecular formula** -- Enter formulas like `H2O`, `SF6` for VSEPR-based geometry
3. **Atom-by-atom** -- Interactive loop to place and bond atoms manually
4. **Preset molecules** -- 10 built-in molecules (ethane, benzene, aspirin, caffeine, etc.)
5. **Peptide builder** -- Build peptides from amino acid sequences with secondary structure

After building, an analysis menu offers functional group detection, bond analysis, SMILES generation, file export, retrosynthesis, process engineering, and 3D visualization.

### SMILES-to-Synthesis Pipeline

```bash
python synthesize.py
```

This standalone script takes a SMILES string and production scale, then runs the complete pipeline:

1. Parses SMILES into a 3D molecule
2. Detects functional groups
3. Checks purchasability against a database of ~200 common chemicals
4. Runs retrosynthetic analysis (beam search, max depth 5)
5. Extracts the best forward synthesis route
6. For each step: selects a reactor, optimizes conditions, plans purification
7. Assesses safety (GHS hazards, PPE, incompatibilities, emergency procedures)
8. Estimates costs (materials, labor, equipment, energy, waste, overhead)
9. Analyzes scale-up (batch sizing, cycle time, annual capacity, capital costs)
10. Generates full ASCII reports

### 3D GUI

```bash
python -c "from molbuilder.gui.app import launch; launch()"
```

A tkinter-based graphical editor with:

- **Element palette** -- Click to select H, C, N, O, F, P, S, Cl, Br, I, or enter a custom element
- **Bond tools** -- Single, double, and triple bond modes
- **3D viewport** -- Embedded matplotlib canvas with CPK-colored atoms, rotatable view
- **Sidebar** -- Molecule info, analysis buttons (functional groups, stereochemistry, retrosynthesis, process engineering), file export (XYZ, MOL/SDF, PDB, JSON, SMILES)
- **File menu** -- Open and save molecules in any supported format
- **Build menu** -- Load from SMILES, formula, or preset molecules

## Python API

### Parse and Write SMILES

```python
from molbuilder.smiles import parse, to_smiles

mol = parse("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
print(f"{len(mol.atoms)} atoms, {len(mol.bonds)} bonds")

canonical = to_smiles(mol)  # canonical SMILES output
```

### Detect Functional Groups

```python
from molbuilder.smiles import parse
from molbuilder.reactions.functional_group_detect import detect_functional_groups

mol = parse("CC(=O)O")  # acetic acid
for fg in detect_functional_groups(mol):
    print(f"  {fg.name} at atoms {fg.atoms}")
# Output: carboxylic_acid at atoms [1, 2, 3]
```

21 functional group detectors: alcohol, aldehyde, ketone, carboxylic acid, ester, amide, amine (primary/secondary/tertiary), alkyl halide, alkene, alkyne, ether, thiol, nitrile, nitro, aromatic ring, epoxide, acid chloride, anhydride, sulfoxide, sulfone, imine.

### Retrosynthetic Analysis

```python
from molbuilder.smiles import parse
from molbuilder.reactions.retrosynthesis import retrosynthesis, format_tree
from molbuilder.reactions.synthesis_route import extract_best_route, format_route

mol = parse("CC(=O)OCC")  # ethyl acetate
tree = retrosynthesis(mol, max_depth=5, beam_width=5)
print(format_tree(tree))

route = extract_best_route(tree)
print(format_route(route))
```

The engine uses beam search to work backwards from the target, matching functional groups against 91 reaction templates. Each disconnection is scored (0-100) based on yield, precursor availability, complexity reduction, and strategic bond preference.

### Process Engineering

```python
from molbuilder.process.reactor import select_reactor
from molbuilder.process.conditions import optimize_conditions
from molbuilder.process.safety import assess_safety
from molbuilder.process.costing import estimate_cost
from molbuilder.process.scale_up import analyze_scale_up

# After obtaining a synthesis route:
for step in route.steps:
    reactor = select_reactor(step.template, scale_kg=100.0)
    conditions = optimize_conditions(step.template, scale_kg=100.0)

safety = assess_safety(route.steps)
cost = estimate_cost(route.steps, scale_kg=100.0)
scaleup = analyze_scale_up(route.steps, target_annual_kg=10000.0)
```

### VSEPR Geometry

```python
from molbuilder.bonding.vsepr import VSEPRMolecule

water = VSEPRMolecule("H2O")
print(water.summary())
# Shows: AXE class, electron geometry, molecular geometry, bond angles, 3D coords
```

### Lewis Structures

```python
from molbuilder.bonding.lewis import LewisStructure

lewis = LewisStructure("CO2")
print(lewis.summary())
# Shows: valence electrons, bonds, lone pairs, formal charges
```

### Quantum Atom

```python
from molbuilder.atomic.quantum_atom import QuantumAtom

fe = QuantumAtom(26)  # iron
print(fe.electron_configuration_string())  # 1s2 2s2 2p6 3s2 3p6 3d6 4s2
print(fe.noble_gas_notation())             # [Ar] 3d6 4s2
print(fe.effective_nuclear_charge(4, 0))   # Z_eff for 4s: ~3.75
```

### Build Peptides

```python
from molbuilder.molecule.amino_acids import build_peptide, AminoAcidType

tripeptide = build_peptide(
    [AminoAcidType.ALA, AminoAcidType.GLY, AminoAcidType.PHE],
    phi=-57, psi=-47  # alpha helix
)
print(f"{len(tripeptide.atoms)} atoms")
```

### File I/O

```python
from molbuilder.smiles import parse
from molbuilder.io.xyz import write_xyz, read_xyz
from molbuilder.io.mol_sdf import write_mol, read_mol
from molbuilder.io.pdb import write_pdb, read_pdb
from molbuilder.io.json_io import write_json, read_json

mol = parse("CCO")

write_xyz(mol, "ethanol.xyz")
write_mol(mol, "ethanol.mol")
write_pdb(mol, "ethanol.pdb")
write_json(mol, "ethanol.json")

mol2 = read_xyz("ethanol.xyz")
mol3 = read_mol("ethanol.mol")
```

Supported formats: XYZ, MOL/SDF (V2000), PDB, JSON, SMILES.

### Report Generation

```python
from molbuilder.reports import (
    generate_molecule_report,
    generate_synthesis_report,
    generate_safety_report,
    generate_cost_report,
)

print(generate_molecule_report(mol))
print(generate_synthesis_report(route))
print(generate_safety_report(safety_assessments))
print(generate_cost_report(cost_estimate))
```

Reports are formatted as clean ASCII tables suitable for terminal display or text file output.

## Reaction Knowledge Base

91 curated reaction templates across 14 categories:

| Category | Examples | Count |
|----------|---------|-------|
| Substitution | SN2 (hydroxide, cyanide, azide), SN1, Williamson ether, Mitsunobu | ~10 |
| Elimination | E2, E1 dehydration, Hofmann, Cope | ~6 |
| Addition | HBr, hydroboration-oxidation, epoxidation, dihydroxylation, ozonolysis | ~10 |
| Oxidation | PCC, Jones, Swern, Dess-Martin, KMnO4, Baeyer-Villiger, Sharpless | ~8 |
| Reduction | NaBH4, LiAlH4, DIBAL-H, Birch, Wolff-Kishner, asymmetric hydrogenation | ~10 |
| Coupling | Grignard, Suzuki, Heck, Sonogashira, Wittig, aldol, C-H activation | ~12 |
| Carbonyl | Fischer esterification, Diels-Alder, Claisen, Michael, Robinson annulation | ~10 |
| Protection | Boc, Fmoc, TBS, benzyl, acetonide (install and remove) | ~12 |
| Rearrangement | Cope, pinacol, Curtius, Beckmann | ~4 |
| Polymerization | ROMP | ~2 |
| Miscellaneous | Appel, Staudinger, olefin metathesis | ~5 |

Each template includes reagents, solvents, catalysts, temperature range, yield range, mechanism description, safety notes, and retrosynthetic transform logic.

## Data Coverage

| Data Set | Entries | Source |
|----------|---------|--------|
| Elements | 118 | IUPAC 2021 |
| Covalent radii | 118 | Cordero et al. 2008, Pyykko & Atsumi 2009 |
| Electronegativity | 103 (Pauling) | CRC Handbook |
| Bond dissociation energies | 49 | CRC Handbook, Luo 2007 |
| Torsion barriers | 16 types | OPLS-AA (Jorgensen et al. 1996) |
| Standard bond lengths | 44 | Experimental averages |
| CPK colors | 56 | Corey-Pauling-Koltun convention |
| GHS hazard statements | 69 | GHS Revision 10 (2023) |
| Purchasable materials | ~200 | Common laboratory chemicals with pricing |
| Reagent database | ~100 reagents, ~30 solvents | CAS numbers, MW, hazards, costs |
| Aufbau exceptions | 22 | Known d-block and f-block anomalies |
| Amino acids | 20 | Standard L-amino acids with full sidechain geometry |

## Architecture

```
molbuilder/
  core/            Constants, elements, geometry, bond data
  atomic/          Bohr model, quantum numbers, wavefunctions
  bonding/         Lewis structures, VSEPR, covalent bond analysis
  molecule/        Molecular graph, conformations, stereochemistry,
                   builders, functional groups, amino acids, peptides
  smiles/          SMILES tokenizer, parser, writer
  io/              File I/O (XYZ, MOL/SDF, PDB, JSON, SMILES)
  reactions/       Reaction templates, reagent database, FG detection,
                   retrosynthesis engine, synthesis route planning
  process/         Reactor selection, conditions, purification, safety,
                   costing, scale-up analysis
  reports/         ASCII report generators (molecule, synthesis, safety, cost)
  visualization/   3D molecule rendering, quantum orbital plots
  gui/             Tkinter-based 3D molecule editor
  cli/             Interactive menu, demos, molecule builder wizard
```

~19,900 lines of source code across 81 Python files in 12 subpackages. No external chemistry dependencies -- only numpy, scipy, and matplotlib.

## Testing

```bash
python -m pytest tests/ -q
```

517 tests covering core chemistry data, atomic models, bonding theory, molecular operations, SMILES round-trips, file I/O, reaction knowledge base, process engineering, and edge cases with scientific validation against known experimental values.

## License

MIT License. See [LICENSE](LICENSE) for details.

Copyright (c) 2025 Taylor C. Powell.
