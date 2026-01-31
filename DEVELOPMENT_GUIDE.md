# Molbuilder Development Guide

Authoritative reference for continuing development after context compactions.
This file should be read at the start of every new session or after context resets.

---

## Project Overview

**molbuilder** is a professional-grade molecular engineering tool restructured from
a flat-file educational chemistry project (~8,000 lines, 13 files) into a proper
Python package (~19,100+ lines, 72 files across 12 subpackages).

- **Pure Python** + numpy + scipy + matplotlib (no RDKit, no OpenBabel)
- **Windows cp1252 compatible** -- no special unicode in source files
- **Python >= 3.11** required (uses `X | Y` union types, `match` statements)
- Entry point: `python -m molbuilder` runs the interactive CLI

---

## Package Structure

```
Molecule_Builder/
  pyproject.toml              # setuptools build config
  requirements.txt            # numpy>=1.24, scipy>=1.10, matplotlib>=3.7
  DEVELOPMENT_GUIDE.md        # THIS FILE
  legacy/                     # Original 13 flat files preserved for reference
  tests/
    __init__.py
  molbuilder/
    __init__.py               # __version__ = "1.0.0"
    __main__.py               # Entry point -> cli.menu.main()
    core/                     # Shared constants, elements, geometry, bond data
    atomic/                   # Bohr model, quantum numbers, wavefunctions
    bonding/                  # Lewis structures, VSEPR, covalent bond analysis
    molecule/                 # Molecule graph, conformations, builders, amino acids
    smiles/                   # SMILES tokenizer, parser, writer
    io/                       # File I/O: XYZ, JSON, MOL/SDF, PDB, SMILES
    reactions/                # Reaction templates, reagent DB, FG detection
    process/                  # [PENDING] Reactor, solvents, costing, safety, scale-up
    visualization/            # Bohr, quantum, and molecule 3D visualizations
    gui/                      # tkinter + matplotlib 3D molecule builder GUI
    reports/                  # [PENDING] Text report generators
    cli/                      # Interactive menu system and demo functions
```

---

## Phase Completion Status

| Phase | Description | Status | Verification |
|-------|-------------|--------|--------------|
| 0 | Infrastructure + core/ | COMPLETE | All core imports work |
| 1 | Migrate existing modules | COMPLETE | `python -m molbuilder all` -- 8 text demos pass |
| 2 | SMILES parser + file I/O | COMPLETE | parse("CCO") works; XYZ/JSON/MOL/PDB/SDF round-trips pass |
| 3 | Reaction knowledge base | COMPLETE | 89 templates, 121 reagents, 32 solvents, 16 FG detectors |
| 4 | Retrosynthetic analysis | COMPLETE | 206 purchasable materials, beam search, forward route extraction |
| 5 | Process engineering | COMPLETE | 7 modules: reactor, solvents, purification, conditions, costing, safety, scale-up |
| 6 | Report generation | COMPLETE | 5 modules: text_formatter, molecule, synthesis, safety, cost reports |
| 7 | CLI wizard | COMPLETE | 5 build flows + 10-option analysis menu, integrated as menu option 11 |
| 8 | 3D GUI | COMPLETE | 6 modules: canvas3d, toolbar, sidebar, dialogs, event_handler, app |
| 9 | Testing + polish | REMAINING | Write formal test suite in tests/ directory |

---

## Key Verification Commands

```bash
# Verify all 8 text demos (Phase 1)
python -m molbuilder all

# Verify SMILES parser (Phase 2)
python -c "from molbuilder.smiles import parse, to_smiles; mol = parse('CCO'); print(to_smiles(mol))"

# Verify file I/O round-trips (Phase 2)
python -c "
from molbuilder.io import write_xyz, read_xyz
from molbuilder.molecule.builders import build_ethane
import tempfile, os
mol = build_ethane()
tf = os.path.join(tempfile.gettempdir(), 'test.xyz')
write_xyz(mol, tf); mol2 = read_xyz(tf)
print(len(mol2.atoms), 'atoms')
"

# Verify functional group detection (Phase 3)
python -c "
from molbuilder.reactions.functional_group_detect import detect_functional_groups
from molbuilder.smiles import parse
mol = parse('CCO')
print([fg.name for fg in detect_functional_groups(mol)])
"

# Verify reaction knowledge base (Phase 3 -- once knowledge_base.py exists)
python -c "
from molbuilder.reactions.knowledge_base import REACTION_TEMPLATES, lookup_by_category
from molbuilder.reactions.reaction_types import ReactionCategory
print(len(REACTION_TEMPLATES), 'templates')
print(len(lookup_by_category(ReactionCategory.OXIDATION)), 'oxidation reactions')
"

# Verify GUI imports (Phase 8 -- headless OK)
python -c "from molbuilder.gui.app import MolBuilderApp; print('GUI OK')"
```

---

## Core Architecture

### The Molecule Class (`molbuilder.molecule.graph.Molecule`)

This is the central data structure. Everything in the system operates on Molecule objects.

```python
class Molecule:
    name: str
    atoms: list[Atom]          # Atom(symbol, position, index, hybridization)
    bonds: list[Bond]          # Bond(atom_i, atom_j, order, rotatable)
    _adj: dict[int, list[int]] # adjacency list (internal)
```

**Key Atom fields:** `atom.symbol` (e.g. "C", "H", "O"), `atom.position` (np.ndarray[3]), `atom.index` (int), `atom.hybridization` (Hybridization enum or None).

**IMPORTANT:** The Atom class uses `.symbol`, NOT `.element`. Any new code must use `atom.symbol`.

**Key Molecule methods:**
- `add_atom(symbol, position, hybridization) -> int` -- returns index
- `add_bond(i, j, order, rotatable) -> Bond`
- `add_atom_bonded(symbol, bonded_to, ...) -> int` -- z-matrix style placement
- `neighbors(idx) -> list[int]` -- indices of bonded atoms
- `get_bond(i, j) -> Bond | None`
- `distance(i, j) -> float`
- `bond_angle(i, j, k) -> float`
- `dihedral_angle(i, j, k, l) -> float`
- `rotate_dihedral(j, k, angle_deg)` -- rotates k-side of bond j-k
- `set_dihedral(i, j, k, l, target_deg)`
- `torsional_energy(j, k) -> StrainEnergy`
- `newman_projection(j, k) -> NewmanProjection`
- `cip_priority(center, nbrs, depth) -> list[int]`
- `assign_rs(center) -> Stereodescriptor`
- `assign_ez(j, k) -> Stereodescriptor`
- `is_in_ring(i, j) -> bool`
- `close_ring(i, j, order)`
- `to_coordinates_dict() -> dict`
- `summary() -> str`

### Enums (all in `molbuilder.molecule.graph`)
- `Hybridization`: SP3, SP2, SP
- `ConformationType`: ECLIPSED, STAGGERED, GAUCHE, ANTI, CUSTOM
- `RingConformation`: CHAIR, BOAT, TWIST_BOAT, HALF_CHAIR, FLAT
- `Stereodescriptor`: R, S, E, Z, NONE

### Bond Data (in `molbuilder.core.bond_data`)
- `STANDARD_BOND_LENGTHS`: dict mapping "(X,Y,order)" tuples to angstroms
- `BDE_TABLE`: bond dissociation energies in kJ/mol
- `TORSION_BARRIERS`: OPLS-AA Fourier coefficients for torsional energy
- `SP3_ANGLE`, `SP2_ANGLE`, `SP_ANGLE`: ideal bond angles in degrees
- `bond_length(sym1, sym2, order) -> float` -- lookup with fallback

---

## Import Patterns

The package uses absolute imports everywhere. Standard patterns:

```python
# Core utilities
from molbuilder.core.constants import PLANCK_CONSTANT, SPEED_OF_LIGHT
from molbuilder.core.elements import ELEMENTS, SYMBOL_TO_Z, from_symbol
from molbuilder.core.element_properties import electronegativity, covalent_radius_pm, cpk_color
from molbuilder.core.geometry import normalize, rotation_matrix, place_atom_zmatrix
from molbuilder.core.bond_data import bond_length, SP3_ANGLE, BDE_TABLE

# Molecule building
from molbuilder.molecule.graph import Molecule, Atom, Bond, Hybridization
from molbuilder.molecule.builders import build_ethane, build_butane, build_cyclohexane
from molbuilder.molecule.conformations import classify_conformation, scan_torsion
from molbuilder.molecule.functional_groups import add_hydroxyl, add_amino_group
from molbuilder.molecule.amino_acids import AminoAcidType, build_peptide

# SMILES
from molbuilder.smiles import parse, to_smiles  # convenience re-exports
from molbuilder.smiles.parser import parse
from molbuilder.smiles.writer import to_smiles

# File I/O (convenience re-exports from io/__init__.py)
from molbuilder.io import write_xyz, read_xyz, write_json, read_json
from molbuilder.io import write_mol, read_mol, write_sdf, read_sdf
from molbuilder.io import write_pdb, read_pdb, write_smiles, read_smiles

# Reactions
from molbuilder.reactions.reaction_types import ReactionCategory, ReactionTemplate
from molbuilder.reactions.reagent_data import REAGENT_DB, SOLVENT_DB
from molbuilder.reactions.functional_group_detect import detect_functional_groups, FunctionalGroup
from molbuilder.reactions.knowledge_base import REACTION_TEMPLATES, lookup_by_category

# Atomic models
from molbuilder.atomic.bohr import BohrAtom
from molbuilder.atomic.quantum_atom import QuantumAtom
from molbuilder.bonding.lewis import LewisStructure, parse_formula
from molbuilder.bonding.vsepr import VSEPRMolecule

# Visualization
from molbuilder.visualization.molecule_viz import visualize_molecule
from molbuilder.visualization.bohr_viz import visualize as visualize_bohr
from molbuilder.visualization.quantum_viz import plot_orbital_3d

# GUI
from molbuilder.gui import launch  # opens tkinter 3D builder window
```

---

## Remaining Work (Phases 4-7, 9)

### Phase 4: Retrosynthetic Analysis

Create two new files in `molbuilder/reactions/`:

**`retrosynthesis.py`** (~500 lines):
- Breadth-first retrosynthetic tree with beam search (top-K=5, max depth=8)
- Algorithm: detect FGs on target -> check purchasability -> lookup matching reactions -> apply reverse transforms -> generate precursor molecules -> score disconnections -> recurse on non-purchasable precursors
- Scoring: strategic bond preference (C-C bonds), symmetry bonus, simplification (fewer rings/stereocenters), precursor availability, reaction yield
- Include a PURCHASABLE_MATERIALS database (~200 common starting materials with name, SMILES, cost_per_kg)
- Key function: `retrosynthesis(mol: Molecule, max_depth=8, beam_width=5) -> RetrosynthesisTree`

**`synthesis_route.py`** (~200 lines):
- `SynthesisStep` dataclass (template, precursors, conditions, expected_yield)
- `SynthesisRoute` dataclass (steps, overall_yield, total_cost_estimate)
- `extract_best_route(tree: RetrosynthesisTree) -> SynthesisRoute` -- convert retro tree to forward plan

Both modules must import Molecule from `molbuilder.molecule.graph` and ReactionTemplate from `molbuilder.reactions.reaction_types`. Use `detect_functional_groups()` from `functional_group_detect.py` and `REACTION_TEMPLATES` / lookup functions from `knowledge_base.py`.

### Phase 5: Process Engineering

Create 7 files in `molbuilder/process/`:

| File | Key exports | Description |
|------|-------------|-------------|
| `reactor.py` | `ReactorType` enum, `ReactorSpec` dataclass, `select_reactor(reaction, scale_kg)` | Decision tree for batch/CSTR/PFR/microreactor selection |
| `solvent_systems.py` | `SolventRecommendation`, `select_solvent(reaction)` | Score solvents by polarity, bp, green chemistry rating |
| `purification.py` | `PurificationMethod` enum, `PurificationStep`, `recommend_purification(product, byproducts, scale)` | Distillation/recrystallization/chromatography/extraction selection |
| `conditions.py` | `ReactionConditions`, `optimize_conditions(reaction, scale)` | Scale-dependent temperature, addition rate, workup adjustments |
| `costing.py` | `CostEstimate`, `estimate_cost(route, scale_kg)` | Sum reagent costs, labor hours, equipment, energy, waste disposal |
| `safety.py` | `SafetyAssessment`, `assess_safety(route)` | GHS hazard lookup per reagent, PPE requirements, emergency procedures |
| `scale_up.py` | `ScaleUpAnalysis`, `analyze_scale_up(route, annual_kg)` | Batch vs continuous, batch sizing, cycle time, capacity |

All process modules consume `SynthesisRoute` from Phase 4 and `REAGENT_DB`/`SOLVENT_DB` from `reagent_data.py`.

### Phase 6: Report Generation

Create 5 files in `molbuilder/reports/`:

| File | Key exports | Description |
|------|-------------|-------------|
| `text_formatter.py` | `ascii_table()`, `section_header()`, `word_wrap()` | ASCII formatting utilities, cp1252 safe |
| `synthesis_report.py` | `generate_synthesis_report(route) -> str` | Full route description with reagents/conditions/yields |
| `safety_report.py` | `generate_safety_report(route) -> str` | GHS hazards per reagent per step, PPE, emergency procedures |
| `cost_report.py` | `generate_cost_report(estimate) -> str` | Materials/labor/equipment/energy/waste breakdown |
| `molecule_report.py` | `generate_molecule_report(mol) -> str` | Formula, MW, atoms, bonds, FGs, stereochemistry |

### Phase 7: CLI Wizard

Create `molbuilder/cli/wizard.py` (~400 lines):

Interactive step-by-step flows:
1. Enter SMILES string -> parse -> analyze
2. Enter formula for simple molecules
3. Atom-by-atom builder (interactive)
4. Choose from preset molecules
5. Build peptide from amino acid sequence

After building: analysis menu (properties, bond analysis, retrosynthesis, process engineering, export to file, visualize 3D).

Integrate into existing `cli/menu.py` as a new menu option.

### Phase 9: Testing + Polish

Create test files in `tests/`:
- `test_core.py` -- constants, elements, geometry, bond_data
- `test_atomic.py` -- BohrAtom, QuantumAtom
- `test_bonding.py` -- Lewis, VSEPR, covalent
- `test_molecule.py` -- Molecule building, conformations, stereochemistry, amino acids
- `test_smiles.py` -- tokenizer, parser, writer round-trips
- `test_io.py` -- file format round-trips
- `test_reactions.py` -- FG detection, knowledge base lookup
- `test_process.py` -- process engineering modules

---

## Known Bug Fixes Applied

1. **`atom.element` vs `atom.symbol`**: The Atom class uses `.symbol`, not `.element`. The `functional_group_detect.py` module was fixed by replacing all `atom.element` with `atom.symbol`.

2. **`mol.bonds_of()` does not exist**: The Molecule class has `mol.neighbors(idx)` returning `list[int]`, not `bonds_of()`. The `_neighbors()` helper in `functional_group_detect.py` was fixed to call `mol.neighbors(idx)` directly.

3. **Import paths in `demos.py`**: Builder functions (`build_ethane`, etc.) are in `molbuilder.molecule.builders`, NOT in `molbuilder.molecule.conformations`. The `Molecule` and `Hybridization` classes are in `molbuilder.molecule.graph`. The `SP3_ANGLE` constant is in `molbuilder.core.bond_data`.

---

## API Gotchas for New Modules

When writing new modules that operate on Molecule objects:

- Use `atom.symbol` (not `.element`) to get element string ("C", "H", "O", etc.)
- Use `mol.neighbors(idx)` to get bonded atom indices (returns `list[int]`)
- Use `mol.get_bond(i, j)` to get a Bond object (returns `Bond | None`)
- Bond order is `bond.order` (int: 1, 2, or 3)
- Bond atoms: `bond.atom_i`, `bond.atom_j` (int indices)
- Atom position: `atom.position` (numpy array of shape (3,))
- Atom index: `atom.index` (int, matches position in `mol.atoms` list)
- Use `mol.add_atom_bonded()` for z-matrix style construction (auto-calculates 3D position)
- Use `mol.close_ring(i, j, order)` to close ring bonds (adds bond without moving atoms)
- The SMILES parser `parse(smiles_str)` returns a Molecule with 3D coordinates and implicit H atoms
- `to_smiles(mol)` generates a canonical SMILES string from a Molecule

---

## Dependency Graph Between Subpackages

```
core/  <-- everything depends on this (no internal dependencies)
  |
  +-- atomic/      (imports core/)
  +-- bonding/     (imports core/)
  +-- molecule/    (imports core/)
  |     |
  |     +-- smiles/  (imports molecule/, core/)
  |     +-- io/      (imports molecule/, smiles/, core/)
  |
  +-- reactions/   (imports molecule/, core/)
  |     |
  |     +-- process/  (imports reactions/, molecule/, core/)
  |
  +-- visualization/  (imports molecule/, atomic/, bonding/, core/)
  +-- gui/            (imports molecule/, smiles/, io/, reactions/, visualization/, core/)
  +-- reports/        (imports reactions/, process/, molecule/, core/)
  +-- cli/            (imports everything above)
```

---

## Data Sizes

| Database | Count | Location |
|----------|-------|----------|
| Elements | 118 | `core/elements.py` ELEMENTS dict |
| Reagents | 121 | `reactions/reagent_data.py` REAGENT_DB |
| Solvents | 32 | `reactions/reagent_data.py` SOLVENT_DB |
| Bond lengths | 27 entries | `core/bond_data.py` STANDARD_BOND_LENGTHS |
| BDE values | 40+ entries | `core/bond_data.py` BDE_TABLE |
| Reaction templates | ~80 | `reactions/knowledge_base.py` REACTION_TEMPLATES |
| Amino acids | 20 | `molecule/amino_acids.py` AminoAcidType enum |
| FG patterns | 16 detectors | `reactions/functional_group_detect.py` |
| Torsion barriers | OPLS-AA params | `core/bond_data.py` TORSION_BARRIERS |

---

## File Encoding

All `.py` files MUST be Windows cp1252 compatible. This means:
- No emojis or unicode mathematical symbols
- Use ASCII approximations: `->` not right-arrow, `>=` not greater-or-equal sign
- Degree symbol: use `deg` or write "degrees" in strings
- Greek letters: spell out (alpha, beta, etc.) or use ASCII approximations

---

## Running the Project

```bash
# From the Molecule_Builder directory:
python -m molbuilder           # Interactive menu
python -m molbuilder all       # Run all 8 text demos
python -m molbuilder 1         # Run demo 1 (Bohr model)
python -m molbuilder 2         # Run demo 2 (Quantum model)
# ... through demo 10

# Launch 3D GUI (requires display):
python -c "from molbuilder.gui import launch; launch()"
```

---

## Plan File

The full approved implementation plan is at:
`C:\Users\taylor\.claude\plans\goofy-coalescing-lemur.md`

It contains the complete migration map, target structure, and detailed specifications
for all 9 phases. Consult it for detailed requirements on any phase.
