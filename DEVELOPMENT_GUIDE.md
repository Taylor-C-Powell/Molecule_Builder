# Molbuilder Development Guide

Authoritative reference for continuing development after context compactions.
Read this file at the start of every new session or after context resets.

**Current state:** All 10 phases (0-9) COMPLETE. All 3 sprints COMPLETE.
517 tests passing (6 warnings). ~19,900 source lines + ~4,100 test lines.
21 FG detectors. 91 reaction templates. 69 GHS codes. 49 BDE entries. 118 covalent radii.

**Sprint 1 completed fixes (P0 critical):**
- P0-1: Dipole moment formula corrected in `bonding/covalent.py`
- P0-2: Partial charge direction verified correct (no fix needed)
- P0-3: SMILES stereochemistry support added (chirality `@`/`@@`, E/Z `/`/`\` bonds)
- P0-4: Bracket atom data preserved in writer (isotope, charge, H count)
- P0-5: FG detection now works with implicit H (valence-based H inference)
- P0-6: Retrosynthesis string transforms validated before use
- P0-7: Lewis structure bond promotion logic fixed for expanded-octet atoms

**Sprint 2 completed fixes (P1 high-priority):**
- P1-SCI-1: Slater's rules verified correct (Fe Z_eff(4s)=3.75, Z_eff(3d)=6.25)
- P1-SCI-2: Bohr model limitation warning for multi-electron atoms
- P1-SCI-3: Reaction time estimation formula corrected
- P1-SCI-4: Safety incompatibility detection expanded (14 checks)
- P1-SCI-5: Aufbau exceptions expanded for d-block and f-block (22 entries)
- P1-SCI-6: 5-membered aromatic ring detection (furan, pyrrole, thiophene)
- P1-SCI-7: Cyclohexane boat conformation puckering fixed
- P1-SCI-8: Replaced hardcoded vacuum permittivity with VACUUM_PERMITTIVITY
- P1-CODE-1: Input validation on all process/ and reports/ public APIs
- P1-CODE-2: Centralized `normalize_reagent_name()` with alias dictionary
- P1-CODE-3: File I/O bounds checking (XYZ, MOL, PDB, SMILES)
- P1-CODE-4: Named tolerance constants in geometry.py
- P1-CODE-5: Convenience re-exports in all __init__.py files
- P1-CODE-6: Unified solvent ratio constant (DEFAULT_SOLVENT_L_PER_KG = 7.0)
- P1-CODE-7: Duplicate reaction templates removed (91 unique)

**Sprint 3 completed improvements (P2 medium):**
- P2-ALGO-1: O(1) name->Z reverse mapping in elements.py
- P2-ALGO-2: LRU cache on most_probable_radius() (wavefunctions.py)
- P2-ALGO-3: Mutable visited set in retrosynthesis (no copying)
- P2-ALGO-4: Zero-distance guard on coulombs_law()
- P2-ROBUST-1: Warnings for silent covalent radius fallback
- P2-ROBUST-2: Parameterized emergency contact numbers in safety_report.py
- P2-ROBUST-3: Steric clash detection method on Molecule class
- P2-DATA-1: Covalent radii extended to all 118 elements
- P2-DATA-2: Expanded BDE table (49 entries total)
- P2-DATA-3: Expanded torsion barriers (sp2-sp3, sp2-sp2 parameters)
- P2-DATA-4: 5 new FG detectors (acid chloride, anhydride, sulfoxide, sulfone, imine)
- P2-DATA-5: 7 new reaction templates (cross-metathesis, C-H activation, ROMP, etc.)
- P2-DATA-6: GHS hazard statements expanded to 69 codes
- P2-TEST-1: Edge case tests (invalid SMILES, empty molecules, boundary conditions)
- P2-TEST-2: Scientific validation tests (Bohr model, bond lengths, VSEPR angles)
- P2-TYPE-1: Type hints added to geometry.py, reports/, quantum_atom.py
- P2-TYPE-2: Protocol classes (SynthesisStepLike, SynthesisRouteLike, etc.)
- Bonus: scipy.special.sph_harm -> sph_harm_y migration (DeprecationWarning fix)

All P0, P1, and P2 items from the development roadmap are COMPLETE.

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Package Structure](#package-structure)
3. [Core Architecture Reference](#core-architecture-reference)
4. [Critical Bugs (P0)](#critical-bugs-p0)
5. [High-Priority Issues (P1)](#high-priority-issues-p1)
6. [Medium-Priority Improvements (P2)](#medium-priority-improvements-p2)
7. [Implementation Roadmap](#implementation-roadmap)
8. [Verification Commands](#verification-commands)
9. [Import Patterns](#import-patterns)
10. [API Gotchas](#api-gotchas)
11. [Encoding and Style Rules](#encoding-and-style-rules)

---

## Project Overview

**molbuilder** is a professional-grade molecular engineering tool. ~19,900 source
lines + ~4,100 test lines across 81 Python files in 12 subpackages.

- **Pure Python** + numpy + scipy + matplotlib (no RDKit, no OpenBabel)
- **Windows cp1252 compatible** -- no special unicode in source files
- **Python >= 3.11** required (uses `X | Y` union types, `match` statements)
- Entry points:
  - `python -m molbuilder` -- interactive CLI menu (11 options)
  - `python synthesize.py` -- SMILES-to-synthesis pipeline
  - `python -m molbuilder.gui` -- 3D GUI molecule builder

---

## Package Structure

```
Molecule_Builder/
  pyproject.toml
  requirements.txt            # numpy>=1.24, scipy>=1.10, matplotlib>=3.7
  synthesize.py               # Standalone synthesis planner script
  DEVELOPMENT_GUIDE.md        # THIS FILE
  legacy/                     # Original 13 flat files preserved
  tests/
    __init__.py
    test_core.py              # 105 tests
    test_atomic.py            # 110 tests
    test_bonding.py           # 72 tests
    test_molecule.py          # 49 tests
    test_smiles.py            # 30 tests
    test_io.py                # 12 tests
    test_reactions.py         # 23 tests
    test_process.py           # 45 tests
    test_edge_cases.py        # 71 tests (edge cases + scientific validation)
  molbuilder/
    __init__.py               # __version__ = "1.0.0"
    __main__.py               # Entry point -> cli.menu.main()
    core/                     # Shared constants, elements, geometry, bond data
      constants.py            # Physical constants, coulombs_law()
      elements.py             # ELEMENTS dict (118), SYMBOL_TO_Z, from_symbol/from_name
      element_properties.py   # Electronegativity, covalent radii, CPK colors
      geometry.py             # normalize, rotation_matrix, place_atom_zmatrix, etc.
      bond_data.py            # STANDARD_BOND_LENGTHS, BDE_TABLE, TORSION_BARRIERS
    atomic/                   # Bohr model, quantum numbers, wavefunctions
      bohr.py                 # BohrAtom class
      quantum_numbers.py      # QuantumState, Subshell, aufbau_order
      quantum_atom.py         # QuantumAtom, slater_zeff
      wavefunctions.py        # Radial/angular wavefunctions, probability densities
    bonding/                  # Lewis structures, VSEPR, covalent bond analysis
      lewis.py                # LewisStructure, parse_formula
      vsepr.py                # VSEPRMolecule, AXEClassification
      covalent.py             # CovalentBond, MolecularBondAnalysis
    molecule/                 # Molecule graph, conformations, builders, amino acids
      graph.py                # Molecule, Atom, Bond, enums (668 lines, central class)
      conformations.py        # classify_conformation, scan_torsion
      stereochemistry.py      # Re-export of Stereodescriptor
      builders.py             # build_ethane, build_butane, build_cyclohexane, etc.
      functional_groups.py    # add_hydroxyl, add_amino, add_phenyl_ring, etc.
      amino_acids.py          # AminoAcidType, build_amino_acid, build_peptide (920 lines)
      peptides.py             # Re-exports from amino_acids.py
    smiles/                   # SMILES tokenizer, parser, writer
      tokenizer.py            # tokenize() -> list[Token]
      parser.py               # parse(smiles) -> Molecule
      writer.py               # to_smiles(mol) -> str
    io/                       # File I/O: XYZ, JSON, MOL/SDF, PDB, SMILES
      xyz.py                  # to/from XYZ format
      json_io.py              # to/from JSON format
      mol_sdf.py              # to/from MOL V2000 / SDF format
      pdb.py                  # to/from PDB format
      smiles_io.py            # SMILES file read/write
    reactions/                # Reaction templates, reagent DB, FG detection
      reaction_types.py       # ReactionCategory enum, ReactionTemplate dataclass
      reagent_data.py         # REAGENT_DB (~100), SOLVENT_DB (~30)
      knowledge_base.py       # 91 reaction templates, lookup functions
      functional_group_detect.py  # 21 FG detectors
      retrosynthesis.py       # Beam search retrosynthesis engine
      synthesis_route.py      # SynthesisStep, SynthesisRoute, extract_best_route
    process/                  # Reactor, solvents, costing, safety, scale-up
      reactor.py              # ReactorType enum, ReactorSpec, select_reactor()
      solvent_systems.py      # SolventRecommendation, select_solvent()
      purification.py         # PurificationMethod enum, recommend_purification()
      conditions.py           # ReactionConditions, optimize_conditions()
      costing.py              # CostBreakdown, CostEstimate, estimate_cost()
      safety.py               # GHS hazards, SafetyAssessment, assess_safety()
      scale_up.py             # ScaleUpAnalysis, analyze_scale_up()
    visualization/            # Bohr, quantum, molecule 3D visualizations
    gui/                      # tkinter + matplotlib 3D molecule builder
    reports/                  # Text report generators
      text_formatter.py       # ascii_table, section_header, word_wrap, etc.
      molecule_report.py      # generate_molecule_report()
      synthesis_report.py     # generate_synthesis_report()
      safety_report.py        # generate_safety_report()
      cost_report.py          # generate_cost_report()
    cli/                      # Interactive menu and wizard
      menu.py                 # Main menu system (11 options)
      demos.py                # Demo functions
      wizard.py               # Step-by-step molecule building wizard
```

---

## Core Architecture Reference

### The Molecule Class (`molbuilder.molecule.graph.Molecule`)

Central data structure. Everything operates on Molecule objects.

```python
class Molecule:
    name: str
    atoms: list[Atom]          # Atom(symbol, position, index, hybridization)
    bonds: list[Bond]          # Bond(atom_i, atom_j, order, rotatable)
    _adj: dict[int, list[int]] # adjacency list (internal)
```

**Key Atom fields:** `atom.symbol` (e.g. "C", "H", "O"), `atom.position`
(np.ndarray[3]), `atom.index` (int), `atom.hybridization` (Hybridization | None).

**CRITICAL:** Atom uses `.symbol`, NOT `.element`. All code must use `atom.symbol`.

**Key Molecule methods:**
- `add_atom(symbol, position, hybridization) -> int`
- `add_bond(i, j, order, rotatable) -> Bond`
- `add_atom_bonded(symbol, bonded_to, ...) -> int` (z-matrix placement)
- `neighbors(idx) -> list[int]`
- `get_bond(i, j) -> Bond | None`
- `distance(i, j) -> float`
- `bond_angle(i, j, k) -> float`
- `dihedral_angle(i, j, k, l) -> float`
- `rotate_dihedral(j, k, angle_deg)`
- `assign_rs(center) -> Stereodescriptor`
- `assign_ez(j, k) -> Stereodescriptor`
- `is_in_ring(i, j) -> bool`
- `close_ring(i, j, order)`

### Enums (in `molbuilder.molecule.graph`)
- `Hybridization`: SP3, SP2, SP
- `ConformationType`: ECLIPSED, STAGGERED, GAUCHE, ANTI, CUSTOM
- `RingConformation`: CHAIR, BOAT, TWIST_BOAT, HALF_CHAIR, FLAT
- `Stereodescriptor`: R, S, E, Z, NONE

---

## Critical Bugs (P0)

These are scientifically incorrect or functionally broken. Fix before any other work.

### P0-1: Dipole Moment Formula Dimensionally Incorrect

**File:** `molbuilder/bonding/covalent.py` lines 229-240
**Impact:** Produces incorrect dipole moment values by orders of magnitude.

**Current (wrong):**
```python
mu = delta_en * bond_length_angstrom * DEBYE_PER_E_ANGSTROM * frac_ionic
```

**Correct formula:**
```python
# mu = partial_charge(e) * bond_length(Angstrom) * 4.8032(D/e*A)
# where partial_charge = percent_ionic_character / 100
mu = (percent_ionic_character / 100.0) * bond_length_angstrom * DEBYE_PER_E_ANGSTROM
```

The electronegativity difference (`delta_en`) is dimensionless and must NOT appear
in the multiplication. The partial charge in electron units comes from the percent
ionic character, not from delta_EN directly.

**Reference:** https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Physical_Properties_of_Matter/Atomic_and_Molecular_Properties/Dipole_Moments

**Verify after fix:**
```python
from molbuilder.bonding.covalent import single_bond
b = single_bond("H", "F")
# HF experimental: ~1.82 D. Should be in range 1.5-2.5 D, not 0.01 or 50.
assert 0.5 < b.dipole_moment_debye < 5.0
```

### P0-2: Verify Partial Charge Assignment Direction

**File:** `molbuilder/bonding/covalent.py` lines 160-171
**Impact:** If inverted, delta+/delta- labels are backwards for every polar bond.

**Expected behavior:**
- The LESS electronegative atom should be `partial_positive` (delta+)
- The MORE electronegative atom should be `partial_negative` (delta-)

**Verify by reading the exact logic:**
```python
from molbuilder.bonding.covalent import single_bond
b = single_bond("C", "O")
# O is more electronegative (3.44) than C (2.55)
assert b.partial_positive == "C"   # C should be delta+
assert b.partial_negative == "O"   # O should be delta-
```

If the test passes, no fix needed. If it fails, swap the logic.

### P0-3: SMILES Stereochemistry Support

**Files:**
- `molbuilder/smiles/tokenizer.py` lines 136-140 (parses @ @@ but discards)
- `molbuilder/smiles/parser.py` (never stores chirality)
- `molbuilder/smiles/writer.py` (never outputs chirality or bracket data)

**Current behavior:** `@` and `@@` are silently discarded. `/` and `\` for E/Z
double bonds are not even tokenized. Bracket atom data (isotopes, charges, H counts)
is parsed by the tokenizer but lost in the writer.

**Required changes:**

1. **Token dataclass** -- add `chirality: str | None` field (values: None, "@", "@@")
2. **Tokenizer** -- store chirality in Token instead of discarding
3. **Atom dataclass** (in `molecule/graph.py`) -- add `chirality: str | None` field
4. **Parser** -- propagate chirality from Token to Atom
5. **Writer** -- output `@`/`@@` for chiral atoms, output bracket atom data
   (isotope, charge, hcount) when present
6. **Tokenizer** -- add `/` and `\` bond stereochemistry tokens
7. **Parser** -- interpret `/`/`\` for E/Z assignment on double bonds
8. **Writer** -- output `/`/`\` for E/Z stereo bonds

**Scope note:** Steps 1-5 are essential. Steps 6-8 are desirable but can be
deferred. At minimum, chirality must survive a parse -> write round-trip.

**Reference:** http://opensmiles.org/opensmiles.html

**Verify after fix:**
```python
from molbuilder.smiles import parse, to_smiles
# Round-trip should preserve chirality
smi = "[C@@H](F)(Cl)Br"
mol = parse(smi)
out = to_smiles(mol)
assert "@" in out  # chirality preserved
```

### P0-4: SMILES Writer Drops Bracket Atom Data

**File:** `molbuilder/smiles/writer.py` lines 186-193

**Current behavior:** Writer only outputs element symbol for bracket atoms.
Isotopes (`[13C]`), charges (`[NH4+]`), and explicit H counts (`[NH2]`) are lost.

**Required:** When writing bracket atoms, include all stored data:
```
[<isotope><symbol><chirality><hcount><charge>]
```
Example: `[13C@@H]`, `[NH4+]`, `[O-]`

**This overlaps with P0-3.** Implement bracket atom output as part of the
stereochemistry work.

### P0-5: Functional Group Detection Assumes Explicit Hydrogens

**File:** `molbuilder/reactions/functional_group_detect.py` lines 267-303

**Bug:** Amine, thiol, and alcohol detection count explicit H atoms via
`_count_element_neighbors(mol, idx, "H")`. If the Molecule was built via SMILES
parser (which adds explicit H), this works. If built from PDB/XYZ (H may be
implicit or missing), detection silently fails.

**Fix options (choose one):**
1. **Valence-based** (preferred): Compute expected H count from valence rules
   instead of counting explicit H neighbors. Formula:
   `implicit_h = standard_valence - sum(bond_orders) - abs(charge)`
2. **Document & validate**: Add a check at the top of `detect_functional_groups()`
   that verifies the molecule has explicit H atoms, and raise a warning if not.

**Affected detectors:** `_detect_alcohols`, `_detect_amines`, `_detect_thiols`,
`_detect_carboxylic_acids`, `_detect_aldehydes`.

**Verify after fix:**
```python
from molbuilder.smiles import parse
from molbuilder.reactions.functional_group_detect import detect_functional_groups
mol = parse("CCO")  # ethanol
fgs = detect_functional_groups(mol)
assert any(fg.name == "alcohol" for fg in fgs)
```

### P0-6: Retrosynthesis String-Replace Transforms Are Fragile

**File:** `molbuilder/reactions/retrosynthesis.py` lines 755-768

**Bug:** Functions like `_replace_oh_with_carbonyl()` do naive string replacement
on SMILES (e.g., `"CO"` -> `"C=O"`). This corrupts molecules where the substring
appears in different context. Example: `"COCO"` (dimethyl ether) becomes `"C=OCO"`.

**Fix:** Replace string-based SMILES transforms with graph-based transforms:
1. Parse the SMILES to a Molecule
2. Identify the target atoms/bonds by local neighborhood (not string position)
3. Modify the Molecule graph (change bond orders, add/remove atoms)
4. Convert back to SMILES via `to_smiles()`

This is a significant refactor. As intermediate fix, add guards:
- Only apply string transforms if the SMILES length decreases (simpler molecule)
- Validate result parses correctly
- Validate result is chemically sensible (correct valence)

### P0-7: Lewis Structure Multiple Bond Promotion Backwards

**File:** `molbuilder/bonding/lewis.py` lines 244-263

**Bug:** The condition `if not can_expand_octet()` prevents bond promotion for
atoms that CAN expand their octet. This is backwards for cases like SO2 where
sulfur needs double bonds.

**Fix:** Change the promotion condition to check whether the central atom's
current electron count is below its target, regardless of expansion ability:
```python
# BEFORE (wrong):
if not can_expand_octet(self.central_symbol):
    # promote bonds...

# AFTER (correct):
electrons_on_central = self._electrons_around(self.central_index)
target = target_electrons(self.central_symbol)
if electrons_on_central < target:
    # promote bonds...
```

**Verify after fix:**
```python
from molbuilder.bonding.lewis import LewisStructure
ls = LewisStructure("SO2")
# S=O double bonds expected
assert any(b.order == 2 for b in ls.bonds)
```

---

## High-Priority Issues (P1)

### P1-SCI: Scientific Accuracy

#### P1-SCI-1: Slater's Rules Implementation

**File:** `molbuilder/atomic/quantum_atom.py` lines 29-83

**Issue:** The grouping and shielding logic is non-standard. Per Slater's rules:
- Electrons are grouped: (1s)(2s,2p)(3s,3p)(3d)(4s,4p)(4d)(4f)...
- s/p electrons in same group: 0.35 each (0.30 for 1s pair)
- For s/p: (n-1) group shields 0.85, all deeper shield 1.00
- For d/f: ALL inner electrons shield 1.00

**Verification:** Iron (Z=26) should give Z_eff(4s) = 3.75, Z_eff(3d) = 6.25.

**Reference:** https://en.wikipedia.org/wiki/Slater's_rules

#### P1-SCI-2: Bohr Model Limitations Warning

**File:** `molbuilder/atomic/bohr.py`

**Issue:** Bohr model is exact only for hydrogen-like atoms (H, He+, Li2+).
For multi-electron atoms, orbital radius, energy levels, and ionization energy
are poor approximations. No warnings are emitted.

**Fix:** Add a `_BOHR_WARNING` constant and emit it in `summary()` and
`ionization_energy()` when `atomic_number > 1 and charge < atomic_number - 1`.

#### P1-SCI-3: Reaction Time Estimation Inverted

**File:** `molbuilder/process/conditions.py` line 117

**Issue:** Formula `1.0 + (100 - yield) * 0.05` makes lower-yield reactions take
longer. In practice, achieving higher conversion requires longer residence time.

**Fix:** `base_hours = 0.5 + (yield_hi / 100.0) * 2.0` (higher yield = longer).
Or use a more nuanced model based on reaction category.

#### P1-SCI-4: Safety Incompatibility Detection Severely Limited

**File:** `molbuilder/process/safety.py` lines 291-325

**Issue:** Only 5 incompatibility pairs checked. Missing critical pairs:
- Hypochlorite + acids -> Cl2 gas
- Permanganate + organics -> fire/explosion
- Alkali metals + water -> violent reaction
- Concentrated acids + concentrated bases -> exothermic
- Nitrates + organics -> explosion risk
- Chlorine/halogens + ammonia -> toxic gases
- Peroxides + metals -> catalytic decomposition
- Strong oxidizers + flammable solvents

**Fix:** Expand `_determine_incompatibilities()` with at least 15 additional pairs.
Each pair needs: reagent keywords (set), hazard description (str), severity (str).

#### P1-SCI-5: Aufbau Exceptions Incomplete

**File:** `molbuilder/atomic/quantum_numbers.py` lines 150-175

**Issue:** Only 13 elements listed. Missing many d-block exceptions and all
f-block exceptions. Lanthanides/actinides will have incorrect configurations.

**Fix:** Add at minimum: Tc(43), Ru(44), Rh(45), Pd(46), Ag(47), Pr(59), Nd(60),
Pm(61), Sm(62), Eu(63), Gd(64), Tb(65), Dy(66), Ho(67), Er(68), Tm(69), Yb(70),
Lu(71), Ac(89), Th(90), Pa(91), U(92), Np(93), Pu(94), Am(95), Cm(96).

Reference data: NIST Atomic Spectra Database.

#### P1-SCI-6: 5-Membered Aromatic Ring Detection Missing

**File:** `molbuilder/reactions/functional_group_detect.py` line 446

**Issue:** `_detect_aromatic_rings()` only searches for 6-membered rings.
Furan, thiophene, pyrrole, imidazole, pyrazole are not detected.

**Fix:** Add `_find_rings_of_size(mol, start, 5)` call alongside the existing
6-membered ring search. Apply aromaticity check (all sp2/aromatic atoms, planar).

#### P1-SCI-7: Cyclohexane Boat Conformation

**File:** `molbuilder/molecule/builders.py` line 133

**Issue:** Puckering pattern `[d, -d*0.5, -d, d, -d, -d*0.5]` is non-standard.

**Fix:** Use symmetric boat: `[d, -d, 0, d, -d, 0]` or `[d, d, -d, d, d, -d]`
depending on which atoms are "flagpole" vs "bowsprit".

#### P1-SCI-8: Hardcoded Vacuum Permittivity in bohr.py

**File:** `molbuilder/atomic/bohr.py` line 126

**Issue:** Uses hardcoded `8.8541878128e-12` instead of importing
`VACUUM_PERMITTIVITY` from `constants.py`.

**Fix:** Replace with `from molbuilder.core.constants import VACUUM_PERMITTIVITY`.

### P1-CODE: Code Quality

#### P1-CODE-1: No Input Validation on Process/Report Public APIs

**Files:** All `process/*.py` and `reports/*.py` public functions

**Issue:** Functions accept `list[Any]` via duck typing with no validation
that `.template` attribute exists. `TypeError` deep in the stack is the only
feedback on bad input.

**Fix:** Add validation at the top of each public function:
```python
def estimate_cost(steps, scale_kg):
    if not steps:
        return CostEstimate(...)  # empty result
    for i, step in enumerate(steps):
        if not hasattr(step, 'template'):
            raise TypeError(
                f"Step {i} must have a 'template' attribute, "
                f"got {type(step).__name__}"
            )
```

Apply to: `estimate_cost()`, `assess_safety()`, `analyze_scale_up()`,
`optimize_conditions()`, `select_reactor()`, `recommend_purification()`,
all `generate_*_report()` functions.

#### P1-CODE-2: Reagent Name Normalization Fragile

**Files:** `process/reactor.py`, `process/conditions.py`, `process/safety.py`,
`process/costing.py`

**Issue:** Each module independently normalizes reagent names with
`r.lower().replace(" ", "_").replace("-", "_")`. Names like "LiAlH4" don't match
"lithium_aluminium_hydride". No canonical lookup exists.

**Fix:** Create a centralized `normalize_reagent_name()` function in
`reactions/reagent_data.py` that:
1. Lowercases
2. Strips whitespace/hyphens/underscores
3. Checks an alias dictionary (e.g., "lialh4" -> "lithium_aluminium_hydride")
4. Returns the canonical key

All modules should import and use this single function.

#### P1-CODE-3: File I/O Lacks Bounds Checking

**Files:** `io/mol_sdf.py`, `io/xyz.py`, `io/pdb.py`, `io/smiles_io.py`

**Issues:**
- `mol_sdf.py`: No check that file has enough lines for atom/bond counts
- `xyz.py`: No validation that atom_count matches actual atom lines
- `pdb.py`: No check that line length >= 78 for element symbol extraction
- `smiles_io.py`: `parse(smi)` unguarded; one bad SMILES crashes entire file read

**Fix for each:**
```python
# xyz.py - validate atom count
if len(lines) < atom_count + 2:
    raise ValueError(f"XYZ file declares {atom_count} atoms but has only "
                     f"{len(lines) - 2} atom lines")

# smiles_io.py - guard parse
try:
    mol = parse(smi)
except (ValueError, IndexError) as e:
    warnings.warn(f"Skipping invalid SMILES '{smi}': {e}")
    continue

# pdb.py - guard line length
element = line[76:78].strip() if len(line) >= 78 else ""

# mol_sdf.py - validate counts
if len(lines) < 4 + n_atoms + n_bonds:
    raise ValueError(f"MOL file truncated: expected {4+n_atoms+n_bonds} lines")
```

#### P1-CODE-4: Tolerance Inconsistencies in geometry.py

**File:** `molbuilder/core/geometry.py`

**Issue:** Three different near-zero thresholds used: `1e-12` (line 20),
`1e-10` (lines 67, 137). No rationale documented.

**Fix:** Define named constants at module level:
```python
_ZERO_VECTOR_TOL = 1e-12     # For zero-length vector detection
_COLLINEAR_TOL = 1e-10       # For collinear atom detection
```
Replace all magic numbers with these constants and add comments explaining
why different values are appropriate (vector magnitude vs cross-product magnitude).

#### P1-CODE-5: Empty __init__.py Files

**Files:** `core/__init__.py`, `bonding/__init__.py`, `atomic/__init__.py`,
`molecule/__init__.py`

**Issue:** No `__all__` definitions, no re-exports. Users must know exact
submodule paths to import anything.

**Fix:** Add `__all__` and convenience re-exports to each `__init__.py`:
```python
# core/__init__.py
"""Core data and utilities for molecular modeling."""
from molbuilder.core.constants import *
from molbuilder.core.elements import ELEMENTS, SYMBOL_TO_Z, from_symbol, from_name
from molbuilder.core.element_properties import electronegativity, covalent_radius_pm
from molbuilder.core.geometry import normalize, rotation_matrix, place_atom_zmatrix
from molbuilder.core.bond_data import bond_length, SP3_ANGLE, SP2_ANGLE, SP_ANGLE
```

Similarly for bonding/, atomic/, molecule/.

#### P1-CODE-6: Solvent Ratio Inconsistency

**Files:** `process/costing.py` vs `process/scale_up.py`

**Issue:** Costing uses `_SOLVENT_L_PER_KG = 8.0`, scale-up uses
`_SOLVENT_RATIO = 6.0` for the same concept.

**Fix:** Define once in `process/__init__.py`:
```python
DEFAULT_SOLVENT_L_PER_KG = 7.0  # Liters of solvent per kg product (compromise)
```
Import from both modules.

#### P1-CODE-7: Duplicate Reaction Templates

**File:** `molbuilder/reactions/knowledge_base.py`

**Issue:** These reactions appear twice:
- Williamson ether synthesis
- Beckmann rearrangement
- Claisen rearrangement

**Fix:** Remove duplicates. Keep the version with more complete data. Verify
test_reactions.py still passes after removal (template count will decrease).

---

## Medium-Priority Improvements (P2)

### P2-TYPE: Type Hints and Protocols

#### P2-TYPE-1: Add Type Hints to ~30% of Functions Missing Them

Priority files (most impact):
- `core/geometry.py`: `add_sp3_hydrogens(mol)` needs `mol: Molecule` type hint
- `reports/*.py`: All `generate_*()` functions use untyped duck typing.
  Define a `SynthesisStepProtocol` in `reports/__init__.py` using `typing.Protocol`
- `io/*.py`: All functions need return type annotations
- `atomic/quantum_atom.py`: `from_symbol()`, `from_name()` need `-> QuantumAtom`

#### P2-TYPE-2: Use Protocol for Duck-Typed Interfaces

The process/ and reports/ modules accept objects with `.template` attribute.
Define:
```python
from typing import Protocol, runtime_checkable

@runtime_checkable
class SynthesisStepLike(Protocol):
    template: ReactionTemplate
    step_number: int
    precursors: list
```

### P2-DATA: Data Completeness

#### P2-DATA-1: Extend Covalent Radii Beyond Element 69

**File:** `core/element_properties.py` COVALENT_RADII_PM

Add radii for elements 70-118. Source: Cordero et al. "Covalent radii revisited"
Dalton Trans. 2008, 2832-2838.

#### P2-DATA-2: Expand BDE Table

**File:** `core/bond_data.py` BDE_TABLE

Add at minimum: C-P, C-Si, N-S, N-Cl, O-F, O-S, S-F, Si-O, P-O, P-N.
Source: CRC Handbook of Chemistry and Physics, Luo "Comprehensive Handbook of
Chemical Bond Energies" (2007).

#### P2-DATA-3: Expand Torsion Barriers

**File:** `core/bond_data.py` TORSION_BARRIERS

Add sp2-sp2, sp2-sp3, aromatic-sp3 torsion parameters.
Source: OPLS-AA force field (Jorgensen et al. JACS 1996, 118, 11225).

#### P2-DATA-4: Add Missing Functional Group Detectors

**File:** `reactions/functional_group_detect.py`

Add detectors for: acid chlorides, anhydrides, sulfones, sulfoxides, imines,
acetals/ketals, silyl ethers, carbamates, hydrazones, oximes, enamines.

Each detector follows the existing pattern:
```python
def _detect_acid_chlorides(mol, results):
    for idx, atom in enumerate(mol.atoms):
        if _element(mol, idx) == "C":
            # Check for C(=O)Cl pattern
            ...
```

#### P2-DATA-5: Add Missing Reaction Types

**File:** `reactions/knowledge_base.py`

Add templates for: cross-metathesis, C-H activation (Pd-catalyzed),
asymmetric hydrogenation, carbene insertion, umpolung (dithiane chemistry),
olefin cross-metathesis (Grubbs), ring-opening metathesis polymerization.

#### P2-DATA-6: Expand GHS Hazard Statements

**File:** `process/safety.py` GHS_HAZARD_STATEMENTS

Add missing codes: H304, H305, H360F, H360D, H360Fd, H360Df, H361, H361f,
H361d, H361fd, H362, H370, H371, H372, H373.
Source: GHS Revision 10 (2023), UN ECE.

### P2-ALGO: Algorithm Improvements

#### P2-ALGO-1: Pre-build Name->Z Reverse Mapping

**File:** `core/elements.py` `from_name()`

Current: O(n) linear scan of ELEMENTS dict.
Fix: Build `_NAME_TO_Z: dict[str, int]` at module load time (like SYMBOL_TO_Z).

#### P2-ALGO-2: Cache most_probable_radius Results

**File:** `atomic/wavefunctions.py` `most_probable_radius()`

Current: Runs scipy.optimize.minimize_scalar on every call.
Fix: Add `@functools.lru_cache(maxsize=256)` decorator.

#### P2-ALGO-3: Retrosynthesis Visited Set Efficiency

**File:** `reactions/retrosynthesis.py` lines 1167-1170

Current: Creates new frozenset at each recursion level (immutable approach).
Fix: Pass mutable `set` by reference. Use `add()` before recursion,
`discard()` after (backtracking).

#### P2-ALGO-4: Add Zero-Distance Guard to coulombs_law

**File:** `core/constants.py` `coulombs_law()`

Add: `if r <= 0: raise ValueError("Distance must be positive")`

### P2-ROBUST: Robustness

#### P2-ROBUST-1: Add Warnings Module Usage

Many functions silently fall back to defaults when data is missing. Add
`import warnings` and emit `UserWarning` when:
- Covalent radius falls back to 150 pm
- BDE lookup returns None
- Reagent not found in REAGENT_DB
- Functional group detection may be unreliable (no explicit H)

#### P2-ROBUST-2: Parameterize Emergency Contact Numbers

**File:** `reports/safety_report.py` lines 263-267

Replace hardcoded US numbers (911, 1-800-222-1222) with configurable
constants at module level.

#### P2-ROBUST-3: Add Steric Clash Detection

**File:** `molecule/graph.py`

Add method `check_steric_clashes(min_distance: float = 0.5) -> list[tuple[int,int]]`
that returns pairs of non-bonded atoms closer than `min_distance` Angstroms.
Useful after building molecules to catch placement errors.

### P2-TEST: Test Coverage Gaps

#### P2-TEST-1: Add Negative/Edge Case Tests

Current tests mostly verify happy paths. Add:
- Invalid SMILES parsing (malformed strings, unclosed brackets)
- Empty molecule operations
- Zero-scale process engineering
- Molecules with no functional groups
- Single-atom molecules
- Very large molecules (100+ atoms) for performance
- All 7 critical bug verifications as regression tests

#### P2-TEST-2: Add Scientific Validation Tests

Compare computed values against known experimental data:
- Bohr model: H ionization energy = 13.6 eV
- Bond lengths: C-C single = 1.54 A, C=C double = 1.34 A, C-C triple = 1.20 A
- VSEPR angles: CH4 = 109.5, H2O = 104.5, NH3 = 107
- Dipole moments: HF = 1.82 D, HCl = 1.08 D, H2O = 1.85 D (after P0-1 fix)

---

## Implementation Roadmap

### Sprint 1: Critical Bug Fixes (P0-1 through P0-7)

**Order matters. Do these in sequence:**

1. **P0-1** (dipole moment) -- single property fix in covalent.py
2. **P0-2** (partial charges) -- verify or fix in covalent.py
3. **P0-7** (Lewis bond promotion) -- condition fix in lewis.py
4. **P0-5** (FG detection H atoms) -- fix detectors or add validation
5. **P0-6** (retro string transforms) -- add guards, plan graph-based refactor
6. **P0-3 + P0-4** (SMILES stereochemistry + bracket atoms) -- largest change,
   touches tokenizer/parser/writer/Atom dataclass

**After Sprint 1:** Run full test suite. Add regression tests for each fix.

### Sprint 2: High-Priority Fixes (P1)

**Parallelizable -- these are independent:**

- P1-SCI-1 (Slater's rules) + P1-SCI-2 (Bohr warning)
- P1-SCI-3 (reaction time) + P1-SCI-4 (safety incompatibilities)
- P1-SCI-5 (Aufbau exceptions) + P1-SCI-6 (5-membered aromatics)
- P1-CODE-1 (input validation) + P1-CODE-2 (reagent normalization)
- P1-CODE-3 (I/O bounds checking) + P1-CODE-4 (tolerance constants)
- P1-CODE-5 (__init__.py exports) + P1-CODE-6 (solvent ratio)
- P1-CODE-7 (duplicate templates) + P1-SCI-7 (boat conformation)

### Sprint 3: Medium Improvements (P2)

**Grouped by theme:**

- **Types:** P2-TYPE-1, P2-TYPE-2
- **Data:** P2-DATA-1 through P2-DATA-6
- **Algorithms:** P2-ALGO-1 through P2-ALGO-4
- **Robustness:** P2-ROBUST-1 through P2-ROBUST-3
- **Tests:** P2-TEST-1, P2-TEST-2

---

## Verification Commands

```bash
# Run full test suite (expect 435 passing before fixes, more after)
python -m pytest tests/ -q

# Run specific test file
python -m pytest tests/test_bonding.py -v

# Verify SMILES round-trip
python -c "
from molbuilder.smiles import parse, to_smiles
for smi in ['CCO', 'c1ccccc1', 'CC(=O)O', 'CC(=O)Oc1ccccc1C(=O)O']:
    mol = parse(smi)
    out = to_smiles(mol)
    print(f'{smi:30s} -> {out}')
"

# Verify functional group detection
python -c "
from molbuilder.smiles import parse
from molbuilder.reactions.functional_group_detect import detect_functional_groups
test_cases = [('CCO', 'alcohol'), ('CC=O', 'aldehyde'), ('CC(=O)O', 'carboxylic_acid')]
for smi, expected in test_cases:
    mol = parse(smi)
    fgs = [fg.name for fg in detect_functional_groups(mol)]
    ok = expected in fgs
    print(f'{smi:15s} expect={expected:20s} found={fgs!r:40s} {\"PASS\" if ok else \"FAIL\"}'  )
"

# Verify dipole moment (after P0-1 fix)
python -c "
from molbuilder.bonding.covalent import single_bond
for a, b, exp_range in [('H','F',(1.0,2.5)), ('H','Cl',(0.5,1.5)), ('C','C',(0,0.01))]:
    bond = single_bond(a, b)
    mu = bond.dipole_moment_debye
    ok = exp_range[0] <= mu <= exp_range[1]
    print(f'{a}-{b}: mu={mu:.3f} D  expected {exp_range}  {\"PASS\" if ok else \"FAIL\"}')
"

# Verify retrosynthesis
python -c "
from molbuilder.smiles import parse, to_smiles
from molbuilder.reactions.retrosynthesis import retrosynthesis, format_tree
mol = parse('CC(=O)Oc1ccccc1C(=O)O')  # aspirin
mol.name = 'aspirin'
tree = retrosynthesis(mol, max_depth=5, beam_width=5)
print(format_tree(tree)[:500])
"

# Verify process engineering pipeline
python -c "
from molbuilder.process.reactor import select_reactor
from molbuilder.process.costing import estimate_cost
from molbuilder.reactions.knowledge_base import REACTION_TEMPLATES
t = REACTION_TEMPLATES[0]
r = select_reactor(t, 10.0)
print(f'Reactor: {r.reactor_type.name}, {r.volume_L:.1f} L')
"

# Run full synthesis pipeline
python synthesize.py
# Enter: CCO (ethanol)
```

---

## Import Patterns

```python
# Core utilities
from molbuilder.core.constants import PLANCK_CONSTANT, SPEED_OF_LIGHT
from molbuilder.core.elements import ELEMENTS, SYMBOL_TO_Z, from_symbol
from molbuilder.core.element_properties import electronegativity, covalent_radius_pm
from molbuilder.core.geometry import normalize, rotation_matrix, place_atom_zmatrix
from molbuilder.core.bond_data import bond_length, SP3_ANGLE, BDE_TABLE

# Molecule building
from molbuilder.molecule.graph import Molecule, Atom, Bond, Hybridization
from molbuilder.molecule.builders import build_ethane, build_butane, build_cyclohexane
from molbuilder.molecule.functional_groups import add_hydroxyl, add_amino

# SMILES
from molbuilder.smiles import parse, to_smiles

# File I/O
from molbuilder.io import write_xyz, read_xyz, write_json, read_json
from molbuilder.io import write_mol, read_mol, write_pdb, read_pdb

# Reactions
from molbuilder.reactions import REACTION_TEMPLATES, detect_functional_groups
from molbuilder.reactions.retrosynthesis import retrosynthesis, format_tree
from molbuilder.reactions.synthesis_route import extract_best_route, format_route

# Process engineering
from molbuilder.process.reactor import select_reactor
from molbuilder.process.costing import estimate_cost
from molbuilder.process.safety import assess_safety

# Reports
from molbuilder.reports import generate_molecule_report
from molbuilder.reports import generate_synthesis_report
```

---

## API Gotchas

When writing new modules that operate on Molecule objects:

- Use `atom.symbol` (NOT `.element`) to get element string ("C", "H", "O")
- Use `mol.neighbors(idx)` to get bonded atom indices (returns `list[int]`)
- Use `mol.get_bond(i, j)` to get a Bond object (returns `Bond | None`)
- Bond order: `bond.order` (int: 1, 2, or 3)
- Bond atoms: `bond.atom_i`, `bond.atom_j` (int indices)
- Atom position: `atom.position` (numpy array of shape (3,))
- Atom index: `atom.index` (int, matches position in `mol.atoms` list)
- `valence_electrons` on QuantumAtom is a `@property`, NOT a method call
- `steric_number()` and `lone_pairs_on_central()` on LewisStructure ARE methods
- VSEPR molecular geometry: `vsepr.axe.molecular_geometry` (not `vsepr.molecular_geometry`)
- `build_2_butene(is_cis=True)` -- parameter is `is_cis`, not `cis`
- AXE notation includes explicit lone pair count: "AX3E1" not "AX3E"
- `from_symbol()` in `core/elements.py` returns `tuple[int, str, str, float]`

---

## Encoding and Style Rules

All `.py` files MUST be Windows cp1252 compatible:
- No emojis or unicode mathematical symbols
- Use ASCII: `->` not right-arrow, `>=` not greater-or-equal sign
- Degree symbol: use `deg` or "degrees" in strings
- Greek letters: spell out (alpha, beta) or use ASCII approximations

---

## Dependency Graph

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
  +-- gui/            (imports everything)
  +-- reports/        (imports reactions/, process/, molecule/, core/)
  +-- cli/            (imports everything)
```

---

## Data Sizes

| Database | Count | Location |
|----------|-------|----------|
| Elements | 118 | `core/elements.py` |
| Covalent radii | 118 | `core/element_properties.py` |
| Reagents | ~100 | `reactions/reagent_data.py` |
| Solvents | ~30 | `reactions/reagent_data.py` |
| Reaction templates | 91 | `reactions/knowledge_base.py` |
| Purchasable materials | ~200 | `reactions/retrosynthesis.py` |
| Bond lengths | 27 entries | `core/bond_data.py` |
| BDE values | 49 entries | `core/bond_data.py` |
| GHS hazard codes | 69 | `process/safety.py` |
| FG detectors | 21 | `reactions/functional_group_detect.py` |
| Amino acids | 20 | `molecule/amino_acids.py` |
| Torsion barriers | 16 types | `core/bond_data.py` |
| Tests | 517 | `tests/` (9 files) |

---

## Plan File

The original approved restructuring plan is at:
`C:\Users\taylor\.claude\plans\goofy-coalescing-lemur.md`

---

## Review Sources

The issues documented above were identified via comprehensive codebase review
using the best-coder and research-scientist analysis skills. Scientific claims
were verified against:

- Slater's Rules: https://en.wikipedia.org/wiki/Slater's_rules
- OpenSMILES Specification: http://opensmiles.org/opensmiles.html
- SMILES Stereochemistry: https://depth-first.com/articles/2020/05/04/stereochemistry-and-atom-parity-in-smiles/
- Dipole Moments: https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Physical_Properties_of_Matter/Atomic_and_Molecular_Properties/Dipole_Moments
- Pauling Ionic Character: https://www.omnicalculator.com/chemistry/percent-ionic-character
- CASP State of Art: https://arxiv.org/abs/2508.01459 (Andronov et al. 2025)
