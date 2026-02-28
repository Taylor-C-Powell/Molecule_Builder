# Changelog

All notable changes to MolBuilder will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0] - 2026-02-28

### Added

#### Core Library
- **SMARTS Pattern Matching Engine** (`smarts/`): Full substructure search with atom primitives (element, charge, H-count, ring membership, aromaticity), bond primitives (order, ring, aromatic), recursive SMARTS, and `SmartsMatcher` class
- **Open Reaction Database Integration** (`data/`): 180 reaction types with pre-computed statistics from ORD; `predict_conditions()` returns empirical temperature, solvent, and atmosphere data with `data_source` attribution ("ORD (n=X)" or "heuristic")
- **Synthetic Accessibility Scoring** (`molecule/sa_score.py`): Ertl-style SA score (1-10 scale) with components: ring_complexity, stereo_penalty, fragment_penalty, size_penalty, sp3_bonus
- **Thermal Hazard Detection** (`process/safety.py`): 9 exothermic reaction patterns (Grignard, diazotization, catalytic hydrogenation, Friedel-Crafts, nitration, sulfonation, ozonolysis, metalation, Wolff-Kishner) with severity ratings and max temperature limits
- **FG SMARTS Cross-Validation** (`reactions/fg_smarts_validation.py`): Validates heuristic FG detection against SMARTS patterns; reports agreement/disagreement with confidence scoring
- **RetroCast Adapter** (`reactions/retrocast_adapter.py`): `tree_to_retrocast_routes()` and `export_retrocast_json()` for benchmarking against 11+ synthesis planners
- **PDF Report Generation** (`reports/pdf_report.py`): ReportLab-based professional reports via `generate_molecule_pdf()` and `generate_process_pdf()` (optional dep: `pip install molbuilder[pdf]`)
- **Condition Prediction Module** (`process/condition_prediction.py`): Substrate-aware template matching with steric/electronic analysis, solvent scoring, and ORD-backed yield adjustment
- **Lipinski Properties** (`molecule/properties.py`): `lipinski_properties()` returning logP, HBD, HBA, TPSA, rotatable bonds, heavy atom count, Ro5 violations
- 3 new FG detectors: acid chlorides, anhydrides, imines (total: 24)
- 94 new reaction templates (total: 185 across 14 categories)
- 71 new reagents with `pricing_tier` field -- commodity/standard/specialty/exotic (total: 171)
- 70+ new purchasable starting materials (total: 270+)
- Solvent-reagent incompatibility checks in `assess_safety()` (NaH+water, n-BuLi+protic solvents, DMSO high-temp decomposition)
- Configurable costing model: `CostParameters` dataclass with region presets, customizable labor/waste/overhead rates
- Example scripts: `examples/smiles_to_manufacturing.py` (full pipeline), `examples/case_study_ibuprofen.py`

#### SaaS API
- **Batch Processing API** (`routers/batch.py`): Submit/poll/list/cancel async batch jobs, tier-limited (free: 10 molecules, pro: 100), threaded worker with 300s timeout
- **Molecule Library API** (`routers/library.py`): Save/get/list/update/delete molecules, tag-based organization, properties caching, import from SMILES
- **File I/O API** (`routers/file_io.py`): Upload XYZ/MOL/SDF/PDB/JSON files, export molecules in any format
- **PDF Reports API** (`routers/reports.py`): POST `/reports/process-pdf` returns PDF bytes
- **Feasibility API** (`routers/feasibility.py`): Quick feasibility analysis for synthesis routes
- **Analytics API** (`routers/analytics.py`): Usage statistics and trends
- **Audit API** (`routers/audit.py`): 21 CFR Part 11 compatible audit trail with HMAC integrity
- **API Versioning** (`middleware_versioning.py`): Version headers, deprecation warnings, sunset dates
- **Billing Integration**: Stripe checkout, customer portal, yearly toggle, tier upgrades
- **Legal Pages** (`routers/legal.py`): Terms of service and privacy policy endpoints
- RequestID middleware (X-Request-ID on every response)
- Security headers middleware (X-Content-Type-Options, X-Frame-Options, X-XSS-Protection, Referrer-Policy, HSTS)
- Usage tracking middleware with rate limiting
- Properties caching via `@lru_cache(maxsize=512)`
- Sentry error monitoring integration
- Admin bootstrap with secure key file generation

#### SDK (v0.2.0)
- Batch processing: `submit_batch()`, `get_batch_status()`, `list_batches()`, `cancel_batch()`, `wait_for_batch()`
- Molecule library: `save_molecule()`, `get_molecule()`, `list_molecules()`, `update_molecule()`, `delete_molecule()`
- File I/O: `import_file()`, `export_file()`
- Reports: `download_report()`
- Async client (`AsyncMolBuilder`) with all methods

#### Frontend
- Library page with molecule management
- Batch processing page with job monitoring
- Process engineering page
- Retrosynthesis page with route visualization
- File upload/download support
- SA score display with traffic-light coloring
- PDF report download buttons

### Fixed
- Broken redox transforms: aldehyde/ketone oxidation and carboxylic acid reduction now produce correct SMILES via graph-based `_replace_carbonyl_with_carboxyl()` and `_replace_carboxyl_with_carbonyl()`
- Retrosynthesis string transforms replaced with graph-based operations to prevent false positives
- Ring detection and FG detection results cached on Molecule instances
- Batch worker race condition: status checks prevent timeout overwrite by daemon threads

### Changed
- Reagent database expanded from ~100 to 171 entries with `pricing_tier` classification
- Safety assessment now includes thermal hazards and solvent-reagent incompatibilities
- `estimate_cost()` uses `CostParameters` dataclass instead of hardcoded values
- `optimize_conditions()` checks ORD empirical data before falling back to heuristics
- `ReactionConditions` and `TemplateMatch` include `data_source` field
- Test suite grew from 585 to 1,510 tests (36 files)
- SaaS test suite: 175 tests (22 files)
- SDK test suite: 64 tests (13 files)

## [1.1.1] - 2026-02-08

### Fixed
- Project URLs in pyproject.toml now point to the correct GitHub repository
- README installation section updated with correct `pip install molbuilder` and clone URL

### Changed
- Top-level `molbuilder` package now exports `Molecule`, `Atom`, `Bond`, `Hybridization`, `parse`, and `to_smiles` for convenient imports
- `molbuilder.core` replaces wildcard import with explicit constant imports
- `molbuilder.io` now exports all read/write functions (`read_xyz`, `write_xyz`, etc.)
- README test count updated from 517 to 585
- Copyright year updated to 2025-2026

### Added
- CONTRIBUTING.md with development setup, testing, and contribution guidelines

## [1.1.0] - 2026-02-05

### Added
- **Molecular Dynamics Engine** (`dynamics/`): Classical force field with Lennard-Jones, Coulomb, harmonic bond/angle, and OPLS-AA torsional terms; Velocity Verlet integrator with Berendsen thermostat; trajectory storage with CubicSpline interpolation for sub-femtosecond time resolution
- **Reaction Mechanism Templates**: Data model for multi-stage reaction mechanisms with electron flow arrows; predefined templates for SN2, E2, radical substitution, and nucleophilic addition to carbonyl
- **Steered-MD Choreography**: Sigmoid-ramped harmonic restraints that guide atoms through mechanism stages, producing smooth bond formation/breaking events with fractional bond orders
- **Extreme Slow-Motion Visualization**: FuncAnimation renderer with configurable slowdown factor (default 1 fs = 1 second), CPK atom rendering, fractional bond-order line styling (dashed for forming, fading for breaking), energy bar overlay, and SI-prefixed time labels
- **Electron Density Rendering**: LCAO-based Monte Carlo point cloud visualization of bonding electron density during bond events, using Slater effective nuclear charges
- **Playback Controls**: Keyboard bindings (Space, arrows, E, L, R, Q) and optional tkinter panel for interactive animation control (play/pause, step, speed, toggle overlays)
- **GUI Integration**: "Simulate" menu in MolBuilder GUI with MD Vibration, Bond Formation, SN2 Mechanism, and Export Animation commands
- **CLI Demo**: `demo_slowmo_interaction()` with interactive options for ethane MD vibration, SN2 mechanism, and C-C bond formation visualization
- **Animation Export**: MP4 (FFMpeg) and GIF (Pillow) export with configurable DPI and frame rate
- **68 new tests**: ForceField parameterization, force computation correctness, Velocity Verlet integration, energy conservation validation (<5% deviation over 1000 steps for harmonic diatomic), trajectory interpolation, mechanism templates, choreography restraints, animation pipeline (headless Agg backend), electron density rendering, GIF export verification

### Scientific Validation
- Energy conservation verified: harmonic diatomic conserves total energy to <5% over 1000 Verlet steps (0.5 fs timestep, NVE ensemble)
- C-C stretch vibration period validated in 15-100 fs range (expected ~33 fs for 1000 cm^-1 stretch)
- UFF Lennard-Jones parameters from Rappe et al., J. Am. Chem. Soc. 1992
- OPLS-AA torsion parameters from Jorgensen et al., J. Am. Chem. Soc. 1996
- Velocity Verlet integrator per Swope et al., J. Chem. Phys. 1982

## [1.0.0] - 2026-02-03

### Added
- **Atomic Physics**: Bohr model, quantum numbers, electron configurations with 22 Aufbau exceptions, Slater's rules for effective nuclear charge, hydrogen-like wavefunctions, orbital probability densities
- **Chemical Bonding**: Lewis structures with octet/expanded octet support, VSEPR geometry prediction (12+ molecular shapes), covalent bond analysis with dipole moments and BDE data
- **Molecular Modeling**: Full molecular graph with 3D coordinate generation, conformational analysis, torsional energy profiles, Newman projections, cyclohexane conformations, R/S and E/Z stereochemistry
- **Biochemistry**: All 20 standard amino acids with L-chirality, peptide bond formation, phi/psi angles, secondary structure templates
- **Cheminformatics**: SMILES parser and writer with chirality, E/Z bond stereochemistry, bracket atoms, aromatic perception, Morgan canonical ordering
- **Retrosynthesis**: Beam-search engine with 91 curated reaction templates across 14 categories, 200+ purchasable starting materials database, scored disconnections
- **Process Engineering**: Reactor selection (batch/CSTR/PFR/microreactor), condition optimization, purification strategy, GHS safety assessment (69 hazard codes), cost estimation, batch-vs-continuous scale-up
- **Interfaces**: Interactive CLI (11 menu options), SMILES-to-synthesis pipeline, Tkinter 3D GUI with molecular editor, Python API
- **File I/O**: XYZ, MOL/SDF (V2000), PDB, JSON, SMILES formats
- **Data**: 118 elements (IUPAC 2021), 118 covalent radii, 49 BDE entries, 16 torsion barrier types, ~100 reagents, ~30 solvents, 69 GHS codes, 20 amino acids
- **Testing**: 517 tests covering core chemistry, atomic models, bonding, molecular operations, SMILES round-trips, file I/O, reactions, process engineering, and edge cases with scientific validation
- **Reports**: ASCII report generators for molecules, synthesis routes, safety assessments, and cost estimates

### Scientific Validation
- Bohr model verified against hydrogen ionization energy (13.6 eV)
- Slater's rules validated (Fe Z_eff(4s) = 3.75, Z_eff(3d) = 6.25)
- Bond lengths validated against experimental averages (C-C 1.54 A, C=C 1.34 A)
- VSEPR angles validated (CH4 109.5, H2O 104.5, NH3 107)
- Dipole moment formula validated against HF experimental (1.82 D)
