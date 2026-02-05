# Changelog

All notable changes to MolBuilder will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
