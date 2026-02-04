# Changelog

All notable changes to MolBuilder will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
