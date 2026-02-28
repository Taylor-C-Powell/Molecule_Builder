# MolBuilder

**From SMILES to manufacturing in pure Python.**

![PyPI](https://img.shields.io/pypi/v/molbuilder) ![Python](https://img.shields.io/pypi/pyversions/molbuilder) ![Tests](https://github.com/Taylor-C-Powell/Molecule_Builder/actions/workflows/test.yml/badge.svg) ![Downloads](https://img.shields.io/pypi/dm/molbuilder) ![License](https://img.shields.io/pypi/l/molbuilder)

```python
from molbuilder.process.condition_prediction import predict_conditions

result = predict_conditions("CCO", reaction_name="oxidation", scale_kg=10.0)

print(result.best_match.template_name)        # TEMPO-mediated oxidation
print(result.best_match.conditions.solvent)    # DCM/water (biphasic)
print(result.best_match.conditions.temperature_C)  # 5.0
print(result.best_match.adjusted_yield_range)  # (85.0, 98.0)
print(result.overall_confidence)               # high
```

MolBuilder is the only open-source Python package that covers the full pipeline from molecular structure through retrosynthesis, reactor selection, safety assessment, cost estimation, and scale-up analysis. One `pip install`, zero C++ dependencies.

```bash
pip install molbuilder
```

## What can it do?

```python
# Parse any SMILES string
from molbuilder.smiles import parse
mol = parse("CC(C)Cc1ccc(cc1)C(C)C(=O)O")  # ibuprofen

# Detect functional groups (24 detectors + SMARTS cross-validation)
from molbuilder.reactions import detect_functional_groups
fgs = detect_functional_groups(mol)
# [carboxylic_acid, aromatic_ring, alcohol, ...]

# Plan a synthesis (185 reaction templates, beam search)
from molbuilder.reactions.retrosynthesis import retrosynthesis
from molbuilder.reactions.synthesis_route import extract_best_route
tree = retrosynthesis(mol, max_depth=5, beam_width=5)
route = extract_best_route(tree)

# Full process engineering for each step
from molbuilder.process.reactor import select_reactor
from molbuilder.process.conditions import optimize_conditions
from molbuilder.process.safety import assess_safety
from molbuilder.process.costing import estimate_cost
from molbuilder.process.scale_up import analyze_scale_up

for step in route.steps:
    reactor = select_reactor(step.template, scale_kg=100.0)   # Batch/CSTR/PFR/microreactor
    conditions = optimize_conditions(step.template, scale_kg=100.0)  # Temp, solvent, time, atmosphere

safety = assess_safety(route.steps)         # GHS hazards, thermal hazards, PPE, incompatibilities
cost = estimate_cost(route.steps, 100.0)    # Materials, labor, equipment, energy, waste
scaleup = analyze_scale_up(route.steps, 10000.0)  # Batch sizing, capex, annual capacity
```

## Why not RDKit?

RDKit is a mature, excellent cheminformatics toolkit. MolBuilder is a different tool for a different job.

| | MolBuilder | RDKit |
|---|---|---|
| **Scope** | Atoms to manufacturing | Cheminformatics (descriptors, fingerprints, substructure) |
| **Process engineering** | Reactor, costing, safety, scale-up | Not available |
| **Retrosynthesis** | 185 templates, beam search | Not included |
| **Source language** | Pure Python (readable, hackable) | C++ with Python bindings |
| **Dependencies** | numpy, scipy, matplotlib | Larger compiled package |

MolBuilder is not a replacement for RDKit's computational accuracy. It's a different tool for a different job: **when you need to go from a molecule to a production plan**, not just compute descriptors. If you need RDKit's 3D coordinate generation, MolBuilder can use it as an optional backend: `pip install molbuilder[rdkit]`.

## Full capabilities

| Layer | What it does |
|-------|-------------|
| **Atomic physics** | Bohr model, quantum numbers, electron configurations, Slater's rules, hydrogen-like wavefunctions |
| **Bonding** | Lewis structures, VSEPR geometry (12+ shapes), covalent bond analysis, dipole moments |
| **Molecular modeling** | 3D coordinate generation (builtin DG+FF or RDKit backend), conformational analysis, Newman projections, R/S and E/Z stereochemistry |
| **Cheminformatics** | SMILES parser/writer with chirality and stereochemistry, Lipinski Ro5, logP, TPSA, pKa prediction, synthetic accessibility scoring |
| **Pattern matching** | SMARTS engine for substructure search, 24 FG detectors with SMARTS cross-validation |
| **Retrosynthesis** | 185 reaction templates across 14 categories, beam-search planner, 270+ purchasable starting materials, scored disconnections |
| **Condition prediction** | Substrate-aware template matching, ORD-backed empirical conditions (180 reaction types), steric/electronic analysis, solvent scoring, yield adjustment |
| **Process engineering** | Reactor selection, condition optimization, purification, GHS safety (69 hazard codes), thermal hazard detection, reagent-solvent incompatibility checks, cost estimation (171 priced reagents), scale-up analysis |
| **Molecular dynamics** | Classical force field (LJ, Coulomb, harmonic, OPLS-AA torsions), Velocity Verlet integrator, reaction mechanism templates (SN2, E2, radical, carbonyl addition), slow-motion visualization |
| **File I/O** | XYZ, MOL/SDF V2000, PDB, JSON, SMILES |
| **Reports** | ASCII and PDF report generators for molecules, synthesis, safety, and cost |
| **Visualization** | 3D ball-and-stick rendering, quantum orbital plots, tkinter GUI editor, animation export (MP4/GIF) |

## SaaS API

MolBuilder is also available as a hosted REST API with 15 router modules:

```bash
# Get an API key at molbuilder-api.up.railway.app
curl -X POST https://molbuilder-api.up.railway.app/api/v1/process/predict-conditions \
  -H "X-API-Key: your-key" \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO", "reaction_name": "oxidation", "scale_kg": 10.0}'
```

**API features:**
- Molecule parsing, properties, 3D coordinates
- Retrosynthesis route planning
- Condition prediction with ORD-backed data
- Process engineering (reactor, safety, costing, scale-up)
- Batch processing (submit/poll async jobs, up to 100 molecules)
- Molecule library (save, search, tag, organize)
- File I/O (upload MOL/SDF/PDB/XYZ, export any format)
- PDF report generation
- Feasibility analysis
- Element lookup (periodic table data)
- Usage analytics and audit trail (21 CFR Part 11 compatible)
- API versioning with deprecation warnings

**Auth:** API key + JWT exchange, RBAC (admin/chemist/viewer)

**Tiers:** Free (100 req/day) | Pro ($49/mo) | Team ($199/mo) | Academic (free with .edu) | Enterprise (custom)

## Python SDK

```bash
pip install molbuilder-client
```

```python
from molbuilder_client import MolBuilder

client = MolBuilder(api_key="mb_...")

# Parse and analyze
mol = client.parse_smiles("CCO")
props = client.get_properties(mol.id)

# Batch processing
job = client.submit_batch(["CCO", "c1ccccc1", "CC(=O)O"], job_type="properties")
result = client.wait_for_batch(job.job_id)

# Molecule library
client.save_molecule(mol.id, tags=["alcohols"])
library = client.list_molecules()

# File I/O
client.import_file("molecule.xyz")
client.export_file("mol_id", "pdb", save_to="output.pdb")

# PDF reports
pdf = client.download_report("CCO", scale_kg=10.0, save_to="report.pdf")
```

Async client also available: `from molbuilder_client import AsyncMolBuilder`

## Web Dashboard

React 19 web interface with:
- Molecule parsing and 3D visualization (3Dmol.js)
- Retrosynthesis route explorer
- Process engineering pipeline
- Batch job submission and monitoring
- Molecule library management
- File upload/download
- SA score display with traffic-light coloring
- PDF report download

## Tutorials

Interactive Jupyter notebooks in [`tutorials/`](tutorials/):

1. **[From SMILES to Production Conditions](tutorials/01_smiles_to_conditions.ipynb)** -- Predict optimal reaction conditions for any substrate
2. **[Retrosynthesis Route Planning](tutorials/02_retrosynthesis_planning.ipynb)** -- Plan multi-step syntheses from purchasable materials
3. **[Process Engineering Pipeline](tutorials/03_process_engineering.ipynb)** -- Reactor selection, safety, costing, and scale-up

Written tutorials in [`docs/tutorials/`](docs/tutorials/):

1. From SMILES to 3D: Building Molecules
2. Drug-Likeness Screening with Lipinski Rule of Five
3. 3D Visualization and Conformational Analysis
4. Quantum Mechanics: Atomic Structure and Orbitals
5. Retrosynthetic Analysis: Planning a Synthesis Route
6. Process Engineering: From Lab to Plant

## Examples

See [`examples/`](examples/) for self-contained demo scripts:

- **`smiles_to_manufacturing.py`** -- Full pipeline: SMILES to FG detection, retrosynthesis, costing, safety, scale-up
- **`case_study_ibuprofen.py`** -- Complete ibuprofen synthesis case study

## Quick reference

```python
# Parse SMILES
from molbuilder.smiles import parse, to_smiles
mol = parse("CC(=O)Oc1ccccc1C(=O)O")  # aspirin

# Functional groups (24 detectors)
from molbuilder.reactions import detect_functional_groups
fgs = detect_functional_groups(mol)

# SMARTS substructure search
from molbuilder.smarts import SmartsMatcher
matcher = SmartsMatcher("[CX3](=O)[OX2H1]")  # carboxylic acid pattern
matches = matcher.match(mol)

# Lipinski properties
from molbuilder.molecule.properties import lipinski_properties
props = lipinski_properties(mol)  # MW, logP, HBD, HBA, TPSA, Ro5

# Synthetic accessibility
from molbuilder.molecule.sa_score import sa_score
result = sa_score(mol)  # 1 (easy) to 10 (hard)

# Condition prediction (ORD-backed)
from molbuilder.process.condition_prediction import predict_conditions
result = predict_conditions("CCO", reaction_name="oxidation")

# Retrosynthesis
from molbuilder.reactions.retrosynthesis import retrosynthesis, format_tree
tree = retrosynthesis(mol, max_depth=5, beam_width=5)
print(format_tree(tree))

# 3D coordinates
from molbuilder.coords import generate_3d
generate_3d(mol)  # modifies positions in-place

# Molecular dynamics
from molbuilder.dynamics import MDSimulation, ForceField
sim = MDSimulation(mol, ForceField())

# File I/O
from molbuilder.io import write_xyz, write_mol, write_pdb
write_xyz(mol, "aspirin.xyz")
write_mol(mol, "aspirin.mol")

# PDF reports (requires pip install molbuilder[pdf])
from molbuilder.reports.pdf_report import generate_molecule_pdf
pdf_bytes = generate_molecule_pdf(mol)

# Interactive CLI
# python -m molbuilder
```

## Data

| Database | Count |
|----------|-------|
| Elements (IUPAC 2021) | 118 |
| Covalent radii | 118 |
| Reaction templates | 185 (14 categories) |
| FG detectors | 24 |
| Reagents (with pricing tiers) | 171 |
| Solvents | 32 |
| Purchasable starting materials | 270+ |
| GHS hazard codes | 69 |
| Thermal hazard patterns | 9 |
| BDE values | 49 |
| Torsion barrier types | 16 |
| Bond length entries | 27 |
| Amino acids | 20 |
| ORD reaction types | 180 |
| SMARTS FG patterns | 24 |

## Testing

1,749 tests across 3 test suites. CI matrix: Python 3.11 / 3.12 / 3.13.

```bash
pip install -e ".[dev]"
python -m pytest tests/ -q          # 1,510 core library tests (36 files)
cd saas && python -m pytest tests/  # 175 SaaS API tests (22 files)
python -m pytest sdk/tests/         # 64 SDK tests (13 files)
```

Coverage gate: 80% minimum enforced in CI.

## Architecture

Five components:

```
molbuilder/         Core Python library (PyPI: molbuilder)
saas/               FastAPI REST API (Railway)
frontend/           React web dashboard (Vercel)
studio/             React 3D molecule editor (Vercel)
sdk/                Python SDK client (PyPI: molbuilder-client)
```

## Future directions

- **Learned retrosynthesis**: Transformer/GNN expansion policy trained on USPTO, with MCTS search
- **ADMET prediction**: hERG, CYP inhibition, Caco-2 permeability, metabolic stability
- **ELN/LIMS integration**: Connectors for Benchling, Signals Notebook
- **Solubility and crystallization**: Prediction models for process development
- **Team management**: Organizations, shared libraries, per-seat billing
- **SSO integration**: Auth0/Okta for enterprise deployments
- **On-premise deployment**: Docker + Kubernetes packaging
- **Patent landscape**: Google Patents API integration for FTO analysis

## License

Apache License 2.0. See [LICENSE](LICENSE) for details.

Copyright 2025-2026 Taylor C. Powell.
