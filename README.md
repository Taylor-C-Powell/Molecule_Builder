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

# Detect functional groups (21 detectors, no SMARTS needed)
from molbuilder.reactions import detect_functional_groups
fgs = detect_functional_groups(mol)
# [carboxylic_acid, aromatic_ring, alcohol, ...]

# Plan a synthesis (91 reaction templates, beam search)
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

safety = assess_safety(route.steps)         # GHS hazards, PPE, emergency procedures
cost = estimate_cost(route.steps, 100.0)    # Materials, labor, equipment, energy, waste
scaleup = analyze_scale_up(route.steps, 10000.0)  # Batch sizing, capex, annual capacity
```

## Why not RDKit?

| | MolBuilder | RDKit |
|---|---|---|
| **Install** | `pip install molbuilder` | Requires conda or C++ compilation |
| **Dependencies** | numpy, scipy, matplotlib | C++ toolchain, Boost, many system libs |
| **Scope** | Atoms to manufacturing | Cheminformatics only |
| **Process engineering** | Reactor, costing, safety, scale-up | Not available |
| **Retrosynthesis** | 91 templates, beam search | Not included |
| **Learning curve** | Pure Python, readable source | C++ bindings, opaque internals |

MolBuilder is not a replacement for RDKit's computational accuracy. It's a different tool for a different job: **when you need to go from a molecule to a production plan**, not just compute descriptors. If you need RDKit's 3D coordinate generation, MolBuilder can use it as an optional backend: `pip install molbuilder[rdkit]`.

## Tutorials

Interactive Jupyter notebooks in [`tutorials/`](tutorials/):

1. **[From SMILES to Production Conditions](tutorials/01_smiles_to_conditions.ipynb)** -- Predict optimal reaction conditions for any substrate
2. **[Retrosynthesis Route Planning](tutorials/02_retrosynthesis_planning.ipynb)** -- Plan multi-step syntheses from purchasable materials
3. **[Process Engineering Pipeline](tutorials/03_process_engineering.ipynb)** -- Reactor selection, safety, costing, and scale-up

## Full capabilities

| Layer | What it does |
|-------|-------------|
| **Atomic physics** | Bohr model, quantum numbers, electron configurations, Slater's rules, hydrogen-like wavefunctions |
| **Bonding** | Lewis structures, VSEPR geometry (12+ shapes), covalent bond analysis, dipole moments |
| **Molecular modeling** | 3D coordinate generation, conformational analysis, Newman projections, R/S and E/Z stereochemistry |
| **Cheminformatics** | SMILES parser/writer with chirality and stereochemistry, Lipinski Ro5, logP, TPSA, pKa prediction |
| **Retrosynthesis** | 91 reaction templates, beam-search planner, 200+ purchasable starting materials, scored disconnections |
| **Condition prediction** | Substrate-aware template matching, steric/electronic analysis, solvent scoring, yield adjustment |
| **Process engineering** | Reactor selection, condition optimization, purification, GHS safety (69 hazard codes), cost estimation, scale-up |
| **File I/O** | XYZ, MOL/SDF V2000, PDB, JSON, SMILES |
| **Visualization** | 3D ball-and-stick rendering, quantum orbital plots, tkinter GUI editor |

## SaaS API

MolBuilder is also available as a hosted REST API:

```bash
# Get an API key at molbuilder-api.up.railway.app
curl -X POST https://molbuilder-api.up.railway.app/api/v1/process/predict-conditions \
  -H "X-API-Key: your-key" \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO", "reaction_name": "oxidation", "scale_kg": 10.0}'
```

Tiers: Free (100 req/day) | Pro ($49/mo) | Team ($199/mo) | Enterprise (custom)

## Quick reference

```python
# Parse SMILES
from molbuilder.smiles import parse, to_smiles
mol = parse("CC(=O)Oc1ccccc1C(=O)O")  # aspirin

# Functional groups
from molbuilder.reactions import detect_functional_groups
fgs = detect_functional_groups(mol)

# Lipinski properties
from molbuilder.molecule.properties import lipinski_properties
props = lipinski_properties(mol)  # MW, logP, HBD, HBA, TPSA, Ro5

# Condition prediction
from molbuilder.process.condition_prediction import predict_conditions
result = predict_conditions("CCO", reaction_name="oxidation")

# Retrosynthesis
from molbuilder.reactions.retrosynthesis import retrosynthesis, format_tree
tree = retrosynthesis(mol, max_depth=5, beam_width=5)
print(format_tree(tree))

# File I/O
from molbuilder.io import write_xyz, write_mol, write_pdb
write_xyz(mol, "aspirin.xyz")
write_mol(mol, "aspirin.mol")

# Interactive CLI
# python -m molbuilder
```

## Testing

1,280+ tests across 23 test files. CI enforces 80% coverage.

```bash
pip install -e ".[dev]"
python -m pytest tests/ -q
```

## License

MIT License. See [LICENSE](LICENSE) for details.

Copyright (c) 2025-2026 Taylor C. Powell.
