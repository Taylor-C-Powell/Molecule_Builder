# Submission Drafts

Ready-to-post text for distribution channels. Copy-paste and submit.

---

## 1. r/cheminformatics Post

**Title:** I built an open-source Python toolkit that goes from SMILES to production conditions -- no RDKit needed

**Body:**

I've been building MolBuilder, a pure-Python molecular engineering toolkit that covers the full pipeline from molecular structure through retrosynthesis, reactor selection, safety assessment, cost estimation, and scale-up analysis.

The newest feature: give it a SMILES string and it predicts optimal reaction conditions:

```python
from molbuilder.process.condition_prediction import predict_conditions

result = predict_conditions("CCO", reaction_name="oxidation", scale_kg=10.0)
print(result.best_match.template_name)     # TEMPO-mediated oxidation
print(result.best_match.conditions.solvent) # DCM/water (biphasic)
print(result.overall_confidence)            # high
```

It analyzes the substrate's steric environment and electronic character, searches 91 reaction templates, scores candidates, and computes optimized conditions for your target scale.

**What makes it different from RDKit:**
- `pip install molbuilder` -- pure Python, no RDKit required
- Goes beyond cheminformatics into process engineering (reactor sizing, GHS safety, cost estimation, scale-up)
- 1,280+ tests, pure Python + numpy/scipy/matplotlib
- 91 reaction templates with retrosynthetic planning
- REST API available for integration

I'd appreciate any feedback from practicing chemists -- especially on whether the condition predictions align with your experience. The tutorial notebooks are in the repo if you want to try it.

- GitHub: https://github.com/Taylor-C-Powell/Molecule_Builder
- PyPI: `pip install molbuilder`
- Tutorials: https://github.com/Taylor-C-Powell/Molecule_Builder/tree/main/tutorials

---

## 2. Show HN Post

**Title:** Show HN: MolBuilder -- SMILES to manufacturing in pure Python (no RDKit)

**URL:** https://github.com/Taylor-C-Powell/Molecule_Builder

**Text:**

MolBuilder is a pure-Python molecular engineering toolkit. Give it a molecule as SMILES, and it can plan a synthesis (91 reaction templates, beam search), predict reaction conditions, select a reactor, assess safety (GHS), estimate costs, and analyze scale-up -- all without RDKit or any C++ dependencies.

Latest feature: substrate-aware condition prediction. It analyzes steric/electronic properties of your molecule and ranks compatible reaction templates with optimized conditions.

The scope is the differentiator: RDKit handles cheminformatics (descriptors, fingerprints, substructure search). MolBuilder covers the downstream pipeline -- retrosynthesis, reactor selection, safety, costing, scale-up. It's not trying to replace RDKit's computational accuracy, it's solving a different problem: going from "I have a molecule" to "here's how to make it at scale."

1,280+ tests, Apache 2.0 license, Python 3.11+. Also available as a REST API.

---

## 3. awesome-cheminformatics PR

**File to edit:** README.md in https://github.com/hsiaoyi0504/awesome-cheminformatics

**Add under "Libraries / Molecular Modeling":**

```markdown
- [MolBuilder](https://github.com/Taylor-C-Powell/Molecule_Builder) - Pure-Python molecular engineering toolkit covering SMILES parsing, retrosynthesis (91 templates), reaction condition prediction, reactor selection, safety assessment, cost estimation, and scale-up analysis. No RDKit dependency.
```

**PR Title:** Add MolBuilder -- pure-Python molecular engineering toolkit

**PR Body:**

MolBuilder is an open-source (Apache 2.0) molecular engineering toolkit in pure Python. It's unique in covering the full pipeline from molecular structure through process engineering:

- SMILES parser/writer with chirality and stereochemistry
- 91 reaction templates with beam-search retrosynthesis
- Substrate-aware reaction condition prediction
- Reactor selection, GHS safety, cost estimation, scale-up
- 1,280+ tests, Python 3.11+, pip-installable

PyPI: https://pypi.org/project/molbuilder/
GitHub: https://github.com/Taylor-C-Powell/Molecule_Builder

---

## 4. awesome-python-chemistry PR

**File to edit:** README.md in https://github.com/lmmentel/awesome-python-chemistry

**Add under "Molecular Modeling":**

```markdown
- [MolBuilder](https://github.com/Taylor-C-Powell/Molecule_Builder) - Pure-Python toolkit for molecular engineering: SMILES parsing, retrosynthesis planning, reaction condition prediction, process engineering (reactor, safety, costing, scale-up). No C++ dependencies.
```

**PR Title:** Add MolBuilder -- pure-Python molecular engineering and process chemistry toolkit

**PR Body:**

Adding MolBuilder, an Apache 2.0-licensed Python package for molecular engineering and process chemistry. Key features:

- Pure Python (numpy/scipy/matplotlib only) -- no RDKit or C++ required
- SMILES parser/writer, 21 functional group detectors
- Retrosynthetic planner with 91 reaction templates
- Process engineering: reactor selection, condition optimization, GHS safety, costing, scale-up
- 1,280+ tests, CI on Python 3.11/3.12/3.13

PyPI: https://pypi.org/project/molbuilder/

---

## 5. LinkedIn Post

I've been building MolBuilder -- an open-source Python toolkit for process chemistry.

The problem: most cheminformatics tools stop at molecular properties. If you're a process chemist, you still have to manually figure out reactor type, reaction conditions, safety assessment, and cost estimation for each synthesis step.

MolBuilder automates that full pipeline. Give it a SMILES string:

1. It plans a synthesis route (91 reaction templates, beam search)
2. Predicts optimal conditions (analyzing substrate sterics and electronics)
3. Selects a reactor (batch/CSTR/PFR/microreactor)
4. Runs GHS safety assessment (69 hazard codes, PPE, emergency procedures)
5. Estimates costs (materials, labor, equipment, energy, waste)
6. Analyzes scale-up (batch sizing, capital costs, annual capacity)

Pure Python, pip-installable, no RDKit needed. Apache 2.0 licensed with 1,280+ tests.

Especially interested in feedback from process development teams at CDMOs and small biotechs -- does this match your workflow needs?

GitHub: https://github.com/Taylor-C-Powell/Molecule_Builder
PyPI: pip install molbuilder

#cheminformatics #processchemistry #drugdiscovery #python #opensource

---

## 6. r/Python Post

**Title:** MolBuilder: pure-Python molecular engineering -- from SMILES to manufacturing plans

**Body:**

I built a Python package that handles the full chemistry pipeline from molecular structure to production planning, without any C++ dependencies.

```bash
pip install molbuilder
```

```python
from molbuilder.process.condition_prediction import predict_conditions

# What's the best way to oxidize ethanol at 10 kg scale?
result = predict_conditions("CCO", reaction_name="oxidation", scale_kg=10.0)
print(result.best_match.template_name)  # TEMPO-mediated oxidation
print(result.best_match.conditions.temperature_C)  # 5.0 C
print(result.overall_confidence)  # high
```

**Why I built it:** RDKit is great for cheminformatics but doesn't cover process engineering. I needed something that could go from "here's a molecule" to "here's how to manufacture it" -- retrosynthesis, reactor selection, safety, costing, and scale-up in one package.

**What it covers:**
- SMILES parsing with chirality/stereochemistry
- 91 reaction templates with retrosynthetic planning
- Substrate-aware condition prediction (new!)
- Reactor selection, GHS safety, cost estimation, scale-up analysis
- Lipinski Ro5, logP, TPSA, pKa prediction
- File I/O: XYZ, MOL/SDF, PDB, JSON

1,280+ tests, Python 3.11+, numpy/scipy/matplotlib only. Apache 2.0 licensed.

Jupyter tutorials are in the repo if you want to try it: https://github.com/Taylor-C-Powell/Molecule_Builder/tree/main/tutorials
