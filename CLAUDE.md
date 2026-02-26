# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MolBuilder is a professional-grade molecular engineering toolkit with five components:

- **`molbuilder/`** -- Core Python library (PyPI: `molbuilder` v1.2.0). Pure Python + numpy/scipy/matplotlib. Optional RDKit backend for 3D coords (`pip install molbuilder[rdkit]`).
- **`saas/`** -- FastAPI REST API (PyPI: `molbuilder-api`). JWT auth, tiered API keys, RBAC, Stripe billing. Deployed to Railway.
- **`frontend/`** -- React SPA web dashboard (Vite + TypeScript + Tailwind CSS 4 + 3Dmol.js). Deployed to Vercel.
- **`studio/`** -- React SPA 3D molecule editor (Vite + TypeScript + Three.js/R3F). Deployed to Vercel.
- **`sdk/`** -- Python SDK client (PyPI: `molbuilder-client`). Uses httpx.

## Build & Install

```bash
# Core library (Python >= 3.11 required)
pip install -e ".[dev]"

# SaaS API
pip install -e "./saas[dev]"

# SDK
pip install -e "./sdk[dev]"

# Frontend (Node 22)
cd frontend && npm ci

# Studio (Node 22)
cd studio && npm ci
```

## Running Tests

```bash
# Core library (1435 tests across 30 files, CI enforces --cov-fail-under=80)
python -m pytest tests/ -v
python -m pytest tests/test_smiles.py -v          # single file
python -m pytest tests/ --cov=molbuilder --cov-report=term-missing

# SaaS API (151 tests across 18 files, uses FastAPI TestClient with temp SQLite DBs)
cd saas && python -m pytest tests/ -v

# SDK (uses respx for HTTP mocking, pytest-asyncio with asyncio_mode=auto)
python -m pytest sdk/tests/ -v

# Frontend / Studio (no unit test runner -- CI checks lint + typecheck + build)
cd frontend && npm run lint && npm run typecheck && npm run build
cd studio && npm run lint && npm run typecheck && npm run build
```

## Linting & Formatting

```bash
# Python (Ruff) -- line length 100, target py311
ruff check molbuilder/ saas/app/ sdk/molbuilder_client/
ruff format molbuilder/

# TypeScript (ESLint 9)
cd frontend && npm run lint
cd studio && npm run lint
```

## Architecture

### Core Library Dependency Graph

```
core/  <-- everything depends on this
  +-- atomic/      (Bohr model, quantum numbers, wavefunctions)
  +-- bonding/     (Lewis structures, VSEPR, covalent bond analysis)
  +-- molecule/    (Molecule graph -- central class)
  |     +-- smiles/  (SMILES tokenizer, parser, writer)
  |     +-- io/      (XYZ, MOL/SDF V2000, PDB, JSON)
  +-- coords/      (3D coordinate generation -- DG+FF builtin or RDKit backend)
  +-- dynamics/    (MD engine: force field, Verlet integrator, reaction mechanisms)
  +-- reactions/   (185 templates, 24 FG detectors, retrosynthesis engine, RetroCast adapter)
  |     +-- process/  (reactor, conditions, costing, safety, scale-up)
  +-- reports/     (ASCII report generators)
  +-- visualization/, gui/, cli/
```

### Central Data Type: `Molecule` (molbuilder/molecule/graph.py)

```python
class Molecule:
    name: str
    atoms: list[Atom]          # Atom(symbol, position, index, hybridization)
    bonds: list[Bond]          # Bond(atom_i, atom_j, order, rotatable)
    _adj: dict[int, list[int]] # adjacency list
```

### SaaS API

FastAPI with Pydantic v2. Auth flow: API key (X-API-Key) -> JWT exchange. Five tiers: free/pro/team/academic/enterprise. RBAC roles: admin/chemist/viewer. SQLite databases for users, audit trail, molecule store.

### Frontend & Studio

React 19, react-router-dom v7, Zustand for state. API proxied via Vite dev server to Railway backend. AuthGuard protects authenticated routes. Lazy-loaded pages with Suspense.

## API Gotchas

- **Use `atom.symbol`, NOT `.element`** -- this is the field name for element strings ("C", "H", "O")
- `mol.neighbors(idx)` returns `list[int]`; `mol.get_bond(i, j)` returns `Bond | None`
- `valence_electrons` on QuantumAtom is a `@property`, NOT a method call
- `steric_number()` and `lone_pairs_on_central()` on LewisStructure ARE method calls
- VSEPR molecular geometry: `vsepr.axe.molecular_geometry` (not `vsepr.molecular_geometry`)
- `build_2_butene(is_cis=True)` -- parameter is `is_cis`, not `cis`
- AXE notation includes explicit lone pair count: "AX3E1" not "AX3E"
- `generate_3d(mol)` modifies positions **in-place** (returns None). Backend: `"auto"` tries RDKit then builtin DG+FF.
- Coverage excludes `gui/`, `cli/menu.py`, `cli/demos.py` (interactive/display code)

## Encoding Rules

All `.py` files MUST be Windows cp1252 compatible -- no emojis, no unicode math symbols. Use ASCII only: `->` not right-arrow, `deg` not degree symbol, spell out Greek letters.

## Key Import Patterns

```python
from molbuilder.molecule.graph import Molecule, Atom, Bond, Hybridization
from molbuilder.smiles import parse, to_smiles
from molbuilder.core.constants import PLANCK_CONSTANT, SPEED_OF_LIGHT
from molbuilder.core.elements import ELEMENTS, SYMBOL_TO_Z, from_symbol
from molbuilder.io import write_xyz, read_xyz, write_mol, read_mol
from molbuilder.coords import generate_3d              # 3D coord generation (auto/builtin/rdkit)
from molbuilder.dynamics import MDSimulation, ForceField  # molecular dynamics
from molbuilder.reactions import REACTION_TEMPLATES, detect_functional_groups
from molbuilder.reactions.retrosynthesis import retrosynthesis, format_tree
from molbuilder.process.reactor import select_reactor
from molbuilder.process.costing import estimate_cost
```

## Entry Points

- `python -m molbuilder` -- interactive CLI menu
- `python synthesize.py` -- SMILES-to-synthesis pipeline
- `uvicorn app.main:app` (from saas/) -- API server
- `npm run dev` (from frontend/) -- dashboard dev server, port 5173
- `npm run dev` (from studio/) -- 3D editor dev server, port 5174

## CI/CD

- **test.yml**: Python 3.11/3.12/3.13 matrix. Ruff lint, core tests (coverage >= 80%), saas tests, sdk tests.
- **frontend-test.yml / studio-test.yml**: Node 22. Lint, typecheck, build.
- **deploy.yml**: Auto-deploys SaaS to Railway on push to main (when saas/ or molbuilder/ change).
- **release-sdk.yml**: Publishes SDK to PyPI on `sdk-v*` tags.

## Reference

See `DEVELOPMENT_GUIDE.md` for detailed architecture, completed sprint history, and verification commands.
