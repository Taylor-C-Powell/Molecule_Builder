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

# SaaS API with PostgreSQL support
pip install -e "./saas[dev,pg]"

# SDK
pip install -e "./sdk[dev]"

# Frontend (Node 22)
cd frontend && npm ci

# Studio (Node 22)
cd studio && npm ci
```

## Running Tests

```bash
# Core library (1510 tests across 36 files, CI enforces --cov-fail-under=80)
python -m pytest tests/ -v
python -m pytest tests/test_smiles.py -v          # single file
python -m pytest tests/ --cov=molbuilder --cov-report=term-missing

# SaaS API (187 tests across 24 files, uses FastAPI TestClient)
cd saas && python -m pytest tests/ -v             # SQLite backend (default)

# SaaS API against PostgreSQL (requires running PG instance)
DATABASE_BACKEND=postgresql DATABASE_URL=postgresql://user:pass@localhost/molbuilder_test \
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
core/  <-- everything depends on this (constants, elements, geometry, bond data)
  +-- atomic/      (Bohr model, quantum numbers, wavefunctions)
  +-- bonding/     (Lewis structures, VSEPR, covalent bond analysis)
  +-- molecule/    (Molecule graph, properties, SA score -- central class)
  |     +-- smiles/  (SMILES tokenizer, parser, writer)
  |     +-- io/      (XYZ, MOL/SDF V2000, PDB, JSON)
  +-- coords/      (3D coordinate generation -- DG+FF builtin or RDKit backend)
  +-- dynamics/    (MD engine: force field, Verlet integrator, reaction mechanisms)
  +-- smarts/      (SMARTS pattern matching engine -- atom/bond primitives, recursive SMARTS)
  +-- data/        (ORD conditions JSON, template-to-ORD mapping)
  +-- reactions/   (185 templates, 24 FG detectors, retrosynthesis, FG SMARTS validation, RetroCast adapter)
  |     +-- process/  (reactor, conditions, condition prediction, costing, safety, scale-up)
  +-- reports/     (ASCII + PDF report generators)
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

FastAPI with Pydantic v2. Auth flow: API key (X-API-Key) -> JWT exchange. Five tiers: free/pro/team/academic/enterprise. RBAC roles: admin/chemist/viewer.

**Database:** Dual-backend via `DatabaseBackend` abstraction (`saas/app/services/database.py`). SQLite for dev/tests (default), PostgreSQL for production. Toggle via `DATABASE_BACKEND` env var (`"sqlite"` or `"postgresql"`). PostgreSQL uses psycopg v3 with `ConnectionPool`. Schema managed by Alembic (`saas/alembic/`).

**5 database tables** (consolidated from 5 separate SQLite files):
- `api_keys` (user_db) -- API key records, Stripe billing info
- `audit_log` (audit_db) -- immutable 21 CFR Part 11 audit trail with HMAC signatures
- `jobs` (job_db) -- batch job state with JSONB input/result
- `library_molecules` (library_db) -- saved molecule library with JSONB tags/properties
- `molecules` (molecule_store) -- SMILES cache with LRU eviction

**15 router modules:** auth, molecule, retrosynthesis, process, batch, library, file_io, reports, billing, analytics, audit, elements, version, feasibility, legal.

**Middleware stack:** SecurityHeaders, RequestID, Versioning, UsageTracking, CORS.

**Services:** molecule_service, retro_service, process_service, feasibility_service, prediction_service, file_io_service, batch_worker, job_db, library_db, molecule_store, user_db, audit_db, stripe_service, usage_tracker.

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

## Database Gotchas

- **Catch `DatabaseIntegrityError`, NOT `sqlite3.IntegrityError`** -- the abstraction layer re-raises driver-specific errors as `DatabaseIntegrityError` (`app.services.database`).
- All DB service classes (`UserDB`, `AuditDB`, `JobDB`, `LibraryDB`, `MoleculeStore`) accept an optional `backend: DatabaseBackend` param. When omitted they fall back to direct SQLite.
- **Audit HMAC signatures** are always computed from the epoch float. PostgreSQL stores both `timestamp` (TIMESTAMPTZ) and `timestamp_epoch` (FLOAT). The `_normalise_record()` method ensures `record["timestamp"]` is always the float for signature verification.
- **JSON fields:** Use `supports_native_json` property to decide between `json.dumps()`/`json.loads()` (SQLite TEXT) and passing dicts directly (PostgreSQL JSONB).
- **Upsert dialect:** SQLite uses `INSERT OR REPLACE`, PostgreSQL uses `INSERT ... ON CONFLICT DO UPDATE`. MoleculeStore handles this in the `put()` method.
- **Tag search dialect:** LibraryDB uses `tags LIKE '%"tag"%'` for SQLite, `tags @> ?::jsonb` for PostgreSQL.
- **SQL placeholders:** All service code uses `?` placeholders. The PostgresBackend auto-translates to `%s` via `_translate_sql()`.
- **Alembic migrations:** Run `cd saas && alembic upgrade head` before first use with PostgreSQL. Schema lives in `saas/alembic/versions/`.

## Encoding Rules

All `.py` files MUST be Windows cp1252 compatible -- no emojis, no unicode math symbols. Use ASCII only: `->` not right-arrow, `deg` not degree symbol, spell out Greek letters.

## Key Import Patterns

```python
# Core library
from molbuilder.molecule.graph import Molecule, Atom, Bond, Hybridization
from molbuilder.smiles import parse, to_smiles
from molbuilder.core.constants import PLANCK_CONSTANT, SPEED_OF_LIGHT
from molbuilder.core.elements import ELEMENTS, SYMBOL_TO_Z, from_symbol
from molbuilder.io import write_xyz, read_xyz, write_mol, read_mol
from molbuilder.coords import generate_3d              # 3D coord generation (auto/builtin/rdkit)
from molbuilder.dynamics import MDSimulation, ForceField  # molecular dynamics
from molbuilder.smarts import SmartsMatcher            # SMARTS substructure search
from molbuilder.reactions import REACTION_TEMPLATES, detect_functional_groups
from molbuilder.reactions.retrosynthesis import retrosynthesis, format_tree
from molbuilder.reactions.fg_smarts_validation import cross_validate_fg  # FG cross-validation
from molbuilder.process.reactor import select_reactor
from molbuilder.process.costing import estimate_cost
from molbuilder.process.condition_prediction import predict_conditions  # ORD-backed
from molbuilder.molecule.properties import lipinski_properties
from molbuilder.molecule.sa_score import sa_score       # synthetic accessibility (1-10)
from molbuilder.reports.pdf_report import generate_molecule_pdf  # optional dep

# SaaS database layer
from app.services.database import DatabaseBackend, SQLiteBackend, PostgresBackend
from app.services.database import DatabaseIntegrityError, get_backend, set_backend
from app.services.user_db import UserDB, get_user_db, set_user_db
from app.services.audit_db import AuditDB, get_audit_db, set_audit_db
from app.services.job_db import JobDB, get_job_db, set_job_db
from app.services.library_db import LibraryDB, get_library_db, set_library_db
from app.services.molecule_store import MoleculeStore, get_molecule_store, set_molecule_store
```

## Entry Points

- `python -m molbuilder` -- interactive CLI menu
- `python synthesize.py` -- SMILES-to-synthesis pipeline
- `uvicorn app.main:app` (from saas/) -- API server
- `npm run dev` (from frontend/) -- dashboard dev server, port 5173
- `npm run dev` (from studio/) -- 3D editor dev server, port 5174
- `cd saas && alembic upgrade head` -- run PostgreSQL schema migrations
- `cd saas && python scripts/migrate_sqlite_to_pg.py` -- one-time SQLite-to-PG data migration
- `cd saas && python scripts/verify_migration.py` -- post-migration row count + signature verification

## CI/CD

- **test.yml**: Python 3.11/3.12/3.13 matrix. Ruff lint, core tests (coverage >= 80%), saas tests (SQLite), sdk tests.
- **test.yml (test-pg job)**: Optional PostgreSQL 16 integration tests. Runs on push to main or PRs labeled `test-pg`. Runs Alembic migrations then full SaaS suite against real PG.
- **frontend-test.yml / studio-test.yml**: Node 22. Lint, typecheck, build.
- **deploy.yml**: Auto-deploys SaaS to Railway on push to main (when saas/ or molbuilder/ change).
- **release-sdk.yml**: Publishes SDK to PyPI on `sdk-v*` tags.

## Reference

See `DEVELOPMENT_GUIDE.md` for detailed architecture, completed sprint history, and verification commands.
