# MolBuilder Launch Posts

Drafted 2026-02-10. Submit to each platform manually.

---

## Hacker News — Show HN

**Submit at:** https://news.ycombinator.com/submit
**URL field:** https://github.com/Taylor-C-Powell/Molecule_Builder
**Title:** Show HN: MolBuilder – Open-source molecular engineering toolkit with REST API

**Comment (post after submitting the link):**

Hi HN, I've been building MolBuilder for the past year -- a Python toolkit for molecular engineering that just got a REST API.

Core features:
- SMILES parsing with 3D coordinate generation
- Retrosynthesis planning (181 reaction templates, beam search algorithm)
- Process engineering evaluation (cost, yield, scale-up analysis)
- Molecular dynamics engine (Velocity Verlet integrator)
- File I/O for MOL, PDB, XYZ, SDF formats

The interesting technical bits: It's pure Python with NumPy/SciPy for the math. The retrosynthesis engine uses a beam search over disconnection transforms, similar to how AlphaZero-based approaches work but much simpler. The API layer adds RBAC, audit trails (21 CFR Part 11 compliance for pharma), and tiered rate limiting.

I'm a solo developer and this started as a way to learn computational chemistry. There are definitely gaps -- the forcefield is basic, the reaction templates are hand-coded, and it's not trying to replace RDKit or other mature tools. But it scratches an itch for programmatic access to synthesis planning without enterprise licensing.

Free tier: 10 req/min. MIT licensed. 663 tests for core package, 107 for API.

- PyPI: `pip install molbuilder`
- GitHub: https://github.com/Taylor-C-Powell/Molecule_Builder
- API: https://www.molbuilder.io
- Docs: https://www.molbuilder.io/docs

Happy to answer questions about the implementation or design choices.

---

## Reddit r/chemistry

**Submit at:** https://www.reddit.com/r/chemistry/submit
**Title:** Open-source toolkit for retrosynthesis planning and molecular engineering -- just added a REST API

**Body:**

I've been working on MolBuilder, a Python package for computational chemistry that might be useful for synthesis planning or teaching.

What it does:
- Parse SMILES and generate 3D structures
- Retrosynthesis planning with 181 reaction templates (beam search)
- Process evaluation: cost estimates, yield predictions, scale-up analysis
- Functional group detection and transformations
- Export to MOL, PDB, XYZ, SDF formats
- Generate GHS safety reports

The retrosynthesis engine is the main feature -- you give it a target molecule and it works backward through reaction templates to find synthesis routes. It ranks routes by cost and feasibility. Not as sophisticated as commercial tools, but it's free and you can extend the template library.

Just deployed a REST API so you can use it without installing Python. Free tier is 10 requests/minute.

Install: `pip install molbuilder`
GitHub: https://github.com/Taylor-C-Powell/Molecule_Builder
API: https://www.molbuilder.io
Docs: https://www.molbuilder.io/docs

MIT licensed. Would love feedback from practicing chemists on what would make this more useful.

---

## Reddit r/cheminformatics

**Submit at:** https://www.reddit.com/r/cheminformatics/submit
**Title:** MolBuilder: Retrosynthesis and process engineering toolkit with REST API (not an RDKit replacement)

**Body:**

I've built an open-source cheminformatics toolkit focused on synthesis planning and process engineering. Just deployed the REST API.

Key capabilities:
- SMILES parsing with 3D coordinate generation
- Retrosynthesis: 181 reaction templates, beam search with cost/feasibility scoring
- Process engineering: yield optimization, cost estimation, scale-up analysis
- Functional group detection (58 SMARTS patterns)
- File I/O: MOL, PDB, XYZ, SDF
- Molecular dynamics with configurable forcefields

Technical approach:
- Pure Python with NumPy/SciPy for numerical work
- Graph-based molecule representation
- Disconnection-transform retrosynthesis (inspired by Corey's logic)
- OPLS-AA parameters for dynamics
- 663 unit tests (core package), 107 integration tests (API)

API features:
- RBAC with role-based access control
- 21 CFR Part 11 audit trails for regulated environments
- Tiered rate limiting (10 req/min free tier)
- OpenAPI/Swagger documentation

Positioning: This is NOT trying to replace RDKit or other established toolkits. It's complementary -- think of it as a higher-level toolkit for synthesis planning workflows. If you need advanced fingerprinting, substructure search, or property prediction, use RDKit. If you need to plan a synthesis route and estimate costs, this might help.

Solo project, MIT licensed.

PyPI: `pip install molbuilder`
GitHub: https://github.com/Taylor-C-Powell/Molecule_Builder
API: https://www.molbuilder.io
Docs: https://www.molbuilder.io/docs

Questions about the implementation or design decisions welcome.

---

## Reddit r/Python

**Submit at:** https://www.reddit.com/r/Python/submit
**Title:** Built a molecular engineering REST API in pure Python -- 663 tests, RBAC, audit trails

**Body:**

I've spent the past year building MolBuilder, an open-source computational chemistry toolkit that just got a production API. Sharing because the engineering might be interesting to this community.

The stack:
- Core: Pure Python with NumPy/SciPy (no C extensions)
- API: FastAPI with Pydantic validation
- Auth: JWT + API key RBAC with role hierarchies
- Storage: SQLite with WAL mode
- Testing: 663 unit tests (pytest), 107 integration tests for API
- Deployment: Docker on Railway with environment-based config

Technical features:
- Graph-based molecule representation with adjacency lists
- Beam search for retrosynthesis planning (181 reaction templates)
- Velocity Verlet integrator for molecular dynamics
- Comprehensive file I/O (MOL, PDB, XYZ, SDF parsers)
- 21 CFR Part 11 audit trail implementation for pharma compliance
- Rate limiting with tiered quotas

Engineering decisions:
- Stuck with pure Python for portability (no compilation needed)
- Dataclasses everywhere for structure
- OpenAPI spec autogenerated from FastAPI
- Modular package design: core, bonding, molecule, smiles, io, reactions, dynamics, cli, gui

The API is free for personal use (10 req/min). MIT licensed so you can fork and extend.

Install: `pip install molbuilder`
GitHub: https://github.com/Taylor-C-Powell/Molecule_Builder
API: https://www.molbuilder.io
Swagger docs: https://www.molbuilder.io/docs

Happy to discuss the architecture or Python-specific implementation choices.
