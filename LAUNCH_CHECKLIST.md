# MolBuilder Launch Checklist

Derived from MolBuilder_GTM_Strategy.md. Track progress by checking items as completed.

---

## Phase 0: Pre-Launch Preparation

### Package and Documentation
- [x] Verify all 585 tests pass: `pytest tests/ -v` (585/585 passing, 2026-02-08)
- [x] Review and clean up public API surface (added top-level exports, explicit core imports, populated io/__init__.py)
- [x] Ensure README.md has: installation, quickstart, feature overview, example code
- [x] Add CHANGELOG.md (covers v1.0.0 and v1.1.0)
- [x] Add CONTRIBUTING.md
- [x] Verify `pyproject.toml` metadata is complete (classifiers, keywords, URLs, PEP 639 license)
- [x] Test local install: verified in clean virtualenv with all subpackage imports
- [x] Test package build: `python -m build` produces molbuilder-1.1.0-py3-none-any.whl (243 KB)

### PyPI Publication
- [x] Create PyPI account (pypi.org)
- [x] Create TestPyPI account (test.pypi.org)
- [x] Upload to TestPyPI first: `twine upload --repository testpypi dist/*`
- [x] Verify TestPyPI install
- [x] Upload to PyPI: `twine upload dist/*`
- [x] Verify PyPI install: `pip install molbuilder` resolves v1.1.0 (2026-02-08)
- [ ] Upload v1.1.1 patch release with corrected URLs and API improvements

### Web Presence
- [ ] Register domain (molbuilder.io or similar)
- [ ] Create landing page with: value proposition, feature highlights, install command, screenshots
- [x] Set up GitHub repository: https://github.com/Taylor-C-Powell/Molecule_Builder
- [ ] Add GitHub topics: chemistry, molecular-modeling, retrosynthesis, cheminformatics, python
- [ ] Add GitHub repo description
- [ ] Create GitHub Discussions or link to community Discord

---

## Phase 1: Soft Launch (Months 1-3)

### API Development
- [ ] Set up FastAPI project in `saas/` directory
- [ ] Implement core endpoints:
  - [ ] `POST /api/v1/molecule/from-smiles` — parse SMILES, return molecular data
  - [ ] `GET /api/v1/molecule/{id}/properties` — molecular properties
  - [ ] `GET /api/v1/molecule/{id}/3d` — 3D coordinates for visualization
  - [ ] `POST /api/v1/retrosynthesis/plan` — run retrosynthesis on target molecule
  - [ ] `POST /api/v1/process/evaluate` — process engineering evaluation
  - [ ] `GET /api/v1/elements/{symbol}` — element data lookup
- [ ] Authentication (API keys for free tier, JWT for Pro+)
- [ ] Rate limiting per tier
- [ ] Basic error handling and input validation
- [ ] API documentation (auto-generated via FastAPI /docs)

### Web Frontend
- [ ] React project scaffold
- [ ] SMILES input component with validation
- [ ] 3D molecule viewer (Three.js or 3Dmol.js)
- [ ] Properties display panel
- [ ] Retrosynthesis tree visualization
- [ ] User account and billing pages

### Billing
- [ ] Stripe account setup
- [ ] Free tier: no payment required, rate-limited API key
- [ ] Pro tier: $49/month or $490/year via Stripe Checkout
- [ ] Team tier: $199/month per seat via Stripe Billing
- [ ] Academic verification: email domain check (.edu, .ac.uk, etc.)

### Community
- [ ] Post on Hacker News (Show HN: MolBuilder — open-source retrosynthesis in Python)
- [ ] Post on Reddit: r/chemistry, r/cheminformatics, r/Python, r/bioinformatics
- [ ] Contact 10 organic chemistry professors with free access offers
- [ ] Write launch blog post for Teapot Commons

---

## Phase 2: Growth (Months 3-6)

### Product
- [ ] Full retrosynthesis visualization in web GUI (interactive tree)
- [ ] Process engineering module in web GUI
- [ ] Cost estimation reports (PDF export)
- [ ] GHS safety reports (PDF export)
- [ ] File upload/download (MOL, PDB, XYZ, SDF)
- [ ] Molecule library (save and organize molecules)

### Marketing
- [ ] Submit paper to Journal of Cheminformatics
- [ ] Attend 1 chemistry conference (ACS, Gordon, or virtual)
- [ ] Publish 4+ technical blog posts on Teapot Commons
- [ ] Record tutorial video series
- [ ] Set up academic licensing page

### Metrics to Track
- [ ] PyPI download count (target: 5,000 by month 6)
- [ ] GitHub stars (target: 500 by month 6)
- [ ] Free registered users (target: 1,000 by month 6)
- [ ] Paying customers (target: 50 by month 6)
- [ ] ARR (target: $50K by month 6)

---

## Phase 3: Scale (Months 6-12)

### Enterprise
- [ ] SSO integration (Auth0/Okta)
- [ ] On-premise deployment option (Docker + Kubernetes)
- [ ] Custom reaction template import
- [ ] ELN/LIMS integration APIs
- [ ] Audit logging
- [ ] SLA documentation

### Sales
- [ ] Institutional licensing program live
- [ ] CRO partnership discussions initiated
- [ ] Biotech accelerator partnerships
- [ ] Enterprise sales collateral (case studies, ROI calculator)

### 12-Month Targets
- [ ] 200+ paying customers
- [ ] $300K+ ARR run rate
- [ ] 15+ academic papers citing MolBuilder
- [ ] 25,000 PyPI downloads
