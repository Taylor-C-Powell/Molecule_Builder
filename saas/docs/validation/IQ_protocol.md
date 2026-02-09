# Installation Qualification (IQ) Protocol

**Document ID**: MB-IQ-001
**System**: MolBuilder SaaS API v0.1.0
**Effective Date**: _______________
**Prepared By**: _______________
**Approved By**: _______________

## 1. Purpose

This protocol establishes the Installation Qualification for the MolBuilder SaaS API,
verifying that the system is installed correctly and meets all specified requirements
per 21 CFR Part 11 and ICH Q2(R1) guidelines.

## 2. Scope

Covers installation of the MolBuilder SaaS API including:
- Python runtime environment
- MolBuilder core library (molbuilder)
- FastAPI web framework and dependencies
- SQLite audit database
- Docker containerization (optional)

## 3. System Requirements

| Component | Requirement |
|-----------|-------------|
| Python | >= 3.11 |
| OS | Linux (production), Windows/macOS (development) |
| RAM | >= 2 GB |
| Disk | >= 1 GB (including SQLite audit DB growth) |
| Network | Port 8000 (configurable) |

### 3.1 Python Dependencies

| Package | Minimum Version |
|---------|----------------|
| molbuilder | 1.1.1 |
| fastapi | 0.109.0 |
| uvicorn | 0.25.0 |
| pydantic | 2.5.0 |
| pydantic-settings | 2.1.0 |
| pyjwt | 2.8.0 |
| python-multipart | 0.0.6 |

## 4. Installation Procedure

### 4.1 Standard Installation

```bash
# Clone repository
git clone https://github.com/Taylor-C-Powell/Molecule_Builder.git
cd Molecule_Builder

# Install core library
pip install .

# Install SaaS dependencies
cd saas
pip install -r requirements.txt  # or: pip install ".[dev]"
```

### 4.2 Docker Installation

```bash
cd Molecule_Builder/saas
docker compose up --build
```

## 5. Verification Steps

### 5.1 Python Version Check

```bash
python --version
# Expected: Python 3.11+
```

| Step | Expected Result | Actual Result | Pass/Fail | Initials | Date |
|------|----------------|---------------|-----------|----------|------|
| 5.1 Python version >= 3.11 | | | | | |

### 5.2 Core Library Installation

```bash
python -c "import molbuilder; print(molbuilder.__version__)"
# Expected: 1.1.1 or higher
```

| Step | Expected Result | Actual Result | Pass/Fail | Initials | Date |
|------|----------------|---------------|-----------|----------|------|
| 5.2 MolBuilder version >= 1.1.1 | | | | | |

### 5.3 Dependency Verification

```bash
pip list | grep -E "fastapi|uvicorn|pydantic|pyjwt"
```

| Step | Expected Result | Actual Result | Pass/Fail | Initials | Date |
|------|----------------|---------------|-----------|----------|------|
| 5.3 All dependencies installed | | | | | |

### 5.4 Application Startup

```bash
cd saas
uvicorn app.main:app --host 0.0.0.0 --port 8000 &
curl http://localhost:8000/health
# Expected: {"status":"ok","version":"1.1.1"}
```

| Step | Expected Result | Actual Result | Pass/Fail | Initials | Date |
|------|----------------|---------------|-----------|----------|------|
| 5.4 Health endpoint returns 200 | | | | | |

### 5.5 Audit Database Initialization

```bash
python -c "
from app.services.audit_db import AuditDB
db = AuditDB(':memory:')
print('Audit DB initialized:', db.count() == 0)
"
# Expected: Audit DB initialized: True
```

| Step | Expected Result | Actual Result | Pass/Fail | Initials | Date |
|------|----------------|---------------|-----------|----------|------|
| 5.5 Audit DB initializes | | | | | |

### 5.6 API Version Header

```bash
curl -I http://localhost:8000/health
# Expected: X-API-Version: 1.0.0
```

| Step | Expected Result | Actual Result | Pass/Fail | Initials | Date |
|------|----------------|---------------|-----------|----------|------|
| 5.6 Version header present | | | | | |

## 6. Deviation Handling

Any deviations from expected results must be documented, investigated, and resolved
before proceeding to Operational Qualification. Deviations require sign-off from
the QA Manager.

## 7. Sign-Off

| Role | Name | Signature | Date |
|------|------|-----------|------|
| QA Manager | | | |
| System Administrator | | | |
| Department Head | | | |
