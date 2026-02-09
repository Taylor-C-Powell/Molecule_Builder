# Operational Qualification (OQ) Protocol

**Document ID**: MB-OQ-001
**System**: MolBuilder SaaS API v0.1.0
**Effective Date**: _______________
**Prepared By**: _______________
**Approved By**: _______________

## 1. Purpose

This protocol verifies that the MolBuilder SaaS API operates correctly within
specified operating ranges per 21 CFR Part 11 and ICH Q2(R1) guidelines. OQ
testing confirms that all features function as designed.

## 2. Prerequisites

- Installation Qualification (MB-IQ-001) completed and approved
- Test environment configured with Python >= 3.11
- pytest and httpx installed (`pip install ".[dev]"`)

## 3. Test Suite Execution

### 3.1 Core Library Tests

```bash
cd Molecule_Builder
python -m pytest tests/ -v --tb=short
```

**Expected**: All tests pass (660+ tests)

| Step | Expected Result | Actual Result | Pass/Fail | Initials | Date |
|------|----------------|---------------|-----------|----------|------|
| 3.1 Core tests: all pass | 660+ passed | | | | |

### 3.2 SaaS API Tests

```bash
cd Molecule_Builder/saas
python -m pytest tests/ -v --tb=short
```

**Expected**: All tests pass (79 tests)

| Step | Expected Result | Actual Result | Pass/Fail | Initials | Date |
|------|----------------|---------------|-----------|----------|------|
| 3.2 SaaS tests: all pass | 79 passed | | | | |

### 3.3 Drug Validation Suite

```bash
cd Molecule_Builder
python -m pytest tests/test_drug_validation.py -v
```

**Expected**: 50 drug molecules processed, all pass or xfail

| Step | Expected Result | Actual Result | Pass/Fail | Initials | Date |
|------|----------------|---------------|-----------|----------|------|
| 3.3 Drug validation: 50 pass | 50 passed | | | | |

### 3.4 Validation Report Generation

```bash
cd Molecule_Builder
python tests/validation_report.py
```

**Expected**: Report generated at `docs/validation/drug_validation_results.md` with >70% route coverage

| Step | Expected Result | Actual Result | Pass/Fail | Initials | Date |
|------|----------------|---------------|-----------|----------|------|
| 3.4 Report generated, >70% coverage | | | | | |

## 4. Functional Verification

### 4.1 Authentication

| Test | Procedure | Expected | Pass/Fail | Initials | Date |
|------|-----------|----------|-----------|----------|------|
| API key registration | POST /api/v1/auth/register | Returns API key | | | |
| Token exchange | POST /api/v1/auth/token | Returns JWT | | | |
| Unauthorized access | GET /api/v1/elements/C (no auth) | 401 | | | |
| Invalid API key | GET with bad X-API-Key | 401 | | | |
| Rate limiting | Exceed tier RPM limit | 429 | | | |

### 4.2 Role-Based Access Control

| Test | Procedure | Expected | Pass/Fail | Initials | Date |
|------|-----------|----------|-----------|----------|------|
| Admin: full access | Admin key on all endpoints | 200 | | | |
| Chemist: GET+POST | Chemist key on computation | 200 | | | |
| Viewer: GET only | Viewer key on POST | 403 | | | |
| Admin: analytics | Admin key on /analytics | 200 | | | |
| Non-admin: analytics | Chemist key on /analytics | 403 | | | |

### 4.3 Audit Trail

| Test | Procedure | Expected | Pass/Fail | Initials | Date |
|------|-----------|----------|-----------|----------|------|
| Record creation | Make authenticated request | Audit record created | | | |
| Signature integrity | GET /audit/log/{id}/verify | is_valid: true | | | |
| Tamper detection | Modify DB record, verify | is_valid: false | | | |
| Admin query | GET /audit/log (admin) | Returns records | | | |
| Non-admin denied | GET /audit/log (chemist) | 403 | | | |
| Export | GET /audit/export | JSON array | | | |

### 4.4 Molecular Operations

| Test | Procedure | Expected | Pass/Fail | Initials | Date |
|------|-----------|----------|-----------|----------|------|
| SMILES parsing | POST /molecule/from-smiles (ethanol) | Molecule data | | | |
| Properties | GET /molecule/{id}/properties | MW, formula | | | |
| 3D structure | GET /molecule/{id}/3d | Coordinates | | | |
| Retrosynthesis | POST /retrosynthesis/plan | Route tree | | | |
| Process evaluation | POST /process/evaluate | Process metrics | | | |

### 4.5 API Versioning

| Test | Procedure | Expected | Pass/Fail | Initials | Date |
|------|-----------|----------|-----------|----------|------|
| Version header | Any request | X-API-Version: 1.0.0 | | | |
| Version endpoint | GET /api/v1/version | Version info + policy | | | |
| Deprecation warning | Send X-API-Version: 0.0.1 | X-Deprecation-Warning | | | |

## 5. Reaction Knowledge Base Verification

| Metric | Requirement | Actual | Pass/Fail | Initials | Date |
|--------|-------------|--------|-----------|----------|------|
| Total reaction templates | >= 175 | | | | |
| Purchasable materials | >= 200 | | | | |
| Functional group detectors | >= 15 | | | | |
| Reaction categories | >= 10 | | | | |
| Named reactions findable | >= 30 | | | | |

## 6. Deviation Handling

Any deviations from expected results must be documented, investigated, and resolved
before proceeding to Performance Qualification. Deviations require sign-off from
the QA Manager.

## 7. Sign-Off

| Role | Name | Signature | Date |
|------|------|-----------|------|
| QA Manager | | | |
| Lead Developer | | | |
| Department Head | | | |
