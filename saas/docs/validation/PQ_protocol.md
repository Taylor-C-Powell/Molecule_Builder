# Performance Qualification (PQ) Protocol

**Document ID**: MB-PQ-001
**System**: MolBuilder SaaS API v0.1.0
**Effective Date**: _______________
**Prepared By**: _______________
**Approved By**: _______________

## 1. Purpose

This protocol verifies that the MolBuilder SaaS API performs consistently and
reliably under expected operating conditions per 21 CFR Part 11 guidelines.
PQ testing confirms the system meets performance requirements in production-like
conditions.

## 2. Prerequisites

- Installation Qualification (MB-IQ-001) completed and approved
- Operational Qualification (MB-OQ-001) completed and approved
- Production-equivalent hardware or Docker environment

## 3. Performance Requirements

| Metric | Requirement |
|--------|-------------|
| Health check response | < 100 ms |
| SMILES parsing (simple) | < 500 ms |
| SMILES parsing (complex drug) | < 2,000 ms |
| Retrosynthesis (simple) | < 5,000 ms |
| Retrosynthesis (complex drug) | < 30,000 ms |
| Element lookup | < 200 ms |
| API startup time | < 10 s |
| Memory usage (idle) | < 512 MB |
| Memory usage (under load) | < 2 GB |
| Concurrent users supported | >= 10 |

## 4. Test Procedures

### 4.1 Response Time Testing

Start the API server and measure response times for key endpoints.

```bash
# Start server
cd Molecule_Builder/saas
uvicorn app.main:app --host 0.0.0.0 --port 8000 &
sleep 5

# Register and get API key
KEY=$(curl -s -X POST http://localhost:8000/api/v1/auth/register \
  -H "Content-Type: application/json" \
  -d '{"email":"perf@test.com","tier":"pro"}' | python -c "import sys,json; print(json.load(sys.stdin)['api_key'])")

# Health check latency (10 iterations)
for i in $(seq 1 10); do
  curl -o /dev/null -s -w "%{time_total}\n" http://localhost:8000/health
done

# SMILES parsing latency
curl -o /dev/null -s -w "%{time_total}\n" \
  -X POST http://localhost:8000/api/v1/molecule/from-smiles \
  -H "X-API-Key: $KEY" \
  -H "Content-Type: application/json" \
  -d '{"smiles":"CCO"}'

# Complex drug parsing
curl -o /dev/null -s -w "%{time_total}\n" \
  -X POST http://localhost:8000/api/v1/molecule/from-smiles \
  -H "X-API-Key: $KEY" \
  -H "Content-Type: application/json" \
  -d '{"smiles":"CC(C)c1c(C(=O)Nc2ccccc2)c(c(c3ccc(F)cc3)n1CC(O)CC(O)CC(O)=O)c4ccccc4"}'
```

| Test | Requirement | Measured | Pass/Fail | Initials | Date |
|------|-------------|----------|-----------|----------|------|
| Health check (avg of 10) | < 100 ms | | | | |
| Simple SMILES parse | < 500 ms | | | | |
| Complex drug parse | < 2,000 ms | | | | |
| Element lookup | < 200 ms | | | | |
| Retrosynthesis (simple) | < 5,000 ms | | | | |
| Retrosynthesis (complex) | < 30,000 ms | | | | |

### 4.2 Memory Usage Testing

```bash
# Measure idle memory
ps aux | grep uvicorn | grep -v grep | awk '{print $6/1024 " MB"}'

# Run 100 SMILES parse requests, then measure
for i in $(seq 1 100); do
  curl -s -X POST http://localhost:8000/api/v1/molecule/from-smiles \
    -H "X-API-Key: $KEY" \
    -H "Content-Type: application/json" \
    -d '{"smiles":"CCO"}' > /dev/null
done
ps aux | grep uvicorn | grep -v grep | awk '{print $6/1024 " MB"}'
```

| Test | Requirement | Measured | Pass/Fail | Initials | Date |
|------|-------------|----------|-----------|----------|------|
| Idle memory | < 512 MB | | | | |
| After 100 requests | < 2 GB | | | | |

### 4.3 Concurrent User Testing

Test with 10 concurrent connections using a simple load script:

```bash
# Run 10 concurrent health checks
for i in $(seq 1 10); do
  curl -s -o /dev/null -w "%{http_code}\n" http://localhost:8000/health &
done
wait
# Expected: all return 200
```

| Test | Requirement | Measured | Pass/Fail | Initials | Date |
|------|-------------|----------|-----------|----------|------|
| 10 concurrent health checks | All 200 | | | | |
| 10 concurrent SMILES parses | All 200 | | | | |

### 4.4 Audit Database Performance

```bash
# Verify audit records accumulate correctly
BEFORE=$(curl -s http://localhost:8000/api/v1/audit/log \
  -H "X-API-Key: $ADMIN_KEY" | python -c "import sys,json; print(json.load(sys.stdin)['total'])")

# Make 50 requests
for i in $(seq 1 50); do
  curl -s http://localhost:8000/health -H "X-API-Key: $KEY" > /dev/null
done

AFTER=$(curl -s http://localhost:8000/api/v1/audit/log \
  -H "X-API-Key: $ADMIN_KEY" | python -c "import sys,json; print(json.load(sys.stdin)['total'])")

echo "Records added: $((AFTER - BEFORE))"
# Expected: 50 new records
```

| Test | Requirement | Measured | Pass/Fail | Initials | Date |
|------|-------------|----------|-----------|----------|------|
| Audit records after 50 requests | 50 new records | | | | |
| Audit query time (1000 records) | < 500 ms | | | | |

### 4.5 Startup Time

```bash
time (uvicorn app.main:app --host 0.0.0.0 --port 8001 &
  sleep 1
  until curl -s http://localhost:8001/health > /dev/null 2>&1; do sleep 0.5; done
  echo "Ready")
```

| Test | Requirement | Measured | Pass/Fail | Initials | Date |
|------|-------------|----------|-----------|----------|------|
| API startup to first healthy response | < 10 s | | | | |

## 5. Stability Testing

Run the API for a sustained period and verify no degradation:

```bash
# Run for 1 hour, making requests every 5 seconds
for i in $(seq 1 720); do
  curl -s -o /dev/null -w "%{time_total}\n" http://localhost:8000/health
  sleep 5
done | awk '{sum+=$1; n++} END {print "Avg:", sum/n*1000, "ms"}'
```

| Test | Requirement | Measured | Pass/Fail | Initials | Date |
|------|-------------|----------|-----------|----------|------|
| 1-hour avg response time | < 100 ms | | | | |
| No 5xx errors in 1 hour | 0 errors | | | | |
| Memory growth in 1 hour | < 50 MB | | | | |

## 6. Docker Performance (if applicable)

```bash
cd Molecule_Builder/saas
docker compose up --build -d
sleep 10
curl http://localhost:8000/health
```

| Test | Requirement | Measured | Pass/Fail | Initials | Date |
|------|-------------|----------|-----------|----------|------|
| Docker container starts | Health returns 200 | | | | |
| Docker memory usage | < 1 GB | | | | |

## 7. Deviation Handling

Any deviations from performance requirements must be documented with root cause
analysis. Acceptable deviations must be approved by QA Manager and documented in
a deviation report.

## 8. Sign-Off

| Role | Name | Signature | Date |
|------|------|-----------|------|
| QA Manager | | | |
| Performance Engineer | | | |
| Department Head | | | |
