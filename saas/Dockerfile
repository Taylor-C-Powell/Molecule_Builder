FROM python:3.13-slim

WORKDIR /app

# Install build deps
RUN pip install --no-cache-dir --upgrade pip setuptools wheel

# Copy the core molbuilder package first (for caching)
COPY molbuilder/ /build/molbuilder/
COPY pyproject.toml README.md /build/
RUN pip install --no-cache-dir /build

# Copy SaaS application
COPY saas/ /app/

# Install SaaS dependencies
RUN pip install --no-cache-dir \
    "fastapi>=0.109.0" \
    "uvicorn[standard]>=0.25.0" \
    "pydantic>=2.5.0" \
    "pydantic-settings>=2.1.0" \
    "pyjwt>=2.8.0" \
    "python-multipart>=0.0.6"

# Create data directory for SQLite audit DB
RUN mkdir -p /data

ENV AUDIT_DB_PATH=/data/molbuilder_audit.db
ENV JWT_SECRET_KEY=change-me-in-production

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=5s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8000/health')" || exit 1

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
