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

# Verify molbuilder imports work at build time
RUN python -c "from molbuilder.smiles.parser import parse; parse('C'); print('molbuilder import OK')"

ENV AUDIT_DB_PATH=/data/molbuilder_audit.db
ENV USER_DB_PATH=/data/molbuilder_users.db

ENV PORT=8000
EXPOSE ${PORT}

HEALTHCHECK --interval=30s --timeout=10s --retries=3 \
    CMD python -c "import os,urllib.request; urllib.request.urlopen(f'http://localhost:{os.environ.get(\"PORT\",8000)}/health')" || exit 1

CMD ["sh", "-c", "uvicorn app.main:app --host 0.0.0.0 --port $PORT --log-level info"]
