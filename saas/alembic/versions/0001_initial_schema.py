"""Initial schema -- all 5 tables migrated from SQLite.

Revision ID: 0001
Revises: None
Create Date: 2026-03-01
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB, TIMESTAMP

revision: str = "0001"
down_revision: Union[str, None] = None
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    # --- api_keys (from user_db) ---
    op.create_table(
        "api_keys",
        sa.Column("id", sa.Integer, primary_key=True, autoincrement=True),
        sa.Column("key_hash", sa.Text, nullable=False, unique=True),
        sa.Column("email", sa.Text, nullable=False),
        sa.Column("tier", sa.Text, nullable=False, server_default="free"),
        sa.Column("role", sa.Text, nullable=False, server_default="chemist"),
        sa.Column(
            "created_at",
            TIMESTAMP(timezone=True),
            nullable=False,
            server_default=sa.func.now(),
        ),
        sa.Column("active", sa.Boolean, nullable=False, server_default=sa.text("true")),
        sa.Column("stripe_customer_id", sa.Text, nullable=True),
        sa.Column("stripe_subscription_id", sa.Text, nullable=True),
        sa.Column("subscription_status", sa.Text, server_default="none"),
    )
    op.create_index("idx_keys_email", "api_keys", ["email"])
    op.create_index("idx_keys_active", "api_keys", ["active"])

    # --- audit_log (from audit_db) ---
    op.create_table(
        "audit_log",
        sa.Column("id", sa.Integer, primary_key=True, autoincrement=True),
        sa.Column(
            "timestamp",
            TIMESTAMP(timezone=True),
            nullable=False,
            server_default=sa.func.now(),
        ),
        sa.Column("timestamp_epoch", sa.Float, nullable=False),
        sa.Column("user_email", sa.Text, nullable=False),
        sa.Column("user_role", sa.Text, nullable=False, server_default="chemist"),
        sa.Column("user_tier", sa.Text, nullable=False, server_default="free"),
        sa.Column("action", sa.Text, nullable=False),
        sa.Column("input_summary", sa.Text, nullable=False, server_default=""),
        sa.Column("output_hash", sa.Text, nullable=False, server_default=""),
        sa.Column("status_code", sa.Integer, nullable=False, server_default=sa.text("0")),
        sa.Column("latency_ms", sa.Float, nullable=False, server_default=sa.text("0")),
        sa.Column("ip_address", sa.Text, nullable=False, server_default=""),
        sa.Column("signature_hash", sa.Text, nullable=False),
    )
    op.create_index("idx_audit_timestamp", "audit_log", ["timestamp"])
    op.create_index("idx_audit_user", "audit_log", ["user_email"])
    op.create_index("idx_audit_action", "audit_log", ["action"])

    # --- jobs (from job_db) ---
    op.create_table(
        "jobs",
        sa.Column("job_id", sa.Text, primary_key=True),
        sa.Column("user_email", sa.Text, nullable=False),
        sa.Column("status", sa.Text, nullable=False, server_default="pending"),
        sa.Column("job_type", sa.Text, nullable=False),
        sa.Column("input_data", JSONB, nullable=False),
        sa.Column("result_data", JSONB, nullable=True),
        sa.Column("error", sa.Text, nullable=True),
        sa.Column(
            "created_at",
            TIMESTAMP(timezone=True),
            nullable=False,
            server_default=sa.func.now(),
        ),
        sa.Column(
            "updated_at",
            TIMESTAMP(timezone=True),
            nullable=False,
            server_default=sa.func.now(),
        ),
        sa.Column("progress_pct", sa.Float, nullable=False, server_default=sa.text("0.0")),
    )
    op.create_index("idx_jobs_user", "jobs", ["user_email"])
    op.create_index("idx_jobs_status", "jobs", ["status"])

    # --- library_molecules (from library_db) ---
    op.create_table(
        "library_molecules",
        sa.Column("id", sa.Integer, primary_key=True, autoincrement=True),
        sa.Column("user_email", sa.Text, nullable=False),
        sa.Column("smiles", sa.Text, nullable=False),
        sa.Column("name", sa.Text, nullable=True),
        sa.Column("tags", JSONB, nullable=False, server_default=sa.text("'[]'::jsonb")),
        sa.Column("notes", sa.Text, nullable=True),
        sa.Column("properties_json", JSONB, nullable=True),
        sa.Column(
            "created_at",
            TIMESTAMP(timezone=True),
            nullable=False,
            server_default=sa.func.now(),
        ),
        sa.Column(
            "updated_at",
            TIMESTAMP(timezone=True),
            nullable=False,
            server_default=sa.func.now(),
        ),
    )
    op.create_index(
        "idx_lib_user_smiles",
        "library_molecules",
        ["user_email", "smiles"],
        unique=True,
    )
    op.create_index("idx_lib_user", "library_molecules", ["user_email"])

    # --- molecules (from molecule_store) ---
    op.create_table(
        "molecules",
        sa.Column("mol_id", sa.Text, primary_key=True),
        sa.Column("smiles", sa.Text, nullable=False),
        sa.Column(
            "created_at",
            TIMESTAMP(timezone=True),
            nullable=False,
            server_default=sa.func.now(),
        ),
        sa.Column(
            "accessed_at",
            TIMESTAMP(timezone=True),
            nullable=False,
            server_default=sa.func.now(),
        ),
    )
    op.create_index("idx_molecules_accessed", "molecules", ["accessed_at"])


def downgrade() -> None:
    op.drop_table("molecules")
    op.drop_table("library_molecules")
    op.drop_table("jobs")
    op.drop_table("audit_log")
    op.drop_table("api_keys")
