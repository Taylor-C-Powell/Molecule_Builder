"""Team management -- teams, members, and team library tables.

Revision ID: 0002
Revises: 0001
Create Date: 2026-03-02
"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB, TIMESTAMP

revision: str = "0002"
down_revision: Union[str, None] = "0001"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    # --- teams ---
    op.create_table(
        "teams",
        sa.Column("id", sa.Integer, primary_key=True, autoincrement=True),
        sa.Column("name", sa.Text, nullable=False),
        sa.Column("slug", sa.Text, nullable=False, unique=True),
        sa.Column("owner_email", sa.Text, nullable=False),
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
    op.create_index("idx_teams_owner", "teams", ["owner_email"])

    # --- team_members ---
    op.create_table(
        "team_members",
        sa.Column("id", sa.Integer, primary_key=True, autoincrement=True),
        sa.Column(
            "team_id", sa.Integer,
            sa.ForeignKey("teams.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column("user_email", sa.Text, nullable=False),
        sa.Column("team_role", sa.Text, nullable=False, server_default="member"),
        sa.Column(
            "joined_at",
            TIMESTAMP(timezone=True),
            nullable=False,
            server_default=sa.func.now(),
        ),
    )
    op.create_index(
        "idx_tm_team_user",
        "team_members",
        ["team_id", "user_email"],
        unique=True,
    )

    # --- team_library_molecules ---
    op.create_table(
        "team_library_molecules",
        sa.Column("id", sa.Integer, primary_key=True, autoincrement=True),
        sa.Column(
            "team_id", sa.Integer,
            sa.ForeignKey("teams.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column("smiles", sa.Text, nullable=False),
        sa.Column("name", sa.Text, nullable=True),
        sa.Column("tags", JSONB, nullable=False, server_default=sa.text("'[]'::jsonb")),
        sa.Column("notes", sa.Text, nullable=True),
        sa.Column("properties_json", JSONB, nullable=True),
        sa.Column("added_by", sa.Text, nullable=False),
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
        "idx_tlib_team_smiles",
        "team_library_molecules",
        ["team_id", "smiles"],
        unique=True,
    )


def downgrade() -> None:
    op.drop_table("team_library_molecules")
    op.drop_table("team_members")
    op.drop_table("teams")
