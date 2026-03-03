"""Team-level role definitions and permission helpers."""

from enum import Enum


class TeamRole(str, Enum):
    OWNER = "owner"
    ADMIN = "admin"
    MEMBER = "member"


def can_manage_members(role: TeamRole) -> bool:
    """Return True if the role can invite/remove members and change roles."""
    return role in (TeamRole.OWNER, TeamRole.ADMIN)


def can_delete_team(role: TeamRole) -> bool:
    """Return True if the role can delete the team."""
    return role == TeamRole.OWNER
