"""Role-based access control: roles and permission checking."""

from enum import Enum


class Role(str, Enum):
    ADMIN = "admin"
    CHEMIST = "chemist"
    VIEWER = "viewer"


# Viewer: GET only; Chemist: GET + POST computation; Admin: all access
_VIEWER_METHODS = {"GET"}
_CHEMIST_METHODS = {"GET", "POST"}

# Paths restricted to admin regardless of method
_ADMIN_ONLY_PREFIXES = (
    "/api/v1/analytics",
    "/api/v1/audit",
    "/api/v1/auth/users",
)

# Self-service paths: any authenticated user, any method
# These endpoints are user-scoped (operations on the user's own data)
_SELF_SERVICE_PREFIXES = (
    "/api/v1/auth/me",
    "/api/v1/library",
    "/api/v1/batch",
)


def check_permission(role: Role, method: str, path: str) -> bool:
    """Return True if the given role is allowed to perform method on path."""
    if role == Role.ADMIN:
        return True

    # Self-service endpoints are open to all authenticated users
    for prefix in _SELF_SERVICE_PREFIXES:
        if path.startswith(prefix):
            return True

    # Admin-only paths
    for prefix in _ADMIN_ONLY_PREFIXES:
        if path.startswith(prefix):
            return False

    if role == Role.CHEMIST:
        return method.upper() in _CHEMIST_METHODS

    if role == Role.VIEWER:
        return method.upper() in _VIEWER_METHODS

    return False
