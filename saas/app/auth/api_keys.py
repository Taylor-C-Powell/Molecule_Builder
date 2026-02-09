"""API key generation, hashing, and validation."""

import hashlib
import secrets
import threading
from app.config import Tier
from app.auth.roles import Role


_PREFIX = {
    Tier.FREE: "mb_free_",
    Tier.PRO: "mb_pro_",
    Tier.TEAM: "mb_team_",
    Tier.ACADEMIC: "mb_acad_",
    Tier.ENTERPRISE: "mb_ent_",
}


class APIKeyRecord:
    __slots__ = ("key_hash", "email", "tier", "role")

    def __init__(self, key_hash: str, email: str, tier: Tier, role: Role = Role.CHEMIST):
        self.key_hash = key_hash
        self.email = email
        self.tier = tier
        self.role = role


class APIKeyStore:
    """Thread-safe in-memory API key store."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._keys: dict[str, APIKeyRecord] = {}  # hash -> record

    @staticmethod
    def _hash(key: str) -> str:
        return hashlib.sha256(key.encode()).hexdigest()

    def create(self, email: str, tier: Tier = Tier.FREE, role: Role = Role.CHEMIST) -> str:
        prefix = _PREFIX[tier]
        raw_key = prefix + secrets.token_urlsafe(32)
        h = self._hash(raw_key)
        record = APIKeyRecord(key_hash=h, email=email, tier=tier, role=role)
        with self._lock:
            self._keys[h] = record
        return raw_key

    def validate(self, key: str) -> APIKeyRecord | None:
        h = self._hash(key)
        with self._lock:
            return self._keys.get(h)

    def find_by_email(self, email: str) -> list[APIKeyRecord]:
        """Return all key records for a given email."""
        with self._lock:
            return [r for r in self._keys.values() if r.email == email]

    def delete_by_email(self, email: str) -> int:
        """Delete all keys for a given email. Returns number deleted."""
        with self._lock:
            to_delete = [h for h, r in self._keys.items() if r.email == email]
            for h in to_delete:
                del self._keys[h]
            return len(to_delete)

    def list_users(self) -> list[dict]:
        """Return a list of unique users with their tier and role."""
        with self._lock:
            seen: dict[str, dict] = {}
            for r in self._keys.values():
                if r.email not in seen:
                    seen[r.email] = {
                        "email": r.email,
                        "tier": r.tier.value,
                        "role": r.role.value,
                    }
            return list(seen.values())


api_key_store = APIKeyStore()
