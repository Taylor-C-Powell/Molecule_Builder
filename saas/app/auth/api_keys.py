"""API key generation, hashing, and validation."""

import hashlib
import hmac
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
    """Thread-safe in-memory API key store with optional SQLite persistence."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._keys: dict[str, APIKeyRecord] = {}  # hash -> record
        self._db = None

    @staticmethod
    def _hash(key: str) -> str:
        """HMAC-SHA256 hash of an API key using the server secret."""
        from app.config import settings
        secret = settings.api_key_hmac_secret
        if secret:
            return hmac.new(secret.encode(), key.encode(), hashlib.sha256).hexdigest()
        # Fallback for tests/bootstrap before secret is set
        return hashlib.sha256(key.encode()).hexdigest()

    @staticmethod
    def _legacy_hash(key: str) -> str:
        """Plain SHA256 for backward compatibility with pre-HMAC keys."""
        return hashlib.sha256(key.encode()).hexdigest()

    def load_from_db(self, db) -> None:
        """Load all active keys from a UserDB and enable write-through."""
        self._db = db
        rows = db.load_all_active()
        with self._lock:
            for row in rows:
                record = APIKeyRecord(
                    key_hash=row["key_hash"],
                    email=row["email"],
                    tier=Tier(row["tier"]),
                    role=Role(row["role"]),
                )
                self._keys[row["key_hash"]] = record

    def create(self, email: str, tier: Tier = Tier.FREE, role: Role = Role.CHEMIST) -> str:
        prefix = _PREFIX[tier]
        raw_key = prefix + secrets.token_urlsafe(32)
        h = self._hash(raw_key)
        record = APIKeyRecord(key_hash=h, email=email, tier=tier, role=role)
        with self._lock:
            self._keys[h] = record
        if self._db is not None:
            self._db.insert_key(h, email, tier.value, role.value)
        return raw_key

    def validate(self, key: str) -> APIKeyRecord | None:
        h = self._hash(key)
        with self._lock:
            record = self._keys.get(h)
            if record is not None:
                return record
            # Backward compatibility: try legacy hash for pre-HMAC keys
            legacy_h = self._legacy_hash(key)
            return self._keys.get(legacy_h)

    def revoke_key(self, key: str) -> bool:
        """Revoke a specific API key. Returns True if found and revoked."""
        h = self._hash(key)
        with self._lock:
            if h in self._keys:
                record = self._keys.pop(h)
                if self._db is not None:
                    self._db.delete_by_hash(h)
                return True
            # Try legacy hash
            legacy_h = self._legacy_hash(key)
            if legacy_h in self._keys:
                self._keys.pop(legacy_h)
                if self._db is not None:
                    self._db.delete_by_hash(legacy_h)
                return True
        return False

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
        if self._db is not None and to_delete:
            self._db.delete_by_email(email)
        return len(to_delete)

    def update_tier(self, email: str, new_tier: Tier) -> int:
        """Update tier for all in-memory keys for an email, write through to DB."""
        count = 0
        with self._lock:
            for r in self._keys.values():
                if r.email == email:
                    r.tier = new_tier
                    count += 1
        if self._db is not None and count > 0:
            self._db.update_tier(email, new_tier.value)
        return count

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
