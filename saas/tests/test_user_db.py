"""Tests for API key persistence layer."""

import hashlib
import pytest
from app.services.user_db import UserDB
from app.auth.api_keys import APIKeyStore
from app.config import Tier
from app.auth.roles import Role


def _hash(key: str) -> str:
    """Compute hash using the same method as APIKeyStore."""
    return APIKeyStore._hash(key)


class TestUserDB:
    """Direct UserDB tests."""

    def test_insert_and_load(self, tmp_path):
        db = UserDB(str(tmp_path / "test.db"))
        db.insert_key("hash1", "alice@example.com", "free", "chemist")
        db.insert_key("hash2", "bob@example.com", "pro", "admin")
        rows = db.load_all_active()
        assert len(rows) == 2
        emails = {r["email"] for r in rows}
        assert emails == {"alice@example.com", "bob@example.com"}
        db.close()

    def test_soft_delete_excludes_from_active(self, tmp_path):
        db = UserDB(str(tmp_path / "test.db"))
        db.insert_key("hash1", "alice@example.com", "free", "chemist")
        db.insert_key("hash2", "alice@example.com", "pro", "chemist")
        deleted = db.delete_by_email("alice@example.com")
        assert deleted == 2
        rows = db.load_all_active()
        assert len(rows) == 0
        db.close()


class TestAPIKeyStoreWriteThrough:
    """Tests for APIKeyStore with DB persistence."""

    def test_create_writes_through(self, tmp_path):
        db = UserDB(str(tmp_path / "test.db"))
        store = APIKeyStore()
        store.load_from_db(db)
        raw_key = store.create("alice@example.com", Tier.FREE)
        rows = db.load_all_active()
        assert len(rows) == 1
        assert rows[0]["email"] == "alice@example.com"
        assert rows[0]["key_hash"] == _hash(raw_key)
        db.close()

    def test_load_from_db_makes_keys_validatable(self, tmp_path):
        db = UserDB(str(tmp_path / "test.db"))
        h = _hash("mb_free_testkey123")
        db.insert_key(h, "alice@example.com", "free", "chemist")
        store = APIKeyStore()
        store.load_from_db(db)
        record = store.validate("mb_free_testkey123")
        assert record is not None
        assert record.email == "alice@example.com"
        assert record.tier == Tier.FREE
        assert record.role == Role.CHEMIST
        db.close()

    def test_delete_writes_through(self, tmp_path):
        db = UserDB(str(tmp_path / "test.db"))
        store = APIKeyStore()
        store.load_from_db(db)
        store.create("alice@example.com", Tier.PRO)
        store.delete_by_email("alice@example.com")
        rows = db.load_all_active()
        assert len(rows) == 0
        db.close()

    def test_no_db_pure_in_memory(self):
        """Without a DB, APIKeyStore works exactly as before."""
        store = APIKeyStore()
        raw_key = store.create("alice@example.com", Tier.FREE)
        record = store.validate(raw_key)
        assert record is not None
        assert record.email == "alice@example.com"
        store.delete_by_email("alice@example.com")
        assert store.validate(raw_key) is None
