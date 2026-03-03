"""Tests for team management endpoints."""

import uuid

import pytest
from fastapi.testclient import TestClient

from app.auth.api_keys import api_key_store
from app.config import Tier
from app.main import app


@pytest.fixture()
def client():
    with TestClient(app) as c:
        yield c


def _make_team_key(email=None):
    email = email or f"team-{uuid.uuid4().hex[:8]}@example.com"
    key = api_key_store.create(email=email, tier=Tier.TEAM)
    return key, email


def _make_free_key(email=None):
    email = email or f"free-{uuid.uuid4().hex[:8]}@example.com"
    key = api_key_store.create(email=email, tier=Tier.FREE)
    return key, email


def _create_team(client, headers, name="Test Team", slug=None):
    slug = slug or f"test-{uuid.uuid4().hex[:6]}"
    resp = client.post(
        "/api/v1/teams/",
        json={"name": name, "slug": slug},
        headers=headers,
    )
    return resp


# ------------------------------------------------------------------ #
# Team CRUD
# ------------------------------------------------------------------ #

class TestTeamCRUD:
    def test_create_team(self, client):
        key, email = _make_team_key()
        resp = _create_team(client, {"X-API-Key": key})
        assert resp.status_code == 201
        data = resp.json()
        assert data["name"] == "Test Team"
        assert data["owner_email"] == email
        assert "id" in data

    def test_create_team_duplicate_slug(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        slug = f"dup-{uuid.uuid4().hex[:6]}"
        resp1 = _create_team(client, headers, slug=slug)
        assert resp1.status_code == 201
        resp2 = _create_team(client, headers, slug=slug)
        assert resp2.status_code == 409

    def test_create_team_free_tier_forbidden(self, client):
        key, _ = _make_free_key()
        resp = _create_team(client, {"X-API-Key": key})
        assert resp.status_code == 403

    def test_list_teams_empty(self, client):
        key, _ = _make_team_key()
        resp = client.get("/api/v1/teams/", headers={"X-API-Key": key})
        assert resp.status_code == 200
        assert resp.json()["teams"] == []

    def test_list_teams_after_create(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        _create_team(client, headers, name="Alpha")
        _create_team(client, headers, name="Beta", slug=f"beta-{uuid.uuid4().hex[:6]}")
        resp = client.get("/api/v1/teams/", headers=headers)
        assert resp.status_code == 200
        assert len(resp.json()["teams"]) == 2

    def test_get_team(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]
        resp = client.get(f"/api/v1/teams/{team_id}", headers=headers)
        assert resp.status_code == 200
        assert resp.json()["id"] == team_id

    def test_get_team_nonmember(self, client):
        key1, _ = _make_team_key()
        key2, _ = _make_team_key()
        create_resp = _create_team(client, {"X-API-Key": key1})
        team_id = create_resp.json()["id"]
        resp = client.get(f"/api/v1/teams/{team_id}", headers={"X-API-Key": key2})
        assert resp.status_code == 403

    def test_update_team(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]
        resp = client.patch(
            f"/api/v1/teams/{team_id}",
            json={"name": "Updated Name"},
            headers=headers,
        )
        assert resp.status_code == 200
        assert resp.json()["name"] == "Updated Name"

    def test_delete_team(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]
        resp = client.delete(f"/api/v1/teams/{team_id}", headers=headers)
        assert resp.status_code == 200
        assert resp.json()["status"] == "deleted"
        # Verify gone
        resp2 = client.get(f"/api/v1/teams/{team_id}", headers=headers)
        assert resp2.status_code == 404

    def test_get_team_not_found(self, client):
        key, _ = _make_team_key()
        resp = client.get("/api/v1/teams/99999", headers={"X-API-Key": key})
        assert resp.status_code == 404


# ------------------------------------------------------------------ #
# Member management
# ------------------------------------------------------------------ #

class TestMembers:
    def test_owner_auto_added(self, client):
        key, email = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]
        resp = client.get(f"/api/v1/teams/{team_id}/members", headers=headers)
        assert resp.status_code == 200
        members = resp.json()["members"]
        assert len(members) == 1
        assert members[0]["user_email"] == email
        assert members[0]["team_role"] == "owner"

    def test_add_member(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        new_email = f"member-{uuid.uuid4().hex[:8]}@example.com"
        resp = client.post(
            f"/api/v1/teams/{team_id}/members",
            json={"email": new_email, "role": "member"},
            headers=headers,
        )
        assert resp.status_code == 201
        assert resp.json()["user_email"] == new_email
        assert resp.json()["team_role"] == "member"

    def test_add_member_duplicate(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        new_email = f"member-{uuid.uuid4().hex[:8]}@example.com"
        client.post(
            f"/api/v1/teams/{team_id}/members",
            json={"email": new_email},
            headers=headers,
        )
        resp = client.post(
            f"/api/v1/teams/{team_id}/members",
            json={"email": new_email},
            headers=headers,
        )
        assert resp.status_code == 409

    def test_add_member_as_owner_rejected(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        resp = client.post(
            f"/api/v1/teams/{team_id}/members",
            json={"email": "x@example.com", "role": "owner"},
            headers=headers,
        )
        assert resp.status_code == 400

    def test_member_cannot_add_members(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        # Add a member
        member_email = f"member-{uuid.uuid4().hex[:8]}@example.com"
        member_key = api_key_store.create(email=member_email, tier=Tier.TEAM)
        client.post(
            f"/api/v1/teams/{team_id}/members",
            json={"email": member_email},
            headers=headers,
        )
        # Member tries to add another
        resp = client.post(
            f"/api/v1/teams/{team_id}/members",
            json={"email": "new@example.com"},
            headers={"X-API-Key": member_key},
        )
        assert resp.status_code == 403

    def test_update_member_role(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        member_email = f"member-{uuid.uuid4().hex[:8]}@example.com"
        client.post(
            f"/api/v1/teams/{team_id}/members",
            json={"email": member_email},
            headers=headers,
        )
        resp = client.patch(
            f"/api/v1/teams/{team_id}/members/{member_email}",
            json={"role": "admin"},
            headers=headers,
        )
        assert resp.status_code == 200
        assert resp.json()["team_role"] == "admin"

    def test_cannot_change_owner_role(self, client):
        key, owner_email = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        resp = client.patch(
            f"/api/v1/teams/{team_id}/members/{owner_email}",
            json={"role": "member"},
            headers=headers,
        )
        assert resp.status_code == 400

    def test_remove_member(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        member_email = f"member-{uuid.uuid4().hex[:8]}@example.com"
        client.post(
            f"/api/v1/teams/{team_id}/members",
            json={"email": member_email},
            headers=headers,
        )
        resp = client.delete(
            f"/api/v1/teams/{team_id}/members/{member_email}",
            headers=headers,
        )
        assert resp.status_code == 200

    def test_cannot_remove_owner(self, client):
        key, owner_email = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        resp = client.delete(
            f"/api/v1/teams/{team_id}/members/{owner_email}",
            headers=headers,
        )
        assert resp.status_code == 400

    def test_self_removal(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        member_email = f"member-{uuid.uuid4().hex[:8]}@example.com"
        member_key = api_key_store.create(email=member_email, tier=Tier.TEAM)
        client.post(
            f"/api/v1/teams/{team_id}/members",
            json={"email": member_email},
            headers=headers,
        )
        # Member removes themselves
        resp = client.delete(
            f"/api/v1/teams/{team_id}/members/{member_email}",
            headers={"X-API-Key": member_key},
        )
        assert resp.status_code == 200

    def test_nonmember_cannot_delete_other(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        member_email = f"member-{uuid.uuid4().hex[:8]}@example.com"
        client.post(
            f"/api/v1/teams/{team_id}/members",
            json={"email": member_email},
            headers=headers,
        )

        # Different non-member tries to remove
        outsider_key, _ = _make_team_key()
        resp = client.delete(
            f"/api/v1/teams/{team_id}/members/{member_email}",
            headers={"X-API-Key": outsider_key},
        )
        assert resp.status_code == 403

    def test_delete_team_cascades_members(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        member_email = f"member-{uuid.uuid4().hex[:8]}@example.com"
        client.post(
            f"/api/v1/teams/{team_id}/members",
            json={"email": member_email},
            headers=headers,
        )
        resp = client.delete(f"/api/v1/teams/{team_id}", headers=headers)
        assert resp.status_code == 200


# ------------------------------------------------------------------ #
# Team library
# ------------------------------------------------------------------ #

class TestTeamLibrary:
    def test_save_molecule(self, client):
        key, email = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        resp = client.post(
            f"/api/v1/teams/{team_id}/library/",
            json={"smiles": "CCO", "name": "Ethanol"},
            headers=headers,
        )
        assert resp.status_code == 200
        data = resp.json()
        assert "C" in data["smiles"]  # canonical form may differ
        assert data["added_by"] == email
        assert data["team_id"] == team_id

    def test_save_duplicate_molecule(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        client.post(
            f"/api/v1/teams/{team_id}/library/",
            json={"smiles": "CCO"},
            headers=headers,
        )
        resp = client.post(
            f"/api/v1/teams/{team_id}/library/",
            json={"smiles": "CCO"},
            headers=headers,
        )
        assert resp.status_code == 409

    def test_list_molecules(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        client.post(
            f"/api/v1/teams/{team_id}/library/",
            json={"smiles": "CCO"},
            headers=headers,
        )
        resp = client.get(
            f"/api/v1/teams/{team_id}/library/",
            headers=headers,
        )
        assert resp.status_code == 200
        data = resp.json()
        assert data["total"] == 1
        assert len(data["molecules"]) == 1

    def test_get_molecule(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        save_resp = client.post(
            f"/api/v1/teams/{team_id}/library/",
            json={"smiles": "CCO", "name": "Ethanol"},
            headers=headers,
        )
        mol_id = save_resp.json()["id"]

        resp = client.get(
            f"/api/v1/teams/{team_id}/library/{mol_id}",
            headers=headers,
        )
        assert resp.status_code == 200
        assert resp.json()["name"] == "Ethanol"

    def test_update_molecule(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        save_resp = client.post(
            f"/api/v1/teams/{team_id}/library/",
            json={"smiles": "CCO"},
            headers=headers,
        )
        mol_id = save_resp.json()["id"]

        resp = client.patch(
            f"/api/v1/teams/{team_id}/library/{mol_id}",
            json={"name": "Updated", "tags": ["alcohol"]},
            headers=headers,
        )
        assert resp.status_code == 200
        assert resp.json()["name"] == "Updated"
        assert resp.json()["tags"] == ["alcohol"]

    def test_delete_molecule(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        save_resp = client.post(
            f"/api/v1/teams/{team_id}/library/",
            json={"smiles": "CCO"},
            headers=headers,
        )
        mol_id = save_resp.json()["id"]

        resp = client.delete(
            f"/api/v1/teams/{team_id}/library/{mol_id}",
            headers=headers,
        )
        assert resp.status_code == 200

        resp2 = client.get(
            f"/api/v1/teams/{team_id}/library/{mol_id}",
            headers=headers,
        )
        assert resp2.status_code == 404

    def test_import_molecules(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        resp = client.post(
            f"/api/v1/teams/{team_id}/library/import",
            json={"smiles_list": ["CCO", "CC", "INVALID!!!"], "tag": "batch"},
            headers=headers,
        )
        assert resp.status_code == 200
        data = resp.json()
        assert data["saved"] == 2
        assert len(data["errors"]) >= 1

    def test_nonmember_cannot_access_library(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        outsider_key, _ = _make_team_key()
        resp = client.get(
            f"/api/v1/teams/{team_id}/library/",
            headers={"X-API-Key": outsider_key},
        )
        assert resp.status_code == 403

    def test_member_can_access_library(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        # Add member
        member_email = f"member-{uuid.uuid4().hex[:8]}@example.com"
        member_key = api_key_store.create(email=member_email, tier=Tier.TEAM)
        client.post(
            f"/api/v1/teams/{team_id}/members",
            json={"email": member_email},
            headers=headers,
        )

        # Member saves a molecule
        resp = client.post(
            f"/api/v1/teams/{team_id}/library/",
            json={"smiles": "CCO"},
            headers={"X-API-Key": member_key},
        )
        assert resp.status_code == 200
        assert resp.json()["added_by"] == member_email

    def test_invalid_smiles(self, client):
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        resp = client.post(
            f"/api/v1/teams/{team_id}/library/",
            json={"smiles": "INVALID!!!"},
            headers=headers,
        )
        assert resp.status_code == 400

    def test_member_delete_not_owner(self, client):
        """Non-owner member cannot delete the team."""
        key, _ = _make_team_key()
        headers = {"X-API-Key": key}
        create_resp = _create_team(client, headers)
        team_id = create_resp.json()["id"]

        member_email = f"member-{uuid.uuid4().hex[:8]}@example.com"
        member_key = api_key_store.create(email=member_email, tier=Tier.TEAM)
        client.post(
            f"/api/v1/teams/{team_id}/members",
            json={"email": member_email},
            headers=headers,
        )

        resp = client.delete(
            f"/api/v1/teams/{team_id}",
            headers={"X-API-Key": member_key},
        )
        assert resp.status_code == 403
