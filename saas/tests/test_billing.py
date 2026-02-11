"""Tests for billing endpoints with mocked Stripe."""

import json
from unittest.mock import MagicMock, patch

import pytest


class TestBillingStatus:
    def test_status_unauthenticated(self, client):
        resp = client.get("/api/v1/billing/status")
        assert resp.status_code == 401

    def test_status_free_user(self, client, auth_headers):
        resp = client.get("/api/v1/billing/status", headers=auth_headers)
        assert resp.status_code == 200
        data = resp.json()
        assert data["tier"] == "free"
        assert data["subscription_status"] == "none"

    def test_status_returns_email(self, client, auth_headers):
        resp = client.get("/api/v1/billing/status", headers=auth_headers)
        assert resp.status_code == 200
        assert "email" in resp.json()


class TestCheckout:
    def test_checkout_requires_auth(self, client):
        resp = client.post(
            "/api/v1/billing/checkout",
            json={"plan": "pro_monthly"},
        )
        assert resp.status_code == 401

    def test_checkout_no_stripe_returns_501(self, client, auth_headers):
        """Without Stripe configured, checkout returns 501."""
        resp = client.post(
            "/api/v1/billing/checkout",
            json={"plan": "pro_monthly"},
            headers=auth_headers,
        )
        assert resp.status_code == 501

    @patch("app.routers.billing.settings")
    @patch("app.services.stripe_service.stripe")
    def test_checkout_with_stripe(self, mock_stripe, mock_settings, client, auth_headers):
        """With Stripe configured, checkout returns a URL."""
        mock_settings.stripe_secret_key = "sk_test_fake"
        mock_settings.stripe_webhook_secret = "whsec_fake"
        mock_settings.stripe_pro_monthly_price_id = "price_monthly"
        mock_settings.stripe_pro_yearly_price_id = "price_yearly"

        # Mock Stripe Customer.list and Customer.create
        mock_stripe.Customer.list.return_value = MagicMock(data=[])
        mock_customer = MagicMock()
        mock_customer.id = "cus_test123"
        mock_stripe.Customer.create.return_value = mock_customer

        # Mock Checkout Session
        mock_session = MagicMock()
        mock_session.url = "https://checkout.stripe.com/test_session"
        mock_stripe.checkout.Session.create.return_value = mock_session

        resp = client.post(
            "/api/v1/billing/checkout",
            json={"plan": "pro_monthly"},
            headers=auth_headers,
        )
        assert resp.status_code == 200
        assert resp.json()["checkout_url"] == "https://checkout.stripe.com/test_session"

    def test_checkout_invalid_plan(self, client, auth_headers):
        resp = client.post(
            "/api/v1/billing/checkout",
            json={"plan": "invalid_plan"},
            headers=auth_headers,
        )
        assert resp.status_code == 422


class TestPortal:
    def test_portal_requires_auth(self, client):
        resp = client.post("/api/v1/billing/portal")
        assert resp.status_code == 401

    def test_portal_no_stripe_returns_501(self, client, auth_headers):
        resp = client.post("/api/v1/billing/portal", headers=auth_headers)
        assert resp.status_code == 501


class TestWebhook:
    def test_webhook_no_stripe_returns_501(self, client):
        resp = client.post(
            "/api/v1/billing/webhook",
            content=b"{}",
            headers={"stripe-signature": "fake"},
        )
        assert resp.status_code == 501

    def test_webhook_bad_signature(self, client):
        """With Stripe configured but bad signature, returns 400."""
        with patch("app.routers.billing.settings") as mock_settings:
            mock_settings.stripe_secret_key = "sk_test_fake"
            mock_settings.stripe_webhook_secret = "whsec_fake"
            resp = client.post(
                "/api/v1/billing/webhook",
                content=b'{"type": "checkout.session.completed"}',
                headers={"stripe-signature": "bad_sig"},
            )
            assert resp.status_code == 400


class TestUpdateTier:
    def test_update_tier_in_memory(self):
        """Verify the api_key_store.update_tier method works."""
        from app.auth.api_keys import api_key_store
        from app.config import Tier

        # Create a free key
        key = api_key_store.create(email="upgrade@test.com", tier=Tier.FREE)
        record = api_key_store.validate(key)
        assert record.tier == Tier.FREE

        # Upgrade to pro
        count = api_key_store.update_tier("upgrade@test.com", Tier.PRO)
        assert count == 1

        # Verify in-memory tier changed
        record = api_key_store.validate(key)
        assert record.tier == Tier.PRO

    def test_update_tier_nonexistent_email(self):
        from app.auth.api_keys import api_key_store
        from app.config import Tier

        count = api_key_store.update_tier("nobody@test.com", Tier.PRO)
        assert count == 0
