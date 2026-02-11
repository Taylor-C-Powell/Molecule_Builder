"""Tests for billing endpoints: create_checkout and billing_status."""

from __future__ import annotations

import httpx
import respx

from molbuilder_client import BillingStatus, CheckoutSession, MolBuilder


def test_create_checkout(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.post("/billing/checkout").mock(
        return_value=httpx.Response(
            200,
            json={"checkout_url": "https://checkout.stripe.com/session_abc"},
        )
    )

    result = client.create_checkout("pro_monthly")

    assert isinstance(result, CheckoutSession)
    assert result.checkout_url == "https://checkout.stripe.com/session_abc"


def test_create_checkout_sends_plan(client: MolBuilder, mock_api: respx.Router) -> None:
    route = mock_api.post("/billing/checkout").mock(
        return_value=httpx.Response(200, json={"checkout_url": "https://example.com"}),
    )

    client.create_checkout("pro_yearly")

    import json
    body = json.loads(route.calls.last.request.content)
    assert body["plan"] == "pro_yearly"


def test_billing_status(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/billing/status").mock(
        return_value=httpx.Response(
            200,
            json={
                "email": "user@example.com",
                "tier": "pro",
                "subscription_status": "active",
                "stripe_customer_id": "cus_abc",
                "stripe_subscription_id": "sub_xyz",
            },
        )
    )

    result = client.billing_status()

    assert isinstance(result, BillingStatus)
    assert result.email == "user@example.com"
    assert result.tier == "pro"
    assert result.subscription_status == "active"
    assert result.stripe_customer_id == "cus_abc"


def test_billing_status_no_subscription(client: MolBuilder, mock_api: respx.Router) -> None:
    mock_api.get("/billing/status").mock(
        return_value=httpx.Response(
            200,
            json={
                "email": "free@example.com",
                "tier": "free",
                "subscription_status": "none",
                "stripe_customer_id": None,
                "stripe_subscription_id": None,
            },
        )
    )

    result = client.billing_status()
    assert result.subscription_status == "none"
    assert result.stripe_customer_id is None
    assert result.stripe_subscription_id is None
