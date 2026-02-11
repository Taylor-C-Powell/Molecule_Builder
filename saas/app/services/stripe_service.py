"""Stripe billing operations."""

import logging

import stripe

from app.auth.api_keys import api_key_store
from app.config import Tier, settings
from app.services.user_db import get_user_db

logger = logging.getLogger("molbuilder.stripe")


def _configure() -> None:
    """Set the Stripe API key from settings."""
    stripe.api_key = settings.stripe_secret_key


def get_or_create_customer(email: str) -> str:
    """Find existing Stripe customer by email, or create one. Returns customer ID."""
    _configure()
    existing = stripe.Customer.list(email=email, limit=1)
    if existing.data:
        return existing.data[0].id
    customer = stripe.Customer.create(email=email)
    db = get_user_db()
    db.update_stripe_customer(email, customer.id)
    return customer.id


def create_checkout_session(
    customer_id: str, price_id: str, success_url: str, cancel_url: str
) -> str:
    """Create a Stripe Checkout session. Returns the checkout URL."""
    _configure()
    session = stripe.checkout.Session.create(
        customer=customer_id,
        payment_method_types=["card"],
        line_items=[{"price": price_id, "quantity": 1}],
        mode="subscription",
        success_url=success_url,
        cancel_url=cancel_url,
    )
    return session.url


def create_portal_session(customer_id: str, return_url: str) -> str:
    """Create a Stripe Customer Portal session. Returns the portal URL."""
    _configure()
    session = stripe.billing_portal.Session.create(
        customer=customer_id,
        return_url=return_url,
    )
    return session.url


def upgrade_user(email: str, tier: Tier, subscription_id: str) -> None:
    """Upgrade a user's tier after successful payment."""
    db = get_user_db()
    db.update_subscription(email, subscription_id, "active", tier.value)
    api_key_store.update_tier(email, tier)
    logger.info("Upgraded %s to %s (subscription %s)", email, tier.value, subscription_id)


def downgrade_user(email: str) -> None:
    """Revert a user to free tier after subscription ends."""
    db = get_user_db()
    db.update_subscription(email, None, "canceled", Tier.FREE.value)
    api_key_store.update_tier(email, Tier.FREE)
    logger.info("Downgraded %s to free tier", email)


def handle_webhook_event(payload: bytes, sig_header: str) -> dict:
    """Verify and dispatch a Stripe webhook event. Returns event summary."""
    event = stripe.Webhook.construct_event(
        payload, sig_header, settings.stripe_webhook_secret
    )

    if event.type == "checkout.session.completed":
        session = event.data.object
        customer_id = session.customer
        subscription_id = session.subscription
        # Look up user by Stripe customer
        db = get_user_db()
        user = db.get_by_stripe_customer(customer_id)
        if user:
            upgrade_user(user["email"], Tier.PRO, subscription_id)
            return {"action": "upgraded", "email": user["email"], "tier": "pro"}
        else:
            logger.warning("Webhook: no user found for customer %s", customer_id)
            return {"action": "skipped", "reason": "unknown_customer"}

    elif event.type == "customer.subscription.deleted":
        subscription = event.data.object
        customer_id = subscription.customer
        db = get_user_db()
        user = db.get_by_stripe_customer(customer_id)
        if user:
            downgrade_user(user["email"])
            return {"action": "downgraded", "email": user["email"]}
        return {"action": "skipped", "reason": "unknown_customer"}

    elif event.type == "invoice.payment_failed":
        invoice = event.data.object
        customer_id = invoice.customer
        db = get_user_db()
        user = db.get_by_stripe_customer(customer_id)
        if user:
            db.update_subscription(
                user["email"],
                user.get("stripe_subscription_id"),
                "past_due",
                user["tier"],
            )
            logger.warning("Payment failed for %s", user["email"])
            return {"action": "payment_failed", "email": user["email"]}
        return {"action": "skipped", "reason": "unknown_customer"}

    elif event.type == "customer.subscription.updated":
        subscription = event.data.object
        customer_id = subscription.customer
        status = subscription.status
        db = get_user_db()
        user = db.get_by_stripe_customer(customer_id)
        if user and status == "active":
            db.update_subscription(
                user["email"], subscription.id, "active", Tier.PRO.value
            )
            api_key_store.update_tier(user["email"], Tier.PRO)
            return {"action": "subscription_updated", "status": status}
        return {"action": "subscription_updated", "status": status}

    return {"action": "ignored", "event_type": event.type}
