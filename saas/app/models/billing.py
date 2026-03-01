"""Billing request/response models."""

from typing import Literal
from pydantic import BaseModel, Field


class CheckoutRequest(BaseModel):
    plan: Literal["pro_monthly", "pro_yearly", "team_monthly", "team_yearly"] = Field(
        ..., description="Subscription plan to purchase"
    )
    success_url: str = Field(
        "https://www.molbuilder.io/docs",
        description="URL to redirect after successful payment",
    )
    cancel_url: str = Field(
        "https://www.molbuilder.io/docs",
        description="URL to redirect if user cancels",
    )


class CheckoutResponse(BaseModel):
    checkout_url: str = Field(..., description="Stripe Checkout URL to redirect user")


class PortalRequest(BaseModel):
    return_url: str = Field(
        "https://www.molbuilder.io/account",
        description="URL to return to after managing subscription",
    )


class PortalResponse(BaseModel):
    portal_url: str = Field(..., description="Stripe Customer Portal URL")


class BillingStatus(BaseModel):
    email: str
    tier: str
    subscription_status: str = Field(
        ..., description="none, active, past_due, or canceled"
    )
    has_billing: bool = Field(
        False, description="Whether user has a Stripe billing account"
    )
    stripe_customer_id: str | None = None
    stripe_subscription_id: str | None = None
