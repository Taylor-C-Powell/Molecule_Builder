"""Billing endpoints: Stripe checkout, portal, webhooks, status."""

from fastapi import APIRouter, Depends, Request
from fastapi.responses import JSONResponse

from app.config import settings
from app.dependencies import UserContext, get_current_user
from app.models.billing import (
    BillingStatus, CheckoutRequest, CheckoutResponse, PortalResponse,
)
from app.services import stripe_service
from app.services.user_db import get_user_db

router = APIRouter(prefix="/api/v1/billing", tags=["billing"])


def _require_stripe() -> None:
    """Raise 501 if Stripe is not configured."""
    if not settings.stripe_secret_key:
        raise _stripe_not_configured()


def _stripe_not_configured():
    from fastapi import HTTPException
    return HTTPException(status_code=501, detail="Billing not configured")


@router.post("/checkout", response_model=CheckoutResponse)
def create_checkout(body: CheckoutRequest, user: UserContext = Depends(get_current_user)):
    """Create a Stripe Checkout session for upgrading to Pro."""
    _require_stripe()

    price_id = (
        settings.stripe_pro_monthly_price_id
        if body.plan == "pro_monthly"
        else settings.stripe_pro_yearly_price_id
    )
    if not price_id:
        from fastapi import HTTPException
        raise HTTPException(status_code=400, detail=f"Price not configured for {body.plan}")

    customer_id = stripe_service.get_or_create_customer(user.email)

    # Store the customer ID in our DB
    db = get_user_db()
    db.update_stripe_customer(user.email, customer_id)

    checkout_url = stripe_service.create_checkout_session(
        customer_id=customer_id,
        price_id=price_id,
        success_url=body.success_url,
        cancel_url=body.cancel_url,
    )
    return CheckoutResponse(checkout_url=checkout_url)


@router.post("/portal", response_model=PortalResponse)
def create_portal(user: UserContext = Depends(get_current_user)):
    """Create a Stripe Customer Portal session for managing subscription."""
    _require_stripe()

    db = get_user_db()
    info = db.get_stripe_info(user.email)
    if not info or not info.get("stripe_customer_id"):
        from fastapi import HTTPException
        raise HTTPException(status_code=400, detail="No billing account found")

    portal_url = stripe_service.create_portal_session(
        customer_id=info["stripe_customer_id"],
        return_url="https://molbuilder-api-production.up.railway.app/docs",
    )
    return PortalResponse(portal_url=portal_url)


@router.post("/webhook")
async def stripe_webhook(request: Request):
    """Stripe webhook handler. Verifies signature and processes events."""
    if not settings.stripe_secret_key or not settings.stripe_webhook_secret:
        return JSONResponse(status_code=501, content={"error": "Billing not configured"})

    payload = await request.body()
    sig_header = request.headers.get("stripe-signature", "")

    try:
        result = stripe_service.handle_webhook_event(payload, sig_header)
        return result
    except Exception as e:
        return JSONResponse(status_code=400, content={"error": str(e)})


@router.get("/status", response_model=BillingStatus)
def billing_status(user: UserContext = Depends(get_current_user)):
    """Get current billing/subscription status."""
    db = get_user_db()
    info = db.get_stripe_info(user.email)

    if not info:
        return BillingStatus(
            email=user.email,
            tier=user.tier.value,
            subscription_status="none",
        )

    return BillingStatus(
        email=user.email,
        tier=info["tier"],
        subscription_status=info.get("subscription_status") or "none",
        stripe_customer_id=info.get("stripe_customer_id"),
        stripe_subscription_id=info.get("stripe_subscription_id"),
    )
