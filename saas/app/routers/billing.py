"""Billing endpoints: Stripe checkout, portal, webhooks, status."""

from urllib.parse import urlparse

from fastapi import APIRouter, Depends, Request
from fastapi.responses import JSONResponse

from app.config import settings
from app.dependencies import UserContext, get_current_user, require_admin
from app.exceptions import MolBuilderAPIError
from app.models.billing import (
    BillingStatus, CheckoutRequest, CheckoutResponse, PortalRequest, PortalResponse,
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


def _validate_redirect_url(url: str) -> str:
    """Validate that a redirect URL points to an allowed host."""
    allowed = {h.strip().lower() for h in settings.allowed_redirect_hosts.split(",")}
    try:
        parsed = urlparse(url)
        if parsed.scheme not in ("https", "http"):
            raise MolBuilderAPIError(422, f"Redirect URL must use https: {url}")
        if parsed.hostname and parsed.hostname.lower() not in allowed:
            raise MolBuilderAPIError(
                422,
                f"Redirect URL host '{parsed.hostname}' not in allowed list. "
                f"Allowed: {', '.join(sorted(allowed))}",
            )
    except ValueError:
        raise MolBuilderAPIError(422, f"Invalid redirect URL: {url}")
    return url


@router.post("/checkout", response_model=CheckoutResponse)
def create_checkout(body: CheckoutRequest, user: UserContext = Depends(get_current_user)):
    """Create a Stripe Checkout session for upgrading to Pro."""
    _require_stripe()

    # Validate redirect URLs to prevent open redirect
    success_url = _validate_redirect_url(body.success_url)
    cancel_url = _validate_redirect_url(body.cancel_url)

    plan_to_price = {
        "pro_monthly": settings.stripe_pro_monthly_price_id,
        "pro_yearly": settings.stripe_pro_yearly_price_id,
        "team_monthly": settings.stripe_team_monthly_price_id,
        "team_yearly": settings.stripe_team_yearly_price_id,
    }
    price_id = plan_to_price.get(body.plan, "")
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
        success_url=success_url,
        cancel_url=cancel_url,
    )
    return CheckoutResponse(checkout_url=checkout_url)


@router.post("/portal", response_model=PortalResponse)
def create_portal(
    body: PortalRequest = PortalRequest(),
    user: UserContext = Depends(get_current_user),
):
    """Create a Stripe Customer Portal session for managing subscription."""
    _require_stripe()

    return_url = _validate_redirect_url(body.return_url)

    db = get_user_db()
    info = db.get_stripe_info(user.email)
    if not info or not info.get("stripe_customer_id"):
        from fastapi import HTTPException
        raise HTTPException(status_code=400, detail="No billing account found")

    portal_url = stripe_service.create_portal_session(
        customer_id=info["stripe_customer_id"],
        return_url=return_url,
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
    except ValueError as e:
        return JSONResponse(status_code=400, content={"error": str(e)})
    except Exception:
        return JSONResponse(status_code=400, content={"error": "Webhook processing failed"})


@router.get("/status", response_model=BillingStatus)
def billing_status(user: UserContext = Depends(get_current_user)):
    """Get current billing/subscription status. Stripe IDs only visible to admins."""
    db = get_user_db()
    info = db.get_stripe_info(user.email)

    is_admin = user.role.value == "admin"

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
        has_billing=bool(info.get("stripe_customer_id")),
        # Only expose Stripe IDs to admin users
        stripe_customer_id=info.get("stripe_customer_id") if is_admin else None,
        stripe_subscription_id=info.get("stripe_subscription_id") if is_admin else None,
    )
