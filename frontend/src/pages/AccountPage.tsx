import { useEffect } from "react";
import { useAuthStore } from "@/stores/auth-store";
import { useBilling } from "@/hooks/useBilling";
import { BillingStatus } from "@/components/billing/BillingStatus";
import { PricingCards } from "@/components/billing/PricingCards";
import { Card, CardHeader, CardTitle } from "@/components/ui/Card";
import { Badge } from "@/components/ui/Badge";
import { Alert } from "@/components/ui/Alert";
import { Skeleton } from "@/components/ui/Skeleton";
import { useSearchParams } from "react-router-dom";

export default function AccountPage() {
  const { email, tier, apiKey } = useAuthStore();
  const { fetchStatus, checkout, openPortal, status, loading, error } = useBilling();
  const [searchParams] = useSearchParams();
  const checkoutResult = searchParams.get("checkout");

  const refreshToken = useAuthStore((s) => s.refreshToken);

  useEffect(() => {
    if (checkoutResult === "success") {
      refreshToken().then(() => fetchStatus());
    } else {
      fetchStatus();
    }
  }, [fetchStatus, checkoutResult, refreshToken]);

  return (
    <div className="space-y-6">
      <div>
        <h1 className="text-2xl font-bold mb-1">Account</h1>
        <p className="text-sm text-text-secondary">
          Manage your API key and subscription.
        </p>
      </div>

      {checkoutResult === "success" && (
        <Alert variant="success">Subscription activated! Your tier has been upgraded.</Alert>
      )}
      {checkoutResult === "cancel" && (
        <Alert variant="warning">Checkout was cancelled. No changes were made.</Alert>
      )}

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* API Key info */}
        <Card>
          <CardHeader>
            <CardTitle>API Key</CardTitle>
          </CardHeader>
          <div className="space-y-3 text-sm">
            <div className="flex justify-between">
              <span className="text-text-secondary">Email</span>
              <span>{email ?? "--"}</span>
            </div>
            <div className="flex justify-between">
              <span className="text-text-secondary">Tier</span>
              <Badge variant={tier === "pro" ? "accent" : "default"}>
                {tier ? tier.charAt(0).toUpperCase() + tier.slice(1) : "Free"}
              </Badge>
            </div>
            <div className="flex justify-between items-start">
              <span className="text-text-secondary">Key</span>
              <code className="text-xs font-mono text-text-muted">
                {apiKey ? `${apiKey.slice(0, 12)}...${apiKey.slice(-4)}` : "--"}
              </code>
            </div>
          </div>
        </Card>

        {/* Billing status */}
        {loading ? (
          <Skeleton className="h-48" />
        ) : error ? (
          <Alert variant="error">{error}</Alert>
        ) : status ? (
          <BillingStatus status={status} onManage={openPortal} />
        ) : null}
      </div>

      <div>
        <h2 className="text-xl font-bold mb-4">Plans</h2>
        <PricingCards
          onCheckout={checkout}
          currentTier={status?.tier ?? tier ?? "free"}
        />
      </div>
    </div>
  );
}
