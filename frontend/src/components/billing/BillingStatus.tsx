import type { BillingStatusResponse } from "@/api/types";
import { Card, CardHeader, CardTitle } from "@/components/ui/Card";
import { Badge } from "@/components/ui/Badge";
import { Button } from "@/components/ui/Button";

interface BillingStatusProps {
  status: BillingStatusResponse;
  onManage: () => void;
}

export function BillingStatus({ status, onManage }: BillingStatusProps) {
  return (
    <Card>
      <CardHeader>
        <CardTitle>Billing</CardTitle>
      </CardHeader>
      <div className="space-y-3 text-sm">
        <div className="flex justify-between">
          <span className="text-text-secondary">Email</span>
          <span>{status.email}</span>
        </div>
        <div className="flex justify-between">
          <span className="text-text-secondary">Tier</span>
          <Badge variant={status.tier === "pro" ? "accent" : "default"}>
            {status.tier.charAt(0).toUpperCase() + status.tier.slice(1)}
          </Badge>
        </div>
        <div className="flex justify-between">
          <span className="text-text-secondary">Subscription</span>
          <Badge variant={status.subscription_status === "active" ? "success" : "default"}>
            {status.subscription_status === "none" ? "None" : status.subscription_status}
          </Badge>
        </div>
      </div>
      {status.stripe_customer_id && (
        <Button variant="secondary" size="sm" className="mt-4" onClick={onManage}>
          Manage Subscription
        </Button>
      )}
    </Card>
  );
}
