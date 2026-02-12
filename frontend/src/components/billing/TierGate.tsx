import { useAuthStore } from "@/stores/auth-store";
import { Card } from "@/components/ui/Card";
import { Button } from "@/components/ui/Button";
import { Link } from "react-router-dom";

interface TierGateProps {
  requiredTier: string;
  children: React.ReactNode;
}

const TIER_ORDER = ["free", "academic", "pro", "team", "enterprise"];

export function TierGate({ requiredTier, children }: TierGateProps) {
  const tier = useAuthStore((s) => s.tier) ?? "free";
  const currentIdx = TIER_ORDER.indexOf(tier);
  const requiredIdx = TIER_ORDER.indexOf(requiredTier);

  if (currentIdx >= requiredIdx) {
    return <>{children}</>;
  }

  return (
    <Card className="flex flex-col items-center justify-center py-12 text-center">
      <h3 className="text-lg font-semibold mb-2">Upgrade Required</h3>
      <p className="text-sm text-text-secondary mb-4">
        This feature requires the {requiredTier.charAt(0).toUpperCase() + requiredTier.slice(1)} plan or higher.
      </p>
      <Link to="/account">
        <Button>View Plans</Button>
      </Link>
    </Card>
  );
}
