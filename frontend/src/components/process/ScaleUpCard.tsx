import type { ScaleUpResponse } from "@/api/types";
import { Card, CardHeader, CardTitle } from "@/components/ui/Card";
import { Badge } from "@/components/ui/Badge";
import { usd, num } from "@/lib/format";

interface ScaleUpCardProps {
  scaleUp: ScaleUpResponse;
}

export function ScaleUpCard({ scaleUp }: ScaleUpCardProps) {
  return (
    <Card>
      <CardHeader>
        <div className="flex items-center justify-between">
          <CardTitle>Scale-Up Analysis</CardTitle>
          <Badge variant="accent">{scaleUp.recommended_mode}</Badge>
        </div>
      </CardHeader>
      <div className="grid grid-cols-2 gap-4 text-sm mb-4">
        <div>
          <span className="text-text-muted block">Target Annual</span>
          <span className="font-bold">{num(scaleUp.target_annual_kg)} kg/yr</span>
        </div>
        <div>
          <span className="text-text-muted block">Annual Capacity</span>
          <span className="font-bold">{num(scaleUp.annual_capacity_kg)} kg/yr</span>
        </div>
        {scaleUp.batch_size_kg != null && (
          <div>
            <span className="text-text-muted block">Batch Size</span>
            <span className="font-bold">{num(scaleUp.batch_size_kg)} kg</span>
          </div>
        )}
        {scaleUp.batches_per_year != null && (
          <div>
            <span className="text-text-muted block">Batches/Year</span>
            <span className="font-bold">{scaleUp.batches_per_year}</span>
          </div>
        )}
        <div>
          <span className="text-text-muted block">Cycle Time</span>
          <span className="font-bold">{num(scaleUp.cycle_time_hours)} h</span>
        </div>
        <div>
          <span className="text-text-muted block">Capital Cost</span>
          <span className="font-bold">{usd(scaleUp.capital_cost_usd)}</span>
        </div>
        <div>
          <span className="text-text-muted block">Operating Cost</span>
          <span className="font-bold">{usd(scaleUp.operating_cost_annual_usd)}/yr</span>
        </div>
      </div>

      {scaleUp.scale_up_risks.length > 0 && (
        <div className="mb-3">
          <h4 className="text-xs font-medium text-text-muted mb-1">Risks</h4>
          <ul className="text-sm text-text-secondary space-y-0.5">
            {scaleUp.scale_up_risks.map((r) => (
              <li key={r} className="flex items-start gap-1.5">
                <span className="text-yellow mt-0.5">!</span>
                {r}
              </li>
            ))}
          </ul>
        </div>
      )}

      {scaleUp.recommendations.length > 0 && (
        <div>
          <h4 className="text-xs font-medium text-text-muted mb-1">Recommendations</h4>
          <ul className="text-sm text-text-secondary space-y-0.5">
            {scaleUp.recommendations.map((r) => (
              <li key={r} className="flex items-start gap-1.5">
                <span className="text-green mt-0.5">-</span>
                {r}
              </li>
            ))}
          </ul>
        </div>
      )}
    </Card>
  );
}
