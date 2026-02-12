import type { CostEstimateResponse } from "@/api/types";
import { Card, CardHeader, CardTitle } from "@/components/ui/Card";
import { usd, num } from "@/lib/format";

interface CostBreakdownProps {
  cost: CostEstimateResponse;
}

export function CostBreakdown({ cost }: CostBreakdownProps) {
  const b = cost.breakdown;
  const items: [string, number][] = [
    ["Raw Materials", b.raw_materials_usd],
    ["Labor", b.labor_usd],
    ["Equipment", b.equipment_usd],
    ["Energy", b.energy_usd],
    ["Waste Disposal", b.waste_disposal_usd],
    ["Overhead", b.overhead_usd],
  ];

  return (
    <Card>
      <CardHeader>
        <div className="flex items-center justify-between">
          <CardTitle>Cost Estimate</CardTitle>
          <span className="text-lg font-bold text-accent">{usd(cost.total_usd)}</span>
        </div>
      </CardHeader>
      <div className="space-y-2 text-sm mb-4">
        <div className="flex justify-between font-medium text-text-secondary border-b border-border pb-1">
          <span>Per kg</span>
          <span>{usd(cost.per_kg_usd)}/kg</span>
        </div>
        <div className="flex justify-between text-text-muted border-b border-border/50 pb-1">
          <span>Scale</span>
          <span>{num(cost.scale_kg)} kg</span>
        </div>
      </div>
      <div className="space-y-1.5">
        {items.map(([label, val]) => {
          const pct = cost.total_usd > 0 ? (val / cost.total_usd) * 100 : 0;
          return (
            <div key={label} className="text-sm">
              <div className="flex justify-between mb-0.5">
                <span className="text-text-secondary">{label}</span>
                <span className="font-mono">{usd(val)}</span>
              </div>
              <div className="h-1.5 bg-border rounded-full overflow-hidden">
                <div
                  className="h-full bg-accent rounded-full"
                  style={{ width: `${pct}%` }}
                />
              </div>
            </div>
          );
        })}
      </div>
      {cost.notes.length > 0 && (
        <div className="mt-4 pt-3 border-t border-border/50">
          {cost.notes.map((n, i) => (
            <p key={i} className="text-xs text-text-muted">{n}</p>
          ))}
        </div>
      )}
    </Card>
  );
}
