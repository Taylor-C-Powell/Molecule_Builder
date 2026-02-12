import type { SynthesisRouteResponse } from "@/api/types";
import { Card, CardHeader, CardTitle } from "@/components/ui/Card";
import { Badge } from "@/components/ui/Badge";
import { pct, usd } from "@/lib/format";

interface RouteStepsProps {
  route: SynthesisRouteResponse;
}

export function RouteSteps({ route }: RouteStepsProps) {
  return (
    <div className="space-y-4">
      {/* Summary */}
      <Card>
        <CardHeader>
          <CardTitle>Best Route Summary</CardTitle>
        </CardHeader>
        <div className="grid grid-cols-2 sm:grid-cols-4 gap-4 text-sm">
          <div>
            <span className="text-text-muted block">Steps</span>
            <span className="font-bold text-lg">{route.total_steps}</span>
          </div>
          <div>
            <span className="text-text-muted block">Overall Yield</span>
            <span className="font-bold text-lg">{pct(route.overall_yield)}</span>
          </div>
          <div>
            <span className="text-text-muted block">Longest Linear</span>
            <span className="font-bold text-lg">{route.longest_linear_sequence}</span>
          </div>
          <div>
            <span className="text-text-muted block">Starting Materials</span>
            <span className="font-bold text-lg">{route.starting_materials.length}</span>
          </div>
        </div>
      </Card>

      {/* Starting materials */}
      <Card>
        <CardHeader>
          <CardTitle>Starting Materials</CardTitle>
        </CardHeader>
        <div className="space-y-2">
          {route.starting_materials.map((sm) => (
            <div key={sm.smiles} className="flex items-center justify-between text-sm py-1 border-b border-border/50 last:border-0">
              <div>
                <span className="font-medium">{sm.name}</span>
                <code className="ml-2 text-xs text-text-muted font-mono">{sm.smiles}</code>
              </div>
              <span className="text-text-secondary">{usd(sm.cost_per_kg)}/kg</span>
            </div>
          ))}
        </div>
      </Card>

      {/* Steps */}
      {route.steps.map((step) => (
        <Card key={step.step_number}>
          <CardHeader>
            <div className="flex items-center justify-between">
              <CardTitle>Step {step.step_number}: {step.reaction_name}</CardTitle>
              <Badge variant={step.expected_yield >= 0.7 ? "success" : step.expected_yield >= 0.4 ? "warning" : "danger"}>
                Yield: {pct(step.expected_yield)}
              </Badge>
            </div>
            {step.named_reaction && (
              <span className="text-xs text-accent">{step.named_reaction}</span>
            )}
          </CardHeader>
          <div className="grid grid-cols-1 sm:grid-cols-2 gap-4 text-sm">
            <div>
              <span className="text-text-muted block mb-1">Precursors</span>
              {step.precursor_smiles.map((s) => (
                <code key={s} className="block text-xs font-mono text-text-secondary">{s}</code>
              ))}
            </div>
            <div>
              <span className="text-text-muted block mb-1">Product</span>
              <span className="font-medium">{step.product_name}</span>
              <code className="block text-xs font-mono text-text-secondary">{step.product_smiles}</code>
            </div>
          </div>
          {step.conditions && (
            <p className="text-xs text-text-muted mt-3">{step.conditions}</p>
          )}
          {step.notes && (
            <p className="text-xs text-text-secondary mt-1">{step.notes}</p>
          )}
        </Card>
      ))}
    </div>
  );
}
