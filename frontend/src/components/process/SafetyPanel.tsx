import type { SafetyAssessmentResponse } from "@/api/types";
import { Card, CardHeader, CardTitle } from "@/components/ui/Card";
import { Badge } from "@/components/ui/Badge";

interface SafetyPanelProps {
  safety: SafetyAssessmentResponse[];
}

function riskBadge(level: string) {
  const l = level.toLowerCase();
  if (l.includes("high") || l.includes("critical")) return "danger" as const;
  if (l.includes("medium") || l.includes("moderate")) return "warning" as const;
  return "success" as const;
}

export function SafetyPanel({ safety }: SafetyPanelProps) {
  if (!safety.length) return null;

  return (
    <Card>
      <CardHeader>
        <CardTitle>Safety Assessment</CardTitle>
      </CardHeader>
      <div className="space-y-4">
        {safety.map((s) => (
          <div key={s.step_number} className="border-b border-border/50 pb-4 last:border-0 last:pb-0">
            <div className="flex items-center justify-between mb-2">
              <span className="font-medium text-sm">Step {s.step_number}: {s.step_name}</span>
              <Badge variant={riskBadge(s.risk_level)}>{s.risk_level}</Badge>
            </div>
            <div className="grid grid-cols-1 sm:grid-cols-2 gap-3 text-sm">
              {s.ppe_required.length > 0 && (
                <div>
                  <span className="text-text-muted text-xs block mb-1">PPE Required</span>
                  <div className="flex flex-wrap gap-1">
                    {s.ppe_required.map((p) => (
                      <Badge key={p} variant="default">{p}</Badge>
                    ))}
                  </div>
                </div>
              )}
              {s.engineering_controls.length > 0 && (
                <div>
                  <span className="text-text-muted text-xs block mb-1">Engineering Controls</span>
                  <ul className="text-xs text-text-secondary space-y-0.5">
                    {s.engineering_controls.map((c) => (
                      <li key={c}>{c}</li>
                    ))}
                  </ul>
                </div>
              )}
            </div>
            {s.waste_classification && (
              <p className="text-xs text-text-muted mt-2">Waste: {s.waste_classification}</p>
            )}
          </div>
        ))}
      </div>
    </Card>
  );
}
