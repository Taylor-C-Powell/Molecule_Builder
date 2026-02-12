import type { StepProcessDetail } from "@/api/types";
import { Card, CardHeader, CardTitle } from "@/components/ui/Card";
import { num, usd } from "@/lib/format";

interface ReactorCardProps {
  step: StepProcessDetail;
}

export function ReactorCard({ step }: ReactorCardProps) {
  const r = step.reactor;
  const c = step.conditions;

  return (
    <Card>
      <CardHeader>
        <CardTitle>Step {step.step_number}: {step.reaction_name}</CardTitle>
      </CardHeader>
      <div className="grid grid-cols-1 sm:grid-cols-2 gap-6 text-sm">
        {/* Reactor */}
        <div>
          <h4 className="font-medium text-text-secondary mb-2">Reactor</h4>
          <div className="space-y-1">
            <Row label="Type" value={r.reactor_type} />
            <Row label="Volume" value={`${num(r.volume_L)} L`} />
            <Row label="Temperature" value={`${num(r.temperature_C)} C`} />
            <Row label="Pressure" value={`${num(r.pressure_atm)} atm`} />
            <Row label="Residence Time" value={`${num(r.residence_time_min)} min`} />
            <Row label="Mixing" value={r.mixing_type} />
            <Row label="Material" value={r.material} />
            <Row label="Est. Cost" value={usd(r.estimated_cost_usd)} />
          </div>
        </div>
        {/* Conditions */}
        <div>
          <h4 className="font-medium text-text-secondary mb-2">Conditions</h4>
          <div className="space-y-1">
            <Row label="Solvent" value={c.solvent} />
            <Row label="Concentration" value={`${num(c.concentration_M)} M`} />
            <Row label="Addition Rate" value={c.addition_rate} />
            <Row label="Reaction Time" value={`${num(c.reaction_time_hours)} h`} />
            <Row label="Atmosphere" value={c.atmosphere} />
          </div>
          {c.workup_procedure && (
            <p className="text-xs text-text-muted mt-2">{c.workup_procedure}</p>
          )}
        </div>
      </div>
      {/* Purification */}
      {step.purification.length > 0 && (
        <div className="mt-4 pt-4 border-t border-border/50">
          <h4 className="font-medium text-text-secondary mb-2 text-sm">Purification</h4>
          {step.purification.map((p, i) => (
            <div key={i} className="text-sm mb-2">
              <span className="font-medium">{p.method}</span>
              <span className="text-text-muted ml-2">
                Recovery: {num(p.estimated_recovery * 100)}% | Purity: {num(p.estimated_purity * 100)}%
              </span>
              <p className="text-xs text-text-muted">{p.description}</p>
            </div>
          ))}
        </div>
      )}
    </Card>
  );
}

function Row({ label, value }: { label: string; value: string }) {
  return (
    <div className="flex justify-between py-0.5">
      <span className="text-text-muted">{label}</span>
      <span className="font-mono">{value}</span>
    </div>
  );
}
