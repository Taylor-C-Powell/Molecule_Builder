import { useState } from "react";
import { PanelContainer } from "./PanelContainer";
import { useProcess } from "@/hooks/useProcess";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { cn } from "@/lib/cn";

function RiskBadge({ level }: { level: string }) {
  const colors: Record<string, string> = {
    HIGH: "bg-red-500/20 text-red-400",
    MEDIUM: "bg-amber-500/20 text-amber-400",
    LOW: "bg-green-500/20 text-green-400",
  };
  return (
    <span
      className={cn(
        "px-1.5 py-0.5 text-[9px] font-bold rounded-[var(--radius-sm)]",
        colors[level] ?? "bg-white/10 text-text-muted",
      )}
    >
      {level}
    </span>
  );
}

export function ProcessPanel() {
  const activeId = useWorkspaceStore((s) => s.activeId);
  const molecules = useWorkspaceStore((s) => s.molecules);
  const entry = activeId ? molecules.get(activeId) : undefined;

  const { evaluate, loading, error, result, clear } = useProcess();
  const [scaleKg, setScaleKg] = useState(1.0);

  const smiles = entry?.smiles;

  async function handleRun() {
    if (!smiles) return;
    await evaluate({ smiles, scale_kg: scaleKg });
  }

  return (
    <PanelContainer title="Process Engineering" defaultOpen={false}>
      {!smiles ? (
        <div className="text-[10px] text-text-muted py-2">
          Load a molecule first
        </div>
      ) : (
        <div className="space-y-2">
          <div className="flex gap-2 items-center">
            <label className="flex items-center gap-1 text-[10px] text-text-muted">
              Scale (kg)
              <input
                type="number"
                min={0.001}
                step={0.1}
                value={scaleKg}
                onChange={(e) => setScaleKg(Number(e.target.value))}
                className="w-16 px-1 py-0.5 bg-bg border border-border rounded-[var(--radius-sm)] text-[10px] text-text-primary text-center"
              />
            </label>
            <button
              className={cn(
                "flex-1 py-1 text-[10px] font-medium rounded-[var(--radius-sm)] transition-colors cursor-pointer",
                loading
                  ? "bg-accent/30 text-accent/60 cursor-wait"
                  : "bg-accent text-white hover:bg-accent/90",
              )}
              disabled={loading}
              onClick={handleRun}
            >
              {loading ? "Evaluating..." : "Evaluate"}
            </button>
          </div>

          {error && (
            <div className="text-[10px] text-red-400 bg-red-400/10 px-2 py-1 rounded-[var(--radius-sm)]">
              {error}
            </div>
          )}

          {result && (
            <div className="space-y-2">
              {!result.route_found ? (
                <div className="text-[10px] text-amber-400">
                  No synthesis route found for this molecule.
                </div>
              ) : (
                <>
                  {/* Summary */}
                  <div className="flex gap-1 flex-wrap">
                    <span className="px-1.5 py-0.5 bg-accent/10 text-accent text-[9px] rounded-[var(--radius-sm)]">
                      {result.total_steps} steps
                    </span>
                    <span className="px-1.5 py-0.5 bg-green-500/10 text-green-400 text-[9px] rounded-[var(--radius-sm)]">
                      {(result.overall_yield * 100).toFixed(1)}% yield
                    </span>
                    <span className="px-1.5 py-0.5 bg-white/5 text-text-muted text-[9px] rounded-[var(--radius-sm)]">
                      {result.scale_kg} kg
                    </span>
                  </div>

                  {/* Step Details */}
                  {result.step_details.map((step) => (
                    <div
                      key={step.step_number}
                      className="bg-white/3 rounded-[var(--radius-sm)] p-2 space-y-1"
                    >
                      <div className="text-[10px] font-bold text-text-primary">
                        Step {step.step_number}: {step.reaction_name}
                      </div>
                      <div className="text-[9px] text-text-muted grid grid-cols-2 gap-x-2 gap-y-0.5">
                        <span>Reactor: {step.reactor.reactor_type}</span>
                        <span>Volume: {step.reactor.volume_L}L</span>
                        <span>Temp: {step.conditions.temperature_C}C</span>
                        <span>Time: {step.conditions.reaction_time_hours}h</span>
                        <span>Solvent: {step.conditions.solvent}</span>
                        <span>Atm: {step.conditions.atmosphere}</span>
                      </div>
                      {step.purification.length > 0 && (
                        <div className="text-[9px] text-text-muted">
                          Purification:{" "}
                          {step.purification.map((p) => p.method).join(", ")}
                        </div>
                      )}
                    </div>
                  ))}

                  {/* Safety */}
                  {result.safety.length > 0 && (
                    <div className="space-y-1">
                      <div className="text-[10px] font-bold text-text-muted">
                        Safety Assessment
                      </div>
                      {result.safety.map((s) => (
                        <div
                          key={s.step_number}
                          className="flex items-start gap-1 text-[9px]"
                        >
                          <RiskBadge level={s.risk_level} />
                          <div>
                            <span className="text-text-secondary">
                              Step {s.step_number}: {s.step_name}
                            </span>
                            {s.ppe_required.length > 0 && (
                              <div className="text-text-muted">
                                PPE: {s.ppe_required.join(", ")}
                              </div>
                            )}
                          </div>
                        </div>
                      ))}
                    </div>
                  )}

                  {/* Cost */}
                  {result.cost && (
                    <div className="bg-white/3 rounded-[var(--radius-sm)] p-2">
                      <div className="text-[10px] font-bold text-text-muted mb-1">
                        Cost Estimate
                      </div>
                      <div className="text-[11px] font-bold text-accent">
                        ${result.cost.total_usd.toFixed(2)} total (${result.cost.per_kg_usd.toFixed(2)}/kg)
                      </div>
                      <div className="text-[9px] text-text-muted grid grid-cols-2 gap-x-2 gap-y-0.5 mt-1">
                        <span>Materials: ${result.cost.breakdown.raw_materials_usd.toFixed(2)}</span>
                        <span>Labor: ${result.cost.breakdown.labor_usd.toFixed(2)}</span>
                        <span>Equipment: ${result.cost.breakdown.equipment_usd.toFixed(2)}</span>
                        <span>Energy: ${result.cost.breakdown.energy_usd.toFixed(2)}</span>
                      </div>
                    </div>
                  )}

                  {/* Scale-up */}
                  {result.scale_up && (
                    <div className="bg-white/3 rounded-[var(--radius-sm)] p-2">
                      <div className="text-[10px] font-bold text-text-muted mb-1">
                        Scale-Up Assessment
                      </div>
                      <div className="text-[9px] text-text-muted grid grid-cols-2 gap-x-2 gap-y-0.5">
                        <span>Mode: {result.scale_up.recommended_mode}</span>
                        <span>Capacity: {result.scale_up.annual_capacity_kg} kg/yr</span>
                        <span>CapEx: ${result.scale_up.capital_cost_usd.toFixed(0)}</span>
                        <span>OpEx: ${result.scale_up.operating_cost_annual_usd.toFixed(0)}/yr</span>
                      </div>
                      {result.scale_up.scale_up_risks.length > 0 && (
                        <div className="mt-1 text-[9px] text-amber-400">
                          Risks: {result.scale_up.scale_up_risks.join("; ")}
                        </div>
                      )}
                    </div>
                  )}
                </>
              )}

              <button
                className="w-full py-1 text-[10px] text-text-muted hover:text-text-primary bg-white/5 hover:bg-white/10 rounded-[var(--radius-sm)] transition-colors cursor-pointer"
                onClick={clear}
              >
                Clear Results
              </button>
            </div>
          )}
        </div>
      )}
    </PanelContainer>
  );
}
