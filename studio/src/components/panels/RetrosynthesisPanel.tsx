import { useState } from "react";
import { PanelContainer } from "./PanelContainer";
import { useRetrosynthesis } from "@/hooks/useRetrosynthesis";
import { useWorkspaceStore } from "@/stores/workspace-store";
import type { RetroNodeResponse } from "@/api/types";
import { cn } from "@/lib/cn";

function RetroNode({ node, depth }: { node: RetroNodeResponse; depth: number }) {
  const [expanded, setExpanded] = useState(depth < 2);
  const hasChildren = node.children.length > 0;

  return (
    <div className="ml-2 border-l border-border/50">
      <div
        className={cn(
          "flex items-start gap-1 px-2 py-1 text-[10px] font-mono",
          hasChildren && "cursor-pointer hover:bg-white/3",
        )}
        onClick={() => hasChildren && setExpanded(!expanded)}
      >
        {hasChildren && (
          <span className="text-text-muted shrink-0 mt-0.5">
            {expanded ? "\u25BC" : "\u25B6"}
          </span>
        )}
        <div className="min-w-0">
          <span
            className={cn(
              "break-all",
              node.is_purchasable ? "text-green-400" : "text-text-secondary",
            )}
          >
            {node.smiles}
          </span>
          {node.is_purchasable && (
            <span className="ml-1 text-green-500 text-[9px]">[buy]</span>
          )}
          {node.best_disconnection && (
            <div className="text-text-muted text-[9px]">
              {node.best_disconnection.reaction_name}
              {node.best_disconnection.named_reaction &&
                ` (${node.best_disconnection.named_reaction})`}
              <span className="ml-1 text-amber-400">
                score: {node.best_disconnection.score.toFixed(2)}
              </span>
            </div>
          )}
        </div>
      </div>
      {expanded &&
        node.children.map((child, i) => (
          <RetroNode key={i} node={child} depth={depth + 1} />
        ))}
    </div>
  );
}

export function RetrosynthesisPanel() {
  const activeId = useWorkspaceStore((s) => s.activeId);
  const molecules = useWorkspaceStore((s) => s.molecules);
  const entry = activeId ? molecules.get(activeId) : undefined;

  const { plan, loading, error, result, clear } = useRetrosynthesis();
  const [maxDepth, setMaxDepth] = useState(3);
  const [beamWidth, setBeamWidth] = useState(3);

  const smiles = entry?.smiles;

  async function handleRun() {
    if (!smiles) return;
    await plan({ smiles, max_depth: maxDepth, beam_width: beamWidth });
  }

  return (
    <PanelContainer title="Retrosynthesis" defaultOpen={false}>
      {!smiles ? (
        <div className="text-[10px] text-text-muted py-2">
          Load a molecule first
        </div>
      ) : (
        <div className="space-y-2">
          {/* Parameters */}
          <div className="flex gap-2">
            <label className="flex items-center gap-1 text-[10px] text-text-muted">
              Depth
              <input
                type="number"
                min={1}
                max={10}
                value={maxDepth}
                onChange={(e) => setMaxDepth(Number(e.target.value))}
                className="w-10 px-1 py-0.5 bg-bg border border-border rounded-[var(--radius-sm)] text-[10px] text-text-primary text-center"
              />
            </label>
            <label className="flex items-center gap-1 text-[10px] text-text-muted">
              Beam
              <input
                type="number"
                min={1}
                max={10}
                value={beamWidth}
                onChange={(e) => setBeamWidth(Number(e.target.value))}
                className="w-10 px-1 py-0.5 bg-bg border border-border rounded-[var(--radius-sm)] text-[10px] text-text-primary text-center"
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
              {loading ? "Planning..." : "Plan Route"}
            </button>
          </div>

          {error && (
            <div className="text-[10px] text-red-400 bg-red-400/10 px-2 py-1 rounded-[var(--radius-sm)]">
              {error}
            </div>
          )}

          {result && (
            <div className="space-y-2">
              {/* Summary badges */}
              <div className="flex gap-1 flex-wrap">
                <span className="px-1.5 py-0.5 bg-accent/10 text-accent text-[9px] rounded-[var(--radius-sm)]">
                  {result.routes_found} routes
                </span>
                <span className="px-1.5 py-0.5 bg-white/5 text-text-muted text-[9px] rounded-[var(--radius-sm)]">
                  depth {result.max_depth}
                </span>
                <span className="px-1.5 py-0.5 bg-white/5 text-text-muted text-[9px] rounded-[var(--radius-sm)]">
                  beam {result.beam_width}
                </span>
              </div>

              {/* Best route summary */}
              {result.best_route && (
                <div className="bg-green-900/20 border border-green-500/20 rounded-[var(--radius-sm)] p-2 space-y-1">
                  <div className="text-[10px] font-bold text-green-400">
                    Best Route: {result.best_route.total_steps} steps,{" "}
                    {(result.best_route.overall_yield * 100).toFixed(1)}% yield
                  </div>
                  <div className="text-[9px] text-text-muted">
                    Starting materials:{" "}
                    {result.best_route.starting_materials.map((m) => m.name).join(", ")}
                  </div>
                  {result.best_route.steps.map((step) => (
                    <div
                      key={step.step_number}
                      className="pl-2 border-l-2 border-green-500/30 text-[9px]"
                    >
                      <span className="text-green-300 font-mono">
                        Step {step.step_number}:
                      </span>{" "}
                      <span className="text-text-secondary">
                        {step.reaction_name}
                      </span>
                      {step.named_reaction && (
                        <span className="text-amber-400 ml-1">
                          ({step.named_reaction})
                        </span>
                      )}
                      <div className="text-text-muted">
                        {step.product_name} | yield:{" "}
                        {(step.expected_yield * 100).toFixed(0)}%
                      </div>
                      {step.conditions && (
                        <div className="text-text-muted italic">
                          {step.conditions}
                        </div>
                      )}
                    </div>
                  ))}
                </div>
              )}

              {/* Tree view */}
              <div className="text-[10px] text-text-muted font-bold mt-1">
                Retrosynthetic Tree
              </div>
              <div className="max-h-[300px] overflow-y-auto bg-bg/50 rounded-[var(--radius-sm)] border border-border/50 py-1">
                <RetroNode node={result.tree} depth={0} />
              </div>

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
