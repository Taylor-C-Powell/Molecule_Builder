import type { DisconnectionResponse } from "@/api/types";

interface TreeNode {
  smiles: string;
  is_purchasable: boolean;
  functional_groups: string[];
  reaction?: string;
  score?: number;
  disconnections?: DisconnectionResponse[];
}

interface NodeDetailProps {
  node: TreeNode;
  onClose: () => void;
}

export function NodeDetail({ node, onClose }: NodeDetailProps) {
  const best = node.disconnections?.find(
    (d) => d.reaction_name === node.reaction,
  );
  const alternatives = node.disconnections?.filter(
    (d) => d.reaction_name !== node.reaction,
  );

  return (
    <div className="absolute right-0 top-0 bottom-0 w-80 bg-bg-card border-l border-border overflow-y-auto p-4 shadow-lg">
      <div className="flex items-center justify-between mb-4">
        <h3 className="text-sm font-semibold text-text-primary">
          Node Detail
        </h3>
        <button
          onClick={onClose}
          className="text-text-secondary hover:text-text-primary text-lg leading-none"
          aria-label="Close detail panel"
        >
          &times;
        </button>
      </div>

      {/* SMILES */}
      <div className="mb-3">
        <span className="text-xs text-text-secondary">SMILES</span>
        <p className="text-xs font-mono text-text-primary break-all mt-0.5">
          {node.smiles}
        </p>
      </div>

      {/* Purchasable status */}
      <div className="mb-3">
        <span className="text-xs text-text-secondary">Status</span>
        <p className="text-xs mt-0.5">
          {node.is_purchasable ? (
            <span className="text-green-400 font-medium">Purchasable</span>
          ) : (
            <span className="text-text-secondary">Intermediate</span>
          )}
        </p>
      </div>

      {/* Functional groups */}
      {node.functional_groups.length > 0 && (
        <div className="mb-3">
          <span className="text-xs text-text-secondary">
            Functional Groups
          </span>
          <div className="flex flex-wrap gap-1 mt-1">
            {node.functional_groups.map((fg) => (
              <span
                key={fg}
                className="text-[10px] px-1.5 py-0.5 rounded bg-white/5 text-text-secondary border border-border"
              >
                {fg}
              </span>
            ))}
          </div>
        </div>
      )}

      {/* Best disconnection */}
      {best && (
        <div className="mb-3">
          <span className="text-xs text-text-secondary">
            Best Disconnection
          </span>
          <div className="mt-1 p-2 rounded border border-border bg-white/[0.02]">
            <div className="flex items-center justify-between">
              <span className="text-xs font-medium text-text-primary">
                {best.reaction_name}
              </span>
              <span
                className="text-xs font-semibold"
                style={{
                  color:
                    best.score >= 70
                      ? "#22c55e"
                      : best.score >= 40
                        ? "#eab308"
                        : "#ef4444",
                }}
              >
                {best.score.toFixed(1)}
              </span>
            </div>
            {best.named_reaction && (
              <p className="text-[10px] text-text-secondary mt-0.5">
                {best.named_reaction}
              </p>
            )}
            <p className="text-[10px] text-text-secondary mt-0.5">
              {best.category}
            </p>
            {best.precursors.length > 0 && (
              <div className="mt-1.5 space-y-0.5">
                {best.precursors.map((p, i) => (
                  <div key={i} className="text-[10px] font-mono text-text-secondary">
                    {p.name || p.smiles}
                    {p.cost_per_kg > 0 && (
                      <span className="text-text-secondary ml-1">
                        ${p.cost_per_kg.toFixed(0)}/kg
                      </span>
                    )}
                  </div>
                ))}
              </div>
            )}
          </div>
        </div>
      )}

      {/* Alternative disconnections */}
      {alternatives && alternatives.length > 0 && (
        <div>
          <span className="text-xs text-text-secondary">
            Alternative Disconnections ({alternatives.length})
          </span>
          <div className="mt-1 space-y-1.5">
            {alternatives.map((alt, i) => (
              <div
                key={i}
                className="p-2 rounded border border-border bg-white/[0.01]"
              >
                <div className="flex items-center justify-between">
                  <span className="text-[11px] text-text-primary">
                    {alt.reaction_name}
                  </span>
                  <span
                    className="text-[10px] font-semibold"
                    style={{
                      color:
                        alt.score >= 70
                          ? "#22c55e"
                          : alt.score >= 40
                            ? "#eab308"
                            : "#ef4444",
                    }}
                  >
                    {alt.score.toFixed(1)}
                  </span>
                </div>
                {alt.named_reaction && (
                  <p className="text-[10px] text-text-secondary mt-0.5">
                    {alt.named_reaction}
                  </p>
                )}
                <p className="text-[10px] text-text-secondary">
                  {alt.category}
                </p>
              </div>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}
