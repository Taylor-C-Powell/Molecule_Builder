import { useWorkspaceStore } from "@/stores/workspace-store";
import { useHistoryStore } from "@/stores/history-store";
import { cn } from "@/lib/cn";

export function WorkspaceTabs() {
  const molecules = useWorkspaceStore((s) => s.molecules);
  const activeId = useWorkspaceStore((s) => s.activeId);
  const setActive = useWorkspaceStore((s) => s.setActive);
  const removeMolecule = useWorkspaceStore((s) => s.removeMolecule);
  const clearHistory = useHistoryStore((s) => s.clear);

  const entries = Array.from(molecules.values());

  if (entries.length <= 1) return null;

  return (
    <div className="flex items-center h-7 bg-bg-toolbar/50 border-b border-border px-1 gap-0.5 overflow-x-auto select-none">
      {entries.map((entry) => {
        const isActive = entry.id === activeId;
        const label = entry.name || entry.smiles;
        const displayLabel = label.length > 20 ? label.slice(0, 18) + "..." : label;

        return (
          <div
            key={entry.id}
            className={cn(
              "flex items-center gap-1 px-2 py-0.5 rounded-t-[var(--radius-sm)] text-[10px] font-mono cursor-pointer transition-colors min-w-0 shrink-0",
              isActive
                ? "bg-bg text-text-primary border-t border-x border-border"
                : "text-text-muted hover:text-text-secondary hover:bg-white/3",
            )}
            onClick={() => {
              if (!isActive) {
                clearHistory();
                setActive(entry.id);
              }
            }}
          >
            <span className="truncate max-w-[120px]" title={label}>
              {displayLabel}
            </span>
            <button
              className="ml-0.5 w-3.5 h-3.5 flex items-center justify-center rounded-[2px] text-[8px] text-text-muted hover:text-red-400 hover:bg-red-400/10 transition-colors cursor-pointer"
              title="Close"
              onClick={(e) => {
                e.stopPropagation();
                removeMolecule(entry.id);
              }}
            >
              x
            </button>
          </div>
        );
      })}
    </div>
  );
}
