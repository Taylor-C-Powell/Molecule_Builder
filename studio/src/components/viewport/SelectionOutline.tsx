import { useEditorStore } from "@/stores/editor-store";
import { useWorkspaceStore } from "@/stores/workspace-store";

export function SelectionOutline() {
  const selectedAtoms = useEditorStore((s) => s.selectedAtoms);
  const activeId = useWorkspaceStore((s) => s.activeId);
  const molecules = useWorkspaceStore((s) => s.molecules);
  const entry = activeId ? molecules.get(activeId) : undefined;
  const structure = entry?.structure;

  if (!structure || selectedAtoms.size === 0) return null;

  return (
    <div className="absolute bottom-2 left-2 bg-bg/80 border border-accent/30 rounded-[var(--radius-sm)] px-3 py-2 pointer-events-none">
      <div className="text-[10px] text-accent font-mono">
        {selectedAtoms.size} atom{selectedAtoms.size !== 1 ? "s" : ""} selected
      </div>
      <div className="text-[10px] text-text-muted font-mono">
        {Array.from(selectedAtoms)
          .slice(0, 8)
          .map((idx) => {
            const atom = structure.atoms[idx];
            return atom ? `${atom.symbol}${idx}` : `#${idx}`;
          })
          .join(", ")}
        {selectedAtoms.size > 8 && "..."}
      </div>
    </div>
  );
}
