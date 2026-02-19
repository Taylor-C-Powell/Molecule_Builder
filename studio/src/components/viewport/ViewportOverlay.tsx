import { useEditorStore } from "@/stores/editor-store";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { Tool } from "@/types/editor";

export function ViewportOverlay() {
  const hoveredAtom = useEditorStore((s) => s.hoveredAtom);
  const tool = useEditorStore((s) => s.tool);
  const bondStartAtom = useEditorStore((s) => s.bondStartAtom);
  const measureAtoms = useEditorStore((s) => s.measureAtoms);
  const activeId = useWorkspaceStore((s) => s.activeId);
  const molecules = useWorkspaceStore((s) => s.molecules);
  const entry = activeId ? molecules.get(activeId) : undefined;
  const structure = entry?.structure;

  return (
    <>
      {/* Hovered atom info */}
      {hoveredAtom != null && structure && (() => {
        const atom = structure.atoms[hoveredAtom];
        if (!atom) return null;
        return (
          <div className="absolute top-2 left-2 bg-bg/80 border border-border rounded-[var(--radius-sm)] px-3 py-2 pointer-events-none">
            <div className="text-xs font-mono text-text-primary">
              <span className="text-accent font-bold">{atom.symbol}</span>
              <span className="text-text-muted ml-1">#{atom.index}</span>
            </div>
            <div className="text-[10px] text-text-muted font-mono">
              ({atom.position[0].toFixed(3)}, {atom.position[1].toFixed(3)},{" "}
              {atom.position[2].toFixed(3)})
            </div>
            {atom.hybridization && (
              <div className="text-[10px] text-text-muted">
                {atom.hybridization}
                {atom.formal_charge !== 0 && (
                  <span className="ml-1">
                    charge: {atom.formal_charge > 0 ? "+" : ""}
                    {atom.formal_charge}
                  </span>
                )}
              </div>
            )}
          </div>
        );
      })()}

      {/* Tool hint */}
      {tool === Tool.ADD_BOND && bondStartAtom != null && structure && (
        <div className="absolute top-2 right-2 bg-green-900/60 border border-green-500/30 rounded-[var(--radius-sm)] px-3 py-1.5 pointer-events-none">
          <div className="text-[10px] text-green-300 font-mono">
            Bond from {structure.atoms[bondStartAtom]?.symbol}#{bondStartAtom}
            {" -> click target atom"}
          </div>
        </div>
      )}

      {/* Measure mode hint */}
      {tool === Tool.MEASURE && (
        <div className="absolute top-2 right-2 bg-amber-900/60 border border-amber-500/30 rounded-[var(--radius-sm)] px-3 py-1.5 pointer-events-none">
          <div className="text-[10px] text-amber-300 font-mono">
            {measureAtoms.length === 0 && "Click atoms to measure (2=dist, 3=angle, 4=dihedral)"}
            {measureAtoms.length === 1 && "Click 2nd atom for distance..."}
            {measureAtoms.length === 2 && "Distance shown. Click 3rd for angle..."}
            {measureAtoms.length === 3 && "Angle shown. Click 4th for dihedral..."}
          </div>
        </div>
      )}
    </>
  );
}
