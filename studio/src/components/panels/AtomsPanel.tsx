import { PanelContainer } from "./PanelContainer";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { useEditorStore } from "@/stores/editor-store";
import { cn } from "@/lib/cn";

export function AtomsPanel() {
  const activeId = useWorkspaceStore((s) => s.activeId);
  const molecules = useWorkspaceStore((s) => s.molecules);
  const entry = activeId ? molecules.get(activeId) : undefined;
  const structure = entry?.structure;
  const selectedAtoms = useEditorStore((s) => s.selectedAtoms);
  const selectAtom = useEditorStore((s) => s.selectAtom);
  const showHydrogens = useEditorStore((s) => s.showHydrogens);

  if (!structure) {
    return (
      <PanelContainer title="Atoms" defaultOpen={false}>
        <p className="text-xs text-text-muted">No molecule loaded</p>
      </PanelContainer>
    );
  }

  const atoms = showHydrogens
    ? structure.atoms
    : structure.atoms.filter((a) => a.symbol !== "H");

  return (
    <PanelContainer title={`Atoms (${atoms.length})`} defaultOpen={false}>
      <div className="max-h-48 overflow-y-auto">
        <table className="w-full text-[10px] font-mono">
          <thead>
            <tr className="text-text-muted">
              <th className="text-left py-0.5 pr-2">#</th>
              <th className="text-left py-0.5 pr-2">Sym</th>
              <th className="text-right py-0.5 pr-2">X</th>
              <th className="text-right py-0.5 pr-2">Y</th>
              <th className="text-right py-0.5">Z</th>
            </tr>
          </thead>
          <tbody>
            {atoms.map((atom) => (
              <tr
                key={atom.index}
                className={cn(
                  "cursor-pointer hover:bg-white/5 transition-colors",
                  selectedAtoms.has(atom.index) && "bg-accent/10",
                )}
                onClick={(e) => selectAtom(atom.index, e.shiftKey)}
              >
                <td className="py-0.5 pr-2 text-text-muted">{atom.index}</td>
                <td className="py-0.5 pr-2 text-text-primary">{atom.symbol}</td>
                <td className="py-0.5 pr-2 text-right text-text-secondary">
                  {atom.position[0].toFixed(3)}
                </td>
                <td className="py-0.5 pr-2 text-right text-text-secondary">
                  {atom.position[1].toFixed(3)}
                </td>
                <td className="py-0.5 text-right text-text-secondary">
                  {atom.position[2].toFixed(3)}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </PanelContainer>
  );
}
