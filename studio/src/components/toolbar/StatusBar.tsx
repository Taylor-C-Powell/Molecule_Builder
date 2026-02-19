import { useWorkspaceStore } from "@/stores/workspace-store";
import { useEditorStore } from "@/stores/editor-store";
import { useHistoryStore } from "@/stores/history-store";
import { Tool } from "@/types/editor";

const TOOL_LABELS: Record<Tool, string> = {
  [Tool.SELECT]: "Select",
  [Tool.ROTATE]: "Rotate",
  [Tool.ADD_ATOM]: "Add Atom",
  [Tool.ADD_BOND]: "Add Bond",
  [Tool.ERASE]: "Erase",
  [Tool.MEASURE]: "Measure",
};

export function StatusBar() {
  const activeId = useWorkspaceStore((s) => s.activeId);
  const molecules = useWorkspaceStore((s) => s.molecules);
  const entry = activeId ? molecules.get(activeId) : undefined;
  const tool = useEditorStore((s) => s.tool);
  const renderStyle = useEditorStore((s) => s.renderStyle);
  const activeElement = useEditorStore((s) => s.activeElement);
  const activeBondOrder = useEditorStore((s) => s.activeBondOrder);
  const selectedAtoms = useEditorStore((s) => s.selectedAtoms);
  const selectedBonds = useEditorStore((s) => s.selectedBonds);
  const measureAtoms = useEditorStore((s) => s.measureAtoms);
  const canUndo = useHistoryStore((s) => s.canUndo);
  const canRedo = useHistoryStore((s) => s.canRedo);

  const bondLabels = ["", "Single", "Double", "Triple"];

  return (
    <div className="flex items-center h-7 px-3 bg-bg-toolbar border-t border-border text-[10px] font-mono text-text-muted gap-4 select-none">
      <span>
        {entry ? entry.name || entry.smiles : "No molecule"}
      </span>
      {entry?.molecule && (
        <>
          <span>{entry.molecule.num_atoms} atoms</span>
          <span>{entry.molecule.num_bonds} bonds</span>
        </>
      )}
      {entry?.properties?.formula && (
        <span>{entry.properties.formula}</span>
      )}

      <div className="flex-1" />

      {/* Selection info */}
      {selectedAtoms.size > 0 && (
        <span className="text-accent">
          {selectedAtoms.size} atom{selectedAtoms.size !== 1 ? "s" : ""} selected
        </span>
      )}
      {selectedBonds.size > 0 && (
        <span className="text-accent">
          {selectedBonds.size} bond{selectedBonds.size !== 1 ? "s" : ""} selected
        </span>
      )}

      {/* Tool context */}
      {tool === Tool.ADD_ATOM && (
        <span className="text-green-400">Element: {activeElement}</span>
      )}
      {tool === Tool.ADD_BOND && (
        <span className="text-blue-400">Bond: {bondLabels[activeBondOrder]}</span>
      )}
      {tool === Tool.MEASURE && measureAtoms.length > 0 && (
        <span className="text-amber-400">
          Measuring: {measureAtoms.length}/4 atoms
        </span>
      )}

      {/* Undo/Redo indicator */}
      {(canUndo || canRedo) && (
        <span>
          {canUndo ? "U" : "-"}/{canRedo ? "R" : "-"}
        </span>
      )}

      <span>Tool: {TOOL_LABELS[tool]}</span>
      <span>Style: {renderStyle}</span>
    </div>
  );
}
