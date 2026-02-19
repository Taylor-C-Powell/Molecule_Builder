import { useEditorStore } from "@/stores/editor-store";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { useHistoryStore } from "@/stores/history-store";
import { Tool } from "@/types/editor";
import { ElementPicker } from "@/components/panels/ElementPicker";
import { BondOrderPicker } from "@/components/panels/BondOrderPicker";
import { FragmentPicker } from "@/components/panels/FragmentPicker";
import type { FragmentDef } from "@/lib/fragments";
import { cn } from "@/lib/cn";

interface ToolButton {
  tool: Tool;
  icon: string;
  title: string;
}

const TOOLS: ToolButton[] = [
  { tool: Tool.SELECT, icon: "\u25B3", title: "Select (V)" },
  { tool: Tool.ROTATE, icon: "\u21BB", title: "Rotate (R)" },
  { tool: Tool.ADD_ATOM, icon: "+", title: "Add Atom (A)" },
  { tool: Tool.ADD_BOND, icon: "\u2014", title: "Add Bond (B)" },
  { tool: Tool.ERASE, icon: "\u2715", title: "Erase (E)" },
  { tool: Tool.MEASURE, icon: "\u2220", title: "Measure (M)" },
];

export function ToolPalette() {
  const activeTool = useEditorStore((s) => s.tool);
  const setTool = useEditorStore((s) => s.setTool);
  const selectedAtoms = useEditorStore((s) => s.selectedAtoms);

  function handleFragmentSelect(fragment: FragmentDef) {
    useEditorStore.getState().setActiveFragment(fragment);
    // If exactly one atom is selected, insert immediately
    if (selectedAtoms.size === 1) {
      const atomIdx = [...selectedAtoms][0]!;
      const entry = useWorkspaceStore.getState().getActive();
      if (entry?.structure) {
        useHistoryStore.getState().pushSnapshot(entry.structure, `Insert ${fragment.shortLabel}`);
        useWorkspaceStore.getState().insertFragment(atomIdx, fragment);
      }
    }
  }

  return (
    <div className="flex flex-col items-center gap-1 py-2 bg-bg-toolbar border-r border-border overflow-y-auto">
      {TOOLS.map((t) => (
        <button
          key={t.tool}
          className={cn(
            "w-9 h-9 flex items-center justify-center rounded-[var(--radius-sm)] text-xs font-bold transition-colors cursor-pointer",
            activeTool === t.tool
              ? "bg-accent text-white"
              : "text-text-muted hover:text-text-primary hover:bg-white/5",
          )}
          title={t.title}
          onClick={() => setTool(t.tool)}
        >
          {t.icon}
        </button>
      ))}

      {/* Contextual sub-tools */}
      {activeTool === Tool.ADD_ATOM && (
        <div className="mt-2 pt-2 border-t border-border px-1">
          <ElementPicker />
        </div>
      )}

      {activeTool === Tool.ADD_BOND && (
        <div className="mt-2 pt-2 border-t border-border px-1 w-full">
          <BondOrderPicker />
        </div>
      )}

      {activeTool === Tool.ADD_ATOM && (
        <div className="mt-1 pt-1 border-t border-border px-1">
          <div className="text-[8px] text-text-muted text-center mb-0.5">
            Fragments
          </div>
          <FragmentPicker onSelect={handleFragmentSelect} />
          {selectedAtoms.size !== 1 && (
            <div className="text-[8px] text-text-muted text-center mt-1">
              Select 1 atom first
            </div>
          )}
        </div>
      )}
    </div>
  );
}
