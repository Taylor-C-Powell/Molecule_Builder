import { useEffect } from "react";
import { useEditorStore } from "@/stores/editor-store";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { useHistoryStore } from "@/stores/history-store";
import { useToastStore } from "@/stores/toast-store";
import { addMissingHydrogens } from "@/lib/hydrogenate";
import { Tool } from "@/types/editor";

interface UseKeyboardShortcutsOptions {
  onShowHelp?: () => void;
  onPasteSmiles?: (smiles: string) => void;
}

export function useKeyboardShortcuts(options?: UseKeyboardShortcutsOptions) {
  const setTool = useEditorStore((s) => s.setTool);
  const clearSelection = useEditorStore((s) => s.clearSelection);
  const clearMeasurements = useEditorStore((s) => s.clearMeasurements);
  const toggleHydrogens = useEditorStore((s) => s.toggleHydrogens);
  const toggleLabels = useEditorStore((s) => s.toggleLabels);

  useEffect(() => {
    function handleKeyDown(e: KeyboardEvent) {
      // Ignore when typing in input fields
      const tag = (e.target as HTMLElement).tagName;
      if (tag === "INPUT" || tag === "TEXTAREA" || tag === "SELECT") return;

      const ctrl = e.ctrlKey || e.metaKey;

      // Add Hydrogens: Ctrl+H
      if (ctrl && e.key === "h") {
        e.preventDefault();
        const entry = useWorkspaceStore.getState().getActive();
        if (entry?.structure) {
          useHistoryStore.getState().pushSnapshot(entry.structure, "Add hydrogens");
          const filled = addMissingHydrogens(entry.structure);
          useWorkspaceStore.getState().updateActiveStructure(filled);
          useToastStore.getState().addToast("Hydrogens added", "success");
        }
        return;
      }

      // Copy SMILES: Ctrl+C
      if (ctrl && e.key === "c" && !e.shiftKey) {
        const entry = useWorkspaceStore.getState().getActive();
        if (entry?.smiles) {
          e.preventDefault();
          navigator.clipboard.writeText(entry.smiles).then(() => {
            useToastStore.getState().addToast("SMILES copied to clipboard", "success");
          });
          return;
        }
      }

      // Paste SMILES: Ctrl+V
      if (ctrl && e.key === "v" && !e.shiftKey) {
        e.preventDefault();
        navigator.clipboard.readText().then((text) => {
          const trimmed = text.trim();
          if (trimmed && options?.onPasteSmiles) {
            options.onPasteSmiles(trimmed);
          }
        });
        return;
      }

      // Undo: Ctrl+Z
      if (ctrl && e.key === "z" && !e.shiftKey) {
        e.preventDefault();
        const entry = useWorkspaceStore.getState().getActive();
        if (!entry?.structure) return;
        const restored = useHistoryStore.getState().undo(entry.structure);
        if (restored) {
          useWorkspaceStore.getState().updateActiveStructure(restored);
        }
        return;
      }

      // Redo: Ctrl+Shift+Z or Ctrl+Y
      if ((ctrl && e.key === "z" && e.shiftKey) || (ctrl && e.key === "y")) {
        e.preventDefault();
        const entry = useWorkspaceStore.getState().getActive();
        if (!entry?.structure) return;
        const restored = useHistoryStore.getState().redo(entry.structure);
        if (restored) {
          useWorkspaceStore.getState().updateActiveStructure(restored);
        }
        return;
      }

      // Select All: Ctrl+A
      if (ctrl && e.key === "a") {
        e.preventDefault();
        const entry = useWorkspaceStore.getState().getActive();
        if (entry?.structure) {
          useEditorStore.getState().selectAllAtoms(entry.structure.atoms.length);
        }
        return;
      }

      // Don't process tool shortcuts when ctrl is held
      if (ctrl) return;

      switch (e.key) {
        case "?":
          options?.onShowHelp?.();
          break;
        case "Escape":
          clearSelection();
          clearMeasurements();
          break;
        case "Delete":
        case "Backspace": {
          e.preventDefault();
          const editor = useEditorStore.getState();
          const ws = useWorkspaceStore.getState();
          const active = ws.getActive();
          if (!active?.structure) break;
          if (editor.selectedAtoms.size > 0) {
            useHistoryStore.getState().pushSnapshot(active.structure, "Delete atoms");
            ws.removeAtoms(editor.selectedAtoms);
            clearSelection();
          } else if (editor.selectedBonds.size > 0) {
            useHistoryStore.getState().pushSnapshot(active.structure, "Delete bonds");
            ws.removeBonds(editor.selectedBonds);
            clearSelection();
          }
          break;
        }
        default:
          // Tool shortcuts use lowercase
          switch (e.key.toLowerCase()) {
            case "v":
              setTool(Tool.SELECT);
              break;
            case "r":
              setTool(Tool.ROTATE);
              break;
            case "a":
              setTool(Tool.ADD_ATOM);
              break;
            case "b":
              setTool(Tool.ADD_BOND);
              break;
            case "e":
              setTool(Tool.ERASE);
              break;
            case "m":
              setTool(Tool.MEASURE);
              break;
            case "d":
              useEditorStore.getState().toggleDistanceLabels();
              break;
            case "g":
              useEditorStore.getState().toggleGrid();
              break;
            case "f": {
              // Fit to view: reset camera via custom event
              const gl = document.querySelector("canvas");
              if (gl) {
                gl.dispatchEvent(
                  new CustomEvent("molbuilder:camera", {
                    detail: { position: [0, 0, 15] },
                  }),
                );
              }
              break;
            }
            case "h":
              toggleHydrogens();
              break;
            case "l":
              toggleLabels();
              break;
          }
      }
    }

    window.addEventListener("keydown", handleKeyDown);
    return () => window.removeEventListener("keydown", handleKeyDown);
  }, [setTool, clearSelection, clearMeasurements, toggleHydrogens, toggleLabels, options]);
}
