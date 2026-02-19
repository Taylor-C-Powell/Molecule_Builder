import { cn } from "@/lib/cn";

interface ShortcutsOverlayProps {
  open: boolean;
  onClose: () => void;
}

interface ShortcutGroup {
  label: string;
  shortcuts: { keys: string; desc: string }[];
}

const GROUPS: ShortcutGroup[] = [
  {
    label: "Tools",
    shortcuts: [
      { keys: "V", desc: "Select tool" },
      { keys: "R", desc: "Rotate (camera orbit)" },
      { keys: "A", desc: "Add Atom tool" },
      { keys: "B", desc: "Add Bond tool" },
      { keys: "E", desc: "Erase tool" },
      { keys: "M", desc: "Measure tool" },
    ],
  },
  {
    label: "Editing",
    shortcuts: [
      { keys: "Ctrl+Z", desc: "Undo" },
      { keys: "Ctrl+Y", desc: "Redo" },
      { keys: "Ctrl+A", desc: "Select all atoms" },
      { keys: "Ctrl+C", desc: "Copy SMILES" },
      { keys: "Ctrl+V", desc: "Paste SMILES" },
      { keys: "Delete", desc: "Delete selected" },
      { keys: "Escape", desc: "Clear selection" },
    ],
  },
  {
    label: "View",
    shortcuts: [
      { keys: "F", desc: "Fit to view" },
      { keys: "H", desc: "Toggle hydrogens" },
      { keys: "L", desc: "Toggle atom labels" },
      { keys: "D", desc: "Toggle distance labels" },
      { keys: "G", desc: "Toggle grid" },
      { keys: "?", desc: "Show this help" },
    ],
  },
  {
    label: "Selection",
    shortcuts: [
      { keys: "Click", desc: "Select atom or bond" },
      { keys: "Shift+Click", desc: "Multi-select" },
      { keys: "Click empty", desc: "Clear selection" },
    ],
  },
];

export function ShortcutsOverlay({ open, onClose }: ShortcutsOverlayProps) {
  if (!open) return null;

  return (
    <div
      className="fixed inset-0 z-[200] flex items-center justify-center bg-black/60"
      onClick={onClose}
    >
      <div
        className={cn(
          "bg-bg-card border border-border rounded-lg shadow-2xl p-6 max-w-lg w-full max-h-[80vh] overflow-y-auto",
          "animate-[fadeIn_0.15s_ease-out]",
        )}
        onClick={(e) => e.stopPropagation()}
      >
        <div className="flex items-center justify-between mb-4">
          <h2 className="text-sm font-bold text-text-primary">
            Keyboard Shortcuts
          </h2>
          <button
            className="text-text-muted hover:text-text-primary text-xs cursor-pointer"
            onClick={onClose}
          >
            Close (Esc)
          </button>
        </div>

        <div className="grid grid-cols-2 gap-4">
          {GROUPS.map((group) => (
            <div key={group.label}>
              <h3 className="text-[10px] font-bold text-accent uppercase tracking-wider mb-2">
                {group.label}
              </h3>
              <div className="space-y-1">
                {group.shortcuts.map((s) => (
                  <div
                    key={s.keys}
                    className="flex items-center justify-between gap-2"
                  >
                    <kbd className="px-1.5 py-0.5 bg-white/5 border border-border rounded text-[10px] font-mono text-text-secondary min-w-[60px] text-center">
                      {s.keys}
                    </kbd>
                    <span className="text-[10px] text-text-muted flex-1">
                      {s.desc}
                    </span>
                  </div>
                ))}
              </div>
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}
