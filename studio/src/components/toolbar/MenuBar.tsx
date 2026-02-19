import { useState, useRef, useEffect } from "react";
import { useAuthStore } from "@/stores/auth-store";
import { useEditorStore } from "@/stores/editor-store";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { useHistoryStore } from "@/stores/history-store";
import { useExport } from "@/hooks/useExport";
import { useImport } from "@/hooks/useImport";
import { useScreenshot } from "@/hooks/useScreenshot";
import { addMissingHydrogens } from "@/lib/hydrogenate";
import { cn } from "@/lib/cn";
import type { RenderStyle } from "@/types/editor";

interface MenuItem {
  label: string;
  shortcut?: string;
  action?: () => void;
  separator?: boolean;
  checked?: boolean;
  disabled?: boolean;
}

interface MenuDef {
  label: string;
  items: MenuItem[];
}

interface MenuBarProps {
  onShowHelp?: () => void;
}

export function MenuBar({ onShowHelp }: MenuBarProps) {
  const [openMenu, setOpenMenu] = useState<string | null>(null);
  const menuRef = useRef<HTMLDivElement>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);

  const { exportMol, exportSdf, exportXyz } = useExport();
  const { importFile } = useImport();
  const { captureScreenshot } = useScreenshot();
  const logout = useAuthStore((s) => s.logout);
  const email = useAuthStore((s) => s.email);
  const toggleHydrogens = useEditorStore((s) => s.toggleHydrogens);
  const toggleLabels = useEditorStore((s) => s.toggleLabels);
  const showHydrogens = useEditorStore((s) => s.showHydrogens);
  const showLabels = useEditorStore((s) => s.showLabels);
  const renderStyle = useEditorStore((s) => s.renderStyle);
  const setRenderStyle = useEditorStore((s) => s.setRenderStyle);
  const clearSelection = useEditorStore((s) => s.clearSelection);
  const clearMeasurements = useEditorStore((s) => s.clearMeasurements);
  const canUndo = useHistoryStore((s) => s.canUndo);
  const canRedo = useHistoryStore((s) => s.canRedo);

  function handleUndo() {
    const entry = useWorkspaceStore.getState().getActive();
    if (!entry?.structure) return;
    const restored = useHistoryStore.getState().undo(entry.structure);
    if (restored) {
      useWorkspaceStore.getState().updateActiveStructure(restored);
    }
  }

  function handleRedo() {
    const entry = useWorkspaceStore.getState().getActive();
    if (!entry?.structure) return;
    const restored = useHistoryStore.getState().redo(entry.structure);
    if (restored) {
      useWorkspaceStore.getState().updateActiveStructure(restored);
    }
  }

  function handleImportClick() {
    fileInputRef.current?.click();
    setOpenMenu(null);
  }

  function handleFileInput(e: React.ChangeEvent<HTMLInputElement>) {
    const file = e.target.files?.[0];
    if (file) importFile(file);
    e.target.value = "";
  }

  const menus: MenuDef[] = [
    {
      label: "File",
      items: [
        { label: "Import File...", shortcut: "Ctrl+O", action: handleImportClick },
        { separator: true, label: "" },
        { label: "Export as .mol", action: exportMol },
        { label: "Export as .sdf", action: exportSdf },
        { label: "Export as .xyz", action: exportXyz },
        { label: "Export as .png", action: captureScreenshot },
        { separator: true, label: "" },
        { label: "Sign Out", action: logout },
      ],
    },
    {
      label: "Edit",
      items: [
        { label: "Undo", shortcut: "Ctrl+Z", action: handleUndo, disabled: !canUndo },
        { label: "Redo", shortcut: "Ctrl+Y", action: handleRedo, disabled: !canRedo },
        { separator: true, label: "" },
        {
          label: "Select All",
          shortcut: "Ctrl+A",
          action: () => {
            const entry = useWorkspaceStore.getState().getActive();
            if (entry?.structure) {
              useEditorStore.getState().selectAllAtoms(entry.structure.atoms.length);
            }
          },
        },
        {
          label: "Invert Selection",
          action: () => {
            const entry = useWorkspaceStore.getState().getActive();
            if (entry?.structure) {
              useEditorStore.getState().invertSelection(entry.structure.atoms.length);
            }
          },
        },
        { label: "Clear Selection", shortcut: "Esc", action: clearSelection },
        { separator: true, label: "" },
        { label: "Clear Measurements", action: clearMeasurements },
        { separator: true, label: "" },
        {
          label: "Add Hydrogens",
          shortcut: "Ctrl+H",
          action: () => {
            const entry = useWorkspaceStore.getState().getActive();
            if (entry?.structure) {
              useHistoryStore.getState().pushSnapshot(entry.structure, "Add hydrogens");
              const filled = addMissingHydrogens(entry.structure);
              useWorkspaceStore.getState().updateActiveStructure(filled);
            }
          },
        },
      ],
    },
    {
      label: "View",
      items: [
        { label: "Show Hydrogens", shortcut: "H", action: toggleHydrogens, checked: showHydrogens },
        { label: "Show Labels", shortcut: "L", action: toggleLabels, checked: showLabels },
        {
          label: "Distance Labels",
          shortcut: "D",
          action: () => useEditorStore.getState().toggleDistanceLabels(),
          checked: useEditorStore.getState().showDistanceLabels,
        },
        {
          label: "Grid",
          shortcut: "G",
          action: () => useEditorStore.getState().toggleGrid(),
          checked: useEditorStore.getState().showGrid,
        },
        {
          label: "Auto Rotate",
          action: () => useEditorStore.getState().toggleAutoRotate(),
          checked: useEditorStore.getState().autoRotate,
        },
        { separator: true, label: "" },
        {
          label: "Ball & Stick",
          action: () => setRenderStyle("ball-and-stick" as RenderStyle),
          checked: renderStyle === "ball-and-stick",
        },
        {
          label: "Stick",
          action: () => setRenderStyle("stick" as RenderStyle),
          checked: renderStyle === "stick",
        },
        {
          label: "Spacefill",
          action: () => setRenderStyle("spacefill" as RenderStyle),
          checked: renderStyle === "spacefill",
        },
      ],
    },
    {
      label: "Help",
      items: [
        { label: "Keyboard Shortcuts", shortcut: "?", action: onShowHelp },
      ],
    },
  ];

  useEffect(() => {
    function handleClickOutside(e: MouseEvent) {
      if (menuRef.current && !menuRef.current.contains(e.target as Node)) {
        setOpenMenu(null);
      }
    }
    document.addEventListener("mousedown", handleClickOutside);
    return () => document.removeEventListener("mousedown", handleClickOutside);
  }, []);

  return (
    <div
      ref={menuRef}
      className="flex items-center h-10 bg-bg-toolbar border-b border-border px-2 select-none"
    >
      <span className="text-xs font-bold text-accent mr-4 tracking-wide">
        MolBuilder Studio
      </span>

      {menus.map((menu) => (
        <div key={menu.label} className="relative">
          <button
            className={cn(
              "px-3 py-1 text-xs text-text-secondary hover:text-text-primary hover:bg-white/5 rounded-[var(--radius-sm)] transition-colors cursor-pointer",
              openMenu === menu.label && "bg-white/5 text-text-primary",
            )}
            onClick={() =>
              setOpenMenu(openMenu === menu.label ? null : menu.label)
            }
            onMouseEnter={() => openMenu && setOpenMenu(menu.label)}
          >
            {menu.label}
          </button>
          {openMenu === menu.label && (
            <div className="absolute top-full left-0 mt-0.5 min-w-[200px] bg-bg-card border border-border rounded-[var(--radius-sm)] py-1 z-50 shadow-lg">
              {menu.items.map((item, i) =>
                item.separator ? (
                  <div key={i} className="h-px bg-border my-1" />
                ) : (
                  <button
                    key={item.label}
                    className={cn(
                      "w-full flex items-center gap-2 px-3 py-1.5 text-xs transition-colors cursor-pointer",
                      item.disabled
                        ? "text-text-muted/40 cursor-not-allowed"
                        : "text-text-secondary hover:text-text-primary hover:bg-white/5",
                    )}
                    disabled={item.disabled}
                    onClick={() => {
                      if (!item.disabled) {
                        item.action?.();
                        setOpenMenu(null);
                      }
                    }}
                  >
                    <span className="w-3 text-center">
                      {item.checked != null && item.checked ? "\u2713" : ""}
                    </span>
                    <span className="flex-1 text-left">{item.label}</span>
                    {item.shortcut && (
                      <span className="text-text-muted text-[10px]">
                        {item.shortcut}
                      </span>
                    )}
                  </button>
                ),
              )}
            </div>
          )}
        </div>
      ))}

      <div className="flex-1" />

      {email && (
        <span className="text-[10px] text-text-muted font-mono mr-2">
          {email}
        </span>
      )}

      <input
        ref={fileInputRef}
        type="file"
        accept=".mol,.sdf,.sd,.xyz"
        className="hidden"
        onChange={handleFileInput}
      />
    </div>
  );
}
