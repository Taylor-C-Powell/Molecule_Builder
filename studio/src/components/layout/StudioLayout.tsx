import { useState, useCallback } from "react";
import { MenuBar } from "@/components/toolbar/MenuBar";
import { WorkspaceTabs } from "@/components/toolbar/WorkspaceTabs";
import { ToolPalette } from "@/components/toolbar/ToolPalette";
import { StatusBar } from "@/components/toolbar/StatusBar";
import { Viewport } from "@/components/viewport/Viewport";
import { ViewportContextMenu } from "@/components/viewport/ViewportContextMenu";
import { SmilesPanel } from "@/components/panels/SmilesPanel";
import { ImportPanel } from "@/components/panels/ImportPanel";
import { TemplatePanel } from "@/components/panels/TemplatePanel";
import { PropertiesPanel } from "@/components/panels/PropertiesPanel";
import { MeasurementsPanel } from "@/components/panels/MeasurementsPanel";
import { RetrosynthesisPanel } from "@/components/panels/RetrosynthesisPanel";
import { ProcessPanel } from "@/components/panels/ProcessPanel";
import { AtomsPanel } from "@/components/panels/AtomsPanel";
import { ExportPanel } from "@/components/panels/ExportPanel";
import { SettingsPanel } from "@/components/panels/SettingsPanel";
import { ShortcutsOverlay } from "@/components/panels/ShortcutsOverlay";
import { ToastContainer } from "@/components/ui/Toast";
import { MobileWarning } from "@/components/ui/MobileWarning";
import { useKeyboardShortcuts } from "@/hooks/useKeyboardShortcuts";
import { useMolecule } from "@/hooks/useMolecule";

export function StudioLayout() {
  const [showHelp, setShowHelp] = useState(false);
  const { parseMolecule } = useMolecule();
  const onPasteSmiles = useCallback(
    (smiles: string) => parseMolecule(smiles),
    [parseMolecule],
  );
  useKeyboardShortcuts({
    onShowHelp: () => setShowHelp(true),
    onPasteSmiles,
  });

  return (
    <div className="h-screen w-screen grid grid-rows-[40px_auto_1fr_28px] grid-cols-[48px_1fr_320px] bg-bg">
      {/* Row 1: Menu bar spans all columns */}
      <div className="col-span-3">
        <MenuBar onShowHelp={() => setShowHelp(true)} />
      </div>

      {/* Row 2: Workspace tabs (conditional, only if >1 molecule) */}
      <div className="col-span-3">
        <WorkspaceTabs />
      </div>

      {/* Row 3 Col 1: Tool palette */}
      <ToolPalette />

      {/* Row 3 Col 2: 3D Viewport */}
      <Viewport />

      {/* Row 3 Col 3: Right sidebar */}
      <div className="bg-bg-panel border-l border-border overflow-y-auto">
        <SmilesPanel />
        <ImportPanel />
        <TemplatePanel />
        <PropertiesPanel />
        <MeasurementsPanel />
        <RetrosynthesisPanel />
        <ProcessPanel />
        <AtomsPanel />
        <ExportPanel />
        <SettingsPanel />
      </div>

      {/* Row 4: Status bar spans all columns */}
      <div className="col-span-3">
        <StatusBar />
      </div>

      {/* Overlays */}
      <ShortcutsOverlay open={showHelp} onClose={() => setShowHelp(false)} />
      <ViewportContextMenu />
      <ToastContainer />
      <MobileWarning />
    </div>
  );
}
