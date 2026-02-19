import { PanelContainer } from "./PanelContainer";
import { Button } from "@/components/ui/Button";
import { useExport } from "@/hooks/useExport";
import { useScreenshot } from "@/hooks/useScreenshot";
import { useWorkspaceStore } from "@/stores/workspace-store";

export function ExportPanel() {
  const { exportMol, exportSdf, exportXyz } = useExport();
  const { captureScreenshot } = useScreenshot();
  const activeId = useWorkspaceStore((s) => s.activeId);
  const molecules = useWorkspaceStore((s) => s.molecules);
  const entry = activeId ? molecules.get(activeId) : undefined;
  const hasStructure = !!entry?.structure;

  return (
    <PanelContainer title="Export" defaultOpen={false}>
      <div className="space-y-1">
        <Button
          variant="secondary"
          size="sm"
          className="w-full justify-start text-xs"
          disabled={!hasStructure}
          onClick={exportMol}
        >
          Download .mol (V2000)
        </Button>
        <Button
          variant="secondary"
          size="sm"
          className="w-full justify-start text-xs"
          disabled={!hasStructure}
          onClick={exportSdf}
        >
          Download .sdf
        </Button>
        <Button
          variant="secondary"
          size="sm"
          className="w-full justify-start text-xs"
          disabled={!hasStructure}
          onClick={exportXyz}
        >
          Download .xyz
        </Button>
        <Button
          variant="secondary"
          size="sm"
          className="w-full justify-start text-xs"
          disabled={!hasStructure}
          onClick={captureScreenshot}
        >
          Screenshot .png
        </Button>
      </div>
    </PanelContainer>
  );
}
