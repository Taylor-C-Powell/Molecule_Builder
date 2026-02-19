import { useCallback } from "react";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { generateMolV2000 } from "@/lib/mol-export";
import { generateSdf } from "@/lib/sdf-export";
import { generateXyz } from "@/lib/xyz-export";

function downloadFile(content: string, filename: string, mimeType: string) {
  const blob = new Blob([content], { type: mimeType });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

export function useExport() {
  const exportMol = useCallback(() => {
    const entry = useWorkspaceStore.getState().getActive();
    if (!entry?.structure) return;
    const mol = generateMolV2000(entry.structure, entry.name);
    downloadFile(mol, `${entry.name || "molecule"}.mol`, "chemical/x-mdl-molfile");
  }, []);

  const exportSdf = useCallback(() => {
    const entry = useWorkspaceStore.getState().getActive();
    if (!entry?.structure) return;
    const sdf = generateSdf(entry.structure, entry.name, entry.properties);
    downloadFile(sdf, `${entry.name || "molecule"}.sdf`, "chemical/x-mdl-sdfile");
  }, []);

  const exportXyz = useCallback(() => {
    const entry = useWorkspaceStore.getState().getActive();
    if (!entry?.structure) return;
    const xyz = generateXyz(entry.structure, entry.name);
    downloadFile(xyz, `${entry.name || "molecule"}.xyz`, "chemical/x-xyz");
  }, []);

  return { exportMol, exportSdf, exportXyz };
}
