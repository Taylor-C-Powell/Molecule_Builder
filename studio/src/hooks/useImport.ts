import { useCallback } from "react";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { useToastStore } from "@/stores/toast-store";
import { parseMolV2000 } from "@/lib/mol-import";
import { parseSdf } from "@/lib/sdf-import";
import { parseXyz } from "@/lib/xyz-import";
import type { MoleculeResponse } from "@/api/types";

type ImportFormat = "mol" | "sdf" | "xyz";

function detectFormat(filename: string): ImportFormat | null {
  const ext = filename.split(".").pop()?.toLowerCase();
  if (ext === "mol") return "mol";
  if (ext === "sdf" || ext === "sd") return "sdf";
  if (ext === "xyz") return "xyz";
  return null;
}

export function useImport() {
  const addToast = useToastStore((s) => s.addToast);

  const importFile = useCallback(
    async (file: File) => {
      const format = detectFormat(file.name);
      if (!format) {
        addToast("Unsupported file format. Use .mol, .sdf, or .xyz", "error");
        return;
      }

      try {
        const content = await file.text();
        let result: { structure: { id: string; atoms: { index: number; symbol: string; position: [number, number, number]; hybridization: string | null; formal_charge: number }[]; bonds: { atom_i: number; atom_j: number; order: number; rotatable: boolean }[] }; name: string };

        switch (format) {
          case "mol":
            result = parseMolV2000(content);
            break;
          case "sdf":
            result = parseSdf(content);
            break;
          case "xyz":
            result = parseXyz(content);
            break;
        }

        // Create a MoleculeResponse wrapper
        const molResponse: MoleculeResponse = {
          id: result.structure.id,
          name: result.name,
          smiles: "",
          num_atoms: result.structure.atoms.length,
          num_bonds: result.structure.bonds.length,
        };

        const ws = useWorkspaceStore.getState();
        ws.addMolecule(result.name, molResponse);
        ws.setStructure(result.structure.id, result.structure);

        addToast(
          `Imported ${result.name} (${result.structure.atoms.length} atoms, ${result.structure.bonds.length} bonds)`,
          "success",
        );
      } catch (err) {
        const msg = err instanceof Error ? err.message : "Failed to parse file";
        addToast(msg, "error");
      }
    },
    [addToast],
  );

  return { importFile };
}
