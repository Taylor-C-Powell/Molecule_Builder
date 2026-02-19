import { useCallback, useState } from "react";
import { useApiClient } from "./useApiClient";
import { parseSmilesApi, getMolecule3dApi, getMoleculePropertiesApi } from "@/api/molecule";
import { useAuthStore } from "@/stores/auth-store";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { useEditorStore } from "@/stores/editor-store";

export function useMolecule() {
  const client = useApiClient();
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const parseMolecule = useCallback(
    async (smiles: string, name?: string) => {
      setLoading(true);
      setError(null);

      try {
        await useAuthStore.getState().getValidToken();

        const mol = await parseSmilesApi(client, { smiles, name });
        useWorkspaceStore.getState().addMolecule(smiles, mol);

        const [structure, properties] = await Promise.all([
          getMolecule3dApi(client, mol.id),
          getMoleculePropertiesApi(client, mol.id),
        ]);

        useWorkspaceStore.getState().setStructure(mol.id, structure);
        useWorkspaceStore.getState().setProperties(mol.id, properties);
        useEditorStore.getState().clearSelection();
      } catch (err) {
        const msg = err instanceof Error ? err.message : "Failed to parse molecule";
        setError(msg);
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  return { parseMolecule, loading, error };
}
