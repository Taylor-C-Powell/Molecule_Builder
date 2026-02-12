import { useState, useCallback } from "react";
import { useApiClient } from "./useApiClient";
import { useMoleculeStore } from "@/stores/molecule-store";
import { parseSmilesApi, getMolecule3dApi, getMoleculePropertiesApi } from "@/api/molecule";
import { useAuthStore } from "@/stores/auth-store";

export function useMolecule() {
  const client = useApiClient();
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Subscribe to store changes via selectors so React re-renders
  const current = useMoleculeStore((s) => s.current);
  const cache = useMoleculeStore((s) => s.cache);
  const currentEntry = current ? cache.get(current) : undefined;

  const parseMolecule = useCallback(
    async (smiles: string) => {
      setLoading(true);
      setError(null);

      const store = useMoleculeStore.getState();

      try {
        await useAuthStore.getState().getValidToken();

        // Check cache first
        const cached = store.getEntry(smiles);
        if (cached?.molecule && cached.structure && cached.properties) {
          store.setCurrent(smiles);
          setLoading(false);
          return;
        }

        // Parse SMILES
        const mol = cached?.molecule ?? await parseSmilesApi(client, { smiles });
        store.setMolecule(smiles, mol);
        store.setCurrent(smiles);

        // Fetch 3D and properties in parallel
        const [structure, properties] = await Promise.all([
          cached?.structure ?? getMolecule3dApi(client, mol.id),
          cached?.properties ?? getMoleculePropertiesApi(client, mol.id),
        ]);

        store.setStructure(smiles, structure);
        store.setProperties(smiles, properties);
      } catch (err) {
        const msg = err instanceof Error ? err.message : "Failed to parse molecule";
        setError(msg);
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  return {
    parseMolecule,
    loading,
    error,
    currentEntry,
  };
}
