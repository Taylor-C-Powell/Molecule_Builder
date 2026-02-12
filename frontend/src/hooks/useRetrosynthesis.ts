import { useState, useCallback } from "react";
import { useApiClient } from "./useApiClient";
import { planRetrosynthesisApi } from "@/api/retrosynthesis";
import { useAuthStore } from "@/stores/auth-store";
import type { RetroResponse } from "@/api/types";

export function useRetrosynthesis() {
  const client = useApiClient();
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<RetroResponse | null>(null);

  const plan = useCallback(
    async (smiles: string, maxDepth = 5, beamWidth = 5) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const res = await planRetrosynthesisApi(client, {
          smiles,
          max_depth: maxDepth,
          beam_width: beamWidth,
        });
        setResult(res);
      } catch (err) {
        setError(err instanceof Error ? err.message : "Retrosynthesis failed");
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  return { plan, loading, error, result };
}
