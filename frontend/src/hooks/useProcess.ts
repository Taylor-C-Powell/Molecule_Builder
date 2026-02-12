import { useState, useCallback } from "react";
import { useApiClient } from "./useApiClient";
import { evaluateProcessApi } from "@/api/process";
import { useAuthStore } from "@/stores/auth-store";
import type { ProcessEvaluateResponse } from "@/api/types";

export function useProcess() {
  const client = useApiClient();
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<ProcessEvaluateResponse | null>(null);

  const evaluate = useCallback(
    async (smiles: string, scaleKg = 1.0, maxDepth = 5, beamWidth = 5) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const res = await evaluateProcessApi(client, {
          smiles,
          scale_kg: scaleKg,
          max_depth: maxDepth,
          beam_width: beamWidth,
        });
        setResult(res);
      } catch (err) {
        setError(err instanceof Error ? err.message : "Process evaluation failed");
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  return { evaluate, loading, error, result };
}
