import { useState, useCallback } from "react";
import { useApiClient } from "./useApiClient";
import { useAuthStore } from "@/stores/auth-store";
import { evaluateProcessApi } from "@/api/process";
import type { ProcessEvaluateRequest, ProcessEvaluateResponse } from "@/api/types";

export function useProcess() {
  const client = useApiClient();
  const getValidToken = useAuthStore((s) => s.getValidToken);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<ProcessEvaluateResponse | null>(null);

  const evaluate = useCallback(
    async (req: ProcessEvaluateRequest) => {
      setLoading(true);
      setError(null);
      try {
        const token = await getValidToken();
        if (!token) throw new Error("Not authenticated");
        const data = await evaluateProcessApi(client, req);
        setResult(data);
      } catch (err) {
        const msg = err instanceof Error ? err.message : "Process evaluation failed";
        setError(msg);
      } finally {
        setLoading(false);
      }
    },
    [client, getValidToken],
  );

  const clear = useCallback(() => {
    setResult(null);
    setError(null);
  }, []);

  return { evaluate, loading, error, result, clear };
}
