import { useState, useCallback } from "react";
import { useApiClient } from "./useApiClient";
import { useAuthStore } from "@/stores/auth-store";
import { planRetrosynthesisApi } from "@/api/retrosynthesis";
import type { RetroRequest, RetroResponse } from "@/api/types";

export function useRetrosynthesis() {
  const client = useApiClient();
  const getValidToken = useAuthStore((s) => s.getValidToken);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<RetroResponse | null>(null);

  const plan = useCallback(
    async (req: RetroRequest) => {
      setLoading(true);
      setError(null);
      try {
        const token = await getValidToken();
        if (!token) throw new Error("Not authenticated");
        const data = await planRetrosynthesisApi(client, req);
        setResult(data);
      } catch (err) {
        const msg = err instanceof Error ? err.message : "Retrosynthesis failed";
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

  return { plan, loading, error, result, clear };
}
