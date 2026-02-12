import { useMemo } from "react";
import { ApiClient } from "@/api/client";
import { API_BASE } from "@/lib/constants";
import { useAuthStore } from "@/stores/auth-store";

export function useApiClient(): ApiClient {
  return useMemo(() => {
    return new ApiClient(API_BASE, () => {
      const { jwt } = useAuthStore.getState();
      return jwt;
    });
  }, []);
}
