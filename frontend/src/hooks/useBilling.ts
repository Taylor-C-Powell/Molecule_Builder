import { useState, useCallback } from "react";
import { useApiClient } from "./useApiClient";
import { getBillingStatusApi, createCheckoutApi, createPortalSessionApi } from "@/api/billing";
import { useAuthStore } from "@/stores/auth-store";
import type { BillingStatusResponse } from "@/api/types";

export function useBilling() {
  const client = useApiClient();
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [status, setStatus] = useState<BillingStatusResponse | null>(null);

  const fetchStatus = useCallback(async () => {
    setLoading(true);
    setError(null);
    try {
      await useAuthStore.getState().getValidToken();
      const res = await getBillingStatusApi(client);
      setStatus(res);
    } catch (err) {
      setError(err instanceof Error ? err.message : "Failed to load billing status");
    } finally {
      setLoading(false);
    }
  }, [client]);

  const checkout = useCallback(
    async (plan: "pro_monthly" | "pro_yearly") => {
      await useAuthStore.getState().getValidToken();
      const res = await createCheckoutApi(client, {
        plan,
        success_url: `${window.location.origin}/account?checkout=success`,
        cancel_url: `${window.location.origin}/account?checkout=cancel`,
      });
      window.location.href = res.checkout_url;
    },
    [client],
  );

  const openPortal = useCallback(async () => {
    await useAuthStore.getState().getValidToken();
    const res = await createPortalSessionApi(client);
    window.location.href = res.portal_url;
  }, [client]);

  return { fetchStatus, checkout, openPortal, status, loading, error };
}
