import { useState, useCallback, useRef, useEffect } from "react";
import { useApiClient } from "./useApiClient";
import { useAuthStore } from "@/stores/auth-store";
import {
  submitBatchApi,
  getBatchStatusApi,
  listBatchJobsApi,
  cancelBatchJobApi,
} from "@/api/batch";
import type { BatchStatusResponse, BatchListResponse } from "@/api/types";

export function useBatch() {
  const client = useApiClient();
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [jobStatus, setJobStatus] = useState<BatchStatusResponse | null>(null);
  const [jobs, setJobs] = useState<BatchListResponse | null>(null);
  const pollRef = useRef<ReturnType<typeof setInterval> | null>(null);

  const stopPolling = useCallback(() => {
    if (pollRef.current) {
      clearInterval(pollRef.current);
      pollRef.current = null;
    }
  }, []);

  useEffect(() => stopPolling, [stopPolling]);

  const submit = useCallback(
    async (
      smilesList: string[],
      jobType: "properties" | "retrosynthesis" | "conditions" | "evaluate",
      params?: Record<string, unknown>,
    ) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const res = await submitBatchApi(client, {
          smiles_list: smilesList,
          job_type: jobType,
          params,
        });
        return res;
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to submit batch job");
        return null;
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  const pollStatus = useCallback(
    (jobId: string, intervalMs = 5000) => {
      stopPolling();
      const poll = async () => {
        try {
          await useAuthStore.getState().getValidToken();
          const status = await getBatchStatusApi(client, jobId);
          setJobStatus(status);
          if (status.status === "completed" || status.status === "failed" || status.status === "cancelled") {
            stopPolling();
          }
        } catch {
          stopPolling();
        }
      };
      void poll();
      pollRef.current = setInterval(poll, intervalMs);
    },
    [client, stopPolling],
  );

  const fetchList = useCallback(
    async (params?: { page?: number; per_page?: number }) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const res = await listBatchJobsApi(client, params);
        setJobs(res);
        return res;
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to list batch jobs");
        return null;
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  const cancel = useCallback(
    async (jobId: string) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        await cancelBatchJobApi(client, jobId);
        stopPolling();
        return true;
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to cancel job");
        return false;
      } finally {
        setLoading(false);
      }
    },
    [client, stopPolling],
  );

  return { submit, pollStatus, stopPolling, fetchList, cancel, loading, error, jobStatus, jobs };
}
