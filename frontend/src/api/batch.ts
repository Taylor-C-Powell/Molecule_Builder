import type { ApiClient } from "./client";
import type {
  BatchSubmitRequest,
  BatchSubmitResponse,
  BatchStatusResponse,
  BatchListResponse,
} from "./types";

export function submitBatchApi(client: ApiClient, req: BatchSubmitRequest) {
  return client.post<BatchSubmitResponse>("/v1/batch/submit", req);
}

export function getBatchStatusApi(client: ApiClient, jobId: string) {
  return client.get<BatchStatusResponse>(`/v1/batch/${jobId}`);
}

export function listBatchJobsApi(
  client: ApiClient,
  params?: { page?: number; per_page?: number },
) {
  const query = new URLSearchParams();
  if (params?.page) query.set("page", String(params.page));
  if (params?.per_page) query.set("per_page", String(params.per_page));
  const qs = query.toString();
  return client.get<BatchListResponse>(`/v1/batch/${qs ? `?${qs}` : ""}`);
}

export function cancelBatchJobApi(client: ApiClient, jobId: string) {
  return client.delete<{ status: string; job_id: string }>(`/v1/batch/${jobId}`);
}
