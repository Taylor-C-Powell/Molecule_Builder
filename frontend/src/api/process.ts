import type { ApiClient } from "./client";
import type { ProcessEvaluateRequest, ProcessEvaluateResponse } from "./types";

export function evaluateProcessApi(
  client: ApiClient,
  req: ProcessEvaluateRequest,
) {
  return client.post<ProcessEvaluateResponse>("/v1/process/evaluate", req);
}

export function downloadProcessPdf(
  client: ApiClient,
  smiles: string,
  scaleKg: number = 1.0,
) {
  const params = new URLSearchParams({ smiles, scale_kg: String(scaleKg) });
  return client.downloadBlob(`/v1/reports/process-pdf?${params}`, "POST");
}
