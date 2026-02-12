import type { ApiClient } from "./client";
import type { ProcessEvaluateRequest, ProcessEvaluateResponse } from "./types";

export function evaluateProcessApi(
  client: ApiClient,
  req: ProcessEvaluateRequest,
) {
  return client.post<ProcessEvaluateResponse>("/v1/process/evaluate", req);
}
