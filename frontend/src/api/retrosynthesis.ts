import type { ApiClient } from "./client";
import type { RetroRequest, RetroResponse } from "./types";

export function planRetrosynthesisApi(client: ApiClient, req: RetroRequest) {
  return client.post<RetroResponse>("/v1/retrosynthesis/plan", req);
}
