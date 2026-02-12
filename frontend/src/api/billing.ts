import type { ApiClient } from "./client";
import type {
  CheckoutRequest,
  CheckoutResponse,
  PortalResponse,
  BillingStatusResponse,
} from "./types";

export function createCheckoutApi(client: ApiClient, req: CheckoutRequest) {
  return client.post<CheckoutResponse>("/v1/billing/checkout", req);
}

export function getBillingStatusApi(client: ApiClient) {
  return client.get<BillingStatusResponse>("/v1/billing/status");
}

export function createPortalSessionApi(client: ApiClient) {
  return client.post<PortalResponse>("/v1/billing/portal");
}
