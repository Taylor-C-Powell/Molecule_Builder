import { API_BASE } from "@/lib/constants";
import type { RegisterRequest, RegisterResponse, TokenResponse } from "./types";

export async function register(email: string): Promise<RegisterResponse> {
  const res = await fetch(`${API_BASE}/v1/auth/register`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ email } satisfies RegisterRequest),
  });
  if (!res.ok) {
    const err = await res.json().catch(() => ({ detail: "Registration failed" }));
    throw new Error(err.detail);
  }
  return res.json();
}

export async function getToken(apiKey: string): Promise<TokenResponse> {
  const res = await fetch(`${API_BASE}/v1/auth/token`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ api_key: apiKey }),
  });
  if (!res.ok) {
    const err = await res.json().catch(() => ({ detail: "Token exchange failed" }));
    throw new Error(err.detail);
  }
  return res.json();
}
