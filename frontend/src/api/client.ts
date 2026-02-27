import { API_BASE } from "@/lib/constants";
import type { ApiError } from "./types";

export class ApiClient {
  private baseUrl: string;
  private getToken: () => string | null;

  constructor(baseUrl: string, getToken: () => string | null) {
    this.baseUrl = baseUrl;
    this.getToken = getToken;
  }

  private async request<T>(
    method: string,
    path: string,
    body?: unknown,
    opts?: { apiKey?: string },
  ): Promise<T> {
    const headers: Record<string, string> = {
      "Content-Type": "application/json",
    };

    if (opts?.apiKey) {
      headers["X-API-Key"] = opts.apiKey;
    } else {
      const token = this.getToken();
      if (token) {
        headers["Authorization"] = `Bearer ${token}`;
      }
    }

    const res = await fetch(`${this.baseUrl}${path}`, {
      method,
      headers,
      body: body ? JSON.stringify(body) : undefined,
    });

    if (!res.ok) {
      let detail = `Request failed: ${res.status}`;
      try {
        const err: ApiError = await res.json();
        detail = err.detail;
      } catch {
        // non-JSON error response
      }
      throw new ApiRequestError(detail, res.status);
    }

    return res.json();
  }

  get<T>(path: string) {
    return this.request<T>("GET", path);
  }

  post<T>(path: string, body?: unknown, opts?: { apiKey?: string }) {
    return this.request<T>("POST", path, body, opts);
  }

  put<T>(path: string, body?: unknown) {
    return this.request<T>("PUT", path, body);
  }

  delete<T>(path: string) {
    return this.request<T>("DELETE", path);
  }

  async postFile<T>(path: string, file: File): Promise<T> {
    const headers: Record<string, string> = {};
    const token = this.getToken();
    if (token) {
      headers["Authorization"] = `Bearer ${token}`;
    }
    const form = new FormData();
    form.append("file", file);

    const res = await fetch(`${this.baseUrl}${path}`, {
      method: "POST",
      headers,
      body: form,
    });

    if (!res.ok) {
      let detail = `Request failed: ${res.status}`;
      try {
        const err: ApiError = await res.json();
        detail = err.detail;
      } catch {
        // non-JSON error
      }
      throw new ApiRequestError(detail, res.status);
    }

    return res.json();
  }

  async downloadBlob(path: string, method: string = "GET"): Promise<Blob> {
    const headers: Record<string, string> = {};
    const token = this.getToken();
    if (token) {
      headers["Authorization"] = `Bearer ${token}`;
    }

    const res = await fetch(`${this.baseUrl}${path}`, {
      method,
      headers,
    });

    if (!res.ok) {
      let detail = `Request failed: ${res.status}`;
      try {
        const err: ApiError = await res.json();
        detail = err.detail;
      } catch {
        // non-JSON error
      }
      throw new ApiRequestError(detail, res.status);
    }

    return res.blob();
  }
}

export class ApiRequestError extends Error {
  status: number;
  constructor(message: string, status: number) {
    super(message);
    this.name = "ApiRequestError";
    this.status = status;
  }
}

let _client: ApiClient | null = null;

export function getClient(getToken: () => string | null): ApiClient {
  if (!_client) {
    _client = new ApiClient(API_BASE, getToken);
  }
  return _client;
}
