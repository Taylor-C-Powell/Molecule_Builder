import type { ApiClient } from "./client";
import type {
  LibrarySaveRequest,
  LibraryUpdateRequest,
  LibraryMoleculeResponse,
  LibraryListResponse,
  LibraryImportRequest,
  LibraryImportResponse,
} from "./types";

export function saveToLibraryApi(client: ApiClient, req: LibrarySaveRequest) {
  return client.post<LibraryMoleculeResponse>("/v1/library/", req);
}

export function getLibraryApi(
  client: ApiClient,
  params?: { tag?: string; search?: string; page?: number; per_page?: number },
) {
  const query = new URLSearchParams();
  if (params?.tag) query.set("tag", params.tag);
  if (params?.search) query.set("search", params.search);
  if (params?.page) query.set("page", String(params.page));
  if (params?.per_page) query.set("per_page", String(params.per_page));
  const qs = query.toString();
  return client.get<LibraryListResponse>(`/v1/library/${qs ? `?${qs}` : ""}`);
}

export function getLibraryMoleculeApi(client: ApiClient, molId: number) {
  return client.get<LibraryMoleculeResponse>(`/v1/library/${molId}`);
}

export function updateLibraryMoleculeApi(
  client: ApiClient,
  molId: number,
  req: LibraryUpdateRequest,
) {
  return client.put<LibraryMoleculeResponse>(`/v1/library/${molId}`, req);
}

export function deleteLibraryMoleculeApi(client: ApiClient, molId: number) {
  return client.delete<{ status: string; id: number }>(`/v1/library/${molId}`);
}

export function importToLibraryApi(client: ApiClient, req: LibraryImportRequest) {
  return client.post<LibraryImportResponse>("/v1/library/import", req);
}
