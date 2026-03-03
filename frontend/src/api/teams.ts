import type { ApiClient } from "./client";
import type {
  TeamCreateRequest,
  TeamUpdateRequest,
  TeamResponse,
  TeamListResponse,
  MemberAddRequest,
  MemberRoleUpdateRequest,
  MemberResponse,
  MemberListResponse,
  TeamLibrarySaveRequest,
  TeamLibraryUpdateRequest,
  TeamLibraryMoleculeResponse,
  TeamLibraryListResponse,
  TeamLibraryImportRequest,
  TeamLibraryImportResponse,
} from "./types";

export function createTeamApi(client: ApiClient, req: TeamCreateRequest) {
  return client.post<TeamResponse>("/v1/teams/", req);
}

export function listTeamsApi(client: ApiClient) {
  return client.get<TeamListResponse>("/v1/teams/");
}

export function getTeamApi(client: ApiClient, teamId: number) {
  return client.get<TeamResponse>(`/v1/teams/${teamId}`);
}

export function updateTeamApi(client: ApiClient, teamId: number, req: TeamUpdateRequest) {
  return client.patch<TeamResponse>(`/v1/teams/${teamId}`, req);
}

export function deleteTeamApi(client: ApiClient, teamId: number) {
  return client.delete<{ status: string }>(`/v1/teams/${teamId}`);
}

// ---- Members ----

export function listMembersApi(client: ApiClient, teamId: number) {
  return client.get<MemberListResponse>(`/v1/teams/${teamId}/members`);
}

export function addMemberApi(client: ApiClient, teamId: number, req: MemberAddRequest) {
  return client.post<MemberResponse>(`/v1/teams/${teamId}/members`, req);
}

export function updateMemberRoleApi(
  client: ApiClient,
  teamId: number,
  email: string,
  req: MemberRoleUpdateRequest,
) {
  return client.patch<MemberResponse>(
    `/v1/teams/${teamId}/members/${encodeURIComponent(email)}`,
    req,
  );
}

export function removeMemberApi(client: ApiClient, teamId: number, email: string) {
  return client.delete<{ status: string }>(
    `/v1/teams/${teamId}/members/${encodeURIComponent(email)}`,
  );
}

// ---- Team Library ----

export function saveTeamMoleculeApi(
  client: ApiClient,
  teamId: number,
  req: TeamLibrarySaveRequest,
) {
  return client.post<TeamLibraryMoleculeResponse>(`/v1/teams/${teamId}/library/`, req);
}

export function getTeamLibraryApi(
  client: ApiClient,
  teamId: number,
  params?: { tag?: string; search?: string; page?: number; per_page?: number },
) {
  const query = new URLSearchParams();
  if (params?.tag) query.set("tag", params.tag);
  if (params?.search) query.set("search", params.search);
  if (params?.page) query.set("page", String(params.page));
  if (params?.per_page) query.set("per_page", String(params.per_page));
  const qs = query.toString();
  return client.get<TeamLibraryListResponse>(
    `/v1/teams/${teamId}/library/${qs ? `?${qs}` : ""}`,
  );
}

export function updateTeamMoleculeApi(
  client: ApiClient,
  teamId: number,
  molId: number,
  req: TeamLibraryUpdateRequest,
) {
  return client.patch<TeamLibraryMoleculeResponse>(
    `/v1/teams/${teamId}/library/${molId}`,
    req,
  );
}

export function deleteTeamMoleculeApi(client: ApiClient, teamId: number, molId: number) {
  return client.delete<{ status: string; id: number }>(
    `/v1/teams/${teamId}/library/${molId}`,
  );
}

export function importTeamMoleculesApi(
  client: ApiClient,
  teamId: number,
  req: TeamLibraryImportRequest,
) {
  return client.post<TeamLibraryImportResponse>(
    `/v1/teams/${teamId}/library/import`,
    req,
  );
}
