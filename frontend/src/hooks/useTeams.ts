import { useState, useCallback } from "react";
import { useApiClient } from "./useApiClient";
import { useAuthStore } from "@/stores/auth-store";
import {
  createTeamApi,
  listTeamsApi,
  getTeamApi,
  updateTeamApi,
  deleteTeamApi,
  listMembersApi,
  addMemberApi,
  updateMemberRoleApi,
  removeMemberApi,
  saveTeamMoleculeApi,
  getTeamLibraryApi,
  updateTeamMoleculeApi,
  deleteTeamMoleculeApi,
  importTeamMoleculesApi,
} from "@/api/teams";
import type {
  TeamResponse,
  MemberResponse,
  TeamLibraryListResponse,
  TeamLibraryImportResponse,
} from "@/api/types";

export function useTeams() {
  const client = useApiClient();
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [teams, setTeams] = useState<TeamResponse[]>([]);
  const [selectedTeam, setSelectedTeam] = useState<TeamResponse | null>(null);
  const [members, setMembers] = useState<MemberResponse[]>([]);
  const [library, setLibrary] = useState<TeamLibraryListResponse | null>(null);

  // ---- Teams CRUD ----

  const fetchTeams = useCallback(async () => {
    setLoading(true);
    setError(null);
    try {
      await useAuthStore.getState().getValidToken();
      const res = await listTeamsApi(client);
      setTeams(res.teams);
      return res.teams;
    } catch (err) {
      setError(err instanceof Error ? err.message : "Failed to fetch teams");
      return null;
    } finally {
      setLoading(false);
    }
  }, [client]);

  const createTeam = useCallback(
    async (name: string, slug: string) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const team = await createTeamApi(client, { name, slug });
        setSelectedTeam(team);
        return team;
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to create team");
        return null;
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  const fetchTeam = useCallback(
    async (teamId: number) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const team = await getTeamApi(client, teamId);
        setSelectedTeam(team);
        return team;
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to fetch team");
        return null;
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  const updateTeam = useCallback(
    async (teamId: number, name: string) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const team = await updateTeamApi(client, teamId, { name });
        setSelectedTeam(team);
        return team;
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to update team");
        return null;
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  const deleteTeam = useCallback(
    async (teamId: number) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        await deleteTeamApi(client, teamId);
        setSelectedTeam(null);
        return true;
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to delete team");
        return false;
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  // ---- Members ----

  const fetchMembers = useCallback(
    async (teamId: number) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const res = await listMembersApi(client, teamId);
        setMembers(res.members);
        return res.members;
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to fetch members");
        return null;
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  const addMember = useCallback(
    async (teamId: number, email: string, role?: string) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const member = await addMemberApi(client, teamId, { email, role });
        return member;
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to add member");
        return null;
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  const updateMemberRole = useCallback(
    async (teamId: number, email: string, role: string) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const member = await updateMemberRoleApi(client, teamId, email, { role });
        return member;
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to update member role");
        return null;
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  const removeMember = useCallback(
    async (teamId: number, email: string) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        await removeMemberApi(client, teamId, email);
        return true;
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to remove member");
        return false;
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  // ---- Team Library ----

  const saveTeamMolecule = useCallback(
    async (
      teamId: number,
      smiles: string,
      name?: string,
      tags?: string[],
      notes?: string,
    ) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const mol = await saveTeamMoleculeApi(client, teamId, {
          smiles,
          name,
          tags,
          notes,
        });
        return mol;
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to save molecule");
        return null;
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  const fetchTeamLibrary = useCallback(
    async (
      teamId: number,
      params?: { tag?: string; search?: string; page?: number; per_page?: number },
    ) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const res = await getTeamLibraryApi(client, teamId, params);
        setLibrary(res);
        return res;
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to fetch team library");
        return null;
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  const updateTeamMolecule = useCallback(
    async (
      teamId: number,
      molId: number,
      data: { name?: string; tags?: string[]; notes?: string },
    ) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const mol = await updateTeamMoleculeApi(client, teamId, molId, data);
        return mol;
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to update molecule");
        return null;
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  const deleteTeamMolecule = useCallback(
    async (teamId: number, molId: number) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        await deleteTeamMoleculeApi(client, teamId, molId);
        return true;
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to delete molecule");
        return false;
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  const importTeamMolecules = useCallback(
    async (
      teamId: number,
      smilesList: string[],
      tag?: string,
    ): Promise<TeamLibraryImportResponse | null> => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const res = await importTeamMoleculesApi(client, teamId, {
          smiles_list: smilesList,
          tag,
        });
        return res;
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to import molecules");
        return null;
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  return {
    loading,
    error,
    teams,
    selectedTeam,
    members,
    library,
    fetchTeams,
    createTeam,
    fetchTeam,
    updateTeam,
    deleteTeam,
    fetchMembers,
    addMember,
    updateMemberRole,
    removeMember,
    saveTeamMolecule,
    fetchTeamLibrary,
    updateTeamMolecule,
    deleteTeamMolecule,
    importTeamMolecules,
  };
}
