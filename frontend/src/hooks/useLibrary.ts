import { useState, useCallback } from "react";
import { useApiClient } from "./useApiClient";
import { useAuthStore } from "@/stores/auth-store";
import {
  saveToLibraryApi,
  getLibraryApi,
  getLibraryMoleculeApi,
  updateLibraryMoleculeApi,
  deleteLibraryMoleculeApi,
  importToLibraryApi,
} from "@/api/library";
import type {
  LibraryMoleculeResponse,
  LibraryListResponse,
  LibraryImportResponse,
} from "@/api/types";

export function useLibrary() {
  const client = useApiClient();
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [list, setList] = useState<LibraryListResponse | null>(null);
  const [selected, setSelected] = useState<LibraryMoleculeResponse | null>(null);

  const save = useCallback(
    async (smiles: string, name?: string, tags?: string[], notes?: string) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const mol = await saveToLibraryApi(client, { smiles, name, tags, notes });
        setSelected(mol);
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

  const fetchList = useCallback(
    async (params?: { tag?: string; search?: string; page?: number; per_page?: number }) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const res = await getLibraryApi(client, params);
        setList(res);
        return res;
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to fetch library");
        return null;
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  const get = useCallback(
    async (molId: number) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const mol = await getLibraryMoleculeApi(client, molId);
        setSelected(mol);
        return mol;
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to get molecule");
        return null;
      } finally {
        setLoading(false);
      }
    },
    [client],
  );

  const update = useCallback(
    async (molId: number, data: { name?: string; tags?: string[]; notes?: string }) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const mol = await updateLibraryMoleculeApi(client, molId, data);
        setSelected(mol);
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

  const remove = useCallback(
    async (molId: number) => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        await deleteLibraryMoleculeApi(client, molId);
        setSelected(null);
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

  const importSmiles = useCallback(
    async (smilesList: string[], tag?: string): Promise<LibraryImportResponse | null> => {
      setLoading(true);
      setError(null);
      try {
        await useAuthStore.getState().getValidToken();
        const res = await importToLibraryApi(client, { smiles_list: smilesList, tag });
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

  return { save, fetchList, get, update, remove, importSmiles, loading, error, list, selected };
}
