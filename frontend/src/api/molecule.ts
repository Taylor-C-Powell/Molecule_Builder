import type { ApiClient } from "./client";
import type {
  ADMETResponse,
  ParseSmilesRequest,
  MoleculeResponse,
  Molecule3DResponse,
  MoleculePropertiesResponse,
  SolubilityResponse,
  FileImportResponse,
} from "./types";

export function parseSmilesApi(client: ApiClient, req: ParseSmilesRequest) {
  return client.post<MoleculeResponse>("/v1/molecule/from-smiles", req);
}

export function getMolecule3dApi(client: ApiClient, id: string) {
  return client.get<Molecule3DResponse>(`/v1/molecule/${id}/3d`);
}

export function getMoleculePropertiesApi(client: ApiClient, id: string) {
  return client.get<MoleculePropertiesResponse>(`/v1/molecule/${id}/properties`);
}

export function getADMETApi(client: ApiClient, id: string) {
  return client.get<ADMETResponse>(`/v1/molecule/${id}/admet`);
}

export function getSolubilityApi(client: ApiClient, id: string) {
  return client.get<SolubilityResponse>(`/v1/molecule/${id}/solubility`);
}

export function importFileApi(client: ApiClient, file: File) {
  return client.postFile<FileImportResponse>("/v1/molecule/import-file", file);
}

export function exportMoleculeApi(client: ApiClient, molId: string, format: string) {
  return client.downloadBlob(`/v1/molecule/${molId}/export/${format}`);
}
