import type { ApiClient } from "./client";
import type {
  ParseSmilesRequest,
  MoleculeResponse,
  Molecule3DResponse,
  MoleculePropertiesResponse,
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
