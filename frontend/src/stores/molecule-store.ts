import { create } from "zustand";
import type {
  MoleculeResponse,
  Molecule3DResponse,
  MoleculePropertiesResponse,
} from "@/api/types";

interface MoleculeEntry {
  smiles: string;
  molecule: MoleculeResponse;
  structure?: Molecule3DResponse;
  properties?: MoleculePropertiesResponse;
}

interface MoleculeState {
  /** SMILES -> cached entry */
  cache: Map<string, MoleculeEntry>;
  /** Currently active SMILES */
  current: string | null;

  setCurrent: (smiles: string | null) => void;
  setMolecule: (smiles: string, mol: MoleculeResponse) => void;
  setStructure: (smiles: string, s: Molecule3DResponse) => void;
  setProperties: (smiles: string, p: MoleculePropertiesResponse) => void;
  getEntry: (smiles: string) => MoleculeEntry | undefined;
  getCurrentEntry: () => MoleculeEntry | undefined;
}

export const useMoleculeStore = create<MoleculeState>((set, get) => ({
  cache: new Map(),
  current: null,

  setCurrent(smiles) {
    set({ current: smiles });
  },

  setMolecule(smiles, mol) {
    const cache = new Map(get().cache);
    const existing = cache.get(smiles);
    cache.set(smiles, { ...existing, smiles, molecule: mol });
    set({ cache });
  },

  setStructure(smiles, s) {
    const cache = new Map(get().cache);
    const existing = cache.get(smiles);
    if (existing) {
      cache.set(smiles, { ...existing, structure: s });
      set({ cache });
    }
  },

  setProperties(smiles, p) {
    const cache = new Map(get().cache);
    const existing = cache.get(smiles);
    if (existing) {
      cache.set(smiles, { ...existing, properties: p });
      set({ cache });
    }
  },

  getEntry(smiles) {
    return get().cache.get(smiles);
  },

  getCurrentEntry() {
    const { current, cache } = get();
    if (!current) return undefined;
    return cache.get(current);
  },
}));
