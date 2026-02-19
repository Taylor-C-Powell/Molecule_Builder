import { create } from "zustand";
import type { Molecule3DResponse } from "@/api/types";

interface HistoryEntry {
  structure: Molecule3DResponse;
  label: string;
}

interface HistoryState {
  past: HistoryEntry[];
  future: HistoryEntry[];
  canUndo: boolean;
  canRedo: boolean;

  pushSnapshot: (structure: Molecule3DResponse, label: string) => void;
  undo: (currentStructure: Molecule3DResponse) => Molecule3DResponse | null;
  redo: (currentStructure: Molecule3DResponse) => Molecule3DResponse | null;
  clear: () => void;
}

const MAX_HISTORY = 50;

function deepCloneStructure(s: Molecule3DResponse): Molecule3DResponse {
  return {
    id: s.id,
    atoms: s.atoms.map((a) => ({
      ...a,
      position: [...a.position] as [number, number, number],
    })),
    bonds: s.bonds.map((b) => ({ ...b })),
  };
}

export const useHistoryStore = create<HistoryState>((set, get) => ({
  past: [],
  future: [],
  canUndo: false,
  canRedo: false,

  pushSnapshot(structure, label) {
    const clone = deepCloneStructure(structure);
    const past = [...get().past, { structure: clone, label }];
    if (past.length > MAX_HISTORY) past.shift();
    set({ past, future: [], canUndo: true, canRedo: false });
  },

  undo(currentStructure) {
    const { past } = get();
    if (past.length === 0) return null;
    const entry = past[past.length - 1]!;
    const newPast = past.slice(0, -1);
    const future = [
      ...get().future,
      { structure: deepCloneStructure(currentStructure), label: "redo" },
    ];
    set({
      past: newPast,
      future,
      canUndo: newPast.length > 0,
      canRedo: true,
    });
    return deepCloneStructure(entry.structure);
  },

  redo(currentStructure) {
    const { future } = get();
    if (future.length === 0) return null;
    const entry = future[future.length - 1]!;
    const newFuture = future.slice(0, -1);
    const past = [
      ...get().past,
      { structure: deepCloneStructure(currentStructure), label: "undo" },
    ];
    set({
      past,
      future: newFuture,
      canUndo: true,
      canRedo: newFuture.length > 0,
    });
    return deepCloneStructure(entry.structure);
  },

  clear() {
    set({ past: [], future: [], canUndo: false, canRedo: false });
  },
}));
