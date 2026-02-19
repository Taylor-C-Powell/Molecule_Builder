import { create } from "zustand";
import type {
  MoleculeResponse,
  Molecule3DResponse,
  MoleculePropertiesResponse,
  AtomResponse,
  BondResponse,
} from "@/api/types";
import type { PlaceableElement, BondOrder } from "@/types/editor";
import type { FragmentDef } from "@/lib/fragments";

export interface MoleculeEntry {
  id: string;
  smiles: string;
  name: string;
  molecule: MoleculeResponse;
  structure?: Molecule3DResponse;
  properties?: MoleculePropertiesResponse;
}

interface WorkspaceState {
  molecules: Map<string, MoleculeEntry>;
  activeId: string | null;

  setActive: (id: string) => void;
  addMolecule: (smiles: string, response: MoleculeResponse) => void;
  setStructure: (id: string, structure: Molecule3DResponse) => void;
  setProperties: (id: string, props: MoleculePropertiesResponse) => void;
  removeMolecule: (id: string) => void;
  getActive: () => MoleculeEntry | null;

  // Phase 2: local editing
  addAtom: (element: PlaceableElement, position: [number, number, number]) => number | null;
  addBond: (atomI: number, atomJ: number, order: BondOrder) => boolean;
  removeAtoms: (indices: Set<number>) => void;
  removeBonds: (keys: Set<string>) => void;
  updateActiveStructure: (structure: Molecule3DResponse) => void;

  // Phase 5: advanced editing
  replaceAtomElement: (index: number, element: PlaceableElement) => void;
  cycleBondOrder: (atomI: number, atomJ: number) => void;

  // Phase 7: atom dragging + fragments
  moveAtom: (index: number, position: [number, number, number]) => void;
  insertFragment: (parentAtomIndex: number, fragment: FragmentDef) => void;
}

export const useWorkspaceStore = create<WorkspaceState>((set, get) => ({
  molecules: new Map(),
  activeId: null,

  setActive(id) {
    set({ activeId: id });
  },

  addMolecule(smiles, response) {
    const molecules = new Map(get().molecules);
    molecules.set(response.id, {
      id: response.id,
      smiles,
      name: response.name,
      molecule: response,
    });
    set({ molecules, activeId: response.id });
  },

  setStructure(id, structure) {
    const molecules = new Map(get().molecules);
    const entry = molecules.get(id);
    if (entry) {
      molecules.set(id, { ...entry, structure });
      set({ molecules });
    }
  },

  setProperties(id, props) {
    const molecules = new Map(get().molecules);
    const entry = molecules.get(id);
    if (entry) {
      molecules.set(id, { ...entry, properties: props });
      set({ molecules });
    }
  },

  removeMolecule(id) {
    const molecules = new Map(get().molecules);
    molecules.delete(id);
    const activeId = get().activeId === id ? null : get().activeId;
    set({ molecules, activeId });
  },

  getActive() {
    const { activeId, molecules } = get();
    if (!activeId) return null;
    return molecules.get(activeId) ?? null;
  },

  addAtom(element, position) {
    const entry = get().getActive();
    if (!entry?.structure) return null;

    const structure = entry.structure;
    const newIndex = structure.atoms.length;
    const newAtom: AtomResponse = {
      index: newIndex,
      symbol: element,
      position: [...position],
      hybridization: null,
      formal_charge: 0,
    };

    const newStructure: Molecule3DResponse = {
      ...structure,
      atoms: [...structure.atoms, newAtom],
    };

    const molecules = new Map(get().molecules);
    molecules.set(entry.id, {
      ...entry,
      structure: newStructure,
      molecule: {
        ...entry.molecule,
        num_atoms: newStructure.atoms.length,
        num_bonds: newStructure.bonds.length,
      },
    });
    set({ molecules });
    return newIndex;
  },

  addBond(atomI, atomJ, order) {
    const entry = get().getActive();
    if (!entry?.structure) return false;

    const structure = entry.structure;
    if (atomI === atomJ) return false;
    if (atomI >= structure.atoms.length || atomJ >= structure.atoms.length) return false;

    // Check if bond already exists
    const exists = structure.bonds.some(
      (b) =>
        (b.atom_i === atomI && b.atom_j === atomJ) ||
        (b.atom_i === atomJ && b.atom_j === atomI),
    );
    if (exists) return false;

    const newBond: BondResponse = {
      atom_i: Math.min(atomI, atomJ),
      atom_j: Math.max(atomI, atomJ),
      order,
      rotatable: order === 1,
    };

    const newStructure: Molecule3DResponse = {
      ...structure,
      bonds: [...structure.bonds, newBond],
    };

    const molecules = new Map(get().molecules);
    molecules.set(entry.id, {
      ...entry,
      structure: newStructure,
      molecule: {
        ...entry.molecule,
        num_bonds: newStructure.bonds.length,
      },
    });
    set({ molecules });
    return true;
  },

  removeAtoms(indices) {
    const entry = get().getActive();
    if (!entry?.structure || indices.size === 0) return;

    const structure = entry.structure;

    // Remove atoms and reindex
    const remaining = structure.atoms.filter((a) => !indices.has(a.index));
    const oldToNew = new Map<number, number>();
    remaining.forEach((a, i) => {
      oldToNew.set(a.index, i);
    });

    const newAtoms: AtomResponse[] = remaining.map((a, i) => ({
      ...a,
      index: i,
      position: [...a.position] as [number, number, number],
    }));

    // Remove bonds connected to deleted atoms and reindex
    const newBonds: BondResponse[] = structure.bonds
      .filter((b) => !indices.has(b.atom_i) && !indices.has(b.atom_j))
      .map((b) => ({
        ...b,
        atom_i: oldToNew.get(b.atom_i)!,
        atom_j: oldToNew.get(b.atom_j)!,
      }));

    const newStructure: Molecule3DResponse = {
      ...structure,
      atoms: newAtoms,
      bonds: newBonds,
    };

    const molecules = new Map(get().molecules);
    molecules.set(entry.id, {
      ...entry,
      structure: newStructure,
      molecule: {
        ...entry.molecule,
        num_atoms: newAtoms.length,
        num_bonds: newBonds.length,
      },
    });
    set({ molecules });
  },

  removeBonds(keys) {
    const entry = get().getActive();
    if (!entry?.structure || keys.size === 0) return;

    const structure = entry.structure;
    const newBonds = structure.bonds.filter((b) => {
      const key1 = `${b.atom_i}-${b.atom_j}`;
      const key2 = `${b.atom_j}-${b.atom_i}`;
      return !keys.has(key1) && !keys.has(key2);
    });

    const newStructure: Molecule3DResponse = {
      ...structure,
      bonds: newBonds,
    };

    const molecules = new Map(get().molecules);
    molecules.set(entry.id, {
      ...entry,
      structure: newStructure,
      molecule: {
        ...entry.molecule,
        num_bonds: newBonds.length,
      },
    });
    set({ molecules });
  },

  updateActiveStructure(structure) {
    const entry = get().getActive();
    if (!entry) return;
    const molecules = new Map(get().molecules);
    molecules.set(entry.id, {
      ...entry,
      structure,
      molecule: {
        ...entry.molecule,
        num_atoms: structure.atoms.length,
        num_bonds: structure.bonds.length,
      },
    });
    set({ molecules });
  },

  replaceAtomElement(index, element) {
    const entry = get().getActive();
    if (!entry?.structure) return;

    const structure = entry.structure;
    if (index < 0 || index >= structure.atoms.length) return;

    const newAtoms: AtomResponse[] = structure.atoms.map((a) =>
      a.index === index
        ? { ...a, symbol: element, position: [...a.position] as [number, number, number] }
        : a,
    );

    const newStructure: Molecule3DResponse = { ...structure, atoms: newAtoms };
    const molecules = new Map(get().molecules);
    molecules.set(entry.id, { ...entry, structure: newStructure });
    set({ molecules });
  },

  moveAtom(index, position) {
    const entry = get().getActive();
    if (!entry?.structure) return;
    if (index < 0 || index >= entry.structure.atoms.length) return;

    const newAtoms = entry.structure.atoms.map((a) =>
      a.index === index ? { ...a, position: [...position] as [number, number, number] } : a,
    );
    const newStructure = { ...entry.structure, atoms: newAtoms };
    const molecules = new Map(get().molecules);
    molecules.set(entry.id, { ...entry, structure: newStructure });
    set({ molecules });
  },

  insertFragment(parentAtomIndex, fragment) {
    const entry = get().getActive();
    if (!entry?.structure) return;
    const structure = entry.structure;
    const parentAtom = structure.atoms[parentAtomIndex];
    if (!parentAtom) return;

    const baseIndex = structure.atoms.length;
    const parentPos = parentAtom.position;

    // Find a direction away from existing bonds for placement
    let dirX = 1, dirY = 0, dirZ = 0;
    const bonded = structure.bonds
      .filter((b) => b.atom_i === parentAtomIndex || b.atom_j === parentAtomIndex)
      .map((b) => b.atom_i === parentAtomIndex ? b.atom_j : b.atom_i);
    if (bonded.length > 0) {
      // Average direction of existing bonds, then negate to go away
      let avgX = 0, avgY = 0, avgZ = 0;
      for (const bi of bonded) {
        const a = structure.atoms[bi];
        if (a) {
          avgX += a.position[0] - parentPos[0];
          avgY += a.position[1] - parentPos[1];
          avgZ += a.position[2] - parentPos[2];
        }
      }
      const len = Math.sqrt(avgX * avgX + avgY * avgY + avgZ * avgZ);
      if (len > 0.01) {
        dirX = -avgX / len;
        dirY = -avgY / len;
        dirZ = -avgZ / len;
      }
    }

    // Create new atoms at offset positions rotated into the placement direction
    const newAtoms: AtomResponse[] = fragment.atoms.map((fa, i) => ({
      index: baseIndex + i,
      symbol: fa.symbol,
      position: [
        parentPos[0] + fa.offset[0] * dirX - fa.offset[2] * dirZ,
        parentPos[1] + fa.offset[0] * dirY + fa.offset[1],
        parentPos[2] + fa.offset[0] * dirZ + fa.offset[2] * dirX,
      ] as [number, number, number],
      hybridization: null,
      formal_charge: 0,
    }));

    // Bond from parent to first fragment atom
    const newBonds: BondResponse[] = [
      { atom_i: parentAtomIndex, atom_j: baseIndex, order: 1, rotatable: true },
    ];

    // Internal fragment bonds
    for (const [fi, fj, order] of fragment.bonds) {
      newBonds.push({
        atom_i: baseIndex + fi,
        atom_j: baseIndex + fj,
        order,
        rotatable: order === 1,
      });
    }

    const newStructure: Molecule3DResponse = {
      ...structure,
      atoms: [...structure.atoms, ...newAtoms],
      bonds: [...structure.bonds, ...newBonds],
    };

    const molecules = new Map(get().molecules);
    molecules.set(entry.id, {
      ...entry,
      structure: newStructure,
      molecule: {
        ...entry.molecule,
        num_atoms: newStructure.atoms.length,
        num_bonds: newStructure.bonds.length,
      },
    });
    set({ molecules });
  },

  cycleBondOrder(atomI, atomJ) {
    const entry = get().getActive();
    if (!entry?.structure) return;

    const structure = entry.structure;
    const bondIdx = structure.bonds.findIndex(
      (b) =>
        (b.atom_i === atomI && b.atom_j === atomJ) ||
        (b.atom_i === atomJ && b.atom_j === atomI),
    );
    if (bondIdx < 0) return;

    const bond = structure.bonds[bondIdx]!;
    const nextOrder = bond.order >= 3 ? 1 : bond.order + 1;

    const newBonds: BondResponse[] = structure.bonds.map((b, i) =>
      i === bondIdx ? { ...b, order: nextOrder, rotatable: nextOrder === 1 } : b,
    );

    const newStructure: Molecule3DResponse = { ...structure, bonds: newBonds };
    const molecules = new Map(get().molecules);
    molecules.set(entry.id, { ...entry, structure: newStructure });
    set({ molecules });
  },
}));
