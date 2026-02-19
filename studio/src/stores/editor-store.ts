import { create } from "zustand";
import { Tool } from "@/types/editor";
import type { RenderStyle, BondOrder, PlaceableElement, Measurement } from "@/types/editor";
import type { FragmentDef } from "@/lib/fragments";

interface EditorState {
  tool: Tool;
  selectedAtoms: Set<number>;
  selectedBonds: Set<string>;
  hoveredAtom: number | null;
  hoveredBond: string | null;
  showHydrogens: boolean;
  showLabels: boolean;
  renderStyle: RenderStyle;

  // Phase 2: editing state
  activeElement: PlaceableElement;
  activeBondOrder: BondOrder;
  bondStartAtom: number | null;
  measurements: Measurement[];
  measureAtoms: number[];
  activeFragment: FragmentDef | null;

  setTool: (t: Tool) => void;
  selectAtom: (idx: number, multi?: boolean) => void;
  selectBond: (key: string, multi?: boolean) => void;
  clearSelection: () => void;
  setHoveredAtom: (idx: number | null) => void;
  setHoveredBond: (key: string | null) => void;
  toggleHydrogens: () => void;
  toggleLabels: () => void;
  setRenderStyle: (s: RenderStyle) => void;

  // Phase 2: editing actions
  setActiveElement: (el: PlaceableElement) => void;
  setActiveBondOrder: (order: BondOrder) => void;
  setBondStartAtom: (idx: number | null) => void;
  addMeasurement: (m: Measurement) => void;
  clearMeasurements: () => void;
  pushMeasureAtom: (idx: number) => void;
  clearMeasureAtoms: () => void;
  setActiveFragment: (f: FragmentDef | null) => void;

  // Phase 5
  showDistanceLabels: boolean;
  autoRotate: boolean;
  showGrid: boolean;
  toggleDistanceLabels: () => void;
  toggleAutoRotate: () => void;
  toggleGrid: () => void;
  selectAllAtoms: (totalAtoms: number) => void;
  selectByElement: (symbol: string, atoms: { index: number; symbol: string }[]) => void;
  invertSelection: (totalAtoms: number) => void;

  // Context menu
  contextMenu: { x: number; y: number; type: "atom" | "bond"; index: number; bondKey?: string } | null;
  openContextMenu: (x: number, y: number, type: "atom" | "bond", index: number, bondKey?: string) => void;
  closeContextMenu: () => void;
}

export const useEditorStore = create<EditorState>((set, get) => ({
  tool: Tool.SELECT,
  selectedAtoms: new Set(),
  selectedBonds: new Set(),
  hoveredAtom: null,
  hoveredBond: null,
  showHydrogens: true,
  showLabels: false,
  renderStyle: "ball-and-stick",

  activeElement: "C",
  activeBondOrder: 1,
  bondStartAtom: null,
  measurements: [],
  measureAtoms: [],
  activeFragment: null,
  showDistanceLabels: false,
  autoRotate: false,
  showGrid: false,
  contextMenu: null,

  setTool(t) {
    set({ tool: t, bondStartAtom: null, measureAtoms: [] });
  },

  selectAtom(idx, multi = false) {
    const current = get().selectedAtoms;
    if (multi) {
      const next = new Set(current);
      if (next.has(idx)) {
        next.delete(idx);
      } else {
        next.add(idx);
      }
      set({ selectedAtoms: next });
    } else {
      set({ selectedAtoms: new Set([idx]), selectedBonds: new Set() });
    }
  },

  selectBond(key, multi = false) {
    const current = get().selectedBonds;
    if (multi) {
      const next = new Set(current);
      if (next.has(key)) {
        next.delete(key);
      } else {
        next.add(key);
      }
      set({ selectedBonds: next });
    } else {
      set({ selectedBonds: new Set([key]), selectedAtoms: new Set() });
    }
  },

  clearSelection() {
    set({ selectedAtoms: new Set(), selectedBonds: new Set() });
  },

  setHoveredAtom(idx) {
    set({ hoveredAtom: idx });
  },

  setHoveredBond(key) {
    set({ hoveredBond: key });
  },

  toggleHydrogens() {
    set({ showHydrogens: !get().showHydrogens });
  },

  toggleLabels() {
    set({ showLabels: !get().showLabels });
  },

  setRenderStyle(s) {
    set({ renderStyle: s });
  },

  setActiveElement(el) {
    set({ activeElement: el });
  },

  setActiveBondOrder(order) {
    set({ activeBondOrder: order });
  },

  setBondStartAtom(idx) {
    set({ bondStartAtom: idx });
  },

  addMeasurement(m) {
    set({ measurements: [...get().measurements, m] });
  },

  clearMeasurements() {
    set({ measurements: [], measureAtoms: [] });
  },

  pushMeasureAtom(idx) {
    const current = get().measureAtoms;
    if (current.includes(idx)) return;
    set({ measureAtoms: [...current, idx] });
  },

  clearMeasureAtoms() {
    set({ measureAtoms: [] });
  },

  setActiveFragment(f) {
    set({ activeFragment: f });
  },

  toggleDistanceLabels() {
    set({ showDistanceLabels: !get().showDistanceLabels });
  },

  toggleAutoRotate() {
    set({ autoRotate: !get().autoRotate });
  },

  toggleGrid() {
    set({ showGrid: !get().showGrid });
  },

  selectAllAtoms(totalAtoms) {
    const all = new Set<number>();
    for (let i = 0; i < totalAtoms; i++) all.add(i);
    set({ selectedAtoms: all, selectedBonds: new Set() });
  },

  selectByElement(symbol, atoms) {
    const matching = new Set(
      atoms.filter((a) => a.symbol === symbol).map((a) => a.index),
    );
    set({ selectedAtoms: matching, selectedBonds: new Set() });
  },

  invertSelection(totalAtoms) {
    const current = get().selectedAtoms;
    const inverted = new Set<number>();
    for (let i = 0; i < totalAtoms; i++) {
      if (!current.has(i)) inverted.add(i);
    }
    set({ selectedAtoms: inverted, selectedBonds: new Set() });
  },

  openContextMenu(x, y, type, index, bondKey) {
    set({ contextMenu: { x, y, type, index, bondKey } });
  },

  closeContextMenu() {
    set({ contextMenu: null });
  },
}));
