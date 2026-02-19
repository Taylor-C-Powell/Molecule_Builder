export enum Tool {
  SELECT = "select",
  ROTATE = "rotate",
  MEASURE = "measure",
  ADD_ATOM = "add-atom",
  ADD_BOND = "add-bond",
  ERASE = "erase",
}

export type RenderStyle = "ball-and-stick" | "stick" | "spacefill";

export type BondOrder = 1 | 2 | 3;

export const PLACEABLE_ELEMENTS = [
  "C", "N", "O", "H", "S", "P", "F", "Cl", "Br", "I",
] as const;

export type PlaceableElement = (typeof PLACEABLE_ELEMENTS)[number];

export interface Measurement {
  type: "distance" | "angle" | "dihedral";
  atomIndices: number[];
  value: number;
  unit: string;
}

export interface ViewportSettings {
  showHydrogens: boolean;
  showLabels: boolean;
  renderStyle: RenderStyle;
}
