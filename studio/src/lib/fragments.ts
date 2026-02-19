/**
 * Fragment definitions for functional group insertion.
 * Each fragment defines atoms (relative positions from attachment point)
 * and internal bonds. The first atom is the attachment point that bonds
 * to the selected atom.
 */

export interface FragmentAtom {
  symbol: string;
  offset: [number, number, number]; // relative to attachment direction
}

export interface FragmentDef {
  name: string;
  shortLabel: string;
  atoms: FragmentAtom[];
  /** Internal bonds as [fromIdx, toIdx, order] within the fragment */
  bonds: [number, number, number][];
}

// Standard C-X bond length in Angstroms for placing fragments
const CC = 1.54;
const CO = 1.43;
const CN = 1.47;
const CH = 1.09;

export const FRAGMENTS: FragmentDef[] = [
  {
    name: "Methyl (-CH3)",
    shortLabel: "-CH3",
    atoms: [
      { symbol: "C", offset: [CC, 0, 0] },
      { symbol: "H", offset: [CC + CH * 0.33, CH * 0.94, 0] },
      { symbol: "H", offset: [CC + CH * 0.33, -CH * 0.47, CH * 0.82] },
      { symbol: "H", offset: [CC + CH * 0.33, -CH * 0.47, -CH * 0.82] },
    ],
    bonds: [[0, 1, 1], [0, 2, 1], [0, 3, 1]],
  },
  {
    name: "Hydroxyl (-OH)",
    shortLabel: "-OH",
    atoms: [
      { symbol: "O", offset: [CO, 0, 0] },
      { symbol: "H", offset: [CO + 0.96, 0, 0] },
    ],
    bonds: [[0, 1, 1]],
  },
  {
    name: "Amino (-NH2)",
    shortLabel: "-NH2",
    atoms: [
      { symbol: "N", offset: [CN, 0, 0] },
      { symbol: "H", offset: [CN + 1.01 * 0.5, 1.01 * 0.87, 0] },
      { symbol: "H", offset: [CN + 1.01 * 0.5, -1.01 * 0.87, 0] },
    ],
    bonds: [[0, 1, 1], [0, 2, 1]],
  },
  {
    name: "Carboxyl (-COOH)",
    shortLabel: "-COOH",
    atoms: [
      { symbol: "C", offset: [CC, 0, 0] },
      { symbol: "O", offset: [CC + 1.23 * 0.5, 1.23 * 0.87, 0] },
      { symbol: "O", offset: [CC + CO, 0, 0] },
      { symbol: "H", offset: [CC + CO + 0.96, 0, 0] },
    ],
    bonds: [[0, 1, 2], [0, 2, 1], [2, 3, 1]],
  },
  {
    name: "Aldehyde (-CHO)",
    shortLabel: "-CHO",
    atoms: [
      { symbol: "C", offset: [CC, 0, 0] },
      { symbol: "O", offset: [CC + 1.23 * 0.5, 1.23 * 0.87, 0] },
      { symbol: "H", offset: [CC + CH, 0, 0] },
    ],
    bonds: [[0, 1, 2], [0, 2, 1]],
  },
  {
    name: "Nitro (-NO2)",
    shortLabel: "-NO2",
    atoms: [
      { symbol: "N", offset: [CN, 0, 0] },
      { symbol: "O", offset: [CN + 1.21 * 0.5, 1.21 * 0.87, 0] },
      { symbol: "O", offset: [CN + 1.21 * 0.5, -1.21 * 0.87, 0] },
    ],
    bonds: [[0, 1, 2], [0, 2, 1]],
  },
  {
    name: "Fluorine (-F)",
    shortLabel: "-F",
    atoms: [{ symbol: "F", offset: [1.35, 0, 0] }],
    bonds: [],
  },
  {
    name: "Chlorine (-Cl)",
    shortLabel: "-Cl",
    atoms: [{ symbol: "Cl", offset: [1.77, 0, 0] }],
    bonds: [],
  },
  {
    name: "Bromine (-Br)",
    shortLabel: "-Br",
    atoms: [{ symbol: "Br", offset: [1.94, 0, 0] }],
    bonds: [],
  },
  {
    name: "Cyano (-CN)",
    shortLabel: "-CN",
    atoms: [
      { symbol: "C", offset: [CC, 0, 0] },
      { symbol: "N", offset: [CC + 1.16, 0, 0] },
    ],
    bonds: [[0, 1, 3]],
  },
];
