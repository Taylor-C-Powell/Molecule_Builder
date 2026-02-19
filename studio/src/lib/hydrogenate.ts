import type { AtomResponse, BondResponse, Molecule3DResponse } from "@/api/types";

/**
 * Typical valences for common elements.
 * Returns the max number of bonds an element normally forms.
 */
const VALENCE: Record<string, number> = {
  H: 1,
  C: 4,
  N: 3,
  O: 2,
  S: 2,
  P: 3,
  F: 1,
  Cl: 1,
  Br: 1,
  I: 1,
  Si: 4,
  B: 3,
};

/**
 * Calculate how many hydrogens are missing for each atom based on
 * typical valence rules minus current bond order sum.
 * Returns a new Molecule3DResponse with hydrogens added.
 */
export function addMissingHydrogens(structure: Molecule3DResponse): Molecule3DResponse {
  // Calculate bond order sum for each atom
  const bondOrderSum = new Map<number, number>();
  for (const atom of structure.atoms) {
    bondOrderSum.set(atom.index, 0);
  }
  for (const bond of structure.bonds) {
    bondOrderSum.set(bond.atom_i, (bondOrderSum.get(bond.atom_i) ?? 0) + bond.order);
    bondOrderSum.set(bond.atom_j, (bondOrderSum.get(bond.atom_j) ?? 0) + bond.order);
  }

  // Build adjacency for direction calculation
  const adj = new Map<number, number[]>();
  for (const bond of structure.bonds) {
    if (!adj.has(bond.atom_i)) adj.set(bond.atom_i, []);
    if (!adj.has(bond.atom_j)) adj.set(bond.atom_j, []);
    adj.get(bond.atom_i)!.push(bond.atom_j);
    adj.get(bond.atom_j)!.push(bond.atom_i);
  }

  const newAtoms: AtomResponse[] = [...structure.atoms];
  const newBonds: BondResponse[] = [...structure.bonds];
  let nextIndex = structure.atoms.length;

  for (const atom of structure.atoms) {
    if (atom.symbol === "H") continue; // don't add H to H

    const maxValence = VALENCE[atom.symbol];
    if (maxValence == null) continue; // unknown element, skip

    const currentSum = bondOrderSum.get(atom.index) ?? 0;
    const formalCharge = atom.formal_charge ?? 0;
    const missing = maxValence - currentSum + formalCharge;

    if (missing <= 0) continue;

    // Calculate placement directions
    const neighbors = adj.get(atom.index) ?? [];
    const hPositions = computeHydrogenPositions(
      atom.position,
      neighbors.map((ni) => structure.atoms[ni]!.position),
      missing,
    );

    for (const pos of hPositions) {
      newAtoms.push({
        index: nextIndex,
        symbol: "H",
        position: pos,
        hybridization: null,
        formal_charge: 0,
      });
      newBonds.push({
        atom_i: atom.index,
        atom_j: nextIndex,
        order: 1,
        rotatable: false,
      });
      nextIndex++;
    }
  }

  return { ...structure, atoms: newAtoms, bonds: newBonds };
}

/**
 * Compute reasonable 3D positions for hydrogen atoms around a parent.
 * Uses existing neighbor directions to place hydrogens in sensible orientations.
 */
function computeHydrogenPositions(
  parentPos: [number, number, number],
  neighborPositions: [number, number, number][],
  count: number,
): [number, number, number][] {
  const H_BOND_LENGTH = 1.09;
  const positions: [number, number, number][] = [];

  if (neighborPositions.length === 0) {
    // No neighbors: place around parent using tetrahedral geometry
    const dirs: [number, number, number][] = [
      [1, 0, 0],
      [-0.33, 0.94, 0],
      [-0.33, -0.47, 0.82],
      [-0.33, -0.47, -0.82],
    ];
    for (let i = 0; i < count && i < dirs.length; i++) {
      const d = dirs[i]!;
      positions.push([
        parentPos[0] + d[0] * H_BOND_LENGTH,
        parentPos[1] + d[1] * H_BOND_LENGTH,
        parentPos[2] + d[2] * H_BOND_LENGTH,
      ]);
    }
    return positions;
  }

  // Average direction of existing neighbors
  let avgX = 0, avgY = 0, avgZ = 0;
  for (const np of neighborPositions) {
    avgX += np[0] - parentPos[0];
    avgY += np[1] - parentPos[1];
    avgZ += np[2] - parentPos[2];
  }
  let len = Math.sqrt(avgX * avgX + avgY * avgY + avgZ * avgZ);
  if (len < 0.001) len = 1;
  avgX /= len; avgY /= len; avgZ /= len;

  // Primary H direction: opposite to average neighbor direction
  const awayX = -avgX, awayY = -avgY, awayZ = -avgZ;

  // Find a perpendicular vector
  let perpX: number, perpY: number, perpZ: number;
  if (Math.abs(awayX) < 0.9) {
    perpX = 0; perpY = -awayZ; perpZ = awayY;
  } else {
    perpX = -awayZ; perpY = 0; perpZ = awayX;
  }
  const pLen = Math.sqrt(perpX * perpX + perpY * perpY + perpZ * perpZ);
  if (pLen > 0.001) {
    perpX /= pLen; perpY /= pLen; perpZ /= pLen;
  }

  if (count === 1) {
    positions.push([
      parentPos[0] + awayX * H_BOND_LENGTH,
      parentPos[1] + awayY * H_BOND_LENGTH,
      parentPos[2] + awayZ * H_BOND_LENGTH,
    ]);
  } else {
    // Spread hydrogens in a cone around the "away" direction
    const angleStep = (2 * Math.PI) / count;
    const coneAngle = 0.6; // ~34 degrees from away axis

    for (let i = 0; i < count; i++) {
      const theta = angleStep * i;
      const sinCone = Math.sin(coneAngle);
      const cosCone = Math.cos(coneAngle);

      // Cross product to get second perpendicular
      const perp2X = awayY * perpZ - awayZ * perpY;
      const perp2Y = awayZ * perpX - awayX * perpZ;
      const perp2Z = awayX * perpY - awayY * perpX;

      const dx = awayX * cosCone + (perpX * Math.cos(theta) + perp2X * Math.sin(theta)) * sinCone;
      const dy = awayY * cosCone + (perpY * Math.cos(theta) + perp2Y * Math.sin(theta)) * sinCone;
      const dz = awayZ * cosCone + (perpZ * Math.cos(theta) + perp2Z * Math.sin(theta)) * sinCone;

      positions.push([
        parentPos[0] + dx * H_BOND_LENGTH,
        parentPos[1] + dy * H_BOND_LENGTH,
        parentPos[2] + dz * H_BOND_LENGTH,
      ]);
    }
  }

  return positions;
}
