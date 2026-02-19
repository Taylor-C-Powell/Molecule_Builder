import type { AtomResponse, BondResponse, Molecule3DResponse } from "@/api/types";

/**
 * Parse a V2000 MOL file string into a Molecule3DResponse.
 * Handles the standard header (3 lines), counts line, atom block, and bond block.
 */
export function parseMolV2000(content: string): { structure: Molecule3DResponse; name: string } {
  const lines = content.split(/\r?\n/);

  if (lines.length < 5) {
    throw new Error("Invalid MOL file: too few lines");
  }

  const name = lines[0]?.trim() || "imported";

  // Counts line is line index 3
  const countsLine = lines[3]!;
  const nAtoms = parseInt(countsLine.substring(0, 3).trim(), 10);
  const nBonds = parseInt(countsLine.substring(3, 6).trim(), 10);

  if (isNaN(nAtoms) || isNaN(nBonds)) {
    throw new Error("Invalid MOL file: cannot parse atom/bond counts");
  }

  const atoms: AtomResponse[] = [];
  for (let i = 0; i < nAtoms; i++) {
    const line = lines[4 + i];
    if (!line) throw new Error(`Invalid MOL file: missing atom line ${i}`);

    const x = parseFloat(line.substring(0, 10).trim());
    const y = parseFloat(line.substring(10, 20).trim());
    const z = parseFloat(line.substring(20, 30).trim());
    const symbol = line.substring(31, 34).trim();

    if (isNaN(x) || isNaN(y) || isNaN(z) || !symbol) {
      throw new Error(`Invalid MOL file: cannot parse atom at line ${4 + i}`);
    }

    atoms.push({
      index: i,
      symbol,
      position: [x, y, z],
      hybridization: null,
      formal_charge: 0,
    });
  }

  const bonds: BondResponse[] = [];
  const bondStart = 4 + nAtoms;
  for (let i = 0; i < nBonds; i++) {
    const line = lines[bondStart + i];
    if (!line) throw new Error(`Invalid MOL file: missing bond line ${i}`);

    const atomI = parseInt(line.substring(0, 3).trim(), 10) - 1; // 1-indexed to 0-indexed
    const atomJ = parseInt(line.substring(3, 6).trim(), 10) - 1;
    const order = parseInt(line.substring(6, 9).trim(), 10);

    if (isNaN(atomI) || isNaN(atomJ) || isNaN(order)) {
      throw new Error(`Invalid MOL file: cannot parse bond at line ${bondStart + i}`);
    }

    bonds.push({
      atom_i: Math.min(atomI, atomJ),
      atom_j: Math.max(atomI, atomJ),
      order: Math.min(Math.max(order, 1), 3),
      rotatable: order === 1,
    });
  }

  const id = crypto.randomUUID();

  return {
    name,
    structure: { id, atoms, bonds },
  };
}
