import type { AtomResponse, Molecule3DResponse } from "@/api/types";

/**
 * Parse an XYZ file string into a Molecule3DResponse.
 * Format: line 1 = atom count, line 2 = comment, lines 3+ = symbol x y z
 * Note: XYZ format has no bond information, so bonds array will be empty.
 */
export function parseXyz(content: string): { structure: Molecule3DResponse; name: string } {
  const lines = content.split(/\r?\n/).filter((l) => l.trim().length > 0);

  if (lines.length < 3) {
    throw new Error("Invalid XYZ file: too few lines");
  }

  const nAtoms = parseInt(lines[0]!.trim(), 10);
  if (isNaN(nAtoms) || nAtoms <= 0) {
    throw new Error("Invalid XYZ file: cannot parse atom count");
  }

  const name = lines[1]?.trim() || "imported";

  const atoms: AtomResponse[] = [];
  for (let i = 0; i < nAtoms; i++) {
    const line = lines[2 + i];
    if (!line) throw new Error(`Invalid XYZ file: missing atom line ${i}`);

    const parts = line.trim().split(/\s+/);
    if (parts.length < 4) {
      throw new Error(`Invalid XYZ file: cannot parse atom at line ${2 + i}`);
    }

    const symbol = parts[0]!;
    const x = parseFloat(parts[1]!);
    const y = parseFloat(parts[2]!);
    const z = parseFloat(parts[3]!);

    if (isNaN(x) || isNaN(y) || isNaN(z)) {
      throw new Error(`Invalid XYZ file: invalid coordinates at line ${2 + i}`);
    }

    atoms.push({
      index: i,
      symbol,
      position: [x, y, z],
      hybridization: null,
      formal_charge: 0,
    });
  }

  const id = crypto.randomUUID();

  return {
    name,
    structure: { id, atoms, bonds: [] },
  };
}
