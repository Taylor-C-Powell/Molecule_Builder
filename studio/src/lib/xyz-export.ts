import type { Molecule3DResponse } from "@/api/types";

/**
 * Generate an XYZ file string from 3D molecule data.
 * Format: atom_count\ncomment\nsymbol x y z (per atom)
 */
export function generateXyz(structure: Molecule3DResponse, name?: string): string {
  const lines: string[] = [];

  lines.push(structure.atoms.length.toString());
  lines.push(name ?? "molecule");

  for (const atom of structure.atoms) {
    const [x, y, z] = atom.position;
    lines.push(`${atom.symbol}  ${x.toFixed(6)}  ${y.toFixed(6)}  ${z.toFixed(6)}`);
  }

  return lines.join("\n");
}
