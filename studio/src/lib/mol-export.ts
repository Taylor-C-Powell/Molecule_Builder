import type { Molecule3DResponse } from "@/api/types";

/**
 * Generate a V2000 MOL file string from 3D molecule data.
 * Spec: https://en.wikipedia.org/wiki/Chemical_table_file#Molfile
 */
export function generateMolV2000(structure: Molecule3DResponse, name?: string): string {
  const lines: string[] = [];

  // Header block (3 lines)
  lines.push(name ?? "molecule");
  lines.push("  MolBuilder3D");
  lines.push("");

  // Counts line
  const nAtoms = structure.atoms.length;
  const nBonds = structure.bonds.length;
  lines.push(
    `${pad(nAtoms, 3)}${pad(nBonds, 3)}  0  0  0  0  0  0  0  0999 V2000`,
  );

  // Atom block
  for (const atom of structure.atoms) {
    const [x, y, z] = atom.position;
    lines.push(
      `${padF(x, 10, 4)}${padF(y, 10, 4)}${padF(z, 10, 4)} ${padR(atom.symbol, 3)}` +
        ` 0  0  0  0  0  0  0  0  0  0  0  0`,
    );
  }

  // Bond block (1-indexed)
  for (const bond of structure.bonds) {
    const order = Math.min(Math.max(Math.round(bond.order), 1), 3);
    lines.push(
      `${pad(bond.atom_i + 1, 3)}${pad(bond.atom_j + 1, 3)}${pad(order, 3)}  0  0  0  0`,
    );
  }

  lines.push("M  END");

  return lines.join("\n");
}

function pad(n: number, width: number): string {
  return n.toString().padStart(width);
}

function padF(n: number, width: number, decimals: number): string {
  return n.toFixed(decimals).padStart(width);
}

function padR(s: string, width: number): string {
  return s.padEnd(width);
}
