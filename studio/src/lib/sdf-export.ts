import type { Molecule3DResponse, MoleculePropertiesResponse } from "@/api/types";
import { generateMolV2000 } from "./mol-export";

/**
 * Generate an SDF file string. Wraps the MOL block with optional property fields.
 */
export function generateSdf(
  structure: Molecule3DResponse,
  name?: string,
  properties?: MoleculePropertiesResponse,
): string {
  const lines: string[] = [];

  lines.push(generateMolV2000(structure, name));

  if (properties) {
    if (properties.formula) {
      lines.push(`> <FORMULA>`);
      lines.push(properties.formula);
      lines.push("");
    }
    if (properties.molecular_weight) {
      lines.push(`> <MOLECULAR_WEIGHT>`);
      lines.push(properties.molecular_weight.toFixed(4));
      lines.push("");
    }
    if (properties.smiles) {
      lines.push(`> <SMILES>`);
      lines.push(properties.smiles);
      lines.push("");
    }
    if (properties.logp != null) {
      lines.push(`> <LOGP>`);
      lines.push(properties.logp.toFixed(2));
      lines.push("");
    }
  }

  lines.push("$$$$");

  return lines.join("\n");
}
