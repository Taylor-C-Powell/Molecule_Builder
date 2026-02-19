import { parseMolV2000 } from "./mol-import";
import type { Molecule3DResponse } from "@/api/types";

/**
 * Parse an SDF file string. Extracts the first MOL block (before $$$$).
 * SDF files contain one or more MOL blocks followed by property data and $$$$ delimiters.
 */
export function parseSdf(content: string): { structure: Molecule3DResponse; name: string } {
  // Find the first $$$$ delimiter
  const delimIndex = content.indexOf("$$$$");
  const molBlock = delimIndex >= 0 ? content.substring(0, delimIndex) : content;

  // Strip trailing M  END and any data fields after it
  const endIndex = molBlock.indexOf("M  END");
  const cleanBlock = endIndex >= 0 ? molBlock.substring(0, endIndex + 6) : molBlock;

  return parseMolV2000(cleanBlock);
}
