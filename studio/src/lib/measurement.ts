import { Vector3 } from "three";
import type { AtomResponse } from "@/api/types";
import type { Measurement } from "@/types/editor";

export function computeDistance(a: AtomResponse, b: AtomResponse): number {
  const va = new Vector3(...a.position);
  const vb = new Vector3(...b.position);
  return va.distanceTo(vb);
}

export function computeAngle(a: AtomResponse, b: AtomResponse, c: AtomResponse): number {
  const va = new Vector3(...a.position);
  const vb = new Vector3(...b.position);
  const vc = new Vector3(...c.position);
  const ba = new Vector3().subVectors(va, vb).normalize();
  const bc = new Vector3().subVectors(vc, vb).normalize();
  const dot = ba.dot(bc);
  return Math.acos(Math.max(-1, Math.min(1, dot))) * (180 / Math.PI);
}

export function computeDihedral(
  a: AtomResponse,
  b: AtomResponse,
  c: AtomResponse,
  d: AtomResponse,
): number {
  const p1 = new Vector3(...a.position);
  const p2 = new Vector3(...b.position);
  const p3 = new Vector3(...c.position);
  const p4 = new Vector3(...d.position);

  const b1 = new Vector3().subVectors(p2, p1);
  const b2 = new Vector3().subVectors(p3, p2);
  const b3 = new Vector3().subVectors(p4, p3);

  const n1 = new Vector3().crossVectors(b1, b2).normalize();
  const n2 = new Vector3().crossVectors(b2, b3).normalize();
  const m = new Vector3().crossVectors(n1, b2.normalize());

  const x = n1.dot(n2);
  const y = m.dot(n2);
  return Math.atan2(y, x) * (180 / Math.PI);
}

export function buildMeasurement(
  atoms: AtomResponse[],
  indices: number[],
): Measurement | null {
  const resolved = indices.map((i) => atoms[i]).filter(Boolean);

  if (resolved.length === 2) {
    return {
      type: "distance",
      atomIndices: indices,
      value: computeDistance(resolved[0]!, resolved[1]!),
      unit: "A",
    };
  }

  if (resolved.length === 3) {
    return {
      type: "angle",
      atomIndices: indices,
      value: computeAngle(resolved[0]!, resolved[1]!, resolved[2]!),
      unit: "deg",
    };
  }

  if (resolved.length === 4) {
    return {
      type: "dihedral",
      atomIndices: indices,
      value: computeDihedral(resolved[0]!, resolved[1]!, resolved[2]!, resolved[3]!),
      unit: "deg",
    };
  }

  return null;
}
