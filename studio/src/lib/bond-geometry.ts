import { Quaternion, Vector3 } from "three";

export interface BondTransform {
  position: [number, number, number];
  quaternion: [number, number, number, number];
  length: number;
}

const _up = new Vector3(0, 1, 0);
const _dir = new Vector3();
const _mid = new Vector3();
const _q = new Quaternion();

/**
 * Compute the transform (position, quaternion, length) for a cylinder
 * between two 3D atom positions. The cylinder is oriented along its Y axis
 * by default (Three.js CylinderGeometry convention).
 */
export function computeBondTransform(
  posA: [number, number, number],
  posB: [number, number, number],
): BondTransform {
  const a = new Vector3(...posA);
  const b = new Vector3(...posB);

  _dir.subVectors(b, a);
  const length = _dir.length();
  _dir.normalize();

  _mid.addVectors(a, b).multiplyScalar(0.5);

  _q.setFromUnitVectors(_up, _dir);

  return {
    position: [_mid.x, _mid.y, _mid.z],
    quaternion: [_q.x, _q.y, _q.z, _q.w],
    length,
  };
}

/**
 * Compute an offset perpendicular to the bond axis for double/triple bonds.
 * Returns a unit vector perpendicular to the bond direction.
 */
export function computeBondOffset(
  posA: [number, number, number],
  posB: [number, number, number],
  offset: number,
): [number, number, number] {
  const a = new Vector3(...posA);
  const b = new Vector3(...posB);
  const dir = new Vector3().subVectors(b, a).normalize();

  // Find a perpendicular vector
  const perp = new Vector3();
  if (Math.abs(dir.x) < 0.9) {
    perp.set(1, 0, 0);
  } else {
    perp.set(0, 1, 0);
  }
  perp.cross(dir).normalize().multiplyScalar(offset);

  return [perp.x, perp.y, perp.z];
}
