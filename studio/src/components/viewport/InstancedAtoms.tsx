import { useRef, useMemo, useEffect } from "react";
import { InstancedMesh, Matrix4, Color, SphereGeometry, MeshStandardMaterial } from "three";
import type { AtomResponse } from "@/api/types";
import { getCpkColor, getAtomRadius } from "@/lib/cpk-colors";
import { useSettingsStore } from "@/stores/settings-store";

interface InstancedAtomsProps {
  atoms: AtomResponse[];
  renderStyle: "ball-and-stick" | "stick" | "spacefill";
}

/**
 * Renders all atoms using a single InstancedMesh per unique geometry config.
 * Much faster than individual meshes for large molecules (>200 atoms).
 * Falls back gracefully -- used alongside AtomMesh components for interaction.
 */
export function InstancedAtoms({ atoms, renderStyle }: InstancedAtomsProps) {
  const meshRef = useRef<InstancedMesh>(null);
  const atomScale = useSettingsStore((s) => s.atomScale);

  const geometry = useMemo(() => new SphereGeometry(1, 16, 16), []);
  const material = useMemo(
    () => new MeshStandardMaterial({ roughness: 0.4, metalness: 0.1 }),
    [],
  );

  const tempMatrix = useMemo(() => new Matrix4(), []);
  const tempColor = useMemo(() => new Color(), []);

  useEffect(() => {
    const mesh = meshRef.current;
    if (!mesh) return;

    for (let i = 0; i < atoms.length; i++) {
      const atom = atoms[i]!;
      const baseRadius = getAtomRadius(atom.symbol);
      const radius =
        renderStyle === "spacefill"
          ? baseRadius * 2.5
          : renderStyle === "stick"
            ? 0.08
            : baseRadius;
      const scaledRadius = radius * atomScale;

      tempMatrix.makeScale(scaledRadius, scaledRadius, scaledRadius);
      tempMatrix.setPosition(atom.position[0], atom.position[1], atom.position[2]);
      mesh.setMatrixAt(i, tempMatrix);

      tempColor.set(getCpkColor(atom.symbol));
      mesh.setColorAt(i, tempColor);
    }

    mesh.instanceMatrix.needsUpdate = true;
    if (mesh.instanceColor) mesh.instanceColor.needsUpdate = true;
  }, [atoms, renderStyle, atomScale, tempMatrix, tempColor]);

  if (atoms.length === 0) return null;

  return (
    <instancedMesh
      ref={meshRef}
      args={[geometry, material, atoms.length]}
      frustumCulled={false}
    />
  );
}
