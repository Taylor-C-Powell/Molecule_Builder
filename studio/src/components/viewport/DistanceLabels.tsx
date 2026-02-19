import { Text } from "@react-three/drei";
import { Vector3 } from "three";
import { useEditorStore } from "@/stores/editor-store";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { useSettingsStore } from "@/stores/settings-store";

/**
 * Renders distance labels at the midpoint of each bond in 3D space.
 * Only visible when showDistanceLabels is enabled in the editor store.
 */
export function DistanceLabels() {
  const showDistanceLabels = useEditorStore((s) => s.showDistanceLabels);
  const showHydrogens = useEditorStore((s) => s.showHydrogens);
  const labelSize = useSettingsStore((s) => s.labelSize);
  const activeId = useWorkspaceStore((s) => s.activeId);
  const molecules = useWorkspaceStore((s) => s.molecules);
  const entry = activeId ? molecules.get(activeId) : undefined;
  const structure = entry?.structure;

  if (!showDistanceLabels || !structure) return null;

  return (
    <group>
      {structure.bonds.map((bond) => {
        const atomA = structure.atoms[bond.atom_i];
        const atomB = structure.atoms[bond.atom_j];
        if (!atomA || !atomB) return null;

        // Skip hydrogen bonds if hydrogens are hidden
        if (!showHydrogens && (atomA.symbol === "H" || atomB.symbol === "H")) return null;

        const posA = new Vector3(...atomA.position);
        const posB = new Vector3(...atomB.position);
        const mid = posA.clone().add(posB).multiplyScalar(0.5);
        const dist = posA.distanceTo(posB);

        // Offset slightly above the bond axis
        mid.y += 0.15;

        return (
          <Text
            key={`${bond.atom_i}-${bond.atom_j}`}
            position={[mid.x, mid.y, mid.z]}
            fontSize={labelSize * 0.6}
            color="#94a3b8"
            anchorX="center"
            anchorY="bottom"
            outlineWidth={0.01}
            outlineColor="#000000"
          >
            {dist.toFixed(2)}
          </Text>
        );
      })}
    </group>
  );
}
