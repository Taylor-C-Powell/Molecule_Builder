import { useRef } from "react";
import { useFrame } from "@react-three/fiber";
import { Line } from "@react-three/drei";
import type { Group } from "three";
import { Vector3 } from "three";
import { useEditorStore } from "@/stores/editor-store";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { Tool } from "@/types/editor";

/**
 * Renders a dashed line from the bond-start atom to the hovered atom
 * while the ADD_BOND tool is active and a start atom is selected.
 */
export function BondPreview() {
  const groupRef = useRef<Group>(null);
  const tool = useEditorStore((s) => s.tool);
  const bondStartAtom = useEditorStore((s) => s.bondStartAtom);
  const hoveredAtom = useEditorStore((s) => s.hoveredAtom);

  const activeId = useWorkspaceStore((s) => s.activeId);
  const molecules = useWorkspaceStore((s) => s.molecules);
  const entry = activeId ? molecules.get(activeId) : undefined;
  const structure = entry?.structure;

  useFrame(() => {
    // Keep the group visible -- the Line component handles positioning
  });

  if (tool !== Tool.ADD_BOND || bondStartAtom == null || !structure) return null;

  const startAtom = structure.atoms[bondStartAtom];
  if (!startAtom) return null;

  // Show preview to hovered atom, or nothing if no hover target
  const endIndex = hoveredAtom;
  if (endIndex == null || endIndex === bondStartAtom) return null;
  const endAtom = structure.atoms[endIndex];
  if (!endAtom) return null;

  const startPos = new Vector3(...startAtom.position);
  const endPos = new Vector3(...endAtom.position);

  return (
    <group ref={groupRef}>
      <Line
        points={[startPos, endPos]}
        color="#3b82f6"
        lineWidth={2}
        dashed
        dashSize={0.2}
        gapSize={0.1}
      />
    </group>
  );
}
