import { useRef } from "react";
import { useThree, type ThreeEvent } from "@react-three/fiber";
import type { Mesh } from "three";
import { Vector3 } from "three";
import { useEditorStore } from "@/stores/editor-store";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { useHistoryStore } from "@/stores/history-store";
import { Tool } from "@/types/editor";

/**
 * Invisible plane that captures clicks for the ADD_ATOM tool.
 * Positioned at z=0 in molecule space, always facing the camera.
 * When clicked, adds a new atom at the intersection point.
 */
export function AddAtomPlane() {
  const planeRef = useRef<Mesh>(null);
  const { camera } = useThree();
  const tool = useEditorStore((s) => s.tool);
  const activeElement = useEditorStore((s) => s.activeElement);

  if (tool !== Tool.ADD_ATOM) return null;

  function handleClick(e: ThreeEvent<MouseEvent>) {
    e.stopPropagation();

    const ws = useWorkspaceStore.getState();
    const entry = ws.getActive();
    if (!entry?.structure) return;

    // Get click point in world space
    const point = e.point;

    // Snap to a reasonable grid (0.1 A resolution)
    const pos: [number, number, number] = [
      Math.round(point.x * 10) / 10,
      Math.round(point.y * 10) / 10,
      Math.round(point.z * 10) / 10,
    ];

    // Check if too close to an existing atom (< 0.5 A) -- reject if so
    const tooClose = entry.structure.atoms.some((a) => {
      const d = new Vector3(...a.position).distanceTo(new Vector3(...pos));
      return d < 0.5;
    });
    if (tooClose) return;

    useHistoryStore.getState().pushSnapshot(entry.structure, `Add ${activeElement}`);
    ws.addAtom(activeElement, pos);
  }

  // Large invisible plane perpendicular to camera
  const dir = new Vector3();
  camera.getWorldDirection(dir);

  return (
    <mesh
      ref={planeRef}
      position={[0, 0, 0]}
      onClick={handleClick}
      visible={false}
    >
      <planeGeometry args={[200, 200]} />
      <meshBasicMaterial transparent opacity={0} depthWrite={false} />
    </mesh>
  );
}
