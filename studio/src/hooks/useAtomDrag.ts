import { useRef, useCallback } from "react";
import { useThree } from "@react-three/fiber";
import { Plane, Raycaster, Vector2, Vector3 } from "three";
import { useEditorStore } from "@/stores/editor-store";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { useHistoryStore } from "@/stores/history-store";
import type { Molecule3DResponse } from "@/api/types";
import { Tool } from "@/types/editor";

interface DragState {
  active: boolean;
  atomIndex: number;
  plane: Plane;
  offset: Vector3;
  snapshotPushed: boolean;
  preSnapshot: Molecule3DResponse | null;
}

const _raycaster = new Raycaster();
const _ndc = new Vector2();
const _target = new Vector3();

/**
 * Hook for dragging atoms in 3D space.
 * Only active when the SELECT tool is chosen.
 * Projects mouse movement onto a plane perpendicular to the camera
 * that passes through the atom's initial position.
 */
export function useAtomDrag(atomIndex: number, centroid: readonly [number, number, number]) {
  const { camera, gl } = useThree();
  const dragRef = useRef<DragState>({
    active: false,
    atomIndex: -1,
    plane: new Plane(),
    offset: new Vector3(),
    snapshotPushed: false,
    preSnapshot: null,
  });

  const rayPlaneHit = useCallback(
    (e: PointerEvent, plane: Plane): Vector3 | null => {
      const rect = gl.domElement.getBoundingClientRect();
      _ndc.set(
        ((e.clientX - rect.left) / rect.width) * 2 - 1,
        -((e.clientY - rect.top) / rect.height) * 2 + 1,
      );
      _raycaster.setFromCamera(_ndc, camera);
      return _raycaster.ray.intersectPlane(plane, _target.clone());
    },
    [camera, gl],
  );

  const onDrag = useCallback(
    (e: PointerEvent) => {
      const state = dragRef.current;
      if (!state.active) return;

      const hit = rayPlaneHit(e, state.plane);
      if (!hit) return;

      const newWorld = hit.sub(state.offset);
      const newPos: [number, number, number] = [
        newWorld.x + centroid[0],
        newWorld.y + centroid[1],
        newWorld.z + centroid[2],
      ];

      if (!state.snapshotPushed && state.preSnapshot) {
        useHistoryStore.getState().pushSnapshot(state.preSnapshot, "Move atom");
        state.snapshotPushed = true;
      }

      useWorkspaceStore.getState().moveAtom(state.atomIndex, newPos);
    },
    [rayPlaneHit, centroid],
  );

  const endDrag = useCallback(
    (e: PointerEvent) => {
      dragRef.current.active = false;
      gl.domElement.releasePointerCapture(e.pointerId);
      gl.domElement.removeEventListener("pointermove", onDrag);
      gl.domElement.removeEventListener("pointerup", endDrag);
    },
    [gl, onDrag],
  );

  const startDrag = useCallback(
    (e: PointerEvent) => {
      const tool = useEditorStore.getState().tool;
      if (tool !== Tool.SELECT) return;

      const entry = useWorkspaceStore.getState().getActive();
      if (!entry?.structure) return;
      const atom = entry.structure.atoms[atomIndex];
      if (!atom) return;

      const atomWorld = new Vector3(
        atom.position[0] - centroid[0],
        atom.position[1] - centroid[1],
        atom.position[2] - centroid[2],
      );
      const cameraDir = new Vector3();
      camera.getWorldDirection(cameraDir);

      const plane = new Plane().setFromNormalAndCoplanarPoint(cameraDir, atomWorld);
      const hit = rayPlaneHit(e, plane);
      if (!hit) return;

      const offset = hit.sub(atomWorld);

      dragRef.current = {
        active: true,
        atomIndex,
        plane,
        offset,
        snapshotPushed: false,
        preSnapshot: entry.structure,
      };

      gl.domElement.setPointerCapture(e.pointerId);
      gl.domElement.addEventListener("pointermove", onDrag);
      gl.domElement.addEventListener("pointerup", endDrag);
    },
    [atomIndex, camera, gl, centroid, rayPlaneHit, onDrag, endDrag],
  );

  return { startDrag };
}
