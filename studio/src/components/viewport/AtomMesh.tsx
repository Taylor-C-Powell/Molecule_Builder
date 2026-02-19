import { useRef, useState, memo } from "react";
import { type ThreeEvent } from "@react-three/fiber";
import { Text } from "@react-three/drei";
import type { Mesh } from "three";
import type { AtomResponse } from "@/api/types";
import { getCpkColor, getAtomRadius } from "@/lib/cpk-colors";
import { useEditorStore } from "@/stores/editor-store";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { useHistoryStore } from "@/stores/history-store";
import { useAtomDrag } from "@/hooks/useAtomDrag";
import { Tool } from "@/types/editor";
import { buildMeasurement } from "@/lib/measurement";

interface AtomMeshProps {
  atom: AtomResponse;
  selected: boolean;
  showLabel: boolean;
  renderStyle: "ball-and-stick" | "stick" | "spacefill";
  centroid: readonly [number, number, number];
}

export const AtomMesh = memo(function AtomMesh({ atom, selected, showLabel, renderStyle, centroid }: AtomMeshProps) {
  const meshRef = useRef<Mesh>(null);
  const [hovered, setHovered] = useState(false);
  const { startDrag } = useAtomDrag(atom.index, centroid);
  const tool = useEditorStore((s) => s.tool);
  const selectAtom = useEditorStore((s) => s.selectAtom);
  const setHoveredAtom = useEditorStore((s) => s.setHoveredAtom);

  const color = getCpkColor(atom.symbol);
  const baseRadius = getAtomRadius(atom.symbol);
  const radius =
    renderStyle === "spacefill"
      ? baseRadius * 2.5
      : renderStyle === "stick"
        ? 0.08
        : baseRadius;
  const scale = hovered ? 1.15 : 1.0;

  // Erase tool: red tint on hover
  const isEraseHover = tool === Tool.ERASE && hovered;
  // Bond tool: highlight start atom, tint targets
  const bondStartAtom = useEditorStore((s) => s.bondStartAtom);
  const isBondStart = tool === Tool.ADD_BOND && bondStartAtom === atom.index;
  const isBondTarget = tool === Tool.ADD_BOND && bondStartAtom != null && bondStartAtom !== atom.index && hovered;

  function handleClick(e: ThreeEvent<MouseEvent>) {
    e.stopPropagation();

    switch (tool) {
      case Tool.SELECT:
      case Tool.ROTATE:
        selectAtom(atom.index, e.shiftKey);
        break;

      case Tool.ADD_ATOM: {
        // Click existing atom in ADD_ATOM mode -> replace element
        const ws = useWorkspaceStore.getState();
        const entry = ws.getActive();
        const activeEl = useEditorStore.getState().activeElement;
        if (entry?.structure && atom.symbol !== activeEl) {
          useHistoryStore.getState().pushSnapshot(entry.structure, `Replace ${atom.symbol} -> ${activeEl}`);
          ws.replaceAtomElement(atom.index, activeEl);
        }
        break;
      }

      case Tool.ADD_BOND: {
        const editor = useEditorStore.getState();
        if (editor.bondStartAtom == null) {
          editor.setBondStartAtom(atom.index);
        } else if (editor.bondStartAtom !== atom.index) {
          const ws = useWorkspaceStore.getState();
          const entry = ws.getActive();
          if (entry?.structure) {
            useHistoryStore.getState().pushSnapshot(entry.structure, "Add bond");
            ws.addBond(editor.bondStartAtom, atom.index, editor.activeBondOrder);
          }
          editor.setBondStartAtom(null);
        }
        break;
      }

      case Tool.ERASE: {
        const ws = useWorkspaceStore.getState();
        const entry = ws.getActive();
        if (entry?.structure) {
          useHistoryStore.getState().pushSnapshot(entry.structure, `Delete ${atom.symbol}${atom.index}`);
          ws.removeAtoms(new Set([atom.index]));
        }
        break;
      }

      case Tool.MEASURE: {
        const editor = useEditorStore.getState();
        editor.pushMeasureAtom(atom.index);
        const chain = useEditorStore.getState().measureAtoms;
        if (chain.length >= 2 && chain.length <= 4) {
          const ws = useWorkspaceStore.getState();
          const entry = ws.getActive();
          if (entry?.structure) {
            const m = buildMeasurement(entry.structure.atoms, chain);
            if (m) {
              editor.addMeasurement(m);
              if (chain.length === 4) {
                editor.clearMeasureAtoms();
              }
            }
          }
        }
        break;
      }

      default:
        break;
    }
  }

  function handlePointerDown(e: ThreeEvent<PointerEvent>) {
    if (e.button !== 0) return; // left click only
    const nativeEvent = e.nativeEvent;
    if (nativeEvent) {
      startDrag(nativeEvent);
    }
  }

  function handleContextMenu(e: ThreeEvent<MouseEvent>) {
    e.stopPropagation();
    const nativeEvent = e.nativeEvent ?? (e as unknown as { originalEvent?: MouseEvent }).originalEvent;
    if (nativeEvent) {
      nativeEvent.preventDefault();
      useEditorStore.getState().openContextMenu(
        nativeEvent.clientX,
        nativeEvent.clientY,
        "atom",
        atom.index,
      );
    }
  }

  function handlePointerOver(e: ThreeEvent<PointerEvent>) {
    e.stopPropagation();
    setHovered(true);
    setHoveredAtom(atom.index);
    document.body.style.cursor = tool === Tool.ERASE ? "crosshair" : "pointer";
  }

  function handlePointerOut() {
    setHovered(false);
    setHoveredAtom(null);
    document.body.style.cursor = "auto";
  }

  // Compute emissive color based on state
  let emissive = "#000000";
  let emissiveIntensity = 0;
  if (isEraseHover) {
    emissive = "#ef4444";
    emissiveIntensity = 0.6;
  } else if (isBondStart) {
    emissive = "#22c55e";
    emissiveIntensity = 0.5;
  } else if (isBondTarget) {
    emissive = "#3b82f6";
    emissiveIntensity = 0.4;
  } else if (selected) {
    emissive = "#3b82f6";
    emissiveIntensity = 0.5;
  } else if (hovered) {
    emissive = "#1e40af";
    emissiveIntensity = 0.2;
  }

  return (
    <group position={atom.position}>
      <mesh
        ref={meshRef}
        scale={scale}
        onClick={handleClick}
        onPointerDown={handlePointerDown}
        onContextMenu={handleContextMenu}
        onPointerOver={handlePointerOver}
        onPointerOut={handlePointerOut}
      >
        <sphereGeometry args={[radius, 32, 32]} />
        <meshStandardMaterial
          color={color}
          emissive={emissive}
          emissiveIntensity={emissiveIntensity}
          roughness={0.4}
          metalness={0.1}
        />
      </mesh>
      {showLabel && (
        <Text
          position={[0, radius + 0.3, 0]}
          fontSize={0.3}
          color="#ededed"
          anchorX="center"
          anchorY="bottom"
        >
          {atom.symbol}
          {atom.index}
        </Text>
      )}
    </group>
  );
});
