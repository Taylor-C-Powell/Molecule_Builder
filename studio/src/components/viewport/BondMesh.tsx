import { useMemo, useState, memo } from "react";
import { type ThreeEvent } from "@react-three/fiber";
import type { AtomResponse, BondResponse } from "@/api/types";
import { getCpkColor } from "@/lib/cpk-colors";
import { computeBondTransform, computeBondOffset } from "@/lib/bond-geometry";
import { useEditorStore } from "@/stores/editor-store";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { useHistoryStore } from "@/stores/history-store";
import { Tool } from "@/types/editor";
import { Quaternion } from "three";

interface BondMeshProps {
  bond: BondResponse;
  atoms: AtomResponse[];
  renderStyle: "ball-and-stick" | "stick" | "spacefill";
}

const BOND_RADIUS = 0.06;
const BOND_SEGMENTS = 8;
const DOUBLE_OFFSET = 0.1;
const TRIPLE_OFFSET = 0.14;

export const BondMesh = memo(function BondMesh({ bond, atoms, renderStyle }: BondMeshProps) {
  const atomA = atoms[bond.atom_i];
  const atomB = atoms[bond.atom_j];
  const [hovered, setHovered] = useState(false);
  const tool = useEditorStore((s) => s.tool);
  const selectedBonds = useEditorStore((s) => s.selectedBonds);

  const bondKey = `${bond.atom_i}-${bond.atom_j}`;
  const isSelected = selectedBonds.has(bondKey) || selectedBonds.has(`${bond.atom_j}-${bond.atom_i}`);
  const isEraseHover = tool === Tool.ERASE && hovered;

  const transforms = useMemo(() => {
    if (!atomA || !atomB) return [];

    const posA = atomA.position;
    const posB = atomB.position;
    const order = Math.round(bond.order);

    if (order === 1) {
      return [computeBondTransform(posA, posB)];
    }

    if (order === 2) {
      const offset = computeBondOffset(posA, posB, DOUBLE_OFFSET);
      const posA1: [number, number, number] = [posA[0] + offset[0], posA[1] + offset[1], posA[2] + offset[2]];
      const posB1: [number, number, number] = [posB[0] + offset[0], posB[1] + offset[1], posB[2] + offset[2]];
      const posA2: [number, number, number] = [posA[0] - offset[0], posA[1] - offset[1], posA[2] - offset[2]];
      const posB2: [number, number, number] = [posB[0] - offset[0], posB[1] - offset[1], posB[2] - offset[2]];
      return [
        computeBondTransform(posA1, posB1),
        computeBondTransform(posA2, posB2),
      ];
    }

    if (order >= 3) {
      const offset = computeBondOffset(posA, posB, TRIPLE_OFFSET);
      const posA1: [number, number, number] = [posA[0] + offset[0], posA[1] + offset[1], posA[2] + offset[2]];
      const posB1: [number, number, number] = [posB[0] + offset[0], posB[1] + offset[1], posB[2] + offset[2]];
      const posA2: [number, number, number] = [posA[0] - offset[0], posA[1] - offset[1], posA[2] - offset[2]];
      const posB2: [number, number, number] = [posB[0] - offset[0], posB[1] - offset[1], posB[2] - offset[2]];
      return [
        computeBondTransform(posA1, posB1),
        computeBondTransform(posA, posB),
        computeBondTransform(posA2, posB2),
      ];
    }

    return [computeBondTransform(posA, posB)];
  }, [atomA, atomB, bond.order]);

  if (!atomA || !atomB || renderStyle === "spacefill") return null;

  const colorA = getCpkColor(atomA.symbol);
  const colorB = getCpkColor(atomB.symbol);
  const radius = renderStyle === "stick" ? 0.08 : BOND_RADIUS;

  function handleClick(e: ThreeEvent<MouseEvent>) {
    e.stopPropagation();

    switch (tool) {
      case Tool.SELECT:
      case Tool.ROTATE:
        useEditorStore.getState().selectBond(bondKey, e.shiftKey);
        break;

      case Tool.ADD_BOND: {
        // Click existing bond in ADD_BOND mode -> cycle order
        const ws = useWorkspaceStore.getState();
        const entry = ws.getActive();
        if (entry?.structure) {
          useHistoryStore.getState().pushSnapshot(entry.structure, "Cycle bond order");
          ws.cycleBondOrder(bond.atom_i, bond.atom_j);
        }
        break;
      }

      case Tool.ERASE: {
        const ws = useWorkspaceStore.getState();
        const entry = ws.getActive();
        if (entry?.structure) {
          useHistoryStore.getState().pushSnapshot(entry.structure, "Delete bond");
          ws.removeBonds(new Set([bondKey]));
        }
        break;
      }

      default:
        break;
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
        "bond",
        bond.atom_i,
        bondKey,
      );
    }
  }

  function handlePointerOver(e: ThreeEvent<PointerEvent>) {
    e.stopPropagation();
    setHovered(true);
    useEditorStore.getState().setHoveredBond(bondKey);
    if (tool === Tool.ERASE) document.body.style.cursor = "crosshair";
  }

  function handlePointerOut() {
    setHovered(false);
    useEditorStore.getState().setHoveredBond(null);
    if (tool === Tool.ERASE) document.body.style.cursor = "auto";
  }

  // Determine highlight
  let emissive = "#000000";
  let emissiveIntensity = 0;
  if (isEraseHover) {
    emissive = "#ef4444";
    emissiveIntensity = 0.6;
  } else if (isSelected) {
    emissive = "#3b82f6";
    emissiveIntensity = 0.4;
  } else if (hovered) {
    emissive = "#1e40af";
    emissiveIntensity = 0.15;
  }

  return (
    <group
      onClick={handleClick}
      onContextMenu={handleContextMenu}
      onPointerOver={handlePointerOver}
      onPointerOut={handlePointerOut}
    >
      {transforms.map((t, i) => {
        const q = new Quaternion(...t.quaternion);
        return (
          <group key={i}>
            <mesh
              position={[
                t.position[0] - (t.length / 4) * (atomB.position[0] - atomA.position[0]) / t.length,
                t.position[1] - (t.length / 4) * (atomB.position[1] - atomA.position[1]) / t.length,
                t.position[2] - (t.length / 4) * (atomB.position[2] - atomA.position[2]) / t.length,
              ]}
              quaternion={q}
            >
              <cylinderGeometry args={[radius, radius, t.length / 2, BOND_SEGMENTS]} />
              <meshStandardMaterial
                color={colorA}
                roughness={0.5}
                metalness={0.1}
                emissive={emissive}
                emissiveIntensity={emissiveIntensity}
              />
            </mesh>
            <mesh
              position={[
                t.position[0] + (t.length / 4) * (atomB.position[0] - atomA.position[0]) / t.length,
                t.position[1] + (t.length / 4) * (atomB.position[1] - atomA.position[1]) / t.length,
                t.position[2] + (t.length / 4) * (atomB.position[2] - atomA.position[2]) / t.length,
              ]}
              quaternion={q}
            >
              <cylinderGeometry args={[radius, radius, t.length / 2, BOND_SEGMENTS]} />
              <meshStandardMaterial
                color={colorB}
                roughness={0.5}
                metalness={0.1}
                emissive={emissive}
                emissiveIntensity={emissiveIntensity}
              />
            </mesh>
          </group>
        );
      })}
    </group>
  );
});
