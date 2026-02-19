import { Line, Text } from "@react-three/drei";
import { Vector3 } from "three";
import { useEditorStore } from "@/stores/editor-store";
import { useWorkspaceStore } from "@/stores/workspace-store";

/**
 * Renders measurement lines and labels in 3D space.
 * Shows active measurements (completed) and the in-progress measure atom chain.
 */
export function MeasurementOverlay() {
  const measurements = useEditorStore((s) => s.measurements);
  const measureAtoms = useEditorStore((s) => s.measureAtoms);
  const activeId = useWorkspaceStore((s) => s.activeId);
  const molecules = useWorkspaceStore((s) => s.molecules);
  const entry = activeId ? molecules.get(activeId) : undefined;
  const structure = entry?.structure;

  if (!structure) return null;

  return (
    <group>
      {/* Completed measurements */}
      {measurements.map((m, i) => {
        const positions = m.atomIndices
          .map((idx) => structure.atoms[idx])
          .filter(Boolean)
          .map((a) => new Vector3(...a!.position));

        if (positions.length < 2) return null;

        // Midpoint for label placement
        const mid = new Vector3();
        for (const p of positions) mid.add(p);
        mid.divideScalar(positions.length);
        mid.y += 0.4;

        const label =
          m.type === "distance"
            ? `${m.value.toFixed(2)} ${m.unit}`
            : `${m.value.toFixed(1)}${m.unit}`;

        return (
          <group key={i}>
            <Line
              points={positions}
              color="#f59e0b"
              lineWidth={1.5}
              dashed
              dashSize={0.15}
              gapSize={0.08}
            />
            <Text
              position={[mid.x, mid.y, mid.z]}
              fontSize={0.25}
              color="#f59e0b"
              anchorX="center"
              anchorY="bottom"
              outlineWidth={0.02}
              outlineColor="#000000"
            >
              {label}
            </Text>
          </group>
        );
      })}

      {/* In-progress measurement chain */}
      {measureAtoms.length > 0 && (() => {
        const positions = measureAtoms
          .map((idx) => structure.atoms[idx])
          .filter(Boolean)
          .map((a) => new Vector3(...a!.position));

        if (positions.length < 2) return null;

        return (
          <Line
            points={positions}
            color="#60a5fa"
            lineWidth={1}
            dashed
            dashSize={0.1}
            gapSize={0.1}
          />
        );
      })()}
    </group>
  );
}
