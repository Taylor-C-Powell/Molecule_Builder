import { Line } from "@react-three/drei";
import { Vector3 } from "three";
import { useEditorStore } from "@/stores/editor-store";

/**
 * Renders a subtle grid and axis lines in the viewport
 * when showGrid is enabled.
 */
export function GridHelper() {
  const showGrid = useEditorStore((s) => s.showGrid);

  if (!showGrid) return null;

  return (
    <group>
      <gridHelper
        args={[40, 40, "#333333", "#1a1a1a"]}
        position={[0, -5, 0]}
      />
      {/* X axis - red */}
      <Line
        points={[new Vector3(-20, -5, 0), new Vector3(20, -5, 0)]}
        color="#ef4444"
        lineWidth={1}
        opacity={0.3}
        transparent
      />
      {/* Z axis - blue */}
      <Line
        points={[new Vector3(0, -5, -20), new Vector3(0, -5, 20)]}
        color="#3b82f6"
        lineWidth={1}
        opacity={0.3}
        transparent
      />
    </group>
  );
}
