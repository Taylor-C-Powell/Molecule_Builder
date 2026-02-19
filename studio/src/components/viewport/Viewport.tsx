import { Canvas } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";
import { MoleculeScene } from "./MoleculeScene";
import { ViewportOverlay } from "./ViewportOverlay";
import { SelectionOutline } from "./SelectionOutline";
import { ViewportGlCapture } from "./ViewportGlCapture";
import { GridHelper } from "./GridHelper";
import { WelcomeScreen } from "./WelcomeScreen";
import { CameraControls } from "./CameraControls";
import { ErrorBoundary } from "@/components/ui/ErrorBoundary";
import { useEditorStore } from "@/stores/editor-store";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { useSettingsStore } from "@/stores/settings-store";

export function Viewport() {
  const clearSelection = useEditorStore((s) => s.clearSelection);
  const autoRotate = useEditorStore((s) => s.autoRotate);
  const activeId = useWorkspaceStore((s) => s.activeId);
  const molecules = useWorkspaceStore((s) => s.molecules);
  const entry = activeId ? molecules.get(activeId) : undefined;
  const hasStructure = !!entry?.structure;

  const bgColor = useSettingsStore((s) => s.bgColor);

  return (
    <div className="relative w-full h-full">
      <ErrorBoundary
        fallback={
          <div className="flex items-center justify-center h-full bg-bg">
            <div className="text-center space-y-2">
              <div className="text-sm text-red-400">WebGL error</div>
              <div className="text-xs text-text-muted">
                The 3D viewport encountered an error. Try refreshing.
              </div>
            </div>
          </div>
        }
      >
        <Canvas
          camera={{ position: [0, 0, 15], fov: 50 }}
          style={{ background: bgColor }}
          gl={{ preserveDrawingBuffer: true }}
          onPointerMissed={() => clearSelection()}
        >
          <ambientLight intensity={0.4} />
          <directionalLight position={[10, 10, 10]} intensity={0.8} />
          <directionalLight position={[-5, -5, -10]} intensity={0.3} />
          <MoleculeScene />
          <GridHelper />
          <OrbitControls
            enableDamping
            dampingFactor={0.1}
            autoRotate={autoRotate}
            autoRotateSpeed={2.0}
          />
          <ViewportGlCapture />
        </Canvas>
      </ErrorBoundary>

      <ViewportOverlay />
      <SelectionOutline />
      <CameraControls />

      {!hasStructure && <WelcomeScreen />}
    </div>
  );
}
