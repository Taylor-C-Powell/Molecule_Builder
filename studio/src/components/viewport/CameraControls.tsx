import { useCallback } from "react";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { useViewportRefStore } from "@/stores/viewport-ref-store";
import { cn } from "@/lib/cn";

interface CameraPreset {
  label: string;
  shortLabel: string;
  position: [number, number, number];
}

const PRESETS: CameraPreset[] = [
  { label: "Front", shortLabel: "F", position: [0, 0, 15] },
  { label: "Back", shortLabel: "B", position: [0, 0, -15] },
  { label: "Top", shortLabel: "T", position: [0, 15, 0] },
  { label: "Right", shortLabel: "R", position: [15, 0, 0] },
  { label: "Left", shortLabel: "L", position: [-15, 0, 0] },
];

export function CameraControls() {
  const activeId = useWorkspaceStore((s) => s.activeId);
  const molecules = useWorkspaceStore((s) => s.molecules);
  const entry = activeId ? molecules.get(activeId) : undefined;
  const hasStructure = !!entry?.structure;

  const setCameraPosition = useCallback(
    (position: [number, number, number]) => {
      const gl = useViewportRefStore.getState().gl;
      if (!gl) return;
      // Access the camera through the gl renderer's domElement parent
      // We dispatch a custom event that the Canvas can pick up
      const event = new CustomEvent("molbuilder:camera", {
        detail: { position },
      });
      gl.domElement.dispatchEvent(event);
    },
    [],
  );

  if (!hasStructure) return null;

  return (
    <div className="absolute bottom-2 right-2 flex gap-0.5">
      {PRESETS.map((p) => (
        <button
          key={p.label}
          title={p.label}
          onClick={() => setCameraPosition(p.position)}
          className={cn(
            "w-6 h-6 flex items-center justify-center text-[9px] font-mono font-bold",
            "bg-bg/70 border border-border rounded-[var(--radius-sm)]",
            "text-text-muted hover:text-text-primary hover:bg-bg/90 transition-colors cursor-pointer",
          )}
        >
          {p.shortLabel}
        </button>
      ))}
    </div>
  );
}
