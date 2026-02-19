import { useCallback } from "react";
import { useViewportRefStore } from "@/stores/viewport-ref-store";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { useToastStore } from "@/stores/toast-store";

export function useScreenshot() {
  const addToast = useToastStore((s) => s.addToast);

  const captureScreenshot = useCallback(() => {
    const gl = useViewportRefStore.getState().gl;
    if (!gl) {
      addToast("Viewport not ready", "error");
      return;
    }

    const entry = useWorkspaceStore.getState().getActive();
    const name = entry?.name || "molecule";

    try {
      const dataUrl = gl.domElement.toDataURL("image/png");
      const a = document.createElement("a");
      a.href = dataUrl;
      a.download = `${name}.png`;
      a.click();
      addToast("Screenshot saved", "success");
    } catch {
      addToast("Failed to capture screenshot", "error");
    }
  }, [addToast]);

  return { captureScreenshot };
}
