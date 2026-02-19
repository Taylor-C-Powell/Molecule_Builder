import { useThree } from "@react-three/fiber";
import { useEffect } from "react";
import { useViewportRefStore } from "@/stores/viewport-ref-store";

/**
 * R3F component that captures the WebGL renderer reference
 * so it can be used outside the Canvas (e.g., for screenshot export).
 */
export function ViewportGlCapture() {
  const { gl } = useThree();
  const setGl = useViewportRefStore((s) => s.setGl);

  useEffect(() => {
    setGl(gl);
  }, [gl, setGl]);

  return null;
}
