import { create } from "zustand";
import type { WebGLRenderer } from "three";

interface ViewportRefState {
  gl: WebGLRenderer | null;
  setGl: (gl: WebGLRenderer) => void;
}

export const useViewportRefStore = create<ViewportRefState>((set) => ({
  gl: null,
  setGl(gl) {
    set({ gl });
  },
}));
