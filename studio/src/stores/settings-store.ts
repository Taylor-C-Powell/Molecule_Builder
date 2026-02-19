import { create } from "zustand";

interface SettingsState {
  bgColor: string;
  atomScale: number;
  bondThickness: number;
  labelSize: number;
  antialias: boolean;

  setBgColor: (color: string) => void;
  setAtomScale: (scale: number) => void;
  setBondThickness: (thickness: number) => void;
  setLabelSize: (size: number) => void;
  setAntialias: (aa: boolean) => void;
}

const STORAGE_KEY = "molbuilder_studio_settings";

function loadSettings(): Partial<SettingsState> {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (raw) return JSON.parse(raw) as Partial<SettingsState>;
  } catch {
    // ignore
  }
  return {};
}

function saveSettings(state: SettingsState) {
  try {
    localStorage.setItem(
      STORAGE_KEY,
      JSON.stringify({
        bgColor: state.bgColor,
        atomScale: state.atomScale,
        bondThickness: state.bondThickness,
        labelSize: state.labelSize,
        antialias: state.antialias,
      }),
    );
  } catch {
    // ignore
  }
}

const defaults = {
  bgColor: "#0a0a0a",
  atomScale: 1.0,
  bondThickness: 1.0,
  labelSize: 0.3,
  antialias: true,
};

const saved = loadSettings();

export const useSettingsStore = create<SettingsState>((set, get) => ({
  bgColor: saved.bgColor ?? defaults.bgColor,
  atomScale: saved.atomScale ?? defaults.atomScale,
  bondThickness: saved.bondThickness ?? defaults.bondThickness,
  labelSize: saved.labelSize ?? defaults.labelSize,
  antialias: saved.antialias ?? defaults.antialias,

  setBgColor(color) {
    set({ bgColor: color });
    saveSettings(get());
  },
  setAtomScale(scale) {
    set({ atomScale: scale });
    saveSettings(get());
  },
  setBondThickness(thickness) {
    set({ bondThickness: thickness });
    saveSettings(get());
  },
  setLabelSize(size) {
    set({ labelSize: size });
    saveSettings(get());
  },
  setAntialias(aa) {
    set({ antialias: aa });
    saveSettings(get());
  },
}));
