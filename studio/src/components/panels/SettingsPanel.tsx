import { PanelContainer } from "./PanelContainer";
import { useSettingsStore } from "@/stores/settings-store";

const BG_PRESETS = [
  { label: "Dark", value: "#0a0a0a" },
  { label: "Charcoal", value: "#1a1a2e" },
  { label: "Midnight", value: "#0d1117" },
  { label: "Navy", value: "#0a192f" },
  { label: "White", value: "#f5f5f5" },
];

export function SettingsPanel() {
  const bgColor = useSettingsStore((s) => s.bgColor);
  const atomScale = useSettingsStore((s) => s.atomScale);
  const bondThickness = useSettingsStore((s) => s.bondThickness);
  const labelSize = useSettingsStore((s) => s.labelSize);
  const setBgColor = useSettingsStore((s) => s.setBgColor);
  const setAtomScale = useSettingsStore((s) => s.setAtomScale);
  const setBondThickness = useSettingsStore((s) => s.setBondThickness);
  const setLabelSize = useSettingsStore((s) => s.setLabelSize);

  return (
    <PanelContainer title="Settings" defaultOpen={false}>
      <div className="space-y-3">
        {/* Background color */}
        <div>
          <label className="text-[10px] text-text-muted block mb-1">
            Background Color
          </label>
          <div className="flex gap-1 flex-wrap">
            {BG_PRESETS.map((p) => (
              <button
                key={p.value}
                className="w-6 h-6 rounded-[var(--radius-sm)] border cursor-pointer transition-transform hover:scale-110"
                style={{
                  backgroundColor: p.value,
                  borderColor: bgColor === p.value ? "var(--color-accent)" : "var(--color-border)",
                  borderWidth: bgColor === p.value ? 2 : 1,
                }}
                title={p.label}
                onClick={() => setBgColor(p.value)}
              />
            ))}
            <input
              type="color"
              value={bgColor}
              onChange={(e) => setBgColor(e.target.value)}
              className="w-6 h-6 rounded-[var(--radius-sm)] border border-border cursor-pointer bg-transparent"
              title="Custom color"
            />
          </div>
        </div>

        {/* Atom scale */}
        <div>
          <label className="text-[10px] text-text-muted flex items-center justify-between mb-1">
            <span>Atom Size</span>
            <span className="font-mono">{atomScale.toFixed(1)}x</span>
          </label>
          <input
            type="range"
            min={0.3}
            max={3.0}
            step={0.1}
            value={atomScale}
            onChange={(e) => setAtomScale(parseFloat(e.target.value))}
            className="w-full h-1 accent-accent"
          />
        </div>

        {/* Bond thickness */}
        <div>
          <label className="text-[10px] text-text-muted flex items-center justify-between mb-1">
            <span>Bond Thickness</span>
            <span className="font-mono">{bondThickness.toFixed(1)}x</span>
          </label>
          <input
            type="range"
            min={0.3}
            max={3.0}
            step={0.1}
            value={bondThickness}
            onChange={(e) => setBondThickness(parseFloat(e.target.value))}
            className="w-full h-1 accent-accent"
          />
        </div>

        {/* Label size */}
        <div>
          <label className="text-[10px] text-text-muted flex items-center justify-between mb-1">
            <span>Label Size</span>
            <span className="font-mono">{labelSize.toFixed(1)}</span>
          </label>
          <input
            type="range"
            min={0.1}
            max={0.8}
            step={0.05}
            value={labelSize}
            onChange={(e) => setLabelSize(parseFloat(e.target.value))}
            className="w-full h-1 accent-accent"
          />
        </div>
      </div>
    </PanelContainer>
  );
}
