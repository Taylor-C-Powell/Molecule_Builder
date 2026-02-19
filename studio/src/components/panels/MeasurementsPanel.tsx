import { useEditorStore } from "@/stores/editor-store";
import { PanelContainer } from "./PanelContainer";

export function MeasurementsPanel() {
  const measurements = useEditorStore((s) => s.measurements);
  const clearMeasurements = useEditorStore((s) => s.clearMeasurements);

  if (measurements.length === 0) return null;

  return (
    <PanelContainer title="Measurements" defaultOpen>
      <div className="space-y-1">
        {measurements.map((m, i) => (
          <div
            key={i}
            className="flex items-center justify-between text-xs font-mono px-2 py-1.5 bg-white/3 rounded-[var(--radius-sm)]"
          >
            <span className="text-text-muted capitalize">{m.type}</span>
            <span className="text-text-secondary">
              [{m.atomIndices.join(", ")}]
            </span>
            <span className="text-amber-400 font-bold">
              {m.type === "distance"
                ? `${m.value.toFixed(3)} ${m.unit}`
                : `${m.value.toFixed(1)}${m.unit}`}
            </span>
          </div>
        ))}
      </div>
      <button
        className="mt-2 w-full py-1.5 text-xs text-text-muted hover:text-text-primary bg-white/5 hover:bg-white/10 rounded-[var(--radius-sm)] transition-colors cursor-pointer"
        onClick={clearMeasurements}
      >
        Clear All
      </button>
    </PanelContainer>
  );
}
