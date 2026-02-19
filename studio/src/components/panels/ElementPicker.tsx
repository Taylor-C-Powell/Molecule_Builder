import { useEditorStore } from "@/stores/editor-store";
import { PLACEABLE_ELEMENTS, type PlaceableElement } from "@/types/editor";
import { getCpkColor } from "@/lib/cpk-colors";
import { cn } from "@/lib/cn";

const ELEMENT_LABELS: Record<PlaceableElement, string> = {
  C: "Carbon",
  N: "Nitrogen",
  O: "Oxygen",
  H: "Hydrogen",
  S: "Sulfur",
  P: "Phosphorus",
  F: "Fluorine",
  Cl: "Chlorine",
  Br: "Bromine",
  I: "Iodine",
};

export function ElementPicker() {
  const activeElement = useEditorStore((s) => s.activeElement);
  const setActiveElement = useEditorStore((s) => s.setActiveElement);

  return (
    <div className="grid grid-cols-5 gap-1">
      {PLACEABLE_ELEMENTS.map((el) => (
        <button
          key={el}
          title={ELEMENT_LABELS[el]}
          className={cn(
            "w-10 h-10 flex items-center justify-center rounded-[var(--radius-sm)] text-xs font-bold transition-all cursor-pointer border",
            activeElement === el
              ? "border-accent ring-1 ring-accent scale-105"
              : "border-transparent hover:border-white/20 hover:bg-white/5",
          )}
          style={{
            color: getCpkColor(el),
          }}
          onClick={() => setActiveElement(el)}
        >
          {el}
        </button>
      ))}
    </div>
  );
}
