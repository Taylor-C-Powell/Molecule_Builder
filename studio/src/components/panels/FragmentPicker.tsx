import { FRAGMENTS, type FragmentDef } from "@/lib/fragments";
import { useEditorStore } from "@/stores/editor-store";
import { cn } from "@/lib/cn";

interface FragmentPickerProps {
  onSelect: (fragment: FragmentDef) => void;
}

export function FragmentPicker({ onSelect }: FragmentPickerProps) {
  const activeFragment = useEditorStore((s) => s.activeFragment);

  return (
    <div className="grid grid-cols-5 gap-0.5 mt-1">
      {FRAGMENTS.map((f) => (
        <button
          key={f.shortLabel}
          title={f.name}
          onClick={() => onSelect(f)}
          className={cn(
            "px-1 py-1 text-[9px] font-mono border rounded-[var(--radius-sm)] transition-colors cursor-pointer",
            activeFragment?.shortLabel === f.shortLabel
              ? "bg-accent/20 border-accent text-accent"
              : "bg-bg border-border text-text-muted hover:text-text-secondary hover:border-border-hover",
          )}
        >
          {f.shortLabel}
        </button>
      ))}
    </div>
  );
}
