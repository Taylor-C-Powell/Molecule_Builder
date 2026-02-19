import { useEffect, useRef } from "react";
import { cn } from "@/lib/cn";

export interface ContextMenuItem {
  label: string;
  action: () => void;
  separator?: boolean;
  disabled?: boolean;
}

interface ContextMenuProps {
  x: number;
  y: number;
  items: ContextMenuItem[];
  onClose: () => void;
}

export function ContextMenu({ x, y, items, onClose }: ContextMenuProps) {
  const ref = useRef<HTMLDivElement>(null);

  useEffect(() => {
    function handleClick(e: MouseEvent) {
      if (ref.current && !ref.current.contains(e.target as Node)) {
        onClose();
      }
    }
    function handleKey(e: KeyboardEvent) {
      if (e.key === "Escape") onClose();
    }
    document.addEventListener("mousedown", handleClick);
    document.addEventListener("keydown", handleKey);
    return () => {
      document.removeEventListener("mousedown", handleClick);
      document.removeEventListener("keydown", handleKey);
    };
  }, [onClose]);

  // Keep menu within viewport bounds
  const style: React.CSSProperties = {
    position: "fixed",
    left: x,
    top: y,
    zIndex: 300,
  };

  return (
    <div ref={ref} style={style}>
      <div className="min-w-[160px] bg-bg-card border border-border rounded-[var(--radius-sm)] py-1 shadow-lg">
        {items.map((item, i) =>
          item.separator ? (
            <div key={i} className="h-px bg-border my-1" />
          ) : (
            <button
              key={item.label}
              disabled={item.disabled}
              className={cn(
                "w-full text-left px-3 py-1.5 text-xs transition-colors cursor-pointer",
                item.disabled
                  ? "text-text-muted/40 cursor-not-allowed"
                  : "text-text-secondary hover:text-text-primary hover:bg-white/5",
              )}
              onClick={() => {
                item.action();
                onClose();
              }}
            >
              {item.label}
            </button>
          ),
        )}
      </div>
    </div>
  );
}
