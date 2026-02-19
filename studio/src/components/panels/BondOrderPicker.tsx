import { useEditorStore } from "@/stores/editor-store";
import type { BondOrder } from "@/types/editor";
import { cn } from "@/lib/cn";

const BOND_ORDERS: { order: BondOrder; label: string; icon: string }[] = [
  { order: 1, label: "Single", icon: "\u2014" },
  { order: 2, label: "Double", icon: "=" },
  { order: 3, label: "Triple", icon: "\u2261" },
];

export function BondOrderPicker() {
  const activeBondOrder = useEditorStore((s) => s.activeBondOrder);
  const setActiveBondOrder = useEditorStore((s) => s.setActiveBondOrder);

  return (
    <div className="flex gap-1">
      {BOND_ORDERS.map(({ order, label, icon }) => (
        <button
          key={order}
          title={`${label} bond`}
          className={cn(
            "flex-1 h-8 flex items-center justify-center rounded-[var(--radius-sm)] text-sm font-bold transition-colors cursor-pointer border",
            activeBondOrder === order
              ? "border-accent text-accent bg-accent/10"
              : "border-transparent text-text-muted hover:text-text-primary hover:bg-white/5",
          )}
          onClick={() => setActiveBondOrder(order)}
        >
          {icon}
        </button>
      ))}
    </div>
  );
}
