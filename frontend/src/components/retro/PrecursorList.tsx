import type { PrecursorResponse } from "@/api/types";
import { Badge } from "@/components/ui/Badge";
import { usd } from "@/lib/format";

interface PrecursorListProps {
  precursors: PrecursorResponse[];
}

export function PrecursorList({ precursors }: PrecursorListProps) {
  if (!precursors.length) return null;

  return (
    <div className="space-y-1">
      {precursors.map((p) => (
        <div key={p.smiles} className="flex items-center gap-2 text-sm">
          <Badge variant="accent">{p.name}</Badge>
          <code className="text-xs text-text-muted font-mono">{p.smiles}</code>
          <span className="ml-auto text-text-secondary">{usd(p.cost_per_kg)}/kg</span>
        </div>
      ))}
    </div>
  );
}
