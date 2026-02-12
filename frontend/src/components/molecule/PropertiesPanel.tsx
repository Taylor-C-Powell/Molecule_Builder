import type { MoleculePropertiesResponse } from "@/api/types";
import { Card, CardHeader, CardTitle } from "@/components/ui/Card";
import { Skeleton } from "@/components/ui/Skeleton";
import { num } from "@/lib/format";
import { LipinskiCard } from "./LipinskiCard";
import { FunctionalGroups } from "./FunctionalGroups";

interface PropertiesPanelProps {
  properties?: MoleculePropertiesResponse;
  loading?: boolean;
}

export function PropertiesPanel({ properties, loading }: PropertiesPanelProps) {
  if (loading) {
    return (
      <div className="space-y-4">
        <Skeleton className="h-48" />
        <Skeleton className="h-32" />
      </div>
    );
  }

  if (!properties) return null;

  const rows: [string, string][] = [
    ["Formula", properties.formula],
    ["Molecular Weight", `${num(properties.molecular_weight)} g/mol`],
    ["Atoms", String(properties.num_atoms)],
    ["Bonds", String(properties.num_bonds)],
    ["Heavy Atoms", properties.heavy_atom_count != null ? String(properties.heavy_atom_count) : "--"],
    ["Rotatable Bonds", properties.rotatable_bonds != null ? String(properties.rotatable_bonds) : "--"],
    ["LogP", properties.logp != null ? num(properties.logp) : "--"],
    ["TPSA", properties.tpsa != null ? `${num(properties.tpsa)} A^2` : "--"],
    ["HBD", properties.hbd != null ? String(properties.hbd) : "--"],
    ["HBA", properties.hba != null ? String(properties.hba) : "--"],
  ];

  return (
    <div className="space-y-4">
      <Card>
        <CardHeader>
          <div className="flex items-center justify-between">
            <CardTitle>Properties</CardTitle>
            <code className="text-xs font-mono text-text-muted">{properties.smiles}</code>
          </div>
        </CardHeader>
        <div className="grid grid-cols-2 gap-x-6 gap-y-2">
          {rows.map(([label, value]) => (
            <div key={label} className="flex justify-between text-sm py-1 border-b border-border/50">
              <span className="text-text-secondary">{label}</span>
              <span className="font-medium font-mono">{value}</span>
            </div>
          ))}
        </div>
      </Card>

      <LipinskiCard properties={properties} />
      <FunctionalGroups groups={properties.functional_groups} />
    </div>
  );
}
