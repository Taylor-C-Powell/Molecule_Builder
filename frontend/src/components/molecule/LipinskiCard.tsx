import type { MoleculePropertiesResponse } from "@/api/types";
import { Card, CardHeader, CardTitle } from "@/components/ui/Card";
import { Badge } from "@/components/ui/Badge";

interface LipinskiCardProps {
  properties: MoleculePropertiesResponse;
}

const RULES = [
  { label: "MW <= 500", key: "molecular_weight" as const, test: (v: number) => v <= 500 },
  { label: "LogP <= 5", key: "logp" as const, test: (v: number | null) => v != null && v <= 5 },
  { label: "HBD <= 5", key: "hbd" as const, test: (v: number | null) => v != null && v <= 5 },
  { label: "HBA <= 10", key: "hba" as const, test: (v: number | null) => v != null && v <= 10 },
];

export function LipinskiCard({ properties }: LipinskiCardProps) {
  if (properties.lipinski_pass == null) return null;

  return (
    <Card>
      <CardHeader>
        <div className="flex items-center justify-between">
          <CardTitle>Lipinski's Rule of Five</CardTitle>
          <Badge variant={properties.lipinski_pass ? "success" : "danger"}>
            {properties.lipinski_pass ? "Pass" : `${properties.lipinski_violations} violation${properties.lipinski_violations !== 1 ? "s" : ""}`}
          </Badge>
        </div>
      </CardHeader>
      <div className="grid grid-cols-2 gap-2">
        {RULES.map((rule) => {
          const val = properties[rule.key];
          const pass = rule.test(val as never);
          return (
            <div
              key={rule.label}
              className="flex items-center gap-2 text-sm py-1"
            >
              <span className={pass ? "text-green" : "text-red"}>
                {pass ? "\u2713" : "\u2717"}
              </span>
              <span className="text-text-secondary">{rule.label}</span>
            </div>
          );
        })}
      </div>
    </Card>
  );
}
