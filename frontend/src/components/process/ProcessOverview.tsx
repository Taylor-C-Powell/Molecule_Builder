import type { ProcessEvaluateResponse } from "@/api/types";
import { Card, CardHeader, CardTitle } from "@/components/ui/Card";
import { Badge } from "@/components/ui/Badge";
import { pct, num } from "@/lib/format";

interface ProcessOverviewProps {
  data: ProcessEvaluateResponse;
}

export function ProcessOverview({ data }: ProcessOverviewProps) {
  return (
    <Card>
      <CardHeader>
        <div className="flex items-center justify-between">
          <CardTitle>Process Overview</CardTitle>
          <Badge variant={data.route_found ? "success" : "danger"}>
            {data.route_found ? "Route Found" : "No Route"}
          </Badge>
        </div>
      </CardHeader>
      <div className="grid grid-cols-2 sm:grid-cols-4 gap-4 text-sm">
        <div>
          <span className="text-text-muted block">Steps</span>
          <span className="font-bold text-lg">{data.total_steps}</span>
        </div>
        <div>
          <span className="text-text-muted block">Overall Yield</span>
          <span className="font-bold text-lg">{pct(data.overall_yield)}</span>
        </div>
        <div>
          <span className="text-text-muted block">Scale</span>
          <span className="font-bold text-lg">{num(data.scale_kg)} kg</span>
        </div>
        <div>
          <span className="text-text-muted block">SMILES</span>
          <code className="text-xs font-mono text-text-secondary break-all">{data.smiles}</code>
        </div>
      </div>
    </Card>
  );
}
