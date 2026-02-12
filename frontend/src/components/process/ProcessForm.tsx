import { useState } from "react";
import { Button } from "@/components/ui/Button";
import { Input } from "@/components/ui/Input";
import { Alert } from "@/components/ui/Alert";

interface ProcessFormProps {
  onSubmit: (smiles: string, scaleKg: number, maxDepth: number, beamWidth: number) => void;
  loading?: boolean;
  error?: string | null;
}

export function ProcessForm({ onSubmit, loading, error }: ProcessFormProps) {
  const [smiles, setSmiles] = useState("");
  const [scaleKg, setScaleKg] = useState(1);
  const [maxDepth, setMaxDepth] = useState(5);
  const [beamWidth, setBeamWidth] = useState(5);

  function handleSubmit(e: React.FormEvent) {
    e.preventDefault();
    if (smiles.trim()) onSubmit(smiles.trim(), scaleKg, maxDepth, beamWidth);
  }

  return (
    <form onSubmit={handleSubmit} className="space-y-4">
      <div className="flex gap-2">
        <input
          type="text"
          value={smiles}
          onChange={(e) => setSmiles(e.target.value)}
          placeholder="Enter SMILES for process evaluation"
          className="flex-1 rounded-[var(--radius-sm)] border border-border bg-bg-card px-4 py-2.5 text-sm text-text-primary font-mono placeholder:text-text-muted focus:outline-none focus:border-accent focus:ring-1 focus:ring-accent transition-colors"
          autoFocus
        />
        <Button type="submit" disabled={loading || !smiles.trim()}>
          {loading ? (
            <span className="flex items-center gap-2">
              <span className="h-4 w-4 animate-spin rounded-full border-2 border-white border-t-transparent" />
              Evaluating
            </span>
          ) : (
            "Evaluate"
          )}
        </Button>
      </div>
      <div className="flex gap-4 flex-wrap">
        <div className="w-32">
          <Input
            label="Scale (kg)"
            type="number"
            min={0.001}
            step={0.1}
            value={scaleKg}
            onChange={(e) => setScaleKg(Number(e.target.value))}
          />
        </div>
        <div className="w-32">
          <Input
            label="Max Depth"
            type="number"
            min={1}
            max={10}
            value={maxDepth}
            onChange={(e) => setMaxDepth(Number(e.target.value))}
          />
        </div>
        <div className="w-32">
          <Input
            label="Beam Width"
            type="number"
            min={1}
            max={10}
            value={beamWidth}
            onChange={(e) => setBeamWidth(Number(e.target.value))}
          />
        </div>
      </div>
      {error && <Alert variant="error">{error}</Alert>}
    </form>
  );
}
