import { useState, useCallback } from "react";
import { Button } from "@/components/ui/Button";
import { Alert } from "@/components/ui/Alert";

interface SmilesInputProps {
  onSubmit: (smiles: string) => void;
  loading?: boolean;
  error?: string | null;
}

const EXAMPLES = [
  { label: "Ethanol", smiles: "CCO" },
  { label: "Caffeine", smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" },
  { label: "Aspirin", smiles: "CC(=O)OC1=CC=CC=C1C(=O)O" },
  { label: "Benzene", smiles: "c1ccccc1" },
  { label: "Ibuprofen", smiles: "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O" },
];

export function SmilesInput({ onSubmit, loading, error }: SmilesInputProps) {
  const [smiles, setSmiles] = useState("");

  const handleSubmit = useCallback(
    (e: React.FormEvent) => {
      e.preventDefault();
      if (smiles.trim()) onSubmit(smiles.trim());
    },
    [smiles, onSubmit],
  );

  return (
    <div className="space-y-3">
      <form onSubmit={handleSubmit} className="flex gap-2">
        <input
          type="text"
          value={smiles}
          onChange={(e) => setSmiles(e.target.value)}
          placeholder="Enter SMILES (e.g., CCO for ethanol)"
          className="flex-1 rounded-[var(--radius-sm)] border border-border bg-bg-card px-4 py-2.5 text-sm text-text-primary font-mono placeholder:text-text-muted focus:outline-none focus:border-accent focus:ring-1 focus:ring-accent transition-colors"
          autoFocus
        />
        <Button type="submit" disabled={loading || !smiles.trim()}>
          {loading ? (
            <span className="flex items-center gap-2">
              <span className="h-4 w-4 animate-spin rounded-full border-2 border-white border-t-transparent" />
              Parsing
            </span>
          ) : (
            "Parse"
          )}
        </Button>
      </form>
      {error && <Alert variant="error">{error}</Alert>}
      <div className="flex flex-wrap gap-2">
        <span className="text-xs text-text-muted pt-1">Examples:</span>
        {EXAMPLES.map((ex) => (
          <button
            key={ex.smiles}
            type="button"
            onClick={() => {
              setSmiles(ex.smiles);
              onSubmit(ex.smiles);
            }}
            className="px-2.5 py-1 text-xs rounded-full border border-border bg-bg-card text-text-secondary hover:border-border-hover hover:text-text-primary transition-colors cursor-pointer"
          >
            {ex.label}
          </button>
        ))}
      </div>
    </div>
  );
}
