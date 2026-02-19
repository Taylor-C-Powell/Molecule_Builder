import { useState } from "react";
import { PanelContainer } from "./PanelContainer";
import { Button } from "@/components/ui/Button";
import { Alert } from "@/components/ui/Alert";
import { useMolecule } from "@/hooks/useMolecule";

const EXAMPLES = [
  { label: "Benzene", smiles: "c1ccccc1" },
  { label: "Aspirin", smiles: "CC(=O)Oc1ccccc1C(=O)O" },
  { label: "Caffeine", smiles: "Cn1c(=O)c2c(ncn2C)n(C)c1=O" },
  { label: "Ethanol", smiles: "CCO" },
  { label: "Ibuprofen", smiles: "CC(C)Cc1ccc(cc1)C(C)C(=O)O" },
  { label: "Paracetamol", smiles: "CC(=O)Nc1ccc(O)cc1" },
];

export function SmilesPanel() {
  const [smiles, setSmiles] = useState("");
  const { parseMolecule, loading, error } = useMolecule();

  function handleSubmit(e: React.FormEvent) {
    e.preventDefault();
    if (smiles.trim()) {
      parseMolecule(smiles.trim());
    }
  }

  function handleExample(s: string) {
    setSmiles(s);
    parseMolecule(s);
  }

  return (
    <PanelContainer title="SMILES Input">
      <form onSubmit={handleSubmit} className="space-y-2">
        <div className="flex gap-1">
          <input
            className="flex-1 px-2 py-1.5 text-xs font-mono bg-bg border border-border rounded-[var(--radius-sm)] text-text-primary outline-none focus:border-accent placeholder:text-text-muted"
            placeholder="Enter SMILES..."
            value={smiles}
            onChange={(e) => setSmiles(e.target.value)}
          />
          <Button type="submit" size="sm" disabled={loading || !smiles.trim()}>
            {loading ? "..." : "Parse"}
          </Button>
        </div>
      </form>

      {error && (
        <Alert variant="error" className="mt-2 text-xs">
          {error}
        </Alert>
      )}

      <div className="mt-2 flex flex-wrap gap-1">
        {EXAMPLES.map((ex) => (
          <button
            key={ex.smiles}
            className="px-2 py-0.5 text-[10px] bg-bg border border-border rounded-full text-text-muted hover:text-text-secondary hover:border-border-hover transition-colors cursor-pointer"
            onClick={() => handleExample(ex.smiles)}
          >
            {ex.label}
          </button>
        ))}
      </div>
    </PanelContainer>
  );
}
