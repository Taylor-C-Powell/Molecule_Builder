import { useMolecule } from "@/hooks/useMolecule";

const QUICK_MOLECULES = [
  { label: "Benzene", smiles: "c1ccccc1", desc: "Aromatic ring" },
  { label: "Aspirin", smiles: "CC(=O)Oc1ccccc1C(=O)O", desc: "Anti-inflammatory" },
  { label: "Caffeine", smiles: "Cn1c(=O)c2c(ncn2C)n(C)c1=O", desc: "Stimulant" },
  { label: "Ethanol", smiles: "CCO", desc: "Simple alcohol" },
  { label: "Glucose", smiles: "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", desc: "Sugar" },
  { label: "Ibuprofen", smiles: "CC(C)Cc1ccc(cc1)C(C)C(=O)O", desc: "NSAID" },
];

export function WelcomeScreen() {
  const { parseMolecule, loading } = useMolecule();

  return (
    <div className="absolute inset-0 flex items-center justify-center pointer-events-none">
      <div className="text-center space-y-6 max-w-md pointer-events-auto">
        <div className="space-y-1">
          <h2 className="text-lg font-bold text-text-primary tracking-tight">
            MolBuilder Studio
          </h2>
          <p className="text-xs text-text-muted">
            3D molecule building and analysis
          </p>
        </div>

        <div className="space-y-3">
          <p className="text-[11px] text-text-secondary">
            Quick start -- click a molecule:
          </p>
          <div className="grid grid-cols-3 gap-2">
            {QUICK_MOLECULES.map((m) => (
              <button
                key={m.smiles}
                disabled={loading}
                onClick={() => parseMolecule(m.smiles, m.label)}
                className="px-3 py-2 bg-bg-card border border-border rounded-[var(--radius-sm)] hover:border-accent/50 hover:bg-accent/5 transition-colors cursor-pointer text-left group"
              >
                <div className="text-xs font-semibold text-text-primary group-hover:text-accent transition-colors">
                  {m.label}
                </div>
                <div className="text-[10px] text-text-muted">{m.desc}</div>
              </button>
            ))}
          </div>
        </div>

        <div className="flex items-center gap-3 justify-center text-[10px] text-text-muted">
          <span>Type SMILES in sidebar</span>
          <span className="text-border">|</span>
          <span>Import .mol / .sdf / .xyz</span>
          <span className="text-border">|</span>
          <span>Use a template</span>
        </div>

        <div className="text-[10px] text-text-muted/60">
          Press <kbd className="px-1 py-0.5 bg-white/5 border border-border rounded text-[9px] font-mono">?</kbd> for keyboard shortcuts
        </div>
      </div>
    </div>
  );
}
