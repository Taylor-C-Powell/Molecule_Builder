import { SmilesInput } from "@/components/molecule/SmilesInput";
import { MoleculeViewer } from "@/components/molecule/MoleculeViewer";
import { PropertiesPanel } from "@/components/molecule/PropertiesPanel";
import { useMolecule } from "@/hooks/useMolecule";

export default function DashboardPage() {
  const { parseMolecule, loading, error, currentEntry } = useMolecule();

  return (
    <div className="space-y-6">
      <div>
        <h1 className="text-2xl font-bold mb-1">Molecule Dashboard</h1>
        <p className="text-sm text-text-secondary">
          Enter a SMILES string to parse, visualize, and analyze a molecule.
        </p>
      </div>

      <SmilesInput onSubmit={parseMolecule} loading={loading} error={error} />

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        <MoleculeViewer
          structure={currentEntry?.structure}
          loading={loading}
        />
        <PropertiesPanel
          properties={currentEntry?.properties}
          loading={loading}
        />
      </div>
    </div>
  );
}
