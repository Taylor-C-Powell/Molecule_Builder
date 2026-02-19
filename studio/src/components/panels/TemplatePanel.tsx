import { useState } from "react";
import { PanelContainer } from "./PanelContainer";
import { useMolecule } from "@/hooks/useMolecule";
import { cn } from "@/lib/cn";

interface Template {
  name: string;
  smiles: string;
}

interface TemplateCategory {
  label: string;
  templates: Template[];
}

const TEMPLATE_CATEGORIES: TemplateCategory[] = [
  {
    label: "Common Organics",
    templates: [
      { name: "Benzene", smiles: "c1ccccc1" },
      { name: "Toluene", smiles: "Cc1ccccc1" },
      { name: "Phenol", smiles: "Oc1ccccc1" },
      { name: "Aniline", smiles: "Nc1ccccc1" },
      { name: "Ethanol", smiles: "CCO" },
      { name: "Acetic Acid", smiles: "CC(O)=O" },
      { name: "Acetone", smiles: "CC(C)=O" },
      { name: "Cyclohexane", smiles: "C1CCCCC1" },
    ],
  },
  {
    label: "Heterocycles",
    templates: [
      { name: "Pyridine", smiles: "c1ccncc1" },
      { name: "Pyrimidine", smiles: "c1ccncc1" },
      { name: "Furan", smiles: "c1ccoc1" },
      { name: "Thiophene", smiles: "c1ccsc1" },
      { name: "Pyrrole", smiles: "c1cc[nH]c1" },
      { name: "Imidazole", smiles: "c1cnc[nH]1" },
      { name: "Indole", smiles: "c1ccc2[nH]ccc2c1" },
      { name: "Quinoline", smiles: "c1ccc2ncccc2c1" },
    ],
  },
  {
    label: "Amino Acids",
    templates: [
      { name: "Glycine", smiles: "NCC(O)=O" },
      { name: "Alanine", smiles: "CC(N)C(O)=O" },
      { name: "Valine", smiles: "CC(C)C(N)C(O)=O" },
      { name: "Leucine", smiles: "CC(C)CC(N)C(O)=O" },
      { name: "Phenylalanine", smiles: "NC(Cc1ccccc1)C(O)=O" },
      { name: "Tryptophan", smiles: "NC(Cc1c[nH]c2ccccc12)C(O)=O" },
      { name: "Serine", smiles: "NC(CO)C(O)=O" },
      { name: "Cysteine", smiles: "NC(CS)C(O)=O" },
    ],
  },
  {
    label: "Drug Molecules",
    templates: [
      { name: "Aspirin", smiles: "CC(=O)Oc1ccccc1C(O)=O" },
      { name: "Caffeine", smiles: "Cn1c(=O)c2c(ncn2C)n(C)c1=O" },
      { name: "Ibuprofen", smiles: "CC(C)Cc1ccc(cc1)C(C)C(O)=O" },
      { name: "Paracetamol", smiles: "CC(=O)Nc1ccc(O)cc1" },
      { name: "Naproxen", smiles: "COc1ccc2cc(ccc2c1)C(C)C(O)=O" },
      { name: "Penicillin G", smiles: "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(O)=O" },
    ],
  },
  {
    label: "Functional Groups",
    templates: [
      { name: "Methanol", smiles: "CO" },
      { name: "Formaldehyde", smiles: "C=O" },
      { name: "Methylamine", smiles: "CN" },
      { name: "Formic Acid", smiles: "OC=O" },
      { name: "Methyl Ester", smiles: "COC=O" },
      { name: "Nitromethane", smiles: "C[N+](=O)[O-]" },
    ],
  },
];

export function TemplatePanel() {
  const { parseMolecule, loading } = useMolecule();
  const [activeCategory, setActiveCategory] = useState(0);
  const [loadingTemplate, setLoadingTemplate] = useState<string | null>(null);

  async function handleInsert(template: Template) {
    setLoadingTemplate(template.smiles);
    await parseMolecule(template.smiles);
    setLoadingTemplate(null);
  }

  const category = TEMPLATE_CATEGORIES[activeCategory]!;

  return (
    <PanelContainer title="Templates" defaultOpen={false}>
      {/* Category tabs */}
      <div className="flex flex-wrap gap-1 mb-2">
        {TEMPLATE_CATEGORIES.map((cat, i) => (
          <button
            key={cat.label}
            className={cn(
              "px-2 py-0.5 text-[10px] rounded-[var(--radius-sm)] transition-colors cursor-pointer",
              i === activeCategory
                ? "bg-accent/20 text-accent"
                : "text-text-muted hover:text-text-secondary hover:bg-white/5",
            )}
            onClick={() => setActiveCategory(i)}
          >
            {cat.label}
          </button>
        ))}
      </div>

      {/* Template grid */}
      <div className="grid grid-cols-2 gap-1">
        {category.templates.map((t) => (
          <button
            key={t.smiles}
            className={cn(
              "px-2 py-1.5 text-left rounded-[var(--radius-sm)] transition-colors cursor-pointer border border-transparent hover:border-white/10 hover:bg-white/3",
              loadingTemplate === t.smiles && "opacity-50",
            )}
            disabled={loading}
            onClick={() => handleInsert(t)}
            title={t.smiles}
          >
            <div className="text-[11px] text-text-primary font-medium truncate">
              {t.name}
            </div>
            <div className="text-[9px] text-text-muted font-mono truncate">
              {t.smiles}
            </div>
          </button>
        ))}
      </div>
    </PanelContainer>
  );
}
