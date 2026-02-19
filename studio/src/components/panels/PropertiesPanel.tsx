import { PanelContainer } from "./PanelContainer";
import { Badge } from "@/components/ui/Badge";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { useToastStore } from "@/stores/toast-store";

export function PropertiesPanel() {
  const activeId = useWorkspaceStore((s) => s.activeId);
  const molecules = useWorkspaceStore((s) => s.molecules);
  const entry = activeId ? molecules.get(activeId) : undefined;
  const props = entry?.properties;

  if (!props) {
    return (
      <PanelContainer title="Properties">
        <p className="text-xs text-text-muted">No molecule loaded</p>
      </PanelContainer>
    );
  }

  const lipinskiPass = props.lipinski_pass;

  function copyText(text: string, label: string) {
    navigator.clipboard.writeText(text).then(() => {
      useToastStore.getState().addToast(`${label} copied`, "success");
    });
  }

  return (
    <PanelContainer title="Properties">
      <div className="space-y-2">
        {entry?.smiles && (
          <div className="flex items-center gap-1 text-xs">
            <span className="text-text-muted">SMILES:</span>
            <span className="text-text-primary font-mono text-[10px] truncate flex-1">{entry.smiles}</span>
            <button
              onClick={() => copyText(entry.smiles, "SMILES")}
              className="text-[9px] text-text-muted hover:text-accent px-1 cursor-pointer"
              title="Copy SMILES"
            >
              copy
            </button>
          </div>
        )}
        <div className="grid grid-cols-2 gap-x-3 gap-y-1 text-xs">
          <Row label="Formula" value={props.formula} copyable onCopy={() => copyText(props.formula, "Formula")} />
          <Row label="MW" value={`${props.molecular_weight.toFixed(2)} Da`} />
          {props.logp != null && <Row label="LogP" value={props.logp.toFixed(2)} />}
          {props.hbd != null && <Row label="HBD" value={props.hbd.toString()} />}
          {props.hba != null && <Row label="HBA" value={props.hba.toString()} />}
          {props.rotatable_bonds != null && (
            <Row label="Rot. Bonds" value={props.rotatable_bonds.toString()} />
          )}
          {props.tpsa != null && (
            <Row label="TPSA" value={`${props.tpsa.toFixed(1)} A^2`} />
          )}
          {props.heavy_atom_count != null && (
            <Row label="Heavy Atoms" value={props.heavy_atom_count.toString()} />
          )}
        </div>

        {lipinskiPass != null && (
          <div className="flex items-center gap-2">
            <span className="text-xs text-text-muted">Lipinski:</span>
            <Badge variant={lipinskiPass ? "success" : "danger"}>
              {lipinskiPass ? "Pass" : "Fail"}
              {props.lipinski_violations != null && ` (${props.lipinski_violations} violations)`}
            </Badge>
          </div>
        )}

        {props.functional_groups.length > 0 && (
          <div>
            <div className="text-xs text-text-muted mb-1">Functional Groups</div>
            <div className="flex flex-wrap gap-1">
              {props.functional_groups.map((fg) => (
                <Badge key={fg} variant="accent">
                  {fg}
                </Badge>
              ))}
            </div>
          </div>
        )}
      </div>
    </PanelContainer>
  );
}

function Row({ label, value, copyable, onCopy }: { label: string; value: string; copyable?: boolean; onCopy?: () => void }) {
  return (
    <>
      <span className="text-text-muted">{label}</span>
      <span className="text-text-primary font-mono">
        {value}
        {copyable && onCopy && (
          <button onClick={onCopy} className="text-[9px] text-text-muted hover:text-accent ml-1 cursor-pointer">
            copy
          </button>
        )}
      </span>
    </>
  );
}
