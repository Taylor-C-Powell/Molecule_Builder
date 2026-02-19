import { useMemo } from "react";
import { useWorkspaceStore } from "@/stores/workspace-store";
import { useEditorStore } from "@/stores/editor-store";
import { AtomMesh } from "./AtomMesh";
import { BondMesh } from "./BondMesh";
import { AddAtomPlane } from "./AddAtomPlane";
import { BondPreview } from "./BondPreview";
import { MeasurementOverlay } from "./MeasurementOverlay";
import { DistanceLabels } from "./DistanceLabels";
import { InstancedAtoms } from "./InstancedAtoms";

const INSTANCED_THRESHOLD = 200;

export function MoleculeScene() {
  const activeId = useWorkspaceStore((s) => s.activeId);
  const molecules = useWorkspaceStore((s) => s.molecules);
  const selectedAtoms = useEditorStore((s) => s.selectedAtoms);
  const showLabels = useEditorStore((s) => s.showLabels);
  const showHydrogens = useEditorStore((s) => s.showHydrogens);
  const renderStyle = useEditorStore((s) => s.renderStyle);

  const entry = activeId ? molecules.get(activeId) : undefined;
  const structure = entry?.structure;

  const centroid = useMemo(() => {
    if (!structure || structure.atoms.length === 0) return [0, 0, 0] as const;
    let cx = 0, cy = 0, cz = 0;
    for (const a of structure.atoms) {
      cx += a.position[0];
      cy += a.position[1];
      cz += a.position[2];
    }
    const n = structure.atoms.length;
    return [cx / n, cy / n, cz / n] as const;
  }, [structure]);

  if (!structure) return null;

  const visibleAtoms = showHydrogens
    ? structure.atoms
    : structure.atoms.filter((a) => a.symbol !== "H");
  const visibleIndices = new Set(visibleAtoms.map((a) => a.index));
  const visibleBonds = structure.bonds.filter(
    (b) => visibleIndices.has(b.atom_i) && visibleIndices.has(b.atom_j),
  );

  // Use instanced rendering for large molecules (background decorative layer)
  const useInstanced = visibleAtoms.length > INSTANCED_THRESHOLD;

  return (
    <>
      <group position={[-centroid[0], -centroid[1], -centroid[2]]}>
        {/* Instanced background atoms for large molecules */}
        {useInstanced && (
          <InstancedAtoms atoms={visibleAtoms} renderStyle={renderStyle} />
        )}

        {/* Interactive atom meshes (always rendered for click handling) */}
        {visibleAtoms.map((atom) => (
          <AtomMesh
            key={atom.index}
            atom={atom}
            selected={selectedAtoms.has(atom.index)}
            showLabel={showLabels}
            renderStyle={renderStyle}
            centroid={centroid}
          />
        ))}
        {visibleBonds.map((bond) => (
          <BondMesh
            key={`${bond.atom_i}-${bond.atom_j}`}
            bond={bond}
            atoms={structure.atoms}
            renderStyle={renderStyle}
          />
        ))}
        <BondPreview />
        <MeasurementOverlay />
        <DistanceLabels />
        <AddAtomPlane />
      </group>
    </>
  );
}
