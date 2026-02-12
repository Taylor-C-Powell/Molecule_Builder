import { useEffect, useRef } from "react";
import type { Molecule3DResponse } from "@/api/types";
import { getCpkColor } from "@/lib/cpk-colors";
import { Card } from "@/components/ui/Card";
import { Skeleton } from "@/components/ui/Skeleton";

interface MoleculeViewerProps {
  structure?: Molecule3DResponse;
  loading?: boolean;
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any
type Viewer = any;

export function MoleculeViewer({ structure, loading }: MoleculeViewerProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<Viewer>(null);

  useEffect(() => {
    if (!containerRef.current || !structure) return;

    let cancelled = false;
    import("3dmol").then(($3Dmol) => {
      if (cancelled || !containerRef.current) return;

      // Clear previous viewer
      if (viewerRef.current) {
        viewerRef.current.clear();
      } else {
        containerRef.current.innerHTML = "";
        viewerRef.current = $3Dmol.createViewer(containerRef.current, {
          backgroundColor: "#0a0a0a",
          antialias: true,
        });
      }

      const viewer = viewerRef.current;
      if (!viewer) return;

      const model = viewer.addModel();

      for (const atom of structure.atoms) {
        model.addAtom({
          elem: atom.symbol,
          x: atom.position[0],
          y: atom.position[1],
          z: atom.position[2],
          serial: atom.index,
        });
      }

      for (const bond of structure.bonds) {
        model.addBond(bond.atom_i, bond.atom_j, bond.order);
      }

      viewer.setStyle({}, {
        stick: { radius: 0.12, colorscheme: "default" },
        sphere: {
          scale: 0.25,
          colorscheme: {
            prop: "elem",
            map: Object.fromEntries(
              structure.atoms.map((a) => [a.symbol, getCpkColor(a.symbol)]),
            ),
          },
        },
      });

      viewer.zoomTo();
      viewer.render();
      viewer.zoom(0.9);
      viewer.render();
    });

    return () => {
      cancelled = true;
    };
  }, [structure]);

  if (loading) {
    return <Skeleton className="w-full h-[400px]" />;
  }

  if (!structure) {
    return (
      <Card className="flex items-center justify-center h-[400px] text-text-muted text-sm">
        Enter a SMILES string to visualize a molecule in 3D
      </Card>
    );
  }

  return (
    <div
      ref={containerRef}
      className="w-full h-[400px] rounded-[var(--radius-md)] border border-border overflow-hidden"
      style={{ position: "relative" }}
    />
  );
}
