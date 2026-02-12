import { useEffect, useRef, useState } from "react";
import type { Molecule3DResponse } from "@/api/types";
import { getCpkColor } from "@/lib/cpk-colors";
import { Skeleton } from "@/components/ui/Skeleton";

interface MoleculeViewerProps {
  structure?: Molecule3DResponse;
  loading?: boolean;
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any
type Mol3D = any;

export function MoleculeViewer({ structure, loading }: MoleculeViewerProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<Mol3D>(null);
  const [lib, setLib] = useState<Mol3D>(null);

  // Load 3Dmol.js once - handle both ESM default and named export patterns
  useEffect(() => {
    import("3dmol").then((mod) => {
      const $3Dmol = mod.default ?? mod;
      setLib($3Dmol);
    });
  }, []);

  // Create viewer once the lib is loaded and the container is mounted
  useEffect(() => {
    if (!lib || !containerRef.current || viewerRef.current) return;
    if (typeof lib.createViewer !== "function") {
      console.error("3Dmol.createViewer not found. Module keys:", Object.keys(lib));
      return;
    }
    containerRef.current.innerHTML = "";
    viewerRef.current = lib.createViewer(containerRef.current, {
      backgroundColor: "#0a0a0a",
      antialias: true,
    });
  }, [lib]);

  // Render molecule when structure changes
  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer || !structure) return;

    viewer.removeAllModels();

    // Build atom specs with bond connectivity (3Dmol.js API)
    // Each atom needs: elem, x, y, z, bonds (neighbor indices), bondOrder
    const atoms = structure.atoms.map((a) => ({
      elem: a.symbol,
      x: a.position[0],
      y: a.position[1],
      z: a.position[2],
      serial: a.index,
      bonds: [] as number[],
      bondOrder: [] as number[],
    }));

    for (const bond of structure.bonds) {
      const a = atoms[bond.atom_i];
      const b = atoms[bond.atom_j];
      if (a && b) {
        a.bonds.push(bond.atom_j);
        a.bondOrder.push(bond.order);
        b.bonds.push(bond.atom_i);
        b.bondOrder.push(bond.order);
      }
    }

    const model = viewer.addModel();
    model.addAtoms(atoms);

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
  }, [structure, lib]);

  if (loading && !structure) {
    return <Skeleton className="w-full h-[400px]" />;
  }

  return (
    <div className="relative w-full h-[400px] rounded-[var(--radius-md)] border border-border overflow-hidden bg-bg">
      <div
        ref={containerRef}
        className="absolute inset-0"
      />
      {!structure && !loading && (
        <div className="absolute inset-0 flex items-center justify-center text-text-muted text-sm pointer-events-none">
          Enter a SMILES string to visualize a molecule in 3D
        </div>
      )}
    </div>
  );
}
