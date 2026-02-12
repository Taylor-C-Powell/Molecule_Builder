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

/** Convert our API response to XYZ format string for 3Dmol.js */
function toXyz(structure: Molecule3DResponse): string {
  const lines = [
    String(structure.atoms.length),
    structure.id,
    ...structure.atoms.map(
      (a) => `${a.symbol} ${a.position[0]} ${a.position[1]} ${a.position[2]}`,
    ),
  ];
  return lines.join("\n");
}

/** Resolve the $3Dmol namespace from the import or window global */
function resolve3Dmol(mod: Mol3D): Mol3D | null {
  // Try module default export (CJS interop)
  if (mod?.default?.createViewer) return mod.default;
  // Try module directly (named exports)
  if (mod?.createViewer) return mod;
  // Fall back to window global (3Dmol.js sets this as a side effect)
  if (typeof window !== "undefined" && (window as Mol3D).$3Dmol?.createViewer) {
    return (window as Mol3D).$3Dmol;
  }
  return null;
}

export function MoleculeViewer({ structure, loading }: MoleculeViewerProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<Mol3D>(null);
  const [lib, setLib] = useState<Mol3D>(null);

  // Load 3Dmol.js once
  useEffect(() => {
    import("3dmol").then((mod) => {
      const $3Dmol = resolve3Dmol(mod);
      if (!$3Dmol) {
        console.error(
          "3Dmol.js loaded but createViewer not found.",
          "Module keys:", Object.keys(mod),
          "default keys:", mod.default ? Object.keys(mod.default).slice(0, 10) : "none",
          "window.$3Dmol:", typeof (window as Mol3D).$3Dmol,
        );
      }
      setLib($3Dmol);
    });
  }, []);

  // Create viewer once the lib is loaded and the container is mounted
  useEffect(() => {
    if (!lib || !containerRef.current || viewerRef.current) return;
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

    // Use XYZ format - the most reliable path in 3Dmol.js
    const xyz = toXyz(structure);
    viewer.addModel(xyz, "xyz");

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
