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

/** Convert our API response to SDF V2000 molfile for 3Dmol.js.
 *  This preserves exact bond topology and bond orders from the API
 *  rather than relying on 3Dmol.js distance-based bond inference. */
function toSdf(structure: Molecule3DResponse): string {
  const nAtoms = structure.atoms.length;
  const nBonds = structure.bonds.length;
  const lines: string[] = [
    structure.id,              // molecule name
    "  MolBuilder  3D",        // program/timestamp
    "",                        // comment
    // counts line: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
    `${String(nAtoms).padStart(3)}${String(nBonds).padStart(3)}  0  0  0  0  0  0  0  0999 V2000`,
  ];
  // Atom block: x(10.4) y(10.4) z(10.4) symbol(3) ...
  for (const a of structure.atoms) {
    const x = a.position[0].toFixed(4).padStart(10);
    const y = a.position[1].toFixed(4).padStart(10);
    const z = a.position[2].toFixed(4).padStart(10);
    const sym = ` ${a.symbol.padEnd(2)}`;
    lines.push(`${x}${y}${z}${sym} 0  0  0  0  0  0  0  0  0  0  0  0`);
  }
  // Bond block: iii(3) jjj(3) ttt(3) ...  (1-indexed atoms)
  for (const b of structure.bonds) {
    const i = String(b.atom_i + 1).padStart(3);
    const j = String(b.atom_j + 1).padStart(3);
    const t = String(Math.min(b.order, 3)).padStart(3);
    lines.push(`${i}${j}${t}  0  0  0  0`);
  }
  lines.push("M  END", "$$$$");
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

    // Use SDF V2000 format to preserve exact bond topology and orders
    const sdf = toSdf(structure);
    viewer.addModel(sdf, "sdf");

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
