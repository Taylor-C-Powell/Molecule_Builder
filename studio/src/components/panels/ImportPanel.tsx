import { useRef, useState, useCallback } from "react";
import { PanelContainer } from "./PanelContainer";
import { useImport } from "@/hooks/useImport";
import { cn } from "@/lib/cn";

export function ImportPanel() {
  const { importFile } = useImport();
  const fileInputRef = useRef<HTMLInputElement>(null);
  const [dragOver, setDragOver] = useState(false);

  const handleFile = useCallback(
    (file: File) => {
      importFile(file);
    },
    [importFile],
  );

  function handleDrop(e: React.DragEvent) {
    e.preventDefault();
    setDragOver(false);
    const file = e.dataTransfer.files[0];
    if (file) handleFile(file);
  }

  function handleDragOver(e: React.DragEvent) {
    e.preventDefault();
    setDragOver(true);
  }

  function handleDragLeave() {
    setDragOver(false);
  }

  function handleInputChange(e: React.ChangeEvent<HTMLInputElement>) {
    const file = e.target.files?.[0];
    if (file) handleFile(file);
    // Reset input so the same file can be re-imported
    e.target.value = "";
  }

  return (
    <PanelContainer title="Import" defaultOpen={false}>
      <div
        className={cn(
          "border-2 border-dashed rounded-[var(--radius-sm)] p-4 text-center transition-colors cursor-pointer",
          dragOver
            ? "border-accent bg-accent/10"
            : "border-border hover:border-white/20 hover:bg-white/3",
        )}
        onDrop={handleDrop}
        onDragOver={handleDragOver}
        onDragLeave={handleDragLeave}
        onClick={() => fileInputRef.current?.click()}
      >
        <div className="text-xs text-text-muted mb-1">
          Drop file here or click to browse
        </div>
        <div className="text-[9px] text-text-muted/60">
          Supports .mol, .sdf, .xyz
        </div>
      </div>
      <input
        ref={fileInputRef}
        type="file"
        accept=".mol,.sdf,.sd,.xyz"
        className="hidden"
        onChange={handleInputChange}
      />
    </PanelContainer>
  );
}
