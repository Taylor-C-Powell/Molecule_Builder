import { useState } from "react";
import { RetroForm } from "@/components/retro/RetroForm";
import { RetroTree } from "@/components/retro/RetroTree";
import { RouteSteps } from "@/components/retro/RouteSteps";
import { useRetrosynthesis } from "@/hooks/useRetrosynthesis";
import { useApiClient } from "@/hooks/useApiClient";
import { downloadProcessPdf } from "@/api/process";
import { Badge } from "@/components/ui/Badge";
import { Button } from "@/components/ui/Button";

export default function RetrosynthesisPage() {
  const { plan, loading, error, result } = useRetrosynthesis();
  const client = useApiClient();
  const [pdfLoading, setPdfLoading] = useState(false);

  const handleDownloadPdf = async () => {
    if (!result?.best_route) return;
    setPdfLoading(true);
    try {
      const blob = await downloadProcessPdf(client, result.best_route.target_smiles);
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `retro_report_${result.best_route.target_smiles.slice(0, 20)}.pdf`;
      a.click();
      URL.revokeObjectURL(url);
    } catch {
      // PDF download failed silently
    } finally {
      setPdfLoading(false);
    }
  };

  return (
    <div className="space-y-6">
      <div>
        <h1 className="text-2xl font-bold mb-1">Retrosynthetic Analysis</h1>
        <p className="text-sm text-text-secondary">
          Plan synthesis routes from purchasable starting materials.
        </p>
      </div>

      <RetroForm onSubmit={plan} loading={loading} error={error} />

      {result && (
        <>
          <div className="flex items-center gap-3">
            <Badge variant="accent">{result.routes_found} route{result.routes_found !== 1 ? "s" : ""} found</Badge>
            <Badge variant="default">Depth: {result.max_depth}</Badge>
            <Badge variant="default">Beam: {result.beam_width}</Badge>
            {result.best_route && (
              <Button onClick={handleDownloadPdf} disabled={pdfLoading} variant="secondary" size="sm">
                {pdfLoading ? "Generating..." : "Download PDF"}
              </Button>
            )}
          </div>

          <RetroTree tree={result.tree} />

          {result.best_route && <RouteSteps route={result.best_route} />}
        </>
      )}
    </div>
  );
}
