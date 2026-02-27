import { useState } from "react";
import { ProcessForm } from "@/components/process/ProcessForm";
import { ProcessOverview } from "@/components/process/ProcessOverview";
import { ReactorCard } from "@/components/process/ReactorCard";
import { SafetyPanel } from "@/components/process/SafetyPanel";
import { CostBreakdown } from "@/components/process/CostBreakdown";
import { ScaleUpCard } from "@/components/process/ScaleUpCard";
import { useProcess } from "@/hooks/useProcess";
import { useApiClient } from "@/hooks/useApiClient";
import { downloadProcessPdf } from "@/api/process";
import { Button } from "@/components/ui/Button";

export default function ProcessPage() {
  const { evaluate, loading, error, result } = useProcess();
  const client = useApiClient();
  const [pdfLoading, setPdfLoading] = useState(false);

  const handleDownloadPdf = async () => {
    if (!result) return;
    setPdfLoading(true);
    try {
      const blob = await downloadProcessPdf(client, result.smiles, result.scale_kg);
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `process_report_${result.smiles.slice(0, 20)}.pdf`;
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
        <h1 className="text-2xl font-bold mb-1">Process Engineering</h1>
        <p className="text-sm text-text-secondary">
          Evaluate manufacturing processes: reactor sizing, safety, cost, and scale-up.
        </p>
      </div>

      <ProcessForm onSubmit={evaluate} loading={loading} error={error} />

      {result && (
        <>
          <div className="flex items-center justify-between">
            <ProcessOverview data={result} />
            <Button onClick={handleDownloadPdf} disabled={pdfLoading} variant="secondary">
              {pdfLoading ? "Generating..." : "Download PDF Report"}
            </Button>
          </div>

          {result.step_details.map((step) => (
            <ReactorCard key={step.step_number} step={step} />
          ))}

          {result.safety.length > 0 && <SafetyPanel safety={result.safety} />}

          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
            {result.cost && <CostBreakdown cost={result.cost} />}
            {result.scale_up && <ScaleUpCard scaleUp={result.scale_up} />}
          </div>
        </>
      )}
    </div>
  );
}
