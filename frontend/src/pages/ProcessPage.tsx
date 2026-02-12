import { ProcessForm } from "@/components/process/ProcessForm";
import { ProcessOverview } from "@/components/process/ProcessOverview";
import { ReactorCard } from "@/components/process/ReactorCard";
import { SafetyPanel } from "@/components/process/SafetyPanel";
import { CostBreakdown } from "@/components/process/CostBreakdown";
import { ScaleUpCard } from "@/components/process/ScaleUpCard";
import { useProcess } from "@/hooks/useProcess";

export default function ProcessPage() {
  const { evaluate, loading, error, result } = useProcess();

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
          <ProcessOverview data={result} />

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
