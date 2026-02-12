import { RetroForm } from "@/components/retro/RetroForm";
import { RetroTree } from "@/components/retro/RetroTree";
import { RouteSteps } from "@/components/retro/RouteSteps";
import { useRetrosynthesis } from "@/hooks/useRetrosynthesis";
import { Badge } from "@/components/ui/Badge";

export default function RetrosynthesisPage() {
  const { plan, loading, error, result } = useRetrosynthesis();

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
          </div>

          <RetroTree tree={result.tree} />

          {result.best_route && <RouteSteps route={result.best_route} />}
        </>
      )}
    </div>
  );
}
