import { useState, useEffect, type FormEvent } from "react";
import { useBatch } from "@/hooks/useBatch";
import { Button } from "@/components/ui/Button";
import { Card } from "@/components/ui/Card";
import { Badge } from "@/components/ui/Badge";
import { Alert } from "@/components/ui/Alert";
import type { BatchJobSummary } from "@/api/types";

const JOB_TYPES = ["properties", "retrosynthesis", "conditions", "evaluate"] as const;

function statusVariant(status: string): "accent" | "default" {
  if (status === "completed") return "accent";
  return "default";
}

function SubmitForm({ onSubmitted }: { onSubmitted: () => void }) {
  const { submit, loading, error } = useBatch();
  const [text, setText] = useState("");
  const [jobType, setJobType] = useState<(typeof JOB_TYPES)[number]>("properties");
  const [paramsStr, setParamsStr] = useState("");

  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();
    const list = text
      .split("\n")
      .map((s) => s.trim())
      .filter(Boolean);
    if (!list.length) return;

    let params: Record<string, unknown> | undefined;
    if (paramsStr.trim()) {
      try {
        params = JSON.parse(paramsStr);
      } catch {
        return;
      }
    }

    const res = await submit(list, jobType, params);
    if (res) {
      setText("");
      setParamsStr("");
      onSubmitted();
    }
  };

  return (
    <Card>
      <h2 className="text-lg font-semibold mb-4">Submit Batch Job</h2>
      <form onSubmit={handleSubmit} className="space-y-3">
        <div className="flex flex-col gap-1.5">
          <label className="text-sm font-medium text-text-secondary">SMILES (one per line)</label>
          <textarea
            value={text}
            onChange={(e) => setText(e.target.value)}
            placeholder={"CCO\nCC(=O)O\nc1ccccc1"}
            rows={5}
            className="w-full rounded-[var(--radius-sm)] border border-border bg-bg-card px-3 py-2 text-sm text-text-primary placeholder:text-text-muted focus:outline-none focus:border-accent focus:ring-1 focus:ring-accent transition-colors resize-y font-mono"
          />
        </div>

        <div className="flex flex-col gap-1.5">
          <label className="text-sm font-medium text-text-secondary">Job Type</label>
          <select
            value={jobType}
            onChange={(e) => setJobType(e.target.value as (typeof JOB_TYPES)[number])}
            className="w-full rounded-[var(--radius-sm)] border border-border bg-bg-card px-3 py-2 text-sm text-text-primary focus:outline-none focus:border-accent focus:ring-1 focus:ring-accent transition-colors"
          >
            {JOB_TYPES.map((t) => (
              <option key={t} value={t}>
                {t.charAt(0).toUpperCase() + t.slice(1)}
              </option>
            ))}
          </select>
        </div>

        <div className="flex flex-col gap-1.5">
          <label className="text-sm font-medium text-text-secondary">Parameters (JSON, optional)</label>
          <textarea
            value={paramsStr}
            onChange={(e) => setParamsStr(e.target.value)}
            placeholder='{"max_depth": 3}'
            rows={2}
            className="w-full rounded-[var(--radius-sm)] border border-border bg-bg-card px-3 py-2 text-sm text-text-primary placeholder:text-text-muted focus:outline-none focus:border-accent focus:ring-1 focus:ring-accent transition-colors resize-y font-mono"
          />
        </div>

        {error && <Alert variant="error">{error}</Alert>}
        <Button type="submit" disabled={loading || !text.trim()}>
          {loading ? "Submitting..." : "Submit Batch"}
        </Button>
      </form>
    </Card>
  );
}

function ProgressBar({ pct }: { pct: number }) {
  return (
    <div className="w-full bg-border rounded-full h-2 overflow-hidden">
      <div
        className="bg-accent h-full rounded-full transition-all duration-300"
        style={{ width: `${Math.min(100, Math.max(0, pct))}%` }}
      />
    </div>
  );
}

function JobRow({
  job,
  onCancel,
  onSelect,
}: {
  job: BatchJobSummary;
  onCancel: (jobId: string) => void;
  onSelect: (jobId: string) => void;
}) {
  const canCancel = job.status === "pending" || job.status === "running";

  return (
    <div className="border border-border rounded-[var(--radius-sm)] bg-bg-card px-4 py-3 space-y-2">
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-3">
          <code className="text-xs font-mono text-text-muted">{job.job_id.slice(0, 12)}</code>
          <Badge variant={statusVariant(job.status)}>{job.status}</Badge>
          <Badge variant="default">{job.job_type}</Badge>
        </div>
        <div className="flex items-center gap-2">
          <button
            onClick={() => onSelect(job.job_id)}
            className="text-xs text-accent hover:underline bg-transparent border-none cursor-pointer"
          >
            Details
          </button>
          {canCancel && (
            <Button size="sm" variant="ghost" onClick={() => onCancel(job.job_id)}>
              Cancel
            </Button>
          )}
        </div>
      </div>
      {(job.status === "running" || job.status === "pending") && (
        <ProgressBar pct={job.progress_pct} />
      )}
      <div className="text-xs text-text-muted">
        Created: {new Date(job.created_at).toLocaleString()}
      </div>
    </div>
  );
}

export default function BatchPage() {
  const { fetchList, cancel, pollStatus, stopPolling, loading, error, jobs, jobStatus } = useBatch();

  const reload = () => {
    fetchList();
  };

  useEffect(() => {
    reload();
    const interval = setInterval(reload, 5000);
    return () => clearInterval(interval);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  const handleSelect = (jobId: string) => {
    pollStatus(jobId);
  };

  const handleCancel = async (jobId: string) => {
    const ok = await cancel(jobId);
    if (ok) reload();
  };

  return (
    <div className="space-y-6">
      <div>
        <h1 className="text-2xl font-bold mb-1">Batch Processing</h1>
        <p className="text-sm text-text-secondary">
          Submit bulk SMILES for parallel analysis.
        </p>
      </div>

      <SubmitForm onSubmitted={reload} />

      {error && <Alert variant="error">{error}</Alert>}

      <Card>
        <div className="flex items-center justify-between mb-4">
          <h2 className="text-lg font-semibold">Jobs</h2>
          <Button size="sm" variant="ghost" onClick={reload} disabled={loading}>
            Refresh
          </Button>
        </div>

        {jobs && jobs.jobs.length === 0 && (
          <p className="text-sm text-text-muted">No batch jobs yet.</p>
        )}

        {jobs && jobs.jobs.length > 0 && (
          <div className="space-y-2">
            {jobs.jobs.map((job) => (
              <JobRow key={job.job_id} job={job} onCancel={handleCancel} onSelect={handleSelect} />
            ))}
          </div>
        )}
      </Card>

      {jobStatus && (
        <Card>
          <h2 className="text-lg font-semibold mb-3">Job Details</h2>
          <div className="space-y-2 text-sm">
            <div>
              <span className="text-text-muted">ID: </span>
              <code className="font-mono">{jobStatus.job_id}</code>
            </div>
            <div>
              <span className="text-text-muted">Status: </span>
              <Badge variant={statusVariant(jobStatus.status)}>{jobStatus.status}</Badge>
            </div>
            <div>
              <span className="text-text-muted">Type: </span>
              <span>{jobStatus.job_type}</span>
            </div>
            <div>
              <span className="text-text-muted">Progress: </span>
              <span>{jobStatus.progress_pct.toFixed(1)}%</span>
            </div>
            <ProgressBar pct={jobStatus.progress_pct} />

            {jobStatus.error && (
              <Alert variant="error">{jobStatus.error}</Alert>
            )}

            {jobStatus.result && (
              <div>
                <span className="text-text-muted block mb-1">Result:</span>
                <pre className="bg-bg rounded-[var(--radius-sm)] border border-border p-3 text-xs font-mono overflow-x-auto max-h-96">
                  {JSON.stringify(jobStatus.result, null, 2)}
                </pre>
              </div>
            )}

            {(jobStatus.status === "pending" || jobStatus.status === "running") && (
              <div className="flex gap-2">
                <Button size="sm" variant="ghost" onClick={stopPolling}>Stop Polling</Button>
                <Button size="sm" variant="ghost" onClick={() => handleCancel(jobStatus.job_id)}>Cancel Job</Button>
              </div>
            )}
          </div>
        </Card>
      )}
    </div>
  );
}
