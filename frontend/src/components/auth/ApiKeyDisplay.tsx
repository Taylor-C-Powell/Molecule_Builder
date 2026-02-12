import { useState } from "react";
import { Button } from "@/components/ui/Button";

interface ApiKeyDisplayProps {
  apiKey: string;
}

export function ApiKeyDisplay({ apiKey }: ApiKeyDisplayProps) {
  const [copied, setCopied] = useState(false);

  async function copy() {
    await navigator.clipboard.writeText(apiKey);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  }

  return (
    <div className="flex items-center gap-2 rounded-[var(--radius-sm)] border border-border bg-bg p-3 font-mono text-sm">
      <code className="flex-1 break-all text-accent">{apiKey}</code>
      <Button variant="secondary" size="sm" onClick={copy}>
        {copied ? "Copied" : "Copy"}
      </Button>
    </div>
  );
}
