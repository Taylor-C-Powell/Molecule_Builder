import { useState, type FormEvent } from "react";
import { Card } from "@/components/ui/Card";
import { Input } from "@/components/ui/Input";
import { Button } from "@/components/ui/Button";
import { Alert } from "@/components/ui/Alert";

function nameToSlug(name: string): string {
  return name
    .toLowerCase()
    .replace(/\s+/g, "-")
    .replace(/[^a-z0-9-]/g, "")
    .replace(/^-+|-+$/g, "")
    .slice(0, 40);
}

export function CreateTeamForm({
  onCreated,
  loading,
  error,
  onCreate,
}: {
  onCreated: () => void;
  loading: boolean;
  error: string | null;
  onCreate: (name: string, slug: string) => Promise<unknown>;
}) {
  const [name, setName] = useState("");
  const [slug, setSlug] = useState("");
  const [autoSlug, setAutoSlug] = useState(true);

  const handleNameChange = (val: string) => {
    setName(val);
    if (autoSlug) {
      setSlug(nameToSlug(val));
    }
  };

  const handleSlugChange = (val: string) => {
    setAutoSlug(false);
    setSlug(val);
  };

  const handleSubmit = async (e: FormEvent) => {
    e.preventDefault();
    const result = await onCreate(name, slug);
    if (result) {
      setName("");
      setSlug("");
      setAutoSlug(true);
      onCreated();
    }
  };

  return (
    <Card>
      <h2 className="text-lg font-semibold mb-4">Create Team</h2>
      <form onSubmit={handleSubmit} className="space-y-3">
        <Input
          label="Team Name"
          value={name}
          onChange={(e) => handleNameChange(e.target.value)}
          placeholder="e.g. Discovery Chemistry"
          required
        />
        <Input
          label="Slug"
          value={slug}
          onChange={(e) => handleSlugChange(e.target.value)}
          placeholder="e.g. discovery-chemistry"
          required
        />
        <p className="text-xs text-text-muted">
          3-40 characters, lowercase letters, numbers, and hyphens only.
        </p>
        {error && <Alert variant="error">{error}</Alert>}
        <Button type="submit" disabled={loading || !name || !slug}>
          {loading ? "Creating..." : "Create Team"}
        </Button>
      </form>
    </Card>
  );
}
