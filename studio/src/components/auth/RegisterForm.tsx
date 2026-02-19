import { useState } from "react";
import { Link } from "react-router-dom";
import { register } from "@/api/auth";
import { useAuthStore } from "@/stores/auth-store";
import { Button } from "@/components/ui/Button";
import { Input } from "@/components/ui/Input";
import { Alert } from "@/components/ui/Alert";

export function RegisterForm() {
  const [email, setEmail] = useState("");
  const [apiKey, setApiKey] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);
  const [copied, setCopied] = useState(false);
  const setApiKeyStore = useAuthStore((s) => s.setApiKey);

  async function handleSubmit(e: React.FormEvent) {
    e.preventDefault();
    if (!email.trim()) return;
    setLoading(true);
    setError(null);

    try {
      const res = await register(email.trim());
      setApiKey(res.api_key);
      setApiKeyStore(res.api_key, res.email, res.tier);
    } catch (err) {
      setError(err instanceof Error ? err.message : "Registration failed");
    } finally {
      setLoading(false);
    }
  }

  function handleCopy() {
    if (!apiKey) return;
    navigator.clipboard.writeText(apiKey);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  }

  if (apiKey) {
    return (
      <div className="space-y-4 max-w-sm mx-auto">
        <h1 className="text-xl font-bold text-text-primary">Registration Complete</h1>
        <Alert variant="success">
          Save your API key -- it won't be shown again.
        </Alert>
        <div className="flex items-center gap-2 bg-bg-card border border-border rounded-[var(--radius-sm)] p-3">
          <code className="flex-1 text-xs text-text-primary font-mono break-all">
            {apiKey}
          </code>
          <Button variant="ghost" size="sm" onClick={handleCopy}>
            {copied ? "Copied" : "Copy"}
          </Button>
        </div>
        <Link to="/login">
          <Button className="w-full">Continue to Login</Button>
        </Link>
      </div>
    );
  }

  return (
    <form onSubmit={handleSubmit} className="space-y-4 max-w-sm mx-auto">
      <h1 className="text-xl font-bold text-text-primary">Create Account</h1>
      <p className="text-sm text-text-secondary">
        Register for a free MolBuilder API key to use the studio.
      </p>

      {error && <Alert variant="error">{error}</Alert>}

      <Input
        label="Email"
        type="email"
        placeholder="you@example.com"
        value={email}
        onChange={(e) => setEmail(e.target.value)}
      />

      <Button type="submit" disabled={loading || !email.trim()} className="w-full">
        {loading ? "Registering..." : "Register"}
      </Button>

      <p className="text-xs text-text-muted text-center">
        Already have an API key?{" "}
        <Link to="/login" className="text-accent hover:underline">
          Sign in
        </Link>
      </p>
    </form>
  );
}
