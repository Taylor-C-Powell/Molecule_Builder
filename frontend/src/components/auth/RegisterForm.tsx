import { useState } from "react";
import { useNavigate } from "react-router-dom";
import { register } from "@/api/auth";
import { useAuthStore } from "@/stores/auth-store";
import { Button } from "@/components/ui/Button";
import { Input } from "@/components/ui/Input";
import { Card, CardHeader, CardTitle } from "@/components/ui/Card";
import { Alert } from "@/components/ui/Alert";
import { ApiKeyDisplay } from "./ApiKeyDisplay";

export function RegisterForm() {
  const [email, setEmail] = useState("");
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [apiKey, setApiKey] = useState<string | null>(null);
  const setApiKeyStore = useAuthStore((s) => s.setApiKey);
  const navigate = useNavigate();

  async function handleSubmit(e: React.FormEvent) {
    e.preventDefault();
    setLoading(true);
    setError(null);
    try {
      const res = await register(email);
      setApiKey(res.api_key);
      setApiKeyStore(res.api_key, res.email, res.tier);
    } catch (err) {
      setError(err instanceof Error ? err.message : "Registration failed");
    } finally {
      setLoading(false);
    }
  }

  if (apiKey) {
    return (
      <Card className="max-w-md mx-auto mt-12">
        <CardHeader>
          <CardTitle>Your API Key</CardTitle>
        </CardHeader>
        <Alert variant="warning" className="mb-4">
          Save this key now -- it will not be shown again.
        </Alert>
        <ApiKeyDisplay apiKey={apiKey} />
        <Button className="w-full mt-4" onClick={() => navigate("/dashboard")}>
          Go to Dashboard
        </Button>
      </Card>
    );
  }

  return (
    <Card className="max-w-md mx-auto mt-12">
      <CardHeader>
        <CardTitle>Get Started Free</CardTitle>
        <p className="text-sm text-text-secondary mt-1">
          Create an account to access the MolBuilder API.
        </p>
      </CardHeader>
      <form onSubmit={handleSubmit} className="flex flex-col gap-4">
        {error && <Alert variant="error">{error}</Alert>}
        <Input
          label="Email"
          type="email"
          value={email}
          onChange={(e) => setEmail(e.target.value)}
          placeholder="you@example.com"
          required
          autoFocus
        />
        <Button type="submit" disabled={loading || !email} className="w-full">
          {loading ? "Creating account..." : "Get API Key"}
        </Button>
        <p className="text-xs text-text-muted text-center">
          Free tier: 10 requests/min. No credit card required.
        </p>
      </form>
    </Card>
  );
}
